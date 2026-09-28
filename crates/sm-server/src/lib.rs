//! The sentromap HTTP server (design §12): one process per index, memory-mapped and
//! pre-faulted; CPU work on the blocking pool and rayon, never on async worker threads.

pub mod api;
pub mod state;

use std::net::SocketAddr;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use anyhow::Result;
use axum::Router;
use axum::routing::{delete, get, post};
use tower_http::compression::CompressionLayer;
use tower_http::cors::CorsLayer;
use tower_http::services::ServeDir;
use tower_http::trace::TraceLayer;
use tracing::info;

pub use state::App;

pub struct Config {
    pub index: PathBuf,
    pub addr: SocketAddr,
    /// Directory of static front-end files served at `/`.
    pub web: Option<PathBuf>,
    pub cache_bytes: usize,
    pub populate: bool,
    pub cors: bool,
    pub calibrate: bool,
}

/// Open the index and annotations, calibrate the cost model.
pub fn load(cfg: &Config) -> Result<App> {
    let t = Instant::now();
    let idx = sm_index::Index::open_with(&cfg.index, sm_index::array::MapOptions { populate: cfg.populate })?;
    let annots = sm_annot::Annotations::open(&cfg.index)?;
    info!(
        secs = t.elapsed().as_secs_f64(),
        tracks = annots.tracks.len(),
        feature_sets = annots.feature_sets.len(),
        "index loaded"
    );
    let model = if cfg.calibrate {
        let t = Instant::now();
        let m = sm_search::CostModel::calibrate(&idx);
        info!(secs = t.elapsed().as_secs_f64(), ?m, "cost model calibrated");
        m
    } else {
        sm_search::CostModel::default()
    };
    Ok(App::new(idx, annots, model, cfg.cache_bytes))
}

pub fn router(app: Arc<App>, web: Option<PathBuf>, cors: bool) -> Router {
    use api::*;
    let api = Router::new()
        .route("/health", get(health))
        .route("/genome", get(genome))
        .route("/tracks", get(tracks))
        .route("/tracks/{name}", get(track_summary))
        .route("/features/{name}", get(features))
        .route("/sequence", get(sequence))
        .route("/kmer", get(kmer_at))
        .route("/queries", post(create_query))
        .route("/queries/{id}", get(query_status))
        .route("/queries/{id}", delete(cancel_query))
        .route("/queries/{id}/raster", get(raster))
        .route("/queries/{id}/histogram", get(histogram))
        .route("/queries/{id}/tree", get(tree))
        .route("/queries/{id}/variants", get(variants))
        .route("/queries/{id}/sites", get(sites))
        .route("/queries/{id}/export", get(export));
    let mut r = Router::new().nest("/v1", api).with_state(app);
    if let Some(dir) = web {
        r = r.fallback_service(ServeDir::new(dir));
    }
    if cors {
        r = r.layer(CorsLayer::permissive());
    }
    r.layer(CompressionLayer::new()).layer(TraceLayer::new_for_http())
}

pub async fn serve(cfg: Config) -> Result<()> {
    let app = Arc::new(load(&cfg)?);
    let r = router(app, cfg.web.clone(), cfg.cors);
    let listener = tokio::net::TcpListener::bind(cfg.addr).await?;
    info!(addr = %cfg.addr, "listening");
    axum::serve(listener, r)
        .with_graceful_shutdown(async {
            tokio::signal::ctrl_c().await.ok();
        })
        .await?;
    Ok(())
}
