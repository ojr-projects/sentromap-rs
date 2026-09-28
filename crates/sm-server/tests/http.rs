//! End-to-end API test on a synthetic genome: query lifecycle, views, orientation and errors.

use std::sync::Arc;

use axum::body::Body;
use axum::http::{Request, StatusCode};
use http_body_util::BodyExt;
use serde_json::{Value, json};
use sm_core::kmer;
use sm_core::synth::{self, SynthSpec};
use tower::ServiceExt;

struct Fixture {
    _dir: tempfile::TempDir,
    router: axum::Router,
    query: String,
    oracle_sites: Vec<Vec<sm_core::Hit>>,
}

fn fixture() -> Fixture {
    let dir = tempfile::tempdir().unwrap();
    let s = synth::generate(&SynthSpec { seed: 9, copies_per_family: 15, ..Default::default() });
    let fa = dir.path().join("g.fa");
    s.genome.write_fasta(&mut std::fs::File::create(&fa).unwrap()).unwrap();
    let idx_dir = dir.path().join("idx");
    sm_index::build(&fa, &idx_dir, &sm_index::BuildOptions::default()).unwrap();
    let cfg = sm_server::Config {
        index: idx_dir,
        addr: "127.0.0.1:0".parse().unwrap(),
        web: None,
        cache_bytes: 1 << 30,
        populate: false,
        cors: false,
        calibrate: false,
    };
    let app = Arc::new(sm_server::load(&cfg).unwrap());
    let query = String::from_utf8(s.family_sources[1][3..34].to_ascii_uppercase()).unwrap();
    let q = kmer::parse(&query).unwrap();
    let oracle_sites = (0..=6).map(|n| sm_core::oracle::search_genome(&s.genome, q, n)).collect();
    Fixture { _dir: dir, router: sm_server::router(app, None, false), query, oracle_sites }
}

async fn call(r: &axum::Router, method: &str, uri: &str, body: Option<Value>) -> (StatusCode, Vec<u8>) {
    let mut req = Request::builder().method(method).uri(uri);
    let body = match body {
        Some(b) => {
            req = req.header("content-type", "application/json");
            Body::from(b.to_string())
        }
        None => Body::empty(),
    };
    let resp = r.clone().oneshot(req.body(body).unwrap()).await.unwrap();
    let status = resp.status();
    (status, resp.into_body().collect().await.unwrap().to_bytes().to_vec())
}

async fn json_call(r: &axum::Router, method: &str, uri: &str, body: Option<Value>) -> (StatusCode, Value) {
    let (s, b) = call(r, method, uri, body).await;
    (s, serde_json::from_slice(&b).unwrap_or(Value::Null))
}

async fn wait_ready(r: &axum::Router, id: u64) -> Value {
    for _ in 0..500 {
        let (_, v) = json_call(r, "GET", &format!("/v1/queries/{id}"), None).await;
        if v["status"] == "ready" {
            return v;
        }
        assert_ne!(v["status"], "failed", "{v}");
        tokio::time::sleep(std::time::Duration::from_millis(10)).await;
    }
    panic!("query never became ready");
}

#[tokio::test(flavor = "multi_thread")]
async fn query_lifecycle() {
    let f = fixture();
    let r = &f.router;
    let (s, g) = json_call(r, "GET", "/v1/genome", None).await;
    assert_eq!(s, StatusCode::OK);
    let len = g["genome_len"].as_u64().unwrap();

    let (s, created) = json_call(r, "POST", "/v1/queries", Some(json!({ "kmer": f.query, "max_n": 6 }))).await;
    assert_eq!(s, StatusCode::OK, "{created}");
    let id = created["id"].as_u64().unwrap();
    let status = wait_ready(r, id).await;
    // Counts per n agree with the oracle.
    for n in 0..=6 {
        assert_eq!(status["counts"][n]["sites"].as_u64().unwrap() as usize, f.oracle_sites[n].len(), "n={n}");
    }

    // Sites endpoint returns exactly the oracle's sites (positions and strands).
    let (_, sites) = json_call(r, "GET", &format!("/v1/queries/{id}/sites?n=4&limit=100000"), None).await;
    let mut got: Vec<(u64, String)> = sites["sites"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| (s["pos"].as_u64().unwrap(), s["strand"].as_str().unwrap().to_string()))
        .collect();
    got.sort();
    let mut want: Vec<(u64, String)> =
        f.oracle_sites[4].iter().map(|h| (h.pos as u64, h.strand.symbol().to_string())).collect();
    want.sort();
    assert_eq!(got, want);

    // Raster covers all sites; binary layout is 4 bytes per cell.
    let (s, bytes) =
        call(r, "GET", &format!("/v1/queries/{id}/raster?n=6&width=300&height=40&start=0&end={len}"), None).await;
    assert_eq!(s, StatusCode::OK);
    assert_eq!(bytes.len(), 300 * 40 * 4);
    let (_, hist) = call(r, "GET", &format!("/v1/queries/{id}/histogram?n=6&width=10"), None).await;
    let total: u64 = hist.as_chunks::<4>().0.iter().map(|b| u32::from_le_bytes(*b) as u64).sum();
    assert_eq!(total as usize, f.oracle_sites[6].len());

    for order in ["per_n", "frozen", "mst"] {
        let (s, t) = json_call(r, "GET", &format!("/v1/queries/{id}/tree?n=5&order={order}"), None).await;
        assert_eq!(s, StatusCode::OK, "{order}: {t}");
        let (_, v) =
            json_call(r, "GET", &format!("/v1/queries/{id}/variants?n=5&order={order}&limit=10000"), None).await;
        assert_eq!(v["variants"].as_array().unwrap().len(), t["rows"].as_u64().unwrap() as usize);
    }

    // Export streams every site within n.
    let (s, tsv) = call(r, "GET", &format!("/v1/queries/{id}/export?n=3&format=tsv"), None).await;
    assert_eq!(s, StatusCode::OK);
    assert_eq!(String::from_utf8(tsv).unwrap().lines().count() - 1, f.oracle_sites[3].len());

    // The reverse complement is served from the cache, mirrored.
    let rc = kmer::to_string(kmer::revcomp(kmer::parse(&f.query).unwrap()));
    let (_, c2) = json_call(r, "POST", "/v1/queries", Some(json!({ "kmer": rc, "max_n": 4 }))).await;
    assert_eq!(c2["cached"], true);
    let id2 = c2["id"].as_u64().unwrap();
    wait_ready(r, id2).await;
    let (_, s2) = json_call(r, "GET", &format!("/v1/queries/{id2}/sites?n=4&limit=100000"), None).await;
    let mut got2: Vec<(u64, String)> = s2["sites"]
        .as_array()
        .unwrap()
        .iter()
        .map(|s| (s["pos"].as_u64().unwrap(), s["strand"].as_str().unwrap().to_string()))
        .collect();
    got2.sort();
    let mut flipped: Vec<(u64, String)> =
        want.iter().map(|(p, s)| (*p, if s == "+" { "-".to_string() } else { "+".to_string() })).collect();
    flipped.sort();
    assert_eq!(got2, flipped);

    // Errors.
    let (s, _) = json_call(r, "POST", "/v1/queries", Some(json!({ "kmer": "ACGT" }))).await;
    assert_eq!(s, StatusCode::BAD_REQUEST);
    let (s, _) = json_call(r, "POST", "/v1/queries", Some(json!({ "kmer": f.query, "max_n": 30 }))).await;
    assert_eq!(s, StatusCode::BAD_REQUEST);
    let (s, _) = json_call(r, "GET", "/v1/queries/999999", None).await;
    assert_eq!(s, StatusCode::NOT_FOUND);
    let (s, _) = json_call(r, "DELETE", &format!("/v1/queries/{id}"), None).await;
    assert_eq!(s, StatusCode::OK);
}
