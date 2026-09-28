#!/usr/bin/env bash
# Fetch reference genomes and annotation tracks into the gitignored data/ directory.
# Idempotent: files already present are skipped. Usage: scripts/fetch-data.sh [genomes|dm6-tracks|all]
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="$ROOT/data"
mkdir -p "$DATA"

fetch() { # url dest
    local url=$1 dest=$2
    if [[ -s $dest ]]; then echo "have  $(basename "$dest")"; return; fi
    echo "fetch $(basename "$dest")"
    mkdir -p "$(dirname "$dest")"
    curl -fsSL --retry 3 -o "$dest.part" "$url"
    mv "$dest.part" "$dest"
}

fetch_gunzip() { # url dest (dest is the decompressed file)
    local url=$1 dest=$2
    if [[ -s $dest ]]; then echo "have  $(basename "$dest")"; return; fi
    echo "fetch $(basename "$dest")"
    curl -fsSL --retry 3 "$url" | gunzip -c > "$dest.part"
    mv "$dest.part" "$dest"
}

NCBI=https://ftp.ncbi.nlm.nih.gov/genomes/all

genomes() {
    local ec=$NCBI/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic
    local sc=$NCBI/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic
    local dm=$NCBI/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic
    fetch_gunzip "$ec.fna.gz" "$DATA/ecoli_k12.fna"
    fetch_gunzip "$ec.gff.gz" "$DATA/ecoli_k12.gff"
    fetch_gunzip "$sc.fna.gz" "$DATA/S_cere_R64.fna"
    fetch_gunzip "$sc.gff.gz" "$DATA/S_cere_R64.gff"
    fetch_gunzip "$dm.fna.gz" "$DATA/d_mel.fna"
    fetch_gunzip "$dm.gff.gz" "$DATA/d_mel.gff"
}

# D. melanogaster dm6 annotation tracks (GC% is computed from the genome by `sentromap annotate`). UCSC and ChIP-Atlas use UCSC chromosome names
# (chr2L ...); dm6.chromAlias.txt maps them to the RefSeq accessions in d_mel.fna.
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/dm6
CHIPATLAS=https://chip-atlas.dbcls.jp/data/dm6/eachData/bw

# name  kind  source-id  description   (ChIP-Atlas: uniformly reprocessed, embryo unless noted)
DM6_CHIPATLAS=(
    "H3K4me3      histone  SRX7175247  active promoters (embryo, GSM4173909)"
    "H3K4me2      histone  SRX481094   promoters/enhancers (14-16h embryo, GSM1339123)"
    "H3K4me1      histone  SRX7175245  enhancers (embryo, GSM4173907)"
    "H3K27ac      histone  SRX7175237  active enhancers/promoters (embryo, GSM4173899)"
    "H3K9ac       histone  SRX7175250  active promoters (embryo, GSM4173912)"
    "H3K36me3     histone  SRX9518343  transcribed gene bodies (embryo st5, GSM4910644)"
    "H4K16ac      histone  SRX6386747  X dosage compensation (male NC12 embryo, GSM3913901)"
    "H4K20me1     histone  SRX287868   gene bodies/X (14-16h embryo, GSM1147376)"
    "H3K27me3     histone  SRX7175239  Polycomb repression (embryo, GSM4173901)"
    "H3K9me2      histone  SRX2548341  heterochromatin (4-8h embryo, GSM2481882)"
    "H3K9me3      histone  SRX7175251  constitutive heterochromatin (embryo, GSM4173913)"
    "HP1a         protein  SRX7175260  Su(var)205, heterochromatin protein (embryo, GSM4173922)"
    "Su_var_3-9   protein  SRX033321   H3K9 methyltransferase (0-12h embryo, GSM636838)"
    "CTCF         protein  SRX11499366 insulator (0-14h embryo, GSM5461683)"
    "PolII        protein  SRX7175265  RNA polymerase II (embryo, GSM4173927)"
    "ATAC         access   SRX13305250 chromatin accessibility (6-10h embryo, GSM5714991)"
    "Input        control  SRX7175240  ChIP input; tracks mappability/copy number (embryo, GSM4173902)"
)

dm6_tracks() {
    local dir=$DATA/tracks/dm6
    mkdir -p "$dir"
    fetch "$UCSC/bigZips/dm6.chromAlias.txt" "$dir/dm6.chromAlias.txt"
    fetch "$UCSC/bigZips/dm6.chrom.sizes" "$dir/dm6.chrom.sizes"
    fetch "$UCSC/bigZips/dm6.fa.out.gz" "$dir/rmsk.fa.out.gz"          # RepeatMasker
    fetch "$UCSC/bigZips/dm6.trf.bed.gz" "$dir/trf.bed.gz"              # tandem repeats
    fetch "$UCSC/phyloP124way/dm6.phyloP124way.bw" "$dir/phyloP124way.bw"
    fetch "$UCSC/phastCons124way/dm6.phastCons124way.bw" "$dir/phastCons124way.bw"

    local manifest=$dir/tracks.tsv
    printf 'name\tkind\tformat\tfile\tsource\tdescription\n' > "$manifest.part"
    printf 'phyloP124way\tconservation\tbigwig\tphyloP124way.bw\tUCSC dm6\tper-base conservation, 124 insects\n' >> "$manifest.part"
    printf 'phastCons124way\tconservation\tbigwig\tphastCons124way.bw\tUCSC dm6\tconserved-element probability, 124 insects\n' >> "$manifest.part"
    printf 'RepeatMasker\trepeats\trmsk-out\trmsk.fa.out.gz\tUCSC dm6\tinterspersed and simple repeats\n' >> "$manifest.part"
    printf 'TRF\trepeats\tbed\ttrf.bed.gz\tUCSC dm6\tTandem Repeats Finder\n' >> "$manifest.part"
    local line name kind id desc
    for line in "${DM6_CHIPATLAS[@]}"; do
        read -r name kind id desc <<< "$line"
        fetch "$CHIPATLAS/$id.bw" "$dir/$name.$id.bw"
        printf '%s\t%s\tbigwig\t%s\tChIP-Atlas %s\t%s\n' "$name" "$kind" "$name.$id.bw" "$id" "$desc" >> "$manifest.part"
    done
    mv "$manifest.part" "$manifest"
}

case "${1:-all}" in
    genomes) genomes ;;
    dm6-tracks) dm6_tracks ;;
    all) genomes; dm6_tracks ;;
    *) echo "usage: $0 [genomes|dm6-tracks|all]" >&2; exit 1 ;;
esac
