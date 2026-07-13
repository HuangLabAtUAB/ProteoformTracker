# Builds the genome-wide exon structure index (R/reference_exon_index.R)
# from a bulk Ensembl GTF. Download the GTF first if not already present:
#
#   mkdir -p data/annotation
#   curl -o data/annotation/Homo_sapiens.GRCh38.116.gtf.gz \
#     http://ftp.ensembl.org/pub/release-116/gtf/homo_sapiens/Homo_sapiens.GRCh38.116.gtf.gz
#
# Validated against CD44 (canonical ENST00000428726, 742 aa; isoform 11
# ENST00000434472, 429 aa) before running genome-wide -- exact match against
# independently cross-checked residue boundaries (see commit history).

source("R/reference_exon_index.R")

gtf_path <- "data/annotation/Homo_sapiens.GRCh38.116.gtf.gz"
if (!file.exists(gtf_path)) {
  stop("No GTF found at ", gtf_path, " -- download it first (see header of this script)")
}

message("Building genome-wide exon index (this reads/parses the full GTF; a few minutes)...")
index <- build_reference_exon_index(gtf_path, prefiltered = FALSE)

out_path <- "data/reference_exon_index.rds"
saveRDS(index, out_path)
message(
  nrow(index), " exon rows across ", length(unique(index$transcript_id)), " transcripts, ",
  length(unique(index$gene_name)), " genes -> ", out_path
)
