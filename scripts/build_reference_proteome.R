# Downloads the human reviewed canonical reference proteome (UniProt
# proteome UP000005640, Swiss-Prot reviewed entries only -- one
# representative sequence per gene) and builds the offline mass index used
# for confounding-protein lookups (R/reference_proteome_index.R).
#
# The UniProtKB REST "stream" endpoint intermittently fails for a result set
# this size ("Error encountered when streaming data"), so this paginates
# through the "search" endpoint instead, following the Link: rel="next"
# cursor header until exhausted.

library(httr)

fasta_path <- "data/UP000005640_9606_reviewed_canonical.fasta"
index_path <- "data/reference_mass_index.rds"

download_reference_proteome_fasta <- function(out_path, page_size = 500) {
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/search",
    "?query=", utils::URLencode("proteome:UP000005640 AND reviewed:true", reserved = TRUE),
    "&format=fasta&size=", page_size
  )

  con <- file(out_path, open = "w")
  on.exit(close(con), add = TRUE)

  page <- 0
  repeat {
    page <- page + 1
    resp <- httr::GET(url)
    httr::stop_for_status(resp)
    writeLines(httr::content(resp, as = "text", encoding = "UTF-8"), con, sep = "")

    message("fetched page ", page)
    link_header <- httr::headers(resp)$link
    next_url <- if (!is.null(link_header)) {
      m <- regmatches(link_header, regexpr('<[^>]*>(?=; rel="next")', link_header, perl = TRUE))
      if (length(m) > 0) sub("^<(.*)>$", "\\1", m) else NA_character_
    } else {
      NA_character_
    }
    if (is.na(next_url)) break
    url <- next_url
  }
  invisible(out_path)
}

if (!file.exists(fasta_path)) {
  message("Downloading human reviewed canonical reference proteome...")
  download_reference_proteome_fasta(fasta_path)
} else {
  message(fasta_path, " already exists, skipping download")
}

source("R/proteoform_schema.R")
source("R/mass_calculation.R")
source("R/reference_proteome_index.R")

message("Building reference mass index...")
index <- build_reference_mass_index(fasta_path, index_path)
message("Reference mass index built: ", nrow(index), " proteins -> ", index_path)
