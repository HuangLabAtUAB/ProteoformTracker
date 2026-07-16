# Builds a shared exon-based coordinate axis across an arbitrary set of
# transcripts of the same gene, so the exon track and MS2 ladder can be
# drawn aligned to one another the way the validated CD44 prototype did --
# but generalized to any gene, using real genomic exon coordinates already
# in the precomputed exon index (R/reference_exon_index.R) instead of a
# hand-built, single-gene alignment.
#
# Genomic coverage (any transcript's exon at all) is merged into "coverage
# super-intervals" purely to decide axis positioning/order -- each
# transcript's own exon still keeps its own true residue width and lands at
# its own precise offset within that super-interval, strand-aware. This
# mirrors R/fasta_pipeline.R's build_exon_alignment_preview() (see that
# function's doc comment for the full story): requiring an EXACT genomic
# boundary match to align two transcripts' exons -- the original design --
# works fine for two independently-annotated transcripts (their shared
# exons really do have identical boundaries), but breaks down for a novel
# long-read sequence's own exon table, whose boundaries come from
# minimap2's raw alignment segments and can legitimately be a large
# un-spliced block that only partially overlaps a known transcript's much
# finer real exon (e.g. a partial read that didn't resolve an intron).
# Exact-match-only silently showed ZERO overlap in that case, even after
# the same comparison's Module 2 exon-alignment preview (built with the
# overlap-aware version) correctly showed it overlapping.

#' Build a shared exon axis from the exon subtable for a set of transcripts.
#'
#' @param exon_table data.frame from get_exon_structure(), already filtered
#'   or filterable to the transcripts of interest (must include columns
#'   transcript_id, start, end, strand, residue_start, residue_end)
#' @param transcript_ids character vector of transcript ids to include
#' @return list(
#'   axis_exons = data.frame(axis_start, axis_end, length) in 5'->3' order,
#'     axis_start/axis_end are 1-based residue positions on the shared axis
#'   segments = named list, transcript_id -> data.frame(local_start,
#'     local_end, axis_start, axis_end) mapping that transcript's own
#'     residue ranges onto the shared axis
#' )
build_shared_exon_axis <- function(exon_table, transcript_ids) {
  sub <- exon_table[exon_table$transcript_id %in% transcript_ids, ]
  if (nrow(sub) == 0) return(NULL)
  # strand is a factor (levels "+","-","*"); as.character() before comparing
  # -- identical(strand, "-") is always FALSE for a factor value even though
  # it *prints* as "-" (a real bug found and fixed earlier: the minus-strand
  # reversal silently never ran without this).
  strand <- as.character(sub$strand[1])

  # Coverage super-intervals: merge OVERLAPPING (not just exactly-matching)
  # exon boundaries across all transcripts.
  all_exons <- unique(sub[, c("start", "end")])
  all_exons <- all_exons[order(all_exons$start), ]
  supers <- list()
  cur_start <- all_exons$start[1]; cur_end <- all_exons$end[1]
  if (nrow(all_exons) > 1) {
    for (i in 2:nrow(all_exons)) {
      if (all_exons$start[i] <= cur_end) {
        cur_end <- max(cur_end, all_exons$end[i])
      } else {
        supers[[length(supers) + 1]] <- data.frame(start = cur_start, end = cur_end)
        cur_start <- all_exons$start[i]; cur_end <- all_exons$end[i]
      }
    }
  }
  supers[[length(supers) + 1]] <- data.frame(start = cur_start, end = cur_end)
  super_df <- do.call(rbind, supers)

  # ceiling(), not floor(): a codon can straddle an exon-exon junction (2nt
  # in one exon + 1nt in the next, say), so nt_length/3 isn't always a whole
  # number. Flooring previously underallocated axis width, letting a
  # residue's mapped position spill into the next slot.
  super_df$length <- ceiling((super_df$end - super_df$start + 1) / 3)
  super_df$length[super_df$length < 1] <- 1

  if (strand == "-") super_df <- super_df[rev(seq_len(nrow(super_df))), ]
  super_df$axis_end <- cumsum(super_df$length)
  super_df$axis_start <- super_df$axis_end - super_df$length + 1

  find_super <- function(s, e) which(super_df$start <= e & super_df$end >= s)[1]

  segments <- list()
  for (tid in transcript_ids) {
    tsub <- sub[sub$transcript_id == tid, ]
    if (nrow(tsub) == 0) next
    tsub <- tsub[order(tsub$residue_start), ]
    rows <- lapply(seq_len(nrow(tsub)), function(i) {
      si <- find_super(tsub$start[i], tsub$end[i])
      if (is.na(si)) return(NULL)
      # This exon's own offset from its super-interval's 5' edge, converted
      # from bp to a residue-axis offset; strand-aware (a "-"-strand
      # super-interval's 5' edge is its HIGH genomic end).
      offset_bp <- if (strand == "-") super_df$end[si] - tsub$end[i] else tsub$start[i] - super_df$start[si]
      axis_start_i <- super_df$axis_start[si] + round(offset_bp / 3)
      axis_start_i <- min(max(axis_start_i, super_df$axis_start[si]), super_df$axis_end[si])
      # Width comes from this exon's OWN true residue count (already
      # correct per-transcript, from its own reading frame), not re-derived
      # from bp/3 here -- avoids compounding rounding error.
      n_res <- tsub$residue_end[i] - tsub$residue_start[i]
      data.frame(
        local_start = tsub$residue_start[i], local_end = tsub$residue_end[i],
        axis_start = axis_start_i, axis_end = axis_start_i + n_res
      )
    })
    rows <- rows[!vapply(rows, is.null, logical(1))]
    if (length(rows) > 0) segments[[tid]] <- do.call(rbind, rows)
  }

  list(
    axis_exons = super_df[, c("axis_start", "axis_end", "length")],
    axis_length = sum(super_df$length),
    segments = segments
  )
}

#' Map one transcript's own local residue position onto the shared axis
#' built by build_shared_exon_axis(), using its segments table.
#'
#' @param segments one transcript's segment data.frame (from the $segments
#'   list returned by build_shared_exon_axis())
#' @param local_pos 1-based residue position in that transcript's own sequence
#' @return integer axis position, or NA if local_pos falls in an exon that
#'   didn't have a genomic-coordinate match on the shared axis (rare; only
#'   for exon-boundary edge residues at an alternative splice site)
map_local_to_axis <- function(segments, local_pos) {
  for (i in seq_len(nrow(segments))) {
    if (local_pos >= segments$local_start[i] && local_pos <= segments$local_end[i]) {
      return(segments$axis_start[i] + (local_pos - segments$local_start[i]))
    }
  }
  NA_integer_
}
