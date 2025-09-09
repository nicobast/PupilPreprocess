#' sparsity filtering
#'
#' This function removes short islands of valid samples that are bounded by large gaps.
#'
#' @param signal numeric vector with NAs where data are invalid.
#' @param timestamp timestamp in the raw data associated with signal. Expects millisecond format timestamps.
#' @param gap_ms gaps >= this length in ms are considered "big" and define split points.
#' @param min_segment_ms valid sections shorter than this in ms are rejected and thus set to NA
#' @return signal cleared of sparsity islands. These are set to NA.
#' @examples
#' #to be done
#' @export
sparsity_filtering <- function(signal,
                            timestamp,
                            gap_ms = 75,
                            min_segment_ms = 50) {

  #stop if vectors do not allign
  stopifnot(length(signal) == length(timestamp))
  if (!length(signal)) return(signal)

  # ms per sample (robust to small timing jitter)
  ms_per_sample <- median(diff(timestamp), na.rm = TRUE)
  if (!is.finite(ms_per_sample) || ms_per_sample <= 0) {
    stop("Cannot compute ms_per_sample from `timestamp`.")
  }

  gap_thresh   <- max(1L, round(gap_ms / ms_per_sample))
  seg_min_len  <- max(1L, round(min_segment_ms / ms_per_sample))

  is_gap <- is.na(signal)
  rl <- rle(is_gap)
  n_runs <- length(rl$lengths)
  run_ends <- cumsum(rl$lengths)
  run_starts <- run_ends - rl$lengths + 1L

  # Identify "big" gap runs (TRUE-runs with length >= gap_thresh)
  big_gap_idx <- which(rl$values & rl$lengths >= gap_thresh)

  # For each valid run (FALSE), if it is between big gaps (or at edge) and too short -> drop
  to_na <- logical(length(signal))

  for (i in seq_len(n_runs)) {
    if (rl$values[i]) next  # skip gap runs; we only process valid segments
    len_i <- rl$lengths[i]
    if (len_i >= seg_min_len) next

    prev_is_big <- (i > 1L)  && ((i - 1L) %in% big_gap_idx)
    next_is_big <- (i < n_runs) && ((i + 1L) %in% big_gap_idx)

    # Treat edges as natural boundaries, mirroring the original guideline
    left_bound  <- prev_is_big || (i == 1L)
    right_bound <- next_is_big || (i == n_runs)

    if (left_bound && right_bound) {
      rng <- run_starts[i]:run_ends[i]
      to_na[rng] <- TRUE
    }
  }

  signal[to_na] <- NA_real_
  return(signal)
}
