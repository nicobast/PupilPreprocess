pd_preprocessing<-function(x,nm){

  #print status
  cat("Now processing:", nm, "\n")

  #TODO:
  #mode as "left", "right", "mean" --> with offset correction

  #define thresholds - settings
  diagnostic_output<-FALSE
  check_raw_validity<-TRUE
  lower_bound_pupil_diameter<-2 #in mm
  upper_bound_pupil_diameter<-8 #in mm
  validity_threshold<-2 #set between 1 (aggressive) to 4 (very conservative) - specific to TOBII TX-300
  MAD_constant<-3 ##--> if change speed is higher than MAD_constant * median change --> values are excluded
  smooth_length<-150 #measured in ms - smoothing
  interpolation_length<-300 # measured in ms - NA inteprolation
  NA_gap_length<-75 #NA sequences longer than X to consider in sparsity filter
  sparsity_filter_island_length<-50 #in ms - islands of data to remove in sparsity filter

  #define variables
  Left_Diameter<-x$left_pupil_diameter
  Right_Diameter<-x$right_pupil_diameter
  RemoteTime<-x$timestamp #timestamp resolution in milliseconds
  if(check_raw_validity){
    Left_Validity<-x$left_pupil_validity
    Right_Validity<-x$right_pupil_validity
  }


  # STEP 1 - exclude invalid data ####
  pl <- exclude_invalid(pupil_diameter=Left_Diameter,
                        lower_bound_pupil_diameter=lower_bound_pupil_diameter,
                        upper_bound_pupil_diameter=upper_bound_pupil_diameter,
                        check_raw_validity=check_raw_validity,
                        validity_data=Left_Validity,
                        validity_threshold=validity_threshold
                        )

  pl <- exclude_invalid(pupil_diameter=Right_Diameter,
                        lower_bound_pupil_diameter=lower_bound_pupil_diameter,
                        upper_bound_pupil_diameter=upper_bound_pupil_diameter,
                        check_raw_validity=check_raw_validity,
                        validity_data=Right_Validity,
                        validity_threshold=validity_threshold
  )


  # STEP 2 - filtering ####
  ## STEP 2A) dilation speed outlier ####
  #take into account time jumps with Remotetimestamps
  #maximum change in pd compared to last and next pd measurement

  #Left speed estimation
  pl.speed1<-diff(pl)/diff(RemoteTime) #compared to last
  pl.speed2<-diff(rev(pl))/diff(rev(RemoteTime)) #compared to next
  pl.speed1<-c(NA,pl.speed1)
  pl.speed2<-c(rev(pl.speed2),NA)
  pl.speed<-pmax(pl.speed1,pl.speed2,na.rm=T)
  rm(pl.speed1,pl.speed2)

  #Right speed estimation
  pr.speed1<-diff(pr)/diff(RemoteTime) #compared to last
  pr.speed2<-diff(rev(pr))/diff(rev(RemoteTime)) #compared to next
  pr.speed1<-c(NA,pr.speed1)
  pr.speed2<-c(rev(pr.speed2),NA)
  pr.speed<-pmax(pr.speed1,pr.speed2,na.rm=T)
  rm(pr.speed1,pr.speed2)

  #median absolute deviation of speed
  #MAD_constant<-3
  pl.speed.med<-median(pl.speed,na.rm=T)
  pl.mad<-median(abs(pl.speed-pl.speed.med),na.rm = T)
  pl.treshold.speed<-pl.speed.med+MAD_constant*pl.mad #treshold.speed units are mm/microsecond
  #plot(abs(pl.speed))+abline(h=pl.treshold.speed)
  pr.speed.med<-median(pr.speed,na.rm=T)
  pr.mad<-median(abs(pr.speed-pr.speed.med),na.rm = T)
  pr.treshold.speed<-pr.speed.med+MAD_constant*pr.mad #treshold.speed units are mm/microsecond
  #plot(abs(pr.speed))+abline(h=pr.treshold.speed)
  #correct pupil dilation for speed outliers
  pl<-ifelse(abs(pl.speed)>pl.treshold.speed,NA,pl)
  pr<-ifelse(abs(pr.speed)>pr.treshold.speed,NA,pr)

  ## STEP 2B) delete data around blinks ####
  #gaps=missing data sections > 75ms; Leonie: also <=250ms, otherwise not likely to be a blink
  #to be excluded: samples within 50 ms of gaps -> +-25 (8 data points) oder 50?
  pl<-blink_cor(pl)
  pr<-blink_cor(pr)

  ## STEP 2C) normalized dilation size - median absolute deviation -SIZE ####
  #applies a two pass approach
  #first pass: exclude deviation from trend line derived from all samples
  #second pass: exclude deviation from trend line derived from samples passing first pass
  #--> reintroduction of sample that might have been falsely excluded due to outliers
  #estimate smooth size based on sampling rate
  #take sampling rate into account (300 vs. 120):
  #smooth.size<-round(smooth.length/mean(diff(RemoteTime)/1000)) #timestamp resolution in microseconds
  smooth.size<-round(smooth_length/median(diff(RemoteTime),na.rm=T)) #timestamp resolution in milliseconds
  is.even<-function(x){x%%2==0}
  smooth.size<-ifelse(is.even(smooth.size)==T,smooth.size+1,smooth.size) #make sure to be odd value (see runmed)
  #Left
  pl.smooth<-na.approx(pl,na.rm=F,rule=2) #impute missing values with interpolation
  #pl.smooth<-runmed(pl.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pl.smooth))!=0){pl.smooth<-runmed(pl.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pl.mad<-median(abs(pl-pl.smooth),na.rm=T)
  #Right
  pr.smooth<-na.approx(pr,na.rm=F,rule=2) #impute missing values with interpolation
  #pr.smooth<-runmed(pr.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pr.smooth))!=0){pr.smooth<-runmed(pr.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pr.mad<-median(abs(pr-pr.smooth),na.rm=T)
  #correct pupil dilation for size outliers - FIRST pass
  pl.pass1<-ifelse((pl>pl.smooth+MAD_constant*pl.mad)|(pl<pl.smooth-MAD_constant*pl.mad),NA,pl)
  pr.pass1<-ifelse((pr>pr.smooth+MAD_constant*pr.mad)|(pr<pr.smooth-MAD_constant*pr.mad),NA,pr)
  #Left
  pl.smooth<-na.approx(pl.pass1,na.rm=F,rule=2) #impute missing values with interpolation
  #pl.smooth<-runmed(pl.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pl.smooth))!=0){pl.smooth<-runmed(pl.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pl.mad<-median(abs(pl-pl.smooth),na.rm=T)
  #Right
  pr.smooth<-na.approx(pr.pass1,na.rm=F,rule=2) #impute missing values with interpolation
  #pr.smooth<-runmed(pr.smooth,k=smooth.size) #smooth algorithm by running median of 15 * 3.3ms
  if(sum(!is.na(pr.smooth))!=0){pr.smooth<-runmed(pr.smooth,k=smooth.size)} #run smooth algo only if not all elements == NA
  pr.mad<-median(abs(pr-pr.smooth),na.rm=T)
  #correct pupil dilation for size outliers - SECOND pass
  pl.pass2<-ifelse((pl>pl.smooth+MAD_constant*pl.mad)|(pl<pl.smooth-MAD_constant*pl.mad),NA,pl)
  pr.pass2<-ifelse((pr>pr.smooth+MAD_constant*pr.mad)|(pr<pr.smooth-MAD_constant*pr.mad),NA,pr)
  pl<-pl.pass2
  pr<-pr.pass2

  if(diagnostic_output){
    #check for sparsity diagnostic check
    pl_before <- pl
    pr_before <- pr
  }

  ## D) sparsity filter ####
  # Remove short islands of valid samples that are bounded by large gaps.
  # signal: numeric vector with NAs where data are invalid
  # time:   numeric timestamp vector (same length as signal), in ms
  # gap_ms: gaps >= this length (ms) are considered "big" and define split points
  # min_segment_ms: valid sections shorter than this (ms) are rejected (set to NA)
  sparsity_filter <- function(signal, time, gap_ms = NA_gap_length, min_segment_ms = sparsity_filter_island_length) {
    stopifnot(length(signal) == length(time))
    if (!length(signal)) return(signal)

    # ms per sample (robust to small timing jitter)
    ms_per_sample <- median(diff(time), na.rm = TRUE)
    if (!is.finite(ms_per_sample) || ms_per_sample <= 0) {
      stop("Cannot compute ms_per_sample from `time`.")
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
    signal
  }

  if(diagnostic_output){
    #diagnostics on sparsity filter
    pct_removed_left  <- mean(is.na(pl)) - mean(is.na(pl_before))
    pct_removed_right <- mean(is.na(pr)) - mean(is.na(pr_before))
    message(sprintf("Sparsity filter removed ~%.1f%% (L) and ~%.1f%% (R) additional samples.",
                    100*pct_removed_left, 100*pct_removed_right))
  }

  # STEP 3 - processing valid samples  ####
  #take offset between left and right into account
  pd.offset<-pl-pr
  pd.offset<-na.approx(pd.offset,rule=2)
  #mean pupil dilation across both eyes
  pl <- ifelse(is.na(pl)==FALSE, pl, pr+pd.offset)
  pr <- ifelse(is.na(pr)==FALSE, pr, pl-pd.offset)

  #interpolation of NA (for <=300ms)
  interpolation_gap<-round(interpolation_length/median(diff(RemoteTime),na.rm=T)) #timestamp resolution in milliseconds
  pl<-na.approx(pl, na.rm=F, maxgap=interpolation_gap, rule=2)
  pr<-na.approx(pr, na.rm=F, maxgap=interpolation_gap, rule=2)

  pd <- (pl+pr)/2
  # end of function --> return ####
  #detach(x)

  x[,'pd']<-pd
  return(x)
}
