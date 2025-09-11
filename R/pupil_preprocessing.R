#' pupil preprocessing
#'
#' wrapper function that takes pupil dilation raw data of both eye and preprocesses into a cleaned pupil size of both eye
#' Includes the following steps that are each avaialble as separate functions:
#' 1. exclude invalid data - exclude_invalid()
#' 2. exclude data around blinks - blink_correciton()
#' 3. exclude dilation speed outlier - exclude_speed_outlier()
#' 4. exclude dilation size outlier - exclude_size_outlier()
#' 5. sparsity filter - exclude islands of data around large gaps of missings
#' 6. interpolate by eye offset - offset_interpolation()
#' 7. interpolate gaps of missing data - missing_interpolation()
#' 8. calculate mean of both eyes
#'
#' @param eye_tracking_raw raw pupil data with estimates for both eyes, timestamps, and optional validity estimates. Expects the following variable names: left_pupil_diameter, right_pupil_diameter, timestamp. Otherwise provide names with provide_variable_names option.
#' @param provide_variable_names provide variable names as string for left_diameter, right_diameter, timestamp as provided by the raw data stream.
#' @param left_diameter_name if provide_variable_names=T, provide name as string
#' @param right_diameter_name if provide_variable_names=T, provide name as string
#' @param timestamp_name if provide_variable_names=T, provide name as string
#' @param diagnostic_output provide diagnostic output that shows excluded data per preprocessing step.
#' @param lower_bound_pupil_diameter smallest pupil size that is considered valid in mm, default = 2.
#' @param upper_bound_pupil_diameter lergest pupil size that is considered valid in mm, default = 8.
#' @param check_raw_validity consider eye tracker validity estimates in excluding invalid data.
#' @param validity_threshold threshold of validity estimates that are considered valid, set between 1 (aggressive) to 4 (very conservative) - specific to TOBII TX-300, default = 2.
#' @param lower_blink_range shortest duration that is considered a blink in milliseconds, default = 75.
#' @param upper_blink_range longest duration that is considered a blink in milliseconds, default = 250.
#' @param cut_blink_data data to cut before and after blinks in ms, default 25ms.
#' @param sampling_rate eye-tracker sampling rate in Hz, default 300Hz.
#' @param MAD_constant Mean Average Deviation constant as threshold of deviation from mean to be used to define outlier. Default 3.
#' @param smooth_length length of linear imputation to be applied before outlier exclusion.
#' @param NA_gap_length gaps >= this length in milliseconds are considered "big" and define split points in the sparsity filtering.
#' @param sparsity_filter_island_length valid sections shorter than this in ms are rejected and thus set to NA.
#' @param interpolation_length gaps shorter than this in millisecond that will be linearly interpolated.
#' @return returns preprocessed pupil size estimate as mean of both eyes.
#' @examples
#' #to be done
#' @export
pupil_preprocessing<-function(
    eye_tracking_raw,
    provide_variable_names=F,
    left_diameter_name=NA,
    right_diameter_name=NA,
    timestamp_name=NA,
    diagnostic_output=T,
    lower_bound_pupil_diameter=2, #in mm
    upper_bound_pupil_diameter=8, #in mm
    check_raw_validity=F,
    validity_threshold=2,
    lower_blink_range=75,
    upper_blink_range=250,
    cut_blink_data=25,
    sampling_rate=300,
    MAD_constant=3,
    smooth_length=150, #measured in ms - linear interpolation before smoothing
    NA_gap_length=75, #NA sequences longer than X to consider in sparsity filter
    sparsity_filter_island_length=50, #in ms - islands of data to remove in sparsity filter
    interpolation_length=300 # measured in ms - NA inteprolation
    ){

  #TODO:
  #diagnostic procedure
  #eye_tracking_raw<-sample_dataframe

  #define variables from raw eye tracking data
  if(provide_variable_names){
    Left_Diameter<-eye_tracking_raw[[left_diameter_name]]
    Right_Diameter<-eye_tracking_raw[[right_diameter_name]]
    RemoteTime<-eye_tracking_raw[[timestamp_name]] #timestamp resolution in milliseconds
  }

  if(!provide_variable_names){
  Left_Diameter<-eye_tracking_raw$left_pupil_diameter
  Right_Diameter<-eye_tracking_raw$right_pupil_diameter
  RemoteTime<-eye_tracking_raw$timestamp #timestamp resolution in milliseconds
  }

  if(check_raw_validity){
    Left_Validity<-eye_tracking_raw$left_pupil_validity
    Right_Validity<-eye_tracking_raw$right_pupil_validity}

  ###CORE PREPROCESSING:

  # STEP 1 - exclude invalid data ####
  pl <- exclude_invalid(pupil_diameter=Left_Diameter,
                        lower_bound_pupil_diameter=lower_bound_pupil_diameter,
                        upper_bound_pupil_diameter=upper_bound_pupil_diameter,
                        check_raw_validity=check_raw_validity,
                        validity_data=Left_Validity,
                        validity_threshold=validity_threshold
                        )

  pr <- exclude_invalid(pupil_diameter=Right_Diameter,
                        lower_bound_pupil_diameter=lower_bound_pupil_diameter,
                        upper_bound_pupil_diameter=upper_bound_pupil_diameter,
                        check_raw_validity=check_raw_validity,
                        validity_data=Right_Validity,
                        validity_threshold=validity_threshold
  )

  # STEP 2) delete data around blinks ####
  #gaps=missing data sections > 75ms and  <=250ms, otherwise not likely to be a blink
  #to be excluded: samples within 50 ms of gaps -> +-25 (8 data points) oder 50?
  pl<-blink_correction(pl,
                       lower_blink_range=lower_blink_range,
                       upper_blink_range=upper_blink_range,
                       cut_blink_data=cut_blink_data,
                       sampling_rate=sampling_rate)

  pr<-blink_correction(pr,
                       lower_blink_range=lower_blink_range,
                       upper_blink_range=upper_blink_range,
                       cut_blink_data=cut_blink_data,
                       sampling_rate=sampling_rate)

  # STEP 3) dilation speed outlier ####
  #take into account time jumps with Remotetimestamps
  #maximum change in pd compared to last and next pd measurement
  pl<-exclude_speed_outlier(pl,timestamp=RemoteTime,MAD_constant=MAD_constant)
  pr<-exclude_speed_outlier(pr,timestamp=RemoteTime,MAD_constant=MAD_constant)

  # STEP 4) dilation size outlier - median absolute deviation ####
  #applies a two pass approach
  #first pass: exclude deviation from trend line derived from all samples
  #second pass: exclude deviation from trend line derived from samples passing first pass
  #--> reintroduction of sample that might have been falsely excluded due to outliers
  #estimate smooth size based on sampling rate
  #take sampling rate into account (300 vs. 120):
  pl<-exclude_size_outlier(pl,timestamp=RemoteTime,
                           MAD_constant=MAD_constant,smooth_length = smooth_length)
  pr<-exclude_size_outlier(pr,timestamp=RemoteTime,
                           MAD_constant=MAD_constant,smooth_length = smooth_length)


  # if(diagnostic_output){
  #   #check for sparsity diagnostic check
  #   pl_before <- pl
  #   pr_before <- pr
  # }

  # STEP 5) sparsity filter ####
  # Remove short islands of valid samples that are bounded by large gaps.
  # signal: numeric vector with NAs where data are invalid
  # time:   numeric timestamp vector (same length as signal), in ms
  # gap_ms: gaps >= this length (ms) are considered "big" and define split points
  # min_segment_ms: valid sections shorter than this (ms) are rejected (set to NA)
  pl<-sparsity_filtering(signal=pl,timestamp=RemoteTime,
                     gap_ms=NA_gap_length,min_segment_ms=sparsity_filter_island_length)
  pr<-sparsity_filtering(signal=pr,timestamp=RemoteTime,
                         gap_ms=NA_gap_length,min_segment_ms=sparsity_filter_island_length)


  # if(diagnostic_output){
  #   #diagnostics on sparsity filter
  #   pct_removed_left  <- mean(is.na(pl)) - mean(is.na(pl_before))
  #   pct_removed_right <- mean(is.na(pr)) - mean(is.na(pr_before))
  #   message(sprintf("Sparsity filter removed ~%.1f%% (L) and ~%.1f%% (R) additional samples.",
  #                   100*pct_removed_left, 100*pct_removed_right))
  # }

  #STEP 6 - offset interpolation ####
  # take offset between left and right into account
  pl<-offset_interpolation(pl,pr,which_eye = 'left')
  pr<-offset_interpolation(pl,pr,which_eye = 'right')

  #STEP 7 - interpolation of NA
  pl<-missing_interpolation(signal=pl,timestamp=RemoteTime,
                            interpolation_length=interpolation_length)
  pr<-missing_interpolation(signal=pr,timestamp=RemoteTime,
                            interpolation_length=interpolation_length)

  #step 8 - mean of both eyes
  pd <- (pl+pr)/2

  eye_tracking_raw['pd']<-pd
  return(eye_tracking_raw)
}
