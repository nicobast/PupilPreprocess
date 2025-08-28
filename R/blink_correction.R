#' blink correction
#'
#' This function corrects a numeric vector of pupil data of potential eye blinks.
#'
#' @param signal signal a numeric vector, typically of raw pupillometry data.
#' @param  lower_blink_range lower bound of blinks in ms. Default 25ms.
#' @param  upper_blink_range upper bound of blinks in ms. Default 250ms.
#' @param  cut_blink_data data to cut before and after blinks in ms. Default 25ms
#' @param  sampling_rate eye-tracker sampling rate in Hz. Default 300Hz
#' @return signal cleared of potential blink data. These are set to NA.
#' @examples
#' #to be done
#' @export
blink_correction <- function(signal,
                      lower_blink_range=75,
                      upper_blink_range=250,
                      cut_blink_data=25,
                      sampling_rate=300) {

  #define settings for blink correction algorithm based on sampling rate
  lower_threshold<-round(lower_blink_range/(1000/sampling_rate))
  upper_threshold<-round(upper_blink_range/(1000/sampling_rate))
  cut_samples<-round(cut_blink_data/(1000/sampling_rate))
  samples_before<-cut_samples
  samples_after<-cut_samples


  #change NA to 999 for rle()-function
  findna <- ifelse(is.na(signal),999,signal)

  #find blinks:
  #output of rle(): how many times values (NA) are repeated
  repets <- rle(findna)
  #stretch to length of PD vector for indexing
  repets <- rep(repets[["lengths"]], times=repets[["lengths"]])
  #difference between two timestamps~3.33ms -> 75/3.333=22.5 -> wenn 23 Reihen PD=NA, dann blink gap
  #if more than 150ms (45 rows) of NA, missing data due to blink unlikely
  #dummy coding of variables (1=at least 23 consecutive repetitions, 0=less than 23 repetitions)
  repets <- ifelse(repets>=lower_threshold & repets<=upper_threshold, 1, 0)
  #exclude cases where other values than NA (999) are repeated >=23 times by changing dummy value to 0:
  repets[findna!=999 & repets==1] <- 0
  #gives out where changes from 0 to 1 (no NA at least 23xNA) or 1 to 0 (23x NA to no NA) appear
  changes <- c(diff(repets),0)
  #define start (interval before blink/missing data)
  changes.start<-which(changes==1) #where NA-sequence starts
  #gives out row numbers of NA (blink) and previous 8 frames
  start.seq<-unlist(lapply(changes.start, function(x) {seq(max(x-(samples_before-1),1), x)}))
  repets[start.seq]<-1
  #define end (interval after blink/missing data)
  changes.end<-which(changes==-1)+1 #where NA.sequence ends
  #gives out row numbers of NA (blink) and subsequent 8 frames
  end.seq<-unlist(lapply(changes.end, function(x) {seq(x, min(x+(samples_before-1),length(repets)))}))
  repets[end.seq]<-1
  #replace PD data in blink interval (start to end) with NA
  signal[repets==1]<-NA
  return(signal)
}
