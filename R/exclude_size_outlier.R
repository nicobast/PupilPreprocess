#' exclude size outlier
#'
#' This function estimates the median absolute deviation of pupil size from a trend line.
#' Applies a two pass approach: first pass excludes deviation from trend line derived from all samples.
#' second pass excludes deviation from trend line derived from samples passing first pass.
#' reintroduction of samples that might have been falsely excluded due to outliers.
#' Detected size outlier are set to NA.
#'
#' @param signal signal as a numeric vector, typically of raw pupillometry data.
#' @param timestamp timestamp in the raw data associated with signal. Expects millisecond format timestamps
#' @param MAD_constant constant as threshold to be used to define outlier. Default 3.
#' @param smooth_length length of linear imputation to be applied before outlier exclusion
#' @return signal cleared of size outliers. These are set to NA.
#' @examples
#' #to be done
#' @importFrom zoo na.approx
#' @export
exclude_size_outlier <- function(
    signal,
    timestamp,
    MAD_constant=3,
    smooth_length=150){

#estimate smooth size in samples of signal based on length in ms and samsignaling rate in timestamp difference
smooth.size<-round(smooth_length/median(diff(timestamp),na.rm=T)) #timestamp resolution in milliseconds
is.even<-function(x){x%%2==0}
smooth.size<-ifelse(is.even(smooth.size)==T,smooth.size+1,smooth.size) #make sure to be odd value (see runmed)

#impute missing values with interpolation
signal.smooth<-na.approx(signal,na.rm=F,rule=2)
#run smooth algo only if not all elements == NA
if(sum(!is.na(signal.smooth))!=0){signal.smooth<-stats::runmed(signal.smooth,k=smooth.size)}
signal.mad<-median(abs(signal-signal.smooth),na.rm=T)
#FIRST pass: correct pupil dilation for size outliers
signal.pass1<-ifelse((signal>signal.smooth+MAD_constant*signal.mad)|(signal<signal.smooth-MAD_constant*signal.mad),NA,signal)
#impute missing values with interpolation
signal.smooth<-na.approx(signal.pass1,na.rm=F,rule=2)
#run smooth algo only if not all elements == NA
if(sum(!is.na(signal.smooth))!=0){signal.smooth<-stats::runmed(signal.smooth,k=smooth.size)}
signal.mad<-median(abs(signal-signal.smooth),na.rm=T)
#SECOND pass: correct pupil dilation for size outliers
signal.pass2<-ifelse((signal>signal.smooth+MAD_constant*signal.mad)|(signal<signal.smooth-MAD_constant*signal.mad),NA,signal)
signal<-signal.pass2
return(signal)

}
