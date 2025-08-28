#' exclude speed outlier
#'
#' This function estimates speed of pupil size change based on timestamp data
#' and excludes data that is outside defined ranges. Detected speed outlier are set to NA.
#' Consider time jumps in timestamp vector. Speed change is estimate forward and backwards.
#'
#' @param signal signal as a numeric vector, typically of raw pupillometry data.
#' @param timestamp timestamp in the raw data associated with signal. Works independently of the resolution of the timestamp
#' @param  MAD_constant constant as threshold to be used to define outlier. Default 3.
#' @return signal cleared of speed outliers. These are set to NA.
#' @examples
#' #to be done
#' @importFrom stats median rnorm
#' @export
exclude_speed_outlier <- function(
  signal,
  timestamp,
  MAD_constant=3
  )
{
#speed estimation
signal.speed1<-diff(signal)/diff(timestamp) #compared to last
signal.speed2<-diff(rev(signal))/diff(rev(timestamp)) #compared to next
signal.speed1<-c(NA,signal.speed1)
signal.speed2<-c(rev(signal.speed2),NA)
signal.speed<-pmax(signal.speed1,signal.speed2,na.rm=T)

#median absolute deviation of speed
signal.speed.med<-median(signal.speed,na.rm=T)
signal.mad<-median(abs(signal.speed-signal.speed.med),na.rm = T)
signal.treshold.speed<-signal.speed.med+MAD_constant*signal.mad
#treshold.speed units are mm/microsecond

#exclude outlier
signal<-ifelse(abs(signal.speed)>signal.treshold.speed,NA,signal)
return(signal)
}

