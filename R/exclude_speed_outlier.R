#' exclude speed outlier
#'
#' This function estimates speed of pupil size change based on timestamp data
#' and excludes data that is outside defined ranges. Detected speed outlier are set to NA.
#'
#' @param signal signal a numeric vector, typically of raw pupillometry data.
#' @param  lower_blink_range lower bound of blinks in ms. Default 25ms.
#' @param  upper_blink_range upper bound of blinks in ms. Default 250ms.
#' @param  cut_blink_data data to cut before and after blinks in ms. Default 25ms
#' @param  sampling_rate eye-tracker sampling rate in Hz. Default 300Hz
#' @return signal cleared of potential blink data. These are set to NA.
#' @examples
#' to be done
#' @export
exclude_speed_outlier <- function()
{}


#take into account time jumps with Remotetimestamps
#maximum change in pd compared to last and next pd measurement

#Left speed estimation
pl.speed1<-diff(pl)/diff(RemoteTime) #compared to last
pl.speed2<-diff(rev(pl))/diff(rev(RemoteTime)) #compared to next
pl.speed1<-c(NA,pl.speed1)
pl.speed2<-c(rev(pl.speed2),NA)
pl.speed<-pmax(pl.speed1,pl.speed2,na.rm=T)
rm(pl.speed1,pl.speed2)
