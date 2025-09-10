#' missing interpolation
#'
#' missing data as NA is linearily interpolated for a specific timeframe
#'
#' @param signal numeric vector of pupil size data.
#' @param  timestamp timestamp in the raw data associated with signal. Expects millisecond format timestamps.
#' @param  interpolation_length gaps of missing data to interpolate in milliseconds
#' @return pupil data with missing interpolated.
#' @examples
#' #to be done
#' @export
missing_interpolation<-function(
    signal,
    timestamp,
    interpolation_length=300
    ){

  #estimate number of samples to be inteprolated based on timestamp resolution in milliseconds
  interpolation_gap<-round(interpolation_length/median(diff(timestamp),na.rm=T))

  #aplly linear interpolation
  signal<-zoo::na.approx(signal, na.rm=F, maxgap=interpolation_gap, rule=2)
  return(signal)

}
