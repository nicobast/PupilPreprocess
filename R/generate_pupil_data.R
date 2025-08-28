#' generate_pupil_data
#'
#' generates a vector in the shape of pupil data.
#'
#' @param sampling_rate sampling rate of the sampled data.
#' @param  phase_length_range range of pupil response durations to sample.
#' @param  number_of_responses number of pupil responses to sample.
#' @param  mean_response_magnitude mean pupil response magnitude in mm
#' @param  introduce_NA simulate missing data
#' @param  NA_proportion proportion of data to set to missing
#' @param  start_value_pupil_size pupil size start value in mm.
#' @param  drift drift in pupil signal per response in mm
#' @return sampled pupil data.
#' @examples
#' to be done
#' @export
generate_pupil_data<-function(
  sampling_rate=300,
  phase_length_range=c(round(sampling_rate/10),round(sampling_rate/5)),
  number_of_responses=100,
  mean_response_magnitude=0.1,
  introduce_NA=T,
  NA_proportion=0.2,
  start_value_pupil_size=4,
  drift=0.01
  ){

  #translate options to required values
  phase_length<-sample(phase_length_range[1]:phase_length_range[2],number_of_responses,replace=T)
  change<-rnorm(number_of_responses+1)*mean_response_magnitude
  difference<-diff(change)
  start_values<-start_value_pupil_size+difference+drift
  end_values<-c(start_values[2:length(start_values)],start_value_pupil_size)

  # Easing function for shaping the ramp
  ease_fn <- function(t, mode = "linear", power = 2L) {
    stopifnot(all(t >= 0 & t <= 1))
    mode <- match.arg(mode, c("linear","ease-in","ease-out","ease-in-out","cosine","sine"))
    switch(
      mode,
      "linear"      = t,
      "ease-in"     = t^power,
      "ease-out"    = 1 - (1 - t)^power,
      "ease-in-out" = ifelse(t < 0.5,
                             0.5 * (2 * t)^power,
                             1 - 0.5 * (2 * (1 - t))^power),
      "cosine"      = (1 - cos(pi * t)) / 2,      # smooth start & end
      "sine"        = sin(pi * t / 2)             # slightly quicker start
    )
  }
  # wave  builder: rising from `from` to `to`, ending in to
  make_wave <- function(from = 3, to = 8, end = to,
                        n_up = 30, n_down = 30,
                        ease_up = "cosine", ease_down = "cosine",
                        power_up = 2, power_down = 2,
                        hold_peak = 0L) {
    stopifnot(n_up >= 2, n_down >= 2, hold_peak >= 0)
    # Parameterize 0..1 for each phase
    t_up   <- seq(0, 1, length.out = n_up)
    t_down <- seq(0, 1, length.out = n_down)

    # Apply easing
    u_up   <- ease_fn(t_up,   ease_up,   power_up)
    u_down <- ease_fn(t_down, ease_down, power_down)

    # Map to values
    up   <- from + (to - from) * u_up
    down <- to   + (end - to) * u_down

    # Build connected sequence: drop the first element of 'down' to avoid duplicating the peak
    c(up, rep(to, hold_peak), down[-1])
  }

  #generate pupil sequence
  test_pupil_data<-unlist(mapply(function(x,y,z){
    make_wave(from = x, to = y, n_up = z, n_down = z)},
    x=start_values,
    y=end_values,
    z=phase_length))

  plot(test_pupil_data)

  #introduce NA
  if(introduce_NA){
  NA_samples<-NA_proportion*length(test_pupil_data)
  min_NA_length<-1
  max_NA_length<-100
  NA_length<-sample(min_NA_length:max_NA_length,round(NA_samples/(max_NA_length-min_NA_length/2)))
  NA_starts<-sample(1:length(test_pupil_data),length(NA_length))
  NA_ends<-NA_starts+NA_length
  NA_ends<-ifelse(NA_ends>length(test_pupil_data),length(test_pupil_data),NA_ends)
  NA_sequences<-unlist(mapply(seq,NA_starts,NA_ends))
  test_pupil_data[NA_sequences]<-NA
  }

  return(test_pupil_data)
}
