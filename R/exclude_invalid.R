#' exclude invalid data
#'
#' This function exclude invalid pupillometry data. Can also check for validity values as provided by Tobii raw data stream
#'
#' @param pupil_diameter vector of pupil data as numeric.
#' @param  lower_bound_pupil_diameter lower bound of plausible pupil diameter.
#' @param  upper_bound_pupil_diameter upper bound of plausible pupil diameter.
#' @param  check_raw_validity TRUE versus FALSE to check pupil data validity based on external validity information.
#' @param  validity_data vector of length pupil diamter that contains external data validity information as integer.
#' @param  validity_threshold integer as threshold to validate pupil data.
#' @return pupil data corrected for invalid data that is outside plausible range (optional: and above validity threshold).
#' @examples
#' #to be done
#' @export
exclude_invalid <- function(
    pupil_diameter,
    lower_bound_pupil_diameter=2,
    upper_bound_pupil_diameter=8,
    check_raw_validity=FALSE,
    validity_data,
    validity_threshold=2 #set between 1 (aggressive) to 4 (very conservative) - specific to TOBII TX-300
){

  # STEP 1 - exclude invalid data ####
  p <- ifelse((pupil_diameter<lower_bound_pupil_diameter|
                 pupil_diameter>upper_bound_pupil_diameter), NA, pupil_diameter)

  # STEP 1B - check validity of raw data
  if(check_raw_validity){
    p <- ifelse(validity_data<=validity_threshold,p,NA)
  }
  return(p)
}



