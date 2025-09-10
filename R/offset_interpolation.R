#' offset interpolation
#'
#' when pupil size data is provided for both eye. The function estimates the
#' size offset between both eyes and interpolates missing data in one eye by
#' size estimates of the other eye plus offset
#'
#' @param left_pupil_data pupil size data of left eye.
#' @param  right_pupil_data pupil size data of right eye.
#' @param  which_eye data of which eye to interpolate by offset: "left" or "right"
#' @return pupil data with missing interpolated for offset.
#' @examples
#' #to be done
#' @export
offset_interpolation<-function(left_pupil_data,
                               right_pupil_data,
                               which_eye
                               ){

  #estimate offset between eyes
  pd.offset<-left_pupil_data-right_pupil_data

  #interpolate offset estimation
  pd.offset<-zoo::na.approx(pd.offset,rule=2)

  #replace missing data in one eye by estimates of other eye + offset
  if(which_eye=='left'){
    signal <- ifelse(is.na(left_pupil_data)==F, left_pupil_data, right_pupil_data+pd.offset)
  }

  if(which_eye=='right'){
    signal <- ifelse(is.na(right_pupil_data)==F, right_pupil_data, left_pupil_data-pd.offset)
  }

  return(signal)

}

