# PupilPreprocess

R package for automated pupil data presprocessing. Provide a wrapper function pupil_preprocessing() that does all preprocessing steps in one go. Also provide single functions for each preprocessing step, which entail:

 1. exclude invalid data - exclude_invalid()
 2. exclude data around blinks - blink_correciton()
 3. exclude dilation speed outlier - exclude_speed_outlier()
 4. exclude dilation size outlier - exclude_size_outlier()
 5. sparsity filter - exclude islands of data around large gaps missings
 6. interpolate by eye offset - offset_interpolation()
 7. interpolate gaps of missing data - missing_interpolation()
 8. calculate mean of both eyes

Package also comes with a generate_pupil_data() function that generates random pupil response data including missings and measurement noise.

See function documentation in R for further insights on setable parameters per function.
---

## Installation

You can install the development version from GitHub with:

```r
# Install remotes if you don't have it
install.packages("remotes")

# Install the package
remotes::install_github("nicobast/PupilPreprocess")
