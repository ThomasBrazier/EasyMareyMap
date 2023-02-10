# Least squares Cross validation
# Lee & Cox 2010
LSCV = function(data = data.frame()) {
  criterion = (1/nrow(data))*sum((data$obs - data$pred)^2, na.rm = TRUE)
  return(criterion)
}
