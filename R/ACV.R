# Absolute Cross-validation
# Lee & Cox 2010
ACV = function(data = data.frame()) {
  criterion = (1/nrow(data))*sum(abs(data$obs - data$pred), na.rm = TRUE)
  return(criterion)
}
