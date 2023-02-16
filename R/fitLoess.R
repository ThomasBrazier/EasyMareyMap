fitLoess = function(x, span = numeric(), degree = 2) {
  fit = loess(gen ~ phys, x, span = span, degree = degree)
  return(fit)
}
