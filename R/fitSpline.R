
fitSpline = function(x, spar = numeric()) {
  fit = smooth.spline(x = x$phys, y = x$gen, spar = spar)
  return(fit)
}
