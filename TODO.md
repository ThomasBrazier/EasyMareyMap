# TODO

## How to test and compile


### During dev

Load all functions with devtools::load_all(), maybe via Ctrl/Cmd + Shift + L.


### Check and Build

```
library(devtools)
devtools::check()
devtools::build()
```

### Tests

Assuming you’re in a package directory, just run usethis::use_test("name") to create a test file, and set up all the other infrastructure you need. If you’re using RStudio, press Cmd/Ctrl + Shift + T (or run devtools::test() if not) to run all the tests in a package.




## TODO List

[ ] Cleaning steps
See ¨The quality of the genetic map was evaluated visually by generating diagnostic plots in polymapR¨
[x] Flip interval/selection/chromosome (flip genetic distances, not genomic, so don't change the mapping)
[x] Mask interval/selection/chromosome
[x] How to treat negative values (set to zero, don't change, NA)


[x] verbose recombinationMap interpolation


[x] plot_marey -> plot.marey et plot.marey.R
[x] summary.marey.R
[x] Add statistics and goodness of fit to summary


[ ] Commenter jeux de données


[x] A map object is only a single chromosome, clean and simple. The unit is the genetic map
[x] Add a slot distances relatives, pour la carte de Marey, et pour la fonction interpolee. NO, keep it simple


[ ] Interpolation
      [ ] Penalized spline


[ ] Quality control and map validation
      [x] Save model outputs: errors...
      [x] Difference entre MareyCI fitted-bootstrapped et true y
      [x] "Ordinary least squares regression indicates that the estimated positions fit well for each chrome (R2 > 0.99) relative to the true genetic map positions" MapFuser
      [x] fitted genetic dist ~ true/data genetic distance
      [x] mesure de noisyness, erreurs de la fonction d'interpolation, dstance fitter/observed

 
[ ] Block-wise cM/Mb
      [ ] Estimate recombination rates in a set of intervals (estimate on the whole Marey map but output a new rec map on these intervals)
      [ ] From Marey map directly + bootstrap
      [ ] From interpolation + bootstrap




[x] Method to get chromosome-wide statistics
      [x] Gini index
      [x] AUC curves and AUC/Lorenz curve
      [x] coefficient of variation
      [x] periphery-bias ratio
      [x] Veller r-bar
[ ] Graphical method to compare different maps
      [ ] Compare statistics from different blocks/intervals + pairwise blocks (e.g. synteny, homeologs)
      [ ] Lorenz curve : prendre une liste de cartes en argument, renvoyer figure ggplot ou data frame
      [ ] Relative maps : representer les cartes de Marey en distances relatives, prendre une liste de cartes en argument, renvoyer figure ggplot ou data frame
      [ ] Exprimer la difference entre deux cartes (comme avec distances en genes) + mesure statistique de la difference chromosome-wide
      [ ] Difference (ratio) entre deux cartes, le long du genomes, est permise par le passage en distances relatives
      


[ ] Comparison in polyploids: Bourke et al. 2015 https://doi.org/10.1534/genetics.115.181008
¨The slopes of the regression lines for each homolog arm were tested for equality in an analysis of covariance by introducing, where necessary, up to three dummy variables (to code for the presence or absence of a homolog) per chromosome arm per parent (Andrade and Estévez-Pérez 2014). ¨


## Tests

[ ]


## Documentation

[ ]


## Diffusion

[ ] Conda recipe
[ ] Pip recipe
[ ] CRAN
[ ] Biocontainer Docker
