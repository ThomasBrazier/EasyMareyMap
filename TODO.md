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
[ ] Flip a segment


[x] A map object is only a single chromosome, clean and simple. The unit is the genetic map

[ ] Manipulate mareyMap objects
      [ ] A method to reduce a mareyMap object to a given interval
      [ ] A method to change x-axis in relative distances/positions
      [ ] A map object support multiple X axes mapping, and can compare between these mappings
      [ ] Different map objects can be compared


[ ] Quality control and map validation
      [ ] Save model outputs: errors...
      [ ] mesure de noisyness, erreurs de la fonction d'interpolation

 
[ ] Block-wise cM/Mb
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


[ ] Comparison in polyploids: Bourke et al. 2015 https://doi.org/10.1534/genetics.115.181008
¨The slopes of the regression lines for each homolog arm were tested for equality in an analysis of covariance by introducing, where necessary, up to three dummy variables (to code for the presence or absence of a homolog) per chromosome arm per parent (Andrade and Estévez-Pérez 2014). ¨



## Diffusion

* Conda recipe
* Pip recipe
* CRAN
* Biocontainer Docker
