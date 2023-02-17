# EasyMareyMap

Estimate local recombination rates with the Marey map method [Brazier and Glémin 2022](https://doi.org/10.1371/journal.pgen.1010141).

This R package provide an easy script solution for the Marey map method, as described in more details in the original MareyMap package. See [Rezvoy et al. 2007](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btm315) and [Siberchicot et al. 2017](https://CRAN.R-project.org/package=MareyMap).

The original MareyMap package is based on a graphical user interface and does not allow to easily script the method. This package is a simple implementation for scripting.

Moreover, it adds new features not included in MareyMap original:
* automatic calibration of the smoothing parameter
* bootstraps to estimate a 95% confidence interval
* graphical zone selection to exclude a bunch of outlier markers


## Installation

Install the package in your R environment with `devtools`.

```
library(devtools)
devtools::install_github("ThomasBrazier/EasyMareyMap")
```


## Usage

The input file is a tab-separated text file with five columns:
* 'set' dataset name
* 'map' chromosome name
* 'mkr' marker name 
* 'phys' physical position on a reference genome (bp)
* 'gen' genetic position on a linkage map (cM)
* 'vld' only markers TRUE are valid and used in analyses


## Example

An example dataset is provided to play with the package. It is a linkage map of *Arabidopsis thaliana* with markers mapped onto a reference genome [Serin et al. 2017](http://journal.frontiersin.org/article/10.3389/fgene.2017.00201/full).



## References

Brazier, Thomas, and Sylvain Glémin. “Diversity and Determinants of Recombination Landscapes in Flowering Plants.” Edited by Ian R. Henderson. PLOS Genetics 18, no. 8 (August 30, 2022): e1010141. https://doi.org/10.1371/journal.pgen.1010141.

Serin, Elise A. R., L. B. Snoek, Harm Nijveen, Leo A. J. Willems, Jose M. Jiménez-Gómez, Henk W. M. Hilhorst, and Wilco Ligterink. “Construction of a High-Density Genetic Map from RNA-Seq Data for an Arabidopsis Bay-0 × Shahdara RIL Population.” Frontiers in Genetics 8 (December 5, 2017): 201. https://doi.org/10.3389/fgene.2017.00201.

