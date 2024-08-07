# EasyMareyMap

Estimate local recombination rates with the Marey map method as described in [Brazier and Glémin 2022](https://doi.org/10.1371/journal.pgen.1010141).

This R package provides an easy command line solution for the Marey map method, which is described in more detail in the original MareyMap package. The original MareyMap package is designed for interactive graphical use only and does not provide scripting capabilities. It does not allow to process many maps at a time. Hence we re-implemented the method in a more machine-friendly approach, in order to use it for batch processing and reproducible scripts. See the original MareyMap package here [Rezvoy et al. 2007](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btm315) and [Siberchicot et al. 2017](https://CRAN.R-project.org/package=MareyMap).


Moreover, it adds new features not included in MareyMap original:
* automatic calibration of the smoothing parameter
* bootstraps to estimate a 95% confidence interval
* graphical zone selection to exclude a bunch of outlier markers
* a statistical framework for comparative genomics (compare marey maps between species, populations or genomes)



## Installation

Install the package in your R environment with `devtools` (development version).

```
library(devtools)
devtools::install_github("ThomasBrazier/EasyMareyMap")
library(EasyMareyMap)
```

## Usage

### Input Data

The input file is a tab-separated text file with five columns:
* 'set' dataset name
* 'map' chromosome name
* 'mkr' marker name 
* 'phys' physical position on a reference genome (bp)
* 'gen' genetic position on a linkage map (cM)
* 'vld' only markers TRUE are valid and used in analyses


### mareyMap objects

A `mareyMap` object must contain only one genetic map (i.e. one chromosome), but the input data can contain multiple chromosome (hence chromosome name must be specified).



## Example

Some example datasets are provided to play with the package in the `data/` directory.
* a linkage map of *Arabidopsis thaliana* with markers mapped onto a reference genome from [Serin et al. 2017](http://journal.frontiersin.org/article/10.3389/fgene.2017.00201/full)
* female, male and sex-averaged maps of Human, from MareyMapOnline [Rezvoy et al. 2007](https://doi.org/10.1093/bioinformatics/btm315)


![Recombination landscape of *Arabidopsis thaliana* chromosome 1](https://github.com/ThomasBrazier/EasyMareyMap/blob/main/inst/extdata/Arabidopsis_thaliana_chromosome1.jpg?raw=true)
*Recombination landscape of *Arabidopsis thaliana* chromosome 1*


## References

Brazier, Thomas, and Sylvain Glémin. “Diversity and Determinants of Recombination Landscapes in Flowering Plants.” Edited by Ian R. Henderson. PLOS Genetics 18, no. 8 (August 30, 2022): e1010141. https://doi.org/10.1371/journal.pgen.1010141.

Rezvoy, C., Charif, D., Gueguen, L., & Marais, G. A. B. (2007). MareyMap : An R-based tool with graphical interface for estimating recombination rates. Bioinformatics, 23(16), Article 16. https://doi.org/10.1093/bioinformatics/btm315


Serin, Elise A. R., L. B. Snoek, Harm Nijveen, Leo A. J. Willems, Jose M. Jiménez-Gómez, Henk W. M. Hilhorst, and Wilco Ligterink. “Construction of a High-Density Genetic Map from RNA-Seq Data for an Arabidopsis Bay-0 × Shahdara RIL Population.” Frontiers in Genetics 8 (December 5, 2017): 201. https://doi.org/10.3389/fgene.2017.00201.

