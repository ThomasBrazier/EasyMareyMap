# EasyMareyMap - Estimate recombination maps with Marey maps

![license](https://badgen.net/badge/license/GPL-3.0/blue)
![release](https://badgen.net/badge/release/0.2.0/blue?icon=github)
[![rworkflows](https://github.com/ThomasBrazier/EasyMareyMap/actions/workflows/rworkflows.yml/badge.svg)](https://github.com/ThomasBrazier/EasyMareyMap/actions/workflows/rworkflows.yml)


Estimate local recombination rates with the Marey map method as described in [Brazier and Glémin 2022](https://doi.org/10.1371/journal.pgen.1010141).



This R package provides an easy command line solution for the Marey map method, which is described in more detail in the original MareyMap package. The original MareyMap package was designed for interactive graphical use only and does not provide scripting capabilities. It does not allow to process many maps at a time. Hence we re-implemented the method in a more machine-friendly approach, in order to use it for batch processing and reproducible scripts. See the original MareyMap package here [Rezvoy et al. 2007](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/btm315) and [Siberchicot et al. 2017](https://CRAN.R-project.org/package=MareyMap).


Moreover, it adds new features not included in MareyMap original:
* automatic calibration of the smoothing parameter
* bootstraps to estimate a 95% confidence interval
* graphical zone selection to exclude a bunch of outlier markers
* a graphical and statistical framework for comparative genomics (compare Marey maps between species, populations or genomes)


Two S3 objects are at the core of the package:
* the `marey_map` object contains all the data and metadata associated with a single map. It contains slots with the Marey map, the recombination map (after estimating it), and additional metadata, such as model goodness of fit.
* the `comparative_marey_map` object is a collection of `marey_map` objects over which you can apply `EasyMareyMap` functions or iterate with your custom function. Hence it is a convenient object for manipulating, computing and comparing mutliple maps (i.e. multiple chromosomes and/or datasets).


`EasyMareyMap` functions allow you to:
* estimate the recombination map from the Marey map, using loess or smooth spline regression
* estimate a bootstrapped confidence interval for recombination rates
* predict genetic positions for a new set of markers
* plot one or many Marey maps, plot one or many recombination maps
* compute chromosome wide summary statistics (mean, median, variance, gini index, periphery-bias ratio...)
* comparative plots (lorenz curve, broken stick)
* apply some corrections to a raw Marey map (e.g. outlier removal)


See the vignette for details.


## Installation

Install the package in your R environment with `devtools` (development version).

```
library(devtools)
devtools::install_github("ThomasBrazier/EasyMareyMap")
library(EasyMareyMap)
```

Installation with CRAN and Conda will be available at the first stable release.


## Usage

### Input Data

The input file is a tab-separated text file (.csv, .tsv) with five mandatory columns:
* 'set' dataset name
* 'map' chromosome name
* 'mkr' marker name 
* 'phys' physical position on a reference genome (bp)
* 'gen' genetic position on a linkage map (cM)
* 'vld' only markers TRUE are valid and used in analyses



## Example

Some example datasets are provided to play with the package in the `inst/extdata/` directory.
* a linkage map of *Arabidopsis thaliana* with markers mapped onto a reference genome from [Serin et al. 2017](http://journal.frontiersin.org/article/10.3389/fgene.2017.00201/full)
* female, male and sex-averaged maps of Human, from MareyMapOnline [Rezvoy et al. 2007](https://doi.org/10.1093/bioinformatics/btm315)
* Marey maps of *Brassica napus*, *Brassica rapa*, *Tritucm aestivum* and *Zea mays* from [Brazier and Glémin 2022](https://doi.org/10.1371/journal.pgen.1010141)


![Recombination landscape of *Arabidopsis thaliana* chromosome 1](https://github.com/ThomasBrazier/EasyMareyMap/blob/main/inst/extdata/Arabidopsis_thaliana_chromosome1.jpg?raw=true)
*Recombination landscape of *Arabidopsis thaliana* chromosome 1*


## References

Brazier, Thomas, and Sylvain Glémin. “Diversity and Determinants of Recombination Landscapes in Flowering Plants.” Edited by Ian R. Henderson. PLOS Genetics 18, no. 8 (August 30, 2022): e1010141. https://doi.org/10.1371/journal.pgen.1010141.

Rezvoy, C., Charif, D., Gueguen, L., & Marais, G. A. B. (2007). MareyMap: An R-based tool with graphical interface for estimating recombination rates. Bioinformatics, 23(16), Article 16. https://doi.org/10.1093/bioinformatics/btm315

Serin, Elise A. R., L. B. Snoek, Harm Nijveen, Leo A. J. Willems, Jose M. Jiménez-Gómez, Henk W. M. Hilhorst, and Wilco Ligterink. “Construction of a High-Density Genetic Map from RNA-Seq Data for an Arabidopsis Bay-0 × Shahdara RIL Population.” Frontiers in Genetics 8 (December 5, 2017): 201. https://doi.org/10.3389/fgene.2017.00201.

Siberchicot, A., Rezvoy, C., Charif, D., Gueguen, L., & Marais, G. (2017). MareyMap: Estimation of meiotic recombination rates using marey maps. https://CRAN.R-project.org/package=MareyMap


