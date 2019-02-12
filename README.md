# Rabbits
This repository contains Python and R code for processing data, running analyses and plotting figures from 
the paper [Rabbits and the Specious Origins of Domestication](https://doi.org/10.1016/j.tree.2017.12.009)

![∂a∂i model](./dadi-model.png?raw=true)

If you reuse any of this code then please cite the paper:
> Irving-Pease, E.K., Frantz, L.A.F., Sykes, N., Callou, C., Larson, G., 2018. 
> Rabbits and the Specious Origins of Domestication. *Trends in Ecology & Evolution* 33, 149–152.
> https://doi.org/10.1016/j.tree.2017.12.009

## Installation

To reproduce the analyses from the paper you will need to install the following dependencies.

### Python

Python ≥ 2.7 with the following modules:

* [dadi](https://bitbucket.org/gutenkunstlab/dadi)
* [luigi](https://github.com/spotify/luigi)
* [matplotlib](https://github.com/matplotlib/matplotlib)
* [numpy](https://github.com/numpy/numpy)

```bash
pip install luigi matplotlib numpy 
```

```bash
pip install git+https://bitbucket.org/gutenkunstlab/dadi.git
```

The full list of Python modules installed in the project environment can be
found in the `requirement.txt` file.

### R

R ≥ 3.4 with the following modules:

* [ape](https://cran.r-project.org/web/packages/ape/)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/)
* [mapdata](https://cran.r-project.org/web/packages/mapdata/)
* [maps](https://cran.r-project.org/web/packages/maps/)
* [reshape2](https://cran.r-project.org/web/packages/reshape2/)
* [scales](https://cran.r-project.org/web/packages/scales/)
* [stringr](https://cran.r-project.org/web/packages/stringr/)

```R
install.packages(c("ape", "ggplot2", "mapdata", "maps", "reshape2", "scales", "splancs", "stringr"))
```

### Other

* [SRA Toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/)
* [Samtools](https://github.com/samtools/samtools)
* [Picard](http://broadinstitute.github.io/picard/)
* [BWA](http://bio-bwa.sourceforge.net/)
* [GATK](https://software.broadinstitute.org/gatk/)
* [Plink](http://www.cog-genomics.org/plink2)
* [ADMIXTURE](http://software.genetics.ucla.edu/admixture/)
* [sNMF](http://membres-timc.imag.fr/Olivier.Francois/snmf/)
* [FlashPCA](https://github.com/gabraham/flashpca)

## Running the pipeline

The pipeline is broken into two [luigi](https://github.com/spotify/luigi) modules.

### Pre-processing

* download the data from the SRA and Ensembl
* align the data with BWA
* sort and merge with Samtools
* deduplicate with Picard
* variant call and filter with GATK
* convert to Plink
* run ADMIXTURE and sNMF
* run FlashPCA
* plot NJ tree 

```bash
luigi --module pipeline_gatk CustomGenomePipeline
```

### Fitting ∂a∂i models 
 
* compute the SFS
* fit ∂a∂i models


```bash
luigi --module pipeline_dadi CustomDadiPipeline
```

## Author

Evan K. Irving-Pease, [PalaeoBARN](https://www.palaeobarn.com/), University of Oxford 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details
