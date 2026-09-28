# PhenoRec

**PhenoRec** is an R package for **phenotype-preserving integration of single-cell data**.

Current single-cell studies commonly profile samples collected from different biological phenotypes or conditions, such as healthy and diseased individuals in case-control studies, patients at different disease stages or before and after treatment in longitudinal studies.  However, phenotype-associated variation can be inadvertently removed during routine single-cell data integration because it is inherently non-identifiable from batch variation. Consequently, while aligning cells across batches, conventional integration also aligns cells across phenotypes, thereby removing biologically meaningful differences between phenotype groups. 

PhenoRec addresses this problem by explicitly modeling phenotype-associated variation within a generalized latent variable model and introducing an identification strategy to disentangle it from batch variation. PhenoRec integrates single-cell data while preserving biologically meaningful phenotype signals. We expect that PhenoRec will provide a useful computational framework for phenotype-oriented single-cell studies.

![Overview of PhenoRec](PhenoRec.png)

## Installation

PhenoRec can be installed using:

```r
devtools::install_github('chaodeng-aca/PhenoRec')
```

Then load the package:

```r
library(PhenoRec)
```

## Tutorials

For a step-by-step tutorial, please refer to the `vignettes` directory in the repository, which provides usage examples.


## Citation

If you use PhenoRec in your research, please cite:

> Deng C, Zhao H, Wang T. **PhenoRec recovers phenotype signals lost during single-cell data integration.**

A complete citation and publication link will be added after publication of the manuscript.


## License

Please see the `LICENSE` file for licensing information.

