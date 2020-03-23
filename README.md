# BiNetEIA

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3724651.svg)](https://doi.org/10.5281/zenodo.3724651)

This is a repository which contains the R-Package *BiNetEIA*. The package has not been uploaded yet to the CRAN repository. In order to use it, you can clone the repository and build it in RStudio.

## How to cite

In case you use de in other software or scientific publications,
please reference this module. It is published and has a DOI. It can be cited
as:

  Schwemmle, R., and Perona, P.: Highlighting action and environmental component interactions using a
  network theory approach, Impact Assessment and Project Appraisal, 1-16,
  [10.1080/14615517.2019.1704113](https://doi.org/10.1080/14615517.2019.1704113), 2020.

## License
This software can be distributed freely under the GPL v3 license. Please read the LICENSE for further information.

© 2019, Robin Schwemmle (<robin.schwemmle@hydrology.uni-freiburg.de>)

## Examples

```r
library(BINetEIA)
```

Here a brief overview is provided how to apply the package. This is exemplified for three real-world cases. The provided code snippets allow to reproduce the figures

### Tallahala
```r
data('tallahala')
vcat <- c("Water", "Land", "Biology", "Socioeconomy", "Infrastructure")
vcat.env.comps <- c(1,1,1,1,2,2,1,1,2,2,3,3,3,3,3,3,3,4,4,4,4)
vcat.actions <- c(1,1,1,5,5,5,1,3,3,5)
visstats(tallahala, negative.values = TRUE, vcat = vcat, vcat.env.comps = vcat.env.comps,
vcat.actions = vcat.actions)
```

### East Karun
```r
data('eastkarun')
vcat <- c("Water", "Land", "Biology", "Socioeconomy")
vcat.env.comps <- c(1,2,2,2,2,1,1,3,3,3,3,3,3,4,4,4,4,4,4,4)
vcat.actions <- c(1,1,1,3,1)
visstats(eastkarun, negative.values = TRUE, vcat = vcat, vcat.env.comps = vcat.env.comps,
vcat.actions = vcat.actions)
```


### Kladovo

```r
data('kladovo')
vcat <- c("Water", "Land", "Biology", "Socioeconomy", "Infrastructure")
vcat.env.comps <- c(1,2,2,2,2,2,3,3,3,3,3,2,2,4,4,4)
vcat.actions <- c(5,5,5,5,5,5,5,5,5)
visstats(kladovo, negative.values = FALSE, vcat = vcat, vcat.env.comps = vcat.env.comps,
vcat.actions = vcat.actions)
```
