### R/shapeQTL: shape QTL mapping experiment with R
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)

[Nicolas Navarro](http://nnavarro.free.fr)

[R/shapeQTL](http://nnavarro.free.fr/programs.html) is an [R](http://www.r-project.org) package to perform 
statistical analyses to map quantitative trait locus for geometric morphometric data

Geometric morphometrics defines shape as a multivariate trait. Gene pleiotropy on landmark coordinates is the 
underlying genetic model by definition and multivariate methods are required to map QTL. The package provides 
Haley-Knott mapping for shape data derived from other R packages: [R/geomorph](https://github.com/cran/geomorph), 
[R/Morpho](https://github.com/zarquon42b/Morpho), [R/Momocs](https://github.com/vbonhomme/Momocs/) or other open source 
softwares for geometric morphometrics : e.g., [Java/MorphoJ](http://flywings.org.uk/MorphoJ_page.htm), and for genotype probabilities derived from specific packages 
depending on the type of crosses (for example calc.genoprob from [R/qtl](http://www.rqtl.org) for inbred strains).

#### Installation

Install R/shapeQTL from its [GitHub repository](https://github.com/nnavarro/shapeQTL).

##### Install prerequisites
Install [devtools](https://github.com/hadley/devtools) package.

```r
install.packages("devtools")
```

Then install R/shapeQTL 

```r
	require(devtools)
	install_github("nnavarro/shapeQTL", local=FALSE)
```
