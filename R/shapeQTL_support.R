#' @name shapeQTL
#' @docType package
#' @aliases shapeQTL
#' @title shape QTL mapping experiment with R
#' @author Nicolas Navarro \email{nicolas.navarro@@ephe.sorbonne.fr}
#'
#' @description Functions in this package allow one to perform statistical analyses
#'  to map quantitative trait locus for geometric morphometric data
#' @details Geometric morphometrics defines shape as a multivariate trait. Gene pleiotropy on landmark coordinates is the underlying genetic model by definition and multivariate methods are required to map QTL. The package provides Haley-Knott mapping for shape data derived from other R packages: \code{\link[geomorph]{geomorph}}, \code{\link[Morpho]{Morpho}}, \code{\link[Momocs]{Momocs}} or other open source softwares for geometric morphometrics : e.g., \code{\link['http://www.flywings.org']{MorphoJ}}, and for genotype probabilities derived from specific packages depending on the type of crosses (for example \code{\link[qtl]{calc.genoprob}} from \code{R/qtl} for inbred strains).
#' 
#' @import qtl 
#' @import MASS
#' @import parallel
#' @useDynLib shapeQTL, .registration = TRUE, .fixes = "C_"
#' @references Haley CS, Knott SA (1992). A simple regression method for mapping quantitative trait loci in line crosses using flanking markers. Heredity 69: 315–324.
#' @references Knott SA, Haley CS (2000). Multitrait least squares for quantitative trait loci detection. Genetics 156: 899–911.
#' @references Klingenberg CP, Leamy LJ, Routman EJ, Cheverud JM (2001). Genetic architecture of mandible shape in mice: effects of quantitative trait loci analyzed by geometric morphometrics. Genetics 157: 785–802.
#' @references Klingenberg CP (2010). Evolution and development of shape: integrating quantitative approaches. Nature Reviews Genetics 11: 623–635.
#' @references Maga AM, Navarro N, Cunningham ML, Cox TC (2015). Quantitative trait loci affecting the 3D skull shape and size in mouse and prioritization of candidate genes in-silico. Front Physiol 6: 1–13. 
#' @examples XXXX
#' 
#' 
NULL

#' Landmark data from simulated drosophila wings
#'
#' @name wings
#' @docType data
#' @author Nicolas Navarro
#' @keywords datasets
NULL

#' Landmark data from a simulated BC of drosophila wings
#'
#' @name wings.bc
#' @docType data
#' @author Nicolas Navarro
#' @keywords datasets
NULL

#' Landmark data from a simulated F2 of drosophila wings
#'
#' @name wings.f2
#' @docType data
#' @author Nicolas Navarro
#' @keywords datasets
NULL