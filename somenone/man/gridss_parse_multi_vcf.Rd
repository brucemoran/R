% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/somenone_gridss.R
\name{gridss_parse_multi_vcf}
\alias{gridss_parse_multi_vcf}
\title{GRIDSS functions}
\usage{
gridss_parse_multi_vcf(vcf, which_genome = NULL)
}
\arguments{
\item{vcf, }{the VCF file path}

\item{which_genome, }{genome assembly used ("hg19", "hg38")}
}
\value{
a list with recurrent SVs and private SVs as elements
}
\description{
Parse multi-sample VCF (recommended output of GRIDSS)
}
