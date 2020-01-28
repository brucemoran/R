##from GWASTools packages
library(tidyverse)
library(GenomicRanges)

UCSC_seqlevel <- paste0("chr", c(1:22,"X","Y"))
left.centromere.hg19 <- c(121535434,92326171,90504854,49660117,46405641,58830166,58054331,43838887,47367679,39254935,51644205,34856694,16000000,16000000,17000000,35335801,22263006,15460898,24681782,26369569,11288129,13000000,58632012,10104553)
right.centromere.hg19 <- c(124535434,95326171,93504854,52660117,49405641,61830166,61054331,46838887,50367679,42254935,54644205,37856694,19000000,19000000,20000000,38335801,25263006,18460898,27681782,29369569,14288129,16000000,61632012,13104553)
left.centromere.hg38 <- c(122026460,92188146,90772459,49708101,46485901,58553889,58169654,44033745,43236168,39686683,51078349,34769408,16000001,16000001,17000001,36311159,22813680,15460900,24498981,26436233,10864561,12954789,58605580,10316945)
right.centromere.hg38 <- c(125184587,94090557,93655574,51743951,50059807,59829934,60828234,45877265,45518558,41593521,54425074,37185252,18051248,18173523,19725254,38280682,26885980,20861206,27190874,30038348,12915808,15054318,62412542,10544039)

hg19_tib <- tibble(UCSC_seqlevel, left.centromere.hg19, right.centromere.hg19)
hg19_chrominfo <- fetchExtendedChromInfoFromUCSC("hg19") %>%
                  dplyr::filter(UCSC_seqlevel %in% hg19_tib$UCSC_seqlevel) %>%
                  left_join(., hg19_tib)
hg19_seqinfo <- Seqinfo(seqnames=hg19_chrominfo$UCSC_seqlevel,
                        seqlengths=hg19_chrominfo$UCSC_seqlength, isCircular=hg19_chrominfo$circular,
                        genome="hg19")
hg38_tib <- tibble(UCSC_seqlevel, left.centromere.hg38, right.centromere.hg38)
hg38_chrominfo <- fetchExtendedChromInfoFromUCSC("hg38") %>%
                  dplyr::filter(UCSC_seqlevel %in% hg38_tib$UCSC_seqlevel) %>%
                  left_join(., hg38_tib)
hg38_seqinfo <- Seqinfo(seqnames=hg38_chrominfo$UCSC_seqlevel,
                        seqlengths=hg38_chrominfo$UCSC_seqlength, isCircular=hg38_chrominfo$circular,
                        genome="hg38")
rm(UCSC_seqlevel, left.centromere.hg19, right.centromere.hg19, left.centromere.hg38, right.centromere.hg38, hg19_tib, hg38_tib)
