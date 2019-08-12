##plotMA, but for ggplot
coef.vec <- c("Int_LIR_vs_RES", "nInt_LIR_vs_RES")
coef.ind <- c(2, 4)

coefv
coefi

plot.tbl <- tibble(ensembl_gene_id = rownames(all.fite$coefficients),
                   Amean = all.fite$Amean,
                   Coefs = all.fite$coefficients[,coefi],
                   Lods = all.fite$lods[,coefi],
                   Pval = all.fite$p.value[,coefi]) %>%
            dplyr:::mutate(mlog10Pval = -log10(Pval)) %>%
            left_join(rhtx2gene, .) %>%
            na.omit()

##status for MA plot
status_topn <- DEList[[coefi]] %>%
               left_join(plot.tbl, .) %>%
               na.omit() %>%
               top_n(n=5) %>%
               dplyr::select(external_gene_name) %>% unlist()
status_botn <- DEList[[coefi]] %>%
               left_join(plot.tbl, .) %>%
               na.omit() %>%
               top_n(n=-5) %>%
               dplyr::select(external_gene_name) %>% unlist()
status_n <- c(status_topn, status_botn)
plot.tbl <- plot.tbl %>% dplyr::mutate(status = ifelse(external_gene_name %in% status_n, 1, 0)) %>%
            dplyr::select(-human_ensembl_gene_id) %>% unique()
plot.tbl$status <- factor(plot.tbl$status)
plot.tbl.s <- plot.tbl %>% left_join(., DEList[[coefi]]) %>% na.omit()
plot.tbl.ss <- plot.tbl.s %>% dplyr::filter(status %in% 1)

#     plotWithHighlights(x=plot.tbl$Amean, y=plot.tbl$Coefs, main=coef.vec[1], status=plot.tbl$status, values="1")
#     text(x=plot.tbl$Amean[plot.tbl$status == 1],
#          y=plot.tbl$Coefs[plot.tbl$status == 1],
#          plot.tbl$external_gene_name[plot.tbl$status == 1],
#          col="red", pos=1)
mapl <- ggplot(plot.tbl, aes(x=Amean, y=Coefs, label=external_gene_name)) +
geom_point(size=0.5) +
geom_point(data=plot.tbl.ss, colour="red", size=0.5) +
geom_text_repel(data=plot.tbl.ss, colour="red", segment.size=0.2, segment.colour="red", fontface = "bold", cex=3) +
labs(x="Average Mean Expression", y="Log Fold Change", title=paste0("MA plot: ", coefv))
ggsave(mapl, filename=paste0(OUTDIR, "/MA.", coefv, ".pdf"))




}
