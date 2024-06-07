library(tidyverse)
library(scater)
library(scran)
library(SingleCellExperiment)
library(gridExtra)
library(DESeq2)


#### Identify over- and underrepresented CaTCH barcodes ####
metadata <- colData(sce) %>% 
            as_tibble() %>% 
            dplyr::select(Sample, Treatment) %>% 
            distinct()
ref.treatment <- levels(metadata$Treatment)[1]

bc.counts <- colData(sce) %>% 
             as_tibble() %>%
             filter(CaTCH.Status == "Singlet") %>%
             dplyr::select(CaTCH.BCs, Sample, BC_ID) %>%
             group_by(CaTCH.BCs, Sample) %>%
             mutate(n = n()) %>%
             ungroup() %>%
             distinct() %>%
             pivot_wider(names_from = Sample, values_from = n, values_fill = 0) %>%
             filter(rowSums(across(starts_with(ref.treatment))) > 0)

idx <- match(metadata$Sample, setdiff(colnames(bc.counts), c("CaTCH.BCs", "BC_ID")))
metadata <- metadata[idx, ]

dds <- DESeqDataSetFromMatrix(countData = bc.counts %>%
                                          dplyr::select(-CaTCH.BCs) %>%
                                          column_to_rownames(var = "BC_ID"),
                              colData = metadata,
                              design= ~ Treatment)
dds <- DESeq(dds)



#### Get the expression values for the genes of interest ####
genes.of.interest <- c("dCas9-VPR", "GFP", "tagBFP", 
                       "KRAS", "TP53", "DUSP6")
GoI.expr <- tibble(BC_ID = sce$BC_ID,
                   Treatment = sce$Treatment)
for (g in genes.of.interest) {
  GoI.expr[[paste0("Expr.", g)]] <- assay(sce, "logcounts")[g,]
}


#### Plot the expression of the gene of interest for each treatment ####
for (t in setdiff(levels(metadata$Treatment), ref.treatment)) {
  print(sprintf("Processing the treatment '%s'...", t))
  r <- results(dds, 
               alpha = 0.05,
               contrast = c("Treatment", t, ref.treatment)) %>%
       as_tibble(rownames = NA) %>% 
       rownames_to_column("BC_ID") %>%
       filter(!is.na(padj), padj <= 0.05) %>%
       mutate(Type = if_else(log2FoldChange > 0, "Enriched", "Depleted")) %>%
       arrange(desc(Type), desc(abs(log2FoldChange)))

  types <- GoI.expr %>%
           left_join(y = (r %>% 
                          dplyr::select(BC_ID, Type)), 
                     by = "BC_ID") %>%
           dplyr::select(Type)

  cbind(GoI.expr, types) %>%
  as_tibble() %>%
  filter(Treatment %in% c(ref.treatment, t)) %>%
  pivot_longer(cols = starts_with("Expr."), names_to = "Gene", values_to = "LogCounts") %>%
  mutate(Gene = str_replace(string = Gene, pattern = "Expr.", replacement = "")) %>%
  ggplot(aes(x = Type, y = LogCounts, fill = Treatment)) + 
    geom_boxplot() + 
    facet_wrap(~Gene)
  ggsave(filename = paste0("~/plots/expr_", t, ".jpeg"), 
         width = 2600, 
         height = 1500, 
         units = "px", 
         device = "jpeg")
}


#enr <- sce$VPR.expr[which(sce$.BCType == "Enriched")]
#dep <- sce$VPR.expr[which(sce$.BCType == "Depleted")]
#t.test(enr, dep)