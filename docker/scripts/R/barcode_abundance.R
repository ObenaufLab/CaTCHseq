library(tidyverse)
library(scater)
library(scran)
library(SingleCellExperiment)
library(dorothea)
library(progeny)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(gridExtra)

load("/data/gcbds/users/nowoshil/Projects/CaTCH2.0/KRASi_ON_vs_OFF/Exp1010/_NF_run_/outputs/sce/sce.rda")
rm(sce.unfiltered)

ref.treatment <- levels(sce$Treatment)[1]

#### Barcode frequencies ####
if (!("BC_ID" %in% names(colData(sce)))) {
  tmp <- colData(sce) %>% 
         as_tibble() %>%
         filter(CaTCH.Status == "Singlet", Treatment == ref.treatment) %>%
         select(CaTCH.BCs, Sample) %>%
         group_by(CaTCH.BCs, Sample) %>%
         mutate(n = n()) %>%
         ungroup() %>%
         distinct() %>%
         pivot_wider(names_from = Sample, values_from = n, values_fill = 0) %>%
         filter(rowSums(across(starts_with(ref.treatment))) > 0) %>% 
         mutate(.Means = rowMeans(across(starts_with(ref.treatment)))) %>%
         arrange(by = desc(.Means)) %>%
         rowid_to_column(".ID") %>%
         mutate(BC_ID = paste0("BC_", .ID)) %>%
         select(-.Means, -.ID) %>%
         relocate(BC_ID, .after = CaTCH.BCs) %>%
         select(CaTCH.BCs, BC_ID)
  
  colData(sce)["BC_ID"] <- colData(sce) %>%
                           as_tibble() %>%
                           left_join(y = tmp, by = "CaTCH.BCs") %>%
                           mutate(BC_ID = factor(BC_ID, levels = str_sort(unique(BC_ID), numeric = TRUE))) %>%
                           select(BC_ID)
}

plots <- list()
for (t in levels(sce$Treatment)) {
  plot.data <- colData(sce) %>%
               as_tibble() %>%
               filter(CaTCH.Status == "Singlet", Treatment == t, !is.na(BC_ID)) %>% 
               select(Sample, BC_ID) %>%
               group_by(Sample, BC_ID) %>%
               mutate(n = n()) %>%
               ungroup() %>%
               arrange(Sample, desc(n)) %>%
               distinct() %>% 
               pivot_wider(names_from = "Sample", values_from = "n") %>%
               mutate(.Mean = round(rowMeans(across(starts_with(t)), na.rm = TRUE))) %>%
               arrange(desc(.Mean)) %>%
               rowid_to_column(var = ".id") %>%
               mutate(Label = if_else(.id <= 10, paste0(BC_ID, ": ", .Mean), NA_character_))
  p.bccounts <- ggplot(plot.data, aes(x = fct_inorder(as.character(BC_ID)), y = .Mean, label = Label)) + 
                 geom_col() + 
                 geom_label_repel(na.rm = TRUE, show.legend = FALSE) +
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       panel.grid.major.x = element_blank(),
                       panel.grid.minor.x = element_blank()) +
                 xlab("CaTCH barcodes") + 
                 ylab("Number of cells") + 
                 ggtitle(t)
  plots[[length(plots) + 1]] <- p.bccounts
  
  highlight.bcs <- plot.data %>%
                   filter(!is.na(Label)) %>%
                   select(BC_ID) %>%
                   pull()
  p.proportions <- colData(sce) %>%
                   as_tibble() %>%
                   filter(Treatment == t) %>%
                   select(BC_ID, CaTCH.Status, Cluster) %>%
                   group_by(Cluster) %>%
                   mutate(ClusterSize = n()) %>%
                   ungroup() %>%
                   filter(CaTCH.Status == "Singlet") %>%
                   select(-CaTCH.Status) %>%
                   filter(BC_ID %in% highlight.bcs) %>%
                   group_by(BC_ID, Cluster) %>%
                   mutate(BC_prop = n()) %>%
                   distinct() %>% 
                   ungroup() %>%
                   complete(Cluster) %>%
                   ggplot(aes(x = Cluster, fill = BC_ID)) + 
                    geom_col(aes(y = (BC_prop / ClusterSize) * 100)) +
                    ylab("Proportion of cells with\nspecified CaTCH barcodes") + 
                    labs(fill = "CaTCH barcode") + 
                    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  plots[[length(plots) + 1]] <- p.proportions
}

jpeg(filename = "~/plots/abundance.jpeg", width = 1200, height = length(levels(sce$Treatment)) * 300)
do.call("grid.arrange", c(plots, ncol = 2))
dev.off()
