mask <- sce.tmp$Treatment == ref.treatment
pw_activities.ref <- pw_activities[, mask]

tmp <- pw_activities.ref %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "Pathway") %>%
  pivot_longer(cols = c(-Pathway), names_to = "CellID", values_to = "Activity") %>%
  left_join(y = colData(sce.tmp) %>%
              as_tibble(rownames = NA) %>%
              rownames_to_column(var = ".CellID") %>%
              select(.CellID, Cluster, Treatment, BarcodeType),
            by = c("CellID" = ".CellID")) %>%
  filter(Treatment == ref.treatment) %>%
  group_by(Cluster) %>%
  filter(n() >= 50) %>%
  ungroup() %>%
  filter(Cluster %in% c(8, 21, 4, 20, 6, 17)) %>%
  filter(Cluster %in% c(4, 8)) %>%
  ggplot(aes(x = Pathway, y = Activity, fill = BarcodeType)) + 
  geom_boxplot() + 
  facet_wrap(~Cluster, nrow = 2)

pw_activities %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "Pathway") %>%
  pivot_longer(cols = c(-Pathway), names_to = "CellID", values_to = "Activity") %>%
  left_join(y = colData(sce.tmp) %>%
              as_tibble(rownames = NA) %>%
              rownames_to_column(var = ".CellID") %>%
              select(.CellID, Cluster, Treatment, BarcodeType),
            by = c("CellID" = ".CellID")) %>%
  #filter(Cluster %in% c(4, 8)) %>%
  filter(BarcodeType != "Invalid barcode", Cluster %in% c(4, 7, 5, 2, 8, 14), Pathway == "TGFb") %>%
  mutate(Class = sprintf("Cluster %d: %s", Cluster, BarcodeType),
         Cluster = factor(Cluster, levels = c(4,7,5,2,8,14))) %>%
  arrange(Cluster, BarcodeType) %>%
  ggplot(aes(x = Cluster, y = Activity, fill = BarcodeType)) + 
  geom_boxplot() 


sce.tmp$TGFb <- pw_activities["TGFb",]
plotReducedDim(sce.tmp, dimred = "TSNE", colour_by = "TGFb", text_by = "Cluster", point_size = 1, point_alpha = 0.7) +
  ggtitle("TGFb activity") 



ggplot(data = tmp %>% filter(Pathway == "TGFb"), aes(x = Pathway, y = Activity, fill = Class)) + 
  geom_boxplot() 






##### Venn diagrams #####
bc.type <- "Depleted"
enr.A <- ff$Amgen %>%
  filter(Type == bc.type) %>%
  select(BC_ID) %>%
  pull()
enr.M <- ff$Mirati %>%
  filter(Type == bc.type) %>%
  select(BC_ID) %>%
  pull()
enr.R <- ff$RevMed %>%
  filter(Type == bc.type) %>%
  select(BC_ID) %>%
  pull()
length(enr.A)
length(enr.M)
length(enr.R)
length(intersect(enr.A, enr.M))
length(intersect(enr.A, enr.R))
length(intersect(enr.M, enr.R))
length(intersect(intersect(enr.A, enr.M), enr.R))
venn.diagram(list("Amgen" = enr.A, "Mirati" = enr.M, "RevMed" = enr.R), filename = "~/venn.tiff")