library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(ggpubr)

rogue.pathway <- readRDS("/rds/general/user/ao225/home/ext/rogue_pathway.rds")

# wide (donors x pathways) -> long, one row per donor
rp <- rogue.pathway |>
  as.data.frame() |>
  rownames_to_column("donor") |>
  pivot_longer(-donor, names_to = "pathway", values_to = "rogue") |>
  filter(!is.na(rogue))

lev <- c("control", "splicing", "sarcomere", "PVneg", "chromatin")
rp$pathway <- factor(rp$pathway, levels = lev)
print(table(rp$pathway))                      

## ---- stats ----
kruskal.test(rogue ~ pathway, data = rp)      

# each DCM group vs control (control is first level -> DCM<control is "greater")
dcm <- setdiff(lev, "control")
vs <- lapply(dcm, function(g)
  wilcox.test(rogue ~ pathway,
              data = droplevels(subset(rp, pathway %in% c("control", g))),
              alternative = "greater"))
print(data.frame(pathway = dcm,
                 p = sapply(vs, `[[`, "p.value"),
                 padj = p.adjust(sapply(vs, `[[`, "p.value"), "BH")))

# LMNA vs all other DCM
rp2 <- subset(rp, pathway != "control")
rp2$grp <- ifelse(rp2$pathway == "chromatin", "LMNA", "otherDCM")
print(wilcox.test(rogue ~ grp, rp2, alternative = "less"))

print(rp |> group_by(pathway) |> summarise(median = median(rogue), n = n()))

dcm <- c("splicing", "sarcomere", "PVneg", "chromatin")
pv  <- sapply(dcm, function(g)
  wilcox.test(rogue ~ pathway,
              data = droplevels(subset(rp, pathway %in% c("control", g))),
              alternative = "greater")$p.value)
padj <- p.adjust(pv, "BH")
star <- cut(padj, c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***","**","*","ns"))
lab  <- data.frame(pathway = factor(dcm, levels = lev),
                   star = as.character(star), y = 0.95)

# --- plot ---
p <- ggplot(rp, aes(pathway, rogue, fill = pathway)) +
  geom_boxplot(width = .55, outlier.shape = NA, alpha = .55) +
  geom_jitter(width = .12, height = 0, size = 2, alpha = .7) +
  geom_text(data = lab, aes(x = pathway, y = y, label = star),
            inherit.aes = FALSE, size = 6) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c("Control","Splicing\n(RBM20)","Sarcomere",
                              "PVneg","Chromatin\n(LMNA)")) +
  labs(x = NULL, y = "ROGUE  (higher = more uniform)",
       title = "Per-donor cardiomyocyte transcriptional uniformity — LV, v3",
       subtitle = "Each point = one donor;  * vs control (one-sided Wilcoxon, BH-adjusted)") +
  theme_classic(base_size = 13) + theme(legend.position = "none")

ggsave("/rds/general/user/ao225/home/ext/rogue_lv_v3.png", p, width = 8, height = 5.5, dpi = 300)
print(p)
