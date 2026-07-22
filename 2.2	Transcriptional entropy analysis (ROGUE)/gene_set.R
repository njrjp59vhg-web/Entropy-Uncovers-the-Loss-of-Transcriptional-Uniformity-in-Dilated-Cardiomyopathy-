library(dplyr)

chromatin = readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs/ent_chromatin.rds")
sarcomere = readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs/ent_sarcomere.rds")
splicing  = readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs/ent_splicing.rds")
PVneg     = readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs/ent_PVneg.rds")
control   = readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs/ent_control.rds")

refset <- Reduce(union, list(chromatin$Gene, sarcomere$Gene,
                             splicing$Gene, PVneg$Gene, control$Gene))

stable_chromatin <- chromatin %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(ds)) %>%
  pull(Gene)

stable_sarcomere <- sarcomere %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(ds)) %>%
  pull(Gene)

stable_splicing <- splicing %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(ds)) %>%
  pull(Gene)

stable_PVneg <- PVneg %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(ds)) %>%
  pull(Gene)

stable_control<- control %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(ds)) %>%
  pull(Gene)

sharedDCM <- Reduce(union, list(stable_chromatin, stable_sarcomere, stable_splicing, stable_PVneg))
DCMonly   <- setdiff(sharedDCM, stable_control)

writeLines(DCMonly, "/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs/gene_set/DCMonly.txt")
writeLines(stable_chromatin, "/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs/gene_set/stable_chromatin.txt")
writeLines(stable_sarcomere, "/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs/gene_set/stable_sarcomere.txt")
writeLines(stable_splicing, "/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs/gene_set/stable_splicing.txt")
writeLines(stable_PVneg, "/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs/gene_set/stable_PVneg.txt")
writeLines(stable_control, "/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs/gene_set/stable_control.txt")
writeLines(refset, "/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs/gene_set/refset.txt")

dir <- "/rds/general/user/ao225/home/CardiaFinal/Results/Science/Rogue_outputs"
for (g in c("chromatin","sarcomere","splicing","PVneg","control")) {
  e <- readRDS(file.path(dir, paste0("ent_", g, ".rds")))
  write.csv(e[, c("Gene","ds","p.adj")],
            file.path(dir, paste0("ds_", g, ".csv")), row.names = FALSE)
}
