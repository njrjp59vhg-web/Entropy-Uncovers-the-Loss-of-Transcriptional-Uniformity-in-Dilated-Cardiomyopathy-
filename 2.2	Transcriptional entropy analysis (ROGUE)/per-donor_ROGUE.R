library(ROGUE)
library(Matrix)
library(tibble)
library(dplyr)

pathway_expr <- readMM("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/Science/expr_lv.mtx")
pathway_cells <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/Science/cells_lv.csv")
pathway_genes <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/Science/genes_lv.csv")
pathway_meta <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/Science/meta_lv.csv")

rownames(pathway_expr) <- pathway_genes[[1]]
colnames(pathway_expr) <- pathway_cells[[1]]

pathway_expr <- as.matrix(pathway_expr)
pathway_expr <- matr.filter(pathway_expr, min.cells = 50, min.genes = 10)
pathway_meta <- pathway_meta[pathway_meta$X %in% colnames(pathway_expr), ]

rogue.pathway <- rogue(
  pathway_expr,
  labels = pathway_meta$pathway,
  samples = pathway_meta$donor_id,
  platform = "UMI",
  filter = FALSE,
  span = 0.75
)

rogue.states <- rogue.states <- rogue(
  pathway_expr ,
  labels = pathway_meta$cell_states,
  samples = pathway_meta$donor_id,
  platform = "UMI",
  filter = FALSE,
  span = 0.75
)

saveRDS(rogue.pathway, "/rds/general/user/ao225/home/CardiaFinal/Results/Science/rogue_pathway.rds")
saveRDS(rogue.states, "/rds/general/user/ao225/home/CardiaFinal/Results/Science/rogue_pathway.rds")

