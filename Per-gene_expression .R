library(ROGUE)
library(Matrix)
library(tibble)
library(dplyr)

expr <- readMM("/rds/general/user/ao225/home/ext/Science_variant/expr_control.mtx")
cells <- read.csv("/rds/general/user/ao225/home/ext/Science_variant/cells_control.csv")
genes <- read.csv("/rds/general/user/ao225/home/ext/Science_variant/genes_control.csv")

rownames(expr) <- genes[[1]]
colnames(expr) <- cells[[1]]

expr <- as.matrix(expr)
expr <- matr.filter(expr, min.cells = 50, min.genes = 10)

ent = SE_fun(expr)
saveRDS(ent,"/rds/general/user/ao225/home/ext/result/ent_sarcomere.rds")

expr <- readMM("/rds/general/user/ao225/home/ext/Science_variant/expr_sarcomere.mtx")
cells <- read.csv("/rds/general/user/ao225/home/ext/Science_variant/cells_sarcomere.csv")
genes <- read.csv("/rds/general/user/ao225/home/ext/Science_variant/genes_sarcomere.csv")

rownames(expr) <- genes[[1]]
colnames(expr) <- cells[[1]]

expr <- as.matrix(expr)
expr <- matr.filter(expr, min.cells = 50, min.genes = 10)

ent = SE_fun(expr)
saveRDS(ent,"/rds/general/user/ao225/home/ext/result/ent_sarcomere.rds")

expr <- readMM("/rds/general/user/ao225/home/ext/Science_variant/expr_chromatin.mtx")
cells <- read.csv("/rds/general/user/ao225/home/ext/Science_variant/cells_chromatin.csv")
genes <- read.csv("/rds/general/user/ao225/home/ext/Science_variant/genes_chromatin.csv")

rownames(expr) <- genes[[1]]
colnames(expr) <- cells[[1]]

expr <- as.matrix(expr)
expr <- matr.filter(expr, min.cells = 50, min.genes = 10)

ent = SE_fun(expr)
saveRDS(ent,"/rds/general/user/ao225/home/ext/result/ent_chromatin.rds")

expr <- readMM("/rds/general/user/ao225/home/ext/Science_variant/expr_splicing.mtx")
cells <- read.csv("/rds/general/user/ao225/home/ext/Science_variant/cells_splicing.csv")
genes <- read.csv("/rds/general/user/ao225/home/ext/Science_variant/genes_splicing.csv")

rownames(expr) <- genes[[1]]
colnames(expr) <- cells[[1]]

expr <- as.matrix(expr)
expr <- matr.filter(expr, min.cells = 50, min.genes = 10)

ent = SE_fun(expr)
saveRDS(ent,"/rds/general/user/ao225/home/ext/result/ent_splicing.rds")

expr <- readMM("/rds/general/user/ao225/home/ext/Science_variant/expr_PVneg.mtx")
cells <- read.csv("/rds/general/user/ao225/home/ext/Science_variant/cells_PVneg.csv")
genes <- read.csv("/rds/general/user/ao225/home/ext/Science_variant/genes_PVneg.csv")

rownames(expr) <- genes[[1]]
colnames(expr) <- cells[[1]]

expr <- as.matrix(expr)
expr <- matr.filter(expr, min.cells = 50, min.genes = 10)

ent = SE_fun(expr)
saveRDS(ent,"/rds/general/user/ao225/home/ext/result/ent_PVneg.rds")