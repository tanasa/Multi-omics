cat("

PANDA is supervised (discriminant), not unsupervised.

PANDA = Pan-omics Discriminant Analysis.

It explicitly uses class labels (phenotypes / subtypes) during training to learn latent components that separate groups.

In the code : PanDAModel <- PanDA(data, labels, numComponents)

* labels = known outcomes (e.g. breast cancer subtypes)

* These labels guide the optimization

* The latent components are discriminant components, not just variance-capturing ones

")

cat("

PANDA jointly learns:

* Shared latent space across omics

* Maximal class separation (supervised signal)

* Low correlation between omics-specific components.

")

cat("

| Method    | Uses labels? | Goal                                         |
| --------- | ------------ | -------------------------------------------- |
| PCA       | ❌ No         | Maximize variance                            |
| MOFA      | ❌ No         | Explain covariance                           |
| SNF       | ❌ No         | Sample similarity                            |
| **PANDA** | ✅ Yes        | **Class separation + cross-omics alignment** |
| PLS-DA    | ✅ Yes        | Discriminant projection                      |


")

cat("Important nuance (very important) :

PANDA is:

Supervised during training

But produces latent representations that you can:

visualize

cluster

use in downstream models (XGBoost, survival, etc.)

")

cat("

What PANDA is based on

PANDA is a supervised, discriminant, projection-based method.

Conceptually, it is closer to:

LDA / PLS-DA

Supervised CCA–style methods

Multi-view discriminant analysis

It learns linear projection matrices that map each omics dataset into a shared, label-informed latent space.

There is no non-negativity constraint in PANDA.

")

cat(" 

PANDA = supervised multi-omics discriminant projection

in contrast :

NMF / IntNMF = unsupervised nonnegative matrix factorization

")

cat("

PanDA’s primary purpose is clinical stratification.

PanDA is designed to identify clinically relevant patient subgroups by integrating multiple omics layers in a supervised way.

It does not discover biology blindly.
It learns discriminant signals that separate known clinical states.

")

cat("In contrast :

PLS = Partial Least Squares

PLS is a supervised dimensionality-reduction method that finds latent components in predictors (X) that best explain an outcome (Y).

PLS finds directions in high-dimensional data that are maximally predictive of an outcome.

PLS learns latent components such that:

They are linear combinations of X

They maximize covariance with Y

They reduce X to a small number of components

Mathematically (conceptually):

𝑋 → 𝑇 = 𝑋 𝑊

such that cov(𝑇,𝑌) is maximized

")

cat(" 

On the side note:  Why PLS exists (the problem it solves)

PLS was designed for situations where:

p ≫ n (many more features than samples)
Features are highly correlated (omics!)
Classical regression fails or overfits

Typical omics scenario:

100 patients
20,000 genes
→ standard regression breaks
→ PLS works

| PCA                     | PLS                             |
| ----------------------- | ------------------------------- |
| Unsupervised            | **Supervised**                  |
| Maximizes variance in X | **Maximizes covariance with Y** |
| Ignores outcome         | **Uses outcome directly**       |
| Discovery               | **Prediction / stratification** |

👉 PCA answers: “What varies most?”
👉 PLS answers: “What predicts the outcome?”

| Method    | Question answered                                                  |
| --------- | ------------------------------------------------------------------ |
| PCA       | What varies most overall?                                          |
| PLS       | What predicts Y best?                                              |
| **PanDA** | **What multi-omics patterns *define clinically distinct groups*?** |

")

## =========================
## Libraries
## =========================
library(corrplot)
library(kernlab)
library(IntNMF)
library(mclust)
library(aricode)
library(ggforce)
library(magrittr)
library(Seurat)
library(ggpubr)
library(pheatmap)
library(ConsensusClusterPlus)

library(PANDA)
library(mixOmics)
library(plotly)


## =========================
## Load data
## =========================
data(breast.TCGA) # TCGA data from mixOmics

data_list <- list()
data_list$mirna   <- t(breast.TCGA$data.train$mirna)
data_list$mrna    <- t(breast.TCGA$data.train$mrna)
data_list$protein <- t(breast.TCGA$data.train$protein)

Y <- breast.TCGA$data.train$subtype
subtype <- factor(Y)

## =========================
## Fit PanDA
## =========================
gnd <- as.numeric(Y)
numComponents <- 10

PanDAModel <- PanDA(data_list, gnd, numComponents, 0.2)

## =========================
## Extract discriminant components (train)
## =========================
mirnaComp   <- as.data.frame(PanDAModel[["PanDAComponents"]][["mirnaComponents"]])
mrnaComp    <- as.data.frame(PanDAModel[["PanDAComponents"]][["mrnaComponents"]])
proteinComp <- as.data.frame(PanDAModel[["PanDAComponents"]][["proteinComponents"]])

head(data_list$mirna)
dim(data_list$mirna)

head(data_list$mrna)
dim(data_list$mrna)

unique(Y)

head(data_list$protein) 
dim(data_list$protein) 

head(mirnaComp)

# Rows = samples (patients / tumors)
# Columns (DC 1 … DC 10) = PanDA discriminant components
# Values = the coordinates of each sample in the learned latent space for the miRNA modality
# Mathematically, this is the projected data:
# 𝑍(miRNA) = 𝑋(miRNA) 𝑊(miRNA)


head(mrnaComp)

# Rows = samples (patients / tumors)
# Columns (DC 1 … DC 10) = PanDA discriminant components
# Values = the coordinates of each sample in the learned latent space for the mRNA modality
# Mathematically, this is the projected data:
# 𝑍(mRNA) = 𝑋(mRNA) 𝑊(mRNA)


head(proteinComp)

# Rows = samples (patients / tumors)
# Columns (DC 1 … DC 10) = PanDA discriminant components
# Values = the coordinates of each sample in the learned latent space for the protein modality
# Mathematically, this is the projected data:
# 𝑍(protein) = 𝑋(protein) 𝑊(protein)


## If you get: object 'DC 1' not found, uncomment this rename block:
# names(mrnaComp)[1:3]    <- c("DC 1","DC 2","DC 3")
# names(mirnaComp)[1:3]   <- c("DC 1","DC 2","DC 3")
# names(proteinComp)[1:3] <- c("DC 1","DC 2","DC 3")

cols <- c('#BF382A', '#0C4B8E', "#fc8d59")

## 3D plots (train)
fig1 <- plot_ly(mrnaComp, x = ~`DC 1`, y = ~`DC 2`, z = ~`DC 3`,
                color = ~subtype, colors = cols) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'DC 1'),
                      yaxis = list(title = 'DC 2'),
                      zaxis = list(title = 'DC 3')))
# fig1

fig2 <- plot_ly(mirnaComp, x = ~`DC 1`, y = ~`DC 2`, z = ~`DC 3`,
                color = ~subtype, colors = cols) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'DC 1'),
                      yaxis = list(title = 'DC 2'),
                      zaxis = list(title = 'DC 3')))
# fig2

fig3 <- plot_ly(proteinComp, x = ~`DC 1`, y = ~`DC 2`, z = ~`DC 3`,
                color = ~subtype, colors = cols) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'DC 1'),
                      yaxis = list(title = 'DC 2'),
                      zaxis = list(title = 'DC 3')))
# fig3


fig1

fig2

fig3

## =========================
## Train-set clustering purity (ConsensusClusterPlus)
## =========================
subgroup <- 3

mRNAres <- ConsensusClusterPlus(
  t(as.matrix(mrnaComp[, c(1,2)])),
  maxK = 5, reps = 100, distance = "euclidean",
  tmyPal = c("white","seagreen"),
  clusterAlg = "hc",
  plot = FALSE
)
mRNAlabel <- unlist(mRNAres[[subgroup]]["consensusClass"])
ClusterPurity(gnd, mRNAlabel)

miRNAres <- ConsensusClusterPlus(
  t(as.matrix(mirnaComp[, c(1,2)])),
  maxK = 5, reps = 100, distance = "euclidean",
  tmyPal = c("white","seagreen"),
  clusterAlg = "hc",
  plot = FALSE
)
miRNAlabel <- unlist(miRNAres[[subgroup]]["consensusClass"])
ClusterPurity(gnd, miRNAlabel)

proteinres <- ConsensusClusterPlus(
  t(as.matrix(proteinComp[, c(1,2)])),
  maxK = 5, reps = 100, distance = "euclidean",
  tmyPal = c("white","seagreen"),
  clusterAlg = "hc",
  plot = FALSE
)
proteinlabel <- unlist(proteinres[[subgroup]]["consensusClass"])
ClusterPurity(gnd, proteinlabel)


## =========================
## Test-set projection + plots
## =========================
testdata <- list()
testdata$mirna <- t(breast.TCGA$data.test$mirna)
testdata$mrna  <- t(breast.TCGA$data.test$mrna)

testlabel <- breast.TCGA[["data.test"]][["subtype"]]
testgnd <- as.numeric(testlabel)

mirnatestComponent <- as.data.frame(t(testdata$mirna) %*% PanDAModel[["projMatrices"]][["Wmirna"]])
mirnatestComponent$subtype <- testgnd

mrnatestComponent <- as.data.frame(t(testdata$mrna) %*% PanDAModel[["projMatrices"]][["Wmrna"]])
mrnatestComponent$subtype <- testgnd

fig4 <- plot_ly(mirnatestComponent, x = ~V1, y = ~V2, z = ~V3,
                color = ~subtype, colors = cols) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'DC 1'),
                      yaxis = list(title = 'DC 2'),
                      zaxis = list(title = 'DC 3')))
# fig4

fig5 <- plot_ly(mrnatestComponent, x = ~V1, y = ~V2, z = ~V3,
                color = ~subtype, colors = cols) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'DC 1'),
                      yaxis = list(title = 'DC 2'),
                      zaxis = list(title = 'DC 3')))
# fig5


fig4

fig5



## =========================
## Test-set clustering purity
## =========================
tmRNAres <- ConsensusClusterPlus(
  t(as.matrix(mrnatestComponent[, c(1,2)])),
  maxK = 5, reps = 100, distance = "euclidean",
  tmyPal = c("white","seagreen"),
  clusterAlg = "hc",
  plot = FALSE
)
tmRNAlabel <- unlist(tmRNAres[[subgroup]]["consensusClass"])
ClusterPurity(testgnd, tmRNAlabel)

tmiRNAres <- ConsensusClusterPlus(
  t(as.matrix(mirnatestComponent[, c(1,2)])),
  maxK = 5, reps = 100, distance = "euclidean",
  tmyPal = c("white","seagreen"),
  clusterAlg = "hc",
  plot = FALSE
)
tmiRNAlabel <- unlist(tmiRNAres[[subgroup]]["consensusClass"])
ClusterPurity(testgnd, tmiRNAlabel)

## =========================
## Correlation among extracted components (within mRNA)
## =========================
M <- cor(mrnaComp[, 1:10])
corrplot.mixed(M)

## =========================
## Cross-omics correlation (first component)
## =========================
combComponents <- cbind(mirnaComp[,1], mrnaComp[,1], proteinComp[,1])
colnames(combComponents) <- c("miRNA","mRNA","Protein")

panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits = 2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8 / strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
upper.panel <- function(x, y){
  points(x, y, pch = 19, col = cols[subtype], cex = 1.5)
}

pairs(combComponents,
      lower.panel = panel.cor,
      upper.panel = upper.panel)

## =========================
## Gene loadings (top markers) for mRNA components
## =========================
refComponents <- 1:3
nfea <- 30

loadingsmRNA <- as.data.frame(PanDAModel[["projMatrices"]][["Wmrna"]][, refComponents])
rownames(loadingsmRNA) <- rownames(data_list[["mrna"]])

tempfeamRNA <- data.frame(matrix(ncol = length(refComponents), nrow = nfea))
colnames(tempfeamRNA) <- paste("DC", refComponents, sep = "_")

genes <- data.frame(matrix(ncol = length(refComponents), nrow = nfea))
colnames(genes) <- paste("DC", refComponents, sep = "_")

signmRNA <- data.frame(matrix(ncol = length(refComponents), nrow = nfea))
colnames(signmRNA) <- paste("DC", refComponents, sep = "_")

for (i in seq_along(refComponents)) {
  ind <- sort(abs(loadingsmRNA[, i]), decreasing = TRUE, index.return = TRUE)$ix
  tempfeamRNA[, i] <- loadingsmRNA[ind[1:nfea], i]
  genes[, i] <- rownames(loadingsmRNA)[ind[1:nfea]]
  for (j in 1:nrow(tempfeamRNA)) {
    signmRNA[j, i] <- ifelse(tempfeamRNA[j, i] > 0, "+", "-")
  }
}

## Plot top genes on component 1
contrib <- "max"
X2 <- breast.TCGA$data.train$mrna

method.group <- list()
cindx <- data.frame(matrix(NA, ncol = 1, nrow = ncol(X2), dimnames = list(colnames(X2), "clusters")))

for (k in 1:ncol(X2)) {
  method.group[[k]] <- tapply(X2[, k], gnd, mean, na.rm = TRUE)
  cindx[k, ] <- which.max(method.group[[k]])
}
cindx[cindx == 1] <- "Basal"
cindx[cindx == 2] <- "Her2"
cindx[cindx == 3] <- "LumA"

ffeamRNA <- data.frame(matrix(ncol = 4, nrow = nfea))
component <- "DC_1"
ffeamRNA[, 1] <- genes[, component]
ffeamRNA[, 2] <- abs(tempfeamRNA[, component])
ffeamRNA[, 3] <- signmRNA[, component]
ffeamRNA[, 4] <- cindx[genes[, component], ]
colnames(ffeamRNA) <- c("genes","weights","sign","clusters")

ggdotchart(
  ffeamRNA, x = "genes", y = "weights",
  color = "clusters",
  palette = cols,
  sorting = "descending",
  add = "segments",
  add.params = list(color = "#999999", size = 3),
  rotate = TRUE,
  dot.size = 6,
  label = ffeamRNA$sign,
  font.label = list(face = "bold", color = "black", size = 12, vjust = 0.4),
  ggtheme = theme_pubr()
) +
  font("x.text", size = 10, color = "black", face = "bold.italic") +
  font("y.text", size = 10, color = "black", face = "bold.italic") +
  font("xy", size = 10, color = "black", face = "bold.italic")

## =========================
## Heatmap of top genes on component 1
## =========================
nredmRNA <- X2[, genes[, component]]

outmRNA <- pheatmap(nredmRNA, cluster_rows = TRUE, scale = "row",
                    show_rownames = FALSE, cutree_rows = 3, cutree_cols = 2)

mrnalabel <- cutree(outmRNA$tree_row, k = 3)[outmRNA$tree_row[["order"]]]

annot_mrna <- data.frame(
  row.names = names(mrnalabel),
  Clusters = as.factor(mrnalabel),
  Groundtruth = Y[outmRNA[["tree_row"]][["order"]]]
)

my_colour <- list(
  Groundtruth = c(Basal = "#BF382A", Her2 = "#0C4B8E", LumA = "#fc8d59"),
  clustering = c("1" = "#CD0BBC", "2" = "#F5C710", "3" = "#28E2E5")
)

pheatmap(nredmRNA, cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, scale = "row", show_colnames = TRUE,
         cutree_rows = 3, annotation_row = annot_mrna,
         annotation_colors = my_colour)

## =========================
## XGBoost classification
## =========================
library(xgboost)
library(caret)

numberOfClasses <- length(unique(gnd))
xgb_params <- list(
  objective = "multi:softprob",
  eval_metric = "mlogloss",
  num_class = numberOfClasses
)
nround <- 50

train_data <- as.matrix(cbind(mrnaComp, mirnaComp))
train_label <- gnd
train_label[train_label == 3] <- 0
train_matrix <- xgb.DMatrix(data = train_data, label = train_label)

test_data <- as.matrix(cbind(
  mrnatestComponent[, 1:(ncol(mrnatestComponent) - 1)],
  mirnatestComponent[, 1:(ncol(mirnatestComponent) - 1)]
))
colnames(test_data) <- colnames(train_data)

test_label <- testgnd
test_label[test_label == 3] <- 0
test_matrix <- xgb.DMatrix(data = test_data, label = test_label)

bst_model <- xgb.train(params = xgb_params,
                       data = train_matrix,
                       nrounds = nround)

test_pred <- predict(bst_model, newdata = test_matrix)

test_prediction <- matrix(test_pred,
                          nrow = numberOfClasses,
                          ncol = length(test_pred) / numberOfClasses) %>%
  t() %>%
  data.frame() %>%
  mutate(label = test_label + 1,
         max_prob = max.col(., "last"))

confusionMatrix(factor(test_prediction$max_prob),
                factor(test_prediction$label),
                mode = "everything")

## =========================
## Test-set clustering purity
## =========================
tmRNAres <- ConsensusClusterPlus(
  t(as.matrix(mrnatestComponent[, c(1,2)])),
  maxK = 5, reps = 100, distance = "euclidean",
  tmyPal = c("white","seagreen"),
  clusterAlg = "hc",
  plot = FALSE
)
tmRNAlabel <- unlist(tmRNAres[[subgroup]]["consensusClass"])
ClusterPurity(testgnd, tmRNAlabel)

tmiRNAres <- ConsensusClusterPlus(
  t(as.matrix(mirnatestComponent[, c(1,2)])),
  maxK = 5, reps = 100, distance = "euclidean",
  tmyPal = c("white","seagreen"),
  clusterAlg = "hc",
  plot = FALSE
)
tmiRNAlabel <- unlist(tmiRNAres[[subgroup]]["consensusClass"])
ClusterPurity(testgnd, tmiRNAlabel)

## =========================
## Correlation among extracted components (within mRNA)
## =========================
M <- cor(mrnaComp[, 1:10])
corrplot.mixed(M)


# mirnaComp[,1]
# mrnaComp[,1]
# proteinComp[,1]

## =========================
## Cross-omics correlation (first component)
## =========================
combComponents <- cbind(mirnaComp[,1], mrnaComp[,1], proteinComp[,1])
colnames(combComponents) <- c("miRNA","mRNA","Protein")

panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits = 2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8 / strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
upper.panel <- function(x, y){
  points(x, y, pch = 19, col = cols[subtype], cex = 1.5)
}

pairs(combComponents,
      lower.panel = panel.cor,
      upper.panel = upper.panel)




## =========================
## Gene loadings (top markers) for mRNA components
## =========================
refComponents <- 1:3
nfea <- 30

loadingsmRNA <- as.data.frame(PanDAModel[["projMatrices"]][["Wmrna"]][, refComponents])
rownames(loadingsmRNA) <- rownames(data_list[["mrna"]])

tempfeamRNA <- data.frame(matrix(ncol = length(refComponents), nrow = nfea))
colnames(tempfeamRNA) <- paste("DC", refComponents, sep = "_")

genes <- data.frame(matrix(ncol = length(refComponents), nrow = nfea))
colnames(genes) <- paste("DC", refComponents, sep = "_")

signmRNA <- data.frame(matrix(ncol = length(refComponents), nrow = nfea))
colnames(signmRNA) <- paste("DC", refComponents, sep = "_")

for (i in seq_along(refComponents)) {
  ind <- sort(abs(loadingsmRNA[, i]), decreasing = TRUE, index.return = TRUE)$ix
  tempfeamRNA[, i] <- loadingsmRNA[ind[1:nfea], i]
  genes[, i] <- rownames(loadingsmRNA)[ind[1:nfea]]
  for (j in 1:nrow(tempfeamRNA)) {
    signmRNA[j, i] <- ifelse(tempfeamRNA[j, i] > 0, "+", "-")
  }
}

genes

## Plot top genes on component 1
contrib <- "max"
X2 <- breast.TCGA$data.train$mrna

method.group <- list()
cindx <- data.frame(matrix(NA, ncol = 1, nrow = ncol(X2), dimnames = list(colnames(X2), "clusters")))

for (k in 1:ncol(X2)) {
  method.group[[k]] <- tapply(X2[, k], gnd, mean, na.rm = TRUE)
  cindx[k, ] <- which.max(method.group[[k]])
}
cindx[cindx == 1] <- "Basal"
cindx[cindx == 2] <- "Her2"
cindx[cindx == 3] <- "LumA"

ffeamRNA <- data.frame(matrix(ncol = 4, nrow = nfea))
component <- "DC_1"
ffeamRNA[, 1] <- genes[, component]
ffeamRNA[, 2] <- abs(tempfeamRNA[, component])
ffeamRNA[, 3] <- signmRNA[, component]
ffeamRNA[, 4] <- cindx[genes[, component], ]
colnames(ffeamRNA) <- c("genes","weights","sign","clusters")

ggdotchart(
  ffeamRNA, x = "genes", y = "weights",
  color = "clusters",
  palette = cols,
  sorting = "descending",
  add = "segments",
  add.params = list(color = "#999999", size = 3),
  rotate = TRUE,
  dot.size = 6,
  label = ffeamRNA$sign,
  font.label = list(face = "bold", color = "black", size = 12, vjust = 0.4),
  ggtheme = theme_pubr()
) +
  font("x.text", size = 10, color = "black", face = "bold.italic") +
  font("y.text", size = 10, color = "black", face = "bold.italic") +
  font("xy", size = 10, color = "black", face = "bold.italic")


## =========================
## Heatmap of top genes on component 1
## =========================
nredmRNA <- X2[, genes[, component]]

outmRNA <- pheatmap(nredmRNA, cluster_rows = TRUE, scale = "row",
                    show_rownames = FALSE, cutree_rows = 3, cutree_cols = 2)

mrnalabel <- cutree(outmRNA$tree_row, k = 3)[outmRNA$tree_row[["order"]]]

annot_mrna <- data.frame(
  row.names = names(mrnalabel),
  Clusters = as.factor(mrnalabel),
  Groundtruth = Y[outmRNA[["tree_row"]][["order"]]]
)

my_colour <- list(
  Groundtruth = c(Basal = "#BF382A", Her2 = "#0C4B8E", LumA = "#fc8d59"),
  clustering = c("1" = "#CD0BBC", "2" = "#F5C710", "3" = "#28E2E5")
)

pheatmap(nredmRNA, cluster_rows = TRUE, cluster_cols = TRUE,
         show_rownames = FALSE, scale = "row", show_colnames = TRUE,
         cutree_rows = 3, annotation_row = annot_mrna,
         annotation_colors = my_colour)



## =========================
## XGBoost classification
## =========================
library(xgboost)
library(caret)

numberOfClasses <- length(unique(gnd))
xgb_params <- list(
  objective = "multi:softprob",
  eval_metric = "mlogloss",
  num_class = numberOfClasses
)
nround <- 50

train_data <- as.matrix(cbind(mrnaComp, mirnaComp))
train_label <- gnd
train_label[train_label == 3] <- 0
train_matrix <- xgb.DMatrix(data = train_data, label = train_label)

test_data <- as.matrix(cbind(
  mrnatestComponent[, 1:(ncol(mrnatestComponent) - 1)],
  mirnatestComponent[, 1:(ncol(mirnatestComponent) - 1)]
))
colnames(test_data) <- colnames(train_data)

test_label <- testgnd
test_label[test_label == 3] <- 0
test_matrix <- xgb.DMatrix(data = test_data, label = test_label)

bst_model <- xgb.train(params = xgb_params,
                       data = train_matrix,
                       nrounds = nround)

test_pred <- predict(bst_model, newdata = test_matrix)

test_prediction <- matrix(test_pred,
                          nrow = numberOfClasses,
                          ncol = length(test_pred) / numberOfClasses) %>%
  t() %>%
  data.frame() %>%
  mutate(label = test_label + 1,
         max_prob = max.col(., "last"))

confusionMatrix(factor(test_prediction$max_prob),
                factor(test_prediction$label),
                mode = "everything")


