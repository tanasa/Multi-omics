# http://compbio.cs.toronto.edu/SNF/SNF/Software.html

# install.packages("heatmap.plus")
# devtools::install_github("maxconway/SNFtool")

# SNF itself is UNSUPERVISED.

# SNF (Similarity Network Fusion) does not use labels.
# It does not learn from outcomes or classes.
# It only fuses sample–sample similarity networks built from multiple data modalities.

cat("Similarity Network Fusion (SNF) is a new computational method for data integration.  
    
Briefly, SNF combines many different types of measurements (such as mRNA expression data, DNA methylation, miRNA expression 
    
and more - clinical data, questionnaires, image data, etc) for a given set of samples (e.g. patients). 
    
SNF first constructs a sample similarity network for each of the data types and then iteratively integrates these networks 
    
using a novel network fusion method. 
    
Working in the sample network space allows SNF to avoid dealing with different scale, 

collection bias and noise in different data types. 
    
Integrating data in a non-linear fashion allows SNF to take advantage of the common as well as complementary information 
    
in different data types.")

cat("

An example on the task of fusing two data types: mRNA expression and DNA methylation for the same cohort of patients. 

1. Construct patient similarity matrices for each data type using pairwise correlation   

2. Patient similarity matrices are equivalent to patient similarity networks where patients 
are nodes and edges represent patients’ pairwise similarities. 

3. Starting with the patient similarity networks run the patient network fusion, 
iteratively updating each of the networks with the information from the other networks, making them more similar with each step. 

4. The final fused network of patients to which the SNF process has converged. 
Edge color indicates which data type has contributed to the given similarity. 

")

# Notes :

# The code uses a function affinityMatrix to convert the Distance Matrix into an Affinity Matrix. 
# The standard mathematical method for this in spectral clustering (which SNF is based on) is the Gaussian Kernel 
# (also known as the Radial Basis Function or Heat Kernel)

library(SNFtool)

# ----------------------------------------------------------
# 1) Parameters
# ----------------------------------------------------------
K     <- 20   # number of neighbors
alpha <- 0.5  # hyperparameter (maps to 'sigma' in this package)
T     <- 10   # SNF iterations

# ----------------------------------------------------------
# 2) Load built-in example data: Data1, Data2
# ----------------------------------------------------------
data(Data1)
data(Data2)

# Verify data dimensions match
stopifnot(nrow(Data1) == nrow(Data2))
n <- nrow(Data1)

# Create Ground Truth labels (Cluster 1 vs Cluster 2)
truelabel <- c(rep(1, n/2), rep(2, n/2))

# ----------------------------------------------------------
# 3) Distance matrices (SQUARED Euclidean)
# ----------------------------------------------------------
# dist2 calculates squared Euclidean distance
Dist1 <- dist2(as.matrix(Data1), as.matrix(Data1))
Dist2 <- dist2(as.matrix(Data2), as.matrix(Data2))

# ----------------------------------------------------------
# 4) Build similarity graphs (affinity matrices)
# ----------------------------------------------------------
# FIX: Use 'sigma = alpha'. The package expects 'sigma'.
W1 <- affinityMatrix(Dist1, K = K, sigma = alpha)
W2 <- affinityMatrix(Dist2, K = K, sigma = alpha)

# Visualize initial clusters
displayClusters(W1, truelabel)
displayClusters(W2, truelabel)

# ----------------------------------------------------------
# 6) Fuse graphs with SNF
# ----------------------------------------------------------
W <- SNF(list(W1, W2), K = K, t = T)

# ----------------------------------------------------------
# 7) Cluster fused network (spectral clustering)
# ----------------------------------------------------------
C <- 2   # We know there are 2 true clusters
group <- spectralClustering(W, K = C)

displayClusters(W, group)

# ----------------------------------------------------------
# 8) Evaluate clustering quality (NMI vs truth)
# ----------------------------------------------------------
SNFNMI <- calNMI(group, truelabel)
cat(sprintf("\nNMI(group, truelabel) = %f\n", SNFNMI))

# ----------------------------------------------------------
# 9) Concordance between fused network and each individual
# ----------------------------------------------------------
# FIX: Added 'C' as the second argument.
ConcordanceMatrix <- concordanceNetworkNMI(list(W, W1, W2), C)
cat("\nConcordance Network NMI matrix:\n")
print(ConcordanceMatrix)

############################################################
# Multi-view example (Modified for missing Digits data)
############################################################
cat("\n--- Running Multi-view Example ---\n")

# FIX : The 'Digits' dataset is missing.
# We will create 'dataL' using Data1 and Data2 to let the code run.

dataL <- list(as.matrix(Data1), as.matrix(Data2))
label <- truelabel
C_multi <- 2  # We are using Data1/Data2, so C is 2, not 3.

# Parameters
K     <- 20
alpha <- 0.5
T     <- 20

# Calculate squared distances for the list
distL <- lapply(dataL, function(x) dist2(x, x))

# FIX: Use 'sigma = alpha'
affinityL <- lapply(distL, function(d) affinityMatrix(d, K = K, sigma = alpha))

# FIX: Passed C_multi (2) instead of 3
Concordance_matrix <- concordanceNetworkNMI(affinityL, C_multi)
cat("\nConcordance_matrix (fused vs each view):\n")
print(Concordance_matrix)

# Fuse the list
W_multi <- SNF(affinityL, K = K, t = T)

# Cluster the fused result
clustering <- spectralClustering(W_multi, K = C_multi)

# Calculate NMI
NMI_multi <- calNMI(clustering, label)
cat(sprintf("\nMulti-view NMI(clustering, label) = %f\n", NMI_multi))

############################################################
# Label propagation example: predict labels for new samples
############################################################
cat("\n--- Running Label Propagation ---\n")

set.seed(1)
                    
# Create train/test split
n_train <- floor(0.8 * length(label))
trainSample <- sample.int(length(label), n_train)

# Split the list data
train  <- lapply(dataL, function(x) x[trainSample, , drop = FALSE])
test   <- lapply(dataL, function(x) x[-trainSample, , drop = FALSE])
groups <- label[trainSample]

# Prediction Parameters
K      <- 20
alpha  <- 0.5
t      <- 20
method <- TRUE

# Note: groupPredict typically uses 'alpha' in most versions, 
# but if it fails, try changing 'alpha =' to 'sigma ='.
                 
newLabel <- groupPredict(train, test, groups, K = K, alpha = alpha, t = t, method = method)

# Calculate accuracy on the test set
# newLabel often includes training labels at the start, so we take the tail
                 
pred_test <- tail(newLabel, length(label) - n_train)
acc <- mean(label[-trainSample] == pred_test)

cat(sprintf("\nLabel propagation accuracy = %f\n", acc))

############################################################
# End
############################################################

head(Data1)

head(Data2)

# They simulate two different omics modalities measured on the same patients.

# Think of them as:

# Data1 → gene expression (RNA-seq–like)
# Data2 → DNA methylation (CpG-like)
# But this is only an analogy — the values are synthetic.

# | Object  | Meaning                                   |
# | ------- | ----------------------------------------- |
# | Rows    | patients / samples                        |
# | Columns | features (genes, methylation sites, etc.) |

# truelabel <- c(rep(1, n/2), rep(2, n/2))

# First half of patients → Cluster 1
# Second half → Cluster 2

# Crucially:
# Each modality alone only partially separates clusters
# SNF fuses them to recover the structure more clearly
# This is exactly what the SNF paper wanted to demonstrate.

# In real use:
# You might have 20,000 genes
# Or 500,000 CpGs
# Or hundreds of proteins

# SNF does not require the modalities to have the same number of features.
# Only sample-level similarities matter; feature count is irrelevant to SNF.

unique(truelabel)

head(Dist1)
head(Dist2) 

# W1 <- affinityMatrix(Dist1, K = K, sigma = alpha)
# W2 <- affinityMatrix(Dist2, K = K, sigma = alpha)

head(W1)
head(W2)

# Visualize initial clusters
# displayClusters(W1, truelabel)
# displayClusters(W2, truelabel)

# SNF

# The input to our algorithm can be feature vectors, pairwise distances, or pairwise similarities. 
# The learned status matrix can then be used for retrieval, clustering, and classification.

displayClustersWithHeatmap(W, group)

library(igraph)

threshold <- quantile(W[upper.tri(W)], 0.75)
W_thresh <- W
W_thresh[W_thresh < threshold] <- 0

g <- graph_from_adjacency_matrix(W_thresh, mode = "undirected", weighted = TRUE, diag = FALSE)

V(g)$color <- c("red", "blue", "green")[group]
E(g)$width <- E(g)$weight * 5

plot(g, vertex.size = 8, vertex.label = NA, layout = layout_with_fr(g), main = "SNF Network")

mds <- cmdscale(1 - W, k = 2)
plot(mds[,1], mds[,2], col = c("red", "blue")[group], pch = 19, main = "SNF Network")

par(mfrow = c(1, 3))
displayClusters(W1, group); title("Network 1")
displayClusters(W2, group); title("Network 2")
displayClusters(W, group); title("Fused")
par(mfrow = c(1, 1))

library(igraph)

# Set up 3 plots side by side
par(mfrow = c(1, 3))

# Network 1 Graph

threshold1 <- quantile(W1[upper.tri(W1)], 0.75)
W1_thresh <- W1
W1_thresh[W1_thresh < threshold1] <- 0
g1 <- graph_from_adjacency_matrix(W1_thresh, mode = "undirected", weighted = TRUE, diag = FALSE)
V(g1)$color <- c("red", "blue", "green")[group]
E(g1)$width <- E(g1)$weight * 5
plot(g1, vertex.size = 8, vertex.label = NA, layout = layout_with_fr(g1), main = "Network 1")

# Network 2 Graph

threshold2 <- quantile(W2[upper.tri(W2)], 0.75)
W2_thresh <- W2
W2_thresh[W2_thresh < threshold2] <- 0
g2 <- graph_from_adjacency_matrix(W2_thresh, mode = "undirected", weighted = TRUE, diag = FALSE)
V(g2)$color <- c("red", "blue", "green")[group]
E(g2)$width <- E(g2)$weight * 5
plot(g2, vertex.size = 8, vertex.label = NA, layout = layout_with_fr(g2), main = "Network 2")

# Fused Network Graph

threshold <- quantile(W[upper.tri(W)], 0.75)
W_thresh <- W
W_thresh[W_thresh < threshold] <- 0
g <- graph_from_adjacency_matrix(W_thresh, mode = "undirected", weighted = TRUE, diag = FALSE)
V(g)$color <- c("red", "blue", "green")[group]
E(g)$width <- E(g)$weight * 5
plot(g, vertex.size = 8, vertex.label = NA, layout = layout_with_fr(g), main = "Fused Network")

# Reset to single plot
par(mfrow = c(1, 1))

library(igraph)

# Create fused network graph
threshold <- quantile(W[upper.tri(W)], 0.75)
W_thresh <- W
W_thresh[W_thresh < threshold] <- 0

g <- graph_from_adjacency_matrix(W_thresh, mode = "undirected", weighted = TRUE, diag = FALSE)

# Node customization
V(g)$color <- c("red", "blue", "green", "orange", "purple")[group]
V(g)$label <- 1:vcount(g)           # Node IDs
V(g)$label.color <- "black"
V(g)$label.cex <- 0.8
V(g)$size <- 15
V(g)$frame.color <- "darkgray"      # Node border color

# Edge customization
E(g)$width <- E(g)$weight * 5       # Thicker = stronger connection
E(g)$color <- rgb(0.5, 0.5, 0.5, alpha = 0.5)  # Semi-transparent gray
E(g)$curved <- 0.1                  # Slight curve to edges

# Plot
plot(g, 
     layout = layout_with_fr(g),
     main = "SNF Fused Network")

# Add legend
legend("topright", 
       legend = paste("Cluster", sort(unique(group))),
       col = c("red", "blue", "green", "orange", "purple")[sort(unique(group))],
       pch = 19,
       pt.cex = 1.5,
       title = "Clusters",
       bty = "n")

# Nodes (circles)

# Each node = one sample/patient/data point
# Red nodes = Cluster 1 (samples 1-100)
# Blue nodes = Cluster 2 (samples 101-200)
# Numbers inside = sample IDs

# Edges (lines connecting nodes)

# Lines connect samples that are similar to each other
# Only showing connections above a threshold (strongest similarities)
# Thicker edges = stronger similarity between samples

library(igraph)

options(repr.plot.width=16, repr.plot.height=12)

threshold <- quantile(W[upper.tri(W)], 0.75)
W_thresh <- W
W_thresh[W_thresh < threshold] <- 0

g <- graph_from_adjacency_matrix(W_thresh, 
                                 mode = "undirected", 
                                 weighted = TRUE, 
                                 diag = FALSE)

V(g)$color <- c("red", "blue", "green")[group]
V(g)$label <- NA  # Remove labels to see edges better

E(g)$width <- E(g)$weight * 3
E(g)$color <- "gray30"  # Darker edges to see better

plot(g, 
     vertex.size = 3,              # SMALL nodes
     layout = layout_with_fr(g),
     main = "SNF Fused Network")

legend("topright", 
       legend = paste("Cluster", unique(group)),
       col = c("red", "blue"),
       pch = 19,
       pt.cex = 2)





cat("Simply another version of the code")

############################################################
# Compatible with CRAN SNFtool
############################################################

library(SNFtool)
set.seed(1)

############################################################
# Basic SNF example (Data1, Data2)
############################################################

# Parameters
K     <- 20    # number of neighbors
T     <- 10    # number of SNF iterations

# Load example data
data(Data1)
data(Data2)

stopifnot(nrow(Data1) == nrow(Data2))
n <- nrow(Data1)

# Ground-truth labels (simulation)
truelabel <- c(rep(1, n/2), rep(2, n/2))

# Normalize (recommended for continuous data)
Data1 <- standardNormalization(Data1)
Data2 <- standardNormalization(Data2)

# Distance matrices (SQUARED Euclidean)
Dist1 <- dist2(as.matrix(Data1), as.matrix(Data1))
Dist2 <- dist2(as.matrix(Data2), as.matrix(Data2))

# Affinity matrices (CRAN-safe: no alpha argument)
W1 <- affinityMatrix(Dist1, K)
W2 <- affinityMatrix(Dist2, K)

# Inspect individual networks
displayClusters(W1, truelabel)
displayClusters(W2, truelabel)

# Similarity Network Fusion
W <- SNF(list(W1, W2), K, T)

# Spectral clustering on fused network
C <- 2
group <- spectralClustering(W, C)

# Inspect fused clustering
displayClusters(W, group)

# Evaluate clustering
SNFNMI <- calNMI(group, truelabel)
cat("SNF NMI:", SNFNMI, "\n")

# Concordance between networks
ConcordanceMatrix <- concordanceNetworkNMI(list(W, W1, W2), C)
print(ConcordanceMatrix)

############################################################
# Multi-view example (Digits)
############################################################

# Load Digits example (provides dataL and label)
data(Digits)

# Parameters
K <- 20
T <- 20

# Normalize each view (optional but recommended)
dataL <- lapply(dataL, standardNormalization)

# Distance matrices
distL <- lapply(dataL, function(x) dist2(x, x))

# Affinity matrices (CRAN-safe)
affinityL <- lapply(distL, function(d) affinityMatrix(d, K))

# Concordance across views
Concordance_matrix <- concordanceNetworkNMI(affinityL, 3)
print(Concordance_matrix)

############################################################
# Subtyping using SNF
############################################################

W_digits <- SNF(affinityL, K, T)
clustering <- spectralClustering(W_digits, 3)

NMI_digits <- calNMI(clustering, label)
cat("Digits NMI:", NMI_digits, "\n")

############################################################
# Label propagation (groupPredict)
############################################################

# Train / test split
n_train <- floor(0.8 * length(label))
trainSample <- sample.int(length(label), n_train)

train  <- lapply(dataL, function(x) x[trainSample, , drop = FALSE])
test   <- lapply(dataL, function(x) x[-trainSample, , drop = FALSE])
groups <- label[trainSample]

# Parameters for prediction
K      <- 20
alpha  <- 0.5
t      <- 20
method <- TRUE

# Predict labels
newLabel <- groupPredict(train, test, groups, K, alpha, t, method)

# Accuracy
pred_test <- tail(newLabel, length(label) - n_train)
accuracy <- mean(label[-trainSample] == pred_test)
cat("Label propagation accuracy:", accuracy, "\n")

############################################################
# End
############################################################






cat("

NMI = Normalized Mutual Information

It's a metric that measures how well your clustering matches the true labels (ground truth).

What NMI tells you:

NMI = 1.0 → Perfect match! Your clusters exactly match the true labels
NMI = 0.0 → No relationship between clusters and true labels (random)
NMI = 0.5-0.8 → Good clustering, captures most of the structure
NMI > 0.8 → Excellent clustering

If you get:

NMI = 0.95  → Your SNF fusion worked great! Almost perfect clustering
NMI = 0.75  → Good! Captured most of the cluster structure
NMI = 0.40  → Poor clustering, missed a lot of the structure
NMI = 0.10  → Failed, clusters don't match reality

")



cat("

ARI = Adjusted Rand Index
It's another metric that measures how well your clustering matches the true labels, similar to NMI but with different properties.
What ARI tells you:

ARI = 1.0 → Perfect match! Clusters exactly match true labels
ARI = 0.0 → Random clustering (no better than chance)
ARI < 0 → Worse than random (yes, it can be negative!)
ARI = 0.5-0.8 → Good clustering
ARI > 0.8 → Excellent clustering

")



# Install if needed
# install.packages("mclust")

library(mclust)

# Calculate ARI
ari_score <- adjustedRandIndex(group, truelabel)
cat("ARI =", ari_score, "\n")



cat("

| Property                       | **NMI (Normalized Mutual Information)** | **ARI (Adjusted Rand Index)**     |
| ------------------------------ | --------------------------------------- | --------------------------------- |
| Value range                    | **0 to 1**                              | **−1 to 1**                       |
| Interpretation                 | Information-theoretic similarity        | Agreement-based similarity        |
| Handles unequal cluster sizes  | **Yes**                                 | **Yes**                           |
| Baseline for random clustering | **0**                                   | **0**                             |
| Can be negative?               | **No**                                  | **Yes**                           |
| Popular in                     | Machine learning, clustering            | Statistics, clustering validation |

")




