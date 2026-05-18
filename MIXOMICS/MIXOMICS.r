library(mixOmics)

data(breast.TCGA)

cat("

The data were divided into a training set with a subset of 150 samples from the mRNA, miRNA and proteomics data, 
and a test set includes 70 samples, but only from the mRNA, miRNA and methylation data (proteomics missing).

The mixOmics TCGA dataset is accessed via `breast.TCGA’ and contains the following:

breast.TCGA$data.train$mirna (continuous matrix): 150 rows and 184 columns. The expression levels of 184 different sections of miRNA.

breast.TCGA$data.train$mrna (continuous matrix): 150 rows and 200 columns. The expression levels of 200 different sections of mRNA.

breast.TCGA$data.train$protein (continuous matrix): 150 rows and 142 columns. The abundance of 142 different proteins

breast.TCGA$data.train$subtype (categorical vector): length of 150. 

Indicates the breast cancer subtype of each subject. Includes Basal, Her2 and LumA.

")

cat("

What is a Multi-Block Model?

A multi-block model integrates multiple datasets (blocks) simultaneously rather than analyzing them separately or in pairs.

Key Differences:

Pairwise approach (lines 26-31 in your script):

Analyzes miRNA vs mRNA
Then miRNA vs proteomics
Then mRNA vs proteomics

3 separate models, each ignoring the third dataset

Multi-block approach (DIABLO, line 86):

Analyzes miRNA + mRNA + proteomics all at once
Finds components that maximize:

Correlation between all datasets (as specified in design matrix)
Discrimination between cancer subtypes

An integrated model that considers all relationships simultaneously :

Why is this better?

Imagine you're trying to understand a disease:

Pairwise: You look at genes, then proteins, then metabolites separately
Multi-block: You see how genes, proteins, AND metabolites work together as a system

Multi-block captures synergistic effects that pairwise analysis misses. For example:

A gene change might only affect cancer subtype when combined with specific protein and miRNA changes
Multi-block finds these coordinated multi-omics signatures

")



head(breast.TCGA$data.train$mirna, 3)
dim(breast.TCGA$data.train$mirna)

head(breast.TCGA$data.train$mrna, 3)
dim(breast.TCGA$data.train$mrna)

head(breast.TCGA$data.train$mirna, 3)
dim(breast.TCGA$data.train$mirna)

head(breast.TCGA$data.train$protein, 3)
dim(breast.TCGA$data.train$protein)

length(breast.TCGA$data.train$subtype)
unique(breast.TCGA$data.train$subtype)

head(breast.TCGA$data.train$subtype)
tail(breast.TCGA$data.train$subtype)

# set a list of all the X dataframes
data = list(miRNA = breast.TCGA$data.train$mirna, 
            mRNA = breast.TCGA$data.train$mrna,
            proteomics = breast.TCGA$data.train$protein)

# head(data)
str(data)

lapply(data, dim) # check their dimensions

Y = breast.TCGA$data.train$subtype # set the response variable as the Y df
summary(Y)

cat("Initial Analysis")

cat("Pairwise PLS Comparisons")

# Specifies to keep 25 features from each dataset for components 1 and 2. This is arbitrary at this stage

list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(25, 25)

print(list.keepX)
print(list.keepY)

# generate three pairwise PLS models
pls1 <- spls(data[["miRNA"]], data[["mRNA"]], 
             keepX = list.keepX, keepY = list.keepY) 

# Creates a sparse PLS model between miRNA and mRNA. Sparse PLS performs feature selection while modeling relationships.

pls2 <- spls(data[["miRNA"]], data[["proteomics"]], 
             keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(data[["mRNA"]], data[["proteomics"]], 
             keepX = list.keepX, keepY = list.keepY)

# Creates two more pairwise PLS models. This helps understand bivariate relationships before the full multi-block model.

pls1

pls2

pls3

# Circle Correlation Plots for pairwise PLS models on the breast TCGA data. 
# Only displays the top 25 features for each dimension, subsetting by those with a correlation above 0.5.

# plot features of first PLS
plotVar(pls1, cutoff = 0.5, title = "(a) miRNA vs mRNA", 
        legend = c("miRNA", "mRNA"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# plot features of second PLS
plotVar(pls2, cutoff = 0.5, title = "(b) miRNA vs proteomics", 
        legend = c("miRNA", "proteomics"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# plot features of third PLS
plotVar(pls3, cutoff = 0.5, title = "(c) mRNA vs proteomics", 
        legend = c("mRNA", "proteomics"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# Below is are correlations between the first component of each dataset for all three PLS models. 
# It can be seen the values are around ~0.8 - 0.9.

# calculate correlation of miRNA and mRNA
cor(pls1$variates$X, pls1$variates$Y) 

# calculate correlation of miRNA and proteins
cor(pls2$variates$X, pls2$variates$Y) 

# calculate correlation of mRNA and proteins
cor(pls3$variates$X, pls3$variates$Y) 

cat("Initial Model")

# for square matrix filled with 0.1s
design = matrix(0.1, 
                ncol = length(data), 
                nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s

design

# Creates a design matrix that specifies how strongly different datasets should be correlated:
# 0.1 = weak connection between all pairs
# Diagonal = 0 (dataset doesn't connect to itself)
# This is a good default starting point.

# An arbitrarily high number of components (ncomp = 5) will be used.

basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 5, design = design) 

# Fits initial DIABLO model with 5 components (arbitrarily high to determine optimal number later). 
# block.splsda = multi-block sparse Partial Least Squares Discriminant Analysis.

cat("Tuning the number of components")

# run component number tuning with repeated CV
perf.diablo = perf(basic.diablo.model, 
                   validation = 'Mfold', 
                   folds = 10, nrepeat = 10) 

plot(perf.diablo) # plot output of tuning

# Performs 10-fold cross-validation repeated 10 times to assess model performance with different numbers of components. 
# This is a rigorous evaluation approach.

# Choosing the number of components in block.plsda using perf() with 10 × 10-fold CV function in the breast.TCGA study. 
# Classification error rates (overall and balanced, see Section 7.3) are represented on the y-axis with respect 
# to the number of components on the x-axis for each prediction distance presented in PLS-DA

# Both  overall and balanced error rate (BER) decrease from 1 to 2 components. 

# set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
# show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote 

# Extracts the optimal number of components based on:

# Overall.BER = Balanced Error Rate (accounts for class imbalance)
# centroids.dist = distance metric used for classification
# WeightedVote = voting scheme for predictions

cat("Tuning the number of features")

# set grid of values for each component to test
test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)), 
                   miRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                   proteomics = c(5:9, seq(10, 18, 2), seq(20,30,5)))

# Creates a grid of candidate values for the number of features to test:
# 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30
# This is a reasonable range for biomarker discovery.

# run the feature selection tuning
tune.TCGA = tune.block.splsda(X = data, 
                              Y = Y, 
                              ncomp = ncomp, 
                              test.keepX = test.keepX, 
                              design = design,
                              validation = 'Mfold', 
                              folds = 10, nrepeat = 1,
                              dist = "centroids.dist")

# This is computationally expensive - tests all combinations of feature numbers for each dataset to find optimal signature. 
# Note: only 1 repeat here (vs 10 earlier) to save time.

list.keepX = tune.TCGA$choice.keepX # set the optimal values of features to retain
list.keepX

cat("Final model")

# set the optimised DIABLO model
final.diablo.model = block.splsda(X = data, 
                                  Y = Y, 
                                  ncomp = ncomp, 
                                  keepX = list.keepX, 
                                  design = design)

# Builds the final model with optimized parameters. This is the model you'll use for interpretation and prediction.

final.diablo.model$design 

# The selected variables can be extracted with the function selectVar(), for example in the mRNA block, 
# as seen below. Note that the stability of selected variables can be extracted from the output of the perf() function.

# the features selected to form the first component
selectVar(final.diablo.model, block = 'mRNA', comp = 1)$mRNA$name 

# Extracts names of selected features from mRNA dataset for component 1. These are your biomarker candidates.

cat("Plots")

# plotDIABLO() is a diagnostic plot to check whether the correlation between components from each data set 
# has been maximised as specified in the design matrix. 
# We specify which dimension to be assessed with the ncomp argument.

plotDiablo(final.diablo.model, ncomp = 1)

# Diagnostic plot showing how well components from each dataset correlate. Also shows sample separation by subtype with confidence ellipses.

# Diagnostic plot from multiblock sPLS-DA applied on the breast.TCGA study. 
# Samples are represented based on the specified component (here ncomp = 1) 
# for each data set (mRNA, miRNA and protein). 
# Samples are coloured by breast cancer subtype and 95% confidence ellipse plots are represented.

# The first components from each data set are highly correlated to each other (indicated by the large numbers in the bottom left). 
# The colours and ellipses related to the sample subtypes indicate the discriminative power of each component to separate 
# the different tumour subtypes. 

# For the first component, the centroids of each subtype are distinct 
# but there is a moderate amount of overlap between each sample group in their confidence ellipses.

# plotIndiv() function projects each sample into the space spanned by the components of each block.

plotIndiv(final.diablo.model, 
          ind.names = FALSE, 
          legend = TRUE, 
          title = 'DIABLO Sample Plots')


# Sample plot - shows where each sample projects in the component space for each omics layer. 
# Helps visualize sample clustering by subtype.

# Sample plot from multiblock sPLS-DA performed on the breast.TCGA study. 
# The samples are plotted according to their scores on the first 2 components for each data set. 
# Samples are coloured by cancer subtype

# Arrow plot - shows agreement between datasets at the sample level:

# Arrow start = centroid across all datasets for a sample
# Arrow tips = sample position in each individual dataset
# Short arrows = good agreement between omics layers



plotArrow(final.diablo.model, 
          ind.names = FALSE, 
          legend = TRUE, 
          title = 'mixomics')

# The start of the arrow indicates the centroid between all data sets for a given sample and the tips of the arrows 
# indicate the location of that sample in each block. 
# Such graphics highlight the agreement between all data sets at the sample level. 

plotArrow(final.diablo.model, 
          ind.names = FALSE, 
          legend = TRUE, 
          title = 'DIABLO')

cat("Variable plots")

# A majority of the miRNA variables are positively correlated with the first component while the mRNA variables seem 
# to separate along this dimension. 
# These first two components correlate highly with the selected variables from the proteomics dataset. 
# From this, the correlation of each selected feature from all three datasets can be evaluated based on their proximity.

plotVar(final.diablo.model, var.names = FALSE, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 17, 15), cex = c(2,2,2), 
        col = c('darkorchid', 'brown1', 'lightgreen'))

# Shows selected features from all datasets in correlation circle.
# Features close together are highly correlated across omics layers.

circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks= c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

# Circos plot - elegant visualization showing correlations > 0.7 between features from different datasets. 
# Each arc connects correlated features.

# Another visualisation of the correlations between the different types of variables is the relevance network, 
# which is also built on the similarity matrix (as is the circos plot). Each colour represents a type of variable.

# Same optimally selected features, but only connections where correlation > 0.7 are drawn.

# Example: A line between miRNA "hsa-mir-200c" and mRNA "ZEB1" means:

# These are both in the final signature
# They're correlated > 0.7 across samples
# Biological interpretation: miR-200c might regulate ZEB1 (known tumor suppressor relationship)

network(final.diablo.model, 
        blocks = c(1,2,3),
        color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)

# Creates network graph where:

# Nodes = features (colored by omics type)
# Edges = correlations > 0.4
# Shows multi-omics regulatory networks.

library(igraph)
my.network = network(final.diablo.model, 
                     blocks = c(1,2,3),
                     color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.4)
write.graph(my.network$gR, file = "myNetwork.gml", format = "gml")

# Features shown: Same selected features, but with lower cutoff (0.4 vs 0.7), so more connections.
# Nodes = individual features (miRNAs, mRNAs, proteins)
# Edges = correlations > 0.4

plotLoadings(final.diablo.model, comp = 1, contrib = 'max', method = 'median')

plotLoadings(final.diablo.model, comp = 2, contrib = 'max', method = 'median')

plotLoadings(final.diablo.model, comp = 3, contrib = 'max', method = 'median')

# Key Insight:

# The number of features decreases through the analysis:

# Start: All features in original data
# miRNA: ~184 features
# mRNA: ~200 features
# Proteomics: ~142 features

# Pairwise PLS: 25 features per dataset (arbitrary)
# After tuning (line 125): Optimal number (might be 8 mRNA, 12 miRNA, 10 proteins)
# Final plots: Only these optimized features appear

# The cimDIABLO() function is a clustered image map specifically implemented 
# to represent the multi-omics molecular signature expression for each sample

cimDiablo(final.diablo.model)

# Clustered heatmap showing expression of selected features across all samples. Helps identify multi-omics signatures for each subtype.

cat("Performance of the model")

# We assess the performance of the model using 10-fold cross-validation repeated 10 times, using the function perf(). 
# The method runs a block.splsda() model on the pre-specified arguments input from the final.diablo.model object, 
# but on cross-validated samples. We then assess the accuracy of the prediction on the left-out samples.

# Since the tune() function was used with the centroid.dist argument, 
# the outputs of the perf() function for that same distance are examined. 

# run repeated CV performance evaluation
perf.diablo = perf(final.diablo.model, 
                   validation = 'Mfold', 
                   M = 10, nrepeat = 10, 
                   dist = 'centroids.dist') 

perf.diablo$MajorityVote.error.rate

# Final rigorous evaluation - 10×10 CV on the optimized model to get reliable performance estimates.

perf.diablo$WeightedVote.error.rate

# Shows classification error rates using two voting schemes:
# MajorityVote - each dataset gets equal vote
# WeightedVote - datasets weighted by performance

auc.splsda = auroc(final.diablo.model, 
                   roc.block = "miRNA", 
                   roc.comp = 2, 
                   print = FALSE)

# Computes AUC-ROC curves for the miRNA block. Good for assessing discriminative power for multi-class problems.

cat("Prediction on an external test set")

# The predict() function predicts the class of samples from a test set. 
# In our specific case, one data set is missing in the test set but the method can still be applied. 
# Make sure the name of the blocks correspond exactly.

data.test.TCGA = list(mRNA = breast.TCGA$data.test$mrna,
                      miRNA = breast.TCGA$data.test$mirna)

predict.diablo = predict(final.diablo.model, newdata = data.test.TCGA)

confusion.mat = get.confusion_matrix(truth = breast.TCGA$data.test$subtype,
                                     predicted = predict.diablo$WeightedVote$centroids.dist[,2])
confusion.mat

get.BER(confusion.mat)



cat("

Do We Want High Correlations Between Omics Layers?

It depends on your goal! This is nuanced:

High correlations (~0.8-0.9) are GOOD when:

Biological coherence: Different omics layers measuring the same biological process should be correlated

Example: mRNA → protein translation

If mRNA levels are high, protein levels should be high too

High correlation = the biology makes sense

Data quality: High correlations suggest:

Your measurements are reliable
Technical variation is low
The datasets are measuring real biological signal

Integration success: The model successfully found:

Common patterns across omics layers
Coordinated molecular changes
Multi-omics signatures

But you DON'T want TOO high correlation because:

Redundancy: If correlation = 1.0, the datasets provide identical information

You're not gaining anything from multi-omics
Just use one dataset

Loss of complementarity: Each omics layer should provide unique information

miRNA regulates genes (regulatory info)
mRNA shows transcription (gene expression)
Proteins show final functional molecules (functional info)
You want ~0.7-0.9 correlation (related but not identical)

")

# design = matrix(0.1, ncol = length(data), nrow = length(data))
# The 0.1 design matrix is actually quite weak. This tells DIABLO:

# "Find components that are somewhat correlated between datasets"
# "But prioritize discrimination of cancer subtypes"

# You could increase this to 0.5 or 0.7 if you want stronger integration, but 0.1 is a good starting point.
# The observed correlations of 0.8-0.9 from the pairwise PLS models (lines 67-73) are actually results, not input. 

# This means:

# Even with weak constraints (0.1), the data naturally shows strong relationships
# This is good - it validates that your omics layers are biologically connected



cat("What Makes Them Biomarkers? 

1. They Discriminate Between Disease States

The DIABLO model selects features that best separate the different breast cancer subtypes. 

Look at what the model does:

final.diablo.model = block.splsda(X = data, Y = Y, ...)

where Y = breast cancer subtype (Basal, Her2, LumA, LumB, Normal)

The model asks: 'Which miRNAs, mRNAs, and proteins differ most between these cancer subtypes?'

Features are selected because:

Basal subtype patients have high expression of certain features
Her2 subtype patients have different expression patterns
LumA/LumB patients have yet different patterns

2. They're Statistically Validated Through Cross-Validation

tune.TCGA = tune.block.splsda(..., validation = 'Mfold', folds = 10, nrepeat = 1)

Features are only kept if they consistently predict cancer subtype across different data splits

This means:

They're not just noise or random fluctuations
They're reproducible predictors
They generalize to new samples

3. They're Part of a Multi-Omics Signature

Unlike single-omics analysis, DIABLO selects features that are:

Biologically coherent across multiple molecular levels:

circosPlot(final.diablo.model, cutoff = 0.7, ...)

If you see connections between:

miRNA 'hsa-mir-21' (purple)
mRNA 'PTEN' (brown)
Protein 'PTEN' (green)

This suggests a regulatory cascade:

miR-21 regulates PTEN mRNA

PTEN mRNA produces PTEN protein

This entire pathway differs between cancer subtypes

")

# This is MORE meaningful than individual markers because:

# It captures biological mechanisms (not just associations)
# It's more robust (multiple correlated measurements)
# It's more likely to be causal (not just correlation)


