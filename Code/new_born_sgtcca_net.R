library(kernlab)
library(cluster)
library(SmCCNet)
library(igraph)
library(RCy3)
# read in all the data
feature_keys <- read.csv('multiomics_feature_key.csv')
feature_table <- read.csv('feature_table.csv')
metabolome_table <- read.csv('metabolome_key.csv')
pheno_data <- read.csv('clinical_targets.csv')
feature_table$ORIG_ID <- NULL
feature_table$GA_AT_SAMPLE <- NULL

# source function for SGTCCA-Net
source('SGTCCA-Source.R')

# separate dataset into multi-view data
lipidome <- feature_table[, which(feature_keys$Dataset == 'Lipidome')]
metabolome <- feature_table[, 1196 + which(feature_keys$Dataset == 'Metabolome')]
proteome <- feature_table[, which(feature_keys$Dataset == 'Proteome') - 4329]



# lipidomics and metabolomics data should be logged
metabolome <- log(metabolome)
lipidome <- log(lipidome)

# define phenotype 
pheno <- as.matrix(pheno_data$SAMPLE_TO_BIRTH)
# center and scale each data
lipidome <- scale(lipidome)
metabolome <- scale(metabolome)
proteome <- scale(proteome)
pheno <- scale(pheno)

# run SGTCCA-Net and obtain canonical weight
data <- list(lipidome,metabolome,proteome, pheno)
# define the type of data for each profile
dtype <- c('lipid', 'metabolite', 'protein', 'pheno')
# run SGTCCA
result <- SGTCCA_Net(data_list = data, prob_subsamp_prop = 0.1, 
          correlation_list = list(c(1,2,3,4), c(1,2,4), c(1,3,4), c(2,3,4), 
                                  c(1,4), c(2,4), c(3,4)), 
          network_exclude_data = 4, 
          data_type = dtype, num_workers = 5,
          pheno = pheno, min_mod_size = 30, max_mod_size = 300,num_iterations = 10, 
          saving_dir = getwd(), 
          common_subsamp_prop = 0.08, distinct_subsamp_prop = 0.02)
# load SGTCCA canonical weight
load('gtcca_networkdata.Rdata')

# get network module
num_lipid <- sum(cc_weight_network[types == 'lipid'] != 0)
num_metabolite <- sum(cc_weight_network[types == 'metabolite'] != 0)
num_protein <- sum(cc_weight_network[types == 'protein'] != 0)

# scale canonical weight so that each molecular feature has the same scale of measurement acros different molecular profiles
cc_weight_network[types == 'lipid'] <- (cc_weight_network[types == 'lipid']/sum(cc_weight_network[types == 'lipid']))
cc_weight_network[types == 'metabolite'] <- (cc_weight_network[types == 'metabolite']/sum(cc_weight_network[types == 'metabolite']))*(num_metabolite/num_lipid)
cc_weight_network[types == 'protein'] <- (cc_weight_network[types == 'protein']/sum(cc_weight_network[types == 'protein']))*(num_protein/num_lipid)

# subset
row.names(cc_weight_network) <- colnames(combined_data)
# examine the distribution of canonical weight
hist(cc_weight_network[cc_weight_network!=0], breaks = 20)
# filter out canonical weight with weaker values
cc_weight_network[cc_weight_network < 0.006] <- 0

# re-scale canonical weight
sum(cc_weight_network[types == 'lipid']!=0)
sum(cc_weight_network[types == 'protein']!=0)
sum(cc_weight_network[types == 'metabolite']!=0)



# subset with only selected samples 
combined_sub <- combined_data[,which(cc_weight_network[,1] != 0)]
type_sub <- types[which(cc_weight_network[,1] != 0)]
lipid_sub <- combined_sub[, which(type_sub == 'lipid')]
metabolite_sub <- combined_sub[, which(type_sub == 'metabolite')]
protein_sub <- combined_sub[, which(type_sub == 'protein')]

# feature keys
feature_keys_lipid <- feature_keys[feature_keys$Feature %in% colnames(lipid_sub),]
feature_keys_protein <- feature_keys[feature_keys$Feature %in% colnames(protein_sub),]
colnames(metabolite_sub) <- sub('X', '', colnames(metabolite_sub))
feature_keys_metabolite <- feature_keys[match(colnames(metabolite_sub), feature_keys$Feature),]
match(feature_keys_metabolite$Feature, colnames(metabolite_sub))
# change unknown metabolite names
feature_keys_metabolite$Name[which(feature_keys_metabolite$Name == 'Unidentified metabolite')] <- feature_keys_metabolite$Feature[which(feature_keys_metabolite$Name == 'Unidentified metabolite')]
colnames(metabolite_sub) <- feature_keys_metabolite$Name

# change feature keys for protein
feature_keys_protein <- feature_keys[match(colnames(protein_sub), feature_keys$Feature),]
match(feature_keys_protein$Feature, colnames(protein_sub))
colnames(protein_sub) <- feature_keys_protein$Name

# change feature keys for lipids
feature_keys_lipid <- feature_keys[match(colnames(lipid_sub), feature_keys$Feature),]
colnames(lipid_sub) <- feature_keys_lipid$Name


# output all lipids that has non-zero canonical weight
write.csv(feature_keys_lipid, 'feature_keys_lipid_global_network.csv', row.names = FALSE)
# output all metabolites that has non-zero canonical weight
metabolite_information <- metabolome_table[metabolome_table$Analyte.Name %in% colnames(metabolite_sub),]
metabolite_information <- metabolite_information[!is.na(metabolite_information$HMDB.ID),]
write.csv(metabolite_information, 'feature_keys_metabolite_global_network.csv', row.names = FALSE)
# output all lipids that has non-zero canonical weight
write.csv(feature_keys_protein, 'feature_keys_protein_global_network.csv', row.names = FALSE)


# concatenate all molecular profiles
colnames(combined_sub) <- c(colnames(lipid_sub), colnames(metabolite_sub), colnames(protein_sub))
# construct higher-order adjacency matrix
abar[[1]] <- cov_comp(list(t(combined_sub), t(combined_sub), t(pheno)), correlation_list = list(c(1,2,3)))[[1]][,,1]
# scale adjacency matrix
abar[[1]] <- abar[[1]]/max(abar[[1]])
Abar <- abar[[1]]
# examine the distribution of edge values
hist(as.vector(Abar), 50)
# edge filtering
Abar[Abar < 0.05] <- 0
# subset the Pearson's correlation matrix
corsub <- cor_mat[which(cc_weight_network != 0), which(cc_weight_network != 0)]
# define data used for clustering and pruning
X_big <- combined_sub
data_type <- type_sub

# Create a function to compute the eigengap heuristic
eigengap_heuristic <- function(L, max_clusters = 10) {
  eigenvalues <- sort(abs(eigen(L, only.values = TRUE)$values), decreasing = TRUE)
  # If there are fewer than max_clusters eigenvalues, use the number available
  num_eigenvalues_to_consider <- min(length(eigenvalues), max_clusters)
  # Calculate gaps only within the first 'max_clusters' eigenvalues
  eigen_gaps <- diff(eigenvalues[1:num_eigenvalues_to_consider])
  # Find the position of the maximum gap
  max_gap_position <- which.max(eigen_gaps)
  # The optimal number of clusters k is the position of the max gap plus 1
  # But we limit it to 'max_clusters'
  optimal_k <- min(max_gap_position + 1, max_clusters)
  return(optimal_k)
}

# Assuming Abar is your adjacency matrix
graph <- graph_from_adjacency_matrix(1-Abar, mode = "undirected", diag = FALSE)
L <- laplacian_matrix(graph)
optimal_k <- eigengap_heuristic(L)
# Perform spectral clustering with the optimal number of clusters
cl <- specc(as.kernelMatrix(Abar), centers = optimal_k)
clusters <- cl@.Data
summary(as.factor(clusters))


for(i in 1:optimal_k)
{
  # define subnetwork
  network_module <- which(cl@.Data == i)
  if (length(network_module) < 30)
    next
  abar_sub <- Abar[network_module, network_module]
  cor_sub <- corsub[network_module,network_module]
  # prune network module
  networkPruning(Abar = abar_sub,CorrMatrix = cor_sub, 
                 type = data_type[network_module], 
                 data = X_big[,network_module],      
                 Pheno = pheno, ModuleIdx = i, min_mod_size = 30, 
                 max_mod_size = 300, method = 'NetSHy', 
                 saving_dir = getwd())
}



################################################# Network visualization: for top 15 molecular faetures only
# load network data
load('size_208_net_1.Rdata')
load('size_60_net_2.Rdata')
load('size_148_net_3.Rdata')
load('size_233_net_4.Rdata')
# make row names unique
row.names(M) <- colnames(M) <- make.unique(row.names(M))
# make adjacency matrix object
net_ppr <- igraph::graph_from_adjacency_matrix(M, weighted = TRUE,
                                               diag = FALSE, mode = "undirected")
# PageRank algorithm: all parameter setups are based on the recommendation
ranking <- igraph::page_rank(net_ppr, directed = FALSE, damping = 0.9, 
                             options = list(niter = 10^5,eps = 1e-06))$vector
# select top 15 molecular features based on PageRank score 
top_15 <- c(names(sort(ranking[sub_type == 'lipid'], decreasing = TRUE))[1:5], names(sort(ranking[sub_type == 'metabolite'], decreasing = TRUE))[1:5], names(sort(ranking[sub_type == 'protein'], decreasing = TRUE))[1:5])
# (only for network module 2)
top_15 <- c(names(sort(ranking, decreasing = TRUE))[1:15])
# subset the adjacency matrix
M_sub <- M[row.names(M) %in% top_15,colnames(M) %in% top_15]
# subset the correlation matrix
cor_sub <- correlation_sub[row.names(M) %in% top_15,colnames(M) %in% top_15]
# subset the type data
type_sub <- sub_type[row.names(M) %in% top_15]
# subset the PageRank score
rank_sub <- ranking[row.names(M) %in% top_15]
# set diagonal to 0
diag(M_sub) <- 0
# create igraph object
graph <- igraph::graph_from_adjacency_matrix(M_sub, mode = 'undirected', weighted = TRUE,
                                             diag = TRUE, add.colnames = NULL, add.rownames = NA)
# need to open cytoscape before visualization
V(graph)$type <- type_sub
V(graph)$type
# create network through Cytoscape
createNetworkFromIgraph(graph,"ptb_network_1")
createNetworkFromIgraph(graph,"ptb_network_2")
createNetworkFromIgraph(graph,"ptb_network_3")
createNetworkFromIgraph(graph,"ptb_network_4")





