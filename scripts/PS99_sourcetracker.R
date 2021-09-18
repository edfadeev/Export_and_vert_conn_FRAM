library("phyloseq")


source('scripts/SourceTracker_nrestart_ndraws_100.R')


PS99_merged<- readRDS("./PS99_phy0_prev.rds")

PS99_phy0_prev.water<-subset_samples(PS99_phy0.prev, Type!="Sediment")
PS99_phy0_prev.water<- prune_taxa(taxa_sums(PS99_phy0_prev.water)>0,PS99_phy0_prev.water)

#ASV table
otus <- as.data.frame(otu_table(PS99_phy0_prev.water, taxa_are_rows = TRUE))
otus <- t(as.matrix(otus))

#metadata
metadata <- sample_data(PS99_phy0.prev)

# extract only those samples in common between the two tables
common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids,]
metadata <- metadata[common.sample.ids,]
# double-check that the mapping file and otu table
# had overlapping samples
if(length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

# extract the source environments and source/sink indices
train.ix <- which(metadata$SourceSink.depth=='source')
test.ix <- which(metadata$SourceSink.depth=='sink')
envs <- metadata$Source.Alt
if(is.element('StationName',colnames(metadata))) desc <- metadata$StationName
if(is.element('StationName',colnames(metadata))) depth <- metadata$Type

#####################################
#Run SourceTracker
#####################################
# tune the alpha values using cross-validation (this is slow!)
#tune.results <- tune.st(otus[train.ix,], envs[train.ix])
#alpha1 <- tune.results$best.alpha1
#alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix], rarefaction_depth=5000)

# Estimate source proportions in test data
ST_results <- predict.sourcetracker(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2, rarefaction_depth=5000, full.results=TRUE)

# Estimate leave-one-out source proportions in training data 
ST_results.train <- predict.sourcetracker(st, alpha1=alpha1, alpha2=alpha2, rarefaction_depth=5000, full.results=TRUE)



# train SourceTracker object on training data
st.PA <- sourcetracker(otus[test.ix,], envs[test.ix], rarefaction_depth=5000)

# Estimate source proportions in test data
ST_results.PA <- predict.sourcetracker(st.PA,otus[train.ix,], alpha1=alpha1, alpha2=alpha2, rarefaction_depth=5000, full.results=TRUE)

# Estimate leave-one-out source proportions in training data 
ST_results.train.PA <- predict.sourcetracker(st.PA, alpha1=alpha1, alpha2=alpha2, rarefaction_depth=5000, full.results=TRUE)

