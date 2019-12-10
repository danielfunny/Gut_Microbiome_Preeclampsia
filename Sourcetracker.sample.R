packages <- c("decontam", "ggplot2", "reshape2", "phyloseq","gridExtra")
sapply(packages, require, character.only = TRUE)
samdf <- read.table("pla_mata0128.metadata.txt", header=TRUE, stringsAsFactors = FALSE,row.names = 1,sep = "\t",comment.char = "")
otu <- read.table("pla.biom.txt",header = T,row.names = 1,sep = "\t")
samdf <-samdf[match(rownames(samdf), colnames(otu)),]

table(samdf$Bodysite)
neg.sites <- c("HO", "SWAB")
plac.sites <- c("placenta")
samdf$neg <- samdf$Bodysite %in% neg.sites
samdf$plac <- samdf$Bodysite %in% plac.sites
samdf$Type <- samdf$Bodysite
samdf$Type[samdf$neg] <- "Negative"
samdf$Type[samdf$plac] <- "Placenta"
st <- t(otu[,match(rownames(samdf), colnames(otu))])
samdf$Reads <- rowSums(st)
tax <- matrix(unlist(strsplit(as.character(otu$taxonomy),";")),ncol = 7,byrow = T)
rownames(tax) <- rownames(otu)
colnames(tax) <- c("kingdom","phylum","class","order","family","genus","species")
ps <- phyloseq(otu_table(st,taxa_are_rows = FALSE),sample_data(samdf),tax_table(tax))


metadata <- sample_data(ps)
metadata$Env <- metadata$Species
metadata$Description <- metadata$ID
metadata$SourceSink <- ifelse(metadata$neg=="TRUE","source","sink")

otus <- otu_table(ps)
otus <- data.frame(otus)

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
train.ix <- which(metadata$SourceSink=='source')
test.ix <- which(metadata$SourceSink=='sink')
envs <- metadata$Env
if(is.element('Description',colnames(metadata))) desc <- metadata$Description


# load SourceTracker package
source('SourceTracker.r')

# tune the alpha values using cross-validation (this is slow!)
#tune.results <- tune.st(otus[train.ix,], envs[train.ix])
# alpha1 <- tune.results$best.alpha1
# alpha2 <- tune.results$best.alpha2
# note: to skip tuning, run this instead:
alpha1 <- alpha2 <- 0.001

# train SourceTracker object on training data
st <- sourcetracker(otus[train.ix,], envs[train.ix],rarefaction_depth = 1000)

# Estimate source proportions in test data
results <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2,full.results = TRUE,rarefaction_depth = 1000)

# Estimate leave-one-out source proportions in training data 
results.train <- predict(st, alpha1=alpha1, alpha2=alpha2,full.results = TRUE,rarefaction_depth = 1000)

# plot results
#labels <- sprintf('%s %s', envs,desc)
#plot(results, labels[test.ix], type='pie')

# other plotting functions
# plot(results, labels[test.ix], type='bar')
# plot(results, labels[test.ix], type='dist')
# plot(results.train, labels[train.ix], type='pie')
# plot(results.train, labels[train.ix], type='bar')
# plot(results.train, labels[train.ix], type='dist')

# plot results with legend
# plot(results, labels[test.ix], type='pie', include.legend=TRUE, env.colors=c('#47697E','#5B7444','#CC6666','#79BEDB','#885588'))


res.mean <- apply(results$full.results,c(2,3,4),mean)

# Get depth of each sample for relative abundance calculation
sample.depths <- apply(results.train$full.results[1,,,,drop=F],4,sum)
#sample.depths <- rep(1073,19)
# create directory to store the results
outdir='All_OTU_Sourcetracker_1073'
filebase='m'
subdir <- paste(outdir,'full_results',sep='/')
dir.create(subdir,showWarnings=FALSE, recursive=TRUE)

# write each environment as a separate file
for(i in 1:length(results$train.envs)){
  env.name <- results$train.envs[i]
  filename.fractions <- sprintf('%s/%s_%s_contributions.txt', subdir, filebase, env.name)
  res.mean.i <- res.mean[i,,]
  # handle the case where there is only one sink sample
  if(is.null(dim(res.mean.i))) res.mean.i <- matrix(res.mean.i,ncol=1)
  
  # make rows be samples, columns be features
  res.mean.i <- t(res.mean.i)
  
  # ensure proper names are retained
  colnames(res.mean.i) <- colnames(otus)
  rownames(res.mean.i) <- results$samplenames
  
  # calculate and save relative abundance
  #res.mean.i.ra <- sweep(res.mean.i,1,sample.depths,'/')
  res.mean.i.ra <- sweep(res.mean.i,1,1073,'/')
  sink(filename.fractions)
  cat('SampleID\t')
  write.table(res.mean.i.ra,quote=F,sep='\t')
}
save.image("All_OTUs_RunSourceTracker.RData")

res.mean.unknown <- res.mean[3,,]
res.mean.unknown <- t(res.mean.unknown)
colnames(res.mean.unknown) <- colnames(otus)
rownames(res.mean.unknown) <- results$samplenames
#res.mean.unknown.ra <- sweep(res.mean.unknown,1,1073,'/')
write.table(res.mean.unknown,file = "All_OTU_Sourcetracker_1073/full_results/Unknown.conts.txt",col.names = NA,quote=F,sep='\t')


data.for.grouped.plot <- data.frame(rbind(results$proportions,results.train$proportions))
data.for.grouped.plot <- merge(data.for.grouped.plot,metadata,by = "row.names")
library(ggplot2);library(reshape2)
data.for.grouped.plot <- subset(data.for.grouped.plot,select = c("HO","SWAB","Unknown","Species"))
grouped.plot.input <- melt(data.for.grouped.plot)
grouped.plot.input$group <- paste(grouped.plot.input$Species,grouped.plot.input$variable,sep = ".")
ggplot(grouped.plot.input,aes(Species,value,group=group,fill=variable)) + geom_boxplot()
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE, conf.interval=.95) {
  library(doBy)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # Collapse the data
  formula <- as.formula(paste(measurevar, paste(groupvars, collapse=" + "), sep=" ~ "))
  datac <- summaryBy(formula, data=data, FUN=c(length2,mean,sd), na.rm=na.rm)
  
  # Rename columns
  names(datac)[ names(datac) == paste(measurevar, ".mean",    sep="") ] <- measurevar
  names(datac)[ names(datac) == paste(measurevar, ".sd",      sep="") ] <- "sd"
  names(datac)[ names(datac) == paste(measurevar, ".length2", sep="") ] <- "N"
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  print(datac$se)
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  print(datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
senew <- summarySE(grouped.plot.input,measurevar = "value",groupvars = c("Species","variable"))
senew$group <- paste(senew$Species,senew$variable,sep = ".")
senew$Species <- factor(senew$Species,levels = c("HO","SWAB","human","mice"))
p <- ggplot(senew, aes(x=Species, y=value, fill=variable))
#if (c == 0) {c = max(senew$value)}
p <- p + #geom_errorbar(aes(ymin=value-se, ymax=value+se),width=.2,fatten =0.4,
          #             position=position_stack()) +
  geom_bar(stat="identity",width = 0.8) +
  theme_classic()+ 
  theme(legend.position = "none") + theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90,size = 18),
        axis.title = element_blank(),
        axis.ticks.x = element_blank())
p
ggsave("sourcetracker.percentage.grouped.pdf",width = 3,height = 2.5)
