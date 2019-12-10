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


contamdf.prev <- isContaminant(subset_samples(ps, neg | plac), neg="neg", 
                               method="prevalence", detailed=TRUE,threshold = 0.5,normalize = TRUE)
table(contamdf.prev$contaminant)
# Make data.frame of prevalence in positive and negative samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(samdf$neg == "TRUE", ps.pa)
ps.pa.pos <- prune_samples(samdf$neg == "FALSE", ps.pa)
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
df.pa.na.omit <- na.omit(df.pa)
ggplot(data=df.pa.na.omit, aes(x=pa.neg, y=pa.pos,color=contaminant)) + geom_point() + theme_classic() + 
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave("decontam.pdf",width = 3.8,height = 3)

otu.contam <- rownames(contamdf.prev)[which(contamdf.prev$contaminant != "TRUE")]
ps.del.Contam <- prune_taxa(taxa = otu.contam,ps)

out <- isNotContaminant(subset_samples(ps.del.Contam, neg | plac), neg="neg", 
                        method="prevalence", detailed=TRUE,threshold = 0.5,normalize = TRUE)
# out$ind <- seq(nrow(out))
# ggplot(data=out, aes(x=ind, y=p)) + geom_point() + scale_y_log10() + geom_hline(yintercept=0.1, color="red") + geom_hline(yintercept=0.01, color="red") + theme_bw()
# summary(out$p)
table(out$not.contaminant)
ps.pa <- transform_sample_counts(ps.del.Contam, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(samdf$neg == "TRUE", ps.pa)
ps.pa.pos <- prune_samples(samdf$neg == "FALSE", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    notcontaminant=out$not.contaminant)
df.pa.na.omit <- na.omit(df.pa)
ggplot(data=df.pa.na.omit, aes(x=pa.neg, y=pa.pos,color=notcontaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

otu.contam2 <- rownames(out)[which(out$not.contaminant != "FALSE")]
ps.del.Contam2 <- prune_taxa(taxa = otu.contam2,ps.del.Contam)

ps.pa <- transform_sample_counts(ps.del.Contam2, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(samdf$neg == "TRUE", ps.pa)
ps.pa.pos <- prune_samples(samdf$neg == "FALSE", ps.pa)
# Make data.frame of prevalence in positive and negative samples
filter.otu <- subset(out,not.contaminant=="TRUE")
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    notcontaminant=filter.otu$not.contaminant)
df.pa.na.omit <- na.omit(df.pa)
ggplot(data=df.pa.na.omit, aes(x=pa.neg, y=pa.pos,color=notcontaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

write.table(filter.otu,file = "decontam.identify.otus.txt",sep = "\t",quote = F,col.names = NA)

### use all otus
all.otus <- read.table("merge.1st.3rd.biom.txt",header = T,row.names = 1,sep = "\t")

only.decomt.filter.otus <- all.otus[rownames(filter.otu),]
write.table(only.decomt.filter.otus,file = "only.decomt.filter.otus.txt",sep = "\t",quote = F,col.names = NA)



### use Sourcetraker filter otus
st.otus <- read.table("sourcetraker.filtered.otus.txt", header=TRUE, stringsAsFactors = FALSE,row.names = 1,sep = "\t",comment.char = "")

combine.filter.otus <- st.otus[rownames(filter.otu),]
write.table(combine.filter.otus,file = "combine.filter.otus.txt",sep = "\t",quote = F,col.names = NA)

