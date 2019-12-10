library(pheatmap)
taxa_shortname <- function(otu,level){
  if("taxonomy" %in% colnames(otu)){
    otu_tax <- otu[,"taxonomy"]
  }else{
    otu_tax <- rownames(otu)
  }
  otu_tax <- gsub(" ","",otu_tax)
  tax.16s <- strsplit(as.character(otu_tax),split = ";",perl = T)
  if(is.numeric(level))
    tax <- data.frame(matrix(unlist(tax.16s), ncol = level, byrow=T),stringsAsFactors=FALSE)
  if(level=="OTU" | level=="otu")
    tax <- data.frame(matrix(unlist(tax.16s), ncol = 7, byrow=T),stringsAsFactors=FALSE)
  tax1 <- apply(tax,2,function(x) gsub(".__","",x,perl = T))
  if(level %in% c("OTU","otu",7)){
    newtax <- apply(tax1,1,function(x) ifelse(all(nchar(x)>0),paste(x[6],x[7],sep = "."),x[which(nchar(x)<1)-1]))
    if(level == 7)
      newtax[which(duplicated(newtax))] <- apply(tax1[which(duplicated(newtax)),],1,function(x) paste(x[which(nchar(x)>1)[length(which(nchar(x)>1))-1]],x[which(nchar(x)>1)[length(which(nchar(x)>1))]],sep="."))
  }
  if(level %in% c(2:6)){
    newtax <- apply(tax1,1,function(x) x[which(nchar(x)>1 & x!="Other")[length(which(nchar(x)>1 & x!="Other"))]])
    if(length(which(duplicated(newtax)))>1){ 
      newtax[which(duplicated(newtax))] <- apply(tax1[which(duplicated(newtax)),],1,function(x) paste(x[which(nchar(x)>1)[length(which(nchar(x)>1))-1]],x[which(nchar(x)>1)[length(which(nchar(x)>1))]],sep="."))
    }else if(length(which(duplicated(newtax)))==1){
      newtax[which(duplicated(newtax))] <- paste(tax1[which(duplicated(newtax)),which(nchar(tax1[which(duplicated(newtax)),])>1)[length(which(nchar(tax1[which(duplicated(newtax)),])>1))-1]],tax1[which(duplicated(newtax)),which(nchar(tax1[which(duplicated(newtax)),])>1)[length(which(nchar(tax1[which(duplicated(newtax)),])>1))]],sep=".")
    }
  }
  if(level=="OTU"| level=="otu")
    newtax <- paste("OTU",rownames(otu),"_",newtax,sep = "")
  return(newtax)
}
read.table.meta = function(filename, ...){
  lines <- readLines(filename)
  n <- grep("^#", lines)
  if(length(n) > 0){start <- n[length(n)]}else{start <- 1}
  end <- length(lines)
  x <- read.table(text=lines[start:end],header=T,sep="\t",comment.char = '',check.names=F,...)
}

all_genus <- read.table.meta("donor_even15980_s1_L6.txt",row.names=1)
rownames(all_genus) <- taxa_shortname(all_genus,6)
all1 <- t(all_genus[c("Clostridium","Dialister","Veillonella","Phascolarctobacterium","Megamonas",
                      "Fusobacterium","Sutterella",
                      "Akkermansia","Faecalibacterium","Bifidobacterium","Streptococcus","Bacteroides","Coprococcus",
                      "Prevotella"),])
map <- read.table.meta("cxhuman0814.txt",row.names = 1)
sid <- intersect(rownames(map),rownames(all1))
dta <- data.frame(group=map[sid,"Statue2"],all1[sid,])
res <- data.frame(t(apply(dta[,-1],2,function(x) tapply(x, dta$group, mean))))
res$h_n_diff <- res$preeclampsia-res$normotension
res$group <- ifelse(res$h_n_diff>0,"preeclampsia","normotension")

donor_genus <- read.table.meta("donor_even15980_s1_L6.txt",row.names=1)
rownames(donor_genus) <- taxa_shortname(donor_genus,6)
donor1 <- t(donor_genus[c("Clostridium","Dialister","Veillonella","Phascolarctobacterium","Megamonas",
                      "Fusobacterium","Sutterella",
                      "Akkermansia","Faecalibacterium","Bifidobacterium","Streptococcus","Bacteroides","Coprococcus",
                      "Prevotella"),])
map <- read.table.meta("donor.txt",row.names = 1)
sid <- intersect(rownames(map),rownames(donor1))
dta <- data.frame(group=map[sid,"group"],donor1[sid,])
dres <- data.frame(t(apply(dta[,-1],2,function(x) tapply(x, dta$group, mean))))
dres$h_n_diff <- dres$preeclampsia-dres$normotension
dres$group <- ifelse(dres$h_n_diff>0,"preeclampsia","normotension")
colnames(dres) <- paste("donor",colnames(dres),sep = ".")
dres$donor.PE <- ifelse(dres$donor.h_n_diff>0,1,-1)
dres$donor.NP <- ifelse(dres$donor.h_n_diff>0,-1,1)

mice_genus <- read.table.meta("recipient_even15980_s1_L6.txt",row.names = 1)
rownames(mice_genus) <- taxa_shortname(mice_genus,level = 6)
mice_map <-  read.table("recipient.txt",header = T,row.names = 1,sep = "\t",comment.char = "")
mice_genus1 <- t(mice_genus[c("Clostridium","Dialister","Veillonella","Phascolarctobacterium","Megamonas",
                              "Fusobacterium","Sutterella",
                              "Akkermansia","Faecalibacterium","Bifidobacterium","Streptococcus","Bacteroides","Coprococcus",
                              "Prevotella"),])
mid <- intersect(rownames(mice_map),rownames(mice_genus1))
mdta <- data.frame(group=mice_map[mid,"group"],mice_genus1[mid,])
mres <- data.frame(t(apply(mdta[,-1],2,function(x) tapply(x, mdta$group, mean))))
mres$pf_nf_diff <- mres$PEFMT-mres$NPFMT
mres$PFMTgroup <- ifelse(mres$pf_nf_diff>0,1,-1)
mres$NFMTgroup <- ifelse(mres$pf_nf_diff>0,-1,1)

t <- cbind(mres,res)
t <- cbind(t,dres)
t1 <- na.omit(t)
t1$res <- ifelse(t1$h_n_diff*t1$pf_nf_diff>0,"overlap","nonoverlap")
annotation_row = data.frame(
  group = t1$res
)
ann_colors = list(group = c(overlap="#009845",nonoverlap="#E5E4F2"))
rownames(annotation_row) = rownames(t1)
pdf(file = "heatmap.donor.recipient.pdf",width = 4.8,height = 2.5)
pheatmap(t1[,c("donor.NP","donor.PE","NFMTgroup","PFMTgroup")],cluster_rows = T,cluster_cols = F,
         annotation_row = annotation_row,annotation_colors = ann_colors)
dev.off()

