library(ggforce);library(ggplot2);library(cowplot);library(reshape)

read.table.meta = function(filename, ...){
  lines <- readLines(filename)
  n <- grep("^#", lines)
  if(length(n) > 0){start <- n[length(n)]}else{start <- 1}
  end <- length(lines)
  x <- read.table(text=lines[start:end],header=T,sep="\t",comment.char = '',check.names=F,...)
}
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

calculate.rate.trans.to.log <- function(x,dta,group){
  test <- aggregate(x,by=list(dta[,group]),function(y) mean(y,na.rm=T))
  log10(test[1,2]/test[2,2])
}
to.transform.to.fold.change.map <- function(data1,data2,map1,map2,groupname1,groupname2){
  data1 <- data.frame(t(data1),check.names = F)
  data1 <- merge(data1,map1,by = "row.names")
  data2 <- data.frame(t(data2),check.names = F)
  data2 <- merge(data2,map2,by = "row.names")
  int.bac <- intersect(colnames(data1),colnames(data2))
  int.data1 <- data1[,c(int.bac,groupname1)];row.names(int.data1) <-int.data1$Row.names;int.data1 <- int.data1[,-1]
  int.data2 <- data2[,c(int.bac,groupname2)];row.names(int.data2) <-int.data2$Row.names;int.data2 <- int.data2[,-1]
  k1 <- apply(int.data1[,1:c(ncol(int.data1)-3)], 2,function(x) calculate.rate.trans.to.log(x = x,dta = int.data1,group = groupname1))
  k2 <- apply(int.data2[,1:c(ncol(int.data2)-3)], 2,function(x) calculate.rate.trans.to.log(x = x,dta = int.data2,group = groupname2))
  k1 <- data.frame(k1);k2 <- data.frame(k2)
  k <- merge(k1,k2,by = "row.names")
  k[sapply(k,is.infinite)]<-NA;k <- na.omit(k)
  k$direction <- ifelse(k$k1*k$k2>0,"Positive","Oppsite")
  colnames(k) <- c("OTU","File1","File2","Direction")
  out <- melt(k,id.vars = c("OTU","Direction"))
  return(out)
}

mice_map <- read.table.meta("recipient.txt",row.names=1)
mice_genus <- read.table.meta("recipient_even15980_s1_L6.txt",row.names=1)
rownames(mice_genus) <- taxa_shortname(mice_genus,6)
mice_map <- subset(mice_map,group %in% c("NPFMT","PEFMT"))
mice_map$group <- factor(mice_map$group,levels = c("PEFMT","NPFMT"))


human_map <- read.table.meta("donor.txt",row.names=1)
human_genus <- read.table.meta("donor_even15980_s1_L6.txt",row.names=1)
rownames(human_genus) <- taxa_shortname(human_genus,6)
human_map$group <- factor(human_map$group,levels = c("preeclampsia","normotension"))

out <- to.transform.to.fold.change.map(data1 = human_genus,data2 = mice_genus,
                                map1 = human_map,map2 = mice_map,
                                groupname1 = "group",groupname2 = "group")

out1 <- cast(data = out,OTU+Direction~variable)
rownames(out1) <- out1$OTU
out1$OTU <- factor(rank(out1$File1))
out2 <- melt(out1,id.vars = c("OTU","Direction"))

out2$variable <- gsub("File1","Human",out2$variable)
out2$variable <- gsub("File2","Mice",out2$variable)
colnames(out2)[2] <- "Tendency"
out2$Tendency <- gsub("Positive","Consistent",out2$Tendency)
out2$Tendency <- gsub("Oppsite","Inconsistent",out2$Tendency)
 

ggplot(out2,aes(x = value,y = OTU,color=Tendency,shape=variable)) +  
  geom_point() + theme_bw() + xlim(-2.2,2.2) + labs(title="Overlap Genus",x="",y="") + 
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),plot.title = element_text(hjust = 0.5))
ggsave("genus.donor.recipient.pdf",width = 2.8,height = 2.5)
