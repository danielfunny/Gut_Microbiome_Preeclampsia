library(psych);library(corrplot);library(ggplot2)
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

adiv <- read.table("adiv.txt",header = T,row.names = 1,sep = "\t",comment.char = "")
all_genus <- read.table("cxhuman0814_even15980_sorted_L6.txt",header=T,row.names=1,sep="\t")
rownames(all_genus) <- taxa_shortname(all_genus,6)
map.all <- read.table("cxhuman0814.txt",header = T,row.names = 1,sep = "\t",comment.char = "")


#------------------
# alpha diversity
#########-----------------
## alpha diversity of PE & NP
m <- merge(adiv,map.all,by = "row.names")
wilcox.test(m$shannon~m$Statue2)
pd <- position_dodge(0.5)
p.sh <- ggplot(m,aes(Statue2,shannon,fill=Statue2)) + 
  geom_jitter(aes(colour=Statue2),size=0.005,width = 0.1) +
  stat_boxplot(geom = "errorbar",width=0.2,position=pd,lwd=0.4,fatten=0.4)+
  geom_boxplot(position=pd,width=0.7,outlier.shape = NA,lwd=0.2) +
  scale_color_manual(values = c("normotension"="#80cbfb","preeclampsia"="#ffa2a3")) +
  scale_fill_manual(values = c("normotension"="#80cbfb","preeclampsia"="#ffa2a3")) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank()) 
ggsave(filename = "adiv.shannon.pdf",plot = p.sh,width = 1,height = 1.2)


wilcox.test(m$PD_whole_tree~m$Statue2)
pd <- position_dodge(0.5)
p.pd <- ggplot(m,aes(Statue2,PD_whole_tree,fill=Statue2)) + 
  geom_jitter(aes(colour=Statue2),size=0.005,width = 0.1) +
  stat_boxplot(geom = "errorbar",width=0.2,position=pd,lwd=0.4,fatten=0.4)+
  geom_boxplot(position=pd,width=0.7,outlier.shape = NA,lwd=0.2) +
  scale_color_manual(values = c("normotension"="#80cbfb","preeclampsia"="#ffa2a3")) +
  scale_fill_manual(values = c("normotension"="#80cbfb","preeclampsia"="#ffa2a3")) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank()) 
ggsave(filename = "adiv.PD_whole_tree.pdf",plot = p.pd,width = 1,height = 1.2)


wilcox.test(m$observed_otus~m$Statue2)
pd <- position_dodge(0.5)
p.ob <- ggplot(m,aes(Statue2,observed_otus,fill=Statue2)) + 
  geom_jitter(aes(colour=Statue2),size=0.005,width = 0.1) +
  stat_boxplot(geom = "errorbar",width=0.2,position=pd,lwd=0.4,fatten=0.4)+
  geom_boxplot(position=pd,width=0.7,outlier.shape = NA,lwd=0.2) +
  scale_color_manual(values = c("normotension"="#80cbfb","preeclampsia"="#ffa2a3")) +
  scale_fill_manual(values = c("normotension"="#80cbfb","preeclampsia"="#ffa2a3")) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank()) 
ggsave(filename = "adiv.observed_otus.pdf",plot = p.ob,width = 1,height = 1.2)


wilcox.test(m$chao1~m$Statue2)
pd <- position_dodge(0.5)
p.ch <- ggplot(m,aes(Statue2,chao1,fill=Statue2)) + 
  geom_jitter(aes(colour=Statue2),size=0.005,width = 0.1) +
  stat_boxplot(geom = "errorbar",width=0.2,position=pd,lwd=0.4,fatten=0.4)+
  geom_boxplot(position=pd,width=0.7,outlier.shape = NA,lwd=0.2) +
  scale_color_manual(values = c("normotension"="#80cbfb","preeclampsia"="#ffa2a3")) +
  scale_fill_manual(values = c("normotension"="#80cbfb","preeclampsia"="#ffa2a3")) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank()) 
ggsave(filename = "adiv.chao1.pdf",plot = p.ch,width = 1,height = 1.2)



########------------------
## alpha diversity of subgroups 
all <- read.table("PE.subgroups.taxa.txt",header = T,row.names = 1,sep = "\t",comment.char = "")
m <- merge(adiv,all,by = "row.names")
wilcox.test(m$shannon~m$Statue3)
pd <- position_dodge(0.5)
p.sh <- ggplot(m,aes(Statue3,shannon,fill=Statue3)) + 
  geom_jitter(aes(colour=Statue3),size=0.005,width = 0.1) +
  stat_boxplot(geom = "errorbar",width=0.2,position=pd,lwd=0.4,fatten=0.4)+
  geom_boxplot(position=pd,width=0.7,outlier.shape = NA,lwd=0.2) +
  scale_color_manual(values = c("mild_preeclampsia"="#ffa2a3","severe_preeclampsia"="#A83434")) +
  scale_fill_manual(values = c("mild_preeclampsia"="#ffa2a3","severe_preeclampsia"="#A83434")) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90,size=12),
        axis.title = element_blank()) 
ggsave(filename = "adiv.shannon.pdf",plot = p.sh,width = 1,height = 1.2)


wilcox.test(m$PD_whole_tree~m$Statue3)
pd <- position_dodge(0.5)
p.pd <- ggplot(m,aes(Statue3,PD_whole_tree,fill=Statue3)) + 
  geom_jitter(aes(colour=Statue3),size=0.005,width = 0.1) +
  stat_boxplot(geom = "errorbar",width=0.2,position=pd,lwd=0.4,fatten=0.4)+
  geom_boxplot(position=pd,width=0.7,outlier.shape = NA,lwd=0.2) +
  scale_color_manual(values = c("mild_preeclampsia"="#ffa2a3","severe_preeclampsia"="#A83434")) +
  scale_fill_manual(values = c("mild_preeclampsia"="#ffa2a3","severe_preeclampsia"="#A83434")) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90,size=12),
        axis.title = element_blank()) 
ggsave(filename = "adiv.PD_whole_tree.pdf",plot = p.pd,width = 1,height = 1.2)


wilcox.test(m$observed_otus~m$Statue3)
pd <- position_dodge(0.5)
p.ob <- ggplot(m,aes(Statue3,observed_otus,fill=Statue3)) + 
  geom_jitter(aes(colour=Statue3),size=0.005,width = 0.1) +
  stat_boxplot(geom = "errorbar",width=0.2,position=pd,lwd=0.4,fatten=0.4)+
  geom_boxplot(position=pd,width=0.7,outlier.shape = NA,lwd=0.2) +
  scale_color_manual(values = c("mild_preeclampsia"="#ffa2a3","severe_preeclampsia"="#A83434")) +
  scale_fill_manual(values = c("mild_preeclampsia"="#ffa2a3","severe_preeclampsia"="#A83434")) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90,size=12),
        axis.title = element_blank()) 
ggsave(filename = "adiv.observed_otus.pdf",plot = p.ob,width = 1,height = 1.2)


wilcox.test(m$chao1~m$Statue3)
pd <- position_dodge(0.5)
p.ch <- ggplot(m,aes(Statue3,chao1,fill=Statue3)) + 
  geom_jitter(aes(colour=Statue3),size=0.005,width = 0.1) +
  stat_boxplot(geom = "errorbar",width=0.2,position=pd,lwd=0.4,fatten=0.4)+
  geom_boxplot(position=pd,width=0.7,outlier.shape = NA,lwd=0.2) +
  scale_color_manual(values = c("mild_preeclampsia"="# ","severe_preeclampsia"="# ")) +
  scale_fill_manual(values = c("mild_preeclampsia"="#ffa2a3","severe_preeclampsia"="#A83434")) +
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y = element_text(angle = 90,size=12),
        axis.title = element_blank()) 
ggsave(filename = "adiv.chao1.pdf",plot = p.ch,width = 1,height = 1.2)






#------------------
# Correlation
####---------------------------
### correlation of all samples
sid <- intersect(colnames(all_genus),rownames(map.all))
all1 <- t(all_genus[c("Clostridium","Dialister","Veillonella","Fusobacterium",
                      "Lachnospira","Akkermansia","Faecalibacterium"),sid])
all_ind <- map.all[sid,c("sbp","dbp","proteinuria","edema","ALT","AST","Cr","TT","ALB","TP","ATIII","neoweight","hemorrhage")]

all_res <- corr.test(all1,all_ind,method = "spearman",adjust = "none")
p <- data.frame(t(all_res$p))
p1 <- sapply(p,p.adjust,method="fdr")
rownames(p1) <- rownames(p)
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
colors <- c("#67001E","#67001D",col2(26),"#053063","#053062")
out.corr <- corrplot(all_res$r,col = rev(colors),
                     p.mat = t(p1), sig.level = .05, insig = "blank",tl.col = "black")

####---------------------------
### correlation of only PE samples
sub_genus <- all_genus[,grepl("PIH",colnames(all_genus))]
#all_genus <- t(all_genus)
ind <- read.table("cxhuman0814.txt",header = T,row.names = 1,sep = "\t",comment.char = "",check.names = F)
ind <- ind[grepl("PIH",rownames(ind)),]
sid <- intersect(colnames(sub_genus),rownames(ind))
all1 <- t(sub_genus[c("Clostridium","Dialister","Veillonella","Fusobacterium",
                      "Lachnospira","Akkermansia","Faecalibacterium"),sid])
all_ind <- ind[sid,c("sbp","dbp","proteinuria","edema","ALT","AST","Cr","TT","ALB","TP","ATIII","neoweight","hemorrhage")]

all_res <- corr.test(all1,all_ind,method = "spearman",adjust = "none")
p <- data.frame(t(all_res$p))
p1 <- sapply(p,p.adjust,method="fdr")
rownames(p1) <- rownames(p)
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
colors <- c("#67001E","#67001D",col2(26),"#053063","#053062")
out.corr <- corrplot(all_res$r,col = rev(colors),
                     p.mat = t(p1), sig.level = .05, insig = "blank",tl.col = "black")





#------------------
# genus compare between PE & NP
tomakeplot <- function(dta,value,variable,c=0){
  pd <- position_dodge(0.8)
  p.out <- wilcox.test(dta[,value]~dta[,variable])
  print(p.out)
  dta[,value] <- dta[,value]*100
  p <- ggplot(dta,aes_string(variable,value,fill=variable)) + 
    geom_jitter(aes_string(colour=variable),size=0.005,position=pd) +
    stat_boxplot(geom = "errorbar",width=0.2,position=pd,lwd=0.4,fatten=0.4)+
    geom_boxplot(position=pd,width=0.7,outlier.shape = NA,lwd=0.2) +
    scale_fill_manual(values = c("normotension"="#80cbfb","preeclampsia"="#ffa2a3")) +
    scale_color_manual(values = c("normotension"="#80cbfb","preeclampsia"="#ffa2a3")) +
    theme(legend.position = "none") + theme_classic()+ ylab(paste(value,"(%)",sep = ""))+
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_text(angle = 90,size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(angle = 90)) +
    scale_y_continuous(limits = c(0,c))
  filename <- paste("two.groups",value,"pdf",sep = ".")
  ggsave(filename = filename,plot = p,width = 2.1,height = 2.4)
}

sid <- intersect(colnames(all_genus),rownames(map.all))
all1 <- t(all_genus[c("Clostridium","Dialister","Veillonella","Fusobacterium","Coprococcus",
                      "Lactococcus","Lachnospira","Akkermansia","Faecalibacterium"),sid])
m <- merge(all1,map.all,by = "row.names")
tomakeplot(dta = m,value = "Clostridium",variable = "Statue2",c=20)
tomakeplot(dta = m,value = "Dialister",variable = "Statue2",c=4)
tomakeplot(dta = m,value = "Veillonella",variable = "Statue2",c=1)
tomakeplot(dta = m,value = "Fusobacterium",variable = "Statue2",c=0.2)
tomakeplot(dta = m,value = "Coprococcus",variable = "Statue2",c=5)
tomakeplot(dta = m,value = "Lachnospira",variable = "Statue2",c=3)
tomakeplot(dta = m,value = "Akkermansia",variable = "Statue2",c=0.5)
tomakeplot(dta = m,value = "Faecalibacterium",variable = "Statue2",c=10)

#------------------
# picrust
description <- read.table("metagenome_prediction.description_table.txt",sep = "\t",header = T,quote = "",row.names = 1,comment.char = "")
description <- description[,-ncol(description)]; description <- data.frame(t(description))
map.for.picrust <- data.frame(group =map.all$Statue2,row.names = rownames(map.all))
m <- merge(map.for.picrust,description,by="row.names")
n=1
out <- data.frame(name=1,p.value=1,mean.NP=1,mean.PE=1)
for (i in 3:ncol(m)){
  wil <- wilcox.test(m[,i]~m$group)
  if (!is.na(wil$p.value)){
    if (sum(m[,i])>1){
      out[n,1]= colnames(m)[i]
      out[n,2]= wil$p.value
      out[n,3]= aggregate(m[,i],by=list(m$group),mean)[1,2]
      out[n,4]= aggregate(m[,i],by=list(m$group),mean)[2,2]
      n= n+1
    }
  }
}
out$p.fdr <- p.adjust(out$p.value,method = "fdr")
