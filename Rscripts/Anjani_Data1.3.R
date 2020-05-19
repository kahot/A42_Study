library(reshape2)
library(ggplot2)
library(cowplot)
library(colorspace)

#Read 1.3ug data file
aa42<-read.csv("Data/a42_mutants1.3.csv", stringsAsFactors = F)
aa42$Mutants<-gsub(" .*","",aa42$Mutants)
aa42[1, 3:4]<-"Control"
aa42[5, 3:4]<-"WT"
#extract the AA and Property info
aa42_2<-aa42[,2:4]

while(length(ind <- which(aa42$Mutated.Amino.Acid == "")) > 0){
    aa42$Mutated.Amino.Acid[ind] <- aa42$Mutated.Amino.Acid[ind -1]
}
while(length(ind <- which(aa42$Properties.of.Amino.Acid == "")) > 0){
    aa42$Properties.of.Amino.Acid[ind] <- aa42$Properties.of.Amino.Acid[ind -1]
}

#extract the AA and Property info
df<-aa42_2[which(aa42_2$Mutated.Amino.Acid!=''),]
#rename the column 
colnames(df)[2:3]<-c("AA","Property")



# calculate the mean and sd for each antibiotics
antibio<-colnames(aa42)[5:8] 
sum<-list()
for (i in 1:length(antibio)){
    df1<-aggregate(aa42[,antibio[i]], list(aa42$Mutants), mean, na.rm=T)
    df2<-aggregate(aa42[,antibio[i]], list(aa42$Mutants), sd, na.rm=T)
    df1$SD<-df2$x
    colnames(df1)[1:2]<-c("Mutants","Mean")
    df1$Antibiotics<-antibio[i]
    sum[[i]]<-df1
}


Summary<-data.frame(do.call(rbind,sum))

Summary<-merge(df, Summary, by="Mutants")
Summary$Mutants<-factor(Summary$Mutants, levels=paste(df$Mutants))

#create colors:
colors<-qualitative_hcl(9, palette="Set2")

#create a plot
for (i in 1:4){
    ab<-antibio[i]
    p<-ggplot(Summary[Summary$Antibiotics==ab,])+
        geom_bar(aes(x=Mutants, y=Mean, fill=Property), position=position_dodge(.9), stat="identity",width=0.8)+
        geom_errorbar(aes(x=Mutants, y=Mean, ymin=Mean-SD, ymax=Mean+SD), position=position_dodge(.9), width=.2, color="gray30")+
        theme_linedraw()+
        ylab("Zone of inhibition (mm)")+
        labs(x="")+  
        theme(axis.text.x = element_text(angle = 90),
              legend.title = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(linetype=3),
              panel.grid.minor.y = element_line(linetype=3))+
        scale_fill_manual(values=colors)
        
    pname<-paste0("Plot_",i)
    assign(pname,p)
}
title <- ggdraw() + draw_label("")
plot_grid(title, title,Plot_1,Plot_2,Plot_3,Plot_4, nrow = 3, ncol=2, labels = c('','',antibio),label_size = 12,
         label_y=1.05, rel_heights = c(0.2, 1,1))
ggsave(filename="Output/Mean.byMutatnts.byAntibio.1.3.pdf",width =12, height =8)




# create a plot with mean by property type:
colnames(aa42)
AA<-melt(aa42, id.vars = "Properties.of.Amino.Acid", measure.vars=c("Amikacin","Kanamycin","Tobramycin","Neomycin" ))

propSummary<-aggregate(.~Properties.of.Amino.Acid+variable, AA, mean, na.rm=T)
propertySD<- aggregate(.~Properties.of.Amino.Acid+variable, AA, sd, na.rm=T)
propSummary$SD<-propertySD$value
colnames(propSummary)[1:3]<-c("Property", "Antibiotics", "Mean")
propSummary$Property<-factor(propSummary$Property, levels=c("Control","WT","Acidic","Aliphatic","Aromatic","Basic","Cyclic","Neutral","Sulfur containing"))

col4<-qualitative_hcl(4, palette="Set2")
ggplot(propSummary)+
    geom_bar(aes(x=Property, y=Mean, fill=Antibiotics ), position=position_dodge(.9), stat="identity",width=0.8)+
    geom_errorbar(aes(x=Property, y=Mean, ymin=Mean-SD, ymax=Mean+SD,fill=Antibiotics), position=position_dodge(.9), width=.2, color="gray30")+
    theme_bw()+
    ylab("Zone of inhibition (mm)")+
    labs(x="")+  
    theme(axis.text.x = element_text(angle = 90, size=12, hjust=1),
          axis.text.y = element_text(size=10, color=1),
          legend.title = element_blank(),
          panel.grid.major.x = element_blank())+
    scale_fill_manual(values=col4)
ggsave(filename="Output/Mean.byProperty.byAntibiotics1.3.pdf",width =7, height =5)

####
wilcox.res<-data.frame(Property=property)
for ( i in 1:7){
    result<-wilcox.test(AA$value[AA$Properties.of.Amino.Acid=="WT"],AA$value[AA$Properties.of.Amino.Acid==property[i]], alternative ="less", paired=FALSE )
    wilcox.res$W[i]<-result[[1]]
    wilcox.res$P_value[i]<-result[[3]]
    wilcox.res$SamplesSie[i]<-length(AA$value[AA$Properties.of.Amino.Acid==property[i]])
}

write.csv(wilcox.res, "Output/Wilcox.results.compareToWT_1.3.csv")

