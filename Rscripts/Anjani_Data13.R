library(reshape2)
library(ggplot2)
library(cowplot)
library(colorspace)
#create colors:
colors<-qualitative_hcl(9, palette="Set2")


## Read 13ug file

a42<-read.csv("Data/A42_mutants13.csv", stringsAsFactors = F)
a42$Mutants<-gsub(" .*","",a42$Mutants)
a42[1, 3:4]<-"Control"
a42[5, 3:4]<-"WT"
#extract the AA and Property info
a42_2<-a42[,2:4]

while(length(ind <- which(a42$Mutated.Amino.Acid == "")) > 0){
    a42$Mutated.Amino.Acid[ind] <- a42$Mutated.Amino.Acid[ind -1]
}
while(length(ind <- which(a42$Properties.of.Amino.Acid == "")) > 0){
    a42$Properties.of.Amino.Acid[ind] <- a42$Properties.of.Amino.Acid[ind -1]
}

#extract the AA and Property info
df<-a42_2[which(a42_2$Mutated.Amino.Acid!=''),]
#rename the column 
colnames(df)[2:3]<-c("AA","Property")



# calculate the mean and sd for each antibiotics
antibio<-colnames(a42)[5:8] 
sum<-list()
for (i in 1:length(antibio)){
    df1<-aggregate(a42[,antibio[i]], list(a42$Mutants), mean, na.rm=T)
    df2<-aggregate(a42[,antibio[i]], list(a42$Mutants), sd, na.rm=T)
    df1$SD<-df2$x
    colnames(df1)[1:2]<-c("Mutants","Mean")
    df1$Antibiotics<-antibio[i]
    sum[[i]]<-df1
}


Summary<-data.frame(do.call(rbind,sum))

Summary<-merge(df, Summary, by="Mutants")
Summary$Mutants<-factor(Summary$Mutants, levels=paste(df$Mutants))


#create a plot
for (i in 1:4){
    ab<-antibio[i]
    p<-ggplot(Summary[Summary$Antibiotics==ab,])+
        geom_bar(aes(x=Mutants, y=Mean, fill=Property, ), position=position_dodge(.9), stat="identity",width=0.8)+
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
ggsave(filename="Output/Mean.byMutatnts.byAntibio.13.pdf",width =12, height =8)




# create a plot with mean by property type:
aa<-melt(a42, id.vars = "Properties.of.Amino.Acid", measure.vars=c("Amikacin","Kanamycin","Tobramycin","Neomycin" ))

propSum<-aggregate(.~Properties.of.Amino.Acid+variable, aa, mean, na.rm=T)
propSD<-aggregate(.~Properties.of.Amino.Acid+variable, aa, sd, na.rm=T)
propSum$SD<-propSD$value
colnames(propSum)[1:3]<-c("Property", "Antibiotics", "Mean")
propSum$Property<-factor(propSum$Property, levels=c("Control","WT","Acidic","Aliphatic","Aromatic","Basic","Cyclic","Neutral","Sulfur containing"))

col4<-qualitative_hcl(4, palette="Set2")
ggplot(propSum)+
    geom_bar(aes(x=Property, y=Mean, fill=Antibiotics ), position=position_dodge(.9), stat="identity",width=0.8)+
    geom_errorbar(aes(x=Property, y=Mean, ymin=Mean-SD, ymax=Mean+SD,fill=Antibiotics), position=position_dodge(.9), width=.2, color="gray30")+
    theme_bw()+
    ylab("Zone of inhibition (mm)")+
    labs(x="")+  
    theme(axis.text.x = element_text(angle = 90, size=12, hjust=1),
          axis.text.y = element_text( size=10, color=1),
          legend.title = element_blank(),
          panel.grid.major.x = element_blank())+
    scale_fill_manual(values=col4)
ggsave(filename="Output/Mean.byProperty.byAntibio13.pdf",width =7, height =5)


######
# Remove neomycin and average over 3 antibiotics:

aa2<-aa[aa$variable!="Neomycin",]

propSum2<-aggregate(aa$value, list(aa$Properties.of.Amino.Acid), mean, na.rm=T)
propSD2<-aggregate(aa$value, list(aa$Properties.of.Amino.Acid), sd, na.rm=T)
propSum2$SD<-propSD2$x
colnames(propSum2)[1:3]<-c("Property", "Mean","SD")

#order 'control', 'wt' and acending order
dt<-propSum2[propSum2$Property!="Control" & propSum2$Property!="WT",]
dt<-dt[order(dt$Mean),]
propSum2$Property<-factor(propSum2$Property, levels=c("Control","WT",paste(dt$Property)))

#col4<-qualitative_hcl(4, palette="Set2")
ggplot(propSum2)+
    geom_bar(aes(x=Property, y=Mean ), stat="identity",width=0.8, fill=colors[6])+
    geom_errorbar(aes(x=Property, y=Mean, ymin=Mean-SD, ymax=Mean+SD), width=.2, color="gray30")+
    theme_bw()+
    ylab("Zone of inhibition (mm)")+
    labs(x="")+  
    theme(axis.text.x = element_text(angle = 90, size=12, hjust=1, color=1),
          axis.text.y = element_text( size=10, color=1),
          legend.title = element_blank(),
          panel.grid.major.x = element_blank())
ggsave(filename="Output/Mean.byProperty_grouped.13.pdf",width =6, height =5)


#######
#run stats on some of them:

#check assumption
#"Shapiro-Wilk normality test".
shapiro.test(aa$value) 
#W = 0.94196, p-value = 4.716e-09

# Levene's Test 
library(Rcmdr) #(testing homegeneity of variance)
leveneTest(value ~ Properties.of.Amino.Acid, data=aa)
#       Df F value   Pr(>F)   
#group   8   3.158 0.001948 **


property<-paste(propSum$Property[1:9])
property<-property[c(-5,-9)]


# 1. WT vs Cylic AA
wilcox.results<-data.frame(Property=property)
for ( i in 1:7){
    result<-wilcox.test(aa2$value[aa2$Properties.of.Amino.Acid=="WT"],aa2$value[aa2$Properties.of.Amino.Acid==property[i]], alternative ="less", paired=FALSE )
    wilcox.results$W[i]<-result[[1]]
    wilcox.results$P_value[i]<-result[[3]]
    }

write.csv(wilcox.results, "Output/Wilcox.results.compareToWT.csv")


