library(tidyverse)

dir <- "/COMMA/results/GWAS/s-LDSC/partitioned_heritability"

x <- read.table(file.path(dir,"partitioned_h2.results"),header=T,check.names=F)
v <- read.csv(file.path(dir,"annotation_validation.csv"))

names <- c("Arterial_Endothelial_Cell","Vascular_Smooth_Muscle_Cell",
           "Capillary_Endothelial_Cell","Pericyte","Venous_Endothelial_Cell","All_OCR")

SD <- v$SD[match(names,v$Annotation)]
stopifnot(nrow(x)==6,all(is.finite(SD)))

M <- 5698384 # Number of LDSC SNPs
h2 <- 0.0157 # Observed h2 from log file

x <- x %>%
  mutate(
    CellType=factor(c("Arterial EC","VSMC","Capillary EC",
                      "Pericyte","Venous EC","All OCR"),
                    levels=c("Arterial EC","VSMC","Capillary EC",
                             "Pericyte","Venous EC","All OCR")),
    SD=SD,
    EnrStat=(h2/M)*(Prop._h2/Prop._SNPs -
                    (1-Prop._h2)/(1-Prop._SNPs)),
    TauStar=Coefficient*SD*M/h2,
    TauStar.se=Coefficient_std_error*SD*M/h2,
    TauStar.p=2*pnorm(-abs(Coefficient/Coefficient_std_error))
  )

bg <- x[6,]

cells <- x[1:5,] %>%
  mutate(
    Enr.score=EnrStat/bg$EnrStat,
    sig=ifelse(TauStar.p<0.05,"*",""),
    xpos=1:5
  )

write.csv(x,file.path(dir,"ldsc_all_metrics.csv"),row.names=F)
write.csv(cells,file.path(dir,"ldsc_final.csv"),row.names=F)

# Bubble plot
p <- ggplot(cells,aes(xpos,1,size=Enr.score,fill=TauStar)) +
  geom_point(shape=21,colour="black",stroke=0.8) +
  geom_text(aes(y=1.35,label=sig),size=6) +
  geom_text(aes(y=0.65,label=sprintf("%.2fx",Enr.score)),size=3.5) +
  scale_size_area(max_size=18,limits=c(0.75,max(2,max(cells$Enr.score))),
                  breaks=c(0.75,1,1.25,1.5,2),name="Relative enrichment") +
  scale_fill_gradient2(low="#2166AC",mid="#F7F7F7",high="#B2182B",
                       midpoint=0,limits=c(-1,1),oob=scales::squish,
                       name=expression(tau^"*")) +
  scale_x_continuous(breaks=1:5,labels=c("Arterial EC","VSMC",
                      "Capillary EC","Pericyte","Venous EC")) +
  scale_y_continuous(limits=c(0.5,1.6),breaks=NULL) +
  coord_fixed(ratio=1) +
  labs(x=NULL,y=NULL,
       title="Vascular cell-type-specific heritability enrichment",
       subtitle="Bubble size: relative enrichment; colour: conditional tau*") +
  theme_bw(base_size=12) +
  theme(panel.grid=element_blank(),
        panel.border=element_rect(colour="black",fill=NA,linewidth=0.8),
        axis.text.x=element_text(angle=35,hjust=1),
        axis.ticks.y=element_blank())

# Figure 6B
ggsave(file.path(dir,"ldsc_bubble.png"),p,width=10,height=5,dpi=300)
ggsave(file.path(dir,"ldsc_bubble.pdf"),p,width=10,height=5)

# Enrichment bar chart
x$Group <- factor(
  ifelse(x$CellType=="All OCR","All vascular OCRs","Cell-type-specific OCRs"),
  levels=c("All vascular OCRs","Cell-type-specific OCRs")
)

# Figure 6A
p2 <- ggplot(x,aes(x=CellType,y=Enrichment,fill=Group)) +
  geom_col(width=0.7,colour="black",linewidth=0.3) +
  geom_errorbar(aes(ymin=pmax(0,Enrichment-Enrichment_std_error),
                    ymax=Enrichment+Enrichment_std_error),
                width=0.2,linewidth=0.5) +
  facet_grid(~Group,scales="free_x",space="free_x") +
  scale_fill_manual(values=c("Cell-type-specific OCRs"="#7393B3",
                             "All vascular OCRs"="#B0B0B0")) +
  labs(x=NULL,y="SNP heritability enrichment (± SE)",
       title="Heritability enrichment across vascular OCR annotations") +
  theme_classic(base_size=12) +
  theme(legend.position="none",
        axis.text.x=element_text(angle=35,hjust=1),
        strip.background=element_blank())

ggsave(file.path(dir,"ldsc_enrichment_bar.png"),p2,width=11,height=5,dpi=300)
ggsave(file.path(dir,"ldsc_enrichment_bar.pdf"),p2,width=11,height=5)

print(p)
print(p2)
print(select(cells,CellType,Enr.score,TauStar,TauStar.se,TauStar.p))