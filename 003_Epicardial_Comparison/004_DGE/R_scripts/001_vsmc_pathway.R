library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)

setwd("/COMMA/results/meta_analysis/epicardial/")

genes <- read.csv("cell_consensus_genes.csv")

cardiac <- genes[genes$Group == "atlas", ]
epicardial <- genes[genes$Group == "epicardial", ]

cardiac_enriched <- enrichGO(
  gene = cardiac$Gene,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE,
  keyType = "SYMBOL"
)

cardiac_reduced <- simplify(cardiac_enriched, cutoff = 0.9, by = "qvalue",
                            select_fun = min, measure = "Wang")

epicardial_enriched <- enrichGO(
  gene = epicardial$Gene,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable = TRUE,
  keyType = "SYMBOL"
)

epicardial_reduced <- simplify(epicardial_enriched, cutoff = 0.9, by = "qvalue",
                               select_fun = min, measure = "Wang")

# Plots
card_df <- as.data.frame(cardiac_reduced) %>% arrange(desc(Count)) %>% slice_head(n=10) %>% mutate(Description=paste0(toupper(substr(Description,1,1)),substr(Description,2,nchar(Description))),Group="Cardiac")
epi_df  <- as.data.frame(epicardial_reduced) %>% arrange(desc(Count)) %>% slice_head(n=10) %>% mutate(Description=paste0(toupper(substr(Description,1,1)),substr(Description,2,nchar(Description))),Group="Epicardial")
max_count <- max(c(card_df$Count,epi_df$Count)); card_lim <- range(card_df$p.adjust); epi_lim <- range(epi_df$p.adjust)

make_dotplot <- function(df, cols) ggplot(df %>% arrange(Count) %>% mutate(Description=factor(str_wrap(Description,30),levels=str_wrap(Description,30))),
                                          aes(Count,Description,colour=p.adjust,size=Count)) +
  geom_point() + cols + scale_size(range=c(3,10),limits=c(0,max_count),name="Count") + facet_wrap(~Group) +
  labs(x="Gene count",y=NULL) + theme_bw() +
  theme(strip.background=element_rect(fill="grey90",colour="black"), strip.text=element_text(face="bold"),
        axis.text.y=element_text(size=10,hjust=1), legend.position="right", plot.margin=margin(8,5.5,5.5,5.5))

# Figure 3E
p1 <- make_dotplot(card_df, scale_colour_gradient(high="#E5F5E0", low="#238B45", limits=card_lim, name="Cardiac adj. p"))
p2 <- make_dotplot(epi_df,  scale_colour_gradient(high="#FEEDDE", low="#E6550D", limits=epi_lim, name="Epicardial adj. p")) + guides(size="none")
leg <- plot_grid(get_legend(p1), get_legend(p2), ncol=1)
final_plot <- plot_grid(plot_grid(p1+theme(legend.position="none"), p2+theme(legend.position="none"), nrow=1, align="hv"),
                        leg, nrow=1, rel_widths=c(1,0.35), align="v")
ggsave("GO_dotplots.pdf", final_plot, width=12, height=7, units="in", device=cairo_pdf)
