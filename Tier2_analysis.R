#Microbiome
source('http://bioconductor.org/biocLite.R')
BiocManager::install('phyloseq')
setwd("~/Documents/project_tracking/Bernstein_Kenneth/DO-7511--08--13--2019")
packageVersion('phyloseq')
## [1] '1.26.1'

library("phyloseq")
library(ggplot2)
#Importing biom-format files
biom_otu_tax <- import_biom("16S_table__Analysis_keep__.biom")
# Add sample data to the dataset using merge
bmsd <- import_qiime_sample_data("map_subset.txt")
class(bmsd)
dim(bmsd)
Bushman <- merge_phyloseq(biom_otu_tax, bmsd)
Bushman
#Interacting with the taxonomic ranks
rank_names(Bushman) #"Rank1" "Rank2" "Rank3" "Rank4" "Rank5" "Rank6" "Rank7"
colnames(tax_table(Bushman)) <- c("Kingdom", "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species")
rank_names(Bushman) #"Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"
get_taxa_unique(Bushman, "Family")

OTU = otu_table(Bushman, taxa_are_rows = TRUE)
TAX = tax_table(Bushman)
physeq = phyloseq(OTU, TAX)
physeq
plot_bar(physeq, fill = "Rank2") # show barplot on Family level
plot_heatmap(physeq)
plot_richness(Bushman, x = "Group") + geom_boxplot()
#### filterby relative abundance #--min_count_fraction 0.00001 (0.001%)
filtered <- import_biom("./19samples/16S_table__Analysis_keep__.biom",parseFunction = parse_taxonomy_greengenes)
bmsd <- import_qiime_sample_data("./19samples/example_map__Analysis_keep__.txt")
Bushman_filtered <- merge_phyloseq(filtered, bmsd)


get_taxa_unique(Bushman_filtered, "Family")
OTU_filtered = otu_table(Bushman_filtered, taxa_are_rows = TRUE)
TAX_filtered = tax_table(Bushman_filtered)
#physeq_filtered = phyloseq(OTU_filtered, TAX_filtered)
physeq_filtered  = Bushman_filtered
physeq_filtered = subset_taxa(physeq_filtered, Family != "NA" )
physeq_filtered = filter_taxa(physeq_filtered, function(x) mean(x) > 0.1, TRUE)


physeq_filtered = transform_sample_counts(physeq_filtered, function(x) x / sum(x) )
physeq_filtered = tax_glom(physeq_filtered, taxrank="Family", NArm=FALSE)
physeq_filtered1 = tax_glom(physeq_filtered, taxrank="Genus", NArm=FALSE)
p <- plot_bar(physeq_filtered,fill = "Family")

#p <- plot_bar(physeq_filtered, x= "Group",fill = "Family")
p  + theme(legend.position="bottom", axis.text.x = element_text(angle = 0)) + guides(fill=guide_legend(nrow=8))

plot_richness(Bushman_filtered, x = "Group", measures=c("Shannon", "Chao1"), color = "Group") + geom_boxplot() + theme(axis.line = element_line(colour = "black"),
                                                                                                                       panel.grid.major = element_blank(),
                                                                                                                       panel.grid.minor = element_blank(),
                                                                                                                       panel.border = element_blank(),
                                                                                                                       panel.background = element_blank(),
                                                                                                                       axis.text.x = element_text(angle = 0)) 

plot_ordination(Bushman_filtered, ordinate(Bushman_filtered, "MDS"), color = "Group") + geom_point(size = 5)
plot_ordination(Bushman, ordinate(Bushman, "MDS"), color = "Group") + geom_point(size = 5)

GP.MDS = ordinate(GP100, method = "MDS", distance = "unifrac") #calculate the unweighted-UniFrac distance for each sample pair in the dataset
GP.NMDS = ordinate(Bushman_filtered, "MDS", "bray")  # perform NMDS on bray-curtis distance
GP.NMDS.wUF.ord = ordinate(GP100, "NMDS", "unifrac", weighted = TRUE)  # weighted-UniFrac

# prune to just the top 100 most abundant OTUs across all samples (crude).
GP100 = prune_taxa(names(sort(taxa_sums(Bushman_filtered), TRUE))[1:100], Bushman_filtered)
library(ggrepel)
p <- plot_ordination(Bushman_filtered, GP.NMDS, type = "samples", color = "Group", title = "PCoA") + scale_color_manual(values = c("deepskyblue3", "seagreen3", "coral3"))
p + theme_bw() + geom_text_repel(mapping = aes(label = X.Sample), size = 4, vjust = 1.5) + theme(text = element_text(size = 16)) + geom_point(size = 3)

plot_heatmap(Bushman_filtered, method ="MDS", distance ="bray" ,taxa.label = "Genus") + xlab("Sample")

################################ DEseq #################################
library("DESeq2")
diagdds = phyloseq_to_deseq2(Bushman_filtered, ~ Group)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")
norm_1=counts(diagdds,normalized = T)
write.csv(norm_1,"Abandance_nrom.csv")
#comp <- "At_vs_Wt"
#pattern <- c("At", "Wt")
#data <- diagdds[,grep(paste(pattern, collapse = "|"), colnames(diagdds),value = TRUE)]

res = results(diagdds, cooksCutoff = FALSE,contrast = c("Group", "At", "Wt"))
res1 = results(diagdds, cooksCutoff = FALSE,contrast = c("Group", "Db", "Wt"))
res2 = results(diagdds, cooksCutoff = FALSE,contrast = c("Group", "At", "Db"))
alpha = 0.05
sigtab = res[which(res1$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Bushman_filtered)[rownames(sigtab), ], "matrix"))
head(sigtab)
dim(sigtab)

write.csv(sigtab,"Sig_OTUs_At_vs_Wt_padj0.05.csv")
write.csv(sigtab,"Sig_OTUs_Db_vs_Wt_padj0.05.csv")
write.csv(sigtab,"Sig_OTUs_At_vs_Db_padj0.05.csv")

library("ggplot2")
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
ggplot(sigtab, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=3) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

