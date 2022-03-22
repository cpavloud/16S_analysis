#Tell me where is the working directory (this command its like pwd in UNIX)
getwd()
#Change the working directory (like cd)
setwd("~/data/MiSeq_SOP/")

# for the actual analysis
library(phyloseq); packageVersion("phyloseq")
# for the beautiful plots
library(ggplot2); packageVersion("ggplot2")
#to create beautiful colour palettes for your plots
library(RColorBrewer); packageVersion("RColorBrewer") 
# Define a default theme for ggplot graphics
theme_set(theme_bw()) 


####################### Import and format your data ##########################
##############################################################################

#import the OTU table (or else biotic data)
biotic <- read.csv("seq_table.csv", sep = ",", header=TRUE, row.names = 1)
#import the taxonomy table
taxonomy <- read.csv("seq_Taxonomy_silva.csv", sep = ",", header=FALSE, row.names = 1, na.strings=c("","NA"))
colnames(taxonomy) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#check where there are NA values in the taxonomy data frame
colSums(is.na(taxonomy))

#fill in Phylum column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Phylum[i])){
      taxonomy$Phylum[i]=taxonomy$Domain[i]
    }
  if (sum(is.na(taxonomy$Phylum))==0) {
    break
  }
}

#fill in Class column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Class[i])){
      taxonomy$Class[i]=taxonomy$Phylum[i]
    }
  if (sum(is.na(taxonomy$Class))==0) {
    break
  }
}

#fill in Order column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Order[i])){
      taxonomy$Order[i]=taxonomy$Class[i]
    }
  if (sum(is.na(taxonomy$Order))==0) {
    break
  }
}

#fill in Family column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Family[i])){
      taxonomy$Family[i]=taxonomy$Order[i]
    }
  if (sum(is.na(taxonomy$Family))==0) {
    break
  }
}

#fill in Genus column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Genus[i])){
      taxonomy$Genus[i]=taxonomy$Family[i]
    }
  if (sum(is.na(taxonomy$Genus))==0) {
    break
  }
}

#fill in Species column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Species[i])){
      taxonomy$Species[i]=taxonomy$Genus[i]
    }
  if (sum(is.na(taxonomy$Species))==0) {
    break
  }
}

#check where there are NA values in the taxonomy data frame
colSums(is.na(taxonomy))

#convert the taxonomy data from data frame to matrix
taxonomy_matrix <- as.matrix(taxonomy)

# prepare the object for the phyloseq object
TAX = tax_table(taxonomy_matrix)
head(TAX)

#convert the biotic data from data frame to matrix
biotic_matrix <- as.matrix(biotic)
# prepare the objects for the phyloseq object
OTU = otu_table(biotic_matrix, taxa_are_rows = TRUE)
head(OTU)

#import the metadata of the samples
metadata_physeq <- read.table("mouse.dpw.metadata", header=TRUE, row.names = 1) 
# prepare the objects for the phyloseq object
META = sample_data(metadata_physeq)
head(META)
#class(META)

###########################PHYLOSEQ analysis##################################
##############################################################################

# combine them all to create the phyloseq object
physeq = phyloseq(OTU, TAX, META)
#check what the phyloseq object contains
physeq

#Check if there are ASVs with no counts and how many there are
any(taxa_sums(physeq) == 0)
sum(taxa_sums(physeq) == 0)
#Remove ASVs/OTUs/taxa that are empty, i.e. that have no counts
physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
# Remove samples that are now empty, i.e. that have no counts
physeq <- prune_samples(sample_sums(physeq) > 0, physeq)

# create a bar plot, for the taxon rank you like, e.g Phylum 
# this bar chart is created with the default colour palette of phyloseq
barchart <- plot_bar(physeq, fill = "Phylum")
pdf("barchart_silva.pdf", width = 12)
print(barchart)
dev.off()

#get the data frame from the phyloseq object
pd <- psmelt(physeq)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd[,c("Phylum")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Set2".
# "Set2" is one of the colorblind friendly palettes
# Another example of a colorblind friendly palette is "Dark2"
# If you want, by running the command display.brewer.all(colorblindFriendly = TRUE)
# you can see all the colorblind friendly palettes of the RColorBrewer package.
## Now the "getPalette" variable, we set it in the if statement of the marker gene.
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
PhylaPalette = getPalette(HowManyPhyla)

#and do the actual plotting
barchart_palette_better <- ggplot(pd, aes(x = Sample, y = Abundance, factor(Phylum), fill = factor(Phylum))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = PhylaPalette) +
  labs(fill = "Phylum") +
  guides(fill=guide_legend(ncol=1)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  ylab("Relative abundance") +
  xlab("Samples") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"), 
        legend.title=element_text(size=14), legend.text=element_text(size=12))
ggsave("barchart_palette_better_silva_custom.png", width = 22, height = 15)

# plot the diversity indices with colour coding by e.g. dpw (info included in the metadata) 
richness <- plot_richness(physeq, measures=c("Observed", "Chao1", "ACE"), color="dpw")
pdf("richness_silva.pdf", width = 14)
print(richness)
dev.off()

# create a heatmap, for the Phylum rank 
heatmap <- plot_heatmap(physeq, taxa.label="Phylum")
pdf("heatmap_Phylum_silva.pdf")
print(heatmap)
dev.off()

# create a heatmap, for the Class rank 
heatmap <- plot_heatmap(physeq, taxa.label="Class")
pdf("heatmap_Class_silva.pdf")
print(heatmap)
dev.off()

# create the nMDS plot colour coded by e.g. dpw (info included in the metadata) 
ord.nmds.bray <- ordinate(physeq, method="NMDS", distance="bray")
p1 <- plot_ordination(physeq, ord.nmds.bray, color="dpw", title="Bray NMDS")
pdf("p1_silva.pdf")
print(p1)
dev.off()

#check if the clustering you see in the nMDS plot is 
#statistically significant by running PERMANOVA
metadata_permanova <- as(sample_data(physeq), "data.frame")
permanova.dpw <- adonis(distance(physeq, method="bray") ~ dpw, data = metadata_permanova)
permanova.dpw

###########################EXPORT ASV ABUN TABLE TO EXCEL#####################
##############################################################################

#merges ASVs that have the same taxonomy at a the Phylum level
physeq_merged_Phylum <- tax_glom(physeq, "Phylum")
#transform counts to percentages
ps0 <- transform_sample_counts(physeq_merged_Phylum, function(x) x / sum(x))
plot_bar(ps0, fill="Phylum")
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "ASV")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "ASV")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "ASV")
#save the final phyto otu table with the taxonomies
write.csv(OTU_TAX_merged, "OTU_TAX_merged_phylum.csv")

############## PHYLOSEQ analysis for specific taxa############################
##############################################################################

# Create subset of the phyloseq object, e.g. 
physeq.Cyanobacteria = subset_taxa(physeq, Phylum=="Cyanobacteria")

