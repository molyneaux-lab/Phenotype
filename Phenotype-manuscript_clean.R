#Code associated with the publication:
# Investigating the lasting effects of SARS-CoV-2 infection and the lung microbiota: 
#   No persistent microbial alterations in recovered 
#   COVID-19 patients with persistent radiological or respiratory abnormalities.
# Author: Nancy Teng
# GitHub: https://github.com/molyneaux-lab/Phenotype
setwd("~/Documents/GitHub/PHENOTYPE/")

path <- getwd()

library("phyloseq")
library("dplyr")
library("tidyverse")
library("vegan")
library("reshape2")
library("qiime2R")
library("decontam")
#BiocManager::install("microbiome")
library("microbiome")
library(readxl)
library(ggrepel)
library(ggplotify)

#### READ ME ####
# This script creates the figures and cleaned microbiota sequencing data produced in 
# the manuscript listed above. 
# BAL samples were collected from patients at various time points and underwent 16S rRNA gene amplicon sequencing.
# This script is shared to promote data transparency and to be used for methodological reproducibility. 

# The script shares the filtering parameters and removal of contaminants, if applicable.

Palette_color = c("#BC3C29FF", "#0072B5FF")
Palette_fill = c("#BC3C29FF", "#0072B5FF")

Palette_10 = c("#003f5a", # deep navy (blue-teal)
               "#0072B5", # blue
               "#6F99AD", # slate blue
               "#6ea6a4", # muted turquoise
               "#20854E", # green-teal
               "#8ab184", # sage
               "#C4644A", # muted pink
               "#BC3C29", # red-brown
               "#de6600", # orange
               "#E18727", # gold
               "#FFDC91", # light gold
               "#fec682" # peach
)
##### Qiime2 clean up #####
# --- Load the metadata file to index the samples to be used ---
#metadata <- readxl::read_xlsx("Clinical_metadata_summarised.xlsx")
#metadata$`PHENOTYPE ID` <- paste0("PH", metadata$`PHENOTYPE ID`)
#metadata$`PHENOTYPE ID` <- paste0(metadata$`PHENOTYPE ID`, ".BAL")
#metadata <- dplyr::rename("sample-id"="PHENOTYPE ID", metadata)

#manifest <- read_tsv("Manifest-map.txt")
#metadata <- full_join(manifest, metadata)
#metadata <- metadata %>%
#  dplyr::select(!2:8)
#sampleIDs <- metadata$`sample-id`

# --- Upload the ddPCR samples ---
#ddPCR <- read_csv("16S_ddPCR_ms.csv")
#ddPCR <- ddPCR %>%
#  filter(`sample-id` %in% sampleIDs) %>%
#  dplyr::select(`sample-id`, Diagnosis, ddPCR)

#setdiff(metadata$`sample-id`, ddPCR$`sample-id`) # 6 samples missing from original sheet
#intersect(metadata$`sample-id`, ddPCR$`sample-id`)

#metadata <- full_join(metadata, ddPCR, by = "sample-id")
#metadata <- metadata %>%
#  dplyr::select(!Diagnosis.y) %>%
#  drop_na(ddPCR)
#metadata <- dplyr::rename("Diagnosis"="Diagnosis.x", metadata)
#metadata$Diagnosis <- gsub("Negatve", "Negative", metadata$Diagnosis)

#metadata$Diagnosis <- gsub("Controls" , "Non-Fibrotic Controls", metadata$Diagnosis)
#metadata$Diagnosis <- gsub("COVID", "COVID-19", metadata$Diagnosis)

write_csv(metadata, "metadata.csv")

#### Loading phyloseq object ####
physeq<-qza_to_phyloseq(
  features="filtered-metadata-table.qza",
  tree="merged-midrooted-tree.qza",
  taxonomy= "taxonomy-merged-final.qza")
physeq 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 372 taxa and 76 samples ]
#tax_table()   Taxonomy Table:    [ 372 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 372 tips and 369 internal nodes ]


mapping = read_csv("metadata.csv")
mapping <- column_to_rownames(mapping, var = "sample-id")

gplots::venn(list(mapping=rownames(mapping), physeq=sample_names(physeq)))
setdiff(sample_names(physeq), rownames(mapping))  

phyloseq::sample_data(physeq)<-mapping
view(sample_data(physeq))

physeq

# --- Looking at library sizes ---
ps<-physeq
summarize_phyloseq(ps)
df <- as.data.frame(sample_data(ps)) 
df$LibrarySize <- sample_sums(ps) # similar to rowSums/colSums but automated
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, colour=Diagnosis)) + geom_point() + 
  ylim(0,70000)

# --- Look at contaminants in low biomass samples ---
# we're using frequency method as we have ddPCR values
sample_data(ps)$is.neg <- sample_data(ps)$Diagnosis=="Negative control" #Only considered reagent control as negative
contamdf.freq <- isContaminant(ps, method="combined", conc="ddPCR", neg="is.neg", threshold=0.1) 
# Gets data frame with a list of potential contaminants

table(contamdf.freq$contaminant)
head(which(contamdf.freq$contaminant)) # not the highest abundant taxa
# [1]  17  65 157 189 196 197

# To make a presence/absence plot
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
# Transforms data to absence or presence 
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Diagnosis =="Negative control", ps.pa)
ps.pa.pos <- prune_samples(!sample_data(ps.pa)$Diagnosis == "Control", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.freq$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") 

plot_frequency(ps, taxa_names(ps)[c(which(contamdf.freq$contaminant))], conc="ddPCR") + 
  xlab("DNA Concentration (ddPCR)")

tax <- as(tax_table(ps), "matrix")
contaminants<-tax[which(contamdf.freq$contaminant),]

#--- Look at contaminats by prevalence ---
sample_data(ps)$is.neg <- sample_data(ps)$Diagnosis=="Negative control"
#Only considered reagent control as negative
view(sample_data(ps))
contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.1) 
# Gets data frame with a list of potential contaminants

table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant)) # not the highest abundant taxa
#[1]  17 196 197 199 217 247

# Comparing prevalence with frequency, 17 is identified as contaminant in both. 
# To make a presence/absence plot
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
# Transforms data to absence or presence 
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Diagnosis =="Negative control", ps.pa)
ps.pa.pos <- prune_samples(!sample_data(ps.pa)$Diagnosis == "Negative control", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") 
# Better split between true and false contaminant

plot_frequency(ps, taxa_names(ps)[c(which(contamdf.prev$contaminant))], conc="ddPCR") + 
  xlab("DNA Concentration (ddPCR)")

tax <- as(tax_table(ps), "matrix")
contaminants<-tax[which(contamdf.prev$contaminant),]
tax_table <- rownames_to_column(data.frame(tax), var="OTU")

#### Decide what are contaminants ####
ps.contam <- prune_taxa(contamdf.prev$contaminant, ps)
plot_bar(ps.contam, facet_grid =Diagnosis~.)
# Surprisingly high number of contaminants in COVID
ps.contam <- prune_taxa(contamdf.freq$contaminant, ps)
plot_bar(ps.contam, facet_grid =Diagnosis~.)
# With frequency method there are a few samples that are clearly contaminated. 

# Identify ASVs above 1000 reads in the list of contaminants. These are big influencers. 
filter <- phyloseq::genefilter_sample(ps.contam, filterfun_sample(function(x) x >= 1000))
ps.contam.1k <- prune_taxa(filter, ps.contam)
otu_table<-as.data.frame(ps.contam@otu_table)
tax_table<-as.data.frame(ps.contam@tax_table) 

library(tibble)
library(dplyr)

# --- Unfiltered data ---
# In the spirit of being transparent we decided not to remove contaminants from our dataset,
# this could influence the relative abundances and skew composition of the microbiota community.
# We instead opt for remove samples that are clearly outliers i.e., most likely fully contaminated. 

summarize_phyloseq(ps)
ps

# Extract taxa information to get the number of unique Phyla
tax <- as(tax_table(ps), "matrix")
tax_df <- as.data.frame(tax)
filterPhyla = unique(tax_df$Phylum)
filterPhyla <- na.omit(filterPhyla)

ps1 = subset_taxa(ps, !(!Phylum %in% filterPhyla))
ps1 # Only keep the Phyla in filterPhyla in the filtered reads dataset
summarize_phyloseq(ps1)

# Check if there are unique Phyla names
unique(as.data.frame(as(tax_table(ps1), "matrix"))$Phylum)

#Filter at 0.005% 
minTotRelAbun = 0.00005
x = taxa_sums(ps1)
keepTaxa = which((x / sum(x)) > minTotRelAbun)
prunedSet = prune_taxa(names(keepTaxa), ps1)
prunedSet # Only keeping phyla that have a minimum relative abundance of 0.005%. 
unique(as.data.frame(as(tax_table(prunedSet), "matrix"))$Phylum)

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 157 taxa and 70 samples ]
#sample_data() Sample Data:       [ 70 samples by 17 sample variables ]
#tax_table()   Taxonomy Table:    [ 157 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 157 tips and 155 internal nodes ]

thrownTaxa = which((x / sum(x)) < minTotRelAbun)
trashSet = prune_taxa(names(thrownTaxa), ps1)
view(tax_table(trashSet)) # just cross check that those that are removed are indeed in low abundance/not normal

sum(taxa_sums(prunedSet)==0) # should be 0

summarize_phyloseq(prunedSet)

#### Set function to normalise samples ####
## Normalising samples i.e., relative abundances
normalizeSample = function(x) {
  x/sum(x)
}

Controls_relative = transformSampleCounts(prunedSet, normalizeSample)
otu_table(Controls_relative)
OTU1 = as(otu_table(Controls_relative), "matrix")
OTUdf = as.data.frame(OTU1)

TAXdf = as(tax_table(Controls_relative), "matrix")
TAXdf = as.data.frame(TAXdf)

Controls_Phylum <- aggregate_taxa(Controls_relative, 'Phylum') #7 phyla most likely a contaminant
Phylum_df<-as.data.frame(Controls_Phylum@otu_table)

Controls_Family <- aggregate_taxa(Controls_relative, 'Family')
Family_df<-as.data.frame(Controls_Family@otu_table)

Controls_Genus <- aggregate_taxa(Controls_relative, 'Genus')
Genus_df<-as.data.frame(Controls_Genus@otu_table)

### Metadata joining ###
Phylum_df <- data.frame(t(Phylum_df))
Phylum_df <- rownames_to_column(Phylum_df, "sample-id")
Phylum_df <- inner_join(metadata, Phylum_df, by = "sample-id")

Family_df <- data.frame(t(Family_df))
Family_df <- rownames_to_column(Family_df, "sample-id")
Family_df <- inner_join(metadata, Family_df, by = "sample-id")

Genus_df <- data.frame(t(Genus_df))
Genus_df <- rownames_to_column(Genus_df, "sample-id")
Genus_df <- inner_join(metadata, Genus_df, by = "sample-id")

### --- Create ps object without negative controls --- ###
physeq

otu_genus = aggregate_taxa(physeq, "Genus")
otu_genus = as.data.frame(t(otu_table(otu_genus)))
otu_genus = rownames_to_column(otu_genus, "sample-id")

physeq_filtered <- subset_samples(physeq, c(!Diagnosis=="Negative control"))

## --- Using the aggregated taxa (%) data set we identify common contaminants --- ##
### - we use the negative controls as our reference group ---
##### Common contaminants #####
df <- read_csv("Genus-normalised-metadata.csv") %>%
  dplyr::select(!Unknown)

negative <- df %>%
  filter(Diagnosis == "Negative control")

metadata <- df %>%
  dplyr::select(`sample-id`, Diagnosis, ddPCR)

# --- Negative control plot #
negative_genera <- negative %>%
  dplyr::select(Actinomyces:ncol(negative))

negative_genera <- negative_genera[,order(colSums(negative_genera),decreasing=TRUE)]
negative_genera <- negative_genera/rowSums(negative_genera)*100 # change to percentages
rowSums(negative_genera)

# only select genera where total abundance exceeds 0.01%
negative_genera <- names(negative_genera)[colSums(negative_genera) > (sum(negative_genera)*0.01)]

df_negative <- df[, negative_genera]
df_negative <- df_negative/rowSums(df_negative)*100
df_negative <- cbind(metadata, df_negative)
df_negative <- df_negative %>%
  filter(rowSums(df_negative[,4:ncol(df_negative)])>0)
df_dropped <- df_negative %>%
  filter(!rowSums(df_negative[,4:ncol(df_negative)])>0)

df_negative_long <- reshape2::melt(df_negative, id.vars = c("sample-id",
                                                  "Diagnosis",
                                                  "ddPCR"),
                         variable.name = "Genus")

df_negative_long_summarised <- df_negative_long %>%
  group_by(Diagnosis, Genus) %>%
  summarise(mean_value=mean(value),
            sd=sd(value))

df_negative_long_summarised$Genus <- as.factor(df_negative_long_summarised$Genus)
genera_contaminant <- df_negative_long_summarised %>% # only keep taxa at a 5% more abundance, more likely to be true contaminant
  filter(mean_value > 5 & Diagnosis == "Negative control")
genera_contaminant <- as.character(genera_contaminant$Genus)

df_negative_long_summarised <- df_negative_long_summarised %>%
  filter(Genus %in% genera_contaminant)

unique(genera_contaminant)

ggplot(df_negative_long_summarised,
       aes(x = mean_value,
           y = Genus,
           fill = Genus)) +
  #Set the width of the bars in the plot
  geom_bar(stat = "identity",
           width = 0.7) +
  geom_errorbar(aes(xmax=mean_value+sd, xmin=mean_value-sd),
                width=.2,
                position=position_dodge(.9)) +
  facet_grid(. ~ Diagnosis,
             scales = "free", 
             drop = TRUE)+ 
  scale_x_continuous(expand = c(0, 0),
                     limits = c(-20, 100))+
  scale_fill_manual(values = Palette_10) +
  labs(title="Relative abundances of likely contaminant taxa present at 5% abundance",
       subtitle = "Negative control specimens: Reagent controls.",
       x = "Mean relative abundance (%)",
       y = "Genus") + theme_classic()+
  theme(#Set the title font size
    plot.title = element_text(size=20),
    plot.subtitle = element_text(size = 18),
    legend.position = "right",
    legend.title = element_text(size=18),
    legend.text = element_text(size=16,
                               face = "italic"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=90,
                               hjust=1,
                               vjust=0.5,
                               size=14),
    axis.title.x = element_text(size=16),
    axis.title = element_text(size=16),
    axis.text.y = element_text(size=14,
                               face = "italic"),
    axis.line = element_line(size = 0.5,
                             linetype = "solid",
                             colour = "black"),
    aspect.ratio = 1.7
  )

# --- Using ASV we will assess beta diversity ---
## --- this plot also identifies outlier samples, clearly contaminated --- 
x = taxa_sums(physeq_filtered)
# Keep taxa seen at least twice in more than 1% of samples.
filteredset = phyloseq::filter_taxa(physeq_filtered, function(x) sum(x > 2) > (0.1*length(x)), TRUE)
filteredset

minTotRelAbun = 0.00005
x = taxa_sums(filteredset)
keepTaxa = which((x / sum(x)) > minTotRelAbun)
prunedSet_beta = prune_taxa(names(keepTaxa), filteredset)

normalizeSample = function(x) {
  x/sum(x)}
tax <- as(tax_table(physeq), "matrix")

Controls_relative_beta = transformSampleCounts(prunedSet_beta, normalizeSample)
Controls_relative_beta <- aggregate_taxa(Controls_relative_beta, "Genus")
otu_table <- as.data.frame(Controls_relative_beta@otu_table)
tax <- as.data.frame(tax)
Controls_relative_beta_genus <- left_join(rownames_to_column(otu_table), (rownames_to_column(tax)))

ordu = ordinate(filteredset, "PCoA", "unifrac", weighted=TRUE)
#Temporary fix to colour bug is to add a "dummy variable" in sample_data:
#sample_data(physeq)[ , 2] <- sample_data(physeq)[ ,1]
p = plot_ordination(filteredset, ordu, color="Diagnosis") + geom_point(size = 3) +
  ggtitle("PCoA on weighted-UniFrac distance") + (scale_colour_brewer(type="qual", palette="Set1"))
print(p)

myplotdiagnosis <- p
myplotdiagnosis + theme_bw() + stat_ellipse() + theme_classic() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggrepel::geom_text_repel(label=rownames(p$data),colour="black", size=3, max.overlaps = 25)+
  scale_colour_manual(values=c(Palette_color))

ggsave("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Manuscript/Figures/PCOA-contaminant.svg",
       width = 12, height = 10)

###### --- Create the top ten taxa barplot ---- #####
## --- helps to identify samples that are highly contaminated ---
df <- read_csv("Genus-normalised-metadata.csv")
## df_filtered includes samples that are most likely contaminated by environmental controls. 
df_filtered <- df %>% filter(Diagnosis=="COVID-19" | Diagnosis=="Non-Fibrotic Controls")

abund_table <- df_filtered %>%
  dplyr::select(1,Actinomyces:ncol(df_filtered))
abund_table <- column_to_rownames(abund_table, var="sample-id")
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100 # change to percentages

meta_table <- df_filtered %>%
  dplyr::select(`sample-id`, Diagnosis)

top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(df_filtered)[2]
taxa_list2 <- colnames(other)[11:N]
N <- length(taxa_list2)
Others <- data.frame(other[,colnames(other) %in% taxa_list2])
Others
# Sum all the other columns into one.
Others <- rowSums(Others)
Others <- as.data.frame(Others)
Others

# Combine the top 10 genus and the other column and check the rows all still add up to 100.
top_other <- cbind(top, Others)
top_other
top_other_rowsums <- rowSums(top_other)
top_other_rowsums <- as.data.frame(top_other_rowsums)
top_other_rowsums

sample_data <- cbind(meta_table, top_other)
sample_data <- sample_data %>%
  dplyr::arrange(Streptococcus, #desc() # Using desc() If you want to arrange in descending order.
  )

sample_data_long <- melt(sample_data, id.vars = c("sample-id", "Diagnosis"),
                         variable.name = "Genus")
sample_data_long
taxa_list

ggplot(sample_data_long,
       aes(x = `sample-id`,
           y = value,
           fill = Genus)) +
  geom_bar(stat = "identity",
           width = 0.7) +
  scale_fill_manual(values = Palette_10) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100.1)) +
  facet_grid(.~Diagnosis,
             scales = "free_x",
             drop=TRUE)+
  labs(title="Ten most abundant genera in COVID patients \nand non-fibrotic controls",
       x = "Sample ID",
       y = "Relative abundance (%)") +
  theme(#Set the title font size
    plot.title = element_text(size=10),
    plot.subtitle = element_text(size = 8),
    axis.text.x = element_text(angle=90,
                               hjust=1,
                               vjust=0.5,
                               size=8),
    axis.title.y = element_text(size = 12)) +
   bbplot::bbc_style() 

# --- Create the plots for the manuscript ---
#### --- Top ten genera, contaminants removed ---#####
df <- read_csv("Genus-normalised-metadata.csv")

df_filtered <- df %>% filter(Diagnosis=="COVID-19" | Diagnosis=="Non-Fibrotic Controls") %>%
  filter(!`sample-id` == "PH1039.BAL" &
           !`sample-id` =="BRU.03853"&
           !`sample-id` =="PH1046.BAL" &
           !`sample-id`=="PH1024.BAL" &
           !`sample-id`=="BRU.03853" &
           !`sample-id` == "BRU.03815" &
           !`sample-id`=="BRU.1032")

abund_table <- df_filtered %>%
  dplyr::select(1,Actinomyces:ncol(df_filtered))
abund_table <- column_to_rownames(abund_table, var="sample-id")
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100 # change to percentages

meta_table <- df_filtered %>%
  dplyr::select(`sample-id`, Diagnosis)

# Make a dataframe "top" that contains the top 10 most abundant genera. You can alter these as you want.
top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(df_filtered)[2]
taxa_list2 <- colnames(other)[11:N]
N <- length(taxa_list2)
Others <- data.frame(other[,colnames(other) %in% taxa_list2])
Others
# Sum all the other columns into one.
Others <- rowSums(Others)
Others <- as.data.frame(Others)
Others

# Combine the top 10 genus and the other column and check the rows all still add up to 100.
top_other <- cbind(top, Others)
top_other
top_other_rowsums <- rowSums(top_other)
top_other_rowsums <- as.data.frame(top_other_rowsums)
top_other_rowsums

sample_data <- cbind(meta_table, top_other)
sample_data <- sample_data %>%
  dplyr::arrange(Streptococcus, #desc() # Using desc() If you want to arrange in descending order.
  )

sample_data_long <- melt(sample_data, id.vars = c("sample-id", "Diagnosis"),
                         variable.name = "Genus")
sample_data_long
taxa_list

ggplot(sample_data_long,
       aes(x = `sample-id`,
           y = value,
           fill = Genus)) +
  geom_bar(stat = "identity",
           width = 0.7) +
  scale_fill_manual(values = Palette_10) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100.1)) +
  facet_grid(.~Diagnosis,
             scales = "free_x",
             drop=TRUE)+
  labs(title="Ten most abundant genera in COVID patients \nand non-fibrotic disease controls",
       subtitle = "Contaminant samples removed",
       x = "Sample ID",
       y = "Relative abundance (%)") +
  theme(#Set the title font size
    plot.title = element_text(size=10),
    plot.subtitle = element_text(size = 6),
    axis.text.x = element_text(angle=90,
                               hjust=1,
                               vjust=0.5,
                               size=8),
    axis.title.y = element_text(size = 12)) +
  bbplot::bbc_style() 

ggsave("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Manuscript/Figures/Barplot-contaminants.svg",
       width = 15, height = 10)

#### --- Permanova analysis --- #####
permanova_diagnosis = adonis2(abund_table ~ Diagnosis,
        data = meta_table, permutations = 999, method = "bray")
permanova_diagnosis = permanova_diagnosis$`Pr(>F)`

meta_table <- df_filtered %>%
  dplyr::select(`sample-id`: ddPCR) %>% 
  filter(Diagnosis == "COVID-19") %>% drop_na()
abund_table <- df_filtered %>%
  filter(`sample-id` %in% meta_table$`sample-id`) %>%
  dplyr::select(Actinomyces:ncol(df_filtered))

adonis2(abund_table ~ Age, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ `CT: Opacified Lung (%)`, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ `BAL Neutrophil (%)`, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ `Severity classification`, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ Microbiology, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ `Blood Neutrophils`, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ `Blood lymphocytes`, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ `Blood CRP`, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ `Blood Ferritin`, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ `Blood Fibrinogen`, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ FEV1, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ `%FEV1`, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ FVC, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ `%FVC`, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ TLCO, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ ddPCR, data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ Sex, data = meta_table, permutations = 999, method = "bray")

dist_matrix <- vegdist(abund_table, method = "bray")
dispersion_model = betadisper(dist_matrix, group = meta_table$Diagnosis)
disperm = permutest(dispersion_model, permutations = 999)
p_value = disperm$tab$`Pr(>F)`

  plot(dispersion_model)

## Dispersion permanova plot ##
points <- data.frame(dispersion_model$vectors[,1:2], 
                     Diagnosis = meta_table$Diagnosis)

centroids <- data.frame(dispersion_model$centroids[,1:2], 
                        Diagnosis = rownames(dispersion_model$centroids))

ggplot(points, aes(x = PCoA1, y = PCoA2, color = Diagnosis)) +
  geom_point(alpha = 0.6, size = 3) +
  stat_ellipse(type = "t", level = 0.95) + # Add 95% confidence ellipses
  geom_point(data = centroids, aes(fill = Diagnosis), size = 3, shape = 13) +
  theme_bw() +
  labs(title = "PCoA of Bray-Curtis Dissimilarity",
       subtitle = glue("PERMDISP p={round(p_value[1],3)}")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=12),
        axis.text = element_text(size=12)) + 
  scale_fill_manual(values = Palette_fill) + 
  scale_color_manual(values = Palette_fill)

ggsave("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Manuscript/Figures/pdisp.svg",
       width = 12, height = 12)

#### ddPCR - bacterial burden ####
# --- Analysis by ddPCR i.e., bacterial burden ---
df <- read_csv("metadata.csv")
meta_table <- df %>%
  dplyr::select(`sample-id`:ddPCR)

meta_ddPCR <- meta_table %>% filter(!Diagnosis == "Negative control") %>% 
  dplyr::select(`sample-id`, Diagnosis, Microbiology, ddPCR, `Severity classification`, Sex)

bacterial_burden_stats <- meta_ddPCR %>%
  group_by(Diagnosis) %>%
  mutate(median_burden = median(ddPCR),
         mean_burden = mean(ddPCR),
         q1 = quantile(ddPCR, 0.25),  # 1st quartile
         q3 = quantile(ddPCR, 0.75), # 3rd quartile
         norm_test = shapiro.test(ddPCR)$p.value) %>% # Get p-value from shapiro.test
  filter(!Diagnosis=="Negative control")

wilcox.test(ddPCR ~ Diagnosis, data = meta_ddPCR) 

p <- ggplot(bacterial_burden_stats, aes(x=Diagnosis, y=ddPCR, fill=Diagnosis)) +
  geom_dotplot(binaxis="y", stackdir = "center", binwidth = 0.1, position="dodge") +
  geom_errorbar(aes(x=Diagnosis, ymin=q1, ymax=q3), width=0.3, color='black', linewidth=1) +
  scale_y_log10() + theme_classic() + 
  labs(title="Bacterial burden across diagnosis groups",
       y="Bacteridal burden (16S rRNA gene/mL of BAL)",
       x="") +
  scale_fill_manual(values=c(Palette_fill)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
p

ggsave("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Manuscript/Figures/ddPCR.svg",
       width = 12, height = 12)

#### Severity Classification ####
# --- By severity classification ---
COVID_ddPCR <- meta_ddPCR %>%
  filter(Diagnosis == "COVID-19")
COVID_ddPCR$`Severity classification` <- factor(COVID_ddPCR$`Severity classification`,
                                                levels=c("Mild", "Moderate", "Severe"))
kruskal.test(ddPCR ~ `Severity classification`, data = COVID_ddPCR) # p=0.8
COVID_ddPCR_stats <- COVID_ddPCR %>%
  group_by(`Severity classification`) %>%
  mutate(median = median(ddPCR),
         q1 = quantile(ddPCR, 0.25), 
         q3 = quantile(ddPCR, 0.75)) 

Palette <- c("Mild"="#588300",
             "Moderate" = "orange",
             "Severe" = "firebrick")
p <- ggplot(COVID_ddPCR_stats, aes(x=`Severity classification`, y=ddPCR, fill=`Severity classification`)) +
  geom_dotplot(binaxis="y", stackdir = "center", binwidth = 0.1, position="dodge") +
  geom_errorbar(aes(x=`Severity classification`, ymin=q1, ymax=q3), width=0.3, color='black', linewidth=1) +
  scale_y_log10() + theme_classic() + 
  labs(title="Bacterial burden across severity groups",
       y="Bacteridal burden (16S rRNA gene/mL of BAL)",
       x="") +
  scale_fill_manual(values=c(Palette)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
p

ggsave("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Manuscript/Figures/ddPCR-classificaiton.svg",
       height = 12, width = 12)
#### Ordered bar plot ####
# --- Ordered bar plot ---
df <- read_csv("Genus-normalised-metadata.csv") 
df <- df %>%
  dplyr::filter(!Diagnosis == "Negative control") %>%
  dplyr::filter(!`sample-id` == "PH1039.BAL" &
           !`sample-id` =="BRU.03853"&
           !`sample-id` =="PH1046.BAL" &
           !`sample-id`=="PH1024.BAL" &
           !`sample-id`=="BRU.03853" &
           !`sample-id` == "BRU.03815" &
           !`sample-id`=="BRU.1032")

meta_table <- df %>%
  dplyr::select(`sample-id`:ddPCR)

abund_table <- df %>%
  dplyr::select(`sample-id`, Actinomyces:ncol(df))
abund_table <- column_to_rownames(abund_table, var = "sample-id")
abund_table <- abund_table/rowSums(abund_table)*100
rowSums(abund_table)

top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

df_genus_top <- cbind(meta_table, top)

library(FSA)
# Assuming you have a data frame df_genus_top, a grouping variable Diagnosis, and a list of taxa taxa_list
for (genus in taxa_list) {
  # Perform the Kruskal-Wallis test
  test_result <- kruskal.test(as.formula(paste(genus, "~ Diagnosis")), data = df_genus_top)
  
  if (test_result$p.value < 0.05) {
    # Print the Kruskal-Wallis result
    print(paste("Kruskal-Wallis test result for", genus, ":"))
    print(test_result)
    
    # Perform the Dunn test
    dunn_result <- dunnTest(df_genus_top[[genus]], df_genus_top$Diagnosis, method = "bonferroni")
    
    # Print the Dunn test result
    print(paste("Dunn test result for", genus, ":"))
    print(dunn_result)
    
    # Create a boxplot for this genus
    boxplot(as.formula(paste(genus, "~ Diagnosis")), data = df_genus_top,
            main = paste("Boxplot for", genus),
            xlab = "Diagnosis",
            ylab = "Relative abundance (%)")
  }
}

dunn_result

# Actinomyces p=0.017 
# Neisseria p=0.012
# Haemophilus p=0.03
# Rothia p=0.02
# Gemella p=0.028

df_genus_top <- df_genus_top %>%
  dplyr::select(`sample-id`, Diagnosis, Streptococcus:Granulicatella)
df_long <- reshape2::melt(df_genus_top, id.vars = "Diagnosis", 
                measure.vars = c("Streptococcus", "Prevotella", 
                                 "Veillonella", "Actinomyces", 
                                 "Neisseria", "Haemophilus", 
                                 "Sphingomonas", "Rothia", 
                                 "Gemella", "Granulicatella"), 
                variable.name = "Genus",
                factorsAsStrings = TRUE, na.rm = TRUE)

df_long_summarised <- df_long %>%
  group_by(Diagnosis, Genus) %>%
  summarise(mean_value=mean(value),
            sd=sd(value),
            median=median(value),
            q1 = quantile(value, 0.25),  # 1st quartile
            q3 = quantile(value, 0.75)) # 3rd quartile

dat_text <- data.frame(
  label=c("*", ""),
  Diagnosis=c("Non-Fibrotic Controls", "COVID-19"),
  x = c(4), 
  y = c(14)
)
dat_text_2 <- data.frame(
  label=c("*", ""),
  Diagnosis=c("Non-Fibrotic Controls", "COVID-19"),
  x = c(5), 
  y = c(14)
)
dat_text_3 <- data.frame(
  label=c("*", ""),
  Diagnosis=c("Non-Fibrotic Controls", "COVID-19"),
  x = c(6), 
  y = c(14)
)
dat_text_4 <- data.frame(
  label=c("*", ""),
  Diagnosis=c("Non-Fibrotic Controls", "COVID-19"),
  x = c(8), 
  y = c(14)
)
dat_text_5 <- data.frame(
  label=c("*", ""),
  Diagnosis=c("Non-Fibrotic Controls", "COVID-19"),
  x = c(9), 
  y = c(14)
)

p <- ggplot(df_long_summarised, aes(x=Genus, y=median)) + 
  geom_bar(aes(y = median, x = Genus, fill = Genus),
           stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=(q1), ymax=(q3)), width=0.3, color='black', linewidth=0.5)+
  labs(x = "Genus", y = "Relative abundance (%)",
       title= "Top ten genera of the PHENOTYPE cohort",
       subtitle = "*p<0.05 Wilcoxon rank sum") +
  guides(fill = guide_legend(title = "Genus")) +
  facet_grid(Diagnosis~.) +bbplot::bbc_style()

df_long_summarised$Diagnosis <- ordered(df_long_summarised$Diagnosis, levels = c("Non-Fibrotic Controls", "COVID-19"))

levels(df_long_summarised$Diagnosis)

p <- p + scale_fill_manual(values = Palette_10)+
  geom_text(data=dat_text, mapping= aes(x=x, y=y, label=label), size = 6) +
  geom_text(data=dat_text_2, mapping= aes(x=x, y=y, label=label), size = 6) +
  geom_text(data=dat_text_3, mapping= aes(x=x, y=y, label=label), size = 6) +
  geom_text(data=dat_text_4, mapping= aes(x=x, y=y, label=label), size = 6) +
  geom_text(data=dat_text_5, mapping= aes(x=x, y=y, label=label), size = 6) +
  theme(axis.text.x = element_blank(),
        legend.text = element_text(face = "italic"),
        axis.title.y = element_text(),
        strip.text = element_text(size=18),
        plot.subtitle = element_text(size=14,
                                     face = "italic"))
p  

ggsave("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Manuscript/Figures/Topten.svg",
       width = 12, height = 10)

#### ISA ####
ISA <- multipatt(x = abund_table,
                 cluster = meta_table$Diagnosis, 
                 duleg= T)
summary(ISA)

#### --- PCoA Plot --- ####
# PCoA plot- Diagnosis
abund_table_filtered_robust.aitchison = vegdist(abund_table, "bray")
RA_pcoa = cmdscale(abund_table_filtered_robust.aitchison, k=2, eig=T)  

RA_pcoa_eig <- RA_pcoa$eig
RA_total_variance <- sum(RA_pcoa_eig[RA_pcoa_eig > 0])

percentage_explained_pco1 <- (RA_pcoa_eig[1] / RA_total_variance) * 100
percentage_explained_pco2 <- (RA_pcoa_eig[2] / RA_total_variance) * 100

RA_pcoa_coord = RA_pcoa$points

colnames(RA_pcoa_coord) = c("PCoA1", "PCoA2")

plot.data <- cbind(meta_table, RA_pcoa_coord)

PcoA_main = ggplot(data = plot.data, aes(x = PCoA1, y = PCoA2)) + 
  geom_point(aes(colour = Diagnosis, shape = Diagnosis, fill = Diagnosis),
             size = 2) + 
  stat_ellipse(aes(colour = Diagnosis)) +
  # scale_"" is used to design the plot
  scale_fill_manual(values = Palette_fill) + 
  scale_colour_manual(values = Palette_color) + 
  scale_shape_manual(values = c(21,21,21)) +
  labs(title = "",
       subtitle = glue("Permanova: p={round(permanova_diagnosis[1], 5)}"),
       x = paste("PCoA axis 1(", round(percentage_explained_pco1, digits = 2), "%)", sep = ""),
       y = paste("PCoA axis 2 (", round(percentage_explained_pco2, digits = 2), "%)", sep = "")) +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18),
    legend.position = c(0.2,0.1), legend.title = element_text(size=12, face = "bold", colour = "firebrick"),
    legend.text = element_text(size=10, face = "italic"),
    legend.background = element_blank(), legend.key = element_blank(),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=90, hjust=0, vjust=0.5,size=14),
    axis.title = element_text(size=16),
    axis.title.y = element_text(size=14), 
    axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"))

pcoa1_bar = ggplot(plot.data) + 
  geom_boxplot(aes(x = Diagnosis, y=PCoA1,
                   fill = Diagnosis),
               show.legend = F) + coord_flip() +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = Palette_fill)

pcoa2_bar = ggplot(plot.data) + 
  geom_boxplot(aes(x = Diagnosis, y=PCoA2,
                   fill = Diagnosis),
               show.legend = F) +
  scale_y_continuous(expand=c(0,0.001)) + 
  labs(x=NULL, y=NULL) + theme_classic() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values = Palette_fill)

PcoA_main = PcoA_main %>%
  aplot::insert_top(pcoa1_bar, height = 0.1) %>%
  aplot::insert_right(pcoa2_bar, width=0.1) %>%
  as.ggplot()

PcoA_main

ggsave("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Manuscript/Figures/PcoA-plot.svg",
       width = 12, height = 10)


#### Alpha diversity ####
##### Species number #####
species_number <- vegan::specnumber(abund_table)
species_number <- as.data.frame(species_number)
species_number <- cbind(meta_table, species_number)
wilcox = wilcox.test(species_number ~ Diagnosis, species_number)
p_value = wilcox$p.value

ggplot(species_number, aes(x = Diagnosis, y = species_number),
                fill = Diagnosis, shape = Diagnosis) +
  geom_boxplot(aes(fill = Diagnosis, shape = Diagnosis, colour = Diagnosis)) + 
  geom_point(position = "identity", aes(fill = Diagnosis, shape = Diagnosis, colour = Diagnosis), 
             size = 2, stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21, 21)) +
  scale_color_manual(values=c(Palette_color)) +
  scale_fill_manual(values=c("white", "white", "white", "white")) +
  labs(#title="Progressive IPF vs. Stable IPF", 
    subtitle = glue("Wilcoxon ranked sum p={round(p_value,2)}"),
    x = "",
    y = "Genus number") +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18),
    legend.position = "none", panel.background = element_blank(),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=14),
    axis.title.x = element_blank(),
    axis.title = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
    aspect.ratio = 1)

ggsave("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Manuscript/Figures/SpecNum.svg",
       width = 10, height = 10)

##### Shannon diversity #####
shannon_diversity <- vegan::diversity(abund_table)
shannon_diversity <- as.data.frame(shannon_diversity)
shannon_diversity <- cbind(meta_table, shannon_diversity)
wilcox = wilcox.test(shannon_diversity ~ Diagnosis, shannon_diversity)
p_value = wilcox$p.value

ggplot(shannon_diversity, aes(x = Diagnosis, y = shannon_diversity),
       fill = Diagnosis, shape = Diagnosis) +
  geom_boxplot(aes(fill = Diagnosis, shape = Diagnosis, colour = Diagnosis)) + 
  geom_point(position = "identity", aes(fill = Diagnosis, shape = Diagnosis, colour = Diagnosis), 
             size = 2, stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21, 21)) +
  scale_color_manual(values=c(Palette_color)) +
  scale_fill_manual(values=c("white", "white", "white", "white")) +
  labs(#title="Progressive IPF vs. Stable IPF", 
    subtitle = glue("Wilcoxon ranked sum p={round(p_value,2)}"),
    x = "",
    y = "Shannon Diversity") +
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18),
    legend.position = "none", panel.background = element_blank(),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=14),
    axis.title.x = element_blank(),
    axis.title = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
    aspect.ratio = 1)

ggsave("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Manuscript/Figures/ShannonDiversity.svg",
       width = 10, height = 10)

#### Correlation matrix ####
# --- Create a correlation matrix with metadata variables ---
# Tutorial on drawing a correlation map using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)

dataframe <- read_csv("metadata.csv")
meta_cor <- dataframe %>%
  dplyr::select(!Microbiology & !`Severity classification` & !Diagnosis & !Sex & !Smoking_status) %>%
  drop_na()
shannon_COVID = shannon_diversity %>%
  dplyr::select(`sample-id`, shannon_diversity)
species_COVID = species_number %>%
  dplyr::select(`sample-id`, species_number)

meta_cor <- left_join(meta_cor, shannon_COVID)
meta_cor <- left_join(meta_cor, species_COVID)

meta_cor <- drop_na(meta_cor)
meta_cor <- column_to_rownames(meta_cor, var="sample-id")

cormat <- round(cor(meta_cor, use="pairwise.complete.obs", method = "spearman"),2)
head(cormat)

library(vegan)
library(reshape2)
melted_cormat <- melt(cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
upper_tri

# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#003f5a", high = "#BC3C29", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "#003f5a", high = "#BC3C29", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 14, hjust = 1),
        axis.text.y = element_text(size=14),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

cor.test(meta_table$species_number, meta_table$`CT: Opacified Lung (%)`, method="spearman")

ggsave("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Manuscript/Figures/corplot.svg",
       width = 12, height = 12)

## with CT cut off of 25%
meta_table = dataframe %>%
  filter(Diagnosis == "COVID-19") %>%
  mutate(CT_group = ifelse(`CT: Opacified Lung (%)` > 25, "High", "Low"))
median(meta_table$`CT: Opacified Lung (%)`)

meta_table <- left_join(meta_table, species_number)

meta_glm = meta_table %>%
  mutate(CT_opacified = `CT: Opacified Lung (%)`/100)

glm_model_prop <- glm(CT_opacified ~ species_number + Age + Sex + Smoking_status,
                      data = meta_glm,
                      family = quasibinomial(link = "logit"))
summary(glm_model_prop)
odds_species = exp(-0.20)
CI = exp(confint(glm_model_prop))

ggplot(meta_table, aes(x = species_number, y=`CT: Opacified Lung (%)`)) + 
  geom_point(aes(color = CT_group)) + geom_smooth(method="lm") + 
  theme_classic() + 
  theme(
    plot.title = element_text(size=20), plot.subtitle = element_text(size=18),
    #legend.position = "none", 
    panel.background = element_blank(),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=45, hjust=1, vjust=1, size=14),
    axis.title = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
    aspect.ratio = 1) +
  labs(y = "Opacified lung on CT (%)", x = "Genus number",
       subtitle = glue("GLM: Odds Ratio = {round(odds_species,2)}, 95% CI = {round(CI[2,1],2)} - {round(CI[2,2],2)}, p-value=0.01")) +
  scale_fill_manual(values=Palette_fill) + scale_color_manual(values=Palette_color)

ggsave("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Manuscript/Figures/glmCT.svg",
       width = 12, height = 12)

abund_COVID = abund_table %>%
  filter(rownames(abund_table) %in% meta_table$`sample-id`)
meta_table = column_to_rownames(meta_table, "sample-id")

Maaslin2(abund_COVID, meta_table, 
         output="Maaslin2-CT",
         fixed_effects = c("CT_group", "Age", "Sex", "Smoking_statusNever"))

#### ISA with CT groups ####
ISA <- multipatt(x = abund_COVID,
                 cluster = meta_table$CT_group,
                 duleg = TRUE)
summary(ISA)
