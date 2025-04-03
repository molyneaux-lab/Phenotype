setwd("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/")

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

#### Qiime2 clean up ####
metadata <- readxl::read_xlsx("Original-files/Clinical_metadata_summarised.xlsx")
metadata$`PHENOTYPE ID` <- paste0("PH", metadata$`PHENOTYPE ID`)
metadata$`PHENOTYPE ID` <- paste0(metadata$`PHENOTYPE ID`, ".BAL")
metadata <- dplyr::rename("sample-id"="PHENOTYPE ID", metadata)

manifest <- read_tsv("Manifest-map.txt")
metadata <- full_join(manifest, metadata)
metadata <- metadata %>%
  dplyr::select(!2:8)
sampleIDs <- metadata$`sample-id`

ddPCR <- read_csv("../NT001_Microbiota-analysis/16S_ddPCR.csv")
ddPCR <- ddPCR %>%
  filter(`sample-id` %in% sampleIDs) %>%
  dplyr::select(`sample-id`, Diagnosis, ddPCR)

setdiff(metadata$`sample-id`, ddPCR$`sample-id`) # 6 samples missing from original sheet
intersect(metadata$`sample-id`, ddPCR$`sample-id`)

metadata <- full_join(metadata, ddPCR, by = "sample-id")
metadata <- metadata %>%
  dplyr::select(!Diagnosis.y) %>%
  drop_na(ddPCR)
# 76 samples needed, 6 samples are dropped (2 no ddPCR data, 4 no sequencing data)
metadata <- dplyr::rename("Diagnosis"="Diagnosis.x", metadata)
metadata$Diagnosis <- gsub("Negatve", "Negative", metadata$Diagnosis)

metadata$Diagnosis <- gsub("Controls" , "Non-Fibrotic Controls", metadata$Diagnosis)
metadata$Diagnosis <- gsub("COVID", "COVID-19", metadata$Diagnosis)

write_csv(metadata, "Original-files/metadata.csv")
metadata <- column_to_rownames(metadata, var = "sample-id")

#### Loading phyloseq object ####
physeq<-qza_to_phyloseq(
  features="Original-files/filtered-metadata-table.qza",
  tree="Original-files/merged-midrooted-tree.qza",
  taxonomy= "Original-files/taxonomy-merged-final.qza")
physeq 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 372 taxa and 76 samples ]
#tax_table()   Taxonomy Table:    [ 372 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 372 tips and 369 internal nodes ]

mapping=metadata
sample_data(physeq)<-mapping
view(sample_data(physeq))
gplots::venn(list(mapping=rownames(mapping), physeq=sample_names(physeq)))
setdiff(rownames(mapping), sample_names(physeq))  

physeq

### Library size ###
ps<-physeq
summarize_phyloseq(ps)
df <- as.data.frame(sample_data(ps)) 
df$LibrarySize <- sample_sums(ps) # similar to rowSums/colSums but automated
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, colour=Diagnosis)) + geom_point() + 
  ylim(0,70000)

### Contaminants by frequency ###
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
# Not a clear segregation between true and false contaminants. This needs a bit more investigation
# There are some contaminants present at high levels in negative controls that aren't identified

plot_frequency(ps, taxa_names(ps)[c(which(contamdf.freq$contaminant))], conc="ddPCR") + 
  xlab("DNA Concentration (ddPCR)")

tax <- as(tax_table(ps), "matrix")
contaminants<-tax[which(contamdf.freq$contaminant),]
write.table(contaminants,file="Original-files/Contaminants-list_freq.txt", col.names=NA, row.names=T,sep="\t")

### Contaminants by prevalence ###
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
write.table(contaminants,file="Original-files/Contaminants-list_prev.txt", col.names=NA, row.names=T,sep="\t")

#### Decide what are contaminants ####
ps.contam <- prune_taxa(contamdf.prev$contaminant, ps)
plot_bar(ps.contam, facet_grid =Diagnosis~.)
# Surprisingly high number of contaminants in COVID
ps.contam <- prune_taxa(contamdf.freq$contaminant, ps)
plot_bar(ps.contam, facet_grid =Diagnosis~.)
# With frequency it's less.

# DECIDE TO USE FREQUENCY OR PREVALENCE

# Identify ASVs above 1000 reads in the list of contaminants. These are big influencers. 
filter <- phyloseq::genefilter_sample(ps.contam, filterfun_sample(function(x) x >= 1000))
ps.contam.1k <- prune_taxa(filter, ps.contam)
otu_table<-as.data.frame(ps.contam@otu_table)
tax_table<-as.data.frame(ps.contam@tax_table) 

library(tibble)
library(dplyr)

### Back to the main UNFILTERED dataset i.e., ps. ###
# Re read the dataset to retain ALL samples, those without ddPCR data won't be removed.
physeq<-qza_to_phyloseq(
  features="Original-files/filtered-metadata-table.qza",
  tree="Original-files/merged-midrooted-tree.qza",
  taxonomy= "Original-files/taxonomy-merged-final.qza")
sample_data(physeq)<-metadata
ps <- physeq

sum(taxa_sums(ps) == 0) # how many taxa aren't present in ANY samples

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
# "1] Min. number of reads = 2367"
# "2] Max. number of reads = 1337757"
# "3] Total number of reads = 3568484"
# "4] Average number of reads = 50978.3428571429"
# "5] Median number of reads = 32291.5"
# "6] Any OTU sum to 1 or less? NO"
# "7] Sparsity = 0.782893539581438"
# "8] Number of singletons = 0"
# "9] Percent of OTUs that are singletons \n 
# "10] Number of sample variables are: 17

#### Set function to normalise samples ####
## Normalising samples i.e., relative abundances
normalizeSample = function(x) {
  x/sum(x)
}

Controls_relative = transformSampleCounts(prunedSet, normalizeSample)
otu_table(Controls_relative)
OTU1 = as(otu_table(Controls_relative), "matrix")
OTUdf = as.data.frame(OTU1)
write.csv(OTUdf, "Original-files/OTUdf.csv")

TAXdf = as(tax_table(Controls_relative), "matrix")
TAXdf = as.data.frame(TAXdf)
write.csv(TAXdf, "Original-files/tax_table.csv")

setwd("Unfiltered")
Controls_Phylum <- aggregate_taxa(Controls_relative, 'Phylum') #7 phyla most likely a contaminant
Phylum_df<-as.data.frame(Controls_Phylum@otu_table)
write.table(Phylum_df,file="Phylum/Phylum-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Family <- aggregate_taxa(Controls_relative, 'Family')
Family_df<-as.data.frame(Controls_Family@otu_table)
write.table(Family_df,file="Family/Family-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Genus <- aggregate_taxa(Controls_relative, 'Genus')
Genus_df<-as.data.frame(Controls_Genus@otu_table)
write.table(Genus_df,file="Genus/Genus-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

### Metadata joining ###
metadata <- rownames_to_column(metadata, var = "sample-id")
Phylum_df <- data.frame(t(Phylum_df))
Phylum_df <- rownames_to_column(Phylum_df, "sample-id")
Phylum_df <- inner_join(metadata, Phylum_df, by = "sample-id")
write_csv(Phylum_df, "Phylum/Phylum-normalised-metadata.csv")

Family_df <- data.frame(t(Family_df))
Family_df <- rownames_to_column(Family_df, "sample-id")
Family_df <- inner_join(metadata, Family_df, by = "sample-id")
write_csv(Family_df, "Family/Family-normalised-metadata.csv")

Genus_df <- data.frame(t(Genus_df))
Genus_df <- rownames_to_column(Genus_df, "sample-id")
Genus_df <- inner_join(metadata, Genus_df, by = "sample-id")
write_csv(Genus_df, "Genus/Genus-normalised-metadata.csv")

#### Create metadata sheet ####
metadata <- read_csv("Original-files/metadata.csv")
metadata <- column_to_rownames(metadata, var = "sample-id")

### Loading phyloseq object ###
physeq<-qza_to_phyloseq(
  features="Original-files/filtered-metadata-table.qza",
  tree="Original-files/merged-midrooted-tree.qza",
  taxonomy= "Original-files/taxonomy-merged-final.qza")

mapping=metadata
sample_data(physeq)<-mapping
view(sample_data(physeq))
gplots::venn(list(mapping=rownames(mapping), physeq=sample_names(physeq)))
setdiff(rownames(mapping), sample_names(physeq))  

physeq
physeq_filtered <- subset_samples(physeq, c(!Diagnosis=="Negative control"))

#### Identification of contaminated samples ####
### Common contaminants ###
setwd("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Unfiltered/Genus")
df <- read_csv("Genus-normalised-metadata.csv") %>%
  dplyr::select(!Unknown)

negative <- df %>%
  filter(Diagnosis == "Negative control")

metadata <- df %>%
  dplyr::select(`sample-id`, Diagnosis, ddPCR)

#### Negative control plot ####
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

df_negative_long <- melt(df_negative, id.vars = c("sample-id",
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
Palette <- hcl.colors(4, palette="Heat2")

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
  scale_fill_manual(values = Palette) +
  labs(title="Relative abundances of likely contaminant taxa present at 5% abundance",
       subtitle = "Negative control specimens: Reagent controls.",
       x = "Mean relative abundance (%)",
       y = "Genus") +
  theme(#Set the title font size
    plot.title = element_text(size=10),
    plot.subtitle = element_text(size = 8),
    legend.position = "right",
    legend.title = element_text(size=8),
    legend.text = element_text(size=8,
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
                               size=8),
    axis.title.x = element_text(size=8),
    axis.title = element_text(size=8),
    axis.text.y = element_text(size=8,
                               face = "italic"),
    axis.line = element_line(size = 0.5,
                             linetype = "solid",
                             colour = "black"),
    aspect.ratio = 1.7
  )

#### Beta-plot ASV ####
x = taxa_sums(physeq_filtered)
# Keep taxa seen at least twice in more than 1% of samples.
filteredset = filter_taxa(physeq_filtered, function(x) sum(x > 2) > (0.1*length(x)), TRUE)
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
Palette <- c("COVID-19" = "#990000",
             "Non-Fibrotic Controls" = "#588300")
myplotdiagnosis <- p
myplotdiagnosis + theme_bw() + stat_ellipse() + theme_classic() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  ggrepel::geom_text_repel(label=rownames(p$data),colour="black", size=3, max.overlaps = 25)+
  scale_colour_manual(values=c(Palette))

#### Top ten genera ####
setwd("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Unfiltered/Genus")
df <- read_csv("Genus-normalised-metadata.csv")
## df_filtered includes samples that are most likely contaminated by environmental controls. 
df_filtered <- df %>% filter(Diagnosis=="COVID-19" | Diagnosis=="Non-Fibrotic Controls")

abund_table <- df_filtered %>%
  select(1,Actinomyces:ncol(df_filtered))
abund_table <- column_to_rownames(abund_table, var="sample-id")
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100 # change to percentages

meta_table <- df_filtered %>%
  select(`sample-id`, Diagnosis)

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

Palette <- hcl.colors(11, palette="Spectral")

ggplot(sample_data_long,
       aes(x = `sample-id`,
           y = value,
           fill = Genus)) +
  geom_bar(stat = "identity",
           width = 0.7) +
  scale_fill_manual(values = Palette) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100.1)) +
  facet_grid(.~Diagnosis,
             scales = "free_x",
             drop=TRUE)+
  labs(title="Ten most abundant genera in COVID patients \nand non-fibrotic disease controls",
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

#### Alpha-plot ASV, contaminants removed ####
setwd("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/")
metadata <- read_csv("Original-files/metadata.csv")

metadata <- metadata %>%
  drop_na(ddPCR) %>%
  filter(!`sample-id` == "PH1039.BAL" &
           !`sample-id` =="BRU.03853"&
           !`sample-id` =="PH1046.BAL" &
           !`sample-id`=="PH1024.BAL" &
           !`sample-id`=="BRU.03853" &
           !`sample-id` == "BRU.03815" &
           !`sample-id`=="BRU.1032")

metadata <- column_to_rownames(metadata, var = "sample-id")
metadata$Diagnosis <- gsub("Negatve", "Negative", metadata$Diagnosis)

physeq<-qza_to_phyloseq(
  features="Original-files/filtered-metadata-table.qza",
  tree="Original-files/merged-midrooted-tree.qza",
  taxonomy= "Original-files/taxonomy-merged-final.qza")

mapping=metadata
sample_data(physeq)<-mapping
view(sample_data(physeq))
gplots::venn(list(mapping=rownames(mapping), physeq=sample_names(physeq)))
setdiff(rownames(mapping), sample_names(physeq))  

physeq
physeq_filtered <- subset_samples(physeq, c(!Diagnosis=="Negative control"))
sample_data(physeq_filtered)$Diagnosis <- factor(sample_data(physeq_filtered)$Diagnosis,
                                                 levels =c("Non-Fibrotic Controls", "COVID-19"))
p = plot_richness(physeq_filtered, x="Diagnosis", color="Diagnosis", measures=c("Observed","Shannon", "Chao1"))

Palette <- c("COVID-19" = "#990000",
             "Non-Fibrotic Controls" = "#588300")

p + geom_boxplot(data = p$data, aes(x = Diagnosis, y = value, color = NULL), alpha = 0.1) +
  scale_colour_manual(values=c(Palette)) + 
  theme_classic()+ theme(axis.text.x=element_blank())

alpha_diversity <- p$data
alpha_diversity_shannon <- alpha_diversity %>% filter(variable=="Shannon")
hist(alpha_diversity_shannon$value)
shapiro.test(alpha_diversity_shannon$value) # p=0.56

qqnorm(alpha_diversity_shannon$value) #normally distributed
qqline(alpha_diversity_shannon$value)

wilcox.test(value ~ Diagnosis, alpha_diversity_shannon) # p=0.2, no evidence to suggest the means are significantly different.

alpha_diversity_chao <- alpha_diversity %>% filter(variable=="Chao1")
hist(alpha_diversity_chao$value)
shapiro.test(alpha_diversity_chao$value) # p<0.05, evidence to reject that the data did not come from normally distribution
qqnorm(alpha_diversity_chao$value)
qqline(alpha_diversity_chao$value)
wilcox.test(value ~ Diagnosis, alpha_diversity_chao) # p=0.35

## No significant differences between the alpha diversity indices between these groups. 

#### Beta-plot ASV, contaminants removed ####
x = taxa_sums(physeq_filtered)
# Keep taxa seen at least twice in more than 1% of samples.
filteredset = filter_taxa(physeq_filtered, function(x) sum(x > 2) > (0.1*length(x)), TRUE)
filteredset

minTotRelAbun = 0.00005
x = taxa_sums(filteredset)
keepTaxa = which((x / sum(x)) > minTotRelAbun)
prunedSet_beta = prune_taxa(names(keepTaxa), filteredset)
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
Palette <- c("COVID-19" = "#990000",
             "Non-Fibrotic Controls" = "#588300")
myplotdiagnosis <- p
myplotdiagnosis + theme_bw() + stat_ellipse() + theme_classic() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  scale_colour_manual(values=c(Palette))

#### Correlation matrix of metadata variables ####
# Tutorial on drawing a correlation map using ggplot2
# by Umer Zeeshan Ijaz (http://userweb.eng.gla.ac.uk/umer.ijaz)

dataframe <- read_csv("Original-files/metadata.csv")
meta_table <- dataframe %>%
  select(!Microbiology & !`Severity classification` & !Diagnosis) %>%
  drop_na()

meta_table <- column_to_rownames(meta_table, var="sample-id")

cormat <- round(cor(meta_table, use="pairwise.complete.obs"),2)
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
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
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
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
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

#### Top ten genera, contaminants removed ####
setwd("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Unfiltered/Genus")
df <- read_csv("Genus-normalised-metadata.csv")
## df_filtered includes samples that are most likely contaminated by environmental controls. 
df_filtered <- df %>% filter(Diagnosis=="COVID-19" | Diagnosis=="Non-Fibrotic Controls") %>%
  # Remove these samples based on PcoA plot (ASV) and the initial bar plot
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

Palette <- hcl.colors(11, palette="Spectral")

ggplot(sample_data_long,
       aes(x = `sample-id`,
           y = value,
           fill = Genus)) +
  geom_bar(stat = "identity",
           width = 0.7) +
  scale_fill_manual(values = Palette) +
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

#### PERMANOVA ####
adonis2(abund_table ~ Diagnosis,
        data = meta_table, permutations = 999, method = "bray")

meta_table <- df_filtered %>%
  dplyr::select(`sample-id`: ddPCR) %>% 
  filter(Diagnosis == "COVID-19") %>% drop_na()
abund_table <- df_filtered %>%
  filter(`sample-id` %in% meta_table$`sample-id`) %>%
  dplyr::select(Actinomyces:ncol(df_filtered))

adonis2(abund_table ~ Age + `BAL Neutrophil (%)` +
          Microbiology + `CT: Opacified Lung (%)` + `Severity classification`,
        data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ `Blood Neutrophils` + `Blood lymphocytes` +
          `Blood CRP` + `Blood Ferritin` + `Blood Fibrinogen`,
        data = meta_table, permutations = 999, method = "bray")
adonis2(abund_table ~ `FEV1` + `%FEV1` + `FVC` + `%FVC` +
          TLCO + `%TLCO` + ddPCR,
        data=meta_table, permutations = 999, method = "bray")

#### ddPCR analysis ####
## Diagnosis ##
setwd("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/")
df <- read_csv("Original-files/metadata.csv")
meta_table <- df %>%
  select(`sample-id`:ddPCR)

meta_ddPCR <- meta_table %>% filter(!Diagnosis == "Negative control") %>% 
  select(`sample-id`, Diagnosis, Microbiology, ddPCR, `Severity classification`)

bacterial_burden_stats <- meta_ddPCR %>%
  group_by(Diagnosis) %>%
  mutate(median_burden = median(ddPCR),
         mean_burden = mean(ddPCR),
         q1 = quantile(ddPCR, 0.25),  # 1st quartile
         q3 = quantile(ddPCR, 0.75), # 3rd quartile
         norm_test = shapiro.test(ddPCR)$p.value) %>% # Get p-value from shapiro.test
  filter(!Diagnosis=="Negative control")

wilcox.test(ddPCR ~ Diagnosis, data = meta_ddPCR) 

Palette <- c("COVID-19" = "#990000",
             "Non-Fibrotic Controls" = "#588300")
p <- ggplot(bacterial_burden_stats, aes(x=Diagnosis, y=ddPCR, fill=Diagnosis)) +
  geom_dotplot(binaxis="y", stackdir = "center", binwidth = 0.1, position="dodge") +
  geom_errorbar(aes(x=Diagnosis, ymin=q1, ymax=q3), width=0.3, color='black', linewidth=1) +
  scale_y_log10() + theme_classic() + 
  labs(title="Bacterial burden across diagnosis groups",
       y="Bacteridal burden (16S rRNA gene/mL of BAL)",
       x="") +
  scale_fill_manual(values=c(Palette)) +
  theme(axis.text.x = element_text(angle=45, hjust=1))
p

## Severity classification ##
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

#### Ordered box plot ####
setwd("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/")

df <- read_csv("Unfiltered/Genus/Genus-normalised-metadata.csv") 
df <- df %>%
  filter(!Diagnosis == "Negative control") %>%
  filter(!`sample-id` == "PH1039.BAL" &
           !`sample-id` =="BRU.03853"&
           !`sample-id` =="PH1046.BAL" &
           !`sample-id`=="PH1024.BAL" &
           !`sample-id`=="BRU.03853" &
           !`sample-id` == "BRU.03815" &
           !`sample-id`=="BRU.1032")

meta_table <- df %>%
  select(`sample-id`:ddPCR)

abund_table <- df %>%
  select(`sample-id`, Actinomyces:ncol(df))
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

df_long <- melt(df_genus_top, id.vars = "Diagnosis", 
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

Palette <- hcl.colors(10, palette="Spectral")

p <- p + scale_fill_manual(values = Palette)+
  geom_text(data=dat_text, mapping= aes(x=x, y=y, label=label), size = 6) +
  geom_text(data=dat_text_2, mapping= aes(x=x, y=y, label=label), size = 6) +
  geom_text(data=dat_text_3, mapping= aes(x=x, y=y, label=label), size = 6) +
  geom_text(data=dat_text_4, mapping= aes(x=x, y=y, label=label), size = 6) +
  geom_text(data=dat_text_5, mapping= aes(x=x, y=y, label=label), size = 6) +
  theme(axis.text.x = element_text(angle=90,
                                   size = 10,
                                   hjust=1),
        legend.text = element_text(face = "italic"),
        axis.title.y = element_text(),
        strip.text = element_text(size=18),
        plot.subtitle = element_text(size=14,
                                     face = "italic"))
p  

#### Stacked bar plot, COVID ####
setwd("~/OneDrive - Imperial College London/Projects/Experiments/NT008_Phenotype/Unfiltered/Phylum/")
df <- read_csv("Phylum-normalised-metadata.csv")
df_filtered <- df %>% filter(Diagnosis=="COVID-19" | Diagnosis=="Non-Fibrotic Controls") %>%
  # Remove these samples based on PcoA plot (ASV) and the initial bar plot
  filter(!`sample-id` == "PH1039.BAL" &
           !`sample-id` =="BRU.03853"&
           !`sample-id` =="PH1046.BAL" &
           !`sample-id`=="PH1024.BAL" &
           !`sample-id`=="BRU.03853" &
           !`sample-id` == "BRU.03815" &
           !`sample-id`=="BRU.1032")

abund_table <- df_filtered %>%
  select(1,Actinobacteria:ncol(df_filtered))
abund_table <- column_to_rownames(abund_table, var="sample-id")
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100 # change to percentages

meta_table <- df_filtered %>%
  select(`sample-id`, Diagnosis)

# Make a dataframe "top" that contains the top 10 most abundant genera. You can alter these as you want.
top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- 5
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

df_phylum_top <- cbind(meta_table, top)
sample_data_long <- melt(df_phylum_top, id.vars = c("sample-id", "Diagnosis"),
                         variable.name = "Phylum")
sample_data_long
taxa_list

kruskal.test(Firmicutes ~ Diagnosis, df_phylum_top) #p=0.7
kruskal.test(Bacteroidetes ~ Diagnosis, df_phylum_top) #p=0.65
kruskal.test(Proteobacteria ~ Diagnosis, df_phylum_top) #p=0.006, contaminant
kruskal.test(Actinobacteria ~ Diagnosis, df_phylum_top) #p=0.32
kruskal.test(Fusobacteria ~ Diagnosis, df_phylum_top) #p=0.09


Palette <- hcl.colors(6, palette= "viridis")

df_long_summarised <- sample_data_long %>%
  group_by(Diagnosis, Phylum) %>%
  summarise(mean_value=mean(value),
            median=median(value),
            sd=sd(value),
            q1 = quantile(value, 0.25),  # 1st quartile
            q3 = quantile(value, 0.75)) # 3rd quartile

ggplot(df_long_summarised, aes(x=Diagnosis, y=mean_value)) + 
  geom_bar(aes(y = mean_value, x = Diagnosis, fill = Phylum),
           stat="identity", position=position_stack()) + 
  scale_fill_manual(values = Pallete) +
  labs(title = "Phyla abundances of PHENOTYPE cohort",
       y = "Relative abundance (%)") + 
  theme(axis.title.y = element_text(size=12)) +
    bbplot::bbc_style() 

#### ISA ####
library(indicspecies)

df <- read_csv("../Genus/Genus-normalised-metadata.csv")
df$Prevotella <- df$Prevotella + df$X.Prevotella. 
df <- df %>%
  dplyr::select(-X.Prevotella., -Unknown)
df_filtered <- df %>% filter(Diagnosis=="COVID-19" | Diagnosis=="Non-Fibrotic Controls") %>%
  # Remove these samples based on PcoA plot (ASV) and the initial bar plot
  filter(!`sample-id` == "PH1039.BAL" &
           !`sample-id` =="BRU.03853"&
           !`sample-id` =="PH1046.BAL" &
           !`sample-id`=="PH1024.BAL" &
           !`sample-id`=="BRU.03853" &
           !`sample-id` == "BRU.03815" &
           !`sample-id`=="BRU.1032")
meta_table <- df_filtered %>%
  dplyr::select(`sample-id`:Diagnosis)
abund_table <- df_filtered %>%
  dplyr::select(Actinomyces:ncol(df_filtered))

ISA <- multipatt(x = abund_table,
                 cluster = meta_table$Diagnosis,
                 duleg = TRUE)
summary(ISA)

#Group COVID-19  #sps.  8 
#                 stat p.value   
#  Halomonas      0.921   0.005 ** #contaminant? 
#  Fusobacterium  0.835   0.030 * 
#  Actinomyces    0.815   0.020 * 
#  Bradyrhizobium 0.792   0.005 ** #contaminant?
#  Catonella      0.792   0.035 * 
#  Pseudomonas    0.739   0.005 ** #contaminant?
#  Oribacterium   0.724   0.050 * 
#  Mycobacterium  0.674   0.005 **
  
#  Group Non-Fibrotic Controls  #sps.  4 
#                   stat p.value   
#  Rothia          0.853   0.020 * 
#  Actinobacillus  0.778   0.005 **
#  Salinibacterium 0.674   0.005 ** #contaminant?
#  Methylibium     0.564   0.005 ** #contaminant?
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

#### Heatmap ####
df <- read_csv("Unfiltered/Genus/Genus-normalised-metadata.csv")
df_filtered <- df %>% filter(Diagnosis=="COVID-19") %>%
  # Remove these samples based on PcoA plot (ASV) and the initial bar plot
  filter(!`sample-id` == "PH1039.BAL" &
           !`sample-id` =="BRU.03853"&
           !`sample-id` =="PH1046.BAL" &
           !`sample-id`=="PH1024.BAL" &
           !`sample-id`=="BRU.03853" &
           !`sample-id` == "BRU.03815" &
           !`sample-id`=="BRU.1032")

abund_table <- df_filtered %>%
  dplyr::select(1,Actinomyces:ncol(df_filtered))
meta_table <- df_filtered %>%
  dplyr::select(1:19)

abund_table <- column_to_rownames(abund_table, var = "sample-id")
abund_table <- abund_table/rowSums(abund_table)*100 # change to percentages

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

# Heatmap will be made using a distance matrix
data_dist_rows <- vegdist(top_other, method = "bray")
row_clustering <- hclust(data_dist_rows, "average")
data_dist_columns <- vegdist(t(top_other), method = "bray")
col_clustering <- hclust(data_dist_columns, "average")

colour_palette <- colorRampPalette(colors=c("white",
                                            "orange",
                                            "darkorange",
                                            "red",
                                            "darkred",
                                            "black"))(100)

# Add annotation to heatmap AFTER distance matrix has been calculated
# the metadata does not influence the distance matrix. 
# only interested in: FVC, DLCO, CT, Fibrosis
annot_df1 <- data.frame(Severity = meta_table$`Severity classification`)
annot_df2 <- data.frame(CT = meta_table$`CT: Opacified Lung (%)`)
annot_df3 <- data.frame(Neutrophil = meta_table$`BAL Neutrophil (%)`)
annot_df4 <- data.frame(ddPCR = meta_table$ddPCR)

col1 = list(Severity = c("Mild"="#588300",
                          "Moderate" = "orange",
                          "Severe" = "firebrick"))
library(circlize)
col2 = list(CT = colorRamp2(c(0,25,50,75,100), 
                            c("lightblue", "white", "orange", "firebrick", "darkred")))
col3 = list(Neutrophil = colorRamp2(c(0,5,10,15), 
                            c("white", "darkolivegreen2", "darkolivegreen3", "darkolivegreen")))
col4 = list(ddPCR = colorRamp2(c(5000,100000,
                                 200000,300000,400000,500000),
                               c("white", "lightblue1", "skyblue2",
                                 "deepskyblue","dodgerblue4", "navyblue")))
library(ComplexHeatmap)
sidebar_annotation1 <- rowAnnotation(df = annot_df1, # Dataframe containing treatment groups
                                     col = col1, # The list of treatment groups and their assigned colours
                                     show_annotation_name = TRUE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Severity Classification", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size
sidebar_annotation2 <- rowAnnotation(df = annot_df2,
                                     col=col2,
                                     show_annotation_name = T,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "CT Opacification (%)", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size## CT annotation
sidebar_annotation3 <- rowAnnotation(df = annot_df3,
                                     col=col3,
                                     show_annotation_name = T,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "BAL Neutrophils (%)", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size## CT annotation
sidebar_annotation4 <- rowAnnotation(df = annot_df4,
                                     col=col4,
                                     show_annotation_name = T,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "ddPCR", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size## CT annotation

heatmap <- Heatmap(as.matrix(top_other), # The dataframe containing the heatmap data
                   name = "Proportion", # Name is used as the title of the heatmap legend if shown
                   col = colour_palette, # The predefined colour palette
                   cluster_rows = row_clustering, # Cluster the rows using the predefined clustering
                   #cluster_columns = col_clustering, # Cluster the columns using the predefined clustering
                   show_row_names = TRUE, # Show or hide the row names, TRUE to show rownames
                   row_names_gp = gpar(fontsize = 6), # Row name font size
                   column_names_gp = gpar(fontsize = 8, # Column name font size
                                          fontface = "italic"), # Column names in italics
                   column_title_gp = gpar(fontsize = 8), # Column title font size
                   row_title = "Samples", # Set row title
                   row_labels = rownames(top_other),
                   row_names_side = c("left"),
                   row_title_gp = gpar(fontsize = 8), # Set row title font size
                   column_title = "", # Set column title
                   column_title_side = "bottom", # Set column title font size
                   heatmap_legend_param = list(title = "Relative\nabundance\n(%)", # Set legend title
                                               at = c(0,20,40,60,80,100), # Set legend scale breaks
                                               labels = c("0","20","40","60","80","100"), # Set legend scale labels
                                               title_gp = gpar(fontsize = 8), # Set legend title font size
                                               labels_gp = gpar(fontsize = 8))) # Set legend label font size

p <- heatmap + sidebar_annotation1 + sidebar_annotation2 + 
  sidebar_annotation3 + sidebar_annotation4
p
