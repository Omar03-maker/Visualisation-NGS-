# INSTALL PACKAGES AND LIBRARY : 
install.packages("plotly")     
install.packages("readxl")    
install.packages("ggplot2") 
install.packages("plotly")
install.packages("sunburstR")
library("phyloseq")
library("ggplot2")      
library("readxl")       
library("dplyr")        
library("tibble")
library("dplyr") 
library("magrittr")
library("plotly")
library("tidyr")



# UPLOAD FILES : 
OTU_MAT <- read_excel("C:/Users/USER/Documents/DIC_G2B/DIC2/DIC2_G2B_2025/VISULAISATION_NGS/Table_OTU.xlsx", sheet = "OTU_Counts")
TAX_MAT <- read_excel("C:/Users/USER/Documents/DIC_G2B/DIC2/DIC2_G2B_2025/VISULAISATION_NGS/Table_taxonomy.xlsx", sheet = "taxonenew")
SAMPLES_DF <- read_excel("C:/Users/USER/Documents/DIC_G2B/DIC2/DIC2_G2B_2025/VISULAISATION_NGS/Table_sample.xlsx", sheet = "sample")

head(OTU_MAT)
head(TAX_MAT)
head(SAMPLES_DF)

OTU_MAT<- OTU_MAT %>%
  tibble::column_to_rownames("OTU") 
TAX_MAT <- TAX_MAT%>% 
  tibble::column_to_rownames("OTU")
SAMPLES_DF<- SAMPLES_DF%>% 
  tibble::column_to_rownames("echantillon") 

OTU_MAT <- as.matrix(OTU_MAT)
TAX_MAT <- as.matrix(TAX_MAT)
# CREATE THE PHYLOSEQ OBJECT 
Otu = otu_table(OTU_MAT, taxa_are_rows = TRUE)
Tax = tax_table(TAX_MAT)
SAMPLES = sample_data(SAMPLES_DF)
Phylo <- phyloseq(Otu, Tax, SAMPLES)
Phylo

# VISUALIZE DATA : 
sample_names(Phylo)
rank_names(Phylo)
sample_variables(Phylo)

# NORMALIZE NUMBER OF READS :
total = median(sample_sums(Phylo))                   # We use median median sequencing depth in each sample 
standf = function(x, t=total) round(t*(x/sum(x)))
Phylo = transform_sample_counts(Phylo, standf)

# BAR GRAPHS 
# Basic bar graph 
# plot_bar(Phylo, fill = "class")

# Plot Without OTU limits : 
plot_bar(Phylo, fill = "kingdom") + 
  geom_bar(aes(color=kingdom, fill=kingdom), stat="identity", position="stack")
plot_bar(Phylo, fill = "phylum") + 
  geom_bar(aes(color=phylum, fill=phylum), stat="identity", position="stack")
plot_bar(Phylo, fill = "order") + 
  geom_bar(aes(color=order, fill=order), stat="identity", position="stack")
plot_bar(Phylo, fill = "family") + 
  geom_bar(aes(color=family, fill=family), stat="identity", position="stack")

# Regroup together PR vs PO samples 
plot_bar(Phylo, x="sample", fill="phylum") +
  geom_bar(aes(color=sample, fill=order), stat="identity", position="stack") +
  theme_minimal() +
  labs(x="Type d'échantillon", y="Abondance")

# HEAT MAP : 
# Basic heat map 
plot_heatmap(Phylo, method = "NMDS", distance = "bray")
# By clustering we can choose only the heat map the consider 20% of the OTU value 
Phylo_ABUND <- filter_taxa(Phylo, function(x) sum(x > total*0.10) > 0, TRUE)
Phylo_ABUND
otu_table(Phylo_ABUND)[1:8, 1:5]
plot_heatmap(Phylo_ABUND, method = "NMDS", distance = "bray")
# Heat map using different distances and different multivaraite methods
plot_heatmap(Phylo_ABUND, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "species", taxa.order = "species", 
             trans=NULL, low="lightblue", high="blue", na.value="beige")
# Different built-in distances
dist_methods <- unlist(distanceMethodList)
print(dist_methods)

# Heatmap for specific taxonomy group : 
plot_heatmap(Phylo_ABUND, method = "NMDS", distance = "bray", 
             taxa.label = "class", taxa.order = "class", 
             low="beige", high="red", na.value="beige")

# ALPHA DIVERSITY : 
# Chao et shannon diversity 
plot_richness(Phylo, measures="Chao1")
plot_richness(Phylo, measures="Shannon")
# Regroup same fraction 
plot_richness(Phylo, measures=c("Chao1", "Shannon"), 
              x="site", color="sample")
plot_richness(Phylo, measures=c("Chao1", "Shannon"), 
              x="pH", color="sample")
# ORDINATION 
print(colnames(sample_data(Phylo)))
Phylo.ORD <- ordinate(Phylo, "NMDS", "bray") # USE NMDS ET BRAY METHODS 
# GLOBAL PLOT 
plot_ordination(Phylo, Phylo.ORD, type="taxa", color="other", shape= "class", 
                title="NMDS DISTANCE phylum ")
# PLOT BY TAXONOMIC DIVISION 
plot_ordination(Phylo, Phylo.ORD, type="taxa", color="phylum", 
                title="NMDS DISTANCE OTHER", label="class") + facet_wrap(~ phylum, 3)
# Display samples 
plot_ordination(Phylo, Phylo.ORD, type="samples", color="sample", 
                shape="site", title="NMDS DISTANCE SAMPLE & OTHER") +  geom_point(size=5)
# Samples and OTU 
plot_ordination(Phylo, Phylo.ORD, type="split", color="phylum", 
                shape="site", title="biplot", label = "site") +  
  geom_point(size=3)

# NETWORK ANALYSIS 
plot_net(Phylo, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="phylum", point_label="kingdom")
plot_net(Phylo.ORD, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="class", point_label="phylum")


