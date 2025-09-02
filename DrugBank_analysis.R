###DRUG BANK analysis
install.packages("remotes")
remotes::install_github("ropensci/dbparser")
library(dbparser)

#create dvobject
#dvobject contains drugs, drug salts, references, products, and cett sections
dv = parseDrugBank(
  db_path = "full database.xml",
  drug_options = drug_node_options(),
  references_options = references_node_options(),
  cett_options = cett_nodes_options(),  # carriers, enzymes, targets, transporters
  parse_products = TRUE,
  parse_salts = TRUE
)
save(dv, file = "dv.Rdata")

#curate gene symbols

# Create a named vector for mapping: names = old symbols, values = new symbols
symbol_map = setNames(Homo_sapiens$V3, Homo_sapiens$V5)

# Replace symbols where a match is found in the reference
polypeptides$gene_name = ifelse(
  polypeptides$gene_name %in% names(symbol_map),
  symbol_map[polypeptides$gene_name],
  polypeptides$gene_name
)


#extract cett information
targets = dv$cett$targets
synonyms = dv$drugs$synonyms
polypeptides = targets$polypeptides$general_information
save(targets, file = "targets.RData")
write_csv(synonyms, "synonyms.csv")
write_csv(polypeptides, "polypeptides.csv")

#get target_id
targets_info = targets$general_information
targets_info = subset(targets_info, tolower(known_action) == "yes")
write_csv(targets_info, "targets_info.csv")
a = targets$general_information

drug_names = dv$drugs[[1]][, c("drugbank_id", "name")]
colnames(drug_names)[2] = "drug_name"
write_csv(drug_names, "drug_names.csv")

#retain only targets with known action
actions = targets$actions

#merge targets info with known actions by target_id
#now it has target_id, organism, name, known_action, action, position, and drugbank id
targets_known = merge(
  targets_info,
  actions,
  by = "target_id"
)
#remove organism and name columns to 
targets_known = targets_known[, c(1, 4:7)]
write_csv(targets_known, "targets_known.csv")

#merge with polypeptides to get gene names
targets_with_genes = merge(
  targets_known,
  polypeptides,
  by.x = "target_id",
  by.y = "target_id",
  all.x = TRUE
)
write_csv(targets_with_genes, "targets_with_genes.csv")

#merge with drug names
drug_targets_full = merge(
  targets_with_genes,
  drug_names,
  by = "drugbank_id",
  all.x = TRUE
)

#merge with synonyms
drug_targets_full = merge(
  drug_targets_full,
  synonyms,
  by = "drugbank_id",
  all.x = TRUE)
write_csv(drug_targets_full, "drug_targets_full.csv")

#deduplicate synonyms 
drug_targets_unique = drug_targets_full[!duplicated(drug_targets_full[, c("drugbank_id", "gene_name")]), ]
write_csv(drug_targets_unique, "drug_targets_unique.csv")

#create drug-gene interactions for blood samples (no covariates)
interactions_blood = merge(
  blood_no_covars_sign,
  drug_targets_unique,
  by.x = "Gene.Symbol",
  by.y = "gene_name"
)
write_csv(interactions_blood, "interactions_blood.csv")

#drug-gene interactions for blood samples with covariates
interactions_blood_covars = merge(
  blood_with_covars_sign,
  drug_targets_unique,
  by.x = "Gene.Symbol",
  by.y = "gene_name"
)
write_csv(interactions_blood_covars, "interactions_blood_covars.csv")


#build drug-gene network (igraph + ggraph)
install.packages("ggraph")
library(igraph)
library(ggraph)
library(ggplot2)  

# Create edge list
edges = interactions_blood[, c("Gene.Symbol", "drug_name", "action")]  # gene -> drug
colnames(edges) = c("from", "to", "action")

genes = unique(interactions_blood[, c("Gene.Symbol", "meta_LFc")])
genes$regulation = ifelse(genes$meta_LFc > 0, "up", "down")
colnames(genes)[1] = "name"  # for igraph

# Create a unique list of drug nodes
drugs = unique(data.frame(name = interactions_blood$drug_name, regulation = "drug"))
drugs$meta_LFc <- NA
# Combine all nodes
nodes = rbind(genes, drugs)

#build i graph
graph = graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)


# Plot with ggraph
ggraph(graph, layout = "fr") +
  geom_edge_link(alpha = 0.8, color = "gray70") +  # edges
  geom_node_point(aes(color = regulation), size = 14) +  #  nodes
  geom_node_text(aes(label = name), size = 5, color = "black") +  # label centered
  scale_color_manual(
    values = c(
      "up" = "#FF9999",      
      "down" = "#9999FF",
      "drug" = "#A8D5BA"
    )
  ) +
  theme_void()


#network for blood samples with covariates
# Create edge list
edges = interactions_blood_covars[, c("Gene.Symbol", "drug_name", "action")]  # gene -> drug
colnames(edges) = c("from", "to", "action")

genes = unique(interactions_blood_covars[, c("Gene.Symbol", "meta_LFc")])
genes$regulation = ifelse(genes$meta_LFc > 0, "up", "down")
colnames(genes)[1] = "name"  # for igraph

# Create a unique list of drug nodes
drugs = unique(data.frame(name = interactions_blood_covars$drug_name, regulation = "drug"))
drugs$meta_LFc <- NA
# Combine all nodes
nodes = rbind(genes, drugs)

#build i graph
graph = graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)


# Plot with ggraph
ggraph(graph, layout = "fr") +
  geom_edge_link(alpha = 0.8, color = "gray70") +  # edges
  geom_node_point(aes(color = regulation), size = 14) +  #  nodes
  geom_node_text(aes(label = name), size = 3, color = "black") +  # label centered
  scale_color_manual(
    values = c(
      "up" = "#FF9999",      
      "down" = "#9999FF",
      "drug" = "#A8D5BA"
    )
  ) +
  theme_void()


#BRAIN SAMPLES
#no covariates
interactions_brain = merge(
  meta_analysis_no_covars_sign,
  drug_targets_unique,
  by.x = "Gene.Symbol",
  by.y = "gene_name"
)
write_csv(interactions_brain, "interactions_brain.csv")

#with covariates
interactions_brain_covars = merge(
  meta_analysis_with_covars_significant,
  drug_targets_unique,
  by.x = "Gene.Symbol",
  by.y = "gene_name"
)
write_csv(interactions_brain_covars, "interactions_brain_covars.csv")

#network for brain samples with covariates
# Create edge list
edges = interactions_brain_covars[, c("Gene.Symbol", "drug_name", "action")]  # gene -> drug
colnames(edges) = c("from", "to", "action")

genes = unique(interactions_brain_covars[, c("Gene.Symbol", "meta_LFc")])
genes$regulation = ifelse(genes$meta_LFc > 0, "up", "down")
colnames(genes)[1] = "name"  # for igraph

# Create a unique list of drug nodes
drugs = unique(data.frame(name = interactions_brain_covars$drug_name, regulation = "drug"))
drugs$meta_LFc <- NA
# Combine all nodes
nodes = rbind(genes, drugs)

#build i graph
graph = graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)


# Plot with ggraph
ggraph(graph, layout = "fr") +
  geom_edge_link(alpha = 0.8, color = "gray70") +  # edges
  geom_node_point(aes(color = regulation), size = 14) +  #  nodes
  geom_node_text(aes(label = name), size = 5, color = "black") +  # label centered
  scale_color_manual(
    values = c(
      "up" = "#FF9999",      
      "down" = "#9999FF",
      "drug" = "#A8D5BA"
    )
  ) +
  theme_void()

#network for brain samples without covariates
# Create edge list
edges = interactions_brain[, c("Gene.Symbol", "drug_name", "action")]  # gene -> drug
colnames(edges) = c("from", "to", "action")

genes = unique(interactions_brain[, c("Gene.Symbol", "meta_LFc")])
genes$regulation = ifelse(genes$meta_LFc > 0, "up", "down")
colnames(genes)[1] = "name"  # for igraph

# Create a unique list of drug nodes
drugs = unique(data.frame(name = interactions_brain$drug_name, regulation = "drug"))
drugs$meta_LFc <- NA
# Combine all nodes
nodes = rbind(genes, drugs)

#build i graph
graph = graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)


# Plot with ggraph
ggraph(graph, layout = "fr") +
  geom_edge_link(alpha = 0.8, color = "gray70") +  # edges
  geom_node_point(aes(color = regulation), size = 14) +  #  nodes
  geom_node_text(aes(label = name), size = 5, color = "black") +  # label centered
  scale_color_manual(
    values = c(
      "up" = "#FF9999",      
      "down" = "#9999FF",
      "drug" = "#A8D5BA"
    )
  ) +
  theme_void()

#BA9
interactions_ba9 = merge(
  ba9_no_covars_sign,
  drug_targets_unique,
  by.x = "Gene.Symbol",
  by.y = "gene_name"
)
write_csv(interactions_ba9, "interactions_ba9.csv")

interactions_ba9_covars = merge(
  ba9_with_covars_sign,
  drug_targets_unique,
  by.x = "Gene.Symbol",
  by.y = "gene_name"
)
write_csv(interactions_ba9_covars, "interactions_ba9_covars.csv")

#network graph for BA9, no covars
# Create edge list
edges = interactions_ba9[, c("Gene.Symbol", "drug_name", "action")]  # gene -> drug
colnames(edges) = c("from", "to", "action")

genes = unique(interactions_ba9[, c("Gene.Symbol", "meta_LFc")])
genes$regulation = ifelse(genes$meta_LFc > 0, "up", "down")
colnames(genes)[1] = "name"  # for igraph

# Create a unique list of drug nodes
drugs = unique(data.frame(name = interactions_ba9$drug_name, regulation = "drug"))
drugs$meta_LFc <- NA
# Combine all nodes
nodes = rbind(genes, drugs)

#build i graph
graph = graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)


# Plot with ggraph
ggraph(graph, layout = "fr") +
  geom_edge_link(alpha = 0.8, color = "gray70") +  # edges
  geom_node_point(aes(color = regulation), size = 14) +  #  nodes
  geom_node_text(aes(label = name), size = 3, color = "black") +  # label centered
  scale_color_manual(
    values = c(
      "up" = "#FF9999",      
      "down" = "#9999FF",
      "drug" = "#A8D5BA"
    )
  ) +
  theme_void()

# Create edge list
edges = interactions_ba9_covars[, c("Gene.Symbol", "drug_name", "action")]  # gene -> drug
colnames(edges) = c("from", "to", "action")

genes = unique(interactions_ba9_covars[, c("Gene.Symbol", "meta_LFc")])
genes$regulation = ifelse(genes$meta_LFc > 0, "up", "down")
colnames(genes)[1] = "name"  # for igraph

# Create a unique list of drug nodes
drugs = unique(data.frame(name = interactions_ba9_covars$drug_name, regulation = "drug"))
drugs$meta_LFc <- NA
# Combine all nodes
nodes = rbind(genes, drugs)

#build i graph
graph = graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)


# Plot with ggraph
ggraph(graph, layout = "fr") +
  geom_edge_link(alpha = 0.8, color = "gray70") +  # edges
  geom_node_point(aes(color = regulation), size = 14) +  #  nodes
  geom_node_text(aes(label = name), size = 3, color = "black") +  # label centered
  scale_color_manual(
    values = c(
      "up" = "#FF9999",      
      "down" = "#9999FF",
      "drug" = "#A8D5BA"
    )
  ) +
  theme_void()

##BA9+BA46
interactions_pfc = merge(
  pfc_no_covars_sign, 
  drug_targets_unique, 
  by.x = "Gene.Symbol", 
  by.y = "gene_name"
)
write_csv(interactions_pfc, "interactions_pfc.csv")

interactions_pfc_covars = merge(
  pfc_with_covars_sign, 
  drug_targets_unique, 
  by.x = "Gene.Symbol", 
  by.y = "gene_name"
)
write_csv(interactions_pfc_covars, "interactions_pfc_covars.csv")

#network graph for BA9+BA46, no covars
# Create edge list
edges = interactions_pfc[, c("Gene.Symbol", "drug_name", "action")]  # gene -> drug
colnames(edges) = c("from", "to", "action")

genes = unique(interactions_pfc[, c("Gene.Symbol", "meta_LFc")])
genes$regulation = ifelse(genes$meta_LFc > 0, "up", "down")
colnames(genes)[1] = "name"  # for igraph

# Create a unique list of drug nodes
drugs = unique(data.frame(name = interactions_pfc$drug_name, regulation = "drug"))
drugs$meta_LFc <- NA
# Combine all nodes
nodes = rbind(genes, drugs)

#build i graph
graph = graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)


# Plot with ggraph
ggraph(graph, layout = "fr") +
  geom_edge_link(alpha = 0.8, color = "gray70") +  # edges
  geom_node_point(aes(color = regulation), size = 14) +  #  nodes
  geom_node_text(aes(label = name), size = 3, color = "black") +  # label centered
  scale_color_manual(
    values = c(
      "up" = "#FF9999",      
      "down" = "#9999FF",
      "drug" = "#A8D5BA"
    )
  ) +
  theme_void()

#with covariates
# Create edge list
edges = interactions_pfc_covars[, c("Gene.Symbol", "drug_name", "action")]  # gene -> drug
colnames(edges) = c("from", "to", "action")

genes = unique(interactions_pfc_covars[, c("Gene.Symbol", "meta_LFc")])
genes$regulation = ifelse(genes$meta_LFc > 0, "up", "down")
colnames(genes)[1] = "name"  # for igraph

# Create a unique list of drug nodes
drugs = unique(data.frame(name = interactions_pfc_covars$drug_name, regulation = "drug"))
drugs$meta_LFc <- NA
# Combine all nodes
nodes = rbind(genes, drugs)

#build i graph
graph = graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)


# Plot with ggraph
ggraph(graph, layout = "fr") +
  geom_edge_link(alpha = 0.8, color = "gray70") +  # edges
  geom_node_point(aes(color = regulation), size = 14) +  #  nodes
  geom_node_text(aes(label = name), size = 3, color = "black") +  # label centered
  scale_color_manual(
    values = c(
      "up" = "#FF9999",      
      "down" = "#9999FF",
      "drug" = "#A8D5BA"
    )
  ) +
  theme_void()

interactions_blood$label = "blood_no_covars"
interactions_blood_covars$label = "blood_with_covars"
interactions_brain$label = "brain_no_covars"
interactions_brain_covars$label = "brain_with_covars"
interactions_ba9$label = "ba9_no_covars"
interactions_ba9_covars$label = "ba9_with_covars"
interactions_pfc$label = "pfc_no_covars"
interactions_pfc_covars$label = "pfc_with_covars"

drugs_genes_interactions = rbind(interactions_blood, interactions_blood_covars,
                                 interactions_brain, interactions_brain_covars,
                                 interactions_ba9, interactions_ba9_covars,
                                 interactions_pfc, interactions_pfc_covars)

drugs_genes_interactions = drugs_genes_interactions[, c("Gene.Symbol", "name", "drug_name", "label")]

write_csv(drugs_genes_interactions, "drugs_genes_interactions.csv")

#classify the drugs
install.packages(c("webchem", "dplyr", "purrr", "readr"))
install.packages("PubChemR")
library(PubChemR)
library(webchem)
library(dplyr)
library(purrr)

#convert drug names into PubChem CIDs
drug_names_dge = drugs_genes_interactions$drug_name
cids = get_cid(drug_names_dge, from = "name", match = "first", domain = "compound")

#remove rows with NA
cids = cids[!is.na(cids$cid),]
cids = cids[!grepl("^N", cids$query), ]
cids = cids[!grepl("^4", cids$query), ]
cids = cids[!grepl("^6", cids$query), ]
cids = cids[!grepl("^3", cids$query), ]
cids = cids[!grepl("^Dihydro-2-thioxo-5", cids$query),]
write(cids$cid, "cids.txt")

category = dv$drugs$categories

drugs_with_categories = merge(
  drug_names,
  category,
  by.x = "drugbank_id",
  by.y = "drugbank_id"
)

drugs_useful = drugs_with_categories[drugs_with_categories$drug_name %in% cids$query,]

#remove rows with missing mesh_id
drugs_useful$mesh_id[drugs_useful$mesh_id == ""] = NA
drugs_useful = drugs_useful[!is.na(drugs_useful$mesh_id),]

#remove rows with non-pharmacological/medical classification
not_category = c("Hormones", "Proteins", "Peptides", 
                 "Enzyme Inhibitors", "Benzazepines", "Benzodiazepinones", "Nervous System",
                 "Peripheral Nervous System Agents", "Alcohols", "Central Nervous System Agents",
                 "Central Nervous System Depressants", "Gastrointestinal Agents", "Chlorohydrins",
                 "Amines", "Autonomic Agents", "Neurotransmitter Agents",
                 "Peripheral Nervous System Agents", "Phenethylamines", "Cardiovascular Agents",
                 "Hematologic Agents", "Purinergic Agents", "Pyridines", "Sulfur Compounds",
                 "Thiophenes", "Ethers", "Methyl Ethers", "Amines", "Oligopeptides",
                 "Pyrimidines", "Pyrimidinones", "Amides", "Hexoses", "Ketoses",
                 "Imidazoles", "Acids, Acyclic", "Acids, Carbocyclic",
                 "Acetates", "Adenine Nucleotides", "Adrenal Cortex Hormones",
                 "Adrenergic Agents", "Amino Acids", "Amino Acids, Peptides, and Proteins",
                 "Anesthetics, Inhalation", "Anesthetics, Intravenous", "Anesthetics, General",
                 "Anemia, Iron-Deficiency", "Analgesics, Non-Narcotic", "Aminobutyrates",
                 "Alkaloids", "Alkenes", "Aniline Compounds", "Benzene Derivatives", "Benzodiazepines and benzodiazepine derivatives",
                 "Benzopyrans", "Biological Factors", "Carbohydrates", "Coenzymes", "Cytochrome P-450 CYP1A2 Inducers",
                 "Heterocyclic Compounds, Fused-Ring", "Cytochrome P-450 CYP1A2 Inhibitors", 
                 "Cytochrome P-450 CYP2B6 Inhibitors", "Cytochrome P-450 CYP2C19 Inhibitors",
                 "Cytochrome P-450 CYP2C8 Inhibitors", "Cytochrome P-450 CYP2C9 Inhibitors",
                 "Cytochrome P-450 CYP2D6 Inhibitors", "Cytochrome P-450 CYP2E1 Inhibitors",
                 "Cytochrome P-450 Enzyme Inhibitors", "GABA Agents", "GABA Modulators",
                 "Cytochrome P-450 CYP3A Inducers", "Cytochrome P-450 CYP3A4 Inducers",
                 "Cytochrome P-450 Enzyme Inducers", "Monosaccharides", "Adrenergic alpha-2 Receptor Agonists",
                 "Hydroxymethylglutaryl-CoA Reductase Inhibitors", "Alkylating Drugs", "Carbamates",
                 "Eicosanoids", "Lipids", "Fatty Acids", "Vasodilation", "Purinergic P2 Receptor Antagonists",
                 "Purinergic Antagonists", "Nitroso Compounds", "Nitrosourea Compounds", "GABA-A Receptor Agonists",
                 "Muscle Relaxants, Centrally Acting Agents", "Fused-Ring Compounds", "Digoxin and derivatives",
                 "Indoles", "Indole Alkaloids", "Purines", "Ribonucleotides", "Purine Nucleotides", "Nucleic Acids, Nucleotides, and Nucleosides",
                 "Protein Kinase Inhibitors", "Tyrosine Kinase Inhibitors", "Iron Compounds",
                 "Depression, Postpartum", "Organometallic compounds", "Minerals", "Neuromuscular Agents",
                 "Delayed-Action Preparations", "Physiological Phenomena", "Diet, Food, and Nutrition", "Pyrazoles",
                 "Octanols", "Nicotinic Acids", "Pyrrolidines", "Tumor Suppressor Proteins", "Starch", 
                 "Corpus Luteum Hormones", "Cyclodextrins", "Dextrins", "Dietary Carbohydrates", 
                 "Glucans", "Polysaccharides", "Pregnenediones", "Pregnenes", "Progesterone Congeners",
                 "Anticholesteremic Agents", "Noxae", "Photosensitizing Agents", "Toxic Actions" , 
                 "Gastrins", "Gastrointestinal Hormones", "Hormones, Hormone Substitutes, and Hormone Antagonists",         
                 "Serotonergic Drugs Shown to Increase Risk of Serotonin Syndrome", "Purinergic P2Y Receptor Antagonists",                            
                 "Thienopyridines", "Nicotinic Antagonists", "Antineoplastic Agents", "Antineoplastic Agents, Alkylating",                              
                 "Adrenergic Agonists", "Adrenergic alpha-Agonists", "Adjuvants, Anesthesia",                                          
                 "GABA Agonists", "Antihypertensive Agents", "Autacoids", "Fatty Acids, Unsaturated",                                       
                 "Prostaglandins I", "Compounds used in a research, industrial, or household setting", 
                 "Digitalis Glycosides","Protective Agents","Piperazines", "Pyrazines", "Biogenic Amines",                                                
                 "Biogenic Monoamines", "Catechols","Ethanolamines","Cytochrome P-450 CYP2C19 Inducers", 
                 "Cytochrome P-450 CYP2C8 Inducers",                               
                 "Cytochrome P-450 CYP2C9 Inducers", "Sleep Aids, Pharmaceutical", "Cytochrome P-450 CYP2B6 Inducers",                               
                 "Cytochrome P-450 CYP3A Inhibitors", "Cytochrome P-450 CYP3A4 Inhibitors", "Antidotes",
                 "Biological Products", "Central Nervous System Stimulants", "Complex Mixtures", "Cyclohexanes",                                                  
                 "Cycloparaffins", "GABA Antagonists", "GABA-A Receptor Antagonists", "Lactones", "Pharmaceutical Preparations",                                    
                 "Plant Extracts", "Plant Preparations", "Sesquiterpenes", "Terpenes", "Toxins, Biological",                                             
                 "Cytosine Nucleotides", "Nucleotides", "Pyrimidine Nucleotides", "Dermatologicals",                                                
                 "Pregnadienes", "Pregnanes", "Steroids, Fluorinated", "Thiobarbiturates",                                                
                 "Naphthalenes", "Hydrocarbons, Halogenated", "Cytochrome P-450 CYP2E1 Inducers",                               
                 "Phenobarbital and similars", "Phenols",                                                        
                 "Cystine Depleting Agents", "Ethylamines", "Mercaptoethylamines","Sulfhydryl Compounds",                                           
                 "Anti-Infective Agents, Local", "Solvents","Cyclooxygenase Inhibitors", "Hydroxybenzoates",                                               
                 "Salicylates",                                                    
                 "Sensory System Agents", "Ethyl Ethers","Heptanoic Acids", "Cytochrome P-450 CYP3A5 Inhibitors",                             
                 "Hydrocarbons, Fluorinated", "Sulfones","Histamine Agents", "Histamine H1 Antagonists",                                       
                 "Histamine H1 Antagonists, Non-Sedating", "Piperidines", "Chemically-Induced Disorders",                                   
                 "Cytochrome P-450 CYP3A5 Inducers",                               
                 "Cytochrome P-450 CYP3A7 Inducers",                               
                 "Excitatory Amino Acid Agents",                                   
                 "Excitatory Amino Acid Antagonists",                             
                 "Aza Compounds",                                                 
                 "Sleep Initiation and Maintenance Disorders",                     
                 "Piperidones",                                                    
                 "Alkanes",                                                        
                 "Alkanesulfonic Acids",                                           
                 "Hydrocarbons, Acyclic",                                          
                 "Sulfonic Acids",                                                 
                 "Sulfur Acids",                                                   
                 "Enzymes and Coenzymes" ,                                         
                 "Prostaglandins D",                                               
                 "Flavins"    ,                                                    
                 "Pigments, Biological"  ,                                         
                 "Pteridines"         ,                                            
                 "Antimetabolites"        ,                                        
                 "Dicarboxylic Acids"     ,                                        
                 "Glutarates"             ,                                        
                 "Fatty Acids, Volatile"  ,                                        
                 "Valerates"              ,                                        
                 "17-Ketosteroids"         ,                                       
                 "Gonadal Steroid Hormones"  ,                                     
                 "Bridged-Ring Compounds" ,                                        
                 "Flavones" ,                                                      
                 "Flavonoids",                                                     
                 "Pyrans" ,                                                        
                 "Nucleosides",                                                    
                 "Purine Nucleosides" ,                                            
                 "Ribonucleosides"    ,                                            
                 "Aminopyridines"   ,                                              
                 "Hydroxy Acids"    ,                                              
                 "Immunologic Factors" ,                                           
                 "Dioxoles"       ,                                                
                 "Phosphodiesterase Inhibitors" ,                                  
                 "Ammonium Compounds"   ,                                          
                 "Nitrogen Compounds"   ,                                          
                 "Onium Compounds"       ,                                         
                 "Surface-Active Agents"  ,                                        
                 "Abortifacient Agents"   ,                                        
                 "Hormonal Contraceptives for Systemic Use"      ,                 
                 "Luteolytic Agents"   ,                                           
                 "Prostaglandins F, Synthetic"     ,                               
                 "Prostaglandins, Synthetic"   ,                                   
                 "Reproductive Control Agents"   ,                                 
                 "Oxazines",                                         
                 "Cytochrome P-450 CYP3A7 Inhibitors"   ,                          
                 "Carbon Radioisotopes",                                 
                 "Fatty Alcohols" ,                                    
                 "Hexanols",                                        
                 "Food",                                        
                 "Micronutrients",                                                                                            
                 "Abortifacient Agents, Nonsteroidal",                          
                 "Prostaglandins F",                                               
                 "Uterotonic agents", "Organometallic Compounds","Aminobenzoates", "Thiazoles", "Hydrazines" )
drugs_useful = drugs_useful[!drugs_useful$category %in% not_category,]

print(unique(drugs_useful$category))
unique(drugs_useful$drug_name)


drugs_useful_grouped = drugs_useful %>%
  group_by(drug_name) %>%
  summarise(
    drugbank_id = first(drugbank_id),  
    category = list(unique(category)), # collect unique categories
    mesh_id = list(unique(mesh_id))    # collect unique MeSH IDs
  )
drugs_useful_grouped = drugs_useful_grouped[, 1:3]

#remove non-drug drug_names
drugs_genes_interactions = drugs_genes_interactions[!grepl("^N", drugs_genes_interactions$drug_name), ]
drugs_genes_interactions = drugs_genes_interactions[!grepl("^4", drugs_genes_interactions$drug_name), ]
drugs_genes_interactions = drugs_genes_interactions[!grepl("^6", drugs_genes_interactions$drug_name), ]
drugs_genes_interactions = drugs_genes_interactions[!grepl("^5", drugs_genes_interactions$drug_name), ]
drugs_genes_interactions = drugs_genes_interactions[!grepl("^3", drugs_genes_interactions$drug_name), ]
drugs_genes_interactions = drugs_genes_interactions[!grepl("^Dihydro-2-thioxo-5", drugs_genes_interactions$drug_name),]
drugs_genes_interactions = drugs_genes_interactions[!grepl("^Human immunoglobulin G", drugs_genes_interactions$drug_name), ]
drugs_genes_interactions = drugs_genes_interactions[!grepl("^Protein S human", drugs_genes_interactions$drug_name), ]
drugs_genes_interactions = drugs_genes_interactions[!grepl("^Human thrombin", drugs_genes_interactions$drug_name), ]


print(unique(drugs_genes_interactions$drug_name))

devtools::install_github("yduan004/drugbankR")
library(drugbankR)

drugs = unique(drugs_genes_interactions$drug_name)