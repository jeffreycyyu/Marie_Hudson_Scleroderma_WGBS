#=====
#THIS IS THE SCRIPT FOR CREATING .TXT FORMAT TABLES OF ALL DMR ANNOTATIONS

#chromosomes involved in this scrpt (for now, place this after the )
chromosomes = c(1:22,"X")

#load all files from PT1
for (chr in chromosomes){
  load(paste0("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Annotation/PT1_Results/PT1:Annotation:Chr", chr, ".RData"))
}


#list of choices for chomosomes stratitfication to do IPA of seperate chromosomes
stratify_by_chr_choices <- list(c(1:22,"X"), c(1:22), "X")
#list of the same length as stratify_by_chr_choices to name resulting output txt files
stratify_labels <- list("Genome", "Autosomal", "ChrX")


#run for all strata choices
for (stratify_index in 1:length(stratify_by_chr_choices)){
  
  #respecify chormosomes
  chromosomes <- stratify_by_chr_choices[[stratify_index]]
  #specify label for this stratification
  stratification_label <- stratify_labels[[stratify_index]]
  
  
  
  #PRIMARY_ANALYSIS_RESULTS======================================================================================================
  #bind all chromosomes together
  genes_primary_analysis_genome <- do.call(rbind, mget(paste0("genes_alpha_1_chr", chromosomes)))
  #order the table of all gene regions by DMR p-values (now the table is for all 23 chromosomes and ordered)
  genes_primary_analysis_genome <- genes_primary_analysis_genome[order(as.numeric(genes_primary_analysis_genome[,"DMR_p.value"])),]
  
  #this is the list of all annotated genes regardless or region
  genes_primary_analysis_genome_list <- unique(as.character(genes_primary_analysis_genome[,"Searched_Gene_Symbol"]))
  #number of genes
  unique_primary = length(genes_primary_analysis_genome_list)
  
  #empty array 
  my_array_primary <- array(NA, c(unique_primary, 4), dimnames = list(NULL, c("Gene", "p-value", "# Assoc. DMRs", "Chr")))
  
  #DO NOT DELETE: genes are in (HGNC Symbol) format
  #for each gene in the gene list...
  for (i in 1:unique_primary){
    
    #find the vector containing row indices of the gene in the annotation dataset
    temp = grep(as.character(genes_primary_analysis_genome_list[i]), genes_primary_analysis_genome[,"Searched_Gene_Symbol"])
    
    #find the DMRs associated with this gene
    unique_dmr = unique(as.character(genes_primary_analysis_genome[temp,"DMR_Start"]))
    #how many DMRs are associated with this gene
    unique_dmr_length = length(unique_dmr)
    
    #information we want to keep
    repl = cbind(as.character(genes_primary_analysis_genome_list[i]), #the gene symbol
                 min(as.numeric(as.character(genes_primary_analysis_genome[temp, "DMR_p.value"]))), #the minimum DMR p-value (NOT the product of all DMR p-values because our regions were defined arbitrarily; by spacing of >=200bp)
                 unique_dmr_length, #number of DMRs associated with this gene
                 genes_primary_analysis_genome[temp[1], "Chr"] #chromosmes
    )
    
    #append to array
    my_array_primary[i, 1:4] <- repl
    
  }
  
  my_array_primary <- na.omit(my_array_primary)
  
  #create table in .txt format for IPA input
  write.table(my_array_primary, file = paste0("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Annotation/PT2_Results/Annotation_Primary_", stratification_label, ".txt"), 
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  
  
  
  
  #SECONDARY_ANALYSIS_RESULTS======================================================================================================
  #bind all chromosomes together
  genes_secondary_analysis_genome <- do.call(rbind, mget(paste0("genes_alpha_2_chr", chromosomes)))
  #order the table of all gene regions by DMR p-values (now the table is for all 23 chromosomes and ordered)
  genes_secondary_analysis_genome <- genes_secondary_analysis_genome[order(as.numeric(genes_secondary_analysis_genome[,"DMR_p.value"])),]
  
  #this is the list of all annotated genes regardless or region
  genes_secondary_analysis_genome_list <- unique(as.character(genes_secondary_analysis_genome[,"Searched_Gene_Symbol"]))
  #number of genes
  unique_secondary = length(genes_secondary_analysis_genome_list)
  
  #empty array 
  my_array_secondary <- array(NA, c(unique_secondary, 4), dimnames = list(NULL, c("Gene", "p-value", "# Assoc. DMRs", "Chr")))
  
  #DO NOT DELETE: genes are in (HGNC Symbol) format
  #for each gene in the gene list...
  for (i in 1:unique_secondary){
    
    #find the vector containing row indices of the gene in the annotation dataset
    temp = grep(as.character(genes_secondary_analysis_genome_list[i]), genes_secondary_analysis_genome[,"Searched_Gene_Symbol"])
    
    #find the DMRs associated with this gene
    unique_dmr = unique(as.character(genes_secondary_analysis_genome[temp,"DMR_Start"]))
    #how many DMRs are associated with this gene
    unique_dmr_length = length(unique_dmr)
    
    #information we want to keep
    repl = cbind(as.character(genes_secondary_analysis_genome_list[i]), #the gene symbol
                 min(as.numeric(as.character(genes_secondary_analysis_genome[temp, "DMR_p.value"]))), #the minimum DMR p-value (NOT the product of all DMR p-values because our regions were defined arbitrarily; by spacing of >=200bp)
                 unique_dmr_length, #number of DMRs associated with this gene
                 genes_secondary_analysis_genome[temp[1], "Chr"] #chromosmes
    )
    
    #append to array
    my_array_secondary[i, 1:4] <- repl
    
  }
  
  my_array_secondary <- na.omit(my_array_secondary)
  
  #create table in .txt format for IPA input
  write.table(my_array_secondary, file = paste0("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Annotation/PT2_Results/Annotation_Secondary_", stratification_label, ".txt"), 
              append = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  
  
  
  #save
  save.image(file = paste0("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Annotation/PT2_Results/PT2:Annotations_", stratification_label, ".RData"))
  
}
