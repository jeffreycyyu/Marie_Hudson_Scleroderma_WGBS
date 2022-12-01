#=====
#THIS IS THE SCRIPT FOR ANNOTING OUR TABLE OF ALL CPG REGIONS WITHIN THE GENOME WITH SIGNIFICANCT DISEASE P-VALUES AS DEFINED BY < 0.05 (USED IN LU ET AL.'S BUMPHUNTER PAPER)


library(annotatr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(stringr)


#load table of CpG regions and associated p-values
load("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Partitioning_Regions_Spacing/PT3_Results/PT3:Constructing_Spacing_DataFrame:GENOME.RData")

#choose annotation shortcuts as specified by the annotatr package; for a full list of possible inputs (not including any custom annotations), type "builtin_annotations()" into console
annotation_shortcuts <- list(head(builtin_annotations()[grepl("hg19", builtin_annotations())], n=11)[c(1:7,10:11)], head(builtin_annotations()[grepl("hg19", builtin_annotations())], n=11)[c(1:7,10:11)])

#choice of reference genome used
my_reference_genome <- 'hg19'

#chromosomes involved in this scrpt
chromosomes = c(1:22,"X")
#=====


#bonferonni p value of 0.05/number of CpG regions
bonferroni_alpha <- 0.05/nrow(Reg_RES_ALL_cg)
#p value of 0.05 used in Tianyuan's study that we want to compare results to
standard_alpha = 0.05
#set the p-value threshold, in this case we use the 0.05 used in Tianyuan's study, and the bonferroni corrected p value for CpG regions
my_alphas <- c(bonferroni_alpha, standard_alpha)

#number of DMRs in the bonferroni_alpha primary analysis
nrow(Reg_RES_ALL_cg[which(Reg_RES_ALL_cg[,"Disease.p.value"] < bonferroni_alpha), ])

#number of DMRs in the standard_alpha secondary analysis
nrow(Reg_RES_ALL_cg[which(Reg_RES_ALL_cg[,"Disease.p.value"] < standard_alpha), ])




#split annotation by chromosome to avoid annotating to a gene that belongs to a different chormosome
for (chr in chromosomes){
  
  #loop through all "region" definitions (note: for comparison to bumphunter, there is only one annotation and it does not extend to any upstream regions)
  for (alpha_use in my_alphas){
    
    
    #this index helps locate annotated regions of interest and label objects in teh enviornment
    alpha_index <- grep(alpha_use, my_alphas)
    
    
    #choose annotation criteria
    annots = annotation_shortcuts[[alpha_index]]
    #build annotations
    annotations = build_annotations(genome = my_reference_genome, annotations = annots)
    region_sites <- as.data.frame(annotations)
    #filter by chromosome
    region_sites <- region_sites[region_sites[,"seqnames"] == paste0('chr',chr),]
    #omit any row that contains NA values (note: each row is representing one annotated region of a gene, not the entire gene; e.g. omiting the row reporting FLT4's promoter will *NOT* omit other regions of FLT4 such as introns, CpG shores, etc.)
    region_sites <- na.omit(region_sites)
    
    
    #order all CpG regions by p-value
    res_order_now <-Reg_RES_ALL_cg[order(Reg_RES_ALL_cg$Disease.p.value),]
    #filter by p-value and exclude all that do not pass the threshold
    res_order_now <-res_order_now[res_order_now$Disease.p.value <= my_alphas[alpha_index],]
    #filter by chromosome
    res_order_per_chr <- res_order_now[res_order_now[, "chr"] == chr,]
    #the number of check to make where each annotated region of the chroomosomes is checked agaisnt each DMR of the chromosomes
    n_check_sites = nrow(region_sites)*nrow(res_order_per_chr)
    
    
    #generate list of all sequences where each sequence represents the positions covered by a DMR, each element of the list represents a unique DMR
    sequences_DMR_regions <- Map(function(num1, num2) seq(as.numeric(num1), as.numeric(num2)), res_order_per_chr[,which(colnames(res_order_per_chr) == "start")], res_order_per_chr[,which(colnames(res_order_per_chr) == "end")])
    #generate list of all sequences where each sequence represents the positions covered by a gene region, each element of the list represents a unique gene region (NOT nessesarily a unique gene region type, since genes can have multiple introns, CpG islands, etc.)
    sequences_gene_regions <- Map(function(num1, num2) seq(as.numeric(num1), as.numeric(num2)), region_sites[,which(colnames(region_sites) == "start")], region_sites[,which(colnames(region_sites) == "end")])
    
    #find the number of intersected positions between pairwise intersection between vector of DMR positions and vector of gene region positions
    #format of this is as a list (# of elements = # of DMRs) of lists (# of elements = # of annotated gene regions) containing positions of intersections between DMR and gene region in question
    intersected_positions <- do.call("cbind", lapply(sequences_gene_regions, function(x) lapply(sequences_DMR_regions, function(y) length(intersect(x,y)))))
    #these are all the matrix indices for values where
    intersection_matrix_indices <- which(intersected_positions != 0, arr.ind = TRUE)
    #label columns
    colnames(intersection_matrix_indices) <- c("DMR index (use these as row indices for res_order_per_chr)", "gene region index (use these as row indices for region_sites)")
    
    
    #create a table to store all genes and thier associated significant DMR's (a significant DMR is required for a gene to be considered)
    gene_check <- cbind(rep(chr, nrow(intersection_matrix_indices)), #chromosome
                        region_sites[unlist(intersection_matrix_indices[,2]), which(colnames(region_sites) == "start")], #gene region start site
                        region_sites[unlist(intersection_matrix_indices[,2]), which(colnames(region_sites) == "end")], #gene region end site
                        region_sites[unlist(intersection_matrix_indices[,2]), which(colnames(region_sites) == "symbol")], #gene symbol
                        region_sites[unlist(intersection_matrix_indices[,2]), which(colnames(region_sites) == "type")], #gene region type
                        res_order_per_chr[unlist(intersection_matrix_indices[,1]), which(colnames(res_order_per_chr) == "start")], #DMR start site
                        res_order_per_chr[unlist(intersection_matrix_indices[,1]), which(colnames(res_order_per_chr) == "end")], #DMR end site
                        res_order_per_chr[unlist(intersection_matrix_indices[,1]), which(colnames(res_order_per_chr) == "Disease.p.value")] #DMR disease p value
    )
    
    
    #name columns
    colnames(gene_check) <- c("Chr", "Searched_Gene_Start", "Searched_Gene_End", "Searched_Gene_Symbol", "Searched_Gene_Region_Type", "DMR_Start", "DMR_End", "DMR_p.value")
    
    #remove rows with NA gene symbol (depends on annotation reference used, does not seem to have any NAs for the hg38 genome)
    rm.id = which(is.na(gene_check[,which(colnames(gene_check) == "Searched_Gene_Symbol")]))
    if(length(rm.id)>0){
      gene_check <- gene_check[-rm.id,]
    }
    
    #name after extraupstream (not used in this comparitive analysis; thus 0bp is reported)
    assign(paste0("genes_alpha_", alpha_index, "_chr", chr), gene_check)
    
  }
  
  #save
  save.image(file = paste0("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Annotation/PT1_Results/PT1:Annotation:Chr", chr, ".RData"))
  
  print(paste0("done ", chr))
}
