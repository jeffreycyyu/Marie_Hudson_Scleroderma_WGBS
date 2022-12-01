require(lattice)
library(ggplot2)
library(annotatr)
library(dplyr)

#chromosomes involved in this script
chromosomes <- c(1:22,"X")

annotation_types_used <- head(builtin_annotations()[grepl("hg19", builtin_annotations())], n=11)[c(1:4, 6, 5, 7, 10)]
annotation_types_names <- c("1-5kb_upstream_TSS", "promoter", "coding_region", paste0("5", "'", "UTR"), "first_exon", "exon", "intron", paste0("3", "'", "UTR"))
annotation_label_colours <- c("yellow", "green", rep("blue",2), "green", "blue", "dodgerblue", "blue")


load("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Manuscript/Bumphunter/Bumphunter_Summary.RData")
master_list_genes <- intersected_gene_list

#chromosomes need to be respecified (unless I clean up my "save" functions in the previous parts)
chromosomes <- c(1:22,"X")



#load all files from PT1
# for (chr in chromosomes){
#  load(paste0("/scratch/greenwood/jeffrey.yu/Marie_Hudson_Scleroderma_WGBS/Annotation/PT1_Results/PT1:Annotation:Chr", chr, ".RData"))
# }
load(paste0("/scratch/jyu//Marie_Hudson_Scleroderma_WGBS/Annotation/PT1_Results/PT1:Annotation:ChrX.RData"))
#load annotation data
#this is annotation table with all information on DMRs and thier annotations (ie. all results are contained in this data table)
#we regenerate this because now it is ordered by position (both from chr 1 to chr X; and ordered by position within each chromosome)
#the plot pdfs generated later, thus, should be in the order they would be found on the genome

#chromosomes need to be respecified (unless I clean up my "save" functions in the previous parts)
chromosomes <- c(1:22,"X")

#list of genes from Bonferoni primary analysis
master_annotation_table <- do.call(rbind, mget(paste0("genes_alpha_2_chr", chromosomes)))





for(gene in master_list_genes){
  print(gene)
  print(as.character(unique(master_annotation_table[master_annotation_table[,"Searched_Gene_Symbol"] == gene, "Chr"])))
}





#for all genes detected in "primary analysis"
#(note: we want only primary analysis beta(t) values since those genes ar elinmited to only gene sof very high significance, the secondary analysis was used for IPA since it was the most comprehensive.)
#(note: having the most comprehensive set of genes detected does not matter as much for beta(t) since we are plotting each gene seperately so we focus on getting only "good" graphs (i.e. "most significant" genes) which is the primary analysis since taht is more restrictive)
gene="SMAD3"

try({
  #all DMRs associated with this gene's start sites
  special_DMR_start <- as.numeric(as.character(unique(master_annotation_table[master_annotation_table[,"Searched_Gene_Symbol"] == gene, "DMR_Start"])))
  #all DMRs associated with this gene's end sites
  special_DMR_end <- as.numeric(as.character(unique(master_annotation_table[master_annotation_table[,"Searched_Gene_Symbol"] == gene, "DMR_End"])))
  #create matrix where rows are DMRs (may only have one row/DMR) and columns indicate DMR start and end sites
  tally_special <- cbind(special_DMR_start, special_DMR_end)
  
  #sanity: reassign chromosome
  chr <- as.character(unique(master_annotation_table[master_annotation_table[,"Searched_Gene_Symbol"] == gene, "Chr"]))
  
  #for each DMR of this gene...
  for (j in 1:nrow(tally_special)){
    
    #load SOMNiBUS analysis results
    load(paste0("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Partitioning_Regions_Spacing/PT2_Results/PT2:Spacing_Info:Chr", chr, ".RData"))
    
    #sequence of all positions within this DMR
    part_elig_pos = seq(tally_special[j,1], tally_special[j,2])
    #intersect all annotated positions with this DMR's positions
    elig_pos <- intersect(part_elig_pos, position)
    
    #start position of the DMR
    i <- as.numeric(tally_special[j, 1])
    
    #sanity check: checks to make sure ampping went correctly and no DMR invades 2 choromosomes
    if (length(which(sub_groups$start == i)) == 1) {
      
      #CpG region/subgroup names (this will be used to filter our origional seuqencing results to only include information relevant to this particular CpG region)
      meth_convert_index <- rownames(sub_groups[as.numeric(which(sub_groups$start == i)),])
      #all CpG nucleotides that are part of that CpG region/subgroup
      methMat_now <- methMat[which(groups == meth_convert_index),]
      
      #get the dinucleotide positions within this CpG region
      naming <- as.character(rownames(methMat_now))
      tmp1 <-as.character(lapply(strsplit(naming,"-"),function(x)x[1]) %>% unlist)
      pos_now <- lapply(strsplit(tmp1,":"),function(x)x[2]) %>% unlist %>% as.numeric
      
      #this is the full beta values
      beta_out <- RES[[which(sub_groups$start == i)]][["out_re_quasi"]][["Beta.out"]]
      #this is the full standard error of betas
      beta_se_out <- RES[[which(sub_groups$start == i)]][["out_re_quasi"]][["SE.out"]]
      
      mean_age = ((52.8*9) + (37.2*4))/(9+4)
      beta_out2 = cbind(exp(beta_out[,"Intercept"] + beta_out[,"Age"]*mean_age + beta_out[,"Disease"])/(1+exp(beta_out[,"Intercept"] + beta_out[,"Age"]*mean_age + beta_out[,"Disease"])),
                        exp(beta_out[,"Intercept"] + beta_out[,"Age"]*mean_age)/(1+exp(beta_out[,"Intercept"] + beta_out[,"Age"]*mean_age)))
      colnames(beta_out2) = c("SSc", "CNTL")
      #methylation proportions
      print(mean(beta_out2[elig_pos %in% bumphunter_bump_sequences,"SSc"]))
      print(mean(beta_out2[elig_pos %in% bumphunter_bump_sequences,"CNTL"]))
      
      
      #this is the disease p-value for the DMR being plotted
      disease_pvalue <- RES[[which(sub_groups$start == i)]]$out_re_quasi$reg.out["Disease", "p-value"]
      
      
      
      #rename master annotation set object to not contain chromosoem number
      assign("genes_alpha_2", as.data.frame(unique(get(paste0("genes_alpha_2_chr", chr)))))
      
      
      #this is a matrix where each row represents a unique annotated gene region with columns represneting the annotated gene region's start, stop, and type
      #each one of thse regions will get thier own coloured rectangle on the final plot spaning from thier start to thier end and labeled with thier region type
      region_types <-unique(genes_alpha_2[grep(i, genes_alpha_2[,"DMR_Start"]), c("Searched_Gene_Region_Type", "Searched_Gene_Start", "Searched_Gene_End", "DMR_Start", "DMR_End", "Searched_Gene_Symbol")])
      
      #retain only annotated regions that correspond to this gene (recall: DMRs can be annotated to multiple genes and we do not want to report annotated regions from other genes in this gene plot)
      region_types <- region_types[which(region_types[, "Searched_Gene_Symbol"] == gene), 1:(ncol(region_types)-1)]
      
      #this is the range of beta values divided into n equidistance parts for n types of annotation we are using
      y_sections <- (max(beta_out[,"Disease"])-min(beta_out[,"Disease"]))/length(annotation_types_names)
      
      #these are all CpG regions after filtering for bumphunter
      bumphunter_full_CpG_regions_per_chr <- bumphunter_full_CpG_regions[which(bumphunter_full_CpG_regions[,"chr"] == paste0("chr", chr)),]
      #only keep DMRs
      bumphunter_full_CpG_regions_per_chr <- bumphunter_full_CpG_regions_per_chr[which(bumphunter_full_CpG_regions_per_chr[, "p.value"] <= 0.05), ]
      #this is the row index in bumphunter_full_CpG_regions_per_chr that corresponds to the CpG regions after filtering of bumphunter overlapping with the DMR of somnibus
      bumphunter_match_index <- which(bumphunter_full_CpG_regions_per_chr[,"start"] >= pos_now[1] & bumphunter_full_CpG_regions_per_chr[,"end"] <= pos_now[length(pos_now)])
      
      
      
      if (length(bumphunter_match_index) >= 1){
        
        
        for (bump_index in 1:length(bumphunter_match_index)){
          #these are positions to plot (overlapping positions, not including posistions non-overlapping with somnibus)
          #note: upstream 1bp for bumphunter DMR since CpG sites in datmeth are coded as is 2 nucleotides and bumphunter is assigning single nucleotide values, this is an intersect() so no new positions are introduced, only sanity coverage for genes such as AP2AT which had an error origionally
          bumphunter_bump_sequences <- intersect(pos_now, seq(bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"start"]-1, bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"end"]))
          
          
          
          #create plot on a pdf
          pdf(paste0("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/final_figures/proportion_SMAD3_Overlapping_Bumphunter_Beta_Values:Chr_", chr, "_Gene_", gene, "_DMR_", j, ".pdf"), width = 10, height = 7)
          #margins
          par(mar=c(6,6,6,14))
          
          #plot in small points since will be replotted to cover error bars
          plot(elig_pos, beta_out[,"Disease"], col = "black",
               pch = 20, cex = 0.75, cex.main = 1.0, cex.axis = 1.25, cex.lab = 1.25,
               xlab = "Position on Chromosome",
               ylab = "Logit[Estimated Disease Effect on Methylation] (i.e., beta values)",
               ylim = c(min(beta_out[,"Disease"] - beta_se_out[,"Disease"], bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"value"]), max(beta_out[,"Disease"] + beta_se_out[,"Disease"], bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"value"] + range(beta_out[,"Disease"] - beta_se_out[,"Disease"])/15)))
          #line at beta(t) = 0
          abline(h = 0, lty = 4)
          #do not start a newplot
          par(new = TRUE)
          #second plot
          plot(elig_pos, beta_out[,"Disease"], pch = 10, cex = 0.25, cex.main = 1.0, cex.axis = 1.25, cex.lab = 1.25, axes = FALSE, bty = "n", xlab = "", ylab = "", ylim = c(min(beta_out[,"Disease"] - beta_se_out[,"Disease"], bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"value"]), max(beta_out[,"Disease"] + beta_se_out[,"Disease"], bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"value"] + range(beta_out[,"Disease"] - beta_se_out[,"Disease"])/15)))
          #add standard errors
          arrows(x0 = elig_pos, y0 = beta_out[,"Disease"] - beta_se_out[,"Disease"], x1 = elig_pos, y1 = beta_out[,"Disease"] + beta_se_out[,"Disease"], code = 3, angle = 90, length = 0.05, col = "#DCDCDC")
          #plot the axis labels
          axis(side=4, at = rep(min(beta_out[,"Disease"])-0.5*y_sections , length(annotation_types_names)) + seq(1:length(annotation_types_names))*y_sections, cex.axis = 1.25, labels = annotation_types_names, las=1)
          #do not start a newplot for replotting of points
          
          par(new = TRUE)
          #BUMPHUNTER OVERLAPPING POINTS in red
          #replot points to front to cover grey error bar
          plot(elig_pos, beta_out[,"Disease"], col = ifelse(elig_pos %in% bumphunter_bump_sequences, "#DA0000", "black"),
               pch = 20, cex = 0.75, cex.main = 1.0, cex.axis = 1.25, cex.lab = 1.25,
               xlab = "Position on Chromosome",
               ylab = "Logit[Estimated Disease Effect on Methylation] (i.e., beta values)",
               ylim = c(min(beta_out[,"Disease"] - beta_se_out[,"Disease"], bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"value"]), max(beta_out[,"Disease"] + beta_se_out[,"Disease"], bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"value"] + range(beta_out[,"Disease"] - beta_se_out[,"Disease"])/15)))
          
          
          
          
          
          
          #for each distinct gene region annotated...
          for (k in 1:nrow(region_types)){
            try({
              #re-state margins because rectangles wont fit otherwise for some reason
              
              #the index for what type of annotated region we are dealing with for this rectangle
              region_type_cat <- grep(region_types[k, "Searched_Gene_Region_Type"], annotation_types_used)
              
              
              par(new = TRUE)
              #add in coloured rectangles to indicate region types annotated
              rect(as.numeric(as.character(region_types[k, "Searched_Gene_Start"])), #left edge
                   min(beta_out[,"Disease"]) + (region_type_cat-1)*y_sections, #bottom edge
                   as.numeric(as.character(region_types[k, "Searched_Gene_End"])), #right edge
                   min(beta_out[,"Disease"]) + region_type_cat*y_sections, #top edge
                   lty = 0, #no box lines
                   col= alpha(annotation_label_colours[region_type_cat], 0.35) # colour of rectangle fill
              )
              
              #this is the center of the rectangle generated
              #determine where the text should lie on the x axis, we want it in the center of the rectangle unless it is outside the graph
              proposed_text_loc <- mean(c(as.numeric(as.character(region_types[k, "Searched_Gene_Start"])), as.numeric(as.character(region_types[k, "Searched_Gene_End"]))))
              if (proposed_text_loc < as.numeric(region_types[k, "DMR_Start"])){proposed_text_loc <- as.numeric(region_types[k, "DMR_Start"])}
              if (proposed_text_loc > as.numeric(region_types[k, "DMR_End"])){proposed_text_loc <- as.numeric(region_types[k, "DMR_End"])}
              
              #generate text label for each box
              text(proposed_text_loc, min(beta_out[,"Disease"]) + (region_type_cat-0.5)*y_sections, labels = annotation_types_names[region_type_cat], cex = 0.5)
              
            })
            
            par(new = TRUE)
            #set new margins for y label since we made custom y ticks
            par(mar=c(6,6,6,14))
            
            #right axis title
            mtext("Annotated Gene Region", side=4, line=12, cex = 1.5)
            
            par(new = TRUE)
            #set new margins for y label since we made custom y ticks
            par(mar=c(6,6,6,14))
            
            #add rectangle displaying where bumphunter DMR is and the average methylation level difference of the DMR
            arrows(x0 = bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"start"]-1, y0 = bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"value"], x1 = bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"end"], y1 = bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"value"], code = 3, angle = 90, length = abs(max(beta_out[,"Disease"])-min(beta_out[,"Disease"]))/30, col = "#DA0000")
            #generate text label for bumphunter DMR
            text(mean(bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"start"]-1, bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index], "end"]), bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"value"] + abs(max(beta_out[,"Disease"])-min(beta_out[,"Disease"]))/10, labels = paste0("bumphunter DMR\n(mean beta: ", signif(bumphunter_full_CpG_regions_per_chr[bumphunter_match_index[bump_index],"value"], digits = 3), ")"), cex = 0.65, col = "#DA0000")
            #legend
            legend("topright", c("Not Identified by bumphunter", "Identified by bumphunter"), col=c("black", "#DA0000"), lwd=4, cex = 1, bty = "o")
            
            
            
          }
          #do not delete me, used at end of plot
          dev.off()
          
          
        }} else {print("Error_Jeff: BUMPHUNTER HAS NO OVERLAPPING DMRS WITH SOMNIBUS")}
    } else {print("Error_Jeff: The associated DMR for the gene does not seem to be located on the genes chromosome")}
    
  }
  
})
