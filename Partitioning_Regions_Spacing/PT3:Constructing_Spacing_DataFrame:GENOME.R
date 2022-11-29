#=====
#THIS IS JUST A FILE TO BIND THE DATA FOR ALL CHROMOSOMES TOGETHER AND TO CALCULATE TOTAL RUNTIMES FOR ALL CHROMOSOMES
#chromosomes involved in this scrpt
chromosomes = c(1:22,"X")
#=====

for (chr in chromosomes){
  
  #load shortened data frame of CpG regions analyzed by SOMNiBUS, along with statistics
  load(paste0("/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Partitioning_Regions_Spacing/PT2_Results/PT2:Spacing_Info:Chr", chr, ".RData"))
  
  #append chromosome to the end of imported objects
  assign(paste0("Reg_res_chr", chr), Reg_res)
  # assign(paste0("sub_groups_", chr), sub_groups)
  assign(paste0("somnibus_time_chr", chr), somnibus_time)
  #NOTE: TIME WAS RUN ON 1 NODE WITH 8 PROCESSORS FOR ALL 23 CHROMOSOMES IN SEQUELTIAL ORDER BUT USED A FOR(CHR IN CHROMOSOMES) LOOP SO SHOULD *NOT* BE REGARDED AS "RUN IN PARRALELE" FOR TEH TIEM CALCULATE
}


#sum individual runtimes to get the genome runtime
#remember to also keep individual runtimes (do not delete: somnibus_time_chr1, somnibus_time_chr2, etc.)
somnibus_time_genome <- sum(unlist(mget(paste0("somnibus_time_chr", chromosomes))))

#create master table with all CpG regions information (nCpG, chromosome, end position, etc.) and staistics (p.values, phi, EDFs, etc.), for all chromosomes
Reg_RES_ALL_cg <- do.call(rbind, mget(paste0("Reg_res_chr", chromosomes)))

#save
save.image(file = "/home/jyu/scratch/Marie_Hudson_Scleroderma_WGBS/Partitioning_Regions_Spacing/PT3_Results/PT3:Constructing_Spacing_DataFrame:GENOME.RData")
