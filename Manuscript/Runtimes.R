#chromosomes involved in this script
chromosomes <- c(1:22,"X")


load("~/scratch/Marie_Hudson_Scleroderma_WGBS/Partitioning_Regions_Spacing/PT3:Constructing_Spacing_DataFrame:GENOME.RData")

plot_me <- matrix(data = NA, nrow = length(chromosomes) + 1, ncol = 2)
colnames(plot_me) <- c("Region", "Runtime")


pdf(paste0("/home/greenwood/jeffrey.yu/scratch/Marie_Hudson_Scleroderma_WGBS/Manuscript/Plots/Runtimes.pdf"), width = 10, height = 7)

for (chr in chromosomes){

  plot_me[grep(chr, chromosomes),"Runtime"] <- get(paste0("somnibus_time_chr", chr))
}

plot_me[,"Region"] <- c(paste0("Chr. ", chromosomes), "Genome")
plot_me[24,"Runtime"] <- get(paste0("somnibus_time_genome"))

plot(plot_me[,"Runtime"]~factor(plot_me[,"Region"]), xaxt="n", type = "p", xlab = "", ylab = "Runtime (minutes)", main = "SOMNiBUS DMR Detection Runtimes")
axis(1, labels = plot_me[,"Region"], at = 1:24, las=2)

dev.off()
