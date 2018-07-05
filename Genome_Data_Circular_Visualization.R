library(RCircos)
round_data <- read.table(file = "/data/Shaurya/RA@NBRI/RCircos_Histogram_Data.bed", header = TRUE)
names(round_data) = c("Chromosome", "ChromStartIndex", "ChromEndIndex", "Occupancy_Log2FC")
chr.exclude <- NULL
cyto.info <- read.table(file = "/data/Shaurya/RA@NBRI/Arabidopsis_thaliana_Chromosomes.txt", header = TRUE)
tracks.inside <- 10
tracks.outside <- 0
RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.cyto <- RCircos.Get.Plot.Ideogram()
rcircos.position <- RCircos.Get.Plot.Positions()
RCircos.List.Plot.Parameters()
rcircos.params <- RCircos.Get.Plot.Parameters()
rcircos.params$base.per.unit <- 3000
RCircos.Reset.Plot.Parameters(rcircos.params)
RCircos.List.Plot.Parameters()
out.file <- "RCircosDemoArabidopsisThalianaGenome.pdf"
pdf(file=out.file, height=8, width=8, compress=TRUE)
RCircos.Set.Plot.Area()

#*******************Try 2****************************************#

install.packages("colorspace")
install.packages("circlize")

library(RCircos)
library(circlize)
bed_file1 = read.table("/data/Shaurya/JCircosPlot/fig2_100/C_vs_CS_Depleted.bed", header = FALSE, sep = "\t")
bed_file2 = read.table("/data/Shaurya/JCircosPlot/fig2_100/C_vs_CS_Enriched.bed", header = FALSE, sep = "\t")
bed_file3 = read.table("/data/Shaurya/JCircosPlot/fig2_100/N_vs_NS_Depleted.bed", header = FALSE, sep = "\t")
bed_file4 = read.table("/data/Shaurya/JCircosPlot/fig2_100/N_vs_NS_Enriched.bed", header = FALSE, sep = "\t")

value1 <- bed_file1[4]
region1 <- bed_file1[,2:3]

value2 <- bed_file2[4]
region2 <- bed_file2[,2:3]

value3 <- bed_file3[4]
region3 <- bed_file3[,2:3]

value4 <- bed_file4[4]
region4 <- bed_file4[,2:3]

par(mar = c(1, 1, 1, 1))
circos.genomicInitialize(as.data.frame(cyto.info))
circos.genomicTrackPlotRegion(bed_file1, numeric.column = 4, track.height = uh (3, "mm"),
                    panel.fun = function(region1, value1, ...) {
                      circos.genomicLines(region1, value1, type = "h", numeric.column = 1, col =c("#FF0000"))
                    })
circos.genomicTrackPlotRegion(bed_file2, numeric.column = 4, track.height = uh (3, "mm"),
                    panel.fun = function(region2, value2, ...) {
                      circos.genomicLines(region2, value2, type = "h", numeric.column = 1, col =c("#008000"))
                    })
circos.genomicTrackPlotRegion(bed_file3, numeric.column = 4, track.height = uh (3, "mm"),
                    panel.fun = function(region3, value3, ...) {
                      circos.genomicLines(region3, value3, type = "h", numeric.column = 1, col =c("#FFCCE5"))
                    })
circos.genomicTrackPlotRegion(bed_file4, numeric.column = 4, track.height = uh (3, "mm"),
                    panel.fun = function(region4, value4, ...) {
                      circos.genomicLines(region4, value4, type = "h", numeric.column = 1, col =c("#00FFFF"))
                    })
legend("topleft", legend = c("C vs CS Depleted", "C vs CS Enriched", "N vs NS Depleted", "N vs NS Enriched"),
       col=c("#FF0000", "#008000", "#FFCCE5", "#00FFFF"), cex = 0.5, fill = c("#FF0000", "#008000", "#FFCCE5", "#00FFFF"))

