# PART 1 — XCMS PROCESSING
# Step 1 — Navigate to Supernatant Folder


library(xcms)
library(CAMERA)

# Set working directory to your Supernatant folder
setwd("C:/Users/WiNDOWS/Desktop/Metabolomics-analysis/Positive/Supernatant")

# Confirm location
getwd()

# Confirm three folders visible
list.files()

# Confirm 6 files inside each folder
list.files("BNC5")
list.files("BNC5RSw")
list.files("RSw")



# Step 2 — Peak Detection

cat("Starting supernatant peak detection:", 
   format(Sys.time(), "%H:%M:%S"), "\n")

xs_sup <- xcmsSet()

cat("Peak detection complete:", 
    format(Sys.time(), "%H:%M:%S"), "\n")

# Confirm immediately
cat("Groups detected:\n")
print(xs_sup@phenoData$class)
cat("Total samples:", length(xs_sup@filepaths), "\n")
cat("Total peaks detected:", nrow(xs_sup@peaks), "\n")



# Step 3 — Round 1 Grouping and RT Correction

cat("Round 1 grouping...\n")
gxs_sup <- group(xs_sup)
cat("Round 1 grouping complete\n")

cat("Round 1 RT correction (span=0.6)...\n")
xs_sup2 <- retcor(gxs_sup, family = "s", 
                  plottype = "m", span = 0.6)
cat("Round 1 RT correction complete\n")
cat("Inspect plot — note y-axis deviation range\n")



# Step 4 — Round 2 Grouping and RT Correction

cat("Round 2 grouping...\n")
gxs_sup2 <- group(xs_sup2)
cat("Round 2 grouping complete\n")

cat("Round 2 RT correction (span=0.6)...\n")
xs_sup3 <- retcor(gxs_sup2, family = "s", 
                  plottype = "m", span = 0.6)
cat("Round 2 RT correction complete\n")
cat("Inspect plot — deviation should be smaller than Round 1\n")



# Step 5 — Final Grouping and Fill Peaks

cat("Final grouping...\n")
gxs_sup3 <- group(xs_sup3, bw = 10)
cat("Feature groups:", nrow(groups(gxs_sup3)), "\n")

cat("Filling missing peaks...\n")
fpgxs_sup3 <- fillPeaks(gxs_sup3)
cat("Features after filling:", nrow(groups(fpgxs_sup3)), "\n")



# Step 6 — CAMERA Annotation

cat("Starting CAMERA annotation...\n")

xsa_sup   <- xsAnnotate(fpgxs_sup3)
xsaF_sup  <- groupFWHM(xsa_sup, perfwhm = 0.6)
xsaC_sup  <- groupCorr(xsaF_sup)
xsaFI_sup <- findIsotopes(xsaC_sup)
xsaFA_sup <- findAdducts(xsaFI_sup, polarity = "positive")

cat("CAMERA annotation complete\n")



# Step 7 — Export Results

# Full CAMERA annotation file
write.csv(getPeaklist(xsaFA_sup),
          file = "BNC5_supernatant_pos_XCMS_CAMERA_full.csv")

# Peak intensity matrix
peakmatrix_sup <- groupval(fpgxs_sup3, value = "into")
rownames(peakmatrix_sup) <- groupnames(fpgxs_sup3, 
                                       mzdec = 4, 
                                       rtdec = 2)

write.csv(peakmatrix_sup,
          file      = "BNC5_supernatant_pos_peakmatrix.csv",
          row.names = TRUE)

cat("Two output files saved:\n")
cat("1. BNC5_supernatant_pos_XCMS_CAMERA_full.csv\n")
cat("2. BNC5_supernatant_pos_peakmatrix.csv\n")



# Step 8 — Save XCMS Session

save(xs_sup, gxs_sup, xs_sup2, gxs_sup2, xs_sup3,
     gxs_sup3, fpgxs_sup3, xsaFA_sup, peakmatrix_sup,
     file = "BNC5_supernatant_pos_XCMS_session.RData")

cat("XCMS session saved\n")
cat("Reload with: load('BNC5_supernatant_pos_XCMS_session.RData')\n")
