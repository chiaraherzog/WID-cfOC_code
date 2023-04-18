# Extract methylation values from bismark output
dir <- "<bismark output dir from script one>"
bismark <- "<path to bismark>"  # Located in 0-data/bismark-0.22.3/
output <- "<path to extraction output>"
files <- list.files(dir, pattern = ".bam")
files <- paste(files, collapse = " ")

# Generate script
line1 <- paste0("cd ", dir)
write(line1, file = "1-preprocessing/3b-bismark-methylation-extr.sh",
      append = F)

line2 <- paste0(bismark, "/bismark_methylation_extractor --comprehensive --multicore 2 -o ", output, " ", files)
write(line2, file = "1-preprocessing/3b-bismark-methylation-extr.sh",
      append = T)

