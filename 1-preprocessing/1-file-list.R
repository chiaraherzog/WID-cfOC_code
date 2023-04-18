# Identify directory for reads and concateante mates for bismark

dir <- "<your data dir>"
files <- list.files(dir)

# mate 1 carries "1_cut"
mates1 <- files[grepl("1_cut", files)]
mates1 <- paste(mates1, collapse = ",")

# mate 2 carries "1_cut"
mates2 <- files[grepl("2_cut", files)]
mates2 <- paste(mates2, collapse = ",")


# Generate script for bismark
bismark <- "<path to bismark>"  # Located in 0-data/bismark-0.22.3/
targets <- "<path to targets file>" # Located in 0-data/targets/
output <- "<path to output dir>" # should end with '/'

line1 <- paste0("cd ", dir)
write(line1, file = "1-preprocessing/1-bismark.sh",
      append = F)
line2 <- paste0(bismark, " --genome ", targets, " -o ", output, " --multicore 4 -1 ", mates1, " -2 ", mates2)
write(line2, file = "1-preprocessing/1-bismark.sh",
      append = T)


