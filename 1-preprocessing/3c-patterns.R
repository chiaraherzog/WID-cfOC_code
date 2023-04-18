# Extract patterns
input <- "<output path from 3b>"
output <- "<path to store output of this script>"
targets <- "<path to targets>"

patterns <- "<path to count-patterns.R file>"
source(patterns)
count_patterns(input,
               output,
               targets)

