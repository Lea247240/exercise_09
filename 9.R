rm(list=ls())

#Cesta:
setwd('V:/MPA_PRG/exercise_09')
setwd('D:/VUT/4-5rocnik/moje/MPA-PRG/exercise_09') # home

# *De Novo* Genome Assembly

## The Greedy Shortest Common Superstring
### Task
# In R, implement a function `GreedySuperstring()` according to the pseudocode.

# Input:
  # `S` A `DNAStringSet` object of strings (reads).

# Output:
  # `S` A `DNAStringSet` object of the shortest common superstring (contig).

# **Hint:** 
#Create also functions:
# `Overlap()` to calculate overlap between two sequences.
#`OverlapMatrix()` to create a matrix of overlaps among all sequences in `S`.

####################################################################
#GreedySuperstring(S)
#1   while length of S > 1
#2     overlapMat <- OverlapMatrix(S)
#3     if max(overlapMat) = 0
#4       return S
#5     else
#6       seq1, seq2 ← Two sequences from S with the longest overlap
#7       Merge seq1 and seq2 and add the new sequence to S
#8       Remove seq1 and seq2 from S
#9   return S
####################################################################
library(Biostrings)
s1 <- DNAStringSet('CATGC')
s2 <-  DNAStringSet('CTAAGT')
S <- c(s1,s2)

# input data
S <- DNAStringSet(c('CATGC','CTAAGT', 'GCTA', 'TTCA','ATGCATC'))

S <- c(DNAStringSet('CATGC'), DNAStringSet('CTAAGT'),DNAStringSet('GCTA'), DNAStringSet('TTCA'), DNAStringSet('ATGCATC'))
S


# experiment to consensus together
read1 <- S[1]
read2 <- S[2]
read1 <- reverse(read1)
read2 <- (read2)
max_number = 0
for (i in (1:width(read1))){
  cons1 <- read1[[1]][i:length(read1)]
  cons2 <- read2[[1]][1:i]
  
  if (cons1 == cons2){
    number <- length(cons2)
    if (number > max_number){
      max_number <- 0
    }
  }

}

read1
read2
read2_rev <- reverse(read2)
read2_rev

read1[[1]][1] # [[]]go to seq and index []


# matrix
n <- length(S)
matrix <- matrix(0, n, n)

names <- c('CATGC','CTAAGT', 'GCTA', 'TTCA','ATGCATC')
rownames(matrix) <- names
colnames(matrix) <- names
matrix

matrix[1,2] <- 1

matrix[row,col] <- Overlap(S)



# function..............................
Overlap <- function(S){
  
}



OverlapMatrix <- function(S){
  
}


GreedySuperstring <- function(S){
  while (length(S)> 1){
    overlapMat <- OverlapMatrix(S)
    if (max(overlapMat) == 0){
      return (S)
    }
    else {
      ...
    }
  }
}
  
  
  
  
  