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

#........................................................................................
# experiment to consensus together
read1 <- S[1]
read2 <- S[4]
read1 <- reverse(read1)

for (i in (1:width(read1))){
  cons1 <- read1[[1]][i:length(read1)]
  cons2 <- read2[[1]][1:i:-1]
  
  if (cons1 == cons2){
    print(cons1)
    print(cons2)
    number <- length(cons2)
  }
  print(number)
}

# function iterative function (correct overlap)
n1 <- width(read1)
n2 <- width(read2)

for (i in (1:min(n1, n2))){ #look for max consensus this is the min consensus length max
  cons1 <- subseq(read1, start = n1 - i + 1, end = n1) #start read1
  cons2 <- subseq(read2, start = 1, end = i)  #end read2
  
  if (as.character(cons1) == as.character(cons2)){
    print(cons1)
    print(cons2)
    number <- i
  }
  print(number)
}

read1
read2

read1[[1]][1] # [[]]go to seq and index []


# matrix create experiment
n <- length(S)
matrix_overlap <- matrix(0,n,n, dimnames = list(as.character(S), as.character(S)))
matrix_overlap

# second result with properties to namec col and row
#matrix <- matrix(0, n, n)
#names <- c('CATGC','CTAAGT', 'GCTA', 'TTCA','ATGCATC')
#rownames(matrix) <- names
#colnames(matrix) <- names
#matrix
n
for (row in (1:n)){
  for (col in (1:n)){
    matrix_overlap[row,col] <- Overlap(S[row], S[col])
  }
}
matrix_overlap
matrix_overlap[row,col] <- Overlap(S)
diag(matrix_overlap) <- -Inf  # nebo NA
matrix_overlap


#experiment with greedy
max_overlap <- max(matrix_overlap)
max_overlap

idx <- which((matrix_overlap == max(matrix_overlap)), arr.ind = TRUE)[1, ]
idx
seq1 <- S[idx[1]] #row
seq2 <- S[idx[2]] # col

seq1
seq2
# merge this 2 sequences
new_seq <- xscat(subseq(seq1, 1, width(seq1)-max_overlap), seq2) # seq1 without overlap # all seq2
new_seq
S
S <- c(S,new_seq)
S <- S[-c(idx[1], idx[2])]

 

# .......................all final function..............................
Overlap <- function(read1, read2){
  n1 <- width(read1)
  n2 <- width(read2)
  number <- 0
  
  for (i in (1:min(n1, n2))){ #look for max consensus this is the min consensus length max
    cons1 <- subseq(read1, start = n1 - i + 1, end = n1) #start read1
    cons2 <- subseq(read2, start = 1, end = i)  #end read2
    
    if (as.character(cons1) == as.character(cons2)){
      #print(cons1)
      #print(cons2)
      number <- i
    }
  }
  return(number)
}

Overlap(S[1],S[5])


OverlapMatrix <- function(S){
  n <- length(S)
  matrix_overlap <- matrix(0,n,n, dimnames = list(as.character(S), as.character(S))) #create matrix with 0

  for (row in (1:n)){
    for (col in (1:n)){
      matrix_overlap[row,col] <- Overlap(S[row], S[col]) #get the number of consensus
    }
  }
  diag(matrix_overlap) <- -Inf  # nebo NA
  return(matrix_overlap)
}
OverlapMatrix(S)

GreedySuperstring <- function(S){
  while (length(S)> 1){
    overlapMat <- OverlapMatrix(S)
    if (max(overlapMat) == 0){
      return (S)
    }
    else {
      max_overlap <- max(overlapMat) # max score from 2 seq
      idx <- which((overlapMat == max(overlapMat)), arr.ind = TRUE)[1, ] # found idx with max score in matrix
      #get> row col
      #      1    1
      seq1 <- S[idx[1]] #row
      seq2 <- S[idx[2]] # col
      
      # merge this 2 sequences
      new_seq <- xscat(subseq(seq1, 1, width(seq1)-max_overlap), seq2) # seq1 without overlap # all seq2
      new_seq
      
      # add new consensus
      S <- c(S,new_seq)
      
      #remove old consensus
      S <- S[-c(idx[1], idx[2])]
    }
  }
  return(S)
}
library(Biostrings)
GreedySuperstring(S <- DNAStringSet(c('CATGC','CTAAGT', 'GCTA', 'TTCA','ATGCATC')))
  
  
  
  