library(GenomicRanges)
library(IRanges)
library(Gviz)
library(Biostrings)
library(Sova)

# Create Gene test
gene <- createGene(
  id = "ENSG00000139618.19",
  symbol = "BRCA2",
  name = "BRCA2 DNA repair associated",
  description = "Plays a role in DNA repair",
  assembly = "hg38",
  chromosome = "chr13",
  start = 32315508,
  end = 32400268,
  strand = "+",
  sequence = "ACGTACCGTCAAA")

# Create Protein Coding Gene

exons <- GRanges(
  seqnames = Rle(c("chr2")),
  ranges = IRanges(end = c(130181602, 130176637, 130175000, 130173656, 130173354, 
                           130172895, 130172677, 130172500, 130167590, 130164445, 
                           130161272, 130157396, 130156675, 130156135, 130155259, 
                           130154482, 130153935, 130153450, 130153171, 130152884),
                   start = c(130181530, 130176554, 130174914, 130173514, 130173279, 
                             130172787, 130172625, 130172349, 130167458, 130164374, 
                             130161186, 130157251, 130156585, 130156035, 130155096, 
                             130154277, 130153702, 130153319, 130153043, 130151408)),
  strand = Rle(c("-")))

SMPD4 <- createProteinCodingGene(
  id = "ENSG00000136699.20",
  symbol = "SMPD4",
  name = "sphingomyelin phosphodiesterase 4, transcript variant 2",
  assembly = "hg38",
  chromosome = "chr2",
  start = 130151408,
  end = 130181602,
  strand = "-",
  transcript_id = "ENST00000680298.1",
  exons = exons, 
  protein_id = "Q9NXE4",
  protein_name = "Sphingomyelin phosphodiesterase 4")