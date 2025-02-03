library(testthat)
library(GenomicRanges)
library(IRanges)
library(Gviz)
library(Biostrings)
library(Sova)

# Create Gene
test_that("createGene creates a valid Gene objects", {
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
  expect_s4_class(gene, "Gene")
  expect_true(gene@start < gene@end)
  })


# Create Protein Coding Gene

test_that("createProteinCodingGene functions correctly", {
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
  expect_s4_class(SMPD4, "ProteinCodingGene")
})

# Create ncRNA
test_that("create_ncRNA functions correctly", {
  ncRNA_miRNA <- create_ncRNA(
    id = "ENSG00000208008",
    symbol = "MIR125A",
    name = "microRNA 125a",
    description = "MIR125A (MicroRNA 125a) is an RNA Gene, and is affiliated with the miRNA class. Diseases associated with MIR125A include Medulloblastoma and Lung Cancer.",
    assembly = "hg38",
    chromosome = "chr19",
    start = 51693254,
    end = 51693339,
    strand = "+",
    sequence = "ACGTACCGTCAAA",         # These sequences are just examples and not real
    transcript_id = "ENST00000385273.1",
    type = 'miRNA',
    RNA_sequence = "UGCCAGUCUC")

  expect_s4_class(ncRNA_miRNA, "ncRNA")
})

# Create tRNA
test_that("create_tRNA works correctly", {
ncRNA_tRNA_Phe <- create_tRNA(
  id = "ENSG00000210049.1",
  name = "Mitochondrial tRNA phenylalanine",
  assembly = "hg38",
  chromosome = "chrM",
  start = 577,
  end = 647,
  strand = "+",
  transcript_id = "ENST00000387314.1",
  type = "tRNA",
  aminoacid = "phenylalanine",
  anticodon = "GAA"
)
expect_s4_class(ncRNA_tRNA_Phe, "tRNA")
})

# create rRNA
test_that("create_rRNA works correctly", {
rRNA_18S <- create_rRNA(
  id = "ENSG00000225840",
  symbol = "RNA18SN5",
  name = "RNA, 18S ribosomal pseudogene",
  description = "45S ribosomal DNA (rDNA) arrays, or clusters, are present on human chromosomes 13, 14, 15...",
  assembly = "hg19",
  chromosome = "chrY",
  start = 10197256,
  end = 10199103,
  strand = "-",
  transcript_id = "ENST00000445125.2",
  type = "tRNA",
  rRNA_type = "18S",
  ribosomal_subunit = "40S"
)
expect_s4_class(rRNA_18S, "rRNA")
})

# create snRNA
test_that("create_snRNA works correctly", {
ncRNA_snRNA <- create_snRNA(
  id = "ENSG00000206652",
  symbol = "RNU1-1",
  name = "RNA, U1 small nuclear 1",
  description = "",
  assembly = "hg38",
  chromosome = "chr1",
  start = 16514122,
  end = 16514285,
  strand = "-",
  sequence = "tttcatacttacctggcagg...",
  RNA_sequence = 'aaaguaug...',
  transcript_id = "ENST00000383925.1",
  snRNA_class = 'U1')
expect_s4_class(ncRNA_snRNA, "snRNA")
})

# create miRNA
test_that("create_miRNA works correctly", {
miRNA125a <- create_miRNA(
  id = "ENSG00000208008",
  symbol = "MIR125A",
  name = "microRNA 125a",
  description = "MIR125A (MicroRNA 125a) is an RNA Gene, and is affiliated with the miRNA class. Diseases associated with MIR125A include Medulloblastoma and Lung Cancer.",
  assembly = "hg38",
  chromosome = "chr19",
  start = 51693254,
  end = 51693339,
  strand = "+",
  sequence = "ACGTACCGTCAAA...",
  transcript_id = "ENST00000385273.1",
  type = 'miRNA',
  RNA_sequence = "UGCCAGUCUC",
  target_genes = c("ENSG00000141736", "ENSG00000065361"),
  seed = 'GAGUCCC',
  mature_sequence = "ucccugagacccuuuaaccuguga")
expect_s4_class(miRNA125a, "miRNA")
})

# Transform class
test_that("transform class works correctly", {
ncRNA <- create_ncRNA(
  id = "ENSG00000210049.1",
  name = "Mitochondrial tRNA phenylalanine",
  assembly = "hg38",
  chromosome = "chrM",
  start = 577,
  end = 647,
  strand = "+",
  transcript_id = "ENST00000387314.1",
  type = "tRNA"
)
expect_s4_class(ncRNA, "ncRNA")

tRNA <- transformClass(ncRNA, "tRNA", list(aminoacid = "phenylalanine", anticodon = "GAA"))
expect_s4_class(tRNA, "tRNA")
})

# Get slot information
test_that("get_slot works correctly", {
  tRNA <- create_tRNA(
    id = "ENSG00000210049.1",
    name = "Mitochondrial tRNA phenylalanine",
    assembly = "hg38",
    chromosome = "chrM",
    start = 577,
    end = 647,
    strand = "+",
    transcript_id = "ENST00000387314.1",
    type = "tRNA",
    aminoacid = "phenylalanine",
    anticodon = "GAA"
  )
anticodon <- get_slot(tRNA, "anticodon")
expect_equal(anticodon, "GAA")
})

# Update slot information
test_that("set_slot works correctly", {
  tRNA <- create_tRNA(
    id = "ENSG00000210049.1",
    name = "Mitochondrial tRNA phenylalanine",
    assembly = "hg38",
    chromosome = "chrM",
    start = 577,
    end = 647,
    strand = "+",
    transcript_id = "ENST00000387314.1",
    type = "tRNA",
    aminoacid = "phenylalanine",
    anticodon = "GAA"
  )
assembly <- get_slot(tRNA, "assembly")

expect_equal(assembly, "hg38")

tRNA <- set_slot(tRNA, "assembly", "hg19")
new_assembly <- get_slot(tRNA, "assembly")

expect_equal(new_assembly, "hg19")
})


# Import from ENSEMBL
test_that("importGeneFromEsnembl works correctly", {
imported <- importGeneFromEnsembl("ENSG00000154975")
expect_s4_class(imported, "Gene")
})

# Visualize Gene
test_that("lengthProduct works correctly", {
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
  expect_silent(visualizeGene(SMPD4, ideogram = "yes"))
})

# Product length
test_that("lengthProduct works correctly", {
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

  len <- lengthProduct(SMPD4)
  expect_type(len, "double")
  expect_gt(len, 0)
})
