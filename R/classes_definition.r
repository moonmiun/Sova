#' @import Biostrings
#' @import biomaRt
#' @import Gviz
#' @import IRanges
#' @importFrom S4Vectors Rle
#' @importFrom GenomicRanges GRanges
NULL

#' Class Gene
#'
#' A virtual class representing a generic gene.
#'
#' @slot id A unique identifier for the gene (ex, Ensembl ID 'ENSG00000139618').
#' @slot symbol Gene symbol (ex, HUGO symbol like 'BRCA2').
#' @slot name Complete name for the gene.
#' @slot description Gene description. Often the summary found on NCBI database.
#' @slot assembly Assembly used to store the position information (better to use
#' ENSEMBL IDs 'hg19' or 'hg38').
#' @slot chromosome Chromosome position.
#' @slot start Starting position (0 based).
#' @slot end Ending position.
#' @slot strand DNA strand ('+' or '-').
#' @slot sequence DNA sequence
#' @keywords classes
#' @export
setClass(
  Class = "Gene",
  slots = c(
    id = "character",
    symbol = "character",
    name = "character",
    description = "character",
    assembly = "character",
    chromosome = "character",
    start = "numeric",
    end = "numeric",
    strand = "character",
    sequence = "DNAString"
  ),
  prototype = list(
    id = NA_character_,
    symbol = NA_character_,
    name = NA_character_,
    description = NA_character_,
    assembly = NA_character_,
    chromosome = NA_character_,
    start = NA_integer_,
    end = NA_integer_,
    strand = NA_character_,
    sequence = Biostrings::DNAString("")
  )
)
setValidity("Gene", function(object) {
  required_fields <- c("id", "assembly", "chromosome", "start", "end", "strand")
  for (field in required_fields) {
    if (is.na(slot(object, field)) || (is.character(slot(object, field))
                                       && slot(object, field) == ""))
    {
      return(paste0("The '", field, "' slot is required and cannot be empty."))
    }
  }

  tryCatch({
    # The function used here is defined in the general_package_functions.r file.
    validateGeneInputs(object@start, object@end, object@strand, object@sequence)
  }, error = function(e) {
    return(e$message)
  })
  TRUE
})


#' Derived gene class: ProteinCodingGene
#'
#' This is a derived class specifically for protein coding genes.
#'
#' @slot transcript_id Transcript ID (either from Ensembl or RefSeq).
#' @slot exons GRanges object with exons positions.
#' @slot protein_id ID of the protein (best if taken from UniProt).
#' @slot protein_name Name of the protein (best if taken from UniProt).
#' @slot protein_sequence Aminoacid sequence.
#' @examples
#' # Example of creating a GRanges object for exons
#' library(GenomicRanges)
#' exons <- GRanges(
#'   seqnames = Rle(c("chr13")),
#'   ranges = IRanges(start = c(32315507,32316421,32319076,32325075,32326100,
#'   32326241,32326498,32329442,32330918,32332271,32336264,32344557,32346826,
#'   32354860,32356427,32357741,32362522,32363178,32370401,32370955,32376669,
#'   32379316,32379749,32380006,32394688,32396897,32398161),
#'                    end = c(32315667,32316527,32319325,32325184,32326150,
#'                    32326282,32326613,32329492,32331030,32333387,32341196,
#'                    32344653,32346896,32355288,32356609,32357929,32362693,
#'                    32363533,32370557,32371100,32376791,32379515,32379913,
#'                    32380145,32394933,32397044,32400268)),
#'   strand = Rle(c("+"))
#' )
#'
#' @export
setClass(
  Class = "ProteinCodingGene",
  contains = "Gene",
  slots = c(
    transcript_id = "character",
    exons = "GRanges",
    protein_id = "character",
    protein_name = "character",
    protein_sequence = "AAString"
  ),
  prototype = list(
    transcript_id = NA_character_,
    exons = GRanges(
      seqnames = Rle(),
      ranges = IRanges(),
      strand = Rle()
    ),
    protein_id = NA_character_,
    protein_name = NA_character_,
    protein_sequence = Biostrings::AAString("")
  )
)
setValidity("ProteinCodingGene", function(object) {
  required_fields <- c("exons", "protein_id")


  for (field in required_fields) {
    if (any(is.na(slot(object, field))) || (is.character(slot(object, field))
                                            && slot(object, field) == "")) {
      return(paste0("The '", field, "' slot is required and cannot be empty."))
    }
  }
  TRUE
})


#' Derived class: ncRNA (Non Coding RNA)
#'
#' A class to define non-coding RNA.
#'
#' @slot transcript_id Transcript ID.
#' @slot type Non coding RNA can be: 'lncRNA' (long non-coding RNA), 'siRNA'
#' (small interfering RNA) and, these ones of which there are also specific
#' classes available, 'miRNA' (microRNA), 'snRNA' (small nuclear RNA),
#' 'tRNA' (transfer RNA), 'rRNA' (ribosomal RNA)
#' @slot RNA_sequence RNA sequence
#' @export
setClass(
  Class = "ncRNA",
  contains = "Gene",
  slots = c(
    transcript_id = "character",
    type = "character",
    RNA_sequence = "RNAString"
  ),
  prototype = list(
    transcript_id = NA_character_,
    type = NA_character_,
    RNA_sequence = Biostrings::RNAString("")
  )
)
setValidity("ncRNA", function(object) {
  if (is.na(object@type) || (is.character(object@type) && object@type == "")) {
    return(paste0("The 'type' slot is required and cannot be empty."))
  }
  TRUE
})


#' Derived class: miRNA
#'
#' Class to define microRNA
#'
#' @slot target_genes vector of Ensembl ID of the gene target of this microRNA.
#' @slot seed Eight-base seed region of the miRNA.
#' @slot mature_sequence Mature miRNA.
#' @examples
#' # Seed region
#' The binding between miRNAs and target genes is not perfect across the whole
#' mature miRNA sequence: in mammals it is dominated by the so-called seed
#' region.
#' This seed region at the 5’ end of the mature miRNA consists of eight
#' nucleotides.
#'
#' # Mature sequence
#' The generation of miRNAs is a multistage process. Briefly, the mature ∼22 nt
#' miRNA sequence is embedded in one strand of an ∼33 bp double-stranded stem
#' characteristic of hairpin structures
#' in primary miRNA (pri-miRNA) transcripts produced by RNA polymerase II or
#' III. The miRNA must therefore be excised during its biogenesis to elicit gene
#'  silencing. To quickly find the mature sequence of the miRNA, the user can
#'  check this database: https://mirbase.org/.
#' @export
setClass(
  Class = "miRNA",
  contains = "ncRNA",
  slots = c(
    target_genes = "character",
    seed = "character",
    mature_sequence = "RNAString"
  ),
  prototype = list(
    target_genes = character(0),
    seed = NA_character_,
    mature_sequence = Biostrings::RNAString("")
  )
)
setValidity("miRNA", function(object) {
  if (length(object@target_genes) == 0 || all(is.na(object@target_genes)))
  {
    return("The 'target_genes' slot is required and cannot be empty.")
  }
  TRUE
})

#' Derived class: tRNA
#'
#' A class to represent tRNA
#'
#' @slot aminoacid Aminoacid transported by the tRNA.
#' @slot anticodon Anticodon sequence.
#' @examples
#' # Correct notation for aminoacids
#' Aminoacid information can be stored in any of the following notations:
#' - full name (ex: 'phenylalanine')
#' - 3 letters notation (ex: 'phe')
#' - 1 letter notation (ex: 'F')
#' None of these notations is required or checked by the package, but it is
#' recommended for the user to only use one of them for all of the objects.
#'
#' # Correct notation for anticodon sequences
#' The standard notation for anticodon sequences is 5'-XXX-3', using RNA
#' nitrogenous bases (which means Thymine is not allowed).
#' For the sake of simplicity, this package only stores information about the
#' basepairs and doesn't allow the user to indicate the 5' and 3' ends: it is
#' good practice, though, to remember to store the anticodon information
#' in the correct order.
#' Anti-codons' sequences can be found here: https://rnacentral.org/.
#' @export
setClass(
  Class = "tRNA",
  contains = "ncRNA",
  slots = c(
    aminoacid = "character",
    anticodon = "character"
  ),
  prototype = list(
    aminoacid = NA_character_,
    anticodon = NA_character_
  )
)
setValidity("tRNA", function(object) {
  if (is.na(object@aminoacid) || (is.character(object@aminoacid))
      && object@aminoacid == "") {
    return(paste0("The 'aminoacid' slot is required and cannot be empty."))
  }

  if (nchar(object@anticodon) != 3) {
    return("The anticodon must be exactly composed by 3 characters.")
  }
  TRUE
})


#' Class rRNA
#'
#' A class that represents genes that code for ribosomal rRNA (rRNA)
#'
#' @slot rRNA_type Type of rRNA subunit. The ones accepted are "18S", "28S",
#' "5S", "5.8S".
#' @slot ribosomal_subunit In which subunit of the ribosome the rRNA is found.
#' @examples
#' # Ribosomal subunit
#' A ribosome is composed by a small and a big subunit. For prokaryotes the
#' small subunit is denoted as 30S and the big one as 50S, while for eukaryotes
#' the former is 40S and the latter is 60S.
#' @export
setClass(
  Class = "rRNA",
  contains = "ncRNA",
  slots = c(
    rRNA_type = "character",
    ribosomal_subunit = "character"
  ),
  prototype = list(
    rRNA_type = NA_character_,
    ribosomal_subunit = NA_character_
  )
)
setValidity("rRNA", function(object) {
  valid_types <- c("18S", "28S", "5S", "5.8S")
  valid_subunits <- c("40S", "60S", "30S", "50S")
  if (!object@rRNA_type %in% valid_types || is.na(object@rRNA_type)) {
    return("Invalid rRNA type. Must be non-empty and one of: 18S, 28S, 5S,
           5.8S.")
  }
  if (!object@ribosomal_subunit %in% valid_subunits) {
    return("Invalid ribosomal subunit. Must be one of: 40S, 60S, 30S, 50S.")
  }
  TRUE
})


#' Class snRNA
#'
#' A class that represents genes that code for small nuclear RNA (snRNA)
#'
#' @slot snRNA_class A short string that defines a specific type of snRNA,
#' usually differentiated by function. (Ex: U1, U2...)
#' @slot subclass Subclass of the snRNA, defined by its localisation: either
#' 'snoRNA' or 'scaRNA'.
#' @slot associated_proteins A list of the ID of the proteins that the snRNA
#' interacts with.
#' @examples
#' # Subclasses
#' By definition, snRNAs are small molecules of RNA that can be found in the
#' nucleus. Two subclasses can be identified depending on the specific
#' localization inside the nucleus:
#' - snoRNA, is the snRNA that if located in the nucleous;
#' - scaRNA, is the snRNA located in the Cajal body.
#' @export
setClass(
  Class = "snRNA",
  contains = "ncRNA",
  slots = c(
    snRNA_class = "character",
    subclass = "character",
    associated_proteins = "character"
  ),
  prototype = list(
    snRNA_class = NA_character_,
    subclass = NA_character_,
    associated_proteins = character(0)
  )
)
setValidity("snRNA", function(object) {
  if (is.na(object@snRNA_class) || (is.character(object@snRNA_class)
                                    && object@snRNA_class == "")) {
    return(paste0("The 'snRNA_class' slot is required and cannot be empty."))
  }
  TRUE
})
