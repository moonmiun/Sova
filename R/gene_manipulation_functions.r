#' @import Biostrings
#' @import biomaRt
#' @import Gviz
#' @import IRanges
#' @importFrom S4Vectors Rle
#' @importFrom GenomicRanges GRanges
NULL

#' Create a Gene object
#'
#' This function creates an instance of the Gene class.
#'
#' @param id Gene ID.
#' @param symbol Gene symbol.
#' @param name Gene name.
#' @param description Description of the gene.
#' @param assembly Assembly used to store the position information.
#' @param chromosome Chromosome on which the gene is located.
#' @param start Start position of the gene on the chromosome (1-based).
#' @param end End position of the gene on the chromosome.
#' @param strand Strand of the gene.
#' @param sequence DNA sequence of the gene
#' @return A Gene object.
#' @examples
#' # 0 and 1 based positions
#' # This package uses 1-based coordinate system since ENSEMBL uses the same
#' # coordinate system. The user is advised to always check the documentation of
#' # any database the information is collected from.
#' # For example, UCSC Genome Browser tables use 0-based coordinate system, which
#' # means that the first base pair is at position 0 and not 1. At the same time,
#' # the UCSC Genome Browser web interface uses the 1-start coordinate system.
#'
#' # Creating a Gene object
#' # gene <- createGene(
#' #  id = "ENSG00000139618.19",
#' #    symbol = "BRCA2",
#' #    name = "BRCA2 DNA repair associated, transcript variant 4",
#' #   description = "Inherited mutations in BRCA1 and this gene, BRCA2, confer
#' #   increased lifetime risk of developing breast or ovarian cancer....",
#' #   assembly = "GRCh38.p14"
#' #   chromosome = "chr13",
#' #   start = 32315508,
#' #   end = 32400268,
#' #   strand = "+",
#' #   sequence = "ACGTACCGTCAAA....."
#' # )
#'
#' @export
createGene <- function(id, symbol = NA_character_, name = NA_character_,
                       description = NA_character_, assembly,
                       chromosome, start, end, strand,
                       sequence = Biostrings::DNAString("")) {
  sequence <- validateGeneInputs(start, end, strand, sequence)
  new("Gene",
      id = id,
      symbol = symbol,
      name = name,
      description = description,
      assembly = assembly,
      chromosome = chromosome,
      start = start,
      end = end,
      strand = strand,
      sequence = sequence)
}


#' Transform an object of a class, in an object of a subclass.
#'
#' This function lets the user create an object of a subclass by simply adding
#' the required parameters to the 'superior' class.
#'
#' @param object The object of the superior class (either 'Gene',
#' 'ProteinCodingGene' or 'ncRNA').
#' @param target_class The target subclass.
#' @param additional_fields A list of the additional fields to add to define a
#' subclass.
#' @return A Gene object.
#' @examples
#'  # Classes hierarchy
#' # It is possible to go from Gene to ProteinCodingGene or ncRNA, from ncRNA to
#' # miRNA, rRNA, tRNA or snRNA and from Gene to one of the subclasses of ncRNA.
#' # When jumping from Gene to one of the subclasses of ncRNA (skipping the
#' # intermediate step on ncRNA), it is important to remember the necessary field
#' # required by ncRNA ("type") too, in addition to subclasses-specific fields.
#'
#' # Transform a general ncRNA object in a tRNA object
#' # my_ncRNA <- create_ncRNA(id ='ENSG00000210049.1', name = 'mitochondrially
#' # encoded tRNA phenylalanine', assembly = 'GRCh38/hg38', chromosome = 'chrM',
#' # start = 577, end = 647, strand = '+', transcript_id = 'ENST00000387314.1',
#' # type = 'tRNA')
#' # my_tRNA <- transformClass(my_ncRNA, 'tRNA', list(aminoacid = 'phenylalanine',
#' # anticodon = 'GAA'))
#'
#' @export
transformClass <- function(object, target_class, additional_fields = list()) {
  if (!is(object, "Gene")) {
    stop("The object must inherit from the 'Gene' class.")
  }

  # Definying class hierarchy
  valid_classes <- c("Gene", "ProteinCodingGene", "ncRNA", "miRNA", "rRNA",
                     "tRNA", "snRNA")
  class_hierarchy <- list(
    Gene = c("ProteinCodingGene", "ncRNA","miRNA", "rRNA", "tRNA", "snRNA"),
    ncRNA = c("miRNA", "rRNA", "tRNA", "snRNA")
  )

  ## Hierarchy check
  # First check the class of the given object
  current_class <- class(object)[1]
  if (!(target_class %in% valid_classes)) {
    stop("The target class is not valid. Must be one of the following: ",
         paste(valid_classes, collapse = ", "))
  }
  if (current_class != "Gene" && !(target_class %in%
                                   class_hierarchy[[current_class]])) {
    stop("The current class", current_class, " can't be transformed in ",
         target_class)
  }

  # Copy the object and add the additional fields
  new_object <- as(object, target_class)

  for (field_name in names(additional_fields)) {
    slot(new_object, field_name) <- additional_fields[[field_name]]
  }

  # Update the class
  validObject(new_object)

  return(new_object)
}

#' Create a Protein Coding Gene object
#'
#' This function creates an instance of the ProteinCodingGene class.
#'
#' @inheritParams createGene
#' @param transcript_id ID of the alternative transcript considered
#' (Ex: Ensembl Transcript ID or RefSeq ID)
#' @param exons GRanges object representing the exons of the gene.
#' @param protein_id Protein ID associated with the gene (preferably from
#' UniProtKB).
#' @param protein_name Protein name (preferably from UniprotKB).
#' @param protein_sequence The aminoacids sequence.
#' @return A ProteinCodingGene object.
#' @examples
#' # gene <- createProteinCodingGene(
#' # id = "ENSG00000139618.19",
#' # symbol = "BRCA2",
#' # name = "BRCA2 DNA repair associated",
#' # description = "Plays a role in DNA repair...",
#' # assembly = "hg38",
#' # chromosome = "chr13",
#' # start = 32315508,
#' # end = 32400268,
#' # strand = "+",
#' # sequence = "ACGTACCGTCAAA...",
#' # transcirpt_id = "ENST00000544455.6",
#' # exons = exons,                    # A Granges() object defined previously
#' # protein_id = "P51587",
#' # protein_name = "Breast cancer type 2 susceptibility protein",
#' # protein_sequence = "mpigskerptffeifktrcn....")
#' @export
createProteinCodingGene <- function(id, symbol = NA_character_,
                                    name = NA_character_,
                                    description = NA_character_, assembly,
                                    chromosome, start, end, strand,
                                    sequence = Biostrings::DNAString(""),
                                    transcript_id = NA_character_,
                                    exons, protein_id,
                                    protein_name = NA_character_,
                                    protein_sequence = Biostrings::AAString(""))
  {
  sequence <- validateGeneInputs(start, end, strand, sequence)
  if (is.character(protein_sequence)) {
    tryCatch({
      protein_sequence <- Biostrings::AAString(protein_sequence)
    }, error = function(e) {
      stop("Invalid input: the provided string is not a valid aminoacid
           sequence. The 1 letter notation should be used.")
    })
  }
  new("ProteinCodingGene",
      id = id,
      symbol = symbol,
      name = name,
      description = description,
      assembly = assembly,
      chromosome = chromosome,
      start = start,
      end = end,
      strand = strand,
      sequence = sequence,
      transcript_id=transcript_id,
      exons = exons,
      protein_id = protein_id,
      protein_name = protein_name,
      protein_sequence = protein_sequence)
}


#' Create a Non-Coding RNA Gene object
#'
#' This function creates an instance of the NonCodingRNA class.
#'
#' @inheritParams createGene
#' @param transcript_id ID of the alternative transcript considered (Ex: Ensembl
#' Transcript ID or RefSeq ID)
#' @param type Type of the non-coding RNA ('lncRNA', 'siRNA','rRNA', 'tRNA'...)
#' @param RNA_sequence RNA sequence.
#' @return A Non Coding RNA (ncRNA) object.
#' @examples
#' # ncRNA_miRNA <- create_ncRNA(
#' # id = "ENSG00000208008",
#' # symbol = "MIR125A",
#' # name = "microRNA 125a",
#' # description = "MIR125A (MicroRNA 125a) is an RNA Gene, and is affiliated
#' # with the miRNA class. Diseases associated with MIR125A include
#' # Medulloblastoma and Lung Cancer.",
#' # assembly = "hg38",
#' # chromosome = "chr19",
#' # start = 51693254,
#' # end = 51693339,
#' # strand = "+",
#' # sequence = "ACGTACCGTCAAA...",
#' # transcript_id = "ENST00000385273.1",
#' # type = 'miRNA',
#' # RNA_sequence = "UGCCAGUCUC....")
#'
#' @export
create_ncRNA <- function(id, symbol = NA_character_, name = NA_character_,
                         description = NA_character_, assembly, chromosome,
                         start, end, strand,
                         sequence = Biostrings::DNAString(""),
                         transcript_id = NA_character_,
                         type, RNA_sequence = Biostrings::RNAString("")) {
  sequence <- validateGeneInputs(start, end, strand, sequence)
  if (is.character(RNA_sequence)) {
    tryCatch({
      RNA_sequence <- Biostrings::RNAString(RNA_sequence)
    }, error = function(e) {
      stop("Invalid input: the provided string is not a valid RNA sequence.")
    })
  }
  new("ncRNA",
      id = id,
      symbol = symbol,
      name = name,
      description = description,
      assembly = assembly,
      chromosome = chromosome,
      start = start,
      end = end,
      strand = strand,
      sequence = sequence,
      transcript_id = transcript_id,
      type = type,
      RNA_sequence = RNA_sequence)
}



#' Create a miRNA Gene object
#'
#' This function creates an instance of the MicroRNA class.
#'
#' @inheritParams create_ncRNA
#' @param target_genes A list containing the IDs of the target genes of the
#' microRNA.
#' @param seed Seed sequence (usually around 8 nucleotides).
#' @param mature_sequence Mature sequence of the miRNA (usally around 22
#' nucleotides).
#' @return A miRNA object.
#' @examples
#' # miRNA125a <- create_miRNA(
#' # id = "ENSG00000208008",
#' # symbol = "MIR125A",
#' # name = "microRNA 125a",
#' # description = "MIR125A (MicroRNA 125a) is an RNA Gene, and is affiliated
#' # with the miRNA class. Diseases associated with MIR125A include
#' # Medulloblastoma and Lung Cancer.",
#' # assembly = "hg38",
#' # chromosome = "chr19",
#' # start = 51693254,
#' # end = 51693339,
#' # strand = "+",
#' # sequence = "ACGTACCGTCAAA",
#' # transcript_id = "ENST00000385273.1",
#' # type = 'miRNA',
#' # RNA_sequence = "UGCCAGUCUC",
#' # target_genes = c("ENSG00000141736", "ENSG00000065361"),
#' # seed = 'GAGUCCC',
#' # mature_sequence = "ucccugagacccuuuaaccuguga")
#'
#' @export
create_miRNA <- function(id, symbol = NA_character_, name = NA_character_,
                         description = NA_character_, assembly, chromosome,
                         start, end, strand,
                         sequence = Biostrings::DNAString(""),
                         transcript_id = NA_character_,
                         type = "miRNA",
                         RNA_sequence = Biostrings::RNAString(""), target_genes,
                         seed = NA_character_,
                         mature_sequence = Biostrings::RNAString("")) {
  sequence <- validateGeneInputs(start, end, strand, sequence)

  if (is.character(RNA_sequence)) {
    tryCatch({
      RNA_sequence <- Biostrings::RNAString(RNA_sequence)
    }, error = function(e) {
      stop("Invalid input: the provided string is not a valid RNA sequence.")
    })
  }

  if (is.character(mature_sequence)) {
    tryCatch({
      mature_sequence <- Biostrings::RNAString(mature_sequence)
    }, error = function(e) {
      stop("Invalid input: the provided string is not a valid RNA sequence.")
    })
  }

  new("miRNA",
      id = id,
      symbol = symbol,
      name = name,
      description = description,
      assembly = assembly,
      chromosome = chromosome,
      start = start,
      end = end,
      strand = strand,
      sequence = sequence,
      transcript_id = transcript_id,
      type = type,
      RNA_sequence = RNA_sequence,
      target_genes = target_genes,
      seed = seed,
      mature_sequence = mature_sequence)
}

#' Function to create a tRNA object
#'
#' Creates an object of class tRNA with the specified parameters
#'
#' @inheritParams create_ncRNA
#' @param aminoacid Name of the aminoacid transported by the tRNA (can be either
#' full name, or 1-3 letters code notation).
#' @param anticodon Anticodon contained in the tRNA (rom 5' to 3').
#' @return An object of class tRNA
#' @examples
#' # ncRNA_tRNA_Phe <- create_tRNA(
#' # id = "ENSG00000210049.1",
#' # name = "Mitochondrial tRNA phenylalanine",
#' # assembly = "hg38",
#' # chromosome = "chrM",
#' # start = 577,
#' # end = 647,
#' # strand = "+",
#' # transcript_id = "ENST00000387314.1",
#' # type = "tRNA",
#' # aminoacid = "phenylalanine",
#' # anticodon = "GAA")
#' @export
create_tRNA <- function(id, symbol = NA_character_, name = NA_character_,
                        description = NA_character_, assembly, chromosome,
                        start, end, strand,
                        sequence = Biostrings::DNAString(""),
                        transcript_id = NA_character_, type = "tRNA",
                        RNA_sequence = Biostrings::RNAString(""), aminoacid,
                        anticodon = NA_character_) {
  sequence <- validateGeneInputs(start, end, strand, sequence)

  if (is.character(RNA_sequence)) {
    tryCatch({
      RNA_sequence <- Biostrings::RNAString(RNA_sequence)
    }, error = function(e) {
      stop("Invalid input: the provided string is not a valid RNA sequence.")
    })
  }

  new("tRNA",
      id = id,
      symbol = symbol,
      name = name,
      description = description,
      assembly = assembly,
      chromosome = chromosome,
      start = start,
      end = end,
      strand = strand,
      sequence = sequence,
      transcript_id = transcript_id,
      type = type,
      RNA_sequence = RNA_sequence,
      aminoacid = aminoacid,
      anticodon = anticodon)
}


#' Function to create a rRNA object
#'
#' Creates an object of class rRNA with the specified parameters
#'
#' @inheritParams create_ncRNA
#' @param rRNA_type Type of rRNA subunit. The ones accepted are "18S", "28S",
#' "5S", "5.8S".
#' @param ribosomal_subunit In which subunit of the ribosome the rRNA is found.
#' @return An object of class rRNA
#' @examples
#' # rRNA_18s <- create_rRNA(
#' # id = "ENSG00000225840",
#' # symbol = "RNA18SN5",
#' # name = "RNA, 18S ribosomal pseudogene",
#' # description = "45S ribosomal DNA (rDNA) arrays, or clusters,
#' # are present on human chromosomes 13, 14, 15...",
#' # assembly = "hg19",
#' # chromosome = "chrY",
#' # start = 10197256,
#' # end = 10199103,
#' # strand = "-",
#' # transcript_id = "ENST00000445125.2",
#' # type = "tRNA",
#' # rRNA_type = "18S",
#' # ribosomal_subunit = "40S")
#'
#' @export
create_rRNA <- function(id, symbol = NA_character_, name = NA_character_,
                        description = NA_character_, assembly, chromosome,
                        start, end, strand,
                        sequence = Biostrings::DNAString(""),
                        transcript_id = NA_character_,
                        type = "rRNA", RNA_sequence = Biostrings::RNAString(""),
                        rRNA_type, ribosomal_subunit = NA_character_) {
  sequence <- validateGeneInputs(start, end, strand, sequence)

  if (is.character(RNA_sequence)) {
    tryCatch({
      RNA_sequence <- Biostrings::RNAString(RNA_sequence)
    }, error = function(e) {
      stop("Invalid input: the provided string is not a valid RNA sequence.")
    })
  }

  new("rRNA",
      id = id,
      symbol = symbol,
      name = name,
      description = description,
      assembly = assembly,
      chromosome = chromosome,
      start = start,
      end = end,
      strand = strand,
      sequence = sequence,
      transcript_id = transcript_id,
      type = type,
      RNA_sequence = RNA_sequence,
      rRNA_type = rRNA_type,
      ribosomal_subunit = ribosomal_subunit)
}



#' Function to create a snRNA object
#'
#' Creates an object of class snRNA with the specified parameters
#'
#' @inheritParams create_ncRNA
#' @param snRNA_class A short string that defines a specific type of snRNA,
#' usually differentiated by function. (Ex: U1, U2...)
#' @param subclass Subclass of the snRNA, defined by its localisation: either
#' 'snoRNA' or 'scaRNA'.
#' @param associated_proteins A list of the ID of the proteins that the snRNA
#' interacts with.
#' @return An object of class snRNA
#' @examples
#' # ncRNA_snRNA <- create_snRNA(
#' # id = "ENSG00000206652",
#' # symbol = "RNU1-1",
#' # name = "RNA, U1 small nuclear 1",
#' # description = "",
#' # assembly = "hg38",
#' # chromosome = "chr1",
#' # start = 16514122,
#' # end = 16514285,
#' # strand = "-",
#' # sequence = "tttcatacttacctggcagg",
#' # RNA_sequence = 'aaaguaug',
#' # transcript_id = "ENST00000383925.1",
#' # snRNA_class = 'U1')
#'
#' @export
create_snRNA <- function(id, symbol = NA_character_, name = NA_character_,
                         description = NA_character_, assembly, chromosome,
                         start, end, strand,
                         sequence = Biostrings::DNAString(""),
                         transcript_id = NA_character_, type = "snRNA",
                         RNA_sequence = Biostrings::RNAString(""), snRNA_class,
                         subclass = NA_character_,
                         associated_proteins = character(0)) {
  sequence <- validateGeneInputs(start, end, strand, sequence)

  if (is.character(RNA_sequence)) {
    tryCatch({
      RNA_sequence <- Biostrings::RNAString(RNA_sequence)
    }, error = function(e) {
      stop("Invalid input: the provided string is not a valid RNA sequence.")
    })
  }

  new("snRNA",
      id = id,
      symbol = symbol,
      name = name,
      description = description,
      assembly = assembly,
      chromosome = chromosome,
      start = start,
      end = end,
      strand = strand,
      sequence = sequence,
      transcript_id = transcript_id,
      type = type,
      RNA_sequence = RNA_sequence,
      snRNA_class = snRNA_class,
      subclass = subclass,
      associated_proteins = associated_proteins)
}

#' Get a slot value
#'
#' Get a slot value without directly having to access it.
#'
#' @param object Object of Gene, and its derivates, class.
#' @param slot_name Name of the slot.
#' @return Slot value
#' @examples
#' # Access the mature sequence of a miRNA
#' # get_slot(my_object, 'mature_sequence')
#' @export
get_slot <- function(object, slot_name) {
  if (!is(object, "Gene")) stop("The object needs to inherit class 'Gene'.")
  if (!slot_name %in% slotNames(object)) stop("This slot name is not present in
                                              this class.")
  slot(object, slot_name)
}


#' Set a specific slot value
#'
#' Set a slot value without directly having to access it; the value can be
#' either added or edited.
#'
#' @param object Object of Gene, and its derivates, class.
#' @param slot_name Name of the slot.
#' @param value The new value to assign to the slot.
#' @return The updated object.
#' @examples
#' # Changing phenylalanine notation for the object tRNA created before.
#' # tRNA <- set_slot(tRNA, 'aminoacid', 'Phe')
#'
#' @export
set_slot <- function(object, slot_name, value) {
  if (!is(object, "Gene")) stop("The object needs to inherit class 'Gene'.")
  if (!slot_name %in% slotNames(object)) stop("This slot name is not present in
                                              this class.")

  slot(object, slot_name) <- value
  validObject(object)
  return(object)
}
