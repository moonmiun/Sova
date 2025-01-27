#' @import Biostrings
#' @import biomaRt
#' @import Gviz
#' @import IRanges
#' @importFrom S4Vectors Rle
#' @importFrom GenomicRanges GRanges
NULL


#' Convert FASTA file.
#'
#' Reads a FASTA file and returns a Biostring::DNAString.
#'
#' @param filepath The FASTA file path.
#' @return A Biostrings::DNAString that contains your DNA sequence.
#' @examples
#' sequence <- readFasta('C:/User/Desktop/file.fa')
#' @export
readFasta <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("The file doesn't exist.")
  }

  # Reads the file line by line
  lines <- readLines(filepath)

  # Removes the FASTA header
  header <- lines[1]
  if (!startsWith(header, ">")) {
    stop("The file doesn't have a valid FASTA header.")
  }

  # Combines all of the lines in a single sequence.
  sequence <- paste(lines[-1], collapse = "")

  # Removes empty spaces
  sequence <- gsub("\\s+", "", sequence)

  # Converts it in a Biostrings::DNAstring element
  sequence <- Biostrings::DNAString(sequence)

  return(sequence)
}

#' Validate sequence and chromosome coordinates
#'
#' This function validates the input sequence and chromosome coordinates.
#' The sequence is returned as a DNAString (from Biostrings) if it complies
#' with the conditions.
#'
#' @param start The position at which the gene starts (1 based).
#' @param end The position at which the gene ends.
#' @param strand The DNA strand where the gene is located.
#' @param sequence The DNA sequence (given as string, DNAString object, FASTA
#' file or path containing a FASTA file.)
#' @return A DNAString object
#' @export
validateGeneInputs <- function(start, end, strand,
                               sequence = Biostrings::DNAString("")) {
  if (inherits(sequence, "DNAString")) {
    return(sequence)
  }
  if (is.na(sequence)) {
    # Case 0: No sequence was provided, return NA or handle accordingly.
    return(Biostrings::DNAString(""))
  }
  if (is.character(sequence)) {
    # Case 1 : a string element was given.
    # Verifying if it's a valid sequence of DNA...
    tryCatch({
      sequence <- Biostrings::DNAString(sequence)
    }, error = function(e) {
      # If it's not a valid sequence, checks if it's a path to a FASTA file
      if (file.exists(sequence) && grepl("\\.fa$", sequence,
                                         ignore.case = TRUE)) {
        sequence <- readFasta(sequence)
      } else {
        stop("Invalid input: the string is neither a valid DNA sequence nor a
             valid FASTA file path.")
      }
    })
  } else if (inherits(sequence, "connection")) {
    # Case 2: the fasta was opened in the RStudio environment.
    # Reads the file line by line
    lines <- readLines(sequence)

    # Removes the FASTA header
    header <- lines[1]
    if (!startsWith(header, ">")) {
      stop("The file doesn't have a valid FASTA header.")
    }

    # Combines all of the lines in a single sequence.
    sequence <- paste(lines[-1], collapse = "")

    # Removes empty spaces
    sequence <- gsub("\\s+", "", sequence)

    # Converts it in a Biostrings::DNAstring element
    sequence <- Biostrings::DNAString(sequence)

  } else if (!inherits(sequence, "DNAString")) {
    # Case 4: the content is not in any valid form.
    stop("Invalid input: must be a DNA sequence string, a DNAString object, a
         valid FASTA file path, or loaded FASTA content.")
  }

  # Validation of other parameters.
  if (start < 0 || end < 0) stop("Start and end positions must be
                                 non-negative.")
  if (start > end) stop("Start position must be less than or equal to end
                        position.")
  if (!strand %in% c("+", "-")) stop("Strand must be either '+' or '-'.")
  return(sequence)
}



#' Import a Gene from ENSEMBL
#'
#' This function imports gene information from ENSEMBL (assembly: GRCh38) and
#' returns an object of class Gene.
#'
#' @param ensembl_id Gene ID, specifically the ENSEMBL one.
#' @param mart An object of type Mart from package biomaRT. It establishes a
#' connection with the biomaRT database to interact with ENSEMBL data.
#' @return An object of class Gene.
#' @export
importGeneFromEnsembl <- function(ensembl_id, mart = NULL) {
  # The default for mart is NULL: if it's not given by the user, the function
  # creates it.
  if (is.null(mart)) {
    mart <- biomaRt::useEnsembl(biomart = "genes",
                                dataset = "hsapiens_gene_ensembl")
  }

  # Obtains information for the gene.
  gene_data <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "description",
                   "chromosome_name",
                   "start_position", "end_position", "strand", "version"),
    filters = "ensembl_gene_id",
    values = ensembl_id,
    mart = mart
  )

  if (nrow(gene_data) == 0) {
    stop("Gene not found in Ensembl.")
  }

  # Converts the strand information (saved as either 1 or 0) in '+' or '-'
  strand <- ifelse(gene_data$strand == 1, "+", "-")

  # Creates a Gene object
  createGene(
    id = gene_data$ensembl_gene_id,
    symbol = gene_data$hgnc_symbol,
    description = gene_data$description,
    assembly = "GRCh38",
    chromosome = paste0("chr", gene_data$chromosome_name),
    start = gene_data$start_position,
    end = gene_data$end_position,
    strand = strand
  )
}



#' Visualize a gene object with Gviz
#'
#' This function lets the user visualize an object of class Gene in a more
#' user-friendly manner.
#'
#' @param gene An object of class Gene.
#' @param ideogram If it's different from 'NA' downloads from ensembl the
#' ideogram of the current chromosome. It requires internet connection and
#' it might require some time.
#' @return A Gviz visualization of the gene object.
#' @export
visualizeGene <- function(gene, ideogram = NA) {
  if (!inherits(gene, "Gene")) {
    stop("The object provided is not a valid Gene object.")
  }

  # Check for subclass (ProteinCodingGene or ncRNA) and uses transcript_id if
  # available
  if (inherits(gene, "ProteinCodingGene") || inherits(gene, "ncRNA")) {
    # Use transcript_id if it's available, otherwise fallback to gene@id
    transcript <- ifelse(!is.null(gene@transcript_id) &&
                           gene@transcript_id != "", gene@transcript_id,
                         gene@id)
  } else {
    # Default case, just use gene@id
    transcript <- gene@id
  }

  gene_range = GRanges(
    seqnames = Rle(c(gene@chromosome)),
    ranges = IRanges(start = c(gene@start),
                     end = c(gene@end)),
    strand = Rle(c(gene@strand)))
  # Uses the function already implemented in Gviz
  geneTrack <- GeneRegionTrack(
    range = gene_range,
    genome = gene@assembly,
    chromosome = gene@chromosome,
    strand = gene@strand,
    transcript = transcript,
    name = "Gene",
    col = "grey",
    fill = "grey",
    lwd = 1
  )

  if (inherits(gene, "ProteinCodingGene")) {
    exonTrack <- AnnotationTrack(
      range = gene@exons,
      name = "Exons",
      genome = gene@assembly,
      col = "black",
      fill = "grey",
      lwd = 1
    )
  }



  # Ensure it doesn't go below 0
  start_position <- max(gene@start - 1000, 0)
  end_position <- gene@end + 1000

  gtrack <- GenomeAxisTrack()

  if (!is.na(ideogram)) {
    itrack <- IdeogramTrack(genome = gene@assembly,
                            chromosome = gene@chromosome)
    if (inherits(gene, "ProteinCodingGene")){
      return(plotTracks(list(itrack, gtrack, geneTrack, exonTrack),
                        from = start_position, to = end_position,
                        main = paste("Gene:", transcript))
      )
    }
    else {
      return(plotTracks(list(itrack, gtrack, geneTrack), from = start_position,
                        to = end_position,
                        main = paste("Gene:", transcript)))
    }
  }

  # Visualize the Gene
  if (inherits(gene, "ProteinCodingGene")){
    return(plotTracks(list(geneTrack, gtrack, exonTrack), from = start_position,
                      to = end_position,
                      main = paste("Gene:", transcript))
    )
  }
  else {
    return(plotTracks(list(gtrack, geneTrack), from = start_position,
                      to = end_position,
                      main = paste("Gene:", transcript)))
  }
}




#' Length of the genetic product function
#'
#' Calculates the length of the genetic product of each class. What is the
#' genetic product differs from class to class.
#'
#' @param object Object of class Gene and its derivates.
#' @return The length of the genetic product.
#' @export
setGeneric("lengthProduct", function(object) standardGeneric("lengthProduct"))

#' @rdname lengthProduct
#' @examples
#' # Class Object
#' This is the case in which there is no information about the gene so the
#' genetic product is the length of its RNA transcription.
#' lengthProduct(Gene)
setMethod("lengthProduct", "Gene", function(object) {
  # 1 is added because the coordinates are inclusive, so that means that the
  # last position is still part of the gene.
  object@end - object@start + 1
})

#' @rdname lengthProduct
#' @examples
#' # Class ProteinCodingGene
#' # Example in which the RNA sequence is given.
#' set_slot(protein_gene, 'protein_sequence', 'MTEYKLVVVG')
#' lengthProduct(protein_gene)
#' > 10
#'
#' # Example in which the RNA sequence is not given so the function uses the
#' length of the exons.
#' lengthProduct(protein_gene)
setMethod("lengthProduct", "ProteinCodingGene", function(object) {
  if (length(object@protein_sequence) > 0
      && nzchar(as.character(object@protein_sequence))) {
    # If the protein sequence is available, then it simply calculates the
    # length of the string.
    len <- nchar(object@protein_sequence)
  } else {
    len <- sum(width(object@exons)) %/% 3
  }
  return(len)
})


#' @rdname lengthProduct
#' @examples
#' # Class ncRNA
#' In this case the genetic product is simply the RNA transcription of the gene,
#' which has the same length as the gene itself.
#' This method is valid also for objects of sub-classes
#' 'tRNA', 'rRNA' and 'snRNA'.
#' lengthProduct(ncRNA)
setMethod("lengthProduct", "ncRNA", function(object) {
  object@end - object@start + 1
})

#' @rdname lengthProduct
#' @examples
#' # Class miRNA
#' set_slot(miRNA_object, 'mature_sequence', 'UGAGGUAGUAGGUUGUAUAGUU')
#' lengthProduct(miRNA_object)
#' > 22
setMethod("lengthProduct", "miRNA", function(object) {
  # If there's the sequence of the mature miRNA, returns that length.
  if (length(object@mature_sequence) > 0) {
    return(nchar(as.character(object@mature_sequence)))
  } else {
    # Altrimenti restituisce la lunghezza del gene
    return(object@end - object@start + 1)
  }
})
