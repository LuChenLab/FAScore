#' AS mapped to isoforms
#'
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param gtf An annotation file was consistent with identifying AS.
#' @param cores The number of threads. (default is 10)
#'
#'
#' @importFrom dplyr `%>%`
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom IRanges relist
#' @importFrom rtracklayer import
#' @importFrom stringr str_split_fixed
#' @importFrom parallel mclapply
#'
#' @return Return a GRange object. The alternative splicing sites were aligned to transcripts by an annotation file.
#'
#'

ASmapIso <- function(object, gtf, cores = 10) {
  gtf_exon <- rtracklayer::import(gtf) %>% subset(., type == "exon")

  gtf_txdb <- makeTxDbFromGFF(gtf)
  gtf_intron_bytrans <- intronsByTranscript(gtf_txdb, use.names = T)

  ### add SJ ID
  tmp <- unlist(gtf_intron_bytrans)
  tmp$ASID <- paste0(seqnames(tmp), ":", start(tmp), "-", end(tmp), ":", strand(tmp))
  tmp$TranscriptID <- names(tmp)
  mt <- gtf_exon[match(tmp$TranscriptID, gtf_exon$transcript_id), ]
  tmp$GeneID <- mt$gene_id
  tmp$GeneName <- mt$gene_name
  gtf_intron_bytrans <- relist(tmp, gtf_intron_bytrans)

  ## mapped
  gtf_intron_gr <- unlist(gtf_intron_bytrans)
  gtf_intron_gr <- subset(gtf_intron_gr, ASID %in% (stringr::str_split_fixed(rownames(assay(object)), "[|]", 2)[, 2]))
  gtf_intron_gr <- split(gtf_intron_gr, gtf_intron_gr$ASID, drop = TRUE)

  MySJmapIso <- mclapply(gtf_intron_gr, function(gr) {
    Iso <- paste0(gr$TranscriptID, collapse = ",")
    gen <- paste0(unique(gr$GeneID), collapse = ",")
    genN <- paste0(unique(gr$GeneName), collapse = ",")

    gr <- unique(gr)
    gr$TranscriptID <- Iso
    gr$GeneID <- gen
    gr$GeneName <- genN

    return(gr)
  }, mc.cores = cores) %>%
    GRangesList(.) %>%
    unlist()

  names(MySJmapIso) <- NULL

  return(MySJmapIso)
}





#' AS mapped to host genes
#'
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param gtf An annotation file was consistent with identifying AS.
#' @param cores The number of threads. (default is 10)
#'
#'
#' @importFrom dplyr `%>%`
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges GRangesList
#' @importFrom rtracklayer import
#' @importFrom stringr str_split_fixed
#' @importFrom parallel mclapply
#'
#'
#' @return Return a GRange object. The alternative splicing sites were aligned to genes by an annotation file.
#' @export

HostGene <- function(object, gtf, cores = 10) {

  gtfFile <- rtracklayer::import(gtf)
  gtfFile <- subset(gtfFile, type == "gene")
  gtfFile <- gtfFile[, c(5, 7)]

  IsoGR <- ASmapIso(object, gtf = gtf, cores = cores)

  candiASID <- base::setdiff(
    str_split_fixed(rownames(assay(object)), "[|]", 2)[, 2], IsoGR$ASID
  ) %>% GRanges()


  ov <- findOverlaps(candiASID, gtfFile, ignore.strand = F, type = c("within"))

  AS <- candiASID[queryHits(ov), ]
  AS$GeneID <- gtfFile[subjectHits(ov), ]$gene_id
  AS$GeneName <- gtfFile[subjectHits(ov), ]$gene_name
  AS$ASID <- paste0(seqnames(AS), ":", start(AS), "-", end(AS), ":", strand(AS))

  AS <- split(AS, AS$ASID, drop = TRUE)

  AS2 <- mclapply(AS, function(gr) {
    gen <- paste0(unique(gr$GeneID), collapse = ",")
    genN <- paste0(unique(gr$GeneName), collapse = ",")

    gr <- unique(gr)
    gr$TranscriptID <- NA
    gr$GeneID <- gen
    gr$GeneName <- genN

    return(gr)
  }, mc.cores = cores) %>%
    GRangesList(.) %>%
    unlist()

  names(AS2) <- NULL
  # object@rowRanges <- c(IsoGR, AS2)
  allrange <- c(IsoGR, AS2)
  rowRanges(object)$ASID <- do.call(rbind, strsplit(names(rowRanges(object)), split = "[|]"))[,2]

  mcols(rowRanges(object)) <- cbind(mcols(rowRanges(object)),mcols(allrange)[match(rowRanges(object)$ASID,allrange$ASID),2:4])

  return(object)
}











