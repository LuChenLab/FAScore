#' AS mapped to isoforms
#'
#' @param object A FAScore object. see the constructor function '\code{\link{FAScoreDataSet}}'.
#' @param gtf An annotation file was consistent with identifying AS.
#' @param AStype The types of AS event sites, including 'exonic' (exonic sites in genomes) and 'intronic' (intronic sites in genomes, eg. splicing sits).
#' @param cores The number of threads. (default:10)
#'
#'
#' @importFrom dplyr `%>%`
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom rtracklayer import
#' @importFrom parallel mclapply
#'
#' @return Return a GRange object. The alternative splicing sites were aligned to transcripts by an annotation file.
#' @export
#'

ASmapIso <- function(object, gtf, AStype = "exonic", cores = 10) {
  if (class(object) != "FAScore") {
    stop("Object must be a FAScore object")
  }

  mygr <- rownames(object) %>% GRanges()
  gtf_exon <- rtracklayer::import(gtf) %>% subset(., type == "exon")
  col <- grep("id|name", colnames(mcols(gtf_exon)), value = T) %>% grep("gene|transcript",., value = T)
  gtf_exon <- gtf_exon[,col]

  if(AStype == "exonic"){
    ov <- findOverlaps(mygr, gtf_exon)

    tmp <- mygr[queryHits(ov),]
    mcols(tmp) <- cbind(mcols(tmp), mcols(gtf_exon[subjectHits(ov),]))
    tmp$ASID <- paste0(seqnames(tmp),":", start(tmp),"-",end(tmp),":", strand(tmp) )

    tmpGR <- split(tmp, tmp$ASID, drop = TRUE)
    MyASmapIso <- mclapply(tmpGR, function(x){
      x$TranscriptID <- paste0(unique(x$transcript_id), collapse = ",")
      x$TranscriptName <- paste0(unique(x$transcript_name), collapse = ",")
      x$GeneID <- paste0(unique(x$gene_id), collapse = ",")
      x$GeneName <- paste0(unique(x$gene_name), collapse = ",")
      return(x)
    }, mc.cores = 10) %>% GRangesList(.) %>%
      unlist() %>%
      .[,c("ASID","TranscriptID","TranscriptName","GeneID","GeneName")] %>% unique

  } else if(AStype == "intronic"){
    gtf_txdb <- makeTxDbFromGFF(gtf)
    tmp <- intronsByTranscript(gtf_txdb, use.names = T) %>% unlist

    tmp$ASID <- paste0(seqnames(tmp), ":", start(tmp), "-", end(tmp), ":", strand(tmp))
    mt <- gtf_exon[match( names(tmp), gtf_exon$transcript_id), ]
    colnames(mcols(mt)) <- c("GeneID", "GeneName", "TranscriptID","TranscriptName")
    mcols(tmp) <- cbind(mcols(tmp), mcols(mt))
    gtf_intron_list <- split(tmp, tmp$ASID, drop = TRUE)


    MyASmapIso <- mclapply(gtf_intron_list, function(x){
      x$TranscriptID <- paste0(unique(x$transcript_id), collapse = ",")
      x$TranscriptName <- paste0(unique(x$transcript_name), collapse = ",")
      x$GeneID <- paste0(unique(x$gene_id), collapse = ",")
      x$GeneName <- paste0(unique(x$gene_name), collapse = ",")
      return(x)
    }, mc.cores = cores) %>% GRangesList(.) %>%
      unlist() %>%
      .[,c("ASID","TranscriptID","TranscriptName","GeneID","GeneName")] %>% unique

    names(MyASmapIso) <- NULL

  } else{
    stop("the AStype must be 'exonic' or 'intronic' ")
  }

  res <- mcols(MyASmapIso)[ match( rownames(object), MyASmapIso$ASID), ]
  mcols(rowRanges(object)) <- res


  return(object)
}




utils::globalVariables(c(".", "type", "ASDyScore", "GeDyScore", "AS", "ASID", "TranscriptID","density","x","y1","y2","y3",
                         "GeneID", "GeneName", "matchTransID", "matchGeneID", "gene_id", "Run","TranscriptName",
                         "Exp", "FAScore", "Color", "value", "CellType", "Type"))



