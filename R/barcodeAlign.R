#' @export
#' @importFrom S4Vectors metadata metadata<- DataFrame
#' @importFrom BiocParallel SerialParam bpstart bpstop bpisup bpmapply
barcodeAlign <- function(sequences, barcodes, gapOpening=5, gapExtension=1, BPPARAM=SerialParam())
# Pulls out the barcodes, aligns them against all possible options,
# and reports the results.
#
# written by Aaron Lun
# created 19 September 2018
{
    by.core <- .parallelize(sequences, BPPARAM)
    current.score <- next.best <- rep(-Inf, length(sequences))
    current.id <- rep(NA_integer_, length(sequences))

    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    outputs <- bpmapply(FUN=.align_BA_internal, barcodes=by.core, 
        MoreArgs=list(reference=barcodes, gap.opening=gapOpening, gap.extension=gapExtension),
        BPPARAM=BPPARAM, SIMPLIFY=FALSE, USE.NAMES=FALSE)
    
    output <- Reduce(rbind, outputs)
    metadata(output) <- list(gapOpening=gapOpening, gapExtension=gapExtension, barcodes=barcodes)
    output
}

#' @importFrom Biostrings quality
.align_BA_internal <- function(barcodes, references, gap.opening, gap.extension) {
    quals <- quality(barcodes)
    
    current.score <- next.best <- rep(-Inf, length(barcodes))
    current.id <- rep(NA_integer_, length(barcodes))
    
    for (r in seq_along(references)) {
        reference <- toupper(as.character(references[r]))
        scores <- .Call(cxx_barcode_align, barcodes, quals, .create_encoding_vector(quals), 
                    gap.opening, gap.extension, reference)
        
        scores <- unlist(scores)
        
        # Need to update both the best and the next best on record.
        keep <- scores > current.score
        second.keep <- !keep & scores > next.best
        
        current.id[keep] <- r
        next.best[keep] <- current.score[keep]
        current.score[keep] <- scores[keep]
        next.best[second.keep] <- scores[second.keep]
    }
    
    DataFrame(barcode=current.id, score=current.score, gap=current.score - next.best)
}
