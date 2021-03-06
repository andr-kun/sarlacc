#' @export
#' @importFrom Biostrings QualityScaledDNAStringSet
#' @importFrom methods is
#' @importFrom BiocParallel bpmapply SerialParam bpstart bpstop bpisup
#' @importFrom ShortRead FastqSampler
tuneAlignment <- function(adaptor1, adaptor2, filepath, tolerance=200, number=10000,
    gapOp.range=c(4, 10), gapExt.range=c(1, 5), qual.type=c("phred", "solexa", "illumina"), 
    BPPARAM=SerialParam()) 
# This function scrambles the input subject sequences and tries to identify the
# alignment parameters that minimize the overlap between the true alignment 
# scores (to the unscrambled subject) and those to the scrambled sequences. 
#
# written by Aaron Lun
# created 26 February 2018
{
    adaptor1 <- toupper(as.character(adaptor1))
    adaptor2 <- toupper(as.character(adaptor2))
    qual.type <- match.arg(qual.type)
    qual.class <- .qual2class(qual.type)

    f <- FastqSampler(filepath, number)
    on.exit(close(f))
    reads <- .FASTQ2QSDS(yield(f),  qual.class)

    if (length(reads)==0L) {
        return(list(parameters=list(gapOpening=NA_integer_, gapExtension=NA_integer_),
                    scores=list(reads=numeric(0), scrambled=numeric(0))))
    }

    # Taking subsequences of the reads for pairwise alignment.
    reads.out <- .get_front_and_back(reads, tolerance)
    reads.start <- reads.out$front
    reads.end <- reads.out$back
    scrambled.start <- .scramble_input(reads.start, TRUE)
    scrambled.end <- .scramble_input(reads.end, TRUE)

    reads.start <- .parallelize(reads.start, BPPARAM)
    reads.end <- .parallelize(reads.end, BPPARAM)
    scrambled.start <- .parallelize(scrambled.start, BPPARAM)
    scrambled.end <- .parallelize(scrambled.end, BPPARAM)

    # Performing a grid search to maximize the separation in scores.
    gapOp.range <- as.integer(cummax(gapOp.range))
    gapExt.range <- as.integer(cummax(gapExt.range))
    max.score <- 0L
    final.gapOp <- final.gapExt <- NA
    final.reads <- final.scrambled <- NULL

    if (!bpisup(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    for (go in seq(gapOp.range[1], gapOp.range[2], by=1)) { 
        for (ge in seq(gapExt.range[1], gapExt.range[2], by=1)) { 
            out <- bpmapply(FUN=.align_TA_internal, reads.start=reads.start, reads.end=reads.end,
                scrambled.start=scrambled.start, scrambled.end=scrambled.end,
                MoreArgs=list(adaptor1=adaptor1, adaptor2=adaptor2, gap.opening=go, gap.extension=ge), 
                SIMPLIFY=FALSE, BPPARAM=BPPARAM, USE.NAMES=FALSE)
            read.scores <- unlist(lapply(out, "[[", i="reads"))
            scrambled.scores <- unlist(lapply(out, "[[", i="scrambled"))

            cur.score <- .tied_overlap(read.scores, scrambled.scores)
            if (max.score < cur.score) {
                max.score <- cur.score
                final.gapOp <- go
                final.gapExt <- ge
                final.reads <- read.scores
                final.scrambled <- scrambled.scores
            }
        }
    }

    list(parameters=list(gapOpening=final.gapOp, gapExtension=final.gapExt),
        scores=list(reads=final.reads, scrambled=final.scrambled))
}

.tied_overlap <- function(real, fake) 
# Calculating overlap (need average of left.open=TRUE/FALSE to handle ties).
{
    fake <- sort(fake)
    upper.bound <- findInterval(real, fake)
    lower.bound <- findInterval(real, fake, left.open=TRUE)
    sum((upper.bound + lower.bound)/2)/(length(real)*length(fake))
}

.align_TA_internal <- function(reads.start, reads.end, scrambled.start, scrambled.end, adaptor1, adaptor2, ...) 
# Wrapper to ensure that the sarlacc namespace is passed along in bpmapply.
{
    read.scores <- .get_alignment_scores(reads.start, reads.end, adaptor1, adaptor2, ...)
    scram.scores <- .get_alignment_scores(scrambled.start, scrambled.end, adaptor1, adaptor2, ...)
    list(
        reads=.resolve_strand(read.scores$START, read.scores$END, read.scores$RSTART, read.scores$REND)$scores,
        scrambled=.resolve_strand(scram.scores$START, scram.scores$END, scram.scores$RSTART, scram.scores$REND)$scores
    )
}

#' @importFrom Biostrings quality 
.get_alignment_scores <- function(reads.start, reads.end, adaptor1, adaptor2, gap.opening, gap.extension) 
# Retrieve all alignment scores for adaptor/read end combinations.
{
    FUN <- function(R, A) {
        Q <- quality(R)
        .Call(cxx_adaptor_align_score_only, R, Q, .create_encoding_vector(Q), gap.opening, gap.extension, A)
    }
    list(
        START=FUN(R=reads.start, A=adaptor1),
        END=FUN(R=reads.end, A=adaptor2),
        RSTART=FUN(R=reads.end, A=adaptor1),
        REND=FUN(R=reads.start, A=adaptor2)
    )
}
