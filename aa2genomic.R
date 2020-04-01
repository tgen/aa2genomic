##amino acid position to genomic coordinates
##by Natalia Briones

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ensembldb")
BiocManager::install("AnnotationHub")
BiocManager::install("Gviz")
library(ensembldb)     
library(AnnotationHub)
library(Gviz)
ah <- AnnotationHub()
qr <- query(ah, c("EnsDb", "Canis lupus familiaris", "98"))
edbx <- qr[[1]]
genes(edbx)

aaHotspots = read.delim("/path/to/aahotspots.txt", header = TRUE)
aaHotspots_df <- as.data.frame(aaHotspots)
aaHotspots_rng <- IRanges(start = aaHotspots_df$aa_start, end = aaHotspots_df$aa_end,
                     names = aaHotspots_df$protein_id)
aaHotspots_gnm <- proteinToGenome(aaHotspots_rng, edbx)
warnings()
aaHotspots_gnm
aaHotspots_gnm_grng <- unlist(GRangesList(aaHotspots_gnm))
write.table(aaHotspots_gnm_grng,file="/path/to/aaHotspotsGenomicCoordinates.txt",sep = "\t", row.names = FALSE)


ABL1_rng <- IRanges(start=252, end=252, names= "ENSCAFP00000029534")
ABL1_gnm <- proteinToGenome(ABL1_rng, edbx)
ABL1_gnm


## Define a genome axis track
gat <- GenomeAxisTrack()

## Get the transcript ID:
txid <- ABL1_gnm[[1]]$tx_id[1]

## Get a GRanges for the transcript
trt <- getGeneRegionTrackForGviz(edbx, filter = TxIdFilter(txid))

## Define a GRanges for the mapped protein domains and add
## metadata columns with the grouping of the ranges and the
## IDs of the corresponding protein domains, so they can be
## identified in the plot
pdoms <- unlist(GRangesList(ABL1_gnm))
pdoms$grp <- rep(1:length(ABL1_rng), lengths(ABL1_gnm))
pdoms$id <- rep(mcols(ABL1_rng)$protein_domain_id, lengths(ABL1_gnm))

## Since we're using Ensembl chromosome names we have to set
options(ucscChromosomeNames = FALSE)

## Plotting the transcript and the mapped protein domains.
plotTracks(list(gat,
                GeneRegionTrack(trt, name = "tx"),
                AnnotationTrack(pdoms, group = pdoms$grp,
                                id = pdoms$id,
                                groupAnnotation = "id",
                                just.group = "above",
                                shape = "box",
                                name = "Protein domains")),
           transcriptAnnotation = "transcript")

## Define the individual tracks:
## - Ideogram
ideo_track <- IdeogramTrack(genome = "canFam3", chromosome = "chr9")
## - Genome axis
gaxis_track <- GenomeAxisTrack()
## - Transcripts
gene_track <- GeneRegionTrack(txs, showId = TRUE, just.group = "right",
                              name = "", geneSymbol = TRUE, size = 0.5)
## - Protein domains
pdom_track <- AnnotationTrack(pdoms_gnm_grng, group = pdoms_gnm_grng$grp,
                              id = pdoms_gnm_grng$id, groupAnnotation = "id",
                              just.group = "right", shape = "box",
                              name = "Protein domains", size = 0.5)

## Generate the plot
plotTracks(list(ideo_track, gaxis_track, gene_track, pdom_track))


rm(list = ls())

