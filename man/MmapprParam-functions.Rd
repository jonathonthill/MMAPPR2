% Manually synchronized with R/param.R.
\name{MmapprParam-functions}
\alias{MmapprParam-functions}
\alias{wtFiles}
\alias{wtFiles<-}
\alias{wtFiles,MmapprParam-method}
\alias{wtFiles<-,MmapprParam-method}
\alias{mutFiles}
\alias{mutFiles<-}
\alias{mutFiles,MmapprParam-method}
\alias{mutFiles<-,MmapprParam-method}
\alias{refFasta}
\alias{refFasta<-}
\alias{refFasta,MmapprParam-method}
\alias{refFasta<-,MmapprParam-method}
\alias{gtf}
\alias{gtf<-}
\alias{gtf,MmapprParam-method}
\alias{gtf<-,MmapprParam-method}
\alias{outputFolder}
\alias{outputFolder<-}
\alias{outputFolder,MmapprParam-method}
\alias{outputFolder<-,MmapprParam-method}
\alias{includeScaffolds}
\alias{includeScaffolds<-}
\alias{includeScaffolds,MmapprParam-method}
\alias{includeScaffolds<-,MmapprParam-method}
\alias{minDepth}
\alias{minDepth<-}
\alias{minDepth,MmapprParam-method}
\alias{minDepth<-,MmapprParam-method}
\alias{homozygoteCutoff}
\alias{homozygoteCutoff<-}
\alias{homozygoteCutoff,MmapprParam-method}
\alias{homozygoteCutoff<-,MmapprParam-method}
\alias{minBaseQuality}
\alias{minBaseQuality<-}
\alias{minBaseQuality,MmapprParam-method}
\alias{minBaseQuality<-,MmapprParam-method}
\alias{minMapQuality}
\alias{minMapQuality<-}
\alias{minMapQuality,MmapprParam-method}
\alias{minMapQuality<-,MmapprParam-method}
\alias{fileAggregation}
\alias{fileAggregation<-}
\alias{fileAggregation,MmapprParam-method}
\alias{fileAggregation<-,MmapprParam-method}
\alias{distancePower}
\alias{distancePower<-}
\alias{distancePower,MmapprParam-method}
\alias{distancePower<-,MmapprParam-method}
\alias{peakIntervalWidth}
\alias{peakIntervalWidth<-}
\alias{peakIntervalWidth,MmapprParam-method}
\alias{peakIntervalWidth<-,MmapprParam-method}
\alias{loessOptResolution}
\alias{loessOptResolution<-}
\alias{loessOptResolution,MmapprParam-method}
\alias{loessOptResolution<-,MmapprParam-method}
\alias{loessOptCutFactor}
\alias{loessOptCutFactor<-}
\alias{loessOptCutFactor,MmapprParam-method}
\alias{loessOptCutFactor<-,MmapprParam-method}
\alias{maxPileupDepth}
\alias{maxPileupDepth<-}
\alias{maxPileupDepth,MmapprParam-method}
\alias{maxPileupDepth<-,MmapprParam-method}
\alias{candidateMinDepth}
\alias{candidateMinDepth<-}
\alias{candidateMinDepth,MmapprParam-method}
\alias{candidateMinDepth<-,MmapprParam-method}
\alias{candidateMinAltDepth}
\alias{candidateMinAltDepth<-}
\alias{candidateMinAltDepth,MmapprParam-method}
\alias{candidateMinAltDepth<-,MmapprParam-method}
\alias{candidateMinAltFreq}
\alias{candidateMinAltFreq<-}
\alias{candidateMinAltFreq,MmapprParam-method}
\alias{candidateMinAltFreq<-,MmapprParam-method}
\alias{candidateMaxWtAltFreq}
\alias{candidateMaxWtAltFreq<-}
\alias{candidateMaxWtAltFreq,MmapprParam-method}
\alias{candidateMaxWtAltFreq<-,MmapprParam-method}
\alias{candidateMinDeltaAF}
\alias{candidateMinDeltaAF<-}
\alias{candidateMinDeltaAF,MmapprParam-method}
\alias{candidateMinDeltaAF<-,MmapprParam-method}
\alias{candidatePoolMode}
\alias{candidatePoolMode<-}
\alias{candidatePoolMode,MmapprParam-method}
\alias{candidatePoolMode<-,MmapprParam-method}
\alias{candidateChunkSize}
\alias{candidateChunkSize<-}
\alias{candidateChunkSize,MmapprParam-method}
\alias{candidateChunkSize<-,MmapprParam-method}
\alias{peakCutoffSd}
\alias{peakCutoffSd<-}
\alias{peakCutoffSd,MmapprParam-method}
\alias{peakCutoffSd<-,MmapprParam-method}
\alias{peakCutoffMethod}
\alias{peakCutoffMethod<-}
\alias{peakCutoffMethod,MmapprParam-method}
\alias{peakCutoffMethod<-,MmapprParam-method}
\alias{peakIntervalMethod}
\alias{peakIntervalMethod<-}
\alias{peakIntervalMethod,MmapprParam-method}
\alias{peakIntervalMethod<-,MmapprParam-method}
\alias{peakResampleIterations}
\alias{peakResampleIterations<-}
\alias{peakResampleIterations,MmapprParam-method}
\alias{peakResampleIterations<-,MmapprParam-method}
\alias{randomSeed}
\alias{randomSeed<-}
\alias{randomSeed,MmapprParam-method}
\alias{randomSeed<-,MmapprParam-method}
\alias{pairedEnd}
\alias{pairedEnd<-}
\alias{pairedEnd,MmapprParam-method}
\alias{pairedEnd<-,MmapprParam-method}
\alias{ignoreStrand}
\alias{ignoreStrand<-}
\alias{ignoreStrand,MmapprParam-method}
\alias{ignoreStrand<-,MmapprParam-method}
\alias{expressionPseudocount}
\alias{expressionPseudocount<-}
\alias{expressionPseudocount,MmapprParam-method}
\alias{expressionPseudocount<-,MmapprParam-method}
\alias{exportAiccPlots}
\alias{exportAiccPlots<-}
\alias{exportAiccPlots,MmapprParam-method}
\alias{exportAiccPlots<-,MmapprParam-method}
\title{MmapprParam getters and validated setters}
\usage{
\S4method{wtFiles}{MmapprParam}(obj)
\S4method{mutFiles}{MmapprParam}(obj)
\S4method{refFasta}{MmapprParam}(obj)
\S4method{gtf}{MmapprParam}(obj)
\S4method{outputFolder}{MmapprParam}(obj)
\S4method{includeScaffolds}{MmapprParam}(obj)
\S4method{minDepth}{MmapprParam}(obj)
\S4method{homozygoteCutoff}{MmapprParam}(obj)
\S4method{minBaseQuality}{MmapprParam}(obj)
\S4method{minMapQuality}{MmapprParam}(obj)
\S4method{fileAggregation}{MmapprParam}(obj)
\S4method{distancePower}{MmapprParam}(obj)
\S4method{peakIntervalWidth}{MmapprParam}(obj)
\S4method{loessOptResolution}{MmapprParam}(obj)
\S4method{loessOptCutFactor}{MmapprParam}(obj)
\S4method{maxPileupDepth}{MmapprParam}(obj)
\S4method{candidateMinDepth}{MmapprParam}(obj)
\S4method{candidateMinAltDepth}{MmapprParam}(obj)
\S4method{candidateMinAltFreq}{MmapprParam}(obj)
\S4method{candidateMaxWtAltFreq}{MmapprParam}(obj)
\S4method{candidateMinDeltaAF}{MmapprParam}(obj)
\S4method{candidatePoolMode}{MmapprParam}(obj)
\S4method{candidateChunkSize}{MmapprParam}(obj)
\S4method{peakCutoffSd}{MmapprParam}(obj)
\S4method{peakCutoffMethod}{MmapprParam}(obj)
\S4method{peakIntervalMethod}{MmapprParam}(obj)
\S4method{peakResampleIterations}{MmapprParam}(obj)
\S4method{randomSeed}{MmapprParam}(obj)
\S4method{pairedEnd}{MmapprParam}(obj)
\S4method{ignoreStrand}{MmapprParam}(obj)
\S4method{expressionPseudocount}{MmapprParam}(obj)
\S4method{exportAiccPlots}{MmapprParam}(obj)
\S4method{wtFiles}{MmapprParam}(obj) <- value
\S4method{mutFiles}{MmapprParam}(obj) <- value
\S4method{refFasta}{MmapprParam}(obj) <- value
\S4method{gtf}{MmapprParam}(obj) <- value
\S4method{outputFolder}{MmapprParam}(obj) <- value
\S4method{includeScaffolds}{MmapprParam}(obj) <- value
\S4method{minDepth}{MmapprParam}(obj) <- value
\S4method{homozygoteCutoff}{MmapprParam}(obj) <- value
\S4method{minBaseQuality}{MmapprParam}(obj) <- value
\S4method{minMapQuality}{MmapprParam}(obj) <- value
\S4method{fileAggregation}{MmapprParam}(obj) <- value
\S4method{distancePower}{MmapprParam}(obj) <- value
\S4method{peakIntervalWidth}{MmapprParam}(obj) <- value
\S4method{loessOptResolution}{MmapprParam}(obj) <- value
\S4method{loessOptCutFactor}{MmapprParam}(obj) <- value
\S4method{maxPileupDepth}{MmapprParam}(obj) <- value
\S4method{candidateMinDepth}{MmapprParam}(obj) <- value
\S4method{candidateMinAltDepth}{MmapprParam}(obj) <- value
\S4method{candidateMinAltFreq}{MmapprParam}(obj) <- value
\S4method{candidateMaxWtAltFreq}{MmapprParam}(obj) <- value
\S4method{candidateMinDeltaAF}{MmapprParam}(obj) <- value
\S4method{candidatePoolMode}{MmapprParam}(obj) <- value
\S4method{candidateChunkSize}{MmapprParam}(obj) <- value
\S4method{peakCutoffSd}{MmapprParam}(obj) <- value
\S4method{peakCutoffMethod}{MmapprParam}(obj) <- value
\S4method{peakIntervalMethod}{MmapprParam}(obj) <- value
\S4method{peakResampleIterations}{MmapprParam}(obj) <- value
\S4method{randomSeed}{MmapprParam}(obj) <- value
\S4method{pairedEnd}{MmapprParam}(obj) <- value
\S4method{ignoreStrand}{MmapprParam}(obj) <- value
\S4method{expressionPseudocount}{MmapprParam}(obj) <- value
\S4method{exportAiccPlots}{MmapprParam}(obj) <- value
}
\arguments{
\item{obj}{A \code{MmapprParam} object.}
\item{value}{Replacement value. The complete object is revalidated after assignment.}
}
\value{The requested parameter value; replacement methods return the updated object.}
\description{Accessors and validated replacement methods for \code{MmapprParam}. BAM replacements are indexed when necessary; replacing \code{refFasta} also refreshes the internal \code{FaFile}.}
\seealso{\code{\link{mmapprParam}}, \code{\linkS4class{MmapprParam}}}
