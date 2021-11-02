#' Checks entropy contributions for all TCRs expressing a specific motif (searchterms) in a specific range (targetrange)
#'
#' @param base CDR3Breakdown object
#' @param searchterms Vector of k-mers to analyze matching the length analyzed in base
#' @param targetrange Two element vector describing the lower and upper bound of CDR3 region analyzed (from 0-1)
#' @param ticks Vector of ticks to analyze, FineBinTicks for single CDR3, FineBinTicksPaired for paired CDR3s
#'
#' @return data.table object describing per-bin entropy of k-mers containing a specific k-mer (searchterms) within a specific region (targetrange)
#' @export
#'
#' @examples \dontrun{SearchIterator(CDR3Breakdown, c("GG","GA"), c(0.4,0.6), FineBinTicks)}
SearchIterator <- function(base, searchterms, targetrange, ticks=FineBinTicks)
{
  BaseEntropy <- base@Output[, mean(WeightEntropy), by=Tick]$V1
  FullEntropyList <- data.table::data.table(Location=ticks, Entropy=BaseEntropy, Source="Base", N=nrow(base@CloneData), Range=paste0(targetrange, collapse="-"))
  for(n in searchterms)
  {
    HitList <- unique(base@Output[Tick>targetrange[1] & Tick<targetrange[2] & Window==n]$Origin)
    SearchData <- base@Output[Origin %in% HitList]
    UniqueWindows <- sapply(ticks, function(tick) unique(SearchData[Tick==tick]$Window))
    ProbabilityData <- SearchData[, list(Probability=sum(WeightProbability)), by=.(Tick, Window)]
    EntropyData <- ProbabilityData[, list(Entropy=sum(sapply(.SD/sum(.SD), function(x) -x*log(x)))), by=Tick, .SDcols="Probability"]
    FullEntropyList <- rbind(FullEntropyList, data.table::data.table(Location=ticks, Entropy=EntropyData$Entropy, Source=I(n), N=length(HitList), Range=paste0(targetrange, collapse="-")))
  }
  FullEntropyList
}


#' Checks entropy contributions for each major motif in each range window
#'
#' @param Data CDR3Breakdown object
#' @param Searches List of vectors containing search ranges, default SearchBoxes for single CDR3 SearchBoxesPaired for paired
#' @param Limit Minimum frequency of k-mers to analyze in each bin, default 0.03
#' @param ticksused Vector of ticks to analyze, FineBinTicks for single CDR3, FineBinTicksPaired for paired CDR3s
#'
#' @return List of repeated SearchIterator outputs for each search range (Searches), and summarized data by k-mers. Summarized data describes essential k-mers and their location.
#' @export
#'
#' @examples \dontrun{EntropyScan(CD3Breakdown, list(c(0.4,0.5), c(0.5,0.6)), 0.03, FineBinTicks)}
EntropyScan <- function(Data, Searches, Limit=0.03, ticksused=FineBinTicks)
{
  EntropySummary <- data.table::data.table()
  EntropyRaw <- data.table::data.table()
  if(!is.null(Data)){
    for(range in Searches)
    {
      hitsinrange <- sum(table(Data@Output[Tick>range[1] & Tick<range[2]]$Window))
      targets <- which(table(Data@Output[Tick>range[1] & Tick<range[2]]$Window)>(hitsinrange*Limit))
      output <- SearchIterator(Data, names(targets), range, ticksused)
      NewEntropySummary <- output[, .(Average=mean(Entropy), Count=unique(N), Range=unique(Range)), by=Source]
      EntropySummary <- rbind(EntropySummary, NewEntropySummary)
      EntropyRaw <- rbind(EntropyRaw, output)
    }
    data.table::set(EntropySummary, j="DeltaAverage", value=max(EntropySummary$Average)-EntropySummary$Average)
  }
  list(EntropySummary, EntropyRaw)
}

