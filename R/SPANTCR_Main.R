# ---- SPANTCR Main Function
# INSTRUCTIONS:
# Prepare data with the following columns: CDR3/gene(TCRA/B)/score/vgene/jgene,id (if paired)
SPANTCR <- function(Data, Filename, FileType,  AAOperation="Hydrophobicity", 
                    TickNumber=100, SlidingWindowTickSize=2, SignificanceCutoff=0.03, 
                    ScoreFunction=ScoreFunctionExponential5, WeightFunction=WeightFunctionLinear)
{
  
  Metadata <- data.table("FileName"=Filename, "FileType"=FileType, "AAMetric"=AAOperation,
                         "Ticks"=TickNumber, "SlidingWindowTickSize"=SlidingWindowTickSize, "SignificanceCutoff"=SignificanceCutoff)
  
  Data[, score := ScoreFunction(score)]
  
  if(FileType=="Paired")
  {
    TRA <- Data[gene=="TRA"]
    colnames(TRA)[colnames(TRA)!="id"] <- paste0("TRA", colnames(TRA)[colnames(TRA)!="id"])
    TRB <- Data[gene=="TRB"]
    colnames(TRB)[colnames(TRB)!="id"] <- paste0("TRB", colnames(TRB)[colnames(TRB)!="id"])
    PairedData <- merge(TRA, TRB, "id")
    
    Clones <- unique(PairedData[, .(TRACDR3, TRAVgene, TRAJgene, TRBCDR3, TRBVgene, TRBJgene)])
    Clones[, Index := .I]
    Clones[, HitIndex := list()]
    Clonecopy <- copy(Clones)
    
    for (num in Clones$Index){
      set(Clones, i=num, "HitIndex", list(which(PairedData$TRACDR3==Clonecopy$TRACDR3[num] & PairedData$TRAVgene==Clonecopy$TRAVgene[num] & PairedData$TRAJgene==Clonecopy$TRAJgene[num] &
                                                  PairedData$TRBCDR3==Clonecopy$TRBCDR3[num] & PairedData$TRBVgene==Clonecopy$TRBVgene[num] & PairedData$TRBJgene==Clonecopy$TRBJgene[num])))
    }
    Clonecopy <- copy(Clones)
    
    for (num in Clones$Index){
      set(Clones, i=num, c("TRAScore","TRBScore"), list(sum(PairedData[Clonecopy$HitIndex[[num]]]$TRAscore),sum(PairedData[Clonecopy$HitIndex[[num]]]$TRBscore)))
    }
    
    Clones[, Score := mean(c(TRAScore, TRBScore)), by=Index]
    AlphaBeta <- Clones[, .(TRACDR3, TRBCDR3)]
    
    CompleteSlidingWindows1 <- apply(AlphaBeta, 1, function(x) 
      c(strsplit(x[1],NULL)[[1]], strsplit(x[2],NULL)[[1]]))
    CompleteSlidingPositions1 <- apply(AlphaBeta, 1, function(x) 
      c(seq(1/nchar(x[1])/2,by=1/nchar(x[1])),
        1+seq(1/nchar(x[2])/2,by=1/nchar(x[2]))))
    CompleteSlidingRanges1 <- apply(AlphaBeta, 1, function(x) 
      c(mapply(c,(seq(1/nchar(x[1]),by=1/nchar(x[1]))-1/nchar(x[1])),
               seq(1/nchar(x[1]),by=1/nchar(x[1])),SIMPLIFY=F),
        mapply(c,(1+seq(1/nchar(x[2]),by=1/nchar(x[2]))-1/nchar(x[2])),
               1+seq(1/nchar(x[2]),by=1/nchar(x[2])),SIMPLIFY=F)))
    
    CompleteSlidingWindows2 <- apply(AlphaBeta, 1, function(x)
      c(paste0(strsplit(x[1],NULL)[[1]], 
               strsplit(x[1],NULL)[[1]][-1])[-nchar(x[1])],
        c(paste0(strsplit(x[2],NULL)[[1]], 
                 strsplit(x[2],NULL)[[1]][-1])[-nchar(x[2])])))
    CompleteSlidingPositions2 <- apply(AlphaBeta, 1, function(x) 
      c(seq(1/(nchar(x[1])-1)/2,by=1/(nchar(x[1])-1))[-nchar(x[1])],
        1+seq(1/(nchar(x[2])-1)/2,by=1/(nchar(x[2])-1))[-nchar(x[2])]))
    CompleteSlidingRanges2 <- apply(AlphaBeta, 1, function(x) 
      c(mapply(c,(seq(1/nchar(x[1]),by=1/nchar(x[1]))-1/nchar(x[1]))[-nchar(x[1])],
               (seq(1/nchar(x[1]),by=1/nchar(x[1]))+1/nchar(x[1]))[-nchar(x[1])],SIMPLIFY=F),
        mapply(c,(1+seq(1/nchar(x[2]),by=1/nchar(x[2]))-1/nchar(x[2]))[-nchar(x[2])],
               (1+seq(1/nchar(x[2]),by=1/nchar(x[2]))+1/nchar(x[2]))[-nchar(x[2])],SIMPLIFY=F)))
    
    CompleteSlidingWindows3 <- apply(AlphaBeta, 1, function(x)
      c(paste0(strsplit(x[1],NULL)[[1]],
               strsplit(x[1],NULL)[[1]][-1],
               strsplit(x[1],NULL)[[1]][-c(1,2)])[c(-nchar(x[1]), -(nchar(x[1])-1))],
        paste0(strsplit(x[2],NULL)[[1]],
               strsplit(x[2],NULL)[[1]][-1],
               strsplit(x[2],NULL)[[1]][-c(1,2)])[c(-nchar(x[2]), -(nchar(x[2])-1))]))
    CompleteSlidingPositions3 <- apply(AlphaBeta, 1, function(x)
      c(seq(1/(nchar(x[1])-2)/2,by=1/(nchar(x[1])-2))[c(-nchar(x[1]), -(nchar(x[1])-1))],
        1+seq(1/(nchar(x[2])-2)/2,by=1/(nchar(x[2])-2))[c(-nchar(x[2]), -(nchar(x[2])-1))]))
    CompleteSlidingRanges3 <- apply(AlphaBeta, 1, function(x)
      c(mapply(c,(seq(1/nchar(x[1]),by=1/nchar(x[1]))-1/nchar(x[1]))[c(-nchar(x[1]), -(nchar(x[1])-1))],
               (seq(1/nchar(x[1]),by=1/nchar(x[1]))+2/nchar(x[1]))[c(-nchar(x[1]), -(nchar(x[1])-1))],SIMPLIFY=F),
        mapply(c,(1+seq(1/nchar(x[2]),by=1/nchar(x[2]))-1/nchar(x[2]))[c(-nchar(x[2]), -(nchar(x[2])-1))],
               (1+seq(1/nchar(x[2]),by=1/nchar(x[2]))+2/nchar(x[2]))[c(-nchar(x[2]), -(nchar(x[2])-1))],SIMPLIFY=F)))
    
    CombinedPairedData <- data.table(CompleteSlidingWindows1=I(CompleteSlidingWindows1), CompleteSlidingPositions1=I(CompleteSlidingPositions1),
                                     CompleteSlidingRanges1=I(CompleteSlidingRanges1),
                                     CompleteSlidingWindows2=I(CompleteSlidingWindows2), CompleteSlidingPositions2=I(CompleteSlidingPositions2), 
                                     CompleteSlidingRanges2=I(CompleteSlidingRanges2),
                                     CompleteSlidingWindows3=I(CompleteSlidingWindows3), CompleteSlidingPositions3=I(CompleteSlidingPositions3), 
                                     CompleteSlidingRanges3=I(CompleteSlidingRanges3))
    Clones <- cbind(Clones, CombinedPairedData)
    
    SlidingWindowInput <- Clones[[paste0("CompleteSlidingWindows",SlidingWindowTickSize)]]
    SlidingWindowList <- unlist(SlidingWindowInput)
    
    SlidingRangesInput <- Clones[[paste0("CompleteSlidingRanges",SlidingWindowTickSize)]]
    
    WindowRange <- unique(SlidingWindowList)
    
    # Prepare single vector of frequency for all windows
    UnlistedIndexes <- c()
    for (n in 1:length(SlidingWindowInput))
    {
      UnlistedIndexes <- c(UnlistedIndexes, rep.int(n, length(SlidingWindowInput[[n]])))
    }
    UnlistedUniqueInputs <- UnlistedFrequencyMap
    UnlistedFrequencyMap <- Clones$Score[UnlistedIndexes]
    UnlistedIndex <- Clones$Index[UnlistedIndexes]
    UnlistedAlphaVMap <- Clones$TRAVgene[UnlistedIndexes]
    UnlistedBetaVMap <- Clones$TRBVgene[UnlistedIndexes]
    UnlistedAlphaJMap <- Clones$TRAJgene[UnlistedIndexes]
    UnlistedBetaJMap <- Clones$TRBJgene[UnlistedIndexes]
    
    
    
    # Determine elements that appear and indexing to reorder
    LogoTickElements <- WindowRange
    # Input operations and indexing here -----
    LogoTickElementsOperation <- unlist(lapply(strsplit(LogoTickElements,NULL),
                                               function(x) sum(as.numeric(AminoAcidFilter[[AAOperation]]
                                                                          [match(x,AminoAcidFilter[["AA"]])]))))
    
    TickOperationIndex <- order(LogoTickElementsOperation)
    LogoTickElementsReindex <- LogoTickElements[TickOperationIndex]
    LogoTickElementsOperationReindex <- LogoTickElementsOperation[TickOperationIndex]
    
    # Prepare bins
    FineBinTicks <- seq(0, 2, by=1/TickNumber)[-(2*TickNumber+1)]+1/TickNumber/2
    FineBinLabels <- paste0("Bin", 1:(2*TickNumber))
    
    OutputTickData <- data.table(Window=factor(), ChainsVgene=factor(), ChainsJgene=factor(),
                                 Probability=numeric(), WeightProbability=numeric(), Tick=numeric(), Entropy=numeric(),
                                 WeightEntropy=numeric(), Max=numeric(), WeightMax=numeric(), Origin=factor())[1:(2*(SlidingWindowTickSize*TickNumber)*nrow(Clones))]
    
    counter <- 1
    for(tick in 1:length(FineBinTicks))
    {
      if(tick<=length(FineBinTicks)/2){
        UnlistedChainsVMap <- UnlistedAlphaVMap
        UnlistedChainsJMap <- UnlistedAlphaJMap
      } else {
        UnlistedChainsVMap <- UnlistedBetaVMap
        UnlistedChainsJMap <- UnlistedBetaJMap
      }
      TempBinningList <- unlist(lapply(SlidingRangesInput, function(rangelist) lapply(rangelist, function(bounds) (FineBinTicks[tick] > bounds[1] & FineBinTicks[tick] < bounds[2]))))
      
      TempBinningWeight <- unlist(lapply(SlidingRangesInput, function(rangelist) lapply(rangelist, function(bounds) WeightFunction((FineBinTicks[tick]-bounds[1])/(bounds[2]-bounds[1])))))
      
      Hits <- which(TempBinningList)
      Range <- counter:(counter+length(Hits)-1)
      
      UniqueWindows <- unique(SlidingWindowList[Hits])
      OutputTickData[Range, Window := SlidingWindowList[Hits]]
      OutputTickData[Range, ChainsVgene := UnlistedChainsVMap[Hits]]
      OutputTickData[Range, ChainsJgene := UnlistedChainsJMap[Hits]]
      HitWeights <- TempBinningWeight[TempBinningList]
      
      HitProbabilities <- UnlistedFrequencyMap[Hits]/sum(UnlistedFrequencyMap[Hits])
      OutputTickData[Range, Probability := HitProbabilities]
      HitWeightProbabilities <- UnlistedFrequencyMap[Hits]*HitWeights/sum(UnlistedFrequencyMap[Hits]*HitWeights)
      OutputTickData[Range, WeightProbability := HitWeightProbabilities]
      OutputTickData[Range, Tick := FineBinTicks[tick]]
      OutputTickData[Range, Entropy := rep(sum(sapply(sapply(unique(UniqueWindows), 
                                                             function(window) sum(HitProbabilities[OutputTickData[Range]$Window==window])), 
                                                      function(p) -p*log(p))), length(Hits))]
      OutputTickData[Range, WeightEntropy := rep(sum(sapply(sapply(unique(UniqueWindows), 
                                                                   function(window) sum(HitWeightProbabilities[OutputTickData[Range]$Window==window])), 
                                                            function(p) -p*log(p))), length(Hits))]
      OutputTickData[Range, Origin := factor(UnlistedIndex[Hits])]
      MaxList <- sapply(UniqueWindows, function(x) sum(UnlistedFrequencyMap[Hits][SlidingWindowList[Hits]==x])/sum(UnlistedFrequencyMap[Hits]))
      MaxWeightList <- sapply(UniqueWindows, function(x) sum((UnlistedFrequencyMap[Hits]*HitWeights)[SlidingWindowList[Hits]==x])/sum(UnlistedFrequencyMap[Hits]*HitWeights))
      OutputTickData[Range, Max := MaxList[match(Window,names(MaxList))]]
      OutputTickData[Range, WeightMax := MaxWeightList[match(Window,names(MaxWeightList))]]
      
      counter <- counter+length(Hits)
    }
    
    
    # Input color scaling variable here -----
    ColorTickScaleInput <- LogoTickElementsOperationReindex
    ColorTickScaleNormalize <- (ColorTickScaleInput-min(ColorTickScaleInput))/(max(ColorTickScaleInput)-min(ColorTickScaleInput))*0.99+0.01
    ColorTickScaleRampPalette <- c("#4D4D4D",colorRampPalette(c("red","white","blue"))(100))
    ColorTickScaleColors <- ColorTickScaleRampPalette[100*round(ColorTickScaleNormalize,2)]
    
    
    OutputTickData$Window <- factor(OutputTickData$Window, 
                                    levels=rev(levels(OutputTickData$Window)[match(LogoTickElementsReindex,
                                                                                   levels(OutputTickData$Window))]))
    OutputTickData$Color <- LogoTickElementsOperationReindex[match(OutputTickData$Window,LogoTickElementsReindex)]
    
    
    OutputTickData[, SignificantOrigin := ifelse(WeightProbability < SignificanceCutoff, 0, Origin)]
    OutputTickData[, SignificantWindow := ifelse(WeightMax < SignificanceCutoff, 0, Window)]
    OutputTickData[, SignificantColor := ifelse(WeightMax < SignificanceCutoff, 0,
                                                100*round(as.numeric((
                                                  (Color-min(ColorTickScaleInput))/
                                                    (max(ColorTickScaleInput)-min(ColorTickScaleInput))*0.99+0.01)),2))]
    
    
    OutputTickData <- OutputTickData[!is.na(Window)]
    PairedCDR3Complete <- CDR3Breakdown(CloneData=Clones, Elements=LogoTickElementsReindex,
                                        Metric=LogoTickElementsOperationReindex, Output=OutputTickData, Metadata=Metadata)
    
  }
  if(FileType %in% c("TRA","TRB"))
  {
    Data <- Data[gene==FileType]
    Clones <- unique(Data[, .(CDR3, Vgene, Jgene)])
    Clones[, Index := .I]
    Clones[, HitIndex := list()]
    Clonecopy <- copy(Clones)
    
    for (num in Clones$Index){
      set(Clones, i=num, "HitIndex", list(which(Data$CDR3==Clonecopy$CDR3[num] & Data$Vgene==Clonecopy$Vgene[num] & Data$Jgene==Clonecopy$Jgene[num])))
    }
    Clonecopy <- copy(Clones)
    
    for (num in Clones$Index){
      set(Clones, i=num, "Score", sum(Data[Clonecopy$HitIndex[[num]]]$score))
    }
    
    CompleteSlidingWindows1 <- sapply(Clones$CDR3, function(x) strsplit(x,NULL)[[1]], USE.NAMES = F)
    CompleteSlidingPositions1 <- sapply(Clones$CDR3, function(x) seq(1/nchar(x)/2,by=1/nchar(x)), USE.NAMES = F)
    CompleteSlidingRanges1 <- sapply(Clones$CDR3, function(x) 
      mapply(c,(seq(1/nchar(x),by=1/nchar(x))-1/nchar(x)),
             seq(1/nchar(x),by=1/nchar(x)),SIMPLIFY=F), USE.NAMES = F)
    
    CompleteSlidingWindows2 <- sapply(Clones$CDR3, function(x)
      c(paste0(strsplit(x,NULL)[[1]], 
               strsplit(x,NULL)[[1]][-1])[-nchar(x)]), USE.NAMES = F)
    CompleteSlidingPositions2 <- sapply(Clones$CDR3, function(x) 
      seq(1/(nchar(x)-1)/2,by=1/(nchar(x)-1))[-nchar(x)], USE.NAMES = F)
    CompleteSlidingRanges2 <- sapply(Clones$CDR3, function(x) 
      mapply(c,(seq(1/nchar(x),by=1/nchar(x))-1/nchar(x))[-nchar(x)],
             (seq(1/nchar(x),by=1/nchar(x))+1/nchar(x))[-nchar(x)],SIMPLIFY=F), USE.NAMES = F)
    
    CompleteSlidingWindows3 <- sapply(Clones$CDR3, function(x) 
      paste0(strsplit(x,NULL)[[1]], 
             strsplit(x,NULL)[[1]][-1],
             strsplit(x,NULL)[[1]][-c(1,2)])[c(-nchar(x), -(nchar(x)-1))], USE.NAMES = F)
    CompleteSlidingPositions3 <- sapply(Clones$CDR3, function(x) 
      seq(1/(nchar(x)-2)/2,by=1/(nchar(x)-2))[c(-nchar(x), -(nchar(x)-1))], USE.NAMES = F)
    CompleteSlidingRanges3 <- sapply(Clones$CDR3, function(x) 
      mapply(c,(seq(1/nchar(x),by=1/nchar(x))-1/nchar(x))[c(-nchar(x), -(nchar(x)-1))],
             (seq(1/nchar(x),by=1/nchar(x))+2/nchar(x))[c(-nchar(x), -(nchar(x)-1))],SIMPLIFY=F), USE.NAMES=F)
    
    CombinedPairedData <- data.table(CompleteSlidingWindows1=I(CompleteSlidingWindows1), CompleteSlidingPositions1=I(CompleteSlidingPositions1),
                                     CompleteSlidingRanges1=I(CompleteSlidingRanges1),
                                     CompleteSlidingWindows2=I(CompleteSlidingWindows2), CompleteSlidingPositions2=I(CompleteSlidingPositions2), 
                                     CompleteSlidingRanges2=I(CompleteSlidingRanges2),
                                     CompleteSlidingWindows3=I(CompleteSlidingWindows3), CompleteSlidingPositions3=I(CompleteSlidingPositions3), 
                                     CompleteSlidingRanges3=I(CompleteSlidingRanges3))
    Clones <- cbind(Clones, CombinedPairedData)
    
    # Input analysis parameter here -----
    SlidingWindowInput <- Clones[[paste0("CompleteSlidingWindows",SlidingWindowTickSize)]]
    SlidingWindowList <- unlist(SlidingWindowInput)
    SlidingRangesInput <- Clones[[paste0("CompleteSlidingRanges",SlidingWindowTickSize)]]
    
    
    WindowRange <- unique(SlidingWindowList)
    
    # Prepare single vector of frequency for all windows
    UnlistedIndexes <- vector("integer",length(unlist(SlidingWindowInput)))
    count <- 1
    for (n in 1:length(SlidingWindowInput))
    {
      UnlistedIndexes[count:(count+length(SlidingWindowInput[[n]])-1)] <- n
      count <- count+length(SlidingWindowInput[[n]])
    }
    UnlistedFrequencyMap <- Clones$Score[UnlistedIndexes]
    UnlistedIndex <- Clones$Index[UnlistedIndexes]
    UnlistedChainsVMap <- Clones$Vgene[UnlistedIndexes]
    UnlistedChainsJMap <- Clones$Jgene[UnlistedIndexes]
    
    # Determine elements that appear and indexing to reorder
    LogoTickElements <- WindowRange
    # Input operations and indexing here -----
    LogoTickElementsOperation <- unlist(lapply(strsplit(LogoTickElements,NULL),
                                               function(x) sum(as.numeric(AminoAcidFilter[[AAOperation]]
                                                                          [match(x,AminoAcidFilter[["AA"]])]))))
    
    TickOperationIndex <- order(LogoTickElementsOperation)
    LogoTickElementsReindex <- LogoTickElements[TickOperationIndex]
    LogoTickElementsOperationReindex <- LogoTickElementsOperation[TickOperationIndex]
    
    
    # Prepare bins
    FineBinTicks <- seq(0, 1, by=1/TickNumber)[-(TickNumber+1)]+1/TickNumber/2
    FineBinLabels <- paste0("Bin", 1:TickNumber)
    
    OutputTickData <- data.table(Window=factor(), ChainsVgene=factor(), ChainsJgene=factor(),
                                 Probability=numeric(), WeightProbability=numeric(), Tick=numeric(), Entropy=numeric(),
                                 WeightEntropy=numeric(), Max=numeric(), WeightMax=numeric(), Origin=factor())[1:((SlidingWindowTickSize*TickNumber)*nrow(Clones))]#1:sum(sapply(FullBinningList, function(x) sum(unlist(x))))]
    
    counter <- 1
    for(tick in 1:length(FineBinTicks))
    {
      TempBinningList <- unlist(lapply(SlidingRangesInput, function(rangelist) lapply(rangelist, function(bounds) (FineBinTicks[tick] > bounds[1] & FineBinTicks[tick] < bounds[2]))))
      # TempBinningList <- unlist(FullBinningList[[tick]])
      TempBinningWeight <- unlist(lapply(SlidingRangesInput, function(rangelist) lapply(rangelist, function(bounds) WeightFunction((FineBinTicks[tick]-bounds[1])/(bounds[2]-bounds[1])))))
      # TempBinningWeight <- unlist(FullBinningWeight[[tick]])
      Hits <- which(TempBinningList)
      Range <- counter:(counter+length(Hits)-1)
      
      UniqueWindows <- unique(SlidingWindowList[Hits])
      OutputTickData[Range, Window := SlidingWindowList[Hits]]
      OutputTickData[Range, ChainsVgene := UnlistedChainsVMap[Hits]]
      OutputTickData[Range, ChainsJgene := UnlistedChainsJMap[Hits]]
      HitWeights <- TempBinningWeight[TempBinningList]
      
      HitProbabilities <- UnlistedFrequencyMap[Hits]/sum(UnlistedFrequencyMap[Hits])
      OutputTickData[Range, Probability := HitProbabilities]
      HitWeightProbabilities <- UnlistedFrequencyMap[Hits]*HitWeights/sum(UnlistedFrequencyMap[Hits]*HitWeights)
      OutputTickData[Range, WeightProbability := HitWeightProbabilities]
      OutputTickData[Range, Tick := FineBinTicks[tick]]
      OutputTickData[Range, Entropy := rep(sum(sapply(sapply(unique(UniqueWindows), 
                                                             function(window) sum(HitProbabilities[OutputTickData[Range]$Window==window])), 
                                                      function(p) -p*log(p))), length(Hits))]
      OutputTickData[Range, WeightEntropy := rep(sum(sapply(sapply(unique(UniqueWindows), 
                                                                   function(window) sum(HitWeightProbabilities[OutputTickData[Range]$Window==window])), 
                                                            function(p) -p*log(p))), length(Hits))]
      OutputTickData[Range, Origin := factor(UnlistedIndex[Hits])]
      MaxList <- sapply(UniqueWindows, function(x) sum(UnlistedFrequencyMap[Hits][SlidingWindowList[Hits]==x])/sum(UnlistedFrequencyMap[Hits]))
      MaxWeightList <- sapply(UniqueWindows, function(x) sum((UnlistedFrequencyMap[Hits]*HitWeights)[SlidingWindowList[Hits]==x])/sum(UnlistedFrequencyMap[Hits]*HitWeights))
      OutputTickData[Range, Max := MaxList[match(Window,names(MaxList))]]
      OutputTickData[Range, WeightMax := MaxWeightList[match(Window,names(MaxWeightList))]]
      
      counter <- counter+length(Hits)
      
    }
    
    OutputTickData <- OutputTickData[!is.na(Window)]
    
    # Input color scaling variable here -----
    ColorTickScaleInput <- LogoTickElementsOperationReindex
    ColorTickScaleNormalize <- (ColorTickScaleInput-min(ColorTickScaleInput))/(max(ColorTickScaleInput)-min(ColorTickScaleInput))*0.99+0.01
    ColorTickScaleRampPalette <- c("#4D4D4D",colorRampPalette(c("red","white","blue"))(100))
    ColorTickScaleColors <- ColorTickScaleRampPalette[100*round(ColorTickScaleNormalize,2)]
    
    
    OutputTickData[, Window := factor(OutputTickData$Window, 
                                      levels=rev(levels(OutputTickData$Window)[match(LogoTickElementsReindex,
                                                                                     levels(OutputTickData$Window))]))]
    OutputTickData[, Color := LogoTickElementsOperationReindex[match(OutputTickData$Window,LogoTickElementsReindex)]]
    
    
    OutputTickData[, SignificantOrigin := ifelse(WeightProbability < SignificanceCutoff, 0, Origin)]
    OutputTickData[, SignificantWindow := ifelse(WeightMax < SignificanceCutoff, 0, Window)]
    OutputTickData[, SignificantColor := ifelse(WeightMax < SignificanceCutoff, 0, 
                                                100*round(as.numeric((
                                                  (Color-min(ColorTickScaleInput))/
                                                    (max(ColorTickScaleInput)-min(ColorTickScaleInput))*0.99+0.01)),2))]
    
    
    
    OutputTickData <- OutputTickData[!is.na(Window)]
    
    PairedCDR3Complete <- CDR3Breakdown(CloneData=Clones, Elements=LogoTickElementsReindex,
                                        Metric=LogoTickElementsOperationReindex,
                                        Output=OutputTickData, Metadata=Metadata)
    
  }
  
  gc()
  PairedCDR3Complete
}
