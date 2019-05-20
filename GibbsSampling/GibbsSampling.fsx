

#r @"..\GibbsSampling\bin\Debug\netstandard2.0\FSharp.Plotly.dll"
#r @"..\GibbsSampling\bin\Debug\netstandard2.0\FSharpAux.dll"
#r @"..\GibbsSampling\bin\Debug\netstandard2.0\FSharpAux.IO.dll"
#r @"..\GibbsSampling\bin\Debug\netstandard2.0\BioFSharp.dll"
#r @"..\GibbsSampling\bin\Debug\netstandard2.0\BioFSharp.Stats.dll"
#r @"..\GibbsSampling\bin\Debug\netstandard2.0\BioFSharp.IO.dll"
#r @"..\GibbsSampling\bin\Debug\netstandard2.0\GibbsSampling.dll"

open System
open FSharp.Plotly
open BioFSharp
open BioFSharp.Stats
open BioFSharp.IO
open BioFSharp.Nucleotides
open FSharpAux
open FSharpAux.IO
open GibbsSampling
open HelperFunctionsAndTypes
open PositionFrequencyMatrix
open PositionWeightMatrix
open FrequencyCompositeVector
open PositionProbabilityMatrix
open ProbabilityCompositeVector


module SiteSampler =

    /// Gives the startPosition and score of the segment with the highest PositionWeightMatrixScore based on the given sequence and PositionWeightMatrix.
    let getBestPWMSsWithBPV (motifLength:int) (alphabet:#IBioItem[]) (source:BioArray.BioArray<#IBioItem>) (pcv:ProbabilityCompositeVector) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        let rec loop (n:int) (highValue:float) (highIndex:int) =
            if n + motifLength > source.Length then log2(highValue), highIndex
            else
                let tmp =
                    let segment =
                        Array.skip n source
                        |> Array.take motifLength
                    let pwMatrix = createPositionWeightMatrix alphabet pcv positionProbabilityMatrix
                    segment
                    |> calculateSegmentScore pwMatrix
                if tmp > highValue then loop (n + 1) tmp n
                else loop (n + 1) highValue highIndex
        loop 0 0. 0

    /// Checks whether downstream of given positions a higher InformationContent is present or not. 
    /// If yes, the new InformationContent and positions are given back, otherwise the old ones.
    let getRightShiftedBestPWMSsWithBPV (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (startPositions:(float*int)[]) =
        let rec loop (n:int) (acc:(float*int)[]) (bestMotif:(float*int)[]) =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> snd item)) = (bestMotif |> Array.map (fun item -> snd item)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    let tmp = Array.append sources.[0..n-1] sources.[n+1..]
                    Array.append bestMotif.[0..n-1] bestMotif.[n+1..]
                    |> Array.map2 (fun (source:BioArray.BioArray<#IBioItem>) (_, position) -> if position <= source.Length - motifLength - 1 then position + 1 else position) tmp
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        (getDefinedSegment motifLength subSequence position) 
                         |> getPositionFrequencyMatrix) unChosenArrays unChosenStartPositions
                    |> fusePositionFrequencyMatrices motifLength
                    |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSsWithBPV motifLength alphabet sources.[n] pcv positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[n] then 
                        acc.[n] <- tmp
                        acc
                     else acc                    
                    )
                    bestMotif
        loop 0 (Array.copy startPositions) (Array.copy startPositions)

    /// Checks whether upstream of given positions a higher PositionWeightMatrixScore is present or not. 
    /// If yes, the new PositionWeightMatrixScore and positions are given back, otherwise the old ones.
    let getLeftShiftedBestPWMSsWithBPV (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (startPositions:(float*int)[]) =
        let rec loop (n:int) (acc:(float*int)[]) (bestMotif:(float*int)[]) =
            if n = sources.Length then
                if (acc |> Array.map (fun item -> snd item)) = (bestMotif |> Array.map (fun item -> snd item)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    Array.append bestMotif.[0..n-1] bestMotif.[n+1..]
                    |> Array.map (fun (_, position) -> if position > 0 then position - 1 else position)
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        (getDefinedSegment motifLength subSequence position)
                         |> getPositionFrequencyMatrix) unChosenArrays unChosenStartPositions
                    |> fusePositionFrequencyMatrices motifLength
                    |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSsWithBPV motifLength alphabet sources.[n] pcv positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[n] then 
                        acc.[n] <- tmp
                        acc
                     else acc                    
                    )
                    bestMotif
        loop 0 (Array.copy startPositions) (Array.copy startPositions)

    /// Checks the given Sequence for the existence of a conserved motif, by scoring each segment based on the given start positions.
    /// The new PositionWeightMatrix is calculated and updated each iteration if segments with higher scores are found until convergence.
    let getBestPWMSsWithStartPositionsWithBPV (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (startPositions:(float*int)[]) =        
        let rec loop (n:int) acc bestMotif =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> snd item)) = (bestMotif |> Array.map (fun item -> snd item)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    Array.append acc.[0..n-1] acc.[n+1..]
                    |> Array.map (fun (_, position) -> position)
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        (getDefinedSegment motifLength subSequence position)
                         |> getPositionFrequencyMatrix) unChosenArrays unChosenStartPositions
                    |> fusePositionFrequencyMatrices motifLength
                    |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSsWithBPV motifLength alphabet sources.[n] pcv positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[n] then 
                        acc.[n] <- tmp
                        acc
                     else acc                    
                    )
                    bestMotif
        loop 0 (Array.copy startPositions) (Array.copy startPositions)

    /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
    /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
    let getPWMOfRandomStartsWithBPV (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) =    
        let rec loop (n:int) acc =
            if n = sources.Length then (List.rev acc) |> List.toArray
            else
                let unChosenArrays =
                    Array.append (sources.[0..n-1]) (sources.[n+1..])
                let randomStartPositions = 
                    unChosenArrays
                    |> Array.map (fun unChosen ->
                        getRandomNumberInSequence motifLength unChosen)
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        (getDefinedSegment motifLength subSequence position)
                         |> getPositionFrequencyMatrix) unChosenArrays randomStartPositions
                    |> fusePositionFrequencyMatrices motifLength
                    |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
                loop (n + 1) (getBestPWMSsWithBPV motifLength alphabet sources.[n] pcv positionProbabilityMatrix::acc)
        loop 0 []

    /// Repeats the search for the best InformationContent of each found PositionWeightMatrix until they converge or a maximum number of repetitions.
    /// Each iteration it is checked if better InformationContent is found or not.
    let getMotifsWithBestInformationContentWithBPV (numberOfRepetitions:int) (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) =
        let rec loop (n:int) (acc:(float*int)[]) (bestPWMS:(float*int)[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            getPWMOfRandomStartsWithBPV motifLength pseudoCount alphabet sources pcv
                            |> getBestPWMSsWithStartPositionsWithBPV motifLength pseudoCount alphabet sources pcv
                            |> getLeftShiftedBestPWMSsWithBPV motifLength pseudoCount alphabet sources pcv
                            |> getRightShiftedBestPWMSsWithBPV motifLength pseudoCount alphabet sources pcv
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|0., 0|]

    /// Gives the startPosition and score of the segment with the highest PositionWeightMatrixScore based on the given sequence and PositionWeightMatrix.
    let getBestPWMSs (motifLength:int) (alphabet:#IBioItem[]) (pseudoCount:float) (source:BioArray.BioArray<#IBioItem>) (bfVector:FrequencyCompositeVector) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        let rec loop (n:int) (highValue:float) (highIndex:int) =
            if n + motifLength > source.Length then log2(highValue), highIndex
            else
                let tmp =
                    let segment =
                        Array.skip n source
                        |> Array.take motifLength
                    let pcv =
                        addInPlaceCountsOfSource source bfVector
                        |> substractCountsOfBFVector segment 
                        |> getNormalizedProbabilityVector alphabet pseudoCount
                    let pwMatrix = createPositionWeightMatrix alphabet pcv positionProbabilityMatrix
                    segment
                    |> calculateSegmentScore pwMatrix
                if tmp > highValue then loop (n + 1) tmp n
                else loop (n + 1) highValue highIndex
        loop 0 0. 0

    /// Checks whether downstream of given positions a higher InformationContent is present or not. 
    /// If yes, the new InformationContent and positions are given back, otherwise the old ones.
    let getRightShiftedBestPWMSs (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (startPositions:(float*int)[]) =
        let rec loop (n:int) (acc:(float*int)[]) (bestMotif:(float*int)[]) =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> snd item)) = (bestMotif |> Array.map (fun item -> snd item)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    let tmp = Array.append sources.[0..n-1] sources.[n+1..]
                    Array.append bestMotif.[0..n-1] bestMotif.[n+1..]
                    |> Array.map2 (fun (source:BioArray.BioArray<#IBioItem>) (_, position) -> if position <= source.Length - motifLength - 1 then position + 1 else position) tmp
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        getFrequencyVectorWithoutSegment motifLength position unchosenArray) unChosenArrays unChosenStartPositions
                    |> fuseFrequencyVectors alphabet
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        getDefinedSegment motifLength subSequence position
                        |> getPositionFrequencyMatrix) unChosenArrays unChosenStartPositions
                    |> fusePositionFrequencyMatrices motifLength
                    |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSs motifLength alphabet pseudoCount sources.[n] frequencyCompositeVector positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[n] then 
                        acc.[n] <- tmp
                        acc
                     else acc                    
                    )
                    bestMotif
        loop 0 (Array.copy startPositions) (Array.copy startPositions)

    /// Checks whether upstream of given positions a higher PositionWeightMatrixScore is present or not. 
    /// If yes, the new PositionWeightMatrixScore and positions are given back, otherwise the old ones.
    let getLeftShiftedBestPWMSs (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (startPositions:(float*int)[]) =
        let rec loop (n:int) (acc:(float*int)[]) (bestMotif:(float*int)[]) =
            if n = sources.Length then
                if (acc |> Array.map (fun item -> snd item)) = (bestMotif |> Array.map (fun item -> snd item)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    Array.append bestMotif.[0..n-1] bestMotif.[n+1..]
                    |> Array.map (fun (_, position) -> if position > 0 then position - 1 else position)
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        getFrequencyVectorWithoutSegment motifLength position unchosenArray) 
                            unChosenArrays unChosenStartPositions
                    |> fuseFrequencyVectors alphabet
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        getDefinedSegment motifLength subSequence position
                        |> getPositionFrequencyMatrix)  unChosenArrays unChosenStartPositions
                        |> fusePositionFrequencyMatrices motifLength
                    |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSs motifLength alphabet pseudoCount sources.[n] frequencyCompositeVector positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[n] then 
                        acc.[n] <- tmp
                        acc
                     else acc                    
                    )
                    bestMotif
        loop 0 (Array.copy startPositions) (Array.copy startPositions)

    /// Checks the given Sequence for the existence of a conserved motif, by scoring each segment based on the given start positions.
    /// The new PositionWeightMatrix is calculated and updated each iteration if segments with higher scores are found until convergence.
    let getBestPWMSsWithStartPositions (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (startPositions:(float*int)[]) =        
        let rec loop (n:int) acc bestMotif =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> snd item)) = (bestMotif |> Array.map (fun item -> snd item)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    Array.append acc.[0..n-1] acc.[n+1..]
                    |> Array.map (fun (_, position) -> position)
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        getFrequencyVectorWithoutSegment motifLength position unchosenArray) 
                            unChosenArrays unChosenStartPositions
                    |> fuseFrequencyVectors alphabet
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        getDefinedSegment motifLength subSequence position
                        |> getPositionFrequencyMatrix) unChosenArrays unChosenStartPositions
                    |> fusePositionFrequencyMatrices motifLength
                    |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
                let tmp = getBestPWMSs motifLength alphabet pseudoCount sources.[n] frequencyCompositeVector positionProbabilityMatrix
                loop 
                    (n + 1) 
                    (if fst tmp > fst acc.[n] then 
                        acc.[n] <- tmp
                        acc
                     else acc                    
                    )
                    bestMotif
        loop 0 (Array.copy startPositions) (Array.copy startPositions)

    /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
    /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
    let getPWMOfRandomStarts (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        let rec loop (n:int) acc =
            if n = sources.Length then (List.rev acc) |> List.toArray
            else
                let unChosenArrays =
                    Array.append (sources.[0..n-1]) (sources.[n+1..])
                let randomStartPositions = 
                    unChosenArrays
                    |> Array.map (fun unChosen ->
                        getRandomNumberInSequence motifLength unChosen)
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        getFrequencyVectorWithoutSegment motifLength position unchosenArray) 
                            unChosenArrays randomStartPositions
                    |> fuseFrequencyVectors alphabet
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence position -> 
                        getDefinedSegment motifLength subSequence position
                        |> getPositionFrequencyMatrix) unChosenArrays randomStartPositions
                    |> fusePositionFrequencyMatrices motifLength
                    |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
                loop (n + 1) (getBestPWMSs motifLength alphabet pseudoCount sources.[n] frequencyCompositeVector positionProbabilityMatrix::acc)
        loop 0 []

    /// Repeats the search for the best InformationContent of each found PositionWeightMatrix until they converge or a maximum number of repetitions.
    /// Each iteration it is checked if better InformationContent is found or not.
    let getMotifsWithBestInformationContent (numberOfRepetitions:int) (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        let rec loop (n:int) (acc:(float*int)[]) (bestPWMS:(float*int)[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            getPWMOfRandomStarts motifLength pseudoCount alphabet sources
                            |> getBestPWMSsWithStartPositions motifLength pseudoCount alphabet sources
                            |> getLeftShiftedBestPWMSs motifLength pseudoCount alphabet sources
                            |> getRightShiftedBestPWMSs motifLength pseudoCount alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|0., 0|]

    /// Creates a random start position for each sequence and calculates a PositionWeightMatrix based on the. 
    /// The PositionWeithMatrix is then used to find the best PositionWeightMatrixScore for each sequence and gives you back the new Positions and PositionWeightMatrixScores.
    let getMotifsWithBestPWMSOfPPM (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        let rec loop (n:int) acc =
            if n = sources.Length then (List.rev acc) |> List.toArray
            else
                let unChosenArrays =
                    Array.append (sources.[0..n-1]) (sources.[n+1..])
                let randomStartPositions = 
                    unChosenArrays
                    |> Array.map (fun unChosen ->
                        getRandomNumberInSequence motifLength unChosen)
                let frequencyCompositeVector =
                    Array.map2 (fun unchosenArray position ->
                        getFrequencyVectorWithoutSegment motifLength position unchosenArray) 
                            unChosenArrays randomStartPositions
                    |> fuseFrequencyVectors alphabet
                loop (n + 1) (getBestPWMSs motifLength alphabet pseudoCount sources.[n] frequencyCompositeVector positionProbabilityMatrix::acc)
        loop 0 []

    /// Repeats the search for the best InformationContent of each found PositionWeightMatrix until they converge or a maximum number of repetitions.
    /// Each iteration it is checked if better InformationContent is found or not.
    let getBestInformationContentOfPPM (numberOfRepetitions:int) (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        let rec loop (n:int) (acc:(float*int)[]) (bestPWMS:(float*int)[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun (pwms, _) -> pwms)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            getMotifsWithBestPWMSOfPPM motifLength pseudoCount alphabet sources positionProbabilityMatrix
                            |> getBestPWMSsWithStartPositions motifLength pseudoCount alphabet sources
                            |> getLeftShiftedBestPWMSs motifLength pseudoCount alphabet sources
                            |> getRightShiftedBestPWMSs motifLength pseudoCount alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|0., 0|]

    let doSiteSamplingWithBPV (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) =
        getPWMOfRandomStartsWithBPV motifLength pseudoCount alphabet sources pcv
        |> getBestPWMSsWithStartPositionsWithBPV motifLength pseudoCount alphabet sources pcv
        |> getLeftShiftedBestPWMSsWithBPV motifLength pseudoCount alphabet sources pcv
        |> getRightShiftedBestPWMSsWithBPV motifLength pseudoCount alphabet sources pcv

    let doSiteSampling (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        getPWMOfRandomStarts motifLength pseudoCount alphabet sources
        |> getBestPWMSsWithStartPositions motifLength pseudoCount alphabet sources
        |> getLeftShiftedBestPWMSs motifLength pseudoCount alphabet sources
        |> getRightShiftedBestPWMSs motifLength pseudoCount alphabet sources

    let doSiteSamplingWithPPM (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (ppM:PositionProbabilityMatrix) =
        getMotifsWithBestPWMSOfPPM motifLength pseudoCount alphabet sources ppM
        |> getBestPWMSsWithStartPositions motifLength pseudoCount alphabet sources
        |> getLeftShiftedBestPWMSs motifLength pseudoCount alphabet sources
        |> getRightShiftedBestPWMSs motifLength pseudoCount alphabet sources

module MotifSampler =

    /// Contains the information of the probability and the start position(s) of the found motifs(s). 
    type MotifIndex =
        {
            PWMS        : float
            Positions   : int list
        }

    /// Creates a MotifIndex based on the probability of the given segement(s) start position(s).
    let createMotifIndex pwms pos =
        {
            MotifIndex.PWMS         = pwms
            MotifIndex.Positions    = pos
        }

    /// Calculates all possible segment combination of m segments, who do not overlap in the given width 
    /// and gives you the probability for it.
    let calculateMultiCombWithWidth (cutOff:float) (width:int) (m:int) (set:list<float*int>) =
        let rec loop prob positions size set = 
            seq
                {
                    match size, set with
                    | n, x::xs ->
                        if n > 0 then
                            if ceckForDistance width (snd x::positions) then
                                if log2(fst x*prob) > cutOff then
                                    yield! loop (fst x*prob) (snd x::positions) (n - 1) xs
                        if n >= 0 then yield! loop prob positions n xs
                    | 0, [] -> yield createMotifIndex (log2(prob)) positions
                    | _, [] -> () 
                }
        loop 1. [] m set
        |> List.ofSeq

    /// Normalizes the probabilities of all MotifMemories to the sum of all probabilities and picks one by random, 
    /// but those with a higher value have a higher chance to get picked.
    let rouletteWheelSelection (pick:float) (items:MotifIndex list) =
        let normalizedItems =
            let sum = List.sum (items |> List.map (fun item -> item.PWMS))
            items
            |> List.map (fun item -> item.PWMS/sum, item)
        let rec loop acc n =
            if acc <= pick && pick <= acc + fst normalizedItems.[n] then items.[n]
            else loop (acc + fst normalizedItems.[n]) (n + 1)
        loop 0. 0       

    /// Calculates the normalized segment score based on the given PositionWeightMatrix and 
    /// BackGroundProbabilityVecor of all segment combinations of a given sequence. 
    /// The amount of combinations is given by the amount of motifs.
    let calculateNormalizedSegmentScores (cutOff:float) (motifAmount:int) (motifLength:int) (source:BioArray.BioArray<#IBioItem>) (pcv:ProbabilityCompositeVector) (pwMatrix:PositionWeightMatrix) =
        let segments =
            let rec loop n acc =
                if n + motifLength = source.Length+1 then List.rev acc
                else
                    let tmp =
                        source
                        |> Array.skip n
                        |> (fun items -> Array.take motifLength items, n)
                    loop (n+1) (tmp::acc)
            loop 0 []
        let segmentScores =
            segments
            |> List.map (fun segment -> 
                calculateSegmentScore pwMatrix (fst segment), (snd segment))
        let backGroundScores =
            segments
            |> List.map (fun segment -> calculateBackGroundSegmentScore pcv (fst segment))
            |> List.map (fun segmentScore -> createMotifIndex segmentScore [])
        let tmp =
            let rec loop n acc =
                if n > motifAmount then List.rev acc
                else loop (n+1) (calculateMultiCombWithWidth cutOff motifLength n segmentScores::acc)
            loop 1 [backGroundScores]
        tmp
        |> List.concat
        //|> List.collect (fun items -> 
        //    items 
        //    |> List.map (fun item -> createMotifIndex (item.PWMS/items.Head.PWMS) item.Positions))

    //let getRightShiftedBestMotifMotifSampler (motifAmount:int) (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (motifMem:MotifIndex[]) =
    //    let rec loop (n:int) (acc:MotifIndex []) (bestMotif:MotifIndex []) =
    //        if n = sources.Length then
    //            if (acc |> Array.map (fun item -> item.Positions)) = (bestMotif |> Array.map (fun item -> item.Positions)) then acc
    //            else loop 0 acc (Array.copy acc)
    //        else
    //            let unChosenStartPositions =
    //                Array.append bestMotif.[0..n-1] bestMotif.[n+1..]
    //                |> Array.map (fun item -> 
    //                    (item.Positions |> Array.ofSeq)
    //                    |> Array.map(fun position -> if position < sources.[n].Length - motifLength then position - 1 else position))
    //            let unChosenArrays =
    //                Array.append sources.[0..n-1] sources.[n+1..]
    //            let backgroundProbabilityVector =
    //                Array.map2 (fun array positions -> 
    //                    positions
    //                    |> Array.map (fun position -> 
    //                        getFrequencyVectorWithoutSegment position motifLength array)
    //                            ) unChosenArrays unChosenStartPositions
    //                |> Array.concat
    //                |> fuseFrequencyVectors alphabet
    //                |> addInPlaceCountsOfSource sources.[n]
    //                |> getNormalizedProbabilityVector alphabet pseudoCount
    //            let positionProbabilityMatrix =
    //                Array.map2 (fun subSequence positions -> 
    //                    positions
    //                    |> Array.map (fun position -> getDefinedSegment motifLength subSequence position)
    //                           ) unChosenArrays unChosenStartPositions
    //                |> Array.concat
    //                |> Array.map getPositionFrequencyMatrix
    //                |> fusePositionFrequencyMatrices motifLength
    //                |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
    //            let positionWeightMatrix = createPositionWeightMatrix alphabet backgroundProbabilityVector positionProbabilityMatrix
    //            let tmp = 
    //                calculateNormalizedSegmentScores motifAmount motifLength sources.[n] backgroundProbabilityVector positionWeightMatrix
    //                |> List.sortByDescending (fun item -> item.PWMS)
    //                |> List.head
    //            loop 
    //                (n + 1) 
    //                (if tmp.PWMS > acc.[n].PWMS then 
    //                    acc.[n] <- tmp
    //                    acc
    //                 else acc                    
    //                )
    //                bestMotif
    //    loop 0 (Array.copy motifMem) (Array.copy motifMem)

    //let getLeftShiftedBestMotifMotifSampler (motifAmount:int) (motifLength:int) (pseudoCount:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (motifMem:MotifIndex[]) =
    //    let rec loop (n:int) (acc:MotifIndex []) (bestMotif:MotifIndex []) =
    //        if n = sources.Length then
    //            if (acc |> Array.map (fun item -> item.Positions)) = (bestMotif |> Array.map (fun item -> item.Positions)) then acc
    //            else loop 0 acc (Array.copy acc)
    //        else
    //            let unChosenStartPositions =
    //                Array.append bestMotif.[0..n-1] bestMotif.[n+1..]
    //                |> Array.map (fun item -> 
    //                    (item.Positions |> Array.ofSeq)
    //                    |> Array.map(fun position -> if position > 0 then position - 1 else position))
    //            let unChosenArrays =
    //                Array.append sources.[0..n-1] sources.[n+1..]
    //            let backgroundProbabilityVector =
    //                Array.map2 (fun array positions -> 
    //                    positions
    //                    |> Array.map (fun position -> 
    //                        getFrequencyVectorWithoutSegment motifLength position array)
    //                            ) unChosenArrays unChosenStartPositions
    //                |> Array.concat
    //                |> fuseFrequencyVectors alphabet
    //                |> addInPlaceCountsOfSource sources.[n]
    //                |> getNormalizedProbabilityVector alphabet pseudoCount
    //            let positionProbabilityMatrix =
    //                Array.map2 (fun subSequence positions -> 
    //                    positions
    //                    |> Array.map (fun position -> getDefinedSegment motifLength subSequence position)
    //                           ) unChosenArrays unChosenStartPositions
    //                |> Array.concat
    //                |> Array.map getPositionFrequencyMatrix
    //                |> fusePositionFrequencyMatrices motifLength
    //                |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
    //            let positionWeightMatrix = createPositionWeightMatrix alphabet backgroundProbabilityVector positionProbabilityMatrix
    //            let tmp = 
    //                calculateNormalizedSegmentScores motifAmount motifLength sources.[n] backgroundProbabilityVector positionWeightMatrix
    //                |> List.sortByDescending (fun item -> item.PWMS)
    //                |> List.head
    //            loop 
    //                (n + 1) 
    //                (if tmp.PWMS > acc.[n].PWMS then 
    //                    acc.[n] <- tmp
    //                    acc
    //                 else acc                    
    //                )
    //                bestMotif
    //    loop 0 (Array.copy motifMem) (Array.copy motifMem)

    /// Checks the given Sequence for the existence of conserved motif(s), by scoring each segment and segment combination based on the given start positions.
    /// The new PositionWeightMatrix is calculated and updated each iteration if segments with higher scores are found until convergence. 
    let getBestPWMSsWithStartPositionsWithBPV (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (motifMem:MotifIndex[]) =
        let rec loop (n:int) acc bestMotif =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> item.Positions)) = (bestMotif |> Array.map (fun item -> item.Positions)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    Array.append acc.[0..n-1] acc.[n+1..]
                    |> Array.map (fun item -> item.Positions |> List.toArray)
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let backgroundProbabilityVector =
                    Array.map2 (fun array positions -> 
                        positions
                        |> Array.map (fun position -> 
                            getFrequencyVectorWithoutSegment motifLength position array)
                                ) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fuseFrequencyVectors alphabet
                    |> addInPlaceCountsOfSource sources.[n]
                    |> getNormalizedProbabilityVector alphabet pseudoCount
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence positions -> 
                        positions
                        |> Array.map (fun position -> 
                            getDefinedSegment motifLength subSequence position
                            |> getPositionFrequencyMatrix)) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fusePositionFrequencyMatrices motifLength
                    |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
                let positionWeightMatrix = createPositionWeightMatrix alphabet backgroundProbabilityVector positionProbabilityMatrix
                let tmp = 
                    calculateNormalizedSegmentScores cutOff motifAmount motifLength sources.[n] backgroundProbabilityVector positionWeightMatrix
                    |> List.sortByDescending (fun item -> item.PWMS)
                    |> List.head
                loop 
                    (n + 1) 
                    (if tmp.PWMS > acc.[n].PWMS then 
                        acc.[n] <- tmp
                        acc
                     else acc                    
                    )
                    bestMotif
        loop 0 (Array.copy motifMem) (Array.copy motifMem)


    /// Checks the given Sequence for the existence of conserved motif(s), by scoring each segment and segment combination based on the given start positions.
    /// The nex segments for the nex PositionWeightMatrix are picked by chance after calculating the segment scores for each possible combination 
    /// but those with higher scores have a higher chance to be picked.
    let getBestPWMSsByHighChanceMotifSamplingWithBPV (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) (motifMem:MotifIndex[]) =
        let rnd = new System.Random()
        let rec loop (n:int) acc =
            if n = sources.Length then (List.rev acc) |> Array.ofList
            else
                let unChosenStartPositions =
                    Array.append motifMem.[0..n-1] motifMem.[n+1..]
                    |> Array.map (fun item -> item.Positions |> List.toArray)
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let backgroundProbabilityVector =
                    Array.map2 (fun array positions -> 
                        positions
                        |> Array.map (fun position ->
                            getFrequencyVectorWithoutSegment motifLength position array)
                                ) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fuseFrequencyVectors alphabet
                    |> addInPlaceCountsOfSource sources.[n]
                    |> getNormalizedProbabilityVector alphabet pseudoCount
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence positions -> 
                        positions
                        |> Array.map (fun position -> 
                            getDefinedSegment motifLength subSequence position
                            |> getPositionFrequencyMatrix)) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fusePositionFrequencyMatrices motifLength
                    |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
                let positionWeightMatrix = createPositionWeightMatrix alphabet backgroundProbabilityVector positionProbabilityMatrix
                let tmp = 
                    calculateNormalizedSegmentScores cutOff motifAmount motifLength sources.[n] backgroundProbabilityVector positionWeightMatrix
                    |> rouletteWheelSelection (rnd.NextDouble())
                loop (n+1) (tmp::acc)
        loop 0 []

    ///Repeats the motif sampler until convergence or until a given number of repetitions is done.
    let getBestPWMSsWithBPV (numberOfRepetitions:int) (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (pcv:ProbabilityCompositeVector) =
        let rec loop (n:int) (acc:MotifIndex[]) (bestPWMS:MotifIndex[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            SiteSampler.getPWMOfRandomStartsWithBPV motifLength pseudoCount alphabet sources pcv
                            |> Array.map (fun (prob, position) -> createMotifIndex prob [position])
                            |> getBestPWMSsByHighChanceMotifSamplingWithBPV motifAmount motifLength pseudoCount cutOff alphabet sources pcv
                            |> getBestPWMSsWithStartPositionsWithBPV motifAmount motifLength pseudoCount cutOff alphabet sources pcv
                            //|> getLeftShiftedBestMotifMotifSampler motifAmount motifLength pseudoCount alphabet sources
                            //|> getRightShiftedBestMotifMotifSampler motifAmount motifLength pseudoCount alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|createMotifIndex 0. []|]

    /// Checks the given Sequence for the existence of conserved motif(s), by scoring each segment and segment combination based on the given start positions.
    /// The new PositionWeightMatrix is calculated and updated each iteration if segments with higher scores are found until convergence. 
    let getBestPWMSsWithStartPositions (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (motifMem:MotifIndex[]) =
        let rec loop (n:int) acc bestMotif =
            if n = sources.Length then 
                if (acc |> Array.map (fun item -> item.Positions)) = (bestMotif |> Array.map (fun item -> item.Positions)) then acc
                else loop 0 acc (Array.copy acc)
            else
                let unChosenStartPositions =
                    Array.append acc.[0..n-1] acc.[n+1..]
                    //|> Array.map (fun item -> item.Positions |> List.toArray)
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let backgroundProbabilityVector =
                    Array.map2 (fun array positions -> 
                        positions.Positions
                        |> List.map (fun position -> 
                            getFrequencyVectorWithoutSegment motifLength position array)
                            |> List.toArray
                                ) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fuseFrequencyVectors alphabet
                    |> addInPlaceCountsOfSource sources.[n]
                    |> getNormalizedProbabilityVector alphabet pseudoCount
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence positions -> 
                        positions.Positions
                        |> List.map (fun position -> 
                            getDefinedSegment motifLength subSequence position
                            |> getPositionFrequencyMatrix)
                            |> List.toArray) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fusePositionFrequencyMatrices motifLength
                    |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
                let positionWeightMatrix = createPositionWeightMatrix alphabet backgroundProbabilityVector positionProbabilityMatrix
                let tmp = 
                    calculateNormalizedSegmentScores cutOff motifAmount motifLength sources.[n] backgroundProbabilityVector positionWeightMatrix
                    |> List.sortByDescending (fun item -> item.PWMS)
                    |> List.head
                loop 
                    (n + 1) 
                    (if tmp.PWMS > acc.[n].PWMS then 
                        acc.[n] <- tmp
                        acc
                     else acc                    
                    )
                    bestMotif
        loop 0 (Array.copy motifMem) (Array.copy motifMem)


    /// Checks the given Sequence for the existence of conserved motif(s), by scoring each segment and segment combination based on the given start positions.
    /// The nex segments for the nex PositionWeightMatrix are picked by chance after calculating the segment scores for each possible combination 
    /// but those with higher scores have a higher chance to be picked.
    let getBestPWMSsByHighChanceMotifSampling (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (motifMem:MotifIndex[]) =
        let rnd = new System.Random()
        let rec loop (n:int) acc =
            if n = sources.Length then (List.rev acc) |> Array.ofList
            else
                let unChosenStartPositions =
                    Array.append motifMem.[0..n-1] motifMem.[n+1..]
                    //|> Array.map (fun item -> item.Positions)
                let unChosenArrays =
                    Array.append sources.[0..n-1] sources.[n+1..]
                let backgroundProbabilityVector =
                    Array.map2 (fun array positions -> 
                        positions.Positions
                        |> List.map (fun position ->
                            getFrequencyVectorWithoutSegment motifLength position array)
                            |> List.toArray
                                ) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fuseFrequencyVectors alphabet
                    |> addInPlaceCountsOfSource sources.[n]
                    |> getNormalizedProbabilityVector alphabet pseudoCount
                let positionProbabilityMatrix =
                    Array.map2 (fun subSequence positions -> 
                        positions.Positions
                        |> List.map (fun position -> 
                            getDefinedSegment motifLength subSequence position
                            |> getPositionFrequencyMatrix)
                            |> List.toArray) unChosenArrays unChosenStartPositions
                    |> Array.concat
                    |> fusePositionFrequencyMatrices motifLength
                    |> getPositionProbabilityMatrix (sources.Length - 1) alphabet pseudoCount
                let positionWeightMatrix = createPositionWeightMatrix alphabet backgroundProbabilityVector positionProbabilityMatrix
                let tmp = 
                    calculateNormalizedSegmentScores cutOff motifAmount motifLength sources.[n] backgroundProbabilityVector positionWeightMatrix
                    |> rouletteWheelSelection (rnd.NextDouble())
                loop (n+1) (tmp::acc)
        loop 0 []

    ///Repeats the motif sampler until convergence or until a given number of repetitions is done.
    let getBestInformationContent (numberOfRepetitions:int) (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        let rec loop (n:int) (acc:MotifIndex[]) (bestPWMS:MotifIndex[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            SiteSampler.getPWMOfRandomStarts motifLength pseudoCount alphabet sources
                            |> Array.map (fun (prob, position) -> createMotifIndex prob [position])
                            |> getBestPWMSsByHighChanceMotifSampling motifAmount motifLength pseudoCount cutOff alphabet sources
                            |> getBestPWMSsWithStartPositions motifAmount motifLength pseudoCount cutOff alphabet sources
                            //|> getLeftShiftedBestMotifMotifSampler motifAmount motifLength pseudoCount alphabet sources
                            //|> getRightShiftedBestMotifMotifSampler motifAmount motifLength pseudoCount alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|createMotifIndex 0. []|]

    ///Repeats the motif sampler until convergence or until a given number of repetitions is done.
    let getBestPWMSsOfPPM (numberOfRepetitions:int) (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        let rec loop (n:int) (acc:MotifIndex[]) (bestPWMS:MotifIndex[]) =
            if n > numberOfRepetitions then
                bestPWMS
            else
                if acc = bestPWMS then 
                    bestPWMS
                else 
                    let informationContentAcc =
                        acc
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    let informationContentBestPWMS =
                        bestPWMS
                        |> Array.map (fun item -> item.PWMS)
                        |> Array.sum
                    if informationContentAcc > informationContentBestPWMS then
                        loop (n + 1) [||] (if Array.isEmpty acc then bestPWMS else acc)
                    else
                        let pwms =
                            SiteSampler.getMotifsWithBestPWMSOfPPM motifLength pseudoCount alphabet sources positionProbabilityMatrix
                            |> Array.map (fun (prob, position) -> createMotifIndex prob [position])
                            |> getBestPWMSsByHighChanceMotifSampling motifAmount motifLength pseudoCount cutOff alphabet sources
                            |> getBestPWMSsWithStartPositions motifAmount motifLength pseudoCount cutOff alphabet sources
                            //|> getLeftShiftedBestMotifMotifSampler motifAmount motifLength pseudoCount alphabet sources
                            //|> getRightShiftedBestMotifMotifSampler motifAmount motifLength pseudoCount alphabet sources
                        loop (n + 1) pwms bestPWMS
        loop 0 [||] [|createMotifIndex 0. []|]

    let doMotifSamplingWithPPM (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) (positionProbabilityMatrix:PositionProbabilityMatrix) =
        SiteSampler.getMotifsWithBestPWMSOfPPM motifLength pseudoCount alphabet sources positionProbabilityMatrix
        |> Array.map (fun (prob, position) -> createMotifIndex prob [position])
        |> getBestPWMSsByHighChanceMotifSampling motifAmount motifLength pseudoCount cutOff alphabet sources
        |> getBestPWMSsWithStartPositions motifAmount motifLength pseudoCount cutOff alphabet sources

    let doMotifSampling (motifAmount:int) (motifLength:int) (pseudoCount:float) (cutOff:float) (alphabet:#IBioItem[]) (sources:BioArray.BioArray<#IBioItem>[]) =
        SiteSampler.getPWMOfRandomStarts motifLength pseudoCount alphabet sources
        |> Array.map (fun (prob, position) -> createMotifIndex prob [position])
        |> getBestPWMSsByHighChanceMotifSampling motifAmount motifLength pseudoCount cutOff alphabet sources
        |> getBestPWMSsWithStartPositions motifAmount motifLength pseudoCount cutOff alphabet sources

open SiteSampler
open MotifSampler

#time

let fileDir1 = "C:/Users/Student/source/repos/CsbScaffold-VP_CSB/GibbsAlgorithm/material/ChlamySeq3.tsv"

let tests =
    [|
        "GTGGCTGCACCACGTGTATGC"
        "ACATCGCATCACGTGACCAGT"
        "CCTCGCACGTGGTGGTACAGT"
        "CTCGTTAGGACCATCACGTGA"
    |]

let profileSequenceForTests =
    [
        "CACGTG"
        "CACGTG"
        "CACGTG"
        "CACGTG"
    ]

let bioTests =
    tests
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let bioTestsWithMultipleSamples = 
    [|
        "GTGGCTGCACCACGTGTATGCCACGTG"
        "ACATCGCATCACGTGACCAGTTAGTTG"
        "CCTCGCACGTGGTGGTACAGTCGTACG"
        "GCATAAAGGACCATCACGTGAAGCTGC"
        "TTTTTTTTTTTTTTTTTTTTTTTTTTT"
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let bioTestsII =
    [|
        "GTAAGTACAGAAAGCCACAGAGTACCATCTAGGAAATTAACATTATACTAACTTTCTACATCGTTGATACTTATGCGTATACATTCATATA"
        "AGACAGAGTCTAAAGATTGCATTACAAGAAAAAAGTTCTCATTACTAACAAGCAAAATGTTTTGTTTCTCCTTTTA"
        "GTATGTTCATGTCTCATTCTCCTTTTCGGCTCCGTTTAGGTGATAAACGTACTATATTGTGAAAGATTATTTACTAACGACACATTGAAG*"
        "GCATGTGTGCTGCCCAAGTTGAGAAGAGATACTAACAAAATGACCGCGGCTCTCAAAAATAATTGACGAGCTTACGGTGATACGCTTACCG"
        "GTATGTTTGACGAGAATTGCTAGTGTGCGGGAAACTTTGCTACCTTTTTTGGTGCGATGCAACAGGTTACTAATATGTAATACTTCAG"
        "TTTCAAGATTAACCACATCTGCTAACTTTCTCCCTATGCTTTTACTAACAAAATTATTCTCACTCCCCGATATTGA"
        "GTAAGTATCCAGATTTTACTTCATATATTTGCCTTTTTCTGTGCTCCGACTTACTAACATTGTATTCTCCCCTTCTTCATTTTAG"
        "GTATGCATAGGCAATAACTTCGGCCTCATACTCAAAGAACACGTTTACTAACATAACTTATTTACATAG"
        "GTATGTAGTAGGGAAATATATCAAAGGAACAAAATGAAAGCTATGTGATTCCGTAATTTACGAAGGCAAATTACTAACATTGAAATACGGG"
        "GTATGTTACTATTTGGAGTTTCATGAGGCTTTTCCCGCCGTAGATCGAACCCAATCTTACTAACAGAGAAAGGGCTTTTTCCCGACCATCA"
        "TATGTAATGATATATTATGAAGTAAGTTCCCCAAAGCCAATTAACTAACCGAATTTTAATCTGCACTCATCATTAG"
        "GTATGTTCATAATGATTTACATCGGAATTCCCTTTGATACAAGAAAACTAACGGGTATCGTACATCAATTTTTGAAAAAAGTCAAGTACTA"
        "GTATGTATATTTTTGACTTTTTGAGTCTCAACTACCGAAGAGAAATAAACTACTAACGTACTTTAATATTTATAG"
        "TTTCGACGCGAATAGACTTTTTCCTTCTTACAGAACGATAATAACTAACATGACTTTAACAG"
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let profileSequenceForTestsII =
    [|"TACTAAC"; "TACTAAT"; "AACTAAC"|]

let hseConsensus =
    [|
        "-----cGTCcaGAAggCGCcaTACg----------"   //HSP70A Promoter
        "gGGAagCTCtgGAAggGCCgcGATg----------"   //HSP70A Promoter
        "tGAAgcTACagGACt--------------------"   //HSP70A Promoter
        "cGAAggGCCgcGACggTTCgaGAAccGACttGAGg"   //HSP70A Promoter
        "cGAAggTTCgaTGAa--------------------"   //HSP70G Promoter
        "gGAAagAACtcAAAgaTTTc---------------"   //HSP70G Promoter
        "aAAAcgTGAtgGAAg--------------------"   //HSP70G Promoter    
        "gGACttATCgcGTGgaTTCtgGACt----------"   //CLPB3 Promoter
        "-----gTTCggGTAcgTTCtgGCAc----------"   //CLPB3 Promoter
        "-----tTGCttGAAgtTCTc---------------"   //CLPB3 Promoter
        "cGAAggTTCtcGAGgtTTCc---------------"   //CLPB3 Promoter
        "-----gTTCcgGTAttGTCc---------------"   //CLPB3 Promoter
        "cACAcaTTCaaGTAc--------------------"   //HSF1 5'UTR
        "-----tTTCccGCAacTTTg---------------"   //HSF1 5'UTR
        "-----tTTCtaGAAgcTTGccAATgcTTTt-----"   //HSF1 5'UTR
        "tGTGgtTTCgaGAAa--------------------"   //HSF1 Promoter
        "gGAAaaTTCgcGAGacTTGc---------------"   //CDJ1 5'UTR
        "tGGAggTTCcgCGAg--------------------"   //CDJ1 5'UTR
        "tGAAtaTGCatGCGtgTTAt---------------"   //CDJ1 Promoter
        "-----gGTGttGAAacTCCtgCTAg----------"   //CDJ1 Promoter
        "aGAAgcTGCgtGGCa--------------------"   //CDJ1 Promoter
        "-----cCTAgaGAGcaTTCc---------------"   //CDJ1 Promoter
        "-----aTGCgtGTAgaTCGttGCAaaTTCa-----"   //LHCBM9 between 5'UTR beginning and ending
        "-----gTTGggGAAcaATTttGAAt----------"   //LHCBM9 between 5'UTR beginning and ending
        "cGAAggTCTtgGCAg--------------------"   //LHCBM9 between 5'UTR beginning and ending
        "cGCActTTCgcCCAa--------------------"   //TIC110 Promoter
        "tTCAcgTACaaGAAg--------------------"   //TIC110 Promoter
        "aGACgcTTCtaGAGaaAGAacTTCagAAAg-----"   //TIC110 Promoter
        "-----cGGCaaGAAagTTCa---------------"   //TIC110 Promoter
        "cTAAaaAGCcaGGAatTTCccAAAc----------"   //TIC40 Promoter
        "-----gTTCcaGAAgtTTGt---------------"   //TIC40 Promoter
        "-----gTTCctGGAaaTAAtgTAAa----------"   //TIC40 Promoter
        "-----cTTGtaGAGacTTCatGTTc----------"   //TIC40 Promoter
        "cGCAgcTTGcaGAAc--------------------"   //FFC
        "-----aTGCcaAAAtgTTCaaCCAccTCCc-----"   //FFC
        "tACGgtTTCaaGAAa--------------------"   //FFC
        "gGAAgtTTTgtGGGc--------------------"   //FFC
        "cGAAcaTGActGAAttGTCc---------------"   //MPA1 5'UTR
        "-----gTTCttGATccAACggGAGg----------"   //MPA1 5'UTR
        "-----cGACcgGAAtgTTCc---------------"   //MPA1 5'UTR
        "-----tTTCgaGTAagTTCg---------------"   //MPA1 5'UTR

        //"-GAA--TTC--GAA--TTC--GAA--TTC-" // HSE-Consensus
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let trimmedHSEConsensus =
    [|
        "cGTCcaGAAggCGCc"   //HSP70A Promoter
        "gCTCtgGAAggGCCg"   //HSP70A Promoter
        "cTACagGACt-----"   //HSP70A Promoter
        "gGCCgcGACggTTCg"   //HSP70A Promoter
        "gTTCgaTGAa-----"   //HSP70G Promoter
        "gAACtcAAAgaTTTc"   //HSP70G Promoter
        "gTGAtgGAAg-----"   //HSP70G Promoter    
        "tATCgcGTGgaTTCt"   //CLPB3 Promoter
        "gTTCggGTAcgTTCt"   //CLPB3 Promoter
        "tTGCttGAAgtTCTc"   //CLPB3 Promoter
        "gTTCtcGAGgtTTCc"   //CLPB3 Promoter
        "gTTCcgGTAttGTCc"   //CLPB3 Promoter
        "aTTCaaGTAc-----"   //HSF1 5'UTR
        "tTTCccGCAacTTTg"   //HSF1 5'UTR
        "tTTCtaGAAgcTTGc"   //HSF1 5'UTR
        "tTTCgaGAAa-----"   //HSF1 Promoter
        "aTTCgcGAGacTTGc"   //CDJ1 5'UTR
        "gTTCcgCGAg-----"   //CDJ1 5'UTR
        "aTGCatGCGtgTTAt"   //CDJ1 Promoter
        "gGTGttGAAacTCCt"   //CDJ1 Promoter
        "cTGCgtGGCa-----"   //CDJ1 Promoter
        "cCTAgaGAGcaTTCc"   //CDJ1 Promoter
        "aTGCgtGTAgaTCGt"   //LHCBM9 between 5'UTR beginning and ending
        "gTTGggGAAcaATTt"   //LHCBM9 between 5'UTR beginning and ending
        "gTCTtgGCAg-----"   //LHCBM9 between 5'UTR beginning and ending
        "tTTCgcCCAa-----"   //TIC110 Promoter
        "gTACaaGAAg-----"   //TIC110 Promoter
        "cTTCtaGAGaaAGAa"   //TIC110 Promoter
        "cGGCaaGAAagTTCa"   //TIC110 Promoter
        "aAGCcaGGAatTTCc"   //TIC40 Promoter
        "gTTCcaGAAgtTTGt"   //TIC40 Promoter
        "gTTCctGGAaaTAAt"   //TIC40 Promoter
        "cTTGtaGAGacTTCa"   //TIC40 Promoter
        "cTTGcaGAAc-----"   //FFC
        "aTGCcaAAAtgTTCa"   //FFC
        "tTTCaaGAAa-----"   //FFC
        "tTTTgtGGGc-----"   //FFC
        "aTGActGAAttGTCc"   //MPA1 5'UTR
        "gTTCttGATccAACg"   //MPA1 5'UTR
        "cGACcgGAAtgTTCc"   //MPA1 5'UTR
        "tTTCgaGTAagTTCg"   //MPA1 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let completeFilledUpHSEConsensus =
    [|
        "CGCGGCGTCCAGAAGGCGCCATACGGCCCGCTGGC"   //HSP70A Promoter
        "gGGAagCTCtgGAAggGCCgcGATgGGGCGCGCGG"   //HSP70A Promoter
        "TGAAGCTACAGGACTGATTTGGCGGGCTATGAGGG"   //HSP70A Promoter
        "cGAAggGCCgcGACggTTCgaGAAccGACttGAGg"   //HSP70A Promoter
        "CGAAGGTTCGATGAACAGGACAACCTGTCCCGTTG"   //HSP70G Promoter
        "GGAAAGAACTCAAAGATTTCGCGGTGACGGCGTGA"   //HSP70G Promoter   //Motif found with HSP70G ending on PromoterSeq
        "AAAACGTGATGGAAGGCCTGCCGTGCTCTGCTGAG"   //HSP70G Promoter    
        "GGACTTATCGCGTGGATTCTGGACTCCCGTGAGGG"   //CLPB3 Promoter
        "AGCGAGTTCGGGTACGTTCTGGCACGCGTGGCCTG"   //CLPB3 Promoter
        "GCTCCTTGCTTGAAGTTCTCGCTCTGTACGACAGA"   //CLPB3 Promoter
        "CGAAGGTTCTCGAGGTTTCCGTGGGAGTGTGAATG"   //CLPB3 Promoter    //Motif found with HSP70G ending on PromoterSeq
        "TACACGTTCCGGTATTGTCCTTATTCGAAGGTTCT"   //CLPB3 Promoter
        "CACACATTCAAGTACTTAACGACTCACTGCAAAGA"   //HSF1 5'UTR
        "CCAGTTTTCCCGCAACTTTGCTTACACTTACTTTA"   //HSF1 5'UTR
        "CGTTCTTTCTAGAAGCTTGCCAATGCTTTTAATGG"   //HSF1 5'UTR        //Motif found with HSP70G ending on PromoterSeq
        "TGTGGTTTCGAGAAATGTAAGCTAGTGTGAATGCA"   //HSF1 Promoter     //Motif found with HSP70G ending on PromoterSeq
        "GGAAAATTCGCGAGACTTGCAGCTCACAACTTCGC"   //CDJ1 5'UTR        //Motif found with HSP70G ending on PromoterSeq
        "TGGAGGTTCCGCGAGGCCATGGGCGGAAAATTCGC"   //CDJ1 5'UTR
        "TGAATATGCATGCGTGTTATGCTTTCAATGCAGCC"   //CDJ1 Promoter
        "TGGCCGGTGTTGAAACTCCTGCTAGCGATGCACCG"   //CDJ1 Promoter     //Motif found with HSP70G ending on PromoterSeq
        "AGAAGCTGCGTGGCATGTTGCTGGCCGGTGTTGAA"   //CDJ1 Promoter
        "AAGTACCTAGAGAGCATTCCCAAGTTGAGTCGGCC"   //CDJ1 Promoter
        "TGTGCATGCGTGTAGATCGTTGCAAATTCAGTGCG"   //LHCBM9 between 5'UTR beginning and ending
        "ACCTAGTTGGGGAACAATTTTGAATGGAAAAGTGA"   //LHCBM9 between 5'UTR beginning and ending
        "CGAAGGTCTTGGCAGTGCGCAACCTTGACGTGGTT"   //LHCBM9 between 5'UTR beginning and ending
        "CGCACTTTCGCCCAAACGTCGACGTGGTCTTTCTT"   //TIC110 Promoter
        "TTCACGTACAAGAAGTGCAAATTCAATTGAAACCT"   //TIC110 Promoter
        "AGACGCTTCTAGAGAAAGAACTTCAGAAAGTCCAG"   //TIC110 Promoter
        "CCGGCCGGCAAGAAAGTTCAGTAGCGTACCTGCGA"   //TIC110 Promoter   //Motif found with HSP70G ending on PromoterSeq
        "CTAAAAAGCCAGGAATTTCCCAAACTCAGAGCTGC"   //TIC40 Promoter
        "AGGGAGTTCCAGAAGTTTGTGTCCAAACCGGTCGC"   //TIC40 Promoter    //Motif found with HSP70G ending on PromoterSeq
        "GTGCAGTTCCTGGAAATAATGTAAACAAGCCCTGG"   //TIC40 Promoter
        "GTATGCTTGTAGAGACTTCATGTTCCGGGCCTTGC"   //TIC40 Promoter
        //"cGCAgcTTGcaGAAc--------------------"   //FFC
        //"-----aTGCcaAAAtgTTCaaCCAccTCCc-----"   //FFC
        //"tACGgtTTCaaGAAa--------------------"   //FFC
        //"gGAAgtTTTgtGGGc--------------------"   //FFC
        "CGAACATGACTGAATTGTCCAATCACAAACAACGA"   //MPA1 5'UTR
        "CCCTGGTTCTTGATCCAACGGGAGGCGCGCTTGAC"   //MPA1 5'UTR
        "TTAACCGACCGGAATGTTCCCTGGTTCTTGATCCA"   //MPA1 5'UTR
        "GCCGCTTTCGAGTAAGTTCGGAGCACAAAACATCA"   //MPA1 5'UTR        //Motif found with HSP70G ending on PromoterSeq

        //"-GAA--TTC--GAA--TTC--GAA--TTC-" // HSE-Consensus
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let hsp70aGene =
    [|
        "GCAGCTCCGCAGGGCACACGTCACGCGAAGGGCCGCGACGGTTCGAGAACCGACTTGAGGGCGCCAAACGAGCCCGAGCCGCCGTTGCGCCA
         GGCGAAACCAGAACCGTAGATTAATGCACTTGAGCTATTCATTGGAGCGATCTGCCGGGGACAGCGGGTCTGGCGTGCGCGCGATTGGAGAT
         CGCAAATTACATATGTCTGCGTGACGGCGGGGAGCTCGCTGAGGCTTGACATGATTGGTGCGTATGTTTGTATGAAGCTACAGGACTGATTT
         GGCGGGCTATGAGGGCGGGGGAAGCTCTGGAAGGGCCGCGATGGGGCGCGCGGCGTCCAGAAGGCGCCATACGGCCCGCTGGCGGCACCCAT
         CCGGTATAAAAGCCCGCGACCCCGAACGGTGA" //400 bases upstream 5' UTR
        //"CCTCCACTTTCAGCGACAAACGAGCACTTATACATACGCGACTATTCTGCCGCTATACATAACCACTCAACTCGCTTAAGAGTCAGTAAAC" //5'UTR
    |]

let hsp70gGene =
    [|
        "GCACTGTAAAACGTGATGGAAGGCCTGCCGTGCTCTGCTGAGTTGCTTGAGGCCATGGGGGACCGCCTCAGGGAACAGGGCCAACTCATGG
         GCCAGGCAGACTGGCATGGCTTCATGTGCAAGCCAGTGCGGTTCCAAGTCAGACACAGCTACGCCCAGGCCCGTGTAGCAAGCATTGCGCG
         GCAGGACGGAAAGAACTCAAAGATTTCGCGGTGACGGCGTGAAGGGGTGCACCCGGGAGACTTGTGAGCGCGAGAGACGTTGTTCGGCACG
         CCACGGCCCTGTGGTGCCAGTGCATACCTGGGCATACGTGTGGGGGGGGGGGGGGGGAAGGGGCTGCCTCCCTGGCTTCGGCCCCGTTAAC
         TCAGCTCGCGTCGTTTGTGCAGGTGTGCAGTTCACCGTCAGTTACCTTGTAACCTAGTGTCATCCTGTCCCCGGTGGGAAAAGGAGGTCGC
         GTGCCCAGGCATGCACCGTTGAATCAACCGAAGGTTCGATGAACAGGACAACCT" // 509 bases upstream 5'UTR
        //"GTCCCGTTGCCACCCTGTCGCCTTACTCTCTGTCCCAAATTCGCTCCTTTTAATACAGCGTCAAATTGTTGCAAATAGGCGTTAAACCAA" //5'UTR
    |]

let clpb3Gene =
    [|
        "CATACACGTTCCGGTATTGTCCTTATTCGAAGGTTCTCGAGGTTTCCGTGGGAGTGTGAATGCCTTGGTCCGCACTTCCTTCCGGCTTGCC
         AGGGCCTAGGGAACTAGTGTGCTGCCGTGCTCGTGTGTGCAGGAGCGCGGAGCGTTAATTGAGTTGGTGCGGAGCATGCATCGCTGCACTG
         GTTGAATCACACCCGAATATGAGCACAGTTCCACGGTGCATGGTGTGCCGTGAAAAGGGATGTCGTGAGGCATGAGGAATGGCCCGGCTCC
         TTGCTTGAAGTTCTCGCTCTGTACGACAGACGTTGCGCTTTGCTGACGGGGGAGGGCGTTCAGGAGGTCATGCCATGACGAGAGCCTGGCC
         CGGCCCACGGAACATGGTCGCAGGTGTTGGGATGCTGGACCCTCGAGATAAGATCCTTCCGCGCCAACGCGCCGTCGAACCCCTCCAACCC
         CCCCGAAAGCATGCAAGCTACCTGTTGCAACCAATGTCTTTGCGTCACATTAGCGAGTTCGGGTACGTTCTGGCACGCGTGGCCTGGACTT
         ATCGCGTGGATTCTGGACTCCCGTGAGGGCCACTCGCAAAACCTTTCAACCCCGGCGCGGGTATATATAAAGATATGTATTCCTCTACGCC
         ATCCAATTAATTA" // 650 bases upstream of the 5'UTR
        //"CGCACACCGCAAACGCTTATGCTTCAGACGCTTCAAAACAGCCTAACAACTTTTGAGCCGCTTATCGCTTCCCAGGCCCTGCCAGCAGCAGCCGAACTCCTTCGCCAAATC" // 5'UTR
    |]

let hsf1Gene =
    [|
        "TGCTCACAGATGAACCGCAGCTCCACGTGCCTCCCGGTCGCCGCCATTCTCTTCGCGGCACGCGCGAGCGCCTTCCGGTCGGCGCTGACAG
         CTGCCGACGGCCGATGGCGCGCGGTTGGTCGGGCCGACAGGCGACAGGGTGGGGTTGGGATGGGGCTTGGCGGGGCCGGAAGGGTTGGAGT
         TGAAATGCACCTGGCCTAGTCTGCCTCCGTATGACCCCCATCACACCCCCGTAAACAGGTCCGGCCACACCTCACCCCAGCCCTCGCCCAT
         GCCTTGAGCATGTGCATGGATTGTGGTTTCGAGAAATGTAAGCTAGTGTGAATGCAAATGCAAGACGCCATCGTGGTTTGACGTGCCGGTA
         GTAGGCCGCACAGCCAGTGTGTCGGCCTAAGGCGTGGGGGCCCCTCCCTCGCGCTGCCGAGCAGTACCTGGTGTTGTCCACGTGCCCACGC
         CATGTTGAGCGGATAGGTTAACCCTAACTGCACGAGCGTGCAAGGACGTGCGGGGGCTAGGCATCGAACCATGCTGCTGAGGCTGGTAGCG
         GTAGGGTATGCGTGTTACTCAAATCGAGCCCAGGCTGGTACTCGTGCAGCATGTGGCGTGCTAAGCCCATAAGTCTAAGCGTGGCAAACGG
         GAGACTGGTCGGCGACATGTTTCGCAATTCCGGTTGCGATTCTGGTTGGCACGTCCCGAGCCC" // 700 bases upstream of the 5'UTR
        "CCACTGCACGGCACTGCGGTGCCTATCTGTATGTGCAGACTACACCTCGCGCGGAGCTAGCTCGCGCTGGCTCTGTGGCGTTCTTTCTAGA
         AGCTTGCCAATGCTTTTAATGGGCACTCGCTCAGAAGCCCGGCTTTTCCAGTTTTCCCGCAACTTTGCTTACACTTACTTTAATCTTTTGT
         TACCTGACACACATTCAAGTACTTAACGACTCACTGCAAAGAATGCGTTGTTTGAGCAAGCGCACGCTCACCGGCTCGCCGGCTTGTTAAT
         CGTCGTCCTTTTAGTCTGCTTAAATCACCTGTTAGTGACATCTTACTTGAGGCCAGGGTCAAAGCCTTTTAGCTTCTTGCATTTTGGGCTC
         GGAGGGGCCTGGCCGGTCGCTGCAACGGAGCAACAAGGAGCT" // 5'UTR
    |]

let cdj1Gene =
    [|
        "AGACTTCATATGCATACTAAATGTATAGCAAAATACGTGTACATGCACTATGCCTTTAGAGCTGTAGTTTTGGAAGTGGAGTCAGGGCACT
         CGGACGGGTTTGCTACCATGTATGCACAACCAAGTACCTAGAGAGCATTCCCAAGTTGAGTCGGCCAGTTCGGCGAGCGCCGACGATGCTG
         ACACGCCACGAACCCGGCATTGATGCCTTTCTTTCCCCCCCTCTCTCTCCCACCCATTCTTCTTCCCGTCCCTAGGCGCTGCAGAAGCTGC
         GTGGCATGTTGCTGGCCGGTGTTGAAACTCCTGCTAGCGATGCACCGCTGCGCCCGGCCCGGCCCGGCACCACGAGTATCCTGCCCGACCT
         GGCCCCCTTGCTTAATGCCCCTGCTGCAGCCGACAGCCTGTTATGCACAACATTTGCCACACCGAACAGCAGCTCATCTTCGCTACACGCC
         GCGCAGGGGCAGGGGAGCTGCAGGACCGCGACGAACCCTCATTTGCGGTCCTTGGGTCGCGGCCAGCTTGGAACCTAGCCCGGCAGCGAGC
         ATGCTCCACCGGCTGAATATGCATGCGTGTTATGCTTTCAATGCAGCCACGGAA" // 600 bases upstream of the 5'UTR
        "CCCATTAGCTGCGCAGAGGTGGAGGTTCCGCGAGGCCATGGGCGGAAAATTCGCGAGACTTGCAGCTCACAACTTCGCGTGGTGTCTTGCT
         GCATACTGTTTTACTTGTGTCAGTCATAATTATATCGCATATCGTTTATTGGGCTGACCAACGTGACCCAGTAGGCGCAAGGGCGCAATCA
         AGTTTGTAGTTTCTTGGTTCACCAGCGAA" // 5'UTR
    |]

let lhcbm9Gene =
    [|
        //"CTGTACTTTTGGGCGAAATCGGTCTCATCCGGGCCCCCGGCGAGCACTGCACAGCGAATGGTGCAACGAATTGTACAGCCCCGCTTTCGCC
        // CCCTCTAGACCTTCCACAGAGTGTTTAGTCCCCGATTATGCCATGGTTATCCCGGGGCGGGAGCGGCGGCGTGTTTGCGGTGGGAGAGGGG
        // GCGAGGGGGGGAGCCGGGGGGAGTGCAGACAGGTTCGCGAGGTCGCCGAGGCCCCGCGCCAACTTGCGCTGCGTAGACCTCGATTAGCAGA
        // TGCTATTTTAGCCTTTTCCAAGCACTTACAAGCTCTGTATTTTCCGTGTACTTTTCGCGTGTGCAATTCGGCGCGAAGCCGAAAATCATGC
        // GATTGACCCTGTACCGTGTACTTTTGCCCCGCGCGATGGCTCCCAGAACCGGCCATTCACTTTGCAAGTCCCTGTAAGCACATACATGCAG
        // GTGAAATACCACAAACAAGCTACCTTGAGGTGGGGGTAGTCCGAA" // 500 bases upstream of the 5'UTR
        //"TCCTAGCGAAATGCCAACCCGCCAAACATGTCGGCAAGCCTGCGCACGCGCCAG"    // 5'UTR beginning
        "GTGCGCGCGTTTTGCCCGGCCGGACTGCTTCAAGTACGCCTGTACGCACAAGTGGATGAAGCAACCGGTACATGCCATTTGTGCAGTCGCC
         ACGCAGGCGCTAGGCACGTCGGCACGTCCCACCTCTCCTGGCAGAGCCGATCGCCCGGCCACCGCGTGGCCAGAGCACGTCGATTTCCGGC
         ACCCAGAGTGCCGCTAATGCATGCTGACTAGTGGGCATCATCACAGCGGGGCGCGTAGGCGCGGCTCACAGGCGCCTGCGCACTGGACCGG
         CGTGCACCACCCCCCAGCAGGTGCATTCCATGACGCCCCACGGAACCTGCTGCACGCCCCCGCACCCGCCACGTCCACGTGGCCCCCAACC
         ACGCCACTACGCACCTTACCCGCTACTGTCAGAACTCAAAAGTATCAGGGGCCCGGGTCAATTTGTTCCCCATTTCTTTCGTCCGCGTCTC
         TGCCGTGACTTGATTGCAGACCATGGCCCCGTTGGCATCCATACATGCCCCGCGCGCGGCAGCAGTTGGCAGCTGACGCTGTGACGCAGCC
         ATTCGATTCCTAGCAGATACGTACCGTAAGGATGTTGGTCGACTGGGAGATAGCGGACGAGCGCTGCCAGCAGATCCTGTGTAGCCTAGTG
         CCCCCGGCCCCTGCGACCAGCAGTGAAACAACCCAGGTTGACCGCCCCCCTATGCGTTCGCACCAAACCGTCTGAGCACATGGTCTCAAGC
         AAGTGCACGGACAAGTGTTCCCCGGTCACGGGAAAAGGGGCTGAGCCCTTCCACAGTCTGAGGCCGCACAACCTCATCGCCCGGACATAGC
         CGGGGTCGGTTTCCACAGGCACATAGCCAAACCAGGTGCCTAGCAACCGCGATGCCCAGCAGAGCCAACCGCAGCCAGACGCAGCACAGTT
         GTCAACAGCCGCACTGCTTGTCATTGCCGCTACGCATGCCCTGCCAGCCCCACCGGAAAAAGGTGCACGGAAAAAGGTATTGTAGGGTAAT
         GGCACAGGGAAAACAACCAGAACCAAGAGGGCACCCCTGCCACAACGCACCGCAAAGCCCAAGCGGTTGGCGGGAGGGGCCGAAGGTCTTG
         GCAGTGCGCAACCTTGACGTGGTTGATGTTACCGTACCTAGTTGGGGAACAATTTTGAATGGAAAAGTGACAGTGGTCCATAGGCCATAGT
         AATGTAGTCGGGAAGCAGTGGGGTTGAGACCTGGCATGCCGCATGCGGCACCGTACGGTACGTAAGCAGGCGCGGGTGTGTGCATGCGTGT
         AGATCGTTGCAAATTCAGTGCGTGCGGCCGGAGGGAAGGGAAGTTAGCGGGATTGGGAGCTAACAATAATGCTAAAGGCGAAATTAAATGG
         GAAGGGGGCCGCGGGTGGAAAAGTGGGGACGGTAGGAGGGGGCCGGAGCTGGCCAGGGCGCAGGAGACGAAGGGGTCTGGCGACCTCGGCC
         TCGGCACCGTGGCCTCTCGCCCTTCATCTCTGGCCTGCTGTCCTTGAAGCCGGCAGCAAGAGATCTTGAGATTCCAGTGTATTTAAGCTGA
         ATAG" // between 5'UTR beginning and ending
        //"GTAGCGGACAAGTTATCTCGCATCTCTGGTCCACTCACAACTCCCACTAGCAAAG"   // 5'UTR ending
    |]

let tic110Gene =
    [|
        "AGCAGCAACATCTGCCGGGCAGCAGGGTTTGCCGCCGCCGGCACGGCACGCCGAAGCACGACGCCAAACCCGAGCTGTCCGGCCGGCAAGA
         AAGTTCAGTAGCGTACCTGCGAAGCGCTTGTCCACCACCACTTTGTTCTTCTGCTTGGGAAAGCGCTGGAAGCGCGGGTCGCTGTGCATGG
         CGGCGAAACGGGGGTCGACCTTGCCCCCGCCCTTGCCTTCCGGTTTACCCATCTTGATGTGTTTATTTGCAATAGATTACAATGTCAAAAA
         CTGTTATTTTGTTCAAGAAGACGCTTCTAGAGAAAGAACTTCAGAAAGTCCAGGTACTTCTCTACGGGGCTCAGTAATTCACGTACAAGAA
         GTGCAAATTCAATTGAAACCTTTGTGATAACAAGAAACACAATTGACCCTGAGACCGGACCCATGCACGATCGCTGATGCGGAGGACCGCA
         CTTTCGCCCAAACGTCGACGTGGTCTTTCTTTTCTAGGGCAGCTA" // 500 bases upstream of the 5'UTR
        //"GGAACCCAGGACTAACTGGGCATCGCCCTTGACTTGACTGCTTTACTTGCACACACGCCGTCAAGCATTTAGCTTAGGTTAGCAGCTCCGTCGCAGCC" //5'UTR
    |]

let tic40Gene =
    [|
        "CACGCCGCGACCGCGCCCATGACCGCCGCAGCCCCCGGATGACTCGTATGCTTGTAGAGACTTCATGTTCCGGGCCTTGCATCTCGGATTG
         ACCTACCGTCCGCGTGGACGACCGGGGGCTGCTCCTCCTGCGTCTCCGCCATGCTCCGGGCTCCTTACAGCGCAGAGCGCAGCTAAATGTG
         TGCCACTGAAATGTTTGACCGAACTTTGTGTCCCGGCTCAGCAGACGGTGCACGACGACGTCAAGTCCGGGAGACCATAGCCCGCCTGTGA
         GGCTCCGCGTCTGTGCAGTTCCTGGAAATAATGTAAACAAGCCCTGGCGTGCAGGTGAAGTACGTCGTTCAGGGAGTTCCAGAAGTTTGTG
         TCCAAACCGGTCGCTTTGTTTTTTTGCTAAGCGTTACAAATGAAACTAAGCGTGCAATTGAATCAGCCTGGAAAGCGGCGCCTAAAAAGCC
         AGGAATTTCCCAAACTCAGAGCTGCATGCAGGCTGAGTGCTTTGTCCCAGCCCCCCACAGGAATTCCCTGGCACTAATTTCTTCTCCCCTC
         TCTA" // 550 bases upstream of the 5'UTR
        //"GTTAGTTGATACTTAGAGGTCAGCTACAATCGTAAGAATGTTGTGTTCCCGCACAACCGCCCCTCGGGGATTCGGGAAGGCAGCCCTTGCT
        // CTGCCCTCCCGCGCTTCCATCAAATGCCGGGCTGTGCAAG" // 5'UTR
    |]

let ffcGene =
    [|
        "Nothing found"
    |]

let mpa1Gene =
    [|
        //"GCGTAGAAACAATTTAAGATCGTAATTTGCCCATCACCAGCACTGTATGCGCCCGTATGCCCGCCTGTGTGTCACCAGCGATTAGGTAA
        // TTTGCCCCTGCATGCCCAAACCTGACGATGGGGGTTTACGCGCGCGCGGGTGAATACCTTTGCGCGAACCCGTTGTCGGCTTCTGAGTA
        // GTATTATTAGCTGCATAACCGC" //200 bases upstream of the 5'UTR
        "GATTCGACGCAGTGGCCAGCAGCGCAAAACCAGGGCATCGCGGCCGCTTTCGAGTAAGTTCGGAGCACAAAACATCAACCGCTTTAACCGA
         CCGGAATGTTCCCTGGTTCTTGATCCAACGGGAGGCGCGCTTGACAATGCATATAAATACGGGCGAACATGACTGAATTGTCCAATCACAA
         ACAACGACTCACCTACTCCGCCCAACACTTAAGCTGTTGGTTTCAGTTGAC" // 5'UTR
    |]
   
let geneCollection = 
    [|hsp70aGene; hsp70gGene; clpb3Gene; hsf1Gene; cdj1Gene; lhcbm9Gene; tic110Gene; tic40Gene; mpa1Gene|]
    |> Array.concat
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let geneCollectionReduced = 
    [|hsp70aGene; hsp70gGene; clpb3Gene; hsf1Gene; cdj1Gene; lhcbm9Gene; tic110Gene; tic40Gene; mpa1Gene|]
    |> Array.concat
    |> Array.map (fun test -> BioArray.ofNucleotideString test)
 
/// Contains a Set of all DNA bases.
let dnaBases = 
    [|Nucleotides.Nucleotide.A; Nucleotides.Nucleotide.T; Nucleotides.Nucleotide.G; Nucleotides.Nucleotide.C; Nucleotides.Nucleotide.Gap|]

/// Contains a Set of all amino acids.
let aminoAcids  =   
    [|
        AminoAcidSymbols.AminoAcidSymbol.Ala; AminoAcidSymbols.AminoAcidSymbol.Arg; AminoAcidSymbols.AminoAcidSymbol.Asn; 
        AminoAcidSymbols.AminoAcidSymbol.Asp; AminoAcidSymbols.AminoAcidSymbol.Asx; AminoAcidSymbols.AminoAcidSymbol.Cys;
        AminoAcidSymbols.AminoAcidSymbol.Xle; AminoAcidSymbols.AminoAcidSymbol.Gln; AminoAcidSymbols.AminoAcidSymbol.Glu;
        AminoAcidSymbols.AminoAcidSymbol.Glx; AminoAcidSymbols.AminoAcidSymbol.Gly; AminoAcidSymbols.AminoAcidSymbol.His;
        AminoAcidSymbols.AminoAcidSymbol.Ile; AminoAcidSymbols.AminoAcidSymbol.Leu; AminoAcidSymbols.AminoAcidSymbol.Lys;
        AminoAcidSymbols.AminoAcidSymbol.Met; AminoAcidSymbols.AminoAcidSymbol.Phe; AminoAcidSymbols.AminoAcidSymbol.Pro;
        AminoAcidSymbols.AminoAcidSymbol.Pyl; AminoAcidSymbols.AminoAcidSymbol.Sel; AminoAcidSymbols.AminoAcidSymbol.Ser;
        AminoAcidSymbols.AminoAcidSymbol.Thr; AminoAcidSymbols.AminoAcidSymbol.Trp; AminoAcidSymbols.AminoAcidSymbol.Val;
    |]

//let testI = Array.init 1 (fun _ -> getMotifsWithBestInformationContent 1 6 0.0001 dnaBases bioTests)

//testI
//|> Array.countBy (fun items -> items |> Array.map (fun item -> snd item))
//|> Array.sortByDescending (fun (_, i) -> i)

//let testII = 
//    Array.init 1 (fun _ -> getMotifsWithBestInformationContent 1 7 0.0001 dnaBases bioTestsII)
//    |> List.ofSeq

//testII
//|> List.map (fun (items) -> items |> Array.map (fun (_, y) -> y))
//|> groupEquals
//|> List.sortByDescending (fun (_, i) -> i)

//for i=0 to testII.[0].Length-1 do
//    printfn "%A" (getDefinedSegment 20 geneCollection.[i] (snd testII.[0].[i]))

//let realTest =
//    fromFileObo fileDir1
//    |> Array.map (fun test -> BioArray.ofAminoAcidSymbolString test)
//    |> Array.filter (fun item -> item.Length <> 0)

//let tmp = Array.init 100 (fun _ -> getBestInformationContent 100 2 6 0.0001 1. dnaBases bioTestsWithMultipleSamples)

//tmp
//|> Array.countBy (fun items -> items |> Array.map (fun item -> item.Positions))
//|> Array.sortByDescending (fun (item, amount) -> amount)

//let dataTest = Array.init 1 (fun _ -> getMotifsWithBestInformationContent 1 10 0.0001 dnaBases geneCollectionReduced)

//dataTest
//|> List.ofArray
//|> groupEquals
//|> List.sortByDescending (fun (item, amount) -> Array.sum (item |> Array.map (fun (prob, _) -> prob)))

//dataTest
//|> List.ofArray
//|> groupEquals
//|> List.sortByDescending (fun (item, amount) -> amount)

//for i=0 to dataTest.[0].Length-1 do
//    printfn "%A" (getDefinedSegment 20 geneCollection.[i] (snd dataTest.[0].[i]))

//let testTMP =
//    Array.init 1 (fun _ ->
//        PSeq.map (fun _ -> doSiteSampling 10 0.0001 dnaBases geneCollection) [1..1]
//        |> Array.ofSeq
//        |> HelperFunctionsAndTypes.getBestInformationContent
//                 )

//testTMP
//|> List.ofArray
//|> List.map (fun (items) -> items |> Array.map (fun (_, y) -> y))
//|> groupEquals
//|> List.sortByDescending (fun (item, amount) -> amount)

//getBestPWMS testTMP

let hseConsensus1stHalf =
    [|
        "cGTCcaGAAggCGCc"   //HSP70A Promoter
        "gCTCtgGAAggGCCg"   //HSP70A Promoter
        "cTACagGACt-----"   //HSP70A Promoter
        "gGCCgcGACggTTCg"   //HSP70A Promoter
        "gTTCgaTGAa-----"   //HSP70G Promoter
        "gAACtcAAAgaTTTc"   //HSP70G Promoter
        "gTGAtgGAAg-----"   //HSP70G Promoter    
        "tATCgcGTGgaTTCt"   //CLPB3 Promoter
        "gTTCggGTAcgTTCt"   //CLPB3 Promoter
        "tTGCttGAAgtTCTc"   //CLPB3 Promoter
        "gTTCtcGAGgtTTCc"   //CLPB3 Promoter
        "gTTCcgGTAttGTCc"   //CLPB3 Promoter
        "aTTCaaGTAc-----"   //HSF1 5'UTR
        "tTTCccGCAacTTTg"   //HSF1 5'UTR
        "tTTCtaGAAgcTTGc"   //HSF1 5'UTR
        "tTTCgaGAAa-----"   //HSF1 Promoter
        "aTTCgcGAGacTTGc"   //CDJ1 5'UTR
        "gTTCcgCGAg-----"   //CDJ1 5'UTR
        "aTGCatGCGtgTTAt"   //CDJ1 Promoter
        "gGTGttGAAacTCCt"   //CDJ1 Promoter
        "cTGCgtGGCa-----"   //CDJ1 Promoter
        "cCTAgaGAGcaTTCc"   //CDJ1 Promoter
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let hseConsensus2ndHalf =
    [|
        "aTGCgtGTAg"   //LHCBM9 between 5'UTR beginning and ending
        "gTTGggGAAc"   //LHCBM9 between 5'UTR beginning and ending
        "gTCTtgGCAg"   //LHCBM9 between 5'UTR beginning and ending
        "tTTCgcCCAa"   //TIC110 Promoter
        "gTACaaGAAg"   //TIC110 Promoter
        "cTTCtaGAGa"   //TIC110 Promoter
        "cGGCaaGAAa"   //TIC110 Promoter
        "aAGCcaGGAa"   //TIC40 Promoter
        "gTTCcaGAAg"   //TIC40 Promoter
        "gTTCctGGAa"   //TIC40 Promoter
        "cTTGtaGAGa"   //TIC40 Promoter
        "cTTGcaGAAc"   //FFC
        "aTGCcaAAAt"   //FFC
        "tTTCaaGAAa"   //FFC
        "tTTTgtGGGc"   //FFC
        "aTGActGAAt"   //MPA1 5'UTR
        "gTTCttGATc"   //MPA1 5'UTR
        "cGACcgGAAt"   //MPA1 5'UTR
        "tTTCgaGTAa"   //MPA1 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let geneCollection1stHalf = 
    [|hsp70aGene; hsp70gGene; clpb3Gene; hsf1Gene; cdj1Gene|]
    |> Array.concat
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let geneCollection2ndHalf = 
    [|lhcbm9Gene; tic110Gene; tic40Gene; mpa1Gene|]
    |> Array.concat
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let ppm =
    hseConsensus2ndHalf
    |> Array.map getPositionFrequencyMatrix
    |> fusePositionFrequencyMatrices hseConsensus2ndHalf.[0].Length
    |> getPositionProbabilityMatrix (hseConsensus2ndHalf.Length) dnaBases (0.0001)

let test = 
    doSiteSamplingWithPPM hseConsensus2ndHalf.[0].Length 0.0001 dnaBases geneCollection1stHalf ppm
    //getBestInformationContentOfPPM 1 hseConsensus2ndHalf.[0].Length 0.0001 dnaBases geneCollection1stHalf ppm

test
|> Array.countBy (fun item -> item)
|> Array.sortByDescending (fun (item, amount) -> amount)

for i=0 to test.Length-1 do
    printfn "%A" (getDefinedSegment 30 geneCollection1stHalf.[i] (snd test.[i]))

//Found motifs size of 15 of ppmOfFirstHalf in genes of 2ndHalf
"ATGCGTGTAGATCGT"
"GTTCAAGAAGACGCT"
"GTTCCAGAAGTTTGT"
"TTTCGAGTAAGTTCG"

//Found motifs size of 10 of ppmOfSecondtHalf in genes of 1fstHalf
"GGTTCGAGAA"
"TTGCTTGAGG"
"TTGCTTGAAG"
"GTTTCGAGAA"
"TTTCTAGAAG"
"TTTTGGAAGT"
"ATTCGCGAGA"

let testPPM = 
    profileSequenceForTests
    |> Array.ofList
    |> Array.map (fun test -> BioArray.ofNucleotideString test)
    |> Array.map getPositionFrequencyMatrix
    |> fusePositionFrequencyMatrices hseConsensus2ndHalf.[0].Length
    |> getPositionProbabilityMatrix (hseConsensus2ndHalf.Length) dnaBases (0.0001)

let testMotifSampling = Array.init 1 (fun _ -> getBestPWMSsOfPPM 1 2 6 0.0001 10. dnaBases bioTestsWithMultipleSamples testPPM)

let au5g5407_t1__Cre14g617400t11 =
    [|
        "GGCGTGGGTGCAAGTCACCGACTCGCGAAGCAGCTCAGACGCCAGACGGCCCATGCCCAGTGAGAAGGGGCTG
         GAGCGGTACACGCGCGCGGGCTTGCAGGCGTAGGGGCTGAGGTAGAAGGCCTGAGGCGTCACGGCGCCCTGCC
         GCGAAGGTCGGACCGGCAGCATCGACCGCATGGCGATGGCGTGAGCGCCGGCTGGCGCGGCGGTCTTGCGGCG
         GGCAGAGGCCGCGCCAGAGGTGGAAGCGACTGACTTCATGAGAGTGGTGGCCATGGCGCTGGCTTAAGCGATA
         GGAAGGAGGCGAGTGTGCTTGAATGTCGCGAGCTCGGCGAACGTCGAATGCGCAAAACCGGTTGGCTTGAGCA
         AGGTGTGCGCTAGGTACTAAAGCGTGTTGGCGTGCGTCGTTTGTATGTATGAACATCAATGGAGGGCGTGGAA
         CTTTTATAGCGGCGATGAGGCGCCTCACGAGAGTCGGGGCCAACGTCTGGAGCGGCCTCGAGCAGCTGGCCCC
         TCACACGTTCTGGAAGATTTGGTTTTGGGCTTTTGGCCTCCTCATCTTCGTCCATTGTATATATCTCAAGCAC
         AACCTTTTATTGGCTT" // 600 bases upstream of the 5'UTR
        "ACTCTAGCCGAGGCCTCGGGAATCCGGGAAACCGAACTTCTAGAACGTGTGAGGGGCCAGCTGCCTGGAGCCG
         CTCCAGACGTTGGCCCCGACTCTCGTAAGGCGCTCCATTGCCAGTATATAAGTCCACGTCCACCTTGCTTGCC
         AGAACATATAGGTAACATACGCAAACACGCTTCTATAGCTTGCGCACACCTTGCTCACGCCAACCGGTTTTGC
         GCATTCGACGTTCGCTCGCCGAACTCGCGACATTCAACCTCGCCTCCTTCCCATCACTTAAGCCAGCGCC" // 5'UTR
     |]
     |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g11124_t1__Cre03g199150t12 =
    [|
        "AGTGAACACGCGTTAATATGCATCAAGTAATCTGCTGCGGAAAATAAGGAGGACAATATGGAGACTCAGAATC
         GCAATTTGAGAGAAATGGCATGATGTACACATACGTTTATTTGGAACGTAATATGCGCAAGCCACCGAAGGCC
         AGTCAGCTGTCAGGGCATGGATCTGGATTCGGGAGTCGGCAGAGCACCGTGGGCCCGTGGCCCAAGCCGCCAA
         ACTTCAACCCGCCAAGGTTCCAGGATGCTATCAGAATGAGACAGAGCGCAATGGACCATGGCGTCAGCGAATA
         CGCAAGCTGGATGCGAGTCAGGTGCGGGGTTGAGGATTGAAGCCTTGGCACGCAGGGGCTGTTGGCACGCCTT
         GGCATGCCCAGGGCACTGTACAGGGCGGACCCGGGGGCACAGATGACCAGCCTTGAATCCGCGCACGACGAAG
         AGCAGTATCGTACATAGATAATCTCATCACAGAGGCGCTGAAATGTGCAGGGACATAGCGATGAGAGAAGCCC
         GCGGGGAGCGGGTGGGGCCGGGCAGGGGGCAAGCGTGTCGGATGGCCATTCACGCGGAGTCGCTCGCATACTC
         GGCTTTCCATAAGCGC" // 600 bases upstream of the 5'UTR
        "ATCGTTTGATTCGCAGCGCACTTAAAACCAATCTGCATGCCAGGCAACCCTACGCTTGCCGGAGCCGCGACTG
         TTGCGCCTACTTCCGTAGTTGTTCTGCGCAAACAGGATTCGCTTCCTTGATTTGATTCACTCGATTGATTGAG
         GCTATTTCTGTTGGCATAAC" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g15219_t1__Cre09g387150t12 =
    [|
        "GCCGATGGCGCGCGGTTGGTCGGGCCGACAGGCGACAGGGTGGGGTTGGGATGGGGCTTGGCGGGGCCGGAAG
         GGTTGGAGTTGAAATGCACCTGGCCTAGTCTGCCTCCGTATGACCCCCATCACACCCCCGTAAACAGGTCCGG
         CCACACCTCACCCCAGCCCTCGCCCATGCCTTGAGCATGTGCATGGATTGTGGTTTCGAGAAATGTAAGCTAG
         TGTGAATGCAAATGCAAGACGCCATCGTGGTTTGACGTGCCGGTAGTAGGCCGCACAGCCAGTGTGTCGGCCT
         AAGGCGTGGGGGCCCCTCCCTCGCGCTGCCGAGCAGTACCTGGTGTTGTCCACGTGCCCACGCCATGTTGAGC
         GGATAGGTTAACCCTAACTGCACGAGCGTGCAAGGACGTGCGGGGGCTAGGCATCGAACCATGCTGCTGAGGC
         TGGTAGCGGTAGGGTATGCGTGTTACTCAAATCGAGCCCAGGCTGGTACTCGTGCAGCATGTGGCGTGCTAAG
         CCCATAAGTCTAAGCGTGGCAAACGGGAGACTGGTCGGCGACATGTTTCGCAATTCCGGTTGCGATTCTGGTT
         GGCACGTCCCGAGCCC" // 600 bases upstream of the 5'UTR
        "CCACTGCACGGCACTGCGGTGCCTATCTGTATGTGCAGACTACACCTCGCGCGGAGCTAGCTCGCGCTGGCTC
         TGTGGCGTTCTTTCTAGAAGCTTGCCAATGCTTTTAATGGGCACTCGCTCAGAAGCCCGGCTTTTCCAGTTTT
         CCCGCAACTTTGCTTACACTTACTTTAATCTTTTGTTACCTGACACACATTCAAGTACTTAACGACTCACTGC
         AAAGAATGCGTTGTTTGAGCAAGCGCACGCTCACCGGCTCGCCGGCTTGTTAATCGTCGTCCTTTTAGTCTGC
         TTAAATCACCTGTTAGTGACATCTTACTTGAGGCCAGGGTCAAAGCCTTTTAGCTTCTTGCATTTTGGGCTCG
         GAGGGGCCTGGCCGGTCGCTGCAACGGAGCAACAAGGAGCT" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g15518_t1__Cre09g402304 =
    [|
        "CTGGGGACAGCGGGCGCCGCAGCGGAAGCGGCTGGGTAGGCGGCGGGGGCCAGCGCTGGTGTGTGCCATAGCC
         AACGTGAACCACGAAGGGGTCGCAGAAGGTGCGCGCGTGTGTGTGCACGGGGGCGGGCGTCCCAGGCGATGTG
         CTAAAGAAGGTGTGACGGCACGTAAGCGCCCATTGGCGCACGTGGCCTGCACAGCGTGGCGGCGCCCCTTCTG
         GCCAAGGGCGAAACCCTATGGTGCCTTCGGCTCTCAGCCCCCGCTACAGTAGTAGCTCTGCAGGGTGCCGTGC
         CTCACACCTTCATCCTCAACCCATTGCAATCACCCACACGCGAACACACTACAGGTTCACGACTGCTCCTTCA
         CCAACAAGTGCACCAATTACTGCCCGCTCACCAACAACAAACAGGGCTGCTGGGCCAGCAAGGACGCCACGCC
         CTGCCAGGTGGATCGCTCGGTGCTTGGGTTTGCTTGCCACGGTGCGCGTGCATCTGTGCACTACGCGCCAGGC
         CTCAGCCATGCATGCCGCGCCTCCACGTACCGGCATGCACCGCTGCCAGCGAGCGACGCTTGCTTGCTTGCTG
         CTCACGAATGAAGACC" // 600 bases upstream of the 5'UTR
        "GCGATTTGCCTACCTATCCTTCACCATGCTGGCTAGCTACTGCTC" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g9850_t1__Cre01g071662t11 =
    [|
        "TTTTAGGGAGAGACAGTCCATGATTTTGCATGTAAAATTCACGGCTGCCGCTGCTGCCACCGTGACGGCTACA
         GCGGCCGCCCCGGCCCCAGTCGATTCCCAGCAGTGGCGGTGGCGGCGGTGGCAGGCACGGACAGGGGGAATCC
         TGGCGGGCTTCGCCTCCTGTAAAGAAGCACATCCGGTCATGCACAGTCACAACGCCCAGAAAATATTGCCTTA
         AGGCAAATTCCTACAGTGTAGATGGAAAATTGAAGAATACGAGCACTGAATACCATCATTACAACAACGGTTG
         GAAATAACTACACACAGCTAGGAGCGCCTGACAAAAGCTGCACAGATGACAGCGGAAAAGACAGCGGAAGGGG
         ACGAAGGGGACTGCTGCGTGCCAACTGCCATGCTAAGGAGGGCCAGGGACATCCAGCTACACCAAGCTGTCCA
         AGTGGCCATGCGTAAGGGAGCCCTAGGAAAGGCTCCTAATCGCCGAAAACGGTACAGTTCTGTCGGATGAAGG
         GGATGGGGGGCCAGTTGGCGGCGGGGGGGGGGGGGGTCGTGGGGGCGCCCAGTCGGAACGACAAGCAGGCATG
         AGCCGCGCCGATGCTT" // 600 bases upstream of the 5'UTR
        "CGCGTCCTTCCCGGCCTAAGCCTGTTACTTAGCACATAACGGTATCGTTCTTTTGCCTCTACAACTGTGTCAG
         CGTGCGTGAGGACAAAGCAAGGTCACGTAATTGCGACGAGTCAAAGCTTGCGACTCGTCACCAGGTCGGCAGC
         TCACCTAAAGTGGTCCATACATTCCATACACAATATGACTGGCGACCATGGCCACCCGGAGGCCGAAGTGTTC
         CCACCCCCTGACTCTATTGTTGCGAAGGCCCACGTTCGTGGGCTAGCGCAATATAAGGAA" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g6062_t1__Cre16g650050t12  =
    [|
        "CACCAGCAGCCACGGCAGCCATGGCGGCCTGCGCGAGCACAAGAAAAGGAAAAATGTTGGTCAGTGCTTGCTG
         TTTGATTTCCGACCCCAAGTCGCGAAGGCCCGGGGACACCAGTAAAGGCCGATGCCACGTCCAATTGATAATG
         CTTGCTGAATCAATGAGCGAGGTGCTAGTCGCCAGGGGCTGGTGAGCGCTGAGGGCCCAGTCGTCGCCGTCCC
         CTGTGGGAGGAATTGCTGGTCCGCTGGTCCCCCACCGGCGCGTATGAGCTGATCGCCTGCTAACGCTAGCTAC
         CCTCTCGCTACTTATGGCCACGGACCTTGGACTCGGCCTTCACAGCCTGGACCTTCACCGCCTGGCGCGAGGC
         CACGCGGGGGCGGCTCACAGCCGAGCGCTGGGCGACGCGGGACACGCACTGCATTTTTGGTTGGTTGGGTGAA
         GTGTAAGGTGGGTGGTGTACACGGGCGGCTTTGACTCGACCGCTATGCAGAGACCGCGAGCGACCGCCCCATC
         GCGGATAAGCCCGATCATGCAAGCCGTGTGAGAGGAAAGCGCCCCTGTTGCCCGCGCCAGGTCATTTGCTTGT
         TGTGATTGACCATTTT" // 600 bases upstream of the 5'UTR
        "CTGCATTTCTGTTTGGTCCTCCTGGCGAGCGCAAAAGTTGTCGACTTCGCAACACCTTTCTGCTGACATCGGC
         CCGCGCGGGTCCATCCGAATAGCGCTTTAGAGCGCGATTGATCTGATTATATATATAAAACAAGGGGAA" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g784_t1__Cre10g457297t11 =
    [|
        "AAAAGTCTGCGGAATTGTGGGATCCAGTGTTTAGTTCGATACTTGTTCTTCTGTTTTTTGCTCCCATATTATT
         ATGTGTTAATGCCTATTGCATGTAAATAGGTGGCAGGGAAAATTCGGGCGGAGCGAGAAGTGCCTTCGCACGT
         CCATGGACAATTGACTTTCGCGAGGTTGCCCGAAAAGCAGCGTCCGGACCAGCCATGGCACCGAAGCGGCGCC
         GAGACGAAGCTGAGAAGGCGGAGGAGGAGAAAGGTGGGGGAGCGATGACCCGCCAGCGCCCGTCCGCTCATCA
         CAGCCCCTCCGCAGCCGCTCGCAATCTGCATGCTGCCAAGCTGTCGCAAGACTGGTAGAGGGGCTCTCCGCCT
         TCCGCCGAAAGCGCTGGGCTGCGTTCATCCAGCGCGACCGCGCGCTCCACCGCGTCGCCAAGCAGCTCACAGG
         CGGCAGGCCCAAGGAGGAGGTGGTGGTGGGGTGGGGTTCGTGGGCCTTCCAGGGAGGGAAGGGCGGCTCCCCC
         ATCTCCGTCAGGGGCGGGCGCGCGCCGACGGGGCGGCTCATCAAGCTGCTCCGTGAGCGCTACGCCAAGCATG
         TGTTCATCATCGATGA" // 600 bases upstream of the 5'UTR
        "ATACAAGACCTCCAAGGTGGGTGGGTTGAGTGGGATTGCAATTGCCGGTGGGTGGGTGGGCCGCTGGGTTTAT
         GCCTGGCTTTAGCGGAAGCCCTTGTGTGGCCCCCCTCCCTCCCACCCTCCCGGCCATGCAGACATGCTACAAC
         TGCGGGTGTCAGGAGATGGCCATCAAGCGGCTTGGGGGGCTGAAGGAGGGGCAGCGGCCCTGGTCGGTCAAAG
         TCTGCAACGACTGCTTGACGACCTGGGTGAGGACCGTGCATGTGCGGGGACGGCTGGGCACTGGGGACGGATT
         GATTGACAGGTTGACGGTTGGGACAGCTCATGTGCGCTCTCCTCTTCCTCCCCGTTCCCCATCCCCAGAACCG
         CGACGTCTCCGCCGCCAACGTGATCCGTGTGCTCCTCCTGCTGAAGCTGATGGGCTTCGAGCGGCCGACCAAG
         CTGCAGCGGCCGCCATGGCCGCCGGCGGCGGCGGGCCGGGCTGAGAGCCTGAACGGCGCTAGCAGGGCGTGGG
         GCTGAGGGTGCACGTGTCGATTGGCGGCGAGTGACGTGACTAGTTTGTTAGCTGCGGGTTAGCACGGACTGTG
         CACCCCACCCAACCGGCCACGTCCGGATTTGCGGGGATGCCAAAGGCCCCCAACATAGAGGCGTGTGCTTAGT
         AGGCGCCCGCGTCAAGGTGGCTGGGTTGATAACGACCCGGGATCAGCCCTTTTCCTGCCATTAGGCAGCCACC
         TCCTTGTTGCTCCCAGGTCGTAGAGTCCAGTGCGCACGCTGCCGCCTCGACCGTATCTAGGCATTTGAACATA
         GAGCAAATC" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g13197_t1__Cre06g289900t11 =
    [|
        "AGGGACTGCTGGATAATGACGGCCAAACGGCAAGGCTCTAAACTGCAGCCCGGAGCCCTTGTCACATTACGTC
         GACAGAAGCGCTATGCGGTCTGGCTAATGTCCATTGTAGCGATTTGAGTTGTCGCGTGCTGCTGTGACACTCA
         CGTAAGGCCTCTTTGCTCCTGCGCTGGGAAACGATGTTTAAAACGATGGCATTTGCGTCTCCAGTGTAACTCC
         TGTAAGGGGGCGGTCAGTTACTCCGCCGGGTCGGAGTACTAGTGCCGAGGCAGGGCTCGTGCGCTCGTCACGC
         AATGGTAGCTCCTGTCACGACAAGTTCTGTCACACGACTTCAGCTAGCTCAGATAGCTTCGGTGCTGTGCGAA
         ATGTTTGGGACGAACCAACGGATCACTTCATAGCCGAGCCTCCGGCTCAAAACTATCCCCAGCAATAGGCGTA
         CTAAAACTGTATCTATAGCAGTGCCAGGGAACGTGTTCGGGGTGCCTATAGACCGTGTCACGATCAGTCGGCA
         GTAGAGGCGGTGAATCGCCGAAGCAGAATGAGTGTGTGCCGCCGAATTCCTATGTTGCTGACTTTTCTGATAG
         TTGCCTGCGACCATGA" // 600 bases upstream of the 5'UTR
        "ATAAAGGGGATGCAGATGGAACAGAGCCTGCGAAGACTGTTCAGAAAGCTAGCACAGAGGAGCCTTCGCGAGC
         GAGACGAGCTGCTGCTGCAAAAGCAGCGAAGACA" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g5992_t1__Cre03g198236t11 =
    [|
        "CGTCAGGATGCGTCAGGGGGCGTCTGGCCGCGTATCGTGACGCGAGTGCTGAGACGCGAACAGTGACGCGAAC
         TACTTTGAAGGGGACCCATACCGTGGTGCGCACGGTGTTTTACACTAGCCTATACTGCTATTGCGACAAATGT
         AGCCTTGCCCGGCAGTATGTGGGCTGTCATGGGCAGATTGGCCCCTGCACCTGCGCCTGCTGCCGCAAGCCCC
         TGCACCCCTGGCCCTGCACCCGGGCGGGCTGGCCAGCCCCGCCTCTGCCCCGTCGCCCCAACCAACCTTTAGG
         CGCAGTTCAGCCAGCCTGCGCCAGCCAGCCGATATGTGCCGATACTGCTTACAACAAGCCTATAAGCCTAGAT
         ACAGTAAGCAAAGGCATTGGCACAGGCGAGCATCGGTGACACCGGTGCAATCGGTCGCCCACCACCAAGTATT
         TGCCCCTTGCCTCGTTTGACCGCTATTCTCTTGAAAACATAAACATGACTGACGAACTCATGAATCACTCATG
         ATACGTCTGCCAACGGGCCTTGTTCCGTTGGTGGCGGCCTGGGGCCTTCAGAGCGCGCTTCACAGCGGGCTTG
         TATGGCTCCAGGCCGA" // 600 bases upstream of the 5'UTR
        "CCAACCCACCTCAGCGCTTACTGCGCACGTAGTGCTGTCAATGTTAACATTAGTGTTAGTGCTGGCGGCTCTG
         GCAACGGTTCACCTGAGCTTTGCTCCTCGGGCTATCAGCCAAAGTTCACAAACGGAAGAGCTACAAATCTTCA
         TTAGGCTTCCTTCAGGGCAACTCGTCACGTTCGTGGGACCTGTCGATGCCACTGTGGGCTCTGTCTTCGACTT
         TATCGCTTCCACGACCCTCTTCCATCTAGATAACAGCCGGATATGGCGCTTGGTCGCGCGCGGCGGCGACACG
         CTGCTGGAGCGTTCGACCGTGCTGGCACGTCGTGGCCTGGCGTCCGGTTCCGAGCTGCAAGTGCTCGCACGGC
         TAATGGGTGGCACTCCAAAGAGAAGGCCATCAGCACCCTCAGGTATAGTTGTTGGGAGTCAAGAGGCTTGGCA
         CGGCAAAGGACTTACAACCACACCGTGTGCACGCAGACAACTGGCTCAGCCTCCACCGCTTCACCCGCTTTCG
         AAACCACTTCGGAGACCCCGTCTACACGTGCACCAACCCTGGATGCCGTTCGCCAACACAGTTGAAGGATAGA
         GCGGCTTGCTCTCGCCACGAAAAATGCTGTGAGGGTGCACCTCGGGACCTCCCTCTCCCTGCCGGTACCGGTG
         GCGAAGGCGGCCTAGGAGACGCTGACTTTCCTGAGGGAGGGAGGGATTCGCCTGCGGGTTCTGTTAGTGAGAG
         CGACAGC" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g8370_t1__Cre02g078226t11 =
    [|
        "TGTTCACGTGGAGGCGGACAGCGCGGCGCCGCCTTCACAGCGGCGCCGGCGCGCAGGTTGAGACCTGCCGGTG
         TTGTCTGGTTAGCGACCTGTGCTGACGAGCATGGCGTAGCGGGCAGTCGGCTGCAGCGGCAGGGGTGTTTTCC
         TTGCTGTCTGGTTGGCCTGCTGCCTGCTGCACTAGTGCTTGATAGCTGGGCGCGACTGGGCACATATGGCGGG
         CGGTGCCTATACAAACCAACCCCGGCTATGTCCGGGAGATGAGGCTAGTATTGGAACCTTCGGCCTCAGACGG
         AGGACGTGGCGTGGGCTCAGCCCACTTTTCCCTTGCAAGGGAGAGCCACCTTTCCTTGTTCATCTTCTGGGCC
         CCCTTTCTGCTCTCCACCCTCCCCGCCCCCTTTCTGCTCTTTACCCTCCCTCCTGTTCCCCACCCTCCTGCCC
         CAATCTCGTTCCCGCCCTCTCCCATAACCTCGCCCCGTCCAGCCCCGTCCAGCCCCGTCCAGCCCCGTCCTCG
         AACAGGGGAACCCATAGCGGTTGGGGACCCAAAGCGGTTGGGGACCCAAGCCATCGTGCGTCGTCTGCGTCGA
         TACATGGATATAATTA" // 600 bases upstream of the 5'UTR
        "TCTAATTGCTCGATTACGTATTAGTTAGTTATCTCAATTGAAATCCCTTATAACACACTCCTTTCGCGACTCA
         GACTACCTGCAACGCGACGACGCAGCGAGAAGCGCTCACACTTACTCAGGGCAAAATTTCAATAGATATGACT
         GTTGTGTAAATTAACTTGAGACAGTTTAGAGGTGTTCGGTGGGATGGGACCTCACGCTTTGAGAACTACTTTT
         TAGCGGTCCATGAGCTTCCTGGCTCGTTTGGACCAATCTACGAGCATGAACTTGTAACACCAATCGATCAATC
         ATTGGCGGGGACGCCCATTAAACTCATTTAACTCATTTTCCTCTATATGATATCCTCATTTTCATTCCTTTGC
         AAAAGCTAGCTGCAGCACATACTAAGCGCCGCAGCGCACAGTTTTTCC" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g5129_t1__Cre13g603550t12 =
    [|
        "AGTTCCCGGCCGTGGTTTGAATTTGATGGGCTCGGGGGCTCCCTTGTGGACCCAGGGACCCGATGAGGATGGG
         CAGGAGACAGAAGCGGCGCCGGCGCGGGAGCCGCCAGGGTGCGGGAGATGGACTGACGGTTGCCATGTGCGGT
         GCACGCGTGCAGACGGGAGGTGGTTGGTCTTTCCTAGAGATGTGGAGGACAAAGTGGCGTGGCGCCAGACGCA
         CATTTGAAGATGGCTAGGTGTTGGTGCAACCCTGCAAAGGATGTGCGTCCATGGGACTAGCCCAGGTACGAAT
         GTCCTACGCGCATTTCACTGTAGGCTTTGCGACTGTTTGTTAGCAGCGCAGTAGTCCCAAAAGTGCCGAAGGT
         GAAGCCGAACCGCGAAGAAGGAAAAGTCCCGGCGCAAATTGCGCTTGAGTTCTTTGGGGGTGATGTTAAAAGA
         TGCTGCTATAAGACCCGGTTTGATGCGTGTGACGAGATTCAAGCAGTGATTTGAAGGCCCATGCTGCAGCCGC
         GCACGGGCCCGAAGGGACTGGCGGGGTCCTGTAACATGATGCACAACTGCAAGGGCATGCGCCCGGGGTTTCA
         GACACGGGGGCACGTT" // 600 bases upstream of the 5'UTR
        "GGGACGGAGAGCTGGCGGAGGTTGGCCAGCTGGGGGGTGCCACGTGGTGCCACTGTGCCAGCAGGCTGCCAGC
         CCCTGCTACCGTGTCTGCTGCCAGGGGTTGGCACTTGGCAACAGCTGCGGACTGCATGCATGCACAGCGTGAC
         TGAGCAATGCGTGTGCAGGCAGGGAGACCTTGCCCAGCAAACTCACTCCAGCCTCTCTCTGGACGGAATGAAG
         CCACGCCGTGTCGTATGTGCCAAGCCGACGTGTGAGGCAGGAACCTCTCTACGGAACAACTGAAGAAGGGTGG
         GTTGTGGGTCCATGGGTCAGAGGCACGCGTACCCCTGACGGTGCGCGCACCCCTGGCATGCCACGTTACATGC
         ATGTGCGGAGCTACGAGTGAGCACCGACTCACCGAGAGCCTGTGACCGTTTGCCGGTTTGCGGCAGCGGCATG
         GGGCCGGAGCTTGAGAACTCGGGAAGCCTCCAGCAGCCCCCCGAACCCCCGAGTGAGCAAGACAGAGCAAGGG
         ACCGAGCCTCAGCCTTTGTGCGGTGCCGGTAGTGCCCCACGTGAGCATGCCCACACGTTGCCATGGGGGTTGA
         GGGAAGGACGGAGGTAGCGAGCCCAGGGGGGGGCCCACGCAGGTACGTGTCGTGGACAGCACGGGACGTGGGC
         TGGGCCGCGTGGGCACATACAGCCACACGCAGGCCGCAGGAAGGCGGGGTGGAACGCGGCAAGAGGGCTGACG
         GGGCGCCAACAGCGGCCCCGTCGACACCAGAGAAAGCTCTACGCAAGAGCCTCAATATGATTTGCTGATGCCA
         CTGATGCAACCCTTGCCACACTACCGTGGCTCGGGCTTTATGGCAGGATGTACCGCGGTCTCGGTAGCCGGTA
         GTCGTCTCTGGAGCAGGTTTTTGAGACCAGTCTTGGTGTGAAGTGTGTATTATCGCTTGACGTTTCCAACATG
         ACCAAGTCGTAGTCATTTTAAAGCGGTGTGACCTTTTGGCAGTCAAGAAGTCGGCCGCGCGCGCGAGAGGCGG
         GGCCCTGCGCTGTGTGCATCATGTTTATACGTTGATTTGTTGCTTTTTGGCCCTAACGTTCGGATTGGACAGT
         CCAACAATGTGCATGGAGGGCCGACGCACGCTTACGTGGAAAGGCTGGTTGCGTGATTTAAGCGCTGCCGCGC
         CTGCGGCGGCCCTATTACAATACAAACAACAAACAAGCAACCTACTTAATAGCCTTCTATTTCCCAAATCTCA
         CATTGCTGAAACAAACTCAGCAAAGCAAGCAACTCCTAACGCCCGCAAAGCACC" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g11925_t1__Cre18g748547t11 =
    [|
        "AGGGTCAGGGCCTTGTGTTGACGGATGAGAATCTGTGCGATCAGCACTCATATACCTGAGCCTAACGACGGGG
         AGTAATGGACCCCACCGGGTCAGCGCTTTTCCCTCACGTCGCTTACGCTATCGTTCAACATCCTGGAGAGGGT
         GCTCTTGCCAGTACTGCTATTGGCCGTATGAGGCAGACCCACAACACACCCTGTCGACTGGTGTAGATATCGC
         TCGCGACCCGCGTCAGCAGGTCCCGCTTGTAACCAATGCAATCTGGTATGCACCTCTCAGAGGCACCCAAGGA
         GGCGCCTTTGGCCCGTGCCTTGCACGTGCGTTATGGCGCGTGCTGCCTCGCATGTGGCCAAACCAGCAGGCGA
         CGCCCAAGGTGACAACAGCAGCGTACTAAAGCGGCGGCTGCTCTTGCGGCGGCGGCGGTGGCGGGAACCCGCC
         ATGTCACAGCGGGTCACAGGGTGCTACTTGCGGTACTGCAAGCTCGACCGCAGCCTGTCAGTGCGCGTGGTGT
         TTGTTGTTGCTGGTGTTGTGGTTGTTTGCTGGCAGGTCAGGATGCTGGCTTTGCCACCCCGCATATTAGCACA
         TGGCGTACGGATGTGC" // 600 bases upstream of the 5'UTR
        "ATTAGGCAGCCCGACTTGACGACTGCCTCGATTCTCGACGCCAAATAATGCTGAACGGCTCAATCCCTGCTGC
         TATAAGGTATTCACATTGGCGCAAGTACACAAACCATCTTGTCTGGAGGCGCGCATATCGAGACTTTGAGTTT
         AGTTGTCATACTATACCTTGTAGTTTGA" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g7722_t1__Cre17g733900t12 =
    [|
        "AAGCTGCATAACAAAGGAGAGAGGCGCAGGTTTGAACACAGTGGGTTTCGCATATTGGCATTCGCACGCATAC
         AGCAAGATGGTCTTGAACGGCGTACGCTCTTAACAGCCCCCCCACTTAGCAGCAGTATCCAGAGGTGCGTCGT
         GTCAGGGGAATCAGCAGCAATACCTGGCCATACCAAGCACCGTGCAAGTGACAGCTCAAACCGATCCCCTACG
         GCGCTCCCAACTACCAGAGCTGCACTATAAGCGATAACAGGTTTCAAAGTGGTGAAACCTAAAACAGCATTGT
         GACGGTAAGATGGATCGTTTTCGGGCTCGTGTACTTGTGCAATTTGAGCTTTTGTTAGTTGTAGCTTTGGATT
         TCTGCTCTGAAGAGGCAAAAGGAATTTTGGTGGTAGCGAGCCGGCGAAGCCGGGGGAGTTGGCATCGCTCGAC
         CGCGGATATACTGTAATTTTAGCGTTCAAGACTGTAATTTTAAGTTTCTTGAGTAATTTAGGTTTCCATGACT
         TGGATCGTAGGTTGCTCACAAACTGGCTGACTTGGTTTTATGTGTGACTCCATTGACGCCGACAGTTTATTTA
         CCCCAATGGCATCGTA" // 600 bases upstream of the 5'UTR
        "AGTATCTGTAGCAAACTCTCAGATGTCGGTTTTGGCTAAGCGAGAAGCGTGTCAGCGTTCTGTGCAAGGATGG
         AACGGCGCCACAATCCTCGAAAACAGGGTCGCTCGGTTTCGGGCACAGTTCTGCCCTTTTTCAGCTGGGGGCT
         CACGGTATTGCACTTCTTCCGTCGCGACGCTTCGAATTCAGACTGCTGCCTTCAGCTTTCAGCAGCAGCAAAA
         TCATATCATGGGTATGCATATGTACAACCGCAACGTTGCATGCAAGTATAAGGAGTTCATAGAAGCAGCTGAA
         CGGCTTACTTGATTTGGCTACTCTCAAAGCATGGAGTTACTGCGATGATGAGCTCATAAATGTGCCTAGACAA
         ACAAAAAGGTGTTGTAAGCGATATCTACACACACGTACGCATGTTCATCATCAGTGCTGCACACCGGTTCTTC
         TCGGCGAAGCCGTGCGGT" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g7896_t1__Cre17g742550t12 =
    [|
        "AACCCCCCAGCGCCGCACCGAACCGACCTTGAAGTGCTCATCCAGCGACGCCTGAATCTTCTCCCTGGCTGTG
         CGGCGGAAGCTCACGTCCAGTACAAGCTCGTTGGTCGCTGGGAAGGCGCGCAGGAAAGCGCGAAGGAATGCGC
         GAGCTTGCTGCGCCCTTTGCTTCCCGACGTCAAGTTGTGCTTTATGGACGTGGAGGCTGAATGCATGCCGCCA
         ACTCGAGCAAACAGCGCGAGCAGAGCTCAGCGCTTGGCAATCCAGGTGGCGCCCAACGTTTAGCAGCAATCCC
         TCTGGCAATCGGTTCGCCCAGTCAATCCAAGAGCTATCGCCATCAAGACCTCTCTGGTGGGGCCTACGGTCGA
         GGTGCAGTCCGGCGAGCGAGTTCGGGCCGTCGATTTCCATACTTCTTTAGTCGTAAGTCATAATTACCACGGT
         GTTTCACTGCCAAATGCATGCAGGTGTGCGGTGGTTTACTCAGCTGAGTTCCAATCGGTGCTTTGGGCGCTGA
         CCCATCGGGTCCACGGTTTCGGCGATTGCTCAAAGCCCTTCTATATCACTCCCCGTGGATACATTATCGGCCC
         CCATGCGTTTTATTTA" // 600 bases upstream of the 5'UTR
        "GCCATCCCCAAACCTGCCTAGTAACTCACACGCACATATAAAAGAAAAGCAGCCACGCTTTACCGCTGTCATT
         ACATTTCGGAGTCTGGCAAAATAGTTGCTGATATCGATCTTACACCACGTAATGCATCATGTTCAGGTGCCGG
         CTTCGGGGTGCCGCTCCCATCTCCGTCACGGTCCCAGCGCGGCCAAGCACAGGTGCAGCCCACGGATGCTAGG
         ATTTTGG" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g7422_t1__Cre17g718950t12 =
    [|
        "GCGTTGCGGAGCTGAGGGACCAGGGACAGGTGGTACGGCCGAGTGAGTTCTTTGGTTGGTTTGGAGGGTCAGG
         CGACTTCGGCGTGGGGCTCTTGAATGGGTGCATGCACGATTGGGTTTCGACGACACACACGACTGGTGATGAC
         GAGCTAGGACTTTGTTGGCCGCACATTGCTGTGATGCCATGCTGCGCGCGCACGGATGACGGCGTACATATGA
         ATGTATGGTGGTATGTAGCAGAGGGGCGGGGGCGAGTCGCTACGTCCAGCGGGCCACTCAGGGAGGTGTTGCA
         CAGGTTGTTGTGTCGGCAGCTAAGCCTGTTGCTTCGTTTGTTGGAATGGTGGAATACATATAAATGAGAAGCA
         CTGGTGTCACGATGGCATCCCCGAGCGGAACGGCGCTTGGCAGAAGTCTGTTACCGGGAATTCCCGAAGGCCC
         AACGCGCCTCCCAACCGCGAGCGTTACAGTCAGATACGTATTCAGATAACCACTGCTAGGCAGTCGTCTTGCT
         TCTCAGGCCGGGGCCGGGGCCGGGGCCGGGGCCGGGGCCGGGGCAAACAGCAGAATCGTGACGAATGGTCCCC
         ATCGCGCCTTCTTTTG" // 600 bases upstream of the 5'UTR
        "GGTCAGGTCGTCGCCAAGCATTGAGGTACTTCCGGGCGCAACAATTACACTGTACACAGCATTACGCCTACCT
         ACTATCAACTCCCCTGCAGTAAACGGCCTCACCAGTGTGGTCCCAACAACCCCAGCATTCCGCCCTCGGGCAA
         CTTAACACCAGCCAGGAATGAGGACCATTCACACGAGTGGCACATGCAGGCTCGTAACCGGCACTCAGCGCTG
         TCGGCGCCCGGTCCCATTGCCCGGGAGAAGCTCATCGATTGTGGCACGCTCGTGGCGCGACTACTTTGCTGAT
         GAGCGAAGCGCTGCCGAGCAGTTCTTCGGTGAGGACTATGAAGAGGAG" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g98_t1__Cre10g422300t12 =
    [|
        "CTGGCCCGTCCACCTGCTTTGAGCCCAGTCCAAAGATGTCTGGGATTACCACCGCAGCCTGTGCGAAGCGCGG
         ACGGAGCTAAGGCGGTGGTCGAGCAAGCGGTGGCATCTACCCATTCCGTGGCCCAGTAGTGCTTCTTCCGGAC
         GGTTGTTATGAACCCTGTTCTGCCCAACACACGCGCACACACATGCGCATACGCACACACATGCACACACACA
         TGCACAAACACGCACACGCAAACGCGTACACACATACATACACACACGCGCACACAAGACGAGCCGTTTGACT
         GGGCGACTCACGAATTGCCCCTTGCCGACGCAGCGAACGTGGGTTCGGCCCGCCCGCACCAGCTCCCCCCCCT
         GGGCGTATAGTGCACCGCGAACAGTGGAGTGGCAGTGTCGCAGTGCATGTGGCTGGTGCCCAATCCAGGGTGC
         GTGTGGTCTTCTGGACTGCTGGTCCAGGAGGCCCAGCCCTTCGCATCCTTGGGGATCCCTGGCACTGCCGGGC
         ACGAAACCCCAAACCCCACAACCCCTGCACGTCGCCGGCATTGAATTCCGCTCCTTGCGACTGCGTCACATTC
         TCATTAAGACTAATTC" // 600 bases upstream of the 5'UTR
        "TACATTTTGCTTGTGTGCATAAAGGAAACAGACCAGGTCAACTCTGTCACTTTGAGCGACGACCAAAATGTTT
         TAGATCTAGCCCTCAATACTTTTAAGACGTGTACAAATTGCAGATGTCCATATGGGACAACATGTTCGGGGCT
         TTCATGCAGAGCACCTGCAACAATGGCGGGATCACAGCTTGCATCTGTGACTCTCTGCCAGTCTTTACAGGAA
         CTGCTGCCCCACAACTTGATATTGTAAGAAGATATATTTGTGCACAGTACACAAG" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g1550_t1__Cre01g028950t12 =
    [|
        "AAGAACTGCGCACGGTGTTTGGCGGAGGCACCCCTGCTGGGCGGAAGGCCTTACTCGCAGGCAGGCAAAGCTG
         GGCGTGTCTGTGCCATTGGCAGTTCATATATATAACATACACCGGGTGTGGCATGGGTCCGATCGGGAAGCTT
         TATCACACTCGGCGTATGAACAACAAGCTTCCCACCACATGACCCGCCTTTTGGCCTGGCGGGTTTGGGGCCC
         AGGCGGGGCCCACCCTGGGTAATCTGCGGCACCCAGCCGGCACGCACCCGAAAATGTCTTACACGCAACCTAA
         CAGCCTACGTGACGGTCAACGTGCCAGGTTGCTGACAAGGCGCTGGAGACATTTTGGCGGTTATGTGCATAAG
         TGCGAGGTGATTGGGGCCATTAGTTGAGACACGTTGGGCCATAATAGCCATGAGTAGTGGAGCTGGATGGGCT
         TCTCCCCTGGAAGACAGCTTGGGCGTGGTTTAGGCCGGGTTGCGGTTCGGGCATCGAAGGTTGTTGTAATCTG
         ACGGAAAGAGTCCAAAAAGCGGCAACCGCTTGCGGGAACTGTGGTGCCTTGTGGTGCGATGTTCGCACCCGTA
         CAACCGATACGACGCC" // 600 bases upstream of the 5'UTR
        "GTCCTCACCAAGACGATGCCAATGCCAAAACGCCGCCCTCAGACGCTGCACGTTTCCGCTATGCTATAAGCTC
         GCACTGCACCTTTCACAACTTGCTGCCAATCGAATACCTGCTGGCCTTGCTGCTGCGAGTCTGCACAGCATAT
         AACCGTCCTCGACGCCTCCGGCGGTACCACTAGTCGCTCAGGAGGCTAAGAACC" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g1053_t1__Cre01g004157t12 =
    [|
        "GCAAGAAAACCATAGTTGAGGCGGATTTCGTCGCATACGGTGGCGCTGGCCTAAATGAATTTGTATAAATGTA
         GAATACTCCATTATAGCCACTCACGTACTTATCGCAATCCTTGCTGGCGGGTGAGCTGCGAGCGCGGGCTTGT
         CTCCTTAATGCAAGCAGCCATGCTGGCAGCCATATCAAGTGCGTCTTTTAGTCGCTGAAGATGCCATCTAAAA
         AGTGCAGACGGTTTCGCCGCAAGCCGACTGTGCCGATCGATCGAACCTCGTGCCGCCTCATGCCCCGCTGTGC
         CCTGCGGATCTTCAACTGGCCATGCTACAGACTGGAAGTTGCTGTTCAATACCTGCAGACCATACGTGCATGC
         CTGAAGAGTAGCGGTGCCCAGTACTCTTCTGAGTGCCGCTGCAGATGTTGATATGTATCCTCTCCTGACGTCA
         ATCCACCCTCCTATCGATCTCATGGCAGTTTCCATAAAATGACGCTTGCGGGGCACGCCCAATGGATCGCTAC
         GCTAGGCGACCTGCTGCTGCTAATCGATTGTTGGCGGGCACGCATTTTGCTTTACTGACACACTCACAGCGTT
         GTTGGACCAGTCATGT" // 600 bases upstream of the 5'UTR
        "TTAGCCTGCCTGCTCGCGCGCAGTGGCCGTCTTGACATCATGACACAGTTTTGGTTGTTTGTCACAGCACCTT
         TAAATCCTGACCTGCGATGTCTAGTATACATATGACAGCATCTTTTATGTGCTGCAACTTACTTGTCTTGAAA
         GCGCTGTGCTCAGTGCTCACTCGTCCGCCACGACTGGATCCCATTCGGTTTCGCAGAGCCTGATGCGCCCAAA
         CGTAAGAGGCCAGCTTATAGTGGCCTAAGTTGTGGATCTAAACTAAAGGGCACATTCCAACATTACGAGGGTG
         CTTACTTGTTCACTGCTCGTTAGTTGCAGGGCCTTGCGTCACTCCTGGAGAGCCTGGCCCGAATAAAGTAGTG
         CAAGTCAGCCTACAGTCTTAAGAACGACGACAACCAGCTTCTTGACGCTCAGCAGGGACCGGCGAGGACTAGG
         GCGAGGGCTCGAAACCCATGCTCTGACAAGGTCAGTCGGTCTCGTTTCAGACGCCCACTTCAGAATCACATGT
         TCGTCACGTAGGTATGATTTAGAAACCACCACGCGATGCTGGTGACCGTGCGCTGCAGCTCTTCGAACACATC
         TCAGGCGAC" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g3913_t1__Cre12g542450t12 =
    [|
        "GGGCCACGCCGGGGAGGCGTGGGCGCTGCTGTATGACGACGCTTGAGGGTGCGACTTATTGACTCCTTACTGC
         CGTGTAGCGTTACAAACCGCCACGGCCCCAAACGACAATCCCAATCTCTCAAACCGACAATAGCCTCCACTCA
         TGCCTCAAGCGGCCTAGCAACTCATTCGTGGCCCTCAGCGGCCTCCTACCTCCGGCCTCGCAGCTCCCGATAA
         CCCCACCAAGTCCGCCGTGCCCGCCCCAGCCCGCCCGTGTTGAGGTTGCACTAGTGGCCGAAAGTGCTGCCAG
         TACTGGGTGTGTCGCATGTATGAAGTGCCTGATAGCAGCAGAGTCCAGACAACCACGCACGCCGCAGCGCCCA
         CGGGTGCCACCACATTAATCCGCGGCGGCACCAGGGGGGGCGGGTGGGTTGTCACCGTCCCGGCAGAGGGACG
         ATCCGAAATACAGTACAGAAGCACAACGGCAGATAAGGCGCCGTGTGCTCCTGACGCGTACAAGACCCAGCTC
         GGTTCGGCCCCATGCACAGGCACGTACCCGAGCGTCCTGCGCCGTGCGTGACTCTAACGCAACACGGCAGTTA
         CGTCGCAATAACTAGA" // 600 bases upstream of the 5'UTR
        "CTTATCTCCACTGCGCTGCGATAAGTCAGCGCATAACAACGGCGTTTGGCGCCCCGGGTGCAAGACCCGTGAC
         GGGACCCGTGACGTGAAGTTGCAATCTTGTCGTACATTGGACATGGGCGCGGTGAAATATCATTCGAAGCGCT
         GGTCTAAGCATTAACAGCTAGGGCTCAGTAGCAGAGGTTCTGGTGGTCAAGCATTGGGCTCGAGTCTGGGCAG
         GTCTGGGCGGCGGGTTGGCTCACAATATCAAGCCATTTGGACAACTCTAGTTTCTATTCCTTAGATCAAGTTG
         GAAATAATCCATCATCGCCATAACACACTCAAAATATTTTTGCGCCTGCAAAGTATCACGA" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g15550_t1__Cre09g403850t11 =
    [|
        "GTAAATACAGTTTGGCAGAATCACCGACGCGTACGACATACTGTACTTGGAGTTGGAGTCTGTCCCCAGCAGT
         GTGTGTAAGCCAAGCATCACTTGGAGTTGGAGTCTGTCCCCAGCAGTGCGTGTAAGCCAAGCATCACATTCAT
         TCGACTGAGCATCCTACCCATCTGGGCTGACAGCTGGAAGCCTGGAACGAAGCAAAGCCCCGGGGCGCTCAAT
         AAGAGGCCTTGGCGCGCAATTCAGACCTGCTTACCGCCAGGCGCAGCAGCAGTGCACTGTATTGGCGTTGAGG
         GGTGGCATACCCGTGGACCAGAGGAGAGAGCCCCGAAGCTTATAGGAAGCCCGCTCTGTCGTCTCTGCCACGC
         GACTGGCTAGGCGTCGAAGGAGAAGCAGTGCAGCGCGCACCGCACTGCGGGCCGGGACGCTGGAGCAGGCTGA
         GCCCGACCCAGCCCTGCGGCACCAAAGAGTGTACGGTAATACGATTGGCCCGCAGCCAGCGACCCAGCGCAGG
         CCTCTGTTTAGTGTGCTGCCATACTGGAGTACACTCATCATCAATTTATGCACATCGATTTGTGATTGCGTCC
         CTCCAGACCCACTCGA" // 600 bases upstream of the 5'UTR
        "GTCGAGTGCCCGAACTTCCAACGTCGCCCCTTACCAGTGACGCGCCCTTGATACCGTTTTACACTATTACTTG
         TCAAACATTTGCACCGTATGCTCATAGCTGCTTTAAGTTTAAGTCGCTGTCACCTCACTCACTTCAACCTGAA
         AGTGCCCCGTAAACAAGAGCAAGGCTGCAGTAGTGTATCTCGTATGCTGAGTTAAACCTCCCAGGGGCAGGGA
         GCGTGGACAAAGCTGCTTCAGGAGTTCGCACCTCACCCGGAGTCCTTCTGCGCCGACGTT" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g3069_t1__Cre12g499850t11 =
    [|
        "TGCATCAGCTCATCGTACGCATGGGGAACATCGCTCATCCCCGGGACCTGCAAGTGATCCAGGTGGTCAGAAC
         TGCGTGCACGTCTCGACGCATCAAGGCAATCTTGCCATTTTACCTTGGGGCTGCGAGTGGGATGCCGGTCAGG
         CGACTGCAGATACACCAGCGGCGTAGCGTTGTATTCGGCCATAGTAACTTCCTTGGTGCCGGTCGCCGGGCGC
         CCGGTGAACGTGTTCAGGTAAAAGAGGCCACCAACCAGCCCCGCAAAGCCAGTGGGGAATGAAAAGTTAGCAC
         GGGTGAAGTAGGCAGAGTAGGACATTTCTCCTGAGCCGCTGGCAAGCAAAGCGCCGAGGATGGTCGACTCCGA
         CAGTGGTGAGAGCGGCAGAGGCTTCAGTGGATTATGTGATGCCTCAGGAGAAATAGTCCTCTCTATATGCAAA
         CACAGCAGCCTGCAGTTCTTCAGCCAGCGACGCAACCGATCCGAGCTCGGTCCCGACTCCCGACCCCATTCAG
         AGCTGTGGCAACTGCCCGTCTGTATGCAAACATACATGCCAATACACGTTTCTGTAAGATGCACGAGATGTGC
         GGTATGCCTTCATAAC" // 600 bases upstream of the 5'UTR
        "GGGCAAGCAGTGACATTTCCTCTCCGCAACACGCTAGCAGCAGCAGCAGATGCCTGGACAATTGTTAATTATA
         GAGTAGGCACTTACTGGTTTTGAGTCAGGAACAAGTGCGAGCGCAAGACAGCCCTCTGGACAGCCCTTGCACC
         TCGTCTGCGGGACTCCTCAAACCTGGACAAATACAATAGCAAAC" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g3976_t1__Cre12g545650t12 =
    [|
        "CCAGCGGCGGCAGCATGTTAGCAGAACCTGCAGGCCGCAGCAGCGGCAGCGCGTCAACAGAAGATGCCTCAGC
         CGCTGCCGGCGCATGTAAGGTGAGGTGGATTCCACTGCCCACAGCAGCGACAACAGCAACACCTAGATCGAGG
         AGGCGGCCTGCCTGGCAGCAGTAGGTTTAGCCAAAGACAAGACACGCTGGAGACGTGCAGCAGCTTGCACCAG
         CAAATCTTCAGTGGTCACCGACCACGTGTGTTTCCTGCAAGGCCTGGGCACCGCCCTGCACCGCCTGCATTCA
         AGGATTGGTGACGTCGCATGGGTCCCCCAGAGACCAAGCTCACCGCACCAATTTCTGTACCAGAAACCCCAAG
         CAGTGTCGGCGCTGTCTACCGGCGGACACTAGGGAGGTCGCATCAGCGGTTAAACACCGCCAGGAAAACCAGG
         GCACCACCCTGTGCCCATGTTCAAGGATTGGCGACGTTGCAAGCGACGGTTGGCACACTTTCTGCCTCACGCC
         ACTGCTCGAGGGGATTCCCAAAGGTGCTTGGGTGTGCCCCAGCTGCACAGCCAAGGGCGTGACGCCGGACAAC
         CTTTGTCAGGCTCAAG" // 600 bases upstream of the 5'UTR 
        "GTTAAGGTCAAGCCGCGTCGTCTACACCGCCAGCCAAGTTCAAGGACCTTCAAGGAGCGCTTGCT" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g13367_t1__Cre06g298600t12 =
    [|
        "GGGCCTCAGGCCCGCGAAAAAGAAAACGGTTAAGGAGCTGTTTTGAAAGCTGAAGAGGCCGCTTAAGCAGGCC
         TCGGAAAAGGGGGAGGCAGGCGCGCTGACCCTTGCTTGAACTCCGGGCTCTGCGATGGGCCGCGATGAGGTAG
         AGAATTTTTTCCCGCCGAGCACCAGGCGCGCCCTTCAAGATTACAACGCAAGGCCTTGCAAAATAACCAAGGA
         CTTATCAACGACAGCGGTGGCAGAGCGCTGAAAAGAGAATTTGTTTTGGCGGGGGCTCGGCCGCACACTTCGG
         CCCTTTGTGCAACCGATGCCTCTATTGACATACGCATAAGTATACCTGCGAAATAGGTCAGCTGCACAGCCGA
         TGGCCCTCCTCGCCGCGCGACTGCAGCCTTCCCTCGACCCGGGGGTCCCGCCCCCCAGCCGCGGGCTTAGCCA
         AATCAATTTAATGTATGCTGCATTTAATCCAGGCCACCGATTTGGTTGGCTGAATTGGATGCGGATGACTGTG
         ACGAAAGCGTCAATGATTCACTAAATCGTAGCGCTCCCTGAAGAACTTCTTTACGTAAACCACAGGTGCATCT
         TTAAATACAGAGTATT" // 600 bases upstream of the 5'UTR 
        "CTGAATAGCAGGACAGACCGCCCGCAGAACTGACAAATACTCCACA" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g4913_t1__Cre13g592551t12 =
    [|
        "GCTGCTTCTGCGGGTCAGAGGGCGGCGGGCCCGGGCAGTGGTGGCGAGCGGGTGGCAGGGTGCCCTTCGGACT
         GTTAACCAGGTCTTGAGGTACCCAGCCCACGCAACTGAACTGTGGCCGGTACGAGCAATGGTATGTTAAGGGC
         AGCTGCTGGTCTCTACCAGCAACACAACCGGGGGCCGCCTCATTCGCATGCGTGGGATTTGGCCAGTGGATAT
         GGTGTGGGTAGTTGTGGGTGCACGCCTGTGCGCCAAGGAAGGCTGGCGTGATGGACGTAGGACAGCGCCACCC
         GGGGGCGAGGCTGCTCAAAATTGCACTTGATTGTCATTGGGTACGATTGTGTATGCGCCGGGGGCGTTTCGTT
         TCCGACGCCACGGCCGTCGACGCCAACCGCCAGCCCTCGCCAGCGTAGTATTAGAGCGAGTACTTGACAGGCG
         CGACAGAACGTGCGGCGCGGCTTCATTGGCAAACCAAAGTTGGCAAACCATTTAATGGTCACAATGCACCGGA
         TTCGAATGAGGTAGGTGCGCCTACCCACATGTCACCTCGACTTTCGATTGCTGCGAGCCCGTTAGCCAGCAGT
         AGGCACATAGTTTAAT" // 600 bases upstream of the 5'UTR 
        "GTACAGAGTATAAAAATCTGATGTGATTGCTGCAAATTTAATCAATCTTGTGAACCTGGCGAGCAACGCCGGT
         GCCCTACCTCGTCTATGCCGGCATCCGAAGACTATCTGATGACAAAACTTGCATTACCTGTGCCTCGCGCTAG
         ATAAAGCCGCTGGAGGGGGCGCAGCTTCGCAATAGCGAGCTCGCGCTGCGAACTCGACCGCGCGTGGAAATCA
         AATGAA" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g5955_t1__Cre17g734564t11 =
    [|
        "GACCTTAGCGTGCCCTTACCGTAGCTCCGGATGAACTGGCCGCCCATGTCATCAACGATGCTTAAGAATGCGG
         AACAGAGAAACCAAGCGAGCACAACGGTTTAGAGGCCTGGACGCCGGCACAGCAGAACAGAAAGCGTGAGCAA
         AGTAGCGACAGCAGTTCAATGACAGTTGCCCGTCGAGATAATTGCTGCGCGGCATCCGGCGTTGTATAGCTCA
         GGACACTCAGCTGGCGCATACGGCGTCACAGGCTGGGCTGGCGCGGGCTAGCCGAAGACATGCGTCACGACAA
         CGCGTACCCCAGAAATGGAACGCAGCCGCATGCACCGCTGAAGAGTGCAGTAACGACGTACTGAAACTTAGCC
         AGTCCATGAAATAATTGTCTAATGGGTGCGGCTTTCTTGGCCCAGCAATGGCGGGGGGTGGGGGTTAAGGGGG
         GTGGTATGGGTGGGCGGCGTGTGGGTGGGGGGTGGGGGGGCGGGAAGGGGGGGGGGGGGAGGCGAGGCGAATG
         AGGGGCAACCGAATGAGGCGAGTGCACAGCCAAACCGGCCAACCCCGCAGCCTACAACCCTGAGACAGCATTT
         GCATGTATAGTATATG" // 600 bases upstream of the 5'UTR 
        "GAAGCACGCTACAGTTGTATGCGCAACGAACAAAGGACACGCGCCACCAGCAGAATTGTCCTCAATAGCCCGA
         CCGCGCTTTTGCACGCTGTTAAGCCTATAGCAGCATTCAACACGGCTAGTTCGTGTGTCAAGTTTTGCGTGCT
         TGGTCTCAGCGGCTGAGCCAGTCTTAACGCTCACCTTATCGTCCGCCTGCGCCTTATCGGTGCCCAGGACCTC
         TTGAGAGCATTGCCGGCCATGCAAATGGGCGCCTAATGCGGCTCACATAACGCTGTGGGCACAGCGTGTGCCC
         GTCTGCCGCCAGCGAGAGCTCAGCAAATTCGGCCATCGTCCACAGCGCTTGGGCCTGAGCA" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g327_t1__Cre10g433950t11 =
    [|
        "CCAGGCCTTTCCCGCCTCTCGCGCGTCTGTTCTGCAATAGCATTCTTGTTGCTGACCCTTTGCTGCACCCGAC
         GACGCTATCAGATACACATTTGCTTCATGTATACAAAAGCCTTAAGCATTTCGCGCGCGCGTAGCAGGACAGG
         GCGCGCACCCCGACCGCGGCAGCTGGCGATGCCGCGCCCCAATAAGGTATCCAAGAACTGCGTCCGCCAAGTG
         GCTTATAGAAAGTTGCATTGTGCATATCCGAGCCCATCGGCCGTGAGGGGAGGCTTGGGCCACAGGAGGGAGG
         AGGGTGCGGACTTGTCGGAGCCTGCCGGCTGGCCCCGCGCTGCTTGCTATGGCTGCCAAGGTCGCTCAATTGA
         TTTAGACATAGAACTAGTCGGGGTAGACGTGGTGCGCGAACTTGAGAATCAATGGGCGCGAGATTTCGGCTAA
         ACATGTGCCGAGTAATTGAAAGAGGTATTTATTTATCAACAAATGCTACGCAGCTCACAGCGGGAGGAATTTC
         GAGTGTTTCAGTGACGCATTTGAGGCAACATGCAGACTCGGGCTGTGCGCAACGGGCTTATAGCTTAAACGTA
         TAACATTGCACTAAAT" // 600 bases upstream of the 5'UTR 
        "GTACTCGTTGCAAGTTTCAGAACTTAGAAAGGTGCAGCGACGCGATTGATCGGTCCTGCCCATGGACCTCTGC
         GGACTTGACCTGGGTGCACGGGCACCCTTGCTTGGTCGGATTGTCCCTTCTTTGTAATATCGCAACCTTAAAG
         AATAGCGACGCGATATACTGTGTCTTGCGCCCGGCGAACTCGGACTTGCCGAATTGATGTTGCCGGGCGCGAG
         CCCACGACAATATGACTTACTTTTAGAGTTAAGTCATTGTTGTAAAGTGCATGAGACAAGGTGATACATAAGG
         AACGCGTGAAATCTGCGACTTGGGACTCACGCAGGAACAGTGACCAGCGGCCCAATCCCGCGAACCCTCCTCC
         CGGTTCAAAACCGACCGCGGCAACGTGCCAAACCCGCGCTGCACTAAAACACGGTTATCCTTTACACCGAACC
         TCTGAACAACCCTGCCAAAA" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g587_t1__Cre10g447300t12 =
    [|
        "GCTTGCGTGTCTCGCGCTGCTGATCGGAATCGGACACAACGGTCGCAACCCAAGCCAATCCCTCGGAGTTTGG
         ACTCGTCTGGGTTAGTAATCAACTGCATCGACAAAAGACTTACACAGCCCGCCACGAGCGTCAGGACATGGTC
         GCAGGGCACAACACCCGGCGCTGCCACAGCCGCTCCACCATTGGAGGCTGAGGCAGCTGGCTGAGTTGGTGAC
         TGCATTTCTGCTGCTACCGGAGGCACGGTCATTCAGACCACTGCGGTCCGCGACAAGCTACTGGAGCGAGCGA
         AGCAGTGTGGGTCAATAAAGCAAGAGCATTTTCATATTAGCAAACACGAGACACAAACGGCTTCAGATGCAGT
         GCGACAAGCGCAACAGTCTTTAATTTTGTTTCTCGGATTTAGGTTTAAGACGACAGAAAGCCGCTGTACAAAT
         GTATACCTGTATAATATCCTAGTTGAAGTATTGCCAGACGGTGTGAACAGCAAATGCAGCTGCCTCGTTGGAA
         TGGCTGAATGGGGCCTCAGCTCATTGCGGGTTTGCGGGAAACGGGTTTCTCCAGTGTCAATAACAAATATACT
         CCCAGCGGACAGCGTA" // 600 bases upstream of the 5'UTR 
        "GACGTCCTACCCATTGGAGCAACGCCCCTGAACTGAAAGAACGAATTTACGCGGGATCGTCTAGGTCTTGGTT
         GTCCCTTCCTCCAGCCCTCCCTTCAATCCTTTTTATCG" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g13365_t1__Cre06g298500t11 =
    [|
        "GGCGTGGGAACATCGCCATTGATGTTCAGGCGCACCGGCCTCGGCTGTGCAGCACGGTACGCCATAGGCGTGC
         GCGGCCTGCTTCAGCGGGTTTACGCCTGTCCCGCACACCCGCCCTCAGTCAAACACCGCTGCCAGCTGCTGTG
         CCAACCGCTGACACAGGCGCTGTGATTCCCAGGTTAATATGAGGCACGGCCAAGTAACAGCATTATCGTAGCG
         CCTGTAGCAGCCGGGGCCCGAACTATCAGAAACAGGTCTGACTTTGAAATTGCTAGCGCGTTCCATCCATCGT
         GCTTACTGGAAGTGCAAGTTACATGTTTATAAACGTGTTGGGCAAAATGCTGTCCTCATTTGCGGCAGAGGGC
         ACCAACTTAAAGTTTGTTTGCGCTGCGACTCGCCAGGTTGGGTTATCGGTCAGGTAGAAGCAAGTTTCTACAT
         AACGGGCACAAGGACAACCTGTTAAAAATGCCTAGCGACGTCAGGGAACATACAAAGCTGCCTGCGACCCATC
         GGCCCATGCAAGTTCCATTCCACATATGCGACCAACGCGTGGGTTCTTGTCTTCAATCCTTTACAATTACTCT
         TCACGGCCTCCGTAAA" // 600 bases upstream of the 5'UTR 
        "TATTTTCGCTGCGATTACTGAAAACAAAGCACAAATACGCAAGTTCTATTGAAACAAAGCGTCGGTGAACACA
         AGCGCGACTTGATTCGTTTTGTGCCTTAGCTACTGTAGAGTATAGGAGATAGCTTAATCCACTATTCAGCCAG
         CGCCTGTTAAGTATGGTTCAGCACAGGACGTTCCTCCTCGCGGCCTTG" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g13460_t1__Cre06g303200t11 =
    [|
        "CAGAGAAAATGCAAACGGAGTGCAGCAAAGCAATGGCCGGGCGACAGTGTGACAAGCAGGGGCGCAAACAGGG
         CGAGGCAAGGCAGAAGAGGGCGCCCGAATGTACGGATGAGTGGCGGGCGAGTCGAACACCGTGTGTGTGTGTG
         TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTTGGCGT
         GGAGCTACTGCCTGGCCGGGCGCCCGAAGCGTGAGGAATGCACATGGGCAAGCTACCAGCCTGCCAGCAAGTC
         AGTTATGTTGCGGTACAGTAACTCCGGCGTTATACCAACCGCGCAATCATACCTTGGGTTGCAAATTGAATGT
         TGGGGTCTGGGCACACGCACCCCCCGCACCACACGCCGACCACAGCTTCTGGTAGGGGTTGCACCGTTTACAT
         TGATACGCTTTGTGTTCAAGAGCGATCGCGGGGAGTGGGGGGAGGGCAGCGGAGAAGGGGGGGAGGCAAATAT
         TACGGTAACGGACCGAACGGGGCTGAGAGCGACACGCCGGTGCGTTGCGCTTTGCTTTCGGGAGCTGTTCAAC
         TCCCCTTTACGCGCTA" // 600 bases upstream of the 5'UTR 
        "GATCCTGCCGCTCGCCGCCTGTACGCTACCATATCGCCCTCAATACTCAGCGCACTGCACCACAAAATCAACC
         TGTCGCGTCAGTCACTCATCCGGAGCGCCGCGCCCGCCGCGTCGCCTTTCTGCGCCCTGGCTGCTGCAATATG
         CTATGAAAACATAGCCCTACCTCAAGAAAGCAGCGCCTGGACAACTTCATGCGCACTTCGCGTGACCTGCTGG
         CCCTGGAGCGATGGGCTGCATGGGGCGGAGCGACACTGACACTTGTCCCGTCCCGGGCTTCAGGCTCACTCAT
         ACCTCTGCTTACACATACAACAAGCGCTTACAAAGAACAATACGAAAACAACTCGCAACACCGCAATATACGT
         CTACATTGACGAGCGAGAACAACGAACACAAGTCGATTCTGTTCAGCACCGACGCGACCGCAGGTGCGCTGGG
         TTCCTCCGAGAGATCGTGCGCTTGTGCACAGGTTCTAGCGGGAGCTCTGACCGAGCGCGCACGCTGAAAACTT
         GTGCCGGCTGCTGACTCGAGCTGTTCATTCCCCGACCTTCGTTGTCCGCAGAAGCCCTTACCAACACCTTCCA
         GGCCTTCCATCTCTATCATACGTGAAGACTGGACCCCCTCGACTGCTGCTGGTGAGTCAGAGGCAGCTCGGGA
         GGCGCTTGCGTGCTTGGCGGAACATCACTGCAACCTGGCGTAAAAGCACTTGGCTCGCGTTGTTTGGGTGCCA
         GCAGAGGCCCATACATATGCTTAGGTTGCGCAATACAGCGTACTGCGGGGCTTCTTCGTGCAGGGCGTAGGGG
         ACTGACTTGAAATTCGTTTGCGAGGGCTAGGGTGAACGCTGGGTCATGCCGCTGGCAGCTCAAGCTATGAGTA
         GTATCAGTGATATGATGTACGCATTCTGGTCGGTGCGCCCTGCGAGATCGAACTCACATGCGGGGGCAACATG
         TGGTACCTCACACTTGGCTCAGCTTCAGGCCTGTTCTGCGCTCTCTCACACGCTTCCATGTGTGTTGATGCTA
         CTTGCAGGTCGACACTATCTGCGTGCCTATAGGGCAACGTCGAAGG" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g9842_t1__Cre11g467721t11 =
    [|
        "GCGTCGCGACTCCAGGAAAGTAGGCGGTATATCACCCAACAACCCCATGCAAGTCCGGTGTGGCATTGGCGCT
         GCGTGGCCAGCTCGGGCAAGCCGGGACACGGCTAGTTGCTCGATCGGCAAGCATTATGAATGCCACCACTTCA
         CCTGTCCCCTCACACCCGGCAAGCAAGTCGACCGCGCCTTGACTCACAATGGTTGCGATATTACAAAGAAGGG
         ACAATCCGACCAAGCAAGGGGACTCACAATGCTGCCCTACCAGGACGTCACGGCATGACCTGTCACGGCATGA
         CTTCTCACTTGTGACCTTTCGCGTTTCTGAATCGAGAAAGCTATAGGCATTTGATCAAGTATATCCATTTCTT
         TCACAACGGGCCCGTACTCGCCGAGATACAGCAATCGCGACCACGCTGCGCCGGGCTCGCTGCCCCCGCTCGC
         CTTCCCTCCGCCCCACCACTGCTGTCCGCCCGCGCCGGGTCCACGCACGGCGCACGCCGGGGCTTTCACGCCG
         TGACGTCCTGGTAGGACATGCCGTTTGTATGTGGCCATGGGTTCTTGCGCGACGTTTTGTGGTCCAGTTTCTC
         CTTGGTATCACATTTT" // 600 bases upstream of the 5'UTR 
        "GGTATCACGTTTGGTGCCCTGGGTAACAGTGGCTGTGACCTGCGACGGCAGTTCACGACAGCTGCGCTCGTTC
         GCAGTACCAGGCACAGCAGTTAGCCTGTGACGGGGCGAAGGCGAGCCGCGGGCCCGACTGCTTCCCACCTACA
         ACGACAGCTTTTGATAGCGATATTTGTAGAGCTTCTAAAGCACTACAAAGCCATACAGGTCTTTCGCGTGGCT
         CCACCACTTCGCCAACGCAACGGATCTAAGCGGCCATGCACGCAAGAGCCTTCGCACATAAGAAGTCGACTGG
         TCATGGTGTACTGTTATCTTCACACAGAATTCTACACATCTTAGCTCAGTCAAACACCAGACGAGTGCACGGA
         ACCAGGAGCTGAGTCTCAAAC" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let au5g1951_t1__Cre01g049132t12  =
    [|
        "TCGCCCAGAGCGCCCACAACAACTCCGCGTTCGCGCGCGCCGCGACCTCCAGCAGCTGCGGCACCAAGCACCGTGGG
         CACACGGGAGGGCCACCGGCGTTGCATGCTATGGCAACAATAGACAGGTTGGCGGGGTGGCGGGGTGGCGGGGAGCC
         GTGCGGCCGCATTACCCCTGGCTGCTGTCTCCAGGCAGGGTACAGGCGGGGCTGCATCAACTGATGGCTGCAGAGGC
         AGCCACGGCGGCTGCGGCGGTGGTGGGTGTGCTGGCGGTGGTGAGGGCAGCGGCAGCAGCAGGAGCGATGGTAGAAG
         GATGGGCACTGCGCCGGCGTTAGGGCACAAGCTGTGCGCCAACCCGCACTTCTCACTGAACAGCGCGTTGATAGCGC
         CCCGCTGCGGCAAGACATGTGTGCGGCGCCGTGCTGCAACATGCAGTGCCTGCCACTACTGGCGGACTGGCAGGGCG
         GCACTCGGGCATCTATAATGGATCCACCGCGCTGGCAACAGGGGCAGTGAAAGCAGTAGTCTAGGCATGAACCAACT
         GAGTTGGCGTACGTTGATACAGTCAGGATGGCTCTGCCCTGCTCAGTGCCCCGCTAAAGTT" // 600 bases upstream of the 5'UTR 
        "GCTACAGCGCTACCCACCGCAGCCCGGCAAGTTTCTGAGACCGCTACTGAGTTTACTAACAACACGCACCCTGGACC
         AGCTTCCTTAGCTGAGCAGGTGCCTGAGTCACTCCATTACTCCGGTTTCGTATCGATTTGTGGCGCTGCGACTTTAG
         AAGGGCGCTTTGAATAGCAAAGCAAGACAAAGTTTTGCAAGGGTCCGACGACCAGGGGGAAATGGATGCGCAATGCG
         CACAAGGTAAAGCCCTCTTGCTGTCGGGACCCGTTGTATTATATTAATATGCGATAACACACACATAGTCGAAAGCA
         ACACAAACAGTAATTCCATAAAAAAA" // 5'UTR
    |]
    |> Array.map (fun test -> BioArray.ofNucleotideString test)

let dataSet =
    [|au5g5407_t1__Cre14g617400t11; au5g11124_t1__Cre03g199150t12; au5g15219_t1__Cre09g387150t12; au5g15518_t1__Cre09g402304; au5g9850_t1__Cre01g071662t11; 
      au5g6062_t1__Cre16g650050t12; au5g784_t1__Cre10g457297t11; au5g13197_t1__Cre06g289900t11; au5g5992_t1__Cre03g198236t11; au5g8370_t1__Cre02g078226t11;
      au5g5129_t1__Cre13g603550t12; au5g11925_t1__Cre18g748547t11; au5g7722_t1__Cre17g733900t12; au5g7896_t1__Cre17g742550t12; au5g7422_t1__Cre17g718950t12; 
      au5g98_t1__Cre10g422300t12; au5g1550_t1__Cre01g028950t12; au5g1053_t1__Cre01g004157t12; au5g3913_t1__Cre12g542450t12; au5g15550_t1__Cre09g403850t11; 
      au5g3069_t1__Cre12g499850t11; au5g3976_t1__Cre12g545650t12; au5g13367_t1__Cre06g298600t12; au5g4913_t1__Cre13g592551t12; au5g5955_t1__Cre17g734564t11;
      au5g327_t1__Cre10g433950t11; au5g587_t1__Cre10g447300t12; au5g13365_t1__Cre06g298500t11; au5g13460_t1__Cre06g303200t11; au5g9842_t1__Cre11g467721t11; 
      au5g1951_t1__Cre01g049132t12
    |]
    |> Array.concat

let hsePPM = 
    trimmedHSEConsensus
    |> Array.map getPositionFrequencyMatrix
    |> fusePositionFrequencyMatrices trimmedHSEConsensus.[0].Length
    |> getPositionProbabilityMatrix (trimmedHSEConsensus.Length) dnaBases (0.0001)

//let hseMotifSampling = doMotifSamplingWithPPM 6 10 0.0001 7.5 dnaBases dataSet.[0..5] hsePPM
//    PSeq.init 10 (fun _ -> doMotifSamplingWithPPM 5 10 0.0001 10. dnaBases dataSet hsePPM)
//    |> Array.ofSeq

//hseMotifSampling
//|> Array.countBy (fun items -> items(* |> Array.map (fun item -> item.Positions)*))
//|> Array.sortByDescending (fun (item, amount) -> amount)


//Positions and probs at cutoff 1., motif amount 1 and motif length 10
let oneMotif =
    [|
        {PWMS = 9.834315053;    Positions = [517];};    {PWMS = 9.332569549;        Positions = [137];}; 
        {PWMS = 7.163690943;    Positions = [298];};    {PWMS = 7.790494176;        Positions = [32];}; 
        {PWMS = 8.755521567;    Positions = [199];};    {PWMS = 9.128686162;        Positions = [143];};
        {PWMS = 9.332569549;    Positions = [467];};    {PWMS = 3.433463156;        Positions = [25];}; 
        {PWMS = 10.4850627;     Positions = [25];};     {PWMS = 8.094380901;        Positions = [234];}; 
        {PWMS = 7.739695868;    Positions = [146];};    {PWMS = 3.956016214;        Positions = [58];};
        {PWMS = 10.4850627;     Positions = [88];};     {PWMS = 10.26109182;        Positions = [707];}; 
        {PWMS = 10.8815023;     Positions = [584];};    {PWMS = 6.164859137;        Positions = [83];}; 
        {PWMS = 7.224061807;    Positions = [464];};    {PWMS = 8.654040537;        Positions = [683];};
        {PWMS = 8.94935387;     Positions = [180];};    {PWMS = 6.197418139;        Positions = [75];}; 
        {PWMS = 8.331297843;    Positions = [238];};    {PWMS = 10.67190145;        Positions = [1242];}; 
        {PWMS = 7.459523032;    Positions = [150];};    {PWMS = 5.974569486;        Positions = [30];};
        {PWMS = 8.887008114;    Positions = [482];};    {PWMS = 12.43053872;        Positions = [254];}; 
        {PWMS = 8.642783719;    Positions = [451];};    {PWMS = 8.552930899;        Positions = [96];}; 
        {PWMS = 9.356206963;    Positions = [100];};    {PWMS = 6.405576351;        Positions = [283];};
        {PWMS = 5.219047088;    Positions = [407];};    {PWMS = 9.466219266;        Positions = [145];}; 
        {PWMS = 8.399229228;    Positions = [497];};    {PWMS = 7.128108484;        Positions = [91];}; 
        {PWMS = 9.396601381;    Positions = [361];};    {PWMS = 6.871814427;        Positions = [5];};
        {PWMS = 9.821240629;    Positions = [316];};    {PWMS = 7.410082936;        Positions = [333];}; 
        {PWMS = 8.521994798;    Positions = [118];};    {PWMS = 6.626232556;        Positions = [5];}; 
        {PWMS = 7.08006631;     Positions = [72];};     {PWMS = 6.064273508;        Positions = [49];};
        {PWMS = 13.084718;      Positions = [249];};    {PWMS = 6.942215207;        Positions = [37];}; 
        {PWMS = 10.13648574;    Positions = [104];};    {PWMS = 2.661333386e-06;    Positions = [];}; 
        {PWMS = 8.882868317;    Positions = [557];};    {PWMS = 11.01557213;        Positions = [25];};
        {PWMS = 8.195406124;    Positions = [582];};    {PWMS = 6.047626663;        Positions = [137];}; 
        {PWMS = 9.378198292;    Positions = [58];};     {PWMS = 9.086529868;        Positions = [265];}; 
        {PWMS = 8.174317492;    Positions = [510];};    {PWMS = 5.337514919;        Positions = [12];};
        {PWMS = 9.064274386;    Positions = [495];};    {PWMS = 7.352533408;        Positions = [13];}; 
        {PWMS = 8.188679756;    Positions = [279];};    {PWMS = 8.243604011;        Positions = [1012];}; 
        {PWMS = 7.940416792;    Positions = [194];};    {PWMS = 7.549793937;        Positions = [254];};
        {PWMS = 10.44684606;    Positions = [98];};     {PWMS = 8.482031534;        Positions = [96];}
    |]

//Positions and probs at cutoff 5., motif amount 2 and motif length 10
let twoMotifs =
    [|
        {PWMS = 38.3420492;     Positions = [306; 7];};     {PWMS = 38.72612873;    Positions = [137; 110];};
        {PWMS = 34.48027315;    Positions = [482; 302];};   {PWMS = 35.13794353;    Positions = [55; 32];};
        {PWMS = 37.89042613;    Positions = [199; 184];};   {PWMS = 36.64160328;    Positions = [345; 194];};
        {PWMS = 37.28132915;    Positions = [518; 467];};   {PWMS = 21.74727942;    Positions = [25; 14];};
        {PWMS = 38.96321105;    Positions = [378; 25];};    {PWMS = 31.79457991;    Positions = [252; 74];};
        {PWMS = 33.82950972;    Positions = [527; 231];};   {PWMS = 26.23483111;    Positions = [98; 0];};
        {PWMS = 37.99685343;    Positions = [171; 88];};    {PWMS = 37.89745549;    Positions = [707; 254];};
        {PWMS = 36.13021741;    Positions = [584; 458];};   {PWMS = 27.23040783;    Positions = [45; 23];};
        {PWMS = 28.24185721;    Positions = [150; 77];};    {PWMS = 34.90672327;    Positions = [346; 310];};
        {PWMS = 34.41779654;    Positions = [324; 64];};    {PWMS = 30.29268647;    Positions = [260; 233];};
        {PWMS = 35.84075551;    Positions = [474; 460];};   {PWMS = 39.33652709;    Positions = [361; 128];};
        {PWMS = 36.38487382;    Positions = [517; 150];};   {PWMS = 25.77972783;    Positions = [111; 30];};
        {PWMS = 37.38815516;    Positions = [482; 187];};   {PWMS = 40.14106608;    Positions = [254; 61];};
        {PWMS = 37.39167426;    Positions = [451; 206];};   {PWMS = 34.03930159;    Positions = [211; 133];};
        {PWMS = 37.27354871;    Positions = [215; 100];};   {PWMS = 29.13407088;    Positions = [235; 189];};
        {PWMS = 35.14878488;    Positions = [434; 407];};   {PWMS = 38.34965131;    Positions = [145; 6];};
        {PWMS = 32.50058922;    Positions = [540; 482];};   {PWMS = 31.21172822;    Positions = [125; 110];};
        {PWMS = 38.09768666;    Positions = [357; 153];};   {PWMS = 36.49169317;    Positions = [546; 362];};
        {PWMS = 34.55341014;    Positions = [554; 285];};   {PWMS = 32.15053505;    Positions = [182; 57];};
        {PWMS = 34.7203189;     Positions = [572; 118];};   {PWMS = 30.00006937;    Positions = [234; 5];};
        {PWMS = 36.82223061;    Positions = [230; 47];};    {PWMS = 34.29639831;    Positions = [107; 49];};
        {PWMS = 39.18861219;    Positions = [474; 249];};   {PWMS = 30.54903627;    Positions = [37; 0];};
        {PWMS = 35.5648808;     Positions = [444; 104];};   {PWMS = 17.92498088;    Positions = [17; 5];};
        {PWMS = 33.42122192;    Positions = [236; 193];};   {PWMS = 30.9687407;     Positions = [121; 87];};
        {PWMS = 37.50182668;    Positions = [582; 397];};   {PWMS = 34.80700288;    Positions = [228; 121];};
        {PWMS = 35.9227677;     Positions = [507; 95];};    {PWMS = 36.30228088;    Positions = [265; 7];};
        {PWMS = 35.81357658;    Positions = [542; 2];};     {PWMS = 27.83200135;    Positions = [32; 12];};
        {PWMS = 35.86964131;    Positions = [495; 268];};   {PWMS = 27.46551886;    Positions = [89; 13];};
        {PWMS = 34.76214864;    Positions = [451; 279];};   {PWMS = 34.95426201;    Positions = [783; 724];};
        {PWMS = 35.53701509;    Positions = [551; 46];};    {PWMS = 33.52076765;    Positions = [356; 254];};
        {PWMS = 36.26114606;    Positions = [430; 98];};    {PWMS = 36.07154359;    Positions = [188; 96];}
    |]

//Positions and probs at cutoff 7.5., motif amount 3 and motif length 10
let threeMotifs =
    [|
        {PWMS = 7.54984697;         Positions = [420; 354; 306];};  {PWMS = 11.11003929;        Positions = [278; 72; 36];};
        {PWMS = 4.110787234e-06;    Positions = [];};               {PWMS = 2.687551998e-06;    Positions = [];}; 
        {PWMS = 9.328505984;        Positions = [398; 199];};       {PWMS = 7.553335323;        Positions = [348; 135; 83];}; 
        {PWMS = 4.576850221e-06;    Positions = [];};               {PWMS = 1.5522762e-06;      Positions = [];}; 
        {PWMS = 5.557683911e-06;    Positions = [];};               {PWMS = 2.695665715e-06;    Positions = [];};
        {PWMS = 4.571262622e-06;    Positions = [];};               {PWMS = 3.347642323e-06;    Positions = [];}; 
        {PWMS = 8.512757513;        Positions = [585; 241; 171];};  {PWMS = 4.232012017e-06;    Positions = [];}; 
        {PWMS = 3.222352909e-06;    Positions = [];};               {PWMS = 2.185835695e-06;    Positions = [];};
        {PWMS = 3.733768213e-06;    Positions = [];};               {PWMS = 10.40337933;        Positions = [723; 683];}; 
        {PWMS = 3.774774615e-06;    Positions = [];};               {PWMS = 3.168722714e-06;    Positions = [];}; 
        {PWMS = 3.756109089e-06;    Positions = [];};               {PWMS = 10.34211951;        Positions = [1279; 1263; 1242];};
        {PWMS = 4.121594181e-06;    Positions = [];};               {PWMS = 2.950767625e-06;    Positions = [];}; 
        {PWMS = 8.81764848;         Positions = [562; 501; 482];};  {PWMS = 2.339212798e-06;    Positions = [];}; 
        {PWMS = 2.966564455e-06;    Positions = [];};               {PWMS = 2.923221664e-06;    Positions = [];};
        {PWMS = 10.21899815;        Positions = [545; 460; 424];};  {PWMS = 3.379880777e-06;    Positions = [];};
        {PWMS = 3.132022113e-06;    Positions = [];};               {PWMS = 1.716230651e-06;    Positions = [];}; 
        {PWMS = 3.728945693e-06;    Positions = [];};               {PWMS = 2.436722694e-06;    Positions = [];}; 
        {PWMS = 9.67350968;         Positions = [537; 372; 361];};  {PWMS = 3.21169744e-06;     Positions = [];}; 
        {PWMS = 4.924680911e-06;    Positions = [];};               {PWMS = 4.435206002e-06;    Positions = [];};
        {PWMS = 4.088074881e-06;    Positions = [];};               {PWMS = 7.509248176;        Positions = [251; 181; 5];};
        {PWMS = 3.69557523e-06;     Positions = [];};               {PWMS = 2.141965518e-06;    Positions = [];}; 
        {PWMS = 8.813923806;        Positions = [506; 287; 249];};  {PWMS = 1.954352506e-06;    Positions = [];}; 
        {PWMS = 11.54414681;        Positions = [247; 153; 104];};  {PWMS = 2.668117849e-06;    Positions = [];}; 
        {PWMS = 4.578733657e-06;    Positions = [];};               {PWMS = 3.588464567e-06;    Positions = [];};
        {PWMS = 5.661502247e-06;    Positions = [];};               {PWMS = 2.468978194e-06;    Positions = [];}; 
        {PWMS = 8.759057212;        Positions = [544; 453; 396];};  {PWMS = 9.931828228;        Positions = [403; 345; 265];}; 
        {PWMS = 8.058623128;        Positions = [582; 558; 542];};  {PWMS = 2.144150008e-06;    Positions = [];}; 
        {PWMS = 3.336132481e-06;    Positions = [];};               {PWMS = 2.186881887e-06;    Positions = [];};
        {PWMS = 9.141948659;        Positions = [515; 472; 451];};  {PWMS = 3.22463516e-06;     Positions = [];};
        {PWMS = 3.393814196e-06;    Positions = [];};               {PWMS = 3.350421745e-06;    Positions = [];}; 
        {PWMS = 3.685508625e-06;    Positions = [];};               {PWMS = 2.437912654e-06;    Positions = [];}
    |]

//Positions and probs at cutoff 10., motif amount 4 and motif length 10
let fourMotifs =
    [|
        {PWMS = 5.11782001e-06;     Positions = [];};                   {PWMS = 2.533653129e-06;    Positions = [];}; 
        {PWMS = 4.511019627e-06;    Positions = [];};                   {PWMS = 2.513731554e-06;    Positions = [];}; 
        {PWMS = 14.98035299;        Positions = [588; 282; 255; 199];}; {PWMS = 18.54194752;        Positions = [253; 229; 194; 83];}; 
        {PWMS = 5.205385669e-06;    Positions = [];};                   {PWMS = 1.749518213e-06;    Positions = [];}; 
        {PWMS = 6.116142372e-06;    Positions = [];};                   {PWMS = 2.506835815e-06;    Positions = [];};
        {PWMS = 5.171315959e-06;    Positions = [];};                   {PWMS = 3.652725442e-06;    Positions = [];}; 
        {PWMS = 5.319608623e-06;    Positions = [];};                   {PWMS = 5.352698287e-06;    Positions = [];}; 
        {PWMS = 3.505121901e-06;    Positions = [];};                   {PWMS = 2.099331634e-06;    Positions = [];};
        {PWMS = 4.174074363e-06;    Positions = [];};                   {PWMS = 4.106563364e-06;    Positions = [];}; 
        {PWMS = 4.465253381e-06;    Positions = [];};                   {PWMS = 2.373601757e-06;    Positions = [];}; 
        {PWMS = 4.346712864e-06;    Positions = [];};                   {PWMS = 5.480890088e-06;    Positions = [];};
        {PWMS = 4.583067287e-06;    Positions = [];};                   {PWMS = 2.585016579e-06;    Positions = [];}; 
        {PWMS = 2.629299288e-06;    Positions = [];};                   {PWMS = 2.666605361e-06;    Positions = [];}; 
        {PWMS = 3.184568966e-06;    Positions = [];};                   {PWMS = 3.257013534e-06;    Positions = [];};
        {PWMS = 6.170658764e-06;    Positions = [];};                   {PWMS = 3.824880691e-06;    Positions = [];}; 
        {PWMS = 3.981656891e-06;    Positions = [];};                   {PWMS = 2.160452815e-06;    Positions = [];}; 
        {PWMS = 4.154174189e-06;    Positions = [];};                   {PWMS = 3.026987893e-06;    Positions = [];};
        {PWMS = 2.930802229e-06;    Positions = [];};                   {PWMS = 2.756276242e-06;    Positions = [];}; 
        {PWMS = 4.846359166e-06;    Positions = [];};                   {PWMS = 4.259600577e-06;    Positions = [];}; 
        {PWMS = 4.350016565e-06;    Positions = [];};                   {PWMS = 2.357472741e-06;    Positions = [];};
        {PWMS = 3.920607354e-06;    Positions = [];};                   {PWMS = 2.601010935e-06;    Positions = [];}; 
        {PWMS = 3.292465851e-06;    Positions = [];};                   {PWMS = 2.555558432e-06;    Positions = [];}; 
        {PWMS = 3.374137512e-06;    Positions = [];};                   {PWMS = 2.261678438e-06;    Positions = [];};
        {PWMS = 5.227205027e-06;    Positions = [];};                   {PWMS = 3.061643894e-06;    Positions = [];}; 
        {PWMS = 6.870328612e-06;    Positions = [];};                   {PWMS = 2.877169248e-06;    Positions = [];}; 
        {PWMS = 3.509160128e-06;    Positions = [];};                   {PWMS = 3.662527476e-06;    Positions = [];};
        {PWMS = 2.877496604e-06;    Positions = [];};                   {PWMS = 2.147478103e-06;    Positions = [];}; 
        {PWMS = 3.496479818e-06;    Positions = [];};                   {PWMS = 2.598386594e-06;    Positions = [];}; 
        {PWMS = 4.538260996e-06;    Positions = [];};                   {PWMS = 3.64165553e-06;     Positions = [];};
        {PWMS = 3.886427377e-06;    Positions = [];};                   {PWMS = 3.607254433e-06;    Positions = [];}; 
        {PWMS = 4.816917586e-06;    Positions = [];};                   {PWMS = 2.463978637e-06;    Positions = [];}
    |]

//Positions and probs at cutoff 10., motif amount 5 and motif length 10
let fiveMotifs =
    [|
        {PWMS = 5.12623009e-06;     Positions = [];};                           {PWMS = 2.61065042e-06; Positions = [];}; 
        {PWMS = 4.621104947e-06;    Positions = [];};                           {PWMS = 2.59628359e-06; Positions = [];}; 
        {PWMS = 12.88455082;        Positions = [488; 388; 333; 220; 199];};    {PWMS = 13.37565129; Positions = [379; 334; 194; 105; 83];}; 
        {PWMS = 5.250899546e-06;    Positions = [];};                           {PWMS = 1.774215574e-06; Positions = [];}; 
        {PWMS = 6.121966285e-06;    Positions = [];};                           {PWMS = 2.577980148e-06; Positions = [];};
        {PWMS = 5.222753761e-06;    Positions = [];};                           {PWMS = 3.779722162e-06; Positions = [];}; 
        {PWMS = 5.345699349e-06;    Positions = [];};                           {PWMS = 5.318299404e-06; Positions = [];}; 
        {PWMS = 3.63190801e-06;     Positions = [];};                           {PWMS = 2.165388185e-06; Positions = [];};
        {PWMS = 4.289388228e-06;    Positions = [];};                           {PWMS = 4.229239143e-06; Positions = [];}; 
        {PWMS = 4.533566591e-06;    Positions = [];};                           {PWMS = 2.521390616e-06; Positions = [];}; 
        {PWMS = 4.437158622e-06;    Positions = [];};                           {PWMS = 5.478911554e-06; Positions = [];};
        {PWMS = 4.680329588e-06;    Positions = [];};                           {PWMS = 2.674659705e-06; Positions = [];}; 
        {PWMS = 2.735733775e-06;    Positions = [];};                           {PWMS = 2.779823234e-06; Positions = [];}; 
        {PWMS = 3.286499845e-06;    Positions = [];};                           {PWMS = 3.392904152e-06; Positions = [];};
        {PWMS = 6.103727148e-06;    Positions = [];};                           {PWMS = 3.956922858e-06; Positions = [];}; 
        {PWMS = 4.040683252e-06;    Positions = [];};                           {PWMS = 2.239937208e-06; Positions = [];}; 
        {PWMS = 4.273196028e-06;    Positions = [];};                           {PWMS = 3.126599954e-06; Positions = [];};
        {PWMS = 3.070887542e-06;    Positions = [];};                           {PWMS = 2.889426819e-06; Positions = [];}; 
        {PWMS = 5.005014391e-06;    Positions = [];};                           {PWMS = 4.464621966e-06; Positions = [];}; 
        {PWMS = 4.484249048e-06;    Positions = [];};                           {PWMS = 2.435479719e-06; Positions = [];};
        {PWMS = 4.074698903e-06;    Positions = [];};                           {PWMS = 2.698648859e-06; Positions = [];}; 
        {PWMS = 3.377482024e-06;    Positions = [];};                           {PWMS = 2.632114296e-06; Positions = [];}; 
        {PWMS = 3.520181828e-06;    Positions = [];};                           {PWMS = 2.321699946e-06; Positions = [];};
        {PWMS = 5.269654163e-06;    Positions = [];};                           {PWMS = 3.182942514e-06; Positions = [];}; 
        {PWMS = 6.736366183e-06;    Positions = [];};                           {PWMS = 2.964365251e-06; Positions = [];}; 
        {PWMS = 3.681833477e-06;    Positions = [];};                           {PWMS = 3.859096609e-06; Positions = [];};
        {PWMS = 3.02484243e-06;     Positions = [];};                           {PWMS = 2.195609271e-06; Positions = [];}; 
        {PWMS = 3.670819555e-06;    Positions = [];};                           {PWMS = 2.733926608e-06; Positions = [];}; 
        {PWMS = 4.42957831e-06;     Positions = [];};                           {PWMS = 3.746268038e-06; Positions = [];};
        {PWMS = 4.003215622e-06;    Positions = [];};                           {PWMS = 3.772126083e-06; Positions = [];}; 
        {PWMS = 4.765353998e-06;    Positions = [];};                           {PWMS = 2.580523595e-06; Positions = [];}
    |]

let expressionRates =
    [|
        au5g5407_t1__Cre14g617400t11,   8.530948,   10
        au5g11124_t1__Cre03g199150t12,  -0.8772487, 3
        au5g15219_t1__Cre09g387150t12,  4.332184,   12
        au5g15518_t1__Cre09g402304,     -0.8688895, 3
        au5g9850_t1__Cre01g071662t11,   -1.853293,  5
        au5g6062_t1__Cre16g650050t12,   2.702928,   6
        au5g784_t1__Cre10g457297t11,    0.7954394,  6
        au5g13197_t1__Cre06g289900t11,  2.652841,   5
        au5g5992_t1__Cre03g198236t11,   -3.708641,  7
        au5g8370_t1__Cre02g078226t11,   0.9523705,  5
        au5g5129_t1__Cre13g603550t12,   -2.943189,  7
        au5g11925_t1__Cre18g748547t11,  -0.8265074, 4
        au5g7722_t1__Cre17g733900t12,   -3.232005,  6
        au5g7896_t1__Cre17g742550t12,   -1.475992,  5
        au5g7422_t1__Cre17g718950t12,   1.947759,   8
        au5g98_t1__Cre10g422300t12,     5.837328,   4
        au5g1550_t1__Cre01g028950t12,   -3.549863,  5
        au5g1053_t1__Cre01g004157t12,   -3.256696,  5
        au5g3913_t1__Cre12g542450t12,   -0.9271605, 6
        au5g15550_t1__Cre09g403850t11,  0.3335226,  6
        au5g3069_t1__Cre12g499850t11,   0.5008245,  5
        au5g3976_t1__Cre12g545650t12,   -0.8303188, 6
        au5g13367_t1__Cre06g298600t12,  3.703181,   5
        au5g4913_t1__Cre13g592551t12,   2.791368,   6
        au5g5955_t1__Cre17g734564t11,   -1.819351,  4
        au5g327_t1__Cre10g433950t11,    -0.05315473,10
        au5g587_t1__Cre10g447300t12,    -0.5116726, 5
        au5g13365_t1__Cre06g298500t11,  1.687511,   4
        au5g13460_t1__Cre06g303200t11,  1.152718,   7
        au5g9842_t1__Cre11g467721t11,   -3.74618,   5
        au5g1951_t1__Cre01g049132t12,   0.7579265,  5
    |]

let expressionRatesOnlyPromoters =
    [|
        au5g5407_t1__Cre14g617400t11,   8.530948,   5
        au5g11124_t1__Cre03g199150t12,  -0.8772487, 2
        au5g15219_t1__Cre09g387150t12,  4.332184,   8
        au5g15518_t1__Cre09g402304,     -0.8688895, 2
        au5g9850_t1__Cre01g071662t11,   -1.853293,  2
        au5g6062_t1__Cre16g650050t12,   2.702928,   3
        au5g784_t1__Cre10g457297t11,    0.7954394,  4
        au5g13197_t1__Cre06g289900t11,  2.652841,   2
        au5g5992_t1__Cre03g198236t11,   -3.708641,  3
        au5g8370_t1__Cre02g078226t11,   0.9523705,  3
        au5g5129_t1__Cre13g603550t12,   -2.943189,  2
        au5g11925_t1__Cre18g748547t11,  -0.8265074, 2
        au5g7722_t1__Cre17g733900t12,   -3.232005,  4
        au5g7896_t1__Cre17g742550t12,   -1.475992,  2
        au5g7422_t1__Cre17g718950t12,   1.947759,   5
        au5g98_t1__Cre10g422300t12,     5.837328,   2
        au5g1550_t1__Cre01g028950t12,   -3.549863,  3
        au5g1053_t1__Cre01g004157t12,   -3.256696,  3
        au5g3913_t1__Cre12g542450t12,   -0.9271605, 3
        au5g15550_t1__Cre09g403850t11,  0.3335226,  2
        au5g3069_t1__Cre12g499850t11,   0.5008245,  2
        au5g3976_t1__Cre12g545650t12,   -0.8303188, 4
        au5g13367_t1__Cre06g298600t12,  3.703181,   4
        au5g4913_t1__Cre13g592551t12,   2.791368,   3
        au5g5955_t1__Cre17g734564t11,   -1.819351,  2
        au5g327_t1__Cre10g433950t11,    -0.05315473,6
        au5g587_t1__Cre10g447300t12,    -0.5116726, 5
        au5g13365_t1__Cre06g298500t11,  1.687511,   2
        au5g13460_t1__Cre06g303200t11,  1.152718,   5
        au5g9842_t1__Cre11g467721t11,   -3.74618,   3
        au5g1951_t1__Cre01g049132t12,   0.7579265,  2
    |]

expressionRatesOnlyPromoters
|> Array.map (fun (_, value, motifCount) -> value/float motifCount)

let getLength (source:MotifIndex[]) =
    let rec loop n acc =
        if n >= source.Length then List.rev acc
        else
            loop (n+2) (source.[n].Positions.Length + source.[n+1].Positions.Length :: acc)
    loop 0 []

let names =
    [|
        "au5g5407_t1__Cre14g617400t11"
        "au5g11124_t1__Cre03g199150t12"
        "au5g15219_t1__Cre09g387150t12"
        "au5g15518_t1__Cre09g402304"
        "au5g9850_t1__Cre01g071662t11"
        "au5g6062_t1__Cre16g650050t12"
        "au5g784_t1__Cre10g457297t11 "
        "au5g13197_t1__Cre06g289900t11"
        "au5g5992_t1__Cre03g198236t11"
        "au5g8370_t1__Cre02g078226t11"
        "au5g5129_t1__Cre13g603550t12"
        "au5g11925_t1__Cre18g748547t11"
        "au5g7722_t1__Cre17g733900t12"
        "au5g7896_t1__Cre17g742550t12"
        "au5g7422_t1__Cre17g718950t12"
        "au5g98_t1__Cre10g422300t12"
        "au5g1550_t1__Cre01g028950t12"
        "au5g1053_t1__Cre01g004157t12"
        "au5g3913_t1__Cre12g542450t12"
        "au5g15550_t1__Cre09g403850t11"
        "au5g3069_t1__Cre12g499850t11"
        "au5g3976_t1__Cre12g545650t12"
        "au5g13367_t1__Cre06g298600t12"
        "au5g4913_t1__Cre13g592551t12"
        "au5g5955_t1__Cre17g734564t11"
        "au5g327_t1__Cre10g433950t11 "
        "au5g587_t1__Cre10g447300t12 "
        "au5g13365_t1__Cre06g298500t11"
        "au5g13460_t1__Cre06g303200t11"
        "au5g9842_t1__Cre11g467721t11"
        "au5g1951_t1__Cre01g049132t12"
    |]

//let motifAmount =
//    expressionRates
//    |> Array.map (fun (_, _, x) -> float x)

//let expressionRate =
//    expressionRates
//    |> Array.map (fun (_, y, _) -> (y))

//let normalizedExpressionRate =
//    expressionRates
//    |> Array.map (fun (_, y, x) -> (y/float x))

//(fun x y -> PearsonCorrelation().Similarity(x, y)) motifAmount expressionRate

//let values = expressionRate
//let keys   = names
//let labels = motifAmount

//Chart.Column(keys,values,Labels=labels,Opacity=0.3,Marker=Marker.init(Color="rgba(222,45,38,0.8)",Size=1)) // Changing the thickness of the bar is not possible at the moment
//|> Chart.Show