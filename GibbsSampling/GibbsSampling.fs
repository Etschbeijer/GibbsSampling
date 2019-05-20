namespace GibbsSampling


open System
open FSharpAux
open FSharpAux.IO
open BioFSharp


module HelperFunctionsAndTypes =

    /// One dimensional array with fixed positions for each element.
    type BaseArray<'a, 'value when 'a :> IBioItem>internal (array:'value []) =
    
        let getIndex (key:'a) =
            (BioItem.symbol key |> int) - 42

        new () =
            let arr:'value [] = Array.zeroCreate 49
            new BaseArray<_, 'value>(arr)

        member this.Array = 
            if (obj.ReferenceEquals(array, null)) then
                raise (ArgumentNullException("array"))
            array
            
        member this.Item
            with get i       = array.[getIndex i]
            and  set i value = array.[getIndex i] <- value

    /// Matrix with fixed positions for nucleotides and amino acids.
    type BaseMatrix<'a, 'value when 'a :> IBioItem>internal(matrix:'value [,]) =

        let getRowArray2DIndex (key:'a) =
            (BioItem.symbol key |> int) - 42

        new (rowLength:int) =
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new BaseMatrix<_, 'value>(arr)

        member this.Matrix = 
            if (obj.ReferenceEquals(matrix, null)) then
                raise (ArgumentNullException("array"))
            matrix
            
        member this.Item
            with get (column, row)       = matrix.[getRowArray2DIndex column, row]
            and  set (column, row) value = matrix.[getRowArray2DIndex column, row] <- value

    let fromFileObo (filePath:string) =
        SeqIO.Seq.fromFileWithSep '\t' filePath
        |> Array.ofSeq
        |> Array.collect (fun item -> Array.tail item)

    ///Checks whether all elements in the list have a wider distance than width or not.
    let ceckForDistance (width:int) (items:int list) =
        if items.Length <= 0 || items.Length = 1 then true
        else
            let rec loop n i =
                if n = items.Length-1 then true
                else
                    if i >= items.Length then loop (n+1) (n+2)
                    else
                        if Operators.abs(items.[n]-items.[i]) > width then
                            loop n (i+1)
                        else false
            loop 0 1

    /// Get an integer which is between 0 and the length of the sequence - segmentLength
    let getRandomNumberInSequence (segmentLength:int) (source:'a[]) =
        let rnd = System.Random()
        Array.init 1 (fun _ -> rnd.Next(0, source.Length-segmentLength+1))
        |> Array.head

    /// Create a specific sub sequence of the source sequence based on the given length and starting position. Do not forget, sequences start counting at 0!
    let getDefinedSegment (subsequenceLength:int) (source:'a[]) (startPoint:int) =
        source 
        |> Array.skip startPoint 
        |> Array.take subsequenceLength 
        |> Array.ofSeq

    ///Finds the best information content in an array of arrays of PWMSs and Positions.
    let getBestInformationContent (item:((float*int)[])[]) =
        let rec loop (n:int) (bestPWMS:(float*int)[]) =
            if n = item.Length then bestPWMS
            else
                let informationContentItem =
                    item.[n]
                    |> Array.fold (fun baseValue (pwms, _) -> pwms + baseValue) 0.
                let informationContentBestPWMS =
                    bestPWMS
                    |> Array.fold (fun baseValue (pwms, _) -> pwms + baseValue) 0.
                if informationContentItem > informationContentBestPWMS then
                    loop (n + 1) item.[n]
                else
                    loop (n + 1) bestPWMS
        loop 0 [||]

open HelperFunctionsAndTypes

module FrequencyCompositeVector =

    /// One dimensional array with fixed positions for each element.
    /// Use to track frequency of elements independent of position in source.
    type FrequencyCompositeVector internal (array:int []) =

        inherit BaseArray<IBioItem, int>(array)

        new() = 
            let arr:'value [] = Array.zeroCreate 49
            new FrequencyCompositeVector(arr)

    /// Increase counter of position by 1.
    let increaseInPlace (bioItem:'a) (frequencyCompositeVector:FrequencyCompositeVector) = 
        frequencyCompositeVector.[bioItem] <- frequencyCompositeVector.[bioItem] + 1
        frequencyCompositeVector

    /// Increase counter of position by n.
    let increaseInPlaceBy (bioItem:'a) (frequencyCompositeVector:FrequencyCompositeVector) n = 
        frequencyCompositeVector.[bioItem] <- frequencyCompositeVector.[bioItem] + n
        frequencyCompositeVector

    /// Create a FrequencyCompositeVector based on BioArrays and exclude the specified segments.
    let getFrequencyVectorOfSources (motifLength:int) (resSources:BioArray.BioArray<#IBioItem>[]) (motifPositions:int[]) =
        let backGroundCounts = new FrequencyCompositeVector()    
        Array.map2 (fun (resSource:BioArray.BioArray<#IBioItem>) (position:int) -> Array.append resSource.[0..(position-1)] resSource.[(position+motifLength)..]) resSources motifPositions
        |> Array.concat
        |> Array.fold (fun bc bioItem -> (increaseInPlace bioItem bc)) backGroundCounts

    /// Create a FrequencyCompositeVector based on BioArrays and exclude the specified segments.
    let getFrequencyVectorOfSource (resSources:BioArray.BioArray<#IBioItem>) =
        let backGroundCounts = new FrequencyCompositeVector()   
        Array.fold (fun bc bioItem -> (increaseInPlace bioItem bc)) backGroundCounts resSources

    /// Create new FrequencyCompositeVector based on the sum of the positions of an array of FrequencyVectors.
    let fuseFrequencyVectors (alphabet:#IBioItem[]) (bfVectors:FrequencyCompositeVector[]) =
        let backgroundFrequencyVector = new FrequencyCompositeVector()
        for bfVector in bfVectors do
            for bioItem in alphabet do
                backgroundFrequencyVector.[bioItem] <- backgroundFrequencyVector.[bioItem] + (bfVector.[bioItem])
        backgroundFrequencyVector

    /// Create FrequencyCompositeVector based on BioArray and excludes the specified segment.
    let getFrequencyVectorWithoutSegment (motifLength:int) (position:int) (resSource:BioArray.BioArray<#IBioItem>) =
        let backGroundCounts = new FrequencyCompositeVector()
        Array.append resSource.[0..(position-1)] resSource.[(position+motifLength)..]
        |> Array.fold (fun bc bioItem -> (increaseInPlace bioItem bc)) backGroundCounts

    /// Create FrequencyCompositeVector based on BioArray.
    let addInPlaceCountsOfSource (resSources:BioArray.BioArray<#IBioItem>) (backGroundCounts:FrequencyCompositeVector) =   
        resSources
        |> Array.fold (fun bc bioItem -> (increaseInPlace bioItem bc)) backGroundCounts

    /// Subtracts the amount of elements in the given source from the FrequencyCompositeVector.
    let substractCountsOfBFVector (source:BioArray.BioArray<#IBioItem>) (bfVector:FrequencyCompositeVector) =
        let bfVec = new FrequencyCompositeVector(bfVector.Array)
        for bioItem in source do
            bfVec.[bioItem] <- (if bfVector.[bioItem] - 1 > 0 then bfVector.[bioItem] - 1 else 0)
        bfVec

/// Includes tools to create and manipulate ProbabilityVectors.
module ProbabilityCompositeVector = 
    
    type ProbabilityCompositeVector internal (array:float []) =
        
        inherit BaseArray<IBioItem,float>(array)

        new() = 
            let arr:'value [] = Array.zeroCreate 49
            new ProbabilityCompositeVector(arr)

    /// Increase counter of position by 1.
    let increaseInPlace (bioItem:'a) (backGroundProbabilityVector:ProbabilityCompositeVector) = 
        backGroundProbabilityVector.[bioItem] <- backGroundProbabilityVector.[bioItem] + 1.
        backGroundProbabilityVector

    /// Increase counter of position by n.
    let increaseInPlaceBy (bioItem:'a) (backGroundProbabilityVector:ProbabilityCompositeVector) n = 
        backGroundProbabilityVector.[bioItem] <- backGroundProbabilityVector.[bioItem] + n
        backGroundProbabilityVector

    /// Create a ProbabilityCompositeVector based on a FrequencyCompositeVector by replacing the integers with floats.
    let getBPVofBFV (caArray:FrequencyCompositeVector.FrequencyCompositeVector) =
        caArray.Array
        |> Array.map (fun item -> float item)
        |> fun item -> new ProbabilityCompositeVector(item)

    /// Create normalized ProbabilityCompositeVector based on FrequencyCompositeVector.
    let getNormalizedProbabilityVector (alphabet:#IBioItem[]) (pseudoCount:float) (frequencyCompositeVector:FrequencyCompositeVector.FrequencyCompositeVector) =
        let backGroundProbabilityVector = getBPVofBFV frequencyCompositeVector
        let sum = (float (Array.sum frequencyCompositeVector.Array)) + ((float alphabet.Length) * pseudoCount)
        for item in alphabet do
            backGroundProbabilityVector.[item] <- (backGroundProbabilityVector.[item] + pseudoCount)/sum
        backGroundProbabilityVector

    /// Calculate the score of the given sequence based on the Probabilities of the ProbabilityCompositeVector.
    let calculateBackGroundSegmentScore (pcv:ProbabilityCompositeVector) (bioItems:BioArray.BioArray<#IBioItem>) =
        Array.fold (fun (value:float) (bios:#IBioItem) -> value * (pcv.[bios])) 1. bioItems

/// Includes tools to create count-matrices in order to create Position-Specific-PWMS-Matrix.
module PositionFrequencyMatrix =

    /// Matrix with fixed positions for nucleotides and amino acids with default value of 0.
    type PositionFrequencyMatrix internal(matrix:int [,]) =

        inherit BaseMatrix<IBioItem, int>(matrix)

        new (rowLength:int) = 
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new PositionFrequencyMatrix(arr)

    /// Increase counter of PositionFrequencyMatrix at fixed position by 1.
    let increaseInPlace (pos:int) (bioItem:'a when 'a :> IBioItem) (positionFrequencyMatrix:PositionFrequencyMatrix) = 
        positionFrequencyMatrix.[bioItem, pos] <- positionFrequencyMatrix.[bioItem, pos] + 1
        positionFrequencyMatrix

    /// Increase counter of PositionFrequencyMatrix at fixed position by n.
    let increaseInPlaceBy (pos:int) (bioItem:'a when 'a :> IBioItem) (positionFrequencyMatrix:PositionFrequencyMatrix) n = 
        positionFrequencyMatrix.[bioItem, pos] <- positionFrequencyMatrix.[bioItem, pos] + n
        positionFrequencyMatrix

    /// Create PositionFrequencyMatrix based on BioArray.
    let getPositionFrequencyMatrix (source:BioArray.BioArray<#IBioItem>) =
        let positionFrequencyMatrix = new PositionFrequencyMatrix(source.Length)
        source
        |> Array.fold (fun (row, cm) column -> row + 1, increaseInPlace row column cm) (0, positionFrequencyMatrix) |> ignore
        positionFrequencyMatrix

    /// Create new PositionFrequencyMatrix based on the sum of the positions of an array of Countmatrices.
    let fusePositionFrequencyMatrices (motifLength:int) (countMatrices:PositionFrequencyMatrix[]) =
        let positionFrequencyMatrix = new PositionFrequencyMatrix(motifLength)
        if countMatrices.Length > 0 then
            for cMatrix in countMatrices do
                for column = 0 to (Array2D.length1 cMatrix.Matrix)-1 do
                        for row = 0 to (Array2D.length2 cMatrix.Matrix)-1 do
                            positionFrequencyMatrix.Matrix.[column, row] <- positionFrequencyMatrix.Matrix.[column, row] + (cMatrix.Matrix.[column, row])
            positionFrequencyMatrix
        else positionFrequencyMatrix

/// Includes tools to work with Weight-Matrices in order to create a Position-Specific-PWMS-Matrix.
module PositionProbabilityMatrix = 
    
    /// Matrix with fixed positions for nucleotides and amino acids with default value of 0. .
    type PositionProbabilityMatrix internal(matrix:float [,]) =

        inherit BaseMatrix<IBioItem, float>(matrix)

        new(rowLength:int) = 
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new PositionProbabilityMatrix(arr)

    /// Increase counter of PositionProbabilityMatrix at fixed position by 1.
    let increaseInPlace (pos:int) (bioItem:'a when 'a :> IBioItem) (positionProbabilityMatrix:PositionProbabilityMatrix) = 
        positionProbabilityMatrix.[bioItem, pos] <- positionProbabilityMatrix.[bioItem, pos] + 1.
        positionProbabilityMatrix

    /// Increase counter of PositionProbabilityMatrix at fixed position by n.
    let increaseInPlaceBy (pos:int) (bioItem:'a when 'a :> IBioItem) (positionProbabilityMatrix:PositionProbabilityMatrix) n = 
        positionProbabilityMatrix.[bioItem, pos] <- positionProbabilityMatrix.[bioItem, pos] + n
        positionProbabilityMatrix

    /// Create new PositionWeightMatrix based on existing PositionFrequencyMatrix. The counts of each position of each element are transformed to floats.
    let transformCountMatrixToPPMatrix (positionFrequencyMatrix:PositionFrequencyMatrix.PositionFrequencyMatrix) =
        positionFrequencyMatrix.Matrix |> Array2D.map (fun item -> float item)
        |> fun item -> new PositionProbabilityMatrix(item)

    /// Create PositionProbabilityMatrix based on PositionFrequencyMatrix.
    let getPositionProbabilityMatrix (sourceCount:int) (alphabet:#IBioItem[]) (pseudoCount:float) (positionFrequencyMatrix:PositionFrequencyMatrix.PositionFrequencyMatrix) =
        let positionProbabilityMatrix = transformCountMatrixToPPMatrix positionFrequencyMatrix
        let sum = (float sourceCount) + ((float alphabet.Length) * pseudoCount)
        for item in alphabet do
            for position = 0 to (Array2D.length2 positionProbabilityMatrix.Matrix) - 1 do
            positionProbabilityMatrix.[item, position] <- (positionProbabilityMatrix.[item, position] + pseudoCount)/sum
        positionProbabilityMatrix

        /// Includes tools to work with Weight-Matrices in order to create a Position-Specific-PWMS-Matrix.
module PositionWeightMatrix = 
    
    type PositionWeightMatrix internal(matrix:float [,]) =

        inherit BaseMatrix<IBioItem, float>(matrix)

        new(rowLength:int) = 
            let arr:'value [,] = Array2D.zeroCreate 49 rowLength 
            new PositionWeightMatrix(arr)

    /// Increase counter of PositionWeightMatrix at fixed position by 1.
    let increaseInPlace (pos:int) (bioItem:'a when 'a :> IBioItem) (positionWeightMatrix:PositionWeightMatrix) = 
        positionWeightMatrix.[bioItem, pos] <- positionWeightMatrix.[bioItem, pos] + 1.
        positionWeightMatrix

    // Increase counter of PositionWeightMatrix at fixed position by n.
    let increaseInPlaceBy (pos:int) (bioItem:'a when 'a :> IBioItem) (positionWeightMatrix:PositionWeightMatrix) n = 
        positionWeightMatrix.[bioItem, pos] <- positionWeightMatrix.[bioItem, pos] + n
        positionWeightMatrix

    /// Create PositionWeightMatrix based on ProbabilityCompositeVector and PositionProbabilityMatrix.
    let createPositionWeightMatrix (alphabet:#IBioItem[]) (pcv:ProbabilityCompositeVector.ProbabilityCompositeVector) (ppMatrix:PositionProbabilityMatrix.PositionProbabilityMatrix) =
        let pwMatrix = new PositionWeightMatrix(Array2D.length2 ppMatrix.Matrix)
        for item in alphabet do
            for position=0 to (Array2D.length2 ppMatrix.Matrix)-1 do
                pwMatrix.[item, position] <- ppMatrix.[item, position]/pcv.[item]
        pwMatrix

    /// Calculate the score of the given sequence based on the Probabilities of the PositionWeightMatrix.
    let calculateSegmentScore (pwMatrix:PositionWeightMatrix) (bioItems:BioArray.BioArray<#IBioItem>) =
        Array.fold (fun (position:int, value:float) (bios:#IBioItem) -> 
            position + 1, value * (pwMatrix.[bios, position])) (0, 1.) bioItems
        |> snd
