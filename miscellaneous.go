//This file contains small functions to help with the main function.

package main

import (
	"fmt"
	"math"
)

//rounding number to x decimal place
func RoundTo(num float64, x int) float64 {
	deciPlace := math.Pow10(x)
	return math.Round(num*deciPlace) / deciPlace
}

//LogLikeliHood returns the LogLikeliHood of a and b. If it's bigger than1,
//then it's more likely to be a than b. Otherwise, it's more likely to be b.
func LogLikeliHood(a, b float64) float64 {
	loglikelihood := math.Log(a / b)
	return loglikelihood
}

//SumOfIntSlice takse a integer slice and return the sum
func SumOfIntSlice(s []int) int {
	sum := 0
	for _, i := range s {
		sum += i
	}
	return sum
}

//MakeMapHeader takes an interger, returns the header of the transition map or emission map
//It's used for the ProfileHMM and creating an empty map
func MakeMapHeader(numAlignment int) []string {
	header := make([]string, (numAlignment+1)*3)
	header[0] = "Start"
	header[1] = "I0"
	header[len(header)-1] = "End"

	for i := 0; i < numAlignment; i++ {
		header[3*i+2] = "M" + fmt.Sprintf("%d", i+1)
		header[3*i+3] = "D" + fmt.Sprintf("%d", i+1)
		header[3*i+4] = "I" + fmt.Sprintf("%d", i+1)
	}

	return header
}

//CreatEmptyMap takes the header of the row and column, creates an empty map of float64
func CreatEmptyMap(rowheader, colheader []string) MtxMap {
	emptyMap := make(map[string]map[string]float64)
	var eachEntry map[string]float64

	for _, r := range rowheader {
		eachEntry = make(map[string]float64)
		for _, c := range colheader {
			eachEntry[c] = 0.0
		}
		emptyMap[r] = eachEntry
	}
	return emptyMap
}

//MapToMtx takes the row header, colheader and the map, turns the map to a matrix
//according to the sequence of the rowheader and colheader.
//Output: rowheader, colheader and the matix aligned with sequences of the headers.
func MapToMtx(rowheader, colheader []string, theMap MtxMap) ([]string, []string, [][]float64) {
	theMtx := make([][]float64, len(theMap))
	
	for r := range rowheader {
		theMtx[r] = make([]float64, len(colheader))
		for c := range colheader {
			theMtx[r][c] = theMap[rowheader[r]][colheader[c]]
		}
	}
	return rowheader, colheader, theMtx
}

//MarkDeletionState takes in the theta value, a series of multiAlign (need to be at same length)
//and then returns a slice of indexes of 0 and 1 with 0 beding the deletion state.
func MarkDeletionState(theta float64, multiAlign []string) []int {
	size := len(multiAlign[0])
	statesInFloat := make([]float64, size)
	statesInInt := make([]int, size)

	for _, string := range multiAlign {
		for a := 0; a < size; a++ {
			if string[a:a+1] == "-" {
				statesInFloat[a] += 1.0 / float64(len(multiAlign))
			}
		}
	}
	for d := range statesInFloat {
		if statesInFloat[d] >= theta {
			statesInInt[d] = 0
		} else {
			statesInInt[d] = 1
		}
	}
	return statesInInt
}
