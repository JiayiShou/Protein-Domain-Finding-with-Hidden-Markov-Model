//This file contains all the functions needed for viterbi decoding.
package main

import (
	"math"
)

type coordinate struct {
	x int
	y int
}

//Input: A string x, Σ, States, transition map, and emission map of an HMM (Σ, States, Transition, Emission).
//The start point is a integer indicating whether there is invisible states at the
//front and end. If it's 1, then it starts from invisible state "Start" and end in "End".
//Otherwise each state have equal opportunity as a starting state.
//output: A path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π.
func ViterbiDecoding(startpoint int, str string, sigma, states []string, trmap, emimap MtxMap) string {
	if len(str) == 0 {
		panic("Can't create a path for string length 0.")
	}
	var viterbi [][]float64      //len(states) x len(str)
	var backtrace [][]coordinate //len(states) x len(str)-1

	if startpoint == 1 { // it starts from the invisible state "Start"
		viterbi, backtrace = FillViterbiWithHiddenStates(str, sigma, states, trmap, emimap)
	} else {
		viterbi, backtrace = FillViterbiNoHiddenStates(str, sigma, states, trmap, emimap)
	}

	endingState := EndingInViterbi(startpoint, str, states, viterbi)
	path := ViterbiTraceBack(endingState, str, states, backtrace) //Finding the path trace back from the ending state.
	return path
}

//If some states in the HMM have no emission (for example, deletion state does
//not emit anything), Then use this function to fill up the viterbi matrix. It
//takes the input string, sigma(emission options) states(for transition), transition
//map and emision map, and return a filled viterbi matrix and a matrix to trace back the path.
func FillViterbiWithHiddenStates(str string, sigma, states []string, trmap, emimap MtxMap) ([][]float64, [][]coordinate) {
	var logTranEmi, max float64
	viterbi := make([][]float64, len(states))      //len(states) x len(str)
	backtrace := make([][]coordinate, len(states)) //len(states) x len(str)-1

	for v1 := range viterbi { //Make the first volum of viterbi and backtrace...
		viterbi[v1] = make([]float64, len(str)+1)
		backtrace[v1] = make([]coordinate, len(str)) //-1
		if v1 == 0 {
			viterbi[0][0] = 1
		} else if states[v1][0:1] == "D" { // if it's a deletions state, its the previous deletion state plus the transition probability.
			viterbi[v1][0] = viterbi[v1-3][0] + math.Log(trmap[states[v1-3]][states[v1]])
		}
	}
	
	//fill in up the viterbi matrix row by row.
	for s := 1; s <= len(str); s++ {
		for v2 := range viterbi { //each entry in v2 is the max{each v1 to v2}
			//if v2 is a deletion state, we travel along column.
			if states[v2][0:1] == "D" || states[v2][0:1] == "E" {
				max = viterbi[0][s] + math.Log(trmap[states[0]][states[v2]])
				for v1 := range viterbi[0:v2] { //v1 states before v2, in viterbi at the same column as v2
					if trmap[states[v1]][states[v2]] != 0 {
						logTranEmi = viterbi[v1][s] + math.Log(trmap[states[v1]][states[v2]]) //no emission for hidden states
					}
					if logTranEmi >= max {
						max = logTranEmi
						backtrace[v2][s-1].x = v1
						backtrace[v2][s-1].y = s - 1
					}
				}
			} else {
				max = viterbi[0][s-1] + math.Log(trmap[states[0]][states[v2]]) + math.Log(emimap[states[v2]][str[s-1:s]])
				for v1 := range viterbi {
					logTranEmi = viterbi[v1][s-1] + math.Log(trmap[states[v1]][states[v2]]) + math.Log(emimap[states[v2]][str[s-1:s]])
					if logTranEmi >= max {
						max = logTranEmi
						backtrace[v2][s-1].x = v1
						backtrace[v2][s-1].y = s - 2
					}
				}
			}
			viterbi[v2][s] = max
		}
	}
	return viterbi, backtrace
}

//If all states in the HMM have emission (for example, there is no deletion), Then
//use this function to fill up the viterbi matrix. It takes the input string, sigma(emission options)
//states(for transition), transition map and emision map, and return a filled viterbi matrix
//and a matrix to trace back the path.
func FillViterbiNoHiddenStates(str string, sigma, states []string, trmap, emimap MtxMap) ([][]float64, [][]coordinate) {
	var logTranEmi, max float64
	viterbi := make([][]float64, len(states))      //len(states) x len(str)
	backtrace := make([][]coordinate, len(states)) //len(states) x len(str)-1

	for v1 := range viterbi {                      // initialize the first colum to eqal probabilities
		viterbi[v1] = make([]float64, len(str))
		backtrace[v1] = make([]coordinate, len(str)-1)
		viterbi[v1][0] = math.Log(1/float64(len(states))) + math.Log(emimap[states[v1]][str[0:1]])
	}
	//fill in up the viterbi matrix row by row.
	for s := 1; s < len(str); s++ {
		for v2 := range viterbi { //each entry in v2 is the max{each v1 to v2}
			max = viterbi[0][s-1] + math.Log(trmap[states[0]][states[v2]]) + math.Log(emimap[states[v2]][str[s:s+1]])
			for v1 := range viterbi {
				logTranEmi = viterbi[v1][s-1] + math.Log(trmap[states[v1]][states[v2]]) + math.Log(emimap[states[v2]][str[s:s+1]])
				if logTranEmi >= max {
					max = logTranEmi
					backtrace[v2][s-1].x = v1
					backtrace[v2][s-1].y = s - 2
				}
			}
			viterbi[v2][s] = max
		}
	}
	return viterbi, backtrace
}

//Find the ending state in Viterbi matrix. The ending state of the viterbi matrix
//is the maximum value of the last column.
func EndingInViterbi(startpoint int, str string, states []string, viterbi [][]float64) coordinate {
	var endingState coordinate

	if startpoint == 1 {
		endingState.x = len(states) - 1
	} else {
		//Finding the max value in the last column of viterbi, as well as the ending state
		max := viterbi[0][len(str)-1]
		for vf := range viterbi {
			if viterbi[vf][len(str)-1] >= max {
				max = viterbi[vf][len(str)-1]
				endingState.x = vf //the vf th states
			}
		}
	}
	return endingState
}

//Trace back viterbi path from the ending state through backtrace mtx.
func ViterbiTraceBack(endingState coordinate, str string, states []string, backtrace [][]coordinate) string {
	path := states[endingState.x] //x is the state, y is the position of the str
	endingState.y = len(backtrace[0]) - 1

	for endingState.x > 0 {
		endingState = backtrace[endingState.x][endingState.y]
		path = states[endingState.x] + " " + path
	}
	return path
}
