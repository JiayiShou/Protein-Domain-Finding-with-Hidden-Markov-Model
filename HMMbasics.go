//This file has functions: PrHiddenPath, PrStrGivenPath, Forward
package main

import (
  "strings"
)

//input: hidden path π, States and transition matrix Transition of an HMM (Σ, States, Transition, Emission).
//output: The probability of this path, Pr(π). Assuming that initial probabilities are equal.
func PrHiddenPath(startpoint int, pathStr string, states []string, trmap MtxMap) float64 {
  path := strings.Fields(pathStr)

  if len(path) == 0 {
		return -1.0
	}
  var pr float64

  if startpoint == 0{
    pr = 1 / float64(len(states))
  } else {
    pr = 1
  }

	previous := path[0]
	var current string
	for p := 1; p < len(path); p++ {
		current = path[p]
		pr = pr * trmap[previous][current]
		previous = current
	}
	return pr
}

//Input: A string x, Σ , States, transition map, and emission map of an HMM (Σ, States, Transition, Emission).
//output: The probability Pr(x) that the HMM emits x.
func Forward(str string, sigma, states []string, trmap, emimap MtxMap) (pr float64) {
	if len(str) == 0 {
		return
	}
	forwardMtx := make([][]float64, len(states))
	for f := range forwardMtx { // initialize the first colum
		forwardMtx[f] = make([]float64, len(str))
		forwardMtx[f][0] = 1 / float64(len(states)) * emimap[states[f]][str[0:1]]
	}

	for s := 1; s < len(str); s++ {
		for f2 := range forwardMtx {
			for f1 := range forwardMtx {
				forwardMtx[f2][s] += forwardMtx[f1][s-1] * trmap[states[f1]][states[f2]] * emimap[states[f2]][str[s:s+1]]
			}
		}
	}
	for f2 := range forwardMtx {
		pr += forwardMtx[f2][len(str)-1]
	}
	return
}

//The null map forms as the null hypothesis for future log likelihood calculation.
//Input: a transition map as a model for the null map.
//Output: a null transition map that has equal possibility to transit to any next states.
func NullTrMap(trmap MtxMap) MtxMap {
	var entries int
	Nullmap := make(MtxMap, len(trmap))

	for row := range trmap {
		entries = 0
		Nullmap[row] = make(map[string]float64, len(trmap[row]))
		for col := range trmap[row] {
			if trmap[row][col] != 0 {
				entries += 1
			}
		}
		for col := range trmap[row] {
			if trmap[row][col] == 0.0 {
				Nullmap[row][col] = 0.0
			} else {
				Nullmap[row][col] = 1 / float64(entries)
			}
		}
	}
	return Nullmap
}

//NullEmiMapProtein takes an emimap and return the emimap with same states but
//the emission rates are background existence rate of each amino acid.
func NullEmiMapProtein(emimap MtxMap) MtxMap {
	Nullmap := make(MtxMap, len(emimap))

	for row := range emimap {
		Nullmap[row] = make(map[string]float64, len(emimap[row]))
		if row[0:1] != "S" && row[0:1] != "E" && row[0:1] != "D" {
			Nullmap[row]["A"] = 0.074
			Nullmap[row]["R"] = 0.042
			Nullmap[row]["N"] = 0.044
			Nullmap[row]["D"] = 0.059
			Nullmap[row]["C"] = 0.033

			Nullmap[row]["E"] = 0.058
			Nullmap[row]["Q"] = 0.037
			Nullmap[row]["G"] = 0.074
			Nullmap[row]["H"] = 0.029
			Nullmap[row]["I"] = 0.038

			Nullmap[row]["L"] = 0.076
			Nullmap[row]["K"] = 0.072
			Nullmap[row]["M"] = 0.018
			Nullmap[row]["F"] = 0.04
			Nullmap[row]["P"] = 0.05

			Nullmap[row]["S"] = 0.08
			Nullmap[row]["T"] = 0.062
			Nullmap[row]["W"] = 0.013
			Nullmap[row]["Y"] = 0.033
			Nullmap[row]["V"] = 0.068
		}
	}
	return Nullmap
}

/*
//input: A string x, Σ, hidden path π, States and emission matrix(map) of an HMM (Σ, States, Transition, Emission).
//output: The conditional probability Pr(x|π) that string x will be emitted by the HMM given the hidden path π.
func PrStrGivenPath(str, path string, sigma, states []string, emimap MtxMap) float64 {
	pr := 1.0
	for p := 0; p < len(path); p++ { //go through path one by one, p is the index
		pr *= emimap[path[p:p+1]][str[p:p+1]]
	}
	return pr
}
*/
