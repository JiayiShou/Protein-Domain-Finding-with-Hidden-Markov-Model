package main

import (
	"fmt"
)

//ProfileTrMap takes the raw transition map (that only has counts of each transition)
//and make it into percentage by devide each position by totalVisits.
func ProfileTrMap(multiAlign, MapHeader []string, eachLength int, md []int) MtxMap {
	trmap, totalVisits := TrmapRaw(multiAlign, MapHeader, eachLength, md)

	for pre := range trmap {
		for after := range trmap[pre] {
			if trmap[pre][after] != 0 {
				trmap[pre][after] = trmap[pre][after] / float64(totalVisits[pre])
			}
		}
	}
	return trmap
}

//TrmapRaw counts each transition and store the counts into the transition map.
//In specific, it takes the multialign slice, map headers, each alignment eachLength
//and a slice of integer md that shows whether that position is a deleiton state.
//It returns a transition matrix that only has the counts but are not normalized.
func TrmapRaw(multiAlign, MapHeader []string, eachLength int, md []int) (MtxMap, map[string]int) {
	trmap := CreatEmptyMap(MapHeader, MapHeader)

	totalVisits := make(map[string]int)
	for _, h := range MapHeader {
		totalVisits[h] = 0.0
	}
	totalVisits["Start"] = len(multiAlign)

	var curr, prev string
	for t := range multiAlign {
		prev = "Start"
		for s := 0; s < eachLength; s++ {
			//Check if the state is a deletion state or not. If it's a deletion state, we
			//only need to count insertions, if its not a deletion state, we count everything.
			if md[s] == 1 { // it's not a deletion position
				if multiAlign[t][s:s+1] != "-" {
					curr = "M" + fmt.Sprintf("%d", SumOfIntSlice(md[0:s+1]))
				} else {
					curr = "D" + fmt.Sprintf("%d", SumOfIntSlice(md[0:s+1]))
				}
				trmap[prev][curr] += 1.0
				totalVisits[curr] += 1
				prev = curr
			} else if md[s] == 0 {
				// in a deletion state, a simbol other than dash is counter as insertion
				if multiAlign[t][s:s+1] != "-" {
					curr = "I" + fmt.Sprintf("%d", SumOfIntSlice(md[0:s+1]))
					trmap[prev][curr] += 1.0
					totalVisits[curr] += 1
					prev = curr
				}
			}

			if s == eachLength-1 {
				trmap[curr]["End"] += 1.0
			}
		}
	}
	return trmap, totalVisits
}

//TrmapPseudoCount takes a transition map that has the counts but not normalized, and add
//the pseudoCount to each transitions that's nor represented in the alignments.
//It returns an unormalized transition map.
func TrmapPseudoCount(MapHeader []string, pseudoCount float64, alignSize int, trmap MtxMap) MtxMap {
	threeStates := []string{"M", "D", "I"}
	beginStates := []string{"Start", "I0"}

	for _, f := range beginStates {
		for _, t := range threeStates {
			if t == "I" {
				t = t + fmt.Sprintf("%d", 0)
			} else {
				t = t + fmt.Sprintf("%d", 1)
			}
			trmap[f][t] += pseudoCount
		}
	}

	for i := 1; i <= alignSize; i++ {
		for _, f := range threeStates {
			f = f + fmt.Sprintf("%d", i)
			for _, t := range threeStates {
				if t == "I" {
					t = t + fmt.Sprintf("%d", i)
					trmap[f][t] += pseudoCount
				} else if i != alignSize {
					t = t + fmt.Sprintf("%d", i+1)
					trmap[f][t] += pseudoCount
				}
			}
			if i == alignSize {
				trmap[f]["End"] += pseudoCount
			}
		}
	}
	return trmap
}
