package main

import (
	"fmt"
)

//ProfileEmimap takes the raw emission map (that only has counts of each emission)
//and normalize it by devide each position by totalVisits.
func ProfileEmimap(multiAlign, sigma, MapHeader []string, eachLength int, md []int) MtxMap {
	emimap, totalVisits := EmimapRaw(multiAlign, sigma, MapHeader, eachLength, md)
	for state := range emimap {
		for letter := range emimap[state] {
			if emimap[state][letter] != 0 {
				emimap[state][letter] = emimap[state][letter] / float64(totalVisits[state])
			}
		}
	}
	return emimap
}

//EmimapRaw counts each emission and store the counts into the emission map.
//It takes in multialignments slice, sigma(emitted symbols, eg amino acid),
//mapheader(states), length of each alignment and md that indicate if a state is
//a deletion state.
//It return a raw emition map and total visits.
func EmimapRaw(multiAlign, sigma, MapHeader []string, eachLength int, md []int) (MtxMap, map[string]int) {
	emimapRaw := CreatEmptyMap(MapHeader, sigma)

	totalVisits := make(map[string]int)
	for _, h := range MapHeader {
		totalVisits[h] = 0.0
	}
	//totalVisits := make([]int, len(md))
	var curr string
	for s := 0; s < eachLength; s++ {
		for t := range multiAlign {

			if md[s] == 1 {
				if multiAlign[t][s:s+1] != "-" {
					curr = "M" + fmt.Sprintf("%d", SumOfIntSlice(md[0:s+1]))
				} else {
					curr = "D" + fmt.Sprintf("%d", SumOfIntSlice(md[0:s+1]))
				}
				emimapRaw[curr][multiAlign[t][s:s+1]] += 1.0
				totalVisits[curr] += 1
			} else if md[s] == 0 {
				if multiAlign[t][s:s+1] != "-" {
					curr = "I" + fmt.Sprintf("%d", SumOfIntSlice(md[0:s+1]))
					emimapRaw[curr][multiAlign[t][s:s+1]] += 1.0
					totalVisits[curr] += 1
				}
			}
		}
	}
	return emimapRaw, totalVisits
}

//EmimapPseudoCount takes a emimap and add the pseudoCount to each emission from
//insertion and matching, but not start, end, nor deletion because at those states,
//there should not be any emission.
func EmimapPseudoCount(pseudoCount float64, emimap MtxMap) MtxMap {

	for row := range emimap {
		if row[0:1] == "I" || row[0:1] == "M" {
			for col := range emimap[row] {
				emimap[row][col] += pseudoCount
			}
		}
	}
	return emimap
}
