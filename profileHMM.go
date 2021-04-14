//This is the main function that builds the profile HMM (transition and emission).
//I have implement the transition and emission matrix as map of maps to have easier
//access of values. A map of map is also easier than matrix to add or delete
//values in case that we need to enlarge or shrink our number of states.

package main

import (
	"math/rand"
	"time"
)

type MtxMap map[string]map[string]float64

//Input: A threshold θ, followed by Σ, followed by a multiple alignment.
//theta indicates whether we should determine the state at that position as deletion
//according to how many real symbols we found at that specific position.
//pseudoCount is a integer boolean with 0 means don't count pseudocount and 1 means
//do count pseudocount. Sigma is the emitted symbols such as amino acid. multiAlign
//is the incoming aligned data tha we are using to construct our transition and emission maps.
//Output: The transition and emission probabilities of the profile HMM HMM(Alignment, θ).
func ProfileHMM(theta, pseudoCount float64, sigma, multiAlign []string) (MapHeader []string, trmap, emimap MtxMap) {
	if len(multiAlign) == 0 || len(multiAlign[0]) == 0 {
		panic("Invalid training data. Failed to construct ProfileHMM.")
	}

	md := MarkDeletionState(theta, multiAlign) //slice of integers {0,1} mark match or deletion state
	alignSize := SumOfIntSlice(md)             //int
	MapHeader = MakeMapHeader(alignSize)       //[]string
	eachLength := len(multiAlign[0])

	trmap = ProfileTrMap(multiAlign, MapHeader, eachLength, md)
	emimap = ProfileEmimap(multiAlign, sigma, MapHeader, eachLength, md)
	//pseudoCount is either 1 or 0 with 1 indicating we want to involve pseudoCount
	//in our matrix and 0 indicating we are not adding pseuroCount in our matrix.
	if pseudoCount != 0.0 {
		trmapWithPseudo := TrmapPseudoCount(MapHeader, pseudoCount, alignSize, trmap)
		emimapWithPseudo := EmimapPseudoCount(pseudoCount, emimap)
		(&trmapWithPseudo).Normalize()
		(&emimapWithPseudo).Normalize()
		return MapHeader, trmapWithPseudo, emimap
	}
	return
}

//Normalize takes a map such as transition map and emission map, collect the sum
//of each row and divide each value of the row to that sum, thus normalize all values (sum of each row adds to 1).
func (anymap *MtxMap) Normalize() { // map[string]map[string]float64
	var sumEachRow float64
	themap := *anymap

	for row := range themap {
		sumEachRow = 0.0
		for col := range themap[row] {
			sumEachRow += themap[row][col]
		}
		for col := range themap[row] {
			if themap[row][col] != 0.0 {
				themap[row][col] = themap[row][col] / sumEachRow
			}
		}
	}
}

//Domain generator path takes in a transition map, it generates fictional domain
//path for DomainSequenceGenerator. It generates the path according to the probabilities
//of the transition map. To have flexibility, we do not choose the transition or Emission
//state according to the highest probablity, instead, we choose randomly according
//to the probability.
func DomainPathGenerator(trmap MtxMap, header []string) []string {
	rand.Seed(time.Now().UnixNano())

	var row, col string
	var fictionalPath []string
	var randNum float64

	for r := 0; r < len(trmap); r++ {
    //since each line of transition rates adds up to 1, rand.Float64() generates a
    //number between 0 and 1, gives each transition a fair proportion of chance.
		randNum = rand.Float64()
		for c := 0; c < len(trmap); c++ {
			row = header[r]
			col = header[c]
			if trmap[row][col] != 0.0 {
				randNum -= trmap[row][col]
				if randNum <= 0 {
					fictionalPath = append(fictionalPath, header[c])
          r = c
					c = len(trmap) - 1

				}
			}
		}
	}
	return fictionalPath
}

//DomainSeqGenerator takes in path that's generated by DomainPathGenerator, sigmas,
//the emission map and generates the fictionalSeq accordingly. It does not pick the
//highest emission rate. Insead, it picks randomly according to emission rates.
func DomainSeqGenerator(path, sigmas []string, emimap MtxMap) string{
  rand.Seed(time.Now().UnixNano())
  var sequence string
  var randNum float64
	
  //since each line of emission rates adds up to 1, rand.Float64() generates a
  //number between 0 and 1, gives each emission result the fair proportion of chance.
  for _, p := range path{
    randNum = rand.Float64()
    for c := 0; c < len(emimap[p]); c++{
      randNum -= emimap[p][sigmas[c]]
      if randNum <=0 {
        sequence = sequence + sigmas[c]
        c = len(emimap[p]) - 1
      }
    }
  }
  return sequence
}

//NonEmissionStateEist checks if there are states in the emission map that does
//not emit any symbols. Deletion states are such non emission states.
func NonEmissionStateExist(emimap MtxMap) bool {
	var sum float64
	for row := range emimap {
		sum = 0
		for col := range emimap[row] {
			sum += emimap[row][col]
		}
		if sum == 0 {
			return true
		}
	}
	return false
}