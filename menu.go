//This file contains the different options to run the program as detailed in the
//main go file. Specifically ,it contains:
//1. Build profile HMM given alignments.
//2. Given profile HMM and a string, the probablity that this string is emmited by the sequence.
//3. Given profile HMM and a string, the most probable path aligning this string to the HMM.
//4. Given profile HMM, generate fictional strings that belong to the group.
package main

import (
	"bufio"
	"fmt"
	"os"
	"strings"
  "strconv"
)

//OPTION1: produce profile HMM
func Option1(amino []string, theta float64) {
	reader := bufio.NewReader(os.Stdin)

	fmt.Println("\nTo produce profile HMM, you just need one alignment file from Pfam or BLAST.")
	fmt.Println(" - If you got the file from Pfam, press P; if from Blast, press B.")
	PorB, err1 := reader.ReadString('\n')
	if err1 != nil {
		panic("PorB read in error.")
	}
	fmt.Println("\nPlease enter the code of domain family.")

	domain, err2 := reader.ReadString('\n')
	if err2 != nil {
		panic("domain code read in error.")
	}
	domain = strings.TrimSuffix(domain, "\n")

	fmt.Println("\nPlease enter file name with path (including .txt).")
	filename, err3 := reader.ReadString('\n')
	if err3 != nil {
		panic("filename read in error.")
	}
	filename = strings.TrimSuffix(filename, "\n")

	file, err4 := os.Open(filename)
	if err4 != nil {
		fmt.Println("Error: something wrong with openning input files.")
	}
	defer file.Close()

	var multiAlign []string
	if PorB == "P\n" || PorB == "p\n" {
		multiAlign = ReadAlignmentsPfam(file)
	} else if PorB == "B\n" || PorB == "b\n" {
		multiAlign = ReadAlignmentsBLAST(file)
	}

	header, trmap, emimap := ProfileHMM(theta, 0.01, amino, multiAlign)
	MapToFile(domain+"EmiMap.txt", header, amino, emimap)
	fmt.Println("Emission matrix of ProfileHMM produced! Find it as " + domain + "EmiMap.txt")
	MapToFile(domain+"TrMap.txt", header, header, trmap)
	fmt.Println("Transition matrix of ProfileHMM produced! Find it as " + domain + "TrMap.txt")
}

//OPTION2: if a sequence belong to a domain family
func Option2() {
	reader := bufio.NewReader(os.Stdin)
	fmt.Println("\nTo check if a sequence belong to a domain family, we would need: ")
	fmt.Println("    1. The sequence itself (with not dashes or space)")
	fmt.Println("    2. The transition and emission matrix")
	fmt.Println("Please enter the sequence(w/o dashes), file name of transition map and emission map, each on a new line. ")

	str, err1 := reader.ReadString('\n')
	if err1 != nil {
		panic("String read in error.")
	}
	str = strings.TrimSuffix(str, "\n")

	trmap, emimap, header, sigma := ReadInStrNMap()

	nulltrmap := NullTrMap(trmap)
	nullemimap := NullEmiMapProtein(emimap)

	Ha := Forward(str, sigma, header, trmap, emimap)
	H0 := Forward(str, sigma, header, nulltrmap, nullemimap)
	likelihood := LogLikeliHood(Ha, H0)

	if likelihood >= 1 {
		fmt.Println("This sequence likely belongs to the domain group compare to our null model.\n LogLikeliHood: ", likelihood)
	} else {
		fmt.Println("This sequence likely belongs to our null model compare to the group .\n LogLikeliHood: ", likelihood)
	}
}

//OPTION3: check the most probably path of a sequence with a given HMM.
func Option3() {
	reader := bufio.NewReader(os.Stdin)
	fmt.Println("\nTo find the most probable path of a sequence aligning to a HMM, we would need: ")
	fmt.Println("    1. The sequence itself (with not dashes or space)")
	fmt.Println("    2. The transition and emission matrix")
	fmt.Println("Please enter the sequence(w/o dashes), file names of transition map and emission map, each on a new line. ")

	str, err1 := reader.ReadString('\n')
	if err1 != nil {
		panic("String read in error.")
	}
	str = strings.TrimSuffix(str, "\n")

	trmap, emimap, header, sigma := ReadInStrNMap()

	NonEmissionTF := 0 //assume no non emission states for now.
	if NonEmissionStateExist(emimap) == true {
		NonEmissionTF = 1
	}

	path := ViterbiDecoding(NonEmissionTF, str, sigma, header, trmap, emimap)
	fmt.Println("The most probable path is: \n\n", path)

  prPath  := PrHiddenPath(NonEmissionTF, path, header, trmap)
  fmt.Println("\nThe probablity of emitting this path from the transition map is: \n\n", prPath)
}

//OPTION4: generate fictional domain sequences with profile HMM.
func Option4() {
  reader := bufio.NewReader(os.Stdin)
	fmt.Println("\nTo generate fictional sequences with an existing HMM, we would need: ")
	fmt.Println("  1. The transition and emission matrix")
	fmt.Println("\nPlease enter how many sequence you want, the file names of transition map and emission map, each on a new line. ")

  numSeqStr, err1 := reader.ReadString('\n')
	if err1 != nil {
		panic("String read in error.")
	}
	numSeq, err2 := strconv.Atoi(strings.TrimSuffix(numSeqStr, "\n"))
  if err2 != nil {
		panic("Problem converting read in value into integer.")
	}

  trmap, emimap, header, sigma := ReadInStrNMap()
  var probablePath []string
  var probablyDomainSeq string

	fmt.Println("\nHere are the fictional domain sequences. You can Blast them and see if you got lucky!!\n")
  for n := 0; n < numSeq; n++{
    probablePath = DomainPathGenerator(trmap, header)
    probablyDomainSeq = DomainSeqGenerator(probablePath, sigma, emimap)
    fmt.Println(probablyDomainSeq)
  }
}
