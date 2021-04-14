package main

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"
)

//Takes the file downloaded from Pfam and collect all the alignments as same
//length strings. Add all the strings into a slice of strings.
func ReadAlignmentsPfam(file io.Reader) []string {
	var alignments []string
	br := bufio.NewReader(file)
	for {
		var ID, strg string
		if _, err := fmt.Fscanln(br, &ID, &strg); err != nil {
			break
		}
		alignments = append(alignments, strg)
	}
	return alignments
}

//Takes the file downloaded from BLAST. First find the length of the query string
//with "-" dashes. Then find the position of the string in the line. Collect all
//the aligned string at the exact positoin, replace the front and end spaces with
//"-" dashes. Finally, add all the aligned strings into a slice of strings.
func ReadAlignmentsBLAST(file io.Reader) []string {
	var alignments []string
	scanner := bufio.NewScanner(file)

	var eachline, q, str string
	var n, length, startIdx int
	for i := 0; i < 6; i++ {
		if scanner.Scan() && i == 5 {
			eachline = scanner.Text()
			fmt.Sscanln(eachline, &q, &n, &str)
			startIdx = strings.Index(eachline, str)
			length = len(str)
			alignments = append(alignments, str)
		}
	}

	for scanner.Scan() {
		eachline = scanner.Text()
		if len(eachline) != 0 {
			alignedStr := eachline[startIdx : startIdx+length]
			alignedStr = strings.ReplaceAll(alignedStr, " ", "-")
			alignments = append(alignments, alignedStr)
		}
	}
	return alignments
}

//Turn the map of map into txt file. The file is formated to have a matrix looking,
//neat style. Each entry has a width of 9 to accommodate 4 decimal spaces.
func MapToFile(outFileName string, rowheader, colheader []string, theMap MtxMap) {
	outFile, err3 := os.Create(outFileName)

	if err3 != nil {
		fmt.Println("Error in creating the file.")
	}
	defer outFile.Close()

	fmt.Fprintf(outFile, "%-9s", " ")
	for _, h := range colheader {
		fmt.Fprintf(outFile, "%-9s", h)
	}
	fmt.Fprintln(outFile, " ")

	_, _, outMap := MapToMtx(rowheader, colheader, theMap)

	for i := range outMap {
		fmt.Fprintf(outFile, "%-9s", rowheader[i])
		for j := range outMap[i] {
			fmt.Fprintf(outFile, "%-9.4f", RoundTo(outMap[i][j], 4))
		}
		fmt.Fprintln(outFile, " ")
	}
}

//FileToMap takes a transition or emission file, returns a map that represents
//the transition or emission matrix(map of map). It also returns a column header.
func FileToMap(file io.Reader) (MtxMap, []string) {
	//OutMap := make(MtxMap, len(rowheader))
	OutMap := make(MtxMap)
	scanner := bufio.NewScanner(file)

	scanner.Scan()
	colheader := strings.Fields(scanner.Text())

	var eachline []string
	var entryStr, rheader string
	var rowEntry map[string]float64
	
	for scanner.Scan() {
		eachline = strings.Fields(scanner.Text())
		rheader = eachline[0]
		rowEntry = make(map[string]float64)
		for c := 1; c < len(eachline); c++ {
			entryStr = eachline[c]
			if entry, err := strconv.ParseFloat(entryStr, 64); err == nil {
				rowEntry[colheader[c-1]] = entry
			}
			OutMap[rheader] = rowEntry
		}
	}
	return OutMap, colheader
}

//This function has assumed user input, as asked from upper level function. It
//takes in name of the files (transition and emission files), and return the
//the contents of the file as map or map for transition and emision matrix. It
//also returns the header (states) and sigmas (amino acid).
func ReadInStrNMap() (MtxMap, MtxMap, []string, []string) {
	reader := bufio.NewReader(os.Stdin)

	trmapName, err2 := reader.ReadString('\n')
	if err2 != nil {
		panic("Trmap name read in error.")
	}
	trmapName = strings.TrimSuffix(trmapName, "\n")

	emimapName, err3 := reader.ReadString('\n')
	if err3 != nil {
		panic("Emimap name read in error.")
	}
	emimapName = strings.TrimSuffix(emimapName, "\n")

	trfile, err4 := os.Open(trmapName)
	if err4 != nil {
		fmt.Println("Error: something wrong with openning input files.")
	}
	defer trfile.Close()

	emifile, err4 := os.Open(emimapName)
	if err4 != nil {
		fmt.Println("Error: something wrong with openning input files.")
	}
	defer emifile.Close()

	trmap, header := FileToMap(trfile)
	emimap, sigma := FileToMap(emifile)
	return trmap, emimap, header, sigma
}
