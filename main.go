//This is the main function of my program involving protein domain and hidden Markov
//Model. This file mainly interact with uses and get in user's request about what
//their requests are. After getting the request, this function calls sub option functions
//in menu.go.

//Programming for Scientist
//Jiayi Shou Dec.11th 2020

package main

import (
	"bufio"
	"fmt"
	"os"
)

func main() {
	amino := []string{"G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T"}
	theta := 0.4 //I set a arbitrary theta value since we are not improving parameters in my program.

	fmt.Println("\n\n\n--------------------------------------<just a casual dividing line>--------------------------------------")
	fmt.Println("\nWelcome to Profile Hidden Markov Model for Protein Domains!!\n")
	reader := bufio.NewReader(os.Stdin)
	fmt.Println(" - To produce profile HMM, please press 1 and enter;")
	fmt.Println(" - To check the chance of a sequece belong to domain family, please press 2 and enter;")
	fmt.Println(" - To see the most probable path of a sequence aligning to a HMM, please press 3 and enter;")
  fmt.Println(" - To generate fictional domain sequences with profile HMM, please press 4 and enter.")

	optionFunction, _ := reader.ReadString('\n')

	if optionFunction == "1\n" {
		//OPTION1: produce profile HMM
		Option1(amino, theta)
	} else if optionFunction == "2\n" {
		//OPTION2: if a sequence belong to a domain family
		Option2()
	} else if optionFunction == "3\n" {
    //OPTION3: check the most probably path of a sequence with a given HMM.
		Option3()
	} else if optionFunction == "4\n" {
    //OPTION4: generate fictional domain sequences with profile HMM.
		Option4()
	} else {
		fmt.Println("\nOooops, didn't match anything. Program end.\n")
	}
	fmt.Println("\n--------------------------------------------------<END>--------------------------------------------------\n")
}
