#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <cmath>
#include <algorithm>

#include "SequenceAlignment.h"

template<typename T>
bool equal(T V1, T V2) { return V1==V2; }

template<typename Ty, Ty Blank>
void printAlignment(AlignedSequence<Ty, Blank> &Result) {
  for (auto &Entry : Result) {
    std::cout << Entry.get(0);
  }
  std::cout << std::endl;
  for (auto &Entry : Result) {
    if (Entry.match())
    std::cout << "|";
    else std::cout << " ";
  }
  std::cout << std::endl;
  for (auto &Entry : Result) {
    std::cout << Entry.get(1);
  }
  std::cout << std::endl;
}

int main(int argc, char** argv) {
  std::string seq1 = "AAAGAATGCAT";
  std::string seq2 = "AAACTCAT";

  // AAA GAA TGCAT
  // | |  |  | |||
  // A A  A CT CAT

  // AAA GAATGCAT
  // |||    | |||
  // AAAC   T CAT

  if (argc>1) {
    seq1 = std::string(argv[1]);
    seq2 = std::string(argv[2]);
  }

  NeedlemanWunschSA<std::string,char,'-'> SA(ScoringSystem(-1,2),equal<char>,seq1,seq2); //could also use `nullptr` instead of `equal<char>`
  AlignedSequence<char,'-'> &Alignment = SA.getResult();
  printAlignment(Alignment);

  return 0;
}
