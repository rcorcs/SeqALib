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
  std::string seq1str = "AAAGAATGCAT";
  std::string seq2str = "AAACTCAT";

  // AAA GAA TGCAT
  // | |  |  | |||
  // A A  A CT CAT

  // AAA GAATGCAT
  // |||    | |||
  // AAAC   T CAT

  if (argc>1) {
    seq1str = std::string(argv[1]);
    seq2str = std::string(argv[2]);
  }

  std::vector<char> seq1;
  std::vector<char> seq2;

  for (char c : seq1str) seq1.push_back(c);
  for (char c : seq2str) seq2.push_back(c);

  NeedlemanWunschSA<std::vector<char>,char,'-'> SA(ScoringSystem(-1,2),equal<char>,seq1,seq2); //could also use `nullptr` instead of `equal<char>`
  AlignedSequence<char,'-'> &Alignment = SA.getResult();
  printAlignment(Alignment);

  return 0;
}
