#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "SequenceAlignment.h"
#include "FASTA.h"

template<typename Ty, Ty Blank>
void printAlignment(AlignedSequence<Ty, Blank> &Result) {
  unsigned CountMismatches = 0;
  unsigned CountMatches = 0;
  for (auto &Entry : Result) {
    std::cout << Entry.get(0);
    if (Entry.match()) {
       CountMatches++;
       std::cout << "-----";
    } else {
       CountMismatches++;
       std::cout << "     ";
    }
    std::cout << Entry.get(1);
    std::cout << std::endl;
  }
  std::cout << "Matches: " << CountMatches << " (" << (float(CountMatches)/float(CountMismatches+CountMatches)*100.f) << "%)" << std::endl;
  std::cout << "Mismatches: " << CountMismatches << " (" << (float(CountMismatches)/float(CountMismatches+CountMatches)*100.f) << "%)" << std::endl;
}

int main(int argc, char** argv) {
  std::vector<char> Seq1;
  std::vector<char> Seq2;
  
  readFASTA(argv[1],Seq1);
  readFASTA(argv[2],Seq2);

  //AlignedSequence<char,'-'> NWAlignment = NeedlemanWunschSA<std::string,char,'-'>().getAlignment(seq1,seq2); 
  AlignedSequence<char,'-'> HAlignment = HirschbergSA<std::vector<char>,char,'-'>(ScoringSystem(-1,1,-1)).getAlignment(Seq1,Seq2);

  std::cout << "# Hirschberg:" << std::endl;
  printAlignment(HAlignment);
  std::cout << "Size #1: " << Seq1.size() << std::endl;
  std::cout << "Size #2: " << Seq2.size() << std::endl;

  return 0;
}
