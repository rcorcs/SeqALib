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
  
  std::string method;
  std::vector<std::string> FNames;
  for (unsigned i = 1; i<argc; i++) {
    if (std::string("-m")==std::string(argv[i-1])) {
      method = argv[i];
    } else if (std::string("-m")!=std::string(argv[i])) {
      FNames.push_back(std::string(argv[i]));
    }
  }
  if (FNames.size()!=2) {
    std::cout << "Usage:\t" << argv[0] << " [Options] <input FAST> <input FAST>\n";
    std::cout << "Options:\n\t-m <Method>\tName of the algorithm used for sequence alignment\n";
    return 0;
  }

  readFASTA(FNames[0],Seq1);
  readFASTA(FNames[1],Seq2);

  ScoringSystem SS(-1,1,-1);

  AlignedSequence<char,'-'> Alignment;
  if (method.empty()) { Alignment = HirschbergSA<std::vector<char>,char,'-'>(SS).getAlignment(Seq1,Seq2); }
  else switch(method[0]) {
  case 'n': Alignment = NeedlemanWunschSA<std::vector<char>,char,'-'>(SS).getAlignment(Seq1,Seq2);
            break;
  case 'h': Alignment = HirschbergSA<std::vector<char>,char,'-'>(SS).getAlignment(Seq1,Seq2);
            break;
  case 's': Alignment = ShortSightedSA<std::vector<char>,char,'-'>(SS,nullptr,4096).getAlignment(Seq1,Seq2);
            break;
  case 'g': Alignment = GlobalGotohSA<std::vector<char>,char,'-'>(SS).getAlignment(Seq1,Seq2);
            break;
  }

  printAlignment(Alignment);
  std::cout << "Size #1: " << Seq1.size() << std::endl;
  std::cout << "Size #2: " << Seq2.size() << std::endl;

  return 0;
}
