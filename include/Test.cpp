#include <iostream>
#include <string>

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

	//std::string seq1 = "AATCG";
	//std::string seq2 = "AACG";
	//std::string seq1 = "AGGATCGGCTAGAGCTAGAGCTAGCTAGTAGC";
	//std::string seq2 = "GAGATCGGCGGATTACAGGCTATCGA";
	std::string seq1 = "GATTACA";
	std::string seq2 = "AACG";

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

  NeedlemanWunschSA<std::string,char,'-'> SA(ScoringSystem(-1,2),equal<char>); //could also use `nullptr` instead of `equal<char>`
  AlignedSequence<char,'-'> NWAlignment = SA.getAlignment(seq1,seq2);

  //SmithWatermanSA<std::string, char, '-'> SA(ScoringSystem(-1,1), equal<char>);
  //AlignedSequence<char,'-'> SWAlignment =  SA.getAlignment(seq1, seq2);

  ////An even simpler way of using it with the default settings
  //AlignedSequence<char,'-'> NWAlignment = NeedlemanWunschSA<std::string,char,'-'>().getAlignment(seq1,seq2);

  std::cout << "# Needleman-Wunsch:" << std::endl;
  printAlignment(NWAlignment);
  std::cout << "" << std::endl;
  
  AlignedSequence<char,'-'> HAlignment = HirschbergSA<std::string,char,'-'>(ScoringSystem(-1,2),equal<char>).getAlignment(seq1,seq2);
  std::cout << "# Hirschberg:" << std::endl;
  printAlignment(HAlignment);
  std::cout << "" << std::endl;

  AlignedSequence<char, '-'> SWAlignment = SmithWatermanSA<std::string, char, '-'>(ScoringSystem(-1, 2, -1, true), equal<char>).getAlignment(seq1, seq2);
  std::cout << "# Smith-Waterman:" << std::endl;
  printAlignment(SWAlignment);
  std::cout << "" << std::endl;

  AlignedSequence<char, '-'> LGAlignment = LocalGotohSA<std::string, char, '-'>(ScoringSystem(-1, 2, -1, true), equal<char>).getAlignment(seq1, seq2);
  std::cout << "# Local Gotoh:" << std::endl;
  printAlignment(LGAlignment);
  std::cout << "" << std::endl;

  AlignedSequence<char, '-'> GGAlignment = GlobalGotohSA<std::string, char, '-'>(ScoringSystem(-1, 2, -1, true), equal<char>).getAlignment(seq1, seq2);
  std::cout << "# Global Gotoh:" << std::endl;
  printAlignment(GGAlignment);
  std::cout << "" << std::endl;

  return 0;
}
