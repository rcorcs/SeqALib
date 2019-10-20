# SeqALib: A Library for Sequence Alignment

SeqALib contains efficient implementation of sequence alignment algorithms from
Bioinformatics.

![Alignment Example](https://raw.githubusercontent.com/rcorcs/SeqALib/master/doc/alignment-example.png)

The algorithms currently provided by SeqALib are:
* [Needleman-Wunsch](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm)
* [Hirschberg](https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm)
* [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
* [Gotoh Local and Global](http://helios.mi.parisdescartes.fr/~lomn/Cours/BI/Material2019/gap-penalty-gotoh.pdf)

AUTHOR: Rodrigo Rocha and Sean Stirling

## Easy to use

See full example in the file: `test/Test.cpp`

```cpp
  std::string seq1 = "AAAGAATGCAT";
  std::string seq2 = "AAACTCAT";

  NeedlemanWunschSA<std::string,char,'-'> SA(ScoringSystem(-1,2));
  AlignedSequence<char,'-'> Alignment = SA.getAlignment(seq1,seq2);

  // The resultng Alignment contains:
  // AAA GAATGCAT
  // |||    | |||
  // AAAC   T CAT
```

