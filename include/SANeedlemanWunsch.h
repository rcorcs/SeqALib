template<typename ContainerType, typename Ty=typename ContainerType::value_type, Ty Blank=Ty(0), typename MatchFnTy=std::function<bool(typename ContainerType::value_type,typename ContainerType::value_type)>>
class NeedlemanWunschSA : public SequenceAligner<ContainerType,Ty,Blank,MatchFnTy> {
private:
  ScoreSystemType *Matrix;
  size_t MatrixRows;
  size_t MatrixCols;
  bool *Matches;
  size_t MatchesRows;
  size_t MatchesCols;


  const static unsigned END = 0;
  const static unsigned DIAGONAL = 1;
  const static unsigned UP = 2;
  const static unsigned LEFT = 3;

  size_t MaxRow;
  size_t MaxCol;

  using BaseType = SequenceAligner<ContainerType,Ty,Blank,MatchFnTy>;

  void cacheAllMatches() {
    if (BaseType::getMatchOperation()==nullptr) {
      Matches = nullptr;
      return;
    }
    auto &Seq1Ref = BaseType::getSequence(0);
    auto &Seq2Ref = BaseType::getSequence(1);
    Ty *Seq1 = &(Seq1Ref[0]);
    Ty *Seq2 = &(Seq2Ref[0]);
    const size_t SizeSeq1 = Seq1Ref.size();
    const size_t SizeSeq2 = Seq2Ref.size();

    MatchesRows = SizeSeq1;
    MatchesCols = SizeSeq2;
    Matches = new bool[SizeSeq1*SizeSeq2];
    for (unsigned i = 0; i < SizeSeq1; i++)
      for (unsigned j = 0; j< SizeSeq2; j++)
        Matches[i*SizeSeq2 + j] = BaseType::match(Seq1[i],Seq2[j]);
  }

  void computeScoreMatrix() {
    const size_t SizeSeq1 = BaseType::getSequence(0).size();
    const size_t SizeSeq2 = BaseType::getSequence(1).size();

    const size_t NumRows = SizeSeq1 + 1;
    const size_t NumCols = SizeSeq2 + 1;
    Matrix = new int[NumRows * NumCols];
    MatrixRows = NumRows;
    MatrixCols = NumCols;

    ScoringSystem &Scoring = BaseType::getScoring();
    const ScoreSystemType Gap = Scoring.getGapPenalty();
    const ScoreSystemType Match = Scoring.getMatchProfit();
    const bool AllowMismatch = Scoring.getAllowMismatch();
    const ScoreSystemType Mismatch = AllowMismatch
                                   ?Scoring.getMismatchPenalty()
                                   :std::numeric_limits<ScoreSystemType>::min();

    for (unsigned i = 0; i < NumRows; i++)
      Matrix[i * NumCols + 0] = i * Gap;
    for (unsigned j = 0; j < NumCols; j++)
      Matrix[0 * NumCols + j] = j * Gap;

    ScoreSystemType MaxScore = std::numeric_limits<ScoreSystemType>::min();
    if (Matches) {
      if (AllowMismatch) {
        for (unsigned i = 1; i < NumRows; i++) {
          for (unsigned j = 1; j < NumCols; j++) {
            ScoreSystemType Similarity = Matches[(i - 1)*MatchesCols + j - 1] ? Match : Mismatch;
            ScoreSystemType Diagonal = Matrix[(i - 1) * NumCols + j - 1] + Similarity;
            ScoreSystemType Upper = Matrix[(i - 1) * NumCols + j] + Gap;
            ScoreSystemType Left = Matrix[i * NumCols + j - 1] + Gap;
            ScoreSystemType Score = std::max(std::max(Diagonal, Upper), Left);
            Matrix[i * NumCols + j] = Score;
            if (Score>=MaxScore) {
              MaxScore = Score;
              MaxRow=i;
              MaxCol=j;
            }
          }
        }
      } else {
        for (unsigned i = 1; i < NumRows; i++) {
          for (unsigned j = 1; j < NumCols; j++) {
            ScoreSystemType Diagonal = Matches[(i - 1)*MatchesCols + j - 1] ? (Matrix[(i - 1) * NumCols + j - 1] + Match) : Mismatch;
            ScoreSystemType Upper = Matrix[(i - 1) * NumCols + j] + Gap;
            ScoreSystemType Left = Matrix[i * NumCols + j - 1] + Gap;
            ScoreSystemType Score = std::max(std::max(Diagonal, Upper), Left);
            Matrix[i * NumCols + j] = Score;
            if (Score>=MaxScore) {
              MaxScore = Score;
              MaxRow=i;
              MaxCol=j;
            }
          }
        }
      }
    } else {
      auto &Seq1Ref = BaseType::getSequence(0);
      auto &Seq2Ref = BaseType::getSequence(1);
      Ty *Seq1 = &(Seq1Ref[0]);
      Ty *Seq2 = &(Seq2Ref[0]);


      if (AllowMismatch) {
        for (unsigned i = 1; i < NumRows; i++) {
          for (unsigned j = 1; j < NumCols; j++) {
            ScoreSystemType Similarity = (Seq1[i - 1]==Seq2[j - 1]) ? Match : Mismatch;
            ScoreSystemType Diagonal = Matrix[(i - 1) * NumCols + j - 1] + Similarity;
            ScoreSystemType Upper = Matrix[(i - 1) * NumCols + j] + Gap;
            ScoreSystemType Left = Matrix[i * NumCols + j - 1] + Gap;
            ScoreSystemType Score = std::max(std::max(Diagonal, Upper), Left);
            Matrix[i * NumCols + j] = Score;
            if (Score>=MaxScore) {
              MaxScore = Score;
              MaxRow=i;
              MaxCol=j;
            }
          }
        }
      } else {
        for (unsigned i = 1; i < NumRows; i++) {
          for (unsigned j = 1; j < NumCols; j++) {
            ScoreSystemType Diagonal = (Seq1[i - 1]==Seq2[j - 1]) ? (Matrix[(i - 1) * NumCols + j - 1] + Match) : Mismatch;
            ScoreSystemType Upper = Matrix[(i - 1) * NumCols + j] + Gap;
            ScoreSystemType Left = Matrix[i * NumCols + j - 1] + Gap;
            ScoreSystemType Score = std::max(std::max(Diagonal, Upper), Left);
            Matrix[i * NumCols + j] = Score;
            if (Score>=MaxScore) {
              MaxScore = Score;
              MaxRow=i;
              MaxCol=j;
            }
          }
        }
      }
    }
  }

  void buildResult() {

    auto &Seq1Ref = BaseType::getSequence(0);
    auto &Seq2Ref = BaseType::getSequence(1);
    Ty *Seq1 = &(Seq1Ref[0]);
    Ty *Seq2 = &(Seq2Ref[0]);

    auto &Data = BaseType::getResult().Data;

    ScoringSystem &Scoring = BaseType::getScoring();
    const ScoreSystemType Gap = Scoring.getGapPenalty();
    const ScoreSystemType Match = Scoring.getMatchProfit();
    const bool AllowMismatch = Scoring.getAllowMismatch();
    const ScoreSystemType Mismatch = AllowMismatch
                                   ?Scoring.getMismatchPenalty()
                                   :std::numeric_limits<ScoreSystemType>::min();

    int i = MatrixRows-1, j = MatrixCols-1;

    while (i>0 || j>0) {
      if (i>0 && j>0) {
       //Diagonal

       bool IsValidMatch = false;

     ScoreSystemType Score = std::numeric_limits<ScoreSystemType>::min();
     if (Matches) {
         IsValidMatch = Matches[(i - 1)*MatchesCols + j - 1];
       } else {
         IsValidMatch = (Seq1[i - 1]==Seq2[j - 1]);
       }

     if (AllowMismatch) {
       Score = Matrix[(i - 1) * MatrixCols + j - 1] + (IsValidMatch ? Match : Mismatch);
     } else {
       Score = IsValidMatch ? (Matrix[(i - 1) * MatrixCols + j - 1] + Match) : Mismatch;
     }

       if (Matrix[i*MatrixCols + j]==Score) {
         if (IsValidMatch || AllowMismatch) {
           Data.push_front(
               typename BaseType::EntryType(Seq1[i - 1], Seq2[j - 1], IsValidMatch));
         } else {
           Data.push_front(
               typename BaseType::EntryType(Seq1[i - 1], Blank, false));
           Data.push_front(
               typename BaseType::EntryType(Blank, Seq2[j - 1], false));
         }

         i--; j--;
       continue;
       }
      }
      if (i>0 && Matrix[i*MatrixCols + j]==(Matrix[(i-1)*MatrixCols + j]+Gap)) {
        //Up
       Data.push_front(
           typename BaseType::EntryType(Seq1[i - 1], Blank, false));
        i--;
      } else {
        //Left
       Data.push_front(
         typename BaseType::EntryType(Blank, Seq2[j - 1], false));
        j--;
      }
    }

  }

  void clearAll() {
    if (Matrix) delete[]Matrix;
    if (Matches) delete[]Matches;
    Matrix = nullptr;
    Matches = nullptr;
  }

public:
  static ScoringSystem getDefaultScoring() {
    return ScoringSystem(-1,2,-1);
  }

  NeedlemanWunschSA(
    ScoringSystem Scoring,
    MatchFnTy Match,
    ContainerType &Seq1,
    ContainerType &Seq2) : SequenceAligner<ContainerType,Ty,Blank,MatchFnTy>(Scoring, Match, Seq1, Seq2),
    Matrix(nullptr), Matches(nullptr) {

    cacheAllMatches();
    //printMatches();
    computeScoreMatrix();
    //printMatrix();
    buildResult();
    clearAll();        
  }

};
