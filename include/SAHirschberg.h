template<typename ContainerType, typename Ty=typename ContainerType::value_type, Ty Blank=Ty(0), typename MatchFnTy=std::function<bool(typename ContainerType::value_type,typename ContainerType::value_type)>>
class HirschbergSA : public SequenceAligner<ContainerType,Ty,Blank,MatchFnTy> {
private:
  ScoreSystemType *FinalScore;
  ScoreSystemType *ScoreAux;
  ScoreSystemType *ScoreCache;  

  using BaseType = SequenceAligner<ContainerType,Ty,Blank,MatchFnTy>;

  template<typename iterator1, typename iterator2>
  void NWScore(iterator1 Begin1, iterator1 End1, iterator2 Begin2, iterator2 End2) {
    const size_t SizeSeq1 = End1-Begin1;
    const size_t SizeSeq2 = End2-Begin2;

    ScoringSystem &Scoring = BaseType::getScoring();
    const ScoreSystemType Gap = Scoring.getGapPenalty();
    const ScoreSystemType Match = Scoring.getMatchProfit();
    const bool AllowMismatch = Scoring.getAllowMismatch();
    const ScoreSystemType Mismatch = AllowMismatch
                                   ?Scoring.getMismatchPenalty()
                                   :std::numeric_limits<ScoreSystemType>::min();

    FinalScore[0] = 0;
    for (size_t j = 1; j<=SizeSeq2; j++) {
      FinalScore[j] = FinalScore[j-1] + Gap; //Ins(F2[j-1]);
    }

    if (BaseType::getMatchOperation()==nullptr) {
      if (AllowMismatch) {
        for (size_t i = 1; i<=SizeSeq1; i++) {
          ScoreAux[0] = FinalScore[0] + Gap; //Del(*(Begin1+(i-1)));
          for (size_t j = 1; j<=SizeSeq2; j++) {
            ScoreSystemType Similarity = (*(Begin1+(i-1))==*(Begin2+(j-1))) ? Match : Mismatch;
            ScoreSystemType ScoreSub = FinalScore[j-1] + Similarity; //Sub(F1[i-1],F2[j-1]);
            ScoreSystemType ScoreDel = FinalScore[j]   + Gap; //Del(F1[i-1]);
            ScoreSystemType ScoreIns = ScoreAux[j-1] + Gap; //Ins(F2[j-1]);
            ScoreAux[j] = std::max(std::max(ScoreSub,ScoreDel),ScoreIns);
          }
          std::swap(FinalScore,ScoreAux);
        }
      } else {
        for (size_t i = 1; i<=SizeSeq1; i++) {
          ScoreAux[0] = FinalScore[0] + Gap; //Del(F1[i-1]);
          for (size_t j = 1; j<=SizeSeq2; j++) {
            ScoreSystemType ScoreSub = (*(Begin1+(i-1))==*(Begin2+(j-1))) ? (FinalScore[j-1] + Match) : Mismatch;
            ScoreSystemType ScoreDel = FinalScore[j]   + Gap; //Del(F1[i-1]);
            ScoreSystemType ScoreIns = ScoreAux[j-1] + Gap; //Ins(F2[j-1]);
            ScoreAux[j] = std::max(std::max(ScoreSub,ScoreDel),ScoreIns);
          }
          std::swap(FinalScore,ScoreAux);
        }
      }
    } else {
      if (AllowMismatch) {
        for (size_t i = 1; i<=SizeSeq1; i++) {
          ScoreAux[0] = FinalScore[0] + Gap; //Del(*(Begin1+(i-1)));
          for (size_t j = 1; j<=SizeSeq2; j++) {
            ScoreSystemType Similarity = BaseType::match(*(Begin1+(i-1)),*(Begin2+(j-1))) ? Match : Mismatch;
            ScoreSystemType ScoreSub = FinalScore[j-1] + Similarity; //Sub(F1[i-1],F2[j-1]);
            ScoreSystemType ScoreDel = FinalScore[j]   + Gap; //Del(F1[i-1]);
            ScoreSystemType ScoreIns = ScoreAux[j-1] + Gap; //Ins(F2[j-1]);
            ScoreAux[j] = std::max(std::max(ScoreSub,ScoreDel),ScoreIns);
          }
          std::swap(FinalScore,ScoreAux);
        }
      } else {
        for (size_t i = 1; i<=SizeSeq1; i++) {
          ScoreAux[0] = FinalScore[0] + Gap; //Del(F1[i-1]);
          for (size_t j = 1; j<=SizeSeq2; j++) {
            ScoreSystemType ScoreSub = BaseType::match(*(Begin1+(i-1)),*(Begin2+(j-1))) ? (FinalScore[j-1] + Match) : Mismatch;
            ScoreSystemType ScoreDel = FinalScore[j]   + Gap; //Del(F1[i-1]);
            ScoreSystemType ScoreIns = ScoreAux[j-1] + Gap; //Ins(F2[j-1]);
            ScoreAux[j] = std::max(std::max(ScoreSub,ScoreDel),ScoreIns);
          }
          std::swap(FinalScore,ScoreAux);
        }
      }
    }
    //last score is in FinalScore
  }

  template<typename ArrayType>
  void HirschbergRec(ArrayType &Seq1, ArrayType &Seq2, AlignedSequence<Ty> &Res) {
//  AlignedSequence<Ty> HirschbergRec(ArrayType &Seq1, ArrayType &Seq2, AlignedSequence<Ty> &Res) {
    //Ty Blank = BaseType::getBlank();
    //AlignedSequence<Ty> Res;
    if (Seq1.size()==0) {
      for (auto Char : Seq2) {
        Res.Data.push_back(typename BaseType::EntryType(Blank,Char,false));
      }
    } else if (Seq2.size()==0) {
      for (auto Char : Seq1) {
        Res.Data.push_back(typename BaseType::EntryType(Char,Blank,false));
      }
    } else if (Seq1.size()==1 || Seq2.size()==1) {
      NeedlemanWunschSA<ArrayView<ContainerType>, Ty, Blank, MatchFnTy> SA(BaseType::getScoring(),
                               BaseType::getMatchOperation(),
                               Seq1, Seq2);
      //Res.Data = std::move(SA.getResult().Data);
      Res.splice(SA.getResult());

    } else {
      int Seq1Mid = Seq1.size()/2;
   
      NWScore(Seq1.begin(),Seq1.begin()+Seq1Mid,Seq2.begin(), Seq2.end());
      std::swap(FinalScore,ScoreCache);
      
      ArrayType SlicedSeq1(Seq1);
      SlicedSeq1.sliceWindow(Seq1Mid,Seq1.size());
      NWScore(SlicedSeq1.rbegin(), SlicedSeq1.rend(), Seq2.rbegin(), Seq2.rend());

      size_t Seq2Mid = 0;
      int MaxScore = std::numeric_limits<int>::min();
      size_t Size2 = Seq2.size();
      for (size_t i = 0; i<Seq2.size(); i++) {
        int S = ScoreCache[i] + FinalScore[Size2-i];
        if (S>=MaxScore) {
          MaxScore = S;
          Seq2Mid = i;
        }
      }

      ArrayType NewSeq1L(Seq1);
      NewSeq1L.sliceWindow(0,Seq1Mid);
      ArrayType NewSeq2L(Seq2);
      NewSeq2L.sliceWindow(0,Seq2Mid);
      //AlignedSequence<Ty> ResL = HirschbergRec(NewSeq1L, NewSeq2L, Res);
      HirschbergRec(NewSeq1L, NewSeq2L, Res); 

      ArrayType NewSeq1R(Seq1);
      NewSeq1R.sliceWindow(Seq1Mid,Seq1.size());
      ArrayType NewSeq2R(Seq2);
      NewSeq2R.sliceWindow(Seq2Mid,Seq2.size());
      //AlignedSequence<Ty> ResR = HirschbergRec(NewSeq1R, NewSeq2R, Res);
      HirschbergRec(NewSeq1R, NewSeq2R, Res);


      //Res.splice(ResL);
      //Res.splice(ResR);
    }

    //return std::move(Res);
  }

public:
  HirschbergSA(
    ScoringSystem Scoring,
    MatchFnTy Match,
    ContainerType &Seq1,
    ContainerType &Seq2) : SequenceAligner<ContainerType,Ty,Blank,MatchFnTy>(Scoring, Match, Seq1, Seq2) {
    ScoreSystemType *ScoreContainer = new ScoreSystemType[3*(Seq2.size()+1)];
    FinalScore = &ScoreContainer[0];
    ScoreAux = &ScoreContainer[Seq2.size()+1];
    ScoreCache = &ScoreContainer[2*(Seq2.size()+1)];
    ArrayView< ContainerType > View1(Seq1);
    ArrayView< ContainerType > View2(Seq2);
    HirschbergRec(View1,View2,BaseType::getResult());
    //BaseType::getResult().Data = std::move(HirschbergRec(View1,View2).Data);
    delete []ScoreContainer;
  }

};
