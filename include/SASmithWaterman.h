template<typename ContainerType, typename Ty=typename ContainerType::value_type, Ty Blank=Ty(0), typename MatchFnTy=std::function<bool(Ty,Ty)>>
class SmithWatermanSA : public SequenceAligner<ContainerType, Ty, Blank, MatchFnTy> {
    private:
        //Matrix is 1D
        //Structured by first row then second row
        ScoreSystemType *Matrix;
        size_t MatrixRows;
        size_t MatrixCols;
        bool *Matches;
        size_t MatchesRows;
        size_t MatchesCols;

        //Enum for where to check for scoring in the similarity matrix
        enum Check { END, DIAGONAL, UP, LEFT};
        //Check check = End;

        size_t MaxRow;
        size_t MaxCol;
        ScoreSystemType MaxScore;

        using BaseType = SequenceAligner<ContainerType, Ty, Blank, MatchFnTy>;

        void cacheAllMatches(ContainerType &Seq1, ContainerType &Seq2) {
            if (BaseType::getMatchOperation() == nullptr) {
                Matches = nullptr;
                return;
            }

            const size_t SizeSeq1 = Seq1.size();
            const size_t SizeSeq2 = Seq2.size();

            MatchesRows = SizeSeq1;
            MatchesCols = SizeSeq2;
            Matches = new bool[SizeSeq1*SizeSeq2];

			MatrixRows = SizeSeq1 + 1;
			MatrixCols = SizeSeq2 + 1;
            for (unsigned i = 0; i < SizeSeq1; i++){
                for (unsigned j = 0; j < SizeSeq2; j++){
                    Matches[i*SizeSeq2 + j] = BaseType::match(Seq1[i], Seq2[j]);
                }
            }
        }

        void computeScoreMatrix(ContainerType &Seq1, ContainerType &Seq2) {
            //Initialise Matrix and scoring elements
            const size_t SizeSeq1 = Seq1.size();
            const size_t SizeSeq2 = Seq2.size();

            const size_t MatrixRows = SizeSeq1 + 1;
            const size_t MatrixCols = SizeSeq2 + 1;
            Matrix = new int[MatrixRows*MatrixCols];

            ScoringSystem &Scoring = BaseType::getScoring();
            const ScoreSystemType Gap = Scoring.getGapPenalty();
            const ScoreSystemType Match = Scoring.getMatchProfit();
            const bool AllowMismatch = Scoring.getAllowMismatch();
            const ScoreSystemType Mismatch = AllowMismatch
                                            ?Scoring.getMismatchPenalty()
                                            :std::numeric_limits<ScoreSystemType>::min();

            //Set up initial matrix scores that we know for sure
            //For Smith-Waterman we know for sure that the first row
            //And first column are straight 0s

            //Set first column to 0
            for (unsigned i = 0; i < MatrixRows; i++){
                Matrix[i*MatrixCols] = 0;
            }

            //Set first row to 0
            for (unsigned j = 0; j < MatrixCols; j++){
                Matrix[j] = 0;
            }

            //No if statement within for loop as this is inefficient and
            //when performed on large code bases can be exponentially slower
            //Compute Similarity Matrix for non trivial cells
            MaxScore = std::numeric_limits<ScoreSystemType>::min();
            //Check for nullptr / If a match was found
            if (Matches) {
                //If we allow for mismatches
                if (AllowMismatch) {
                    for (unsigned i = 1; i < MatrixRows; i++){
                        for (unsigned j = 1; j < MatrixCols; j++){

                            //If there is a match set Similarity to a match
                            ScoreSystemType Similarity = Matches[(i-1)*MatchesCols + j -1] ? Match : Mismatch;

                            //Calculate scores from diagonal, upper and left
                            ScoreSystemType Diagonal = Matrix[(i-1) * MatrixCols + j -1] + Similarity;
                            ScoreSystemType Upper = Matrix[(i-1) * MatrixCols + j] + Gap;
                            ScoreSystemType Left = Matrix[i * MatrixCols + j -1] + Gap;
                            ScoreSystemType Zero = 0;
                            ScoreSystemType Score = std::max({Diagonal, Upper, Left, Zero});

                            //Set the score in the similarity matrix
                            Matrix[i*MatrixCols + j] = Score;
                            
                            //Save the max score found and save it's position in the matrix
                            if (Score>=MaxScore) {
                                MaxScore = Score;
                                MaxRow = i;
                                MaxCol = j;
                            }

                        }
                    }
                }
                //Mismatch not allowed
                else {
                    for (unsigned i = 1; i < MatrixRows; i++){
                        for (unsigned j = 1; j < MatrixCols; j++){

                            //Calculate scores from diagonal, upper and left
                            ScoreSystemType Diagonal = Matches[(i-1) * MatchesCols + j -1] ? (Matrix[(i-1)*MatrixCols + j -1] + Match) : Mismatch;
                            ScoreSystemType Upper = Matrix[(i-1) * MatrixCols + j] + Gap;
                            ScoreSystemType Left = Matrix[i * MatrixCols + j -1] + Gap;
                            ScoreSystemType Zero = 0;
                            ScoreSystemType Score = std::max({Diagonal, Upper, Left, Zero});

                            //Set the score in the similarity matrix
                            Matrix[i*MatrixCols + j] = Score;
                            
                            //Save the max score found and save it's position in the matrix
                            if (Score>=MaxScore) {
                                MaxScore = Score;
                                MaxRow = i;
                                MaxCol = j;
                            }

                        }
                    }

                }
            }
            //If Matches is null / If a match was not found
            else {
                if (AllowMismatch) {
                    for (unsigned int i = 1; i < MatrixRows; i++){
                        for (unsigned j = 1; j < MatrixCols; j++){

                            //If the point in the sequences match then match otherwise mismatch
                            ScoreSystemType Similarity = (Seq1[i-1] == Seq2[j-1]) ? Match : Mismatch;

                            //Calculate scores from diagonal, upper and left
                            ScoreSystemType Diagonal = Matrix[(i-1) * MatrixCols + j -1] + Similarity;
                            ScoreSystemType Upper = Matrix[(i-1) * MatrixCols + j] + Gap;
                            ScoreSystemType Left = Matrix[i * MatrixCols + j -1] + Gap;
                            ScoreSystemType Zero = 0;
                            ScoreSystemType Score = std::max({Diagonal, Upper, Left, Zero});

                            //Set the score in the similarity matrix
                            Matrix[i*MatrixCols + j] = Score;

                            //Save the max score found and save it's position in the matrix
                            if (Score>=MaxScore) {
                                MaxScore = Score;
                                MaxRow = i;
                                MaxCol = j;
                            }

                        }
                    }
                }
                else {
                    for (unsigned int i = 0; i < MatrixRows; i++){
                        for (unsigned int j = 0; j < MatrixCols; j++){

                            //If the point in the sequences match then match otherwise mismatch
                            ScoreSystemType Diagonal = (Seq1[i-1] == Seq2[j-1]) ? (Matrix[(i-1) * MatrixCols + j - 1] + Match) : Mismatch;
                            ScoreSystemType Upper = Matrix[(i-1) * MatrixCols + j] + Gap;
                            ScoreSystemType Left = Matrix[i * MatrixCols + j -1] + Gap;
                            ScoreSystemType Zero = 0;
                            ScoreSystemType Score = std::max({Diagonal, Upper, Left, Zero});

                            //Set the score in the similarity matrix
                            Matrix[i*MatrixCols + j] = Score;
                            
                            //Save the max score found and save it's position in the matrix
                            if (Score>=MaxScore) {
                                MaxScore = Score;
                                MaxRow = i;
                                MaxCol = j;
                            }

                        }
                    }
                }

            }            
        }

        //Build the resulting aligned sequence
        void buildResult(ContainerType &Seq1, ContainerType &Seq2, AlignedSequence<Ty, Blank> &Result){

            auto &Data = Result.Data;

            ScoringSystem &Scoring = BaseType::getScoring();
            const ScoreSystemType Gap = Scoring.getGapPenalty();
            const ScoreSystemType Match = Scoring.getMatchProfit();
            const bool AllowMismatch = Scoring.getAllowMismatch();
            const ScoreSystemType Mismatch = AllowMismatch ? Scoring.getMismatchPenalty() : std::numeric_limits<ScoreSystemType>::min();

            int i = MaxRow, j = MaxCol;

            while (i>0 || j>0) {
				if (i == 0 || j == 0) {
					break;
				}
                if (i>0 && j>0){
                    //Diagonal
                    
                    bool IsValidMatch = false;

                    ScoreSystemType Score = std::numeric_limits<ScoreSystemType>::min();

                    if (Matches) {
                        IsValidMatch = Matches[(i-1)*MatchesCols + j - 1];
                    }
                    else {
                        IsValidMatch = (Seq1[i-1] == Seq2[j-1]);
                    }

                    if (AllowMismatch) {
                        //Score diagonal + if valid match found then the match score otherwise mismatch score
                        Score = Matrix[(i-1) * MatrixCols + j - 1] + (IsValidMatch ? Match : Mismatch);
						Score = std::max(Score, 0);
                    }
                    else {
                        Score = IsValidMatch ? (Matrix[(i-1) * MatrixCols + j - 1] + Match) : Mismatch;
						Score = std::max(Score, 0);
                    }

                    if (Matrix[i*MatrixCols + j] == Score) {
                        //If the end of the highest alignment has been reached
                        if (Score == 0) {
                            break;
                        }

                        if (IsValidMatch || AllowMismatch){
                            Data.push_front(
                                typename BaseType::EntryType(Seq1[i-1], Seq2[j-1], IsValidMatch)
                            );
                        }
                        else {
                            Data.push_front(
                                typename BaseType::EntryType(Seq1[i-1], Blank, false)
                            );
                            Data.push_front(
                                typename BaseType::EntryType(Blank, Seq2[j-1], false)
                            );
                        }

                        i--;
                        j--;
                        //Go to start of loop
                        continue;
                    }
                }

				if (i > 0 && Matrix[i * MatrixCols + j] == (Matrix[(i - 1) * MatrixCols + j] + Gap)) {

					//If the end of the highest alignment has been reached
					if ((Matrix[(i - 1) * MatrixCols + j] + Gap) <= 0) {
						break;
					}
					//Up
					Data.push_front(
						typename BaseType::EntryType(Seq1[i - 1], Blank, false)
					);
					i--;
				}
				else {

					//If the end of the highest alignment has been reached
					if ((Matrix[(i)*MatrixCols + j - 1] + Gap) <= 0) {
						break;
					}

					//Left
					Data.push_front(
						typename BaseType::EntryType(Blank, Seq2[j - 1], false)
					);
					j--;
				}
            }
        }

        void clearAll(){
            if (Matrix) delete[]Matrix;
            if (Matches) delete[]Matches;
            Matrix = nullptr;
            Matches = nullptr;
        }
    
    

    public:
        static ScoringSystem getDefaultScoring(){
            return ScoringSystem(-1,1,-1);
        }

        SmithWatermanSA() : BaseType(getDefaultScoring(), nullptr), Matrix(nullptr), Matches(nullptr) {}

        SmithWatermanSA(ScoringSystem Scoring, MatchFnTy Match = nullptr)
         : BaseType(Scoring, Match),
           Matrix(nullptr), Matches(nullptr) {}

        virtual AlignedSequence<Ty, Blank> getAlignment(ContainerType &Seq1, ContainerType &Seq2){
            AlignedSequence<Ty, Blank> Result;
            cacheAllMatches(Seq1, Seq2);
            computeScoreMatrix(Seq1, Seq2);
            buildResult(Seq1, Seq2, Result);
            clearAll();
            return Result;
        }

};