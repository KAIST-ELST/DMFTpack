#ifndef model_INCLUDE
#define model_INCLUDE

extern int * isOrbitalCorr,  *isOrbitalHartr, *CorrIndex, *HartrIndex;
extern int *CorrIndex, *HartrIndex;
void setCorrelatedSpaceIndex( std::vector<int> HartreeOrbital_idx , int NumCorrAtom   );
extern std::vector<Eigen::MatrixXcd> DF_CorrBase;
extern double muDFT;


#endif
