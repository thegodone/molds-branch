//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   //
// Copyright (C) 2012-2014 Michihiro Okuyama
//                                                                        // 
// This file is part of MolDS.                                            // 
//                                                                        // 
// MolDS is free software: you can redistribute it and/or modify          // 
// it under the terms of the GNU General Public License as published by   // 
// the Free Software Foundation, either version 3 of the License, or      // 
// (at your option) any later version.                                    // 
//                                                                        // 
// MolDS is distributed in the hope that it will be useful,               // 
// but WITHOUT ANY WARRANTY; without even the implied warranty of         // 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          // 
// GNU General Public License for more details.                           // 
//                                                                        // 
// You should have received a copy of the GNU General Public License      // 
// along with MolDS.  If not, see <http://www.gnu.org/licenses/>.         // 
//************************************************************************//
#ifndef INCLUDED_ZINDOS
#define INCLUDED_ZINDOS
namespace MolDS_zindo{

/***
 *  Main Reference for Zindo is [RZ_1973]
 */
class ZindoS : public MolDS_cndo::Cndo2{
public:
   ZindoS();
   virtual ~ZindoS();
   virtual void SetMolecule(MolDS_base::Molecule* molecule);
   void DoCIS();
   void OutputCISResults() const;
   void CalcOverlapSingletSDsWithAnotherElectronicStructure(double** overlapSingletSDs, 
                                                            double const* const* overlapMOs) const;
   void CalcOverlapESsWithAnotherElectronicStructure(double** overlapESs, 
                                                     double const* const* overlapSingletSDs,
                                                     const MolDS_base::ElectronicStructure& lhsElectronicStructure) const;
protected:
   std::string errorMessageDavidsonNotConverged;
   std::string errorMessageCalcCISMatrix;
   std::string errorMessageCalcZMatrixForceEtaNull;
   std::string messageStartCIS;
   std::string messageDoneCIS;
   std::string messageDavidsonConverge;
   std::string messageStartCalcCISMatrix;
   std::string messageOmpElapsedTimeCalcCISMarix;
   std::string messageOmpElapsedTimeCIS;
   std::string messageDoneCalcCISMatrix;
   double*** zMatrixForce;
   double*** etaMatrixForce;
   struct MoIndexPair{int moI; int moJ; bool isMoICIMO; bool isMoJCIMO;};
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual void CalcCISProperties();
   virtual void CalcElectronicTransitionDipoleMoment(double* transitionDipoleMoment,
                                                     int to, int from,
                                                     double const* const* fockMatrix,
                                                     double const* const* matrixCIS,
                                                     double const* const* const* cartesianMatrix,
                                                     const MolDS_base::Molecule& molecule, 
                                                     double const* const* orbitalElectronPopulation,
                                                     double const* const* overlapAOs,
                                                     double const* groundStateDipole) const;
   virtual void CalcGammaAB(double** gammaAB, const MolDS_base::Molecule& molecule) const;
   virtual double GetFockDiagElement(const MolDS_base_atoms::Atom& atomA, 
                                     int indexAtomA, 
                                     int mu, 
                                     const MolDS_base::Molecule& molecule, 
                                     double const* const* gammaAB,
                                     double const* const* orbitalElectronPopulation, 
                                     double const* atomicElectronPopulation,
                                     double const* const* const* const* const* const* twoElecsTwoAtomCores,
                                     bool isGuess) const;
   virtual double GetFockOffDiagElement(const MolDS_base_atoms::Atom& atomA, 
                                        const MolDS_base_atoms::Atom& atomB, 
                                        int indexAtomA, 
                                        int indexAtomB, 
                                        int mu, int nu, 
                                        const MolDS_base::Molecule& molecule, 
                                        double const* const* gammaAB, 
                                        double const* const* overelap,
                                        double const* const* orbitalElectronPopulation, 
                                        double const* const* const* const* const* const* twoElecsTwoAtomCores,
                                        bool isGuess) const;
   virtual void CalcDiatomicOverlapAOsInDiatomicFrame(double** diatomicOverlapAOs, 
                                                      const MolDS_base_atoms::Atom& atomA, 
                                                      const MolDS_base_atoms::Atom& atomB) const;
   virtual void CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(double** diatomicOverlapAOsDeri, 
                                                                   const MolDS_base_atoms::Atom& atomA, 
                                                                   const MolDS_base_atoms::Atom& atomB) const;
   virtual void CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(double** diatomicOverlapAOs2ndDeri, 
                                                                   const MolDS_base_atoms::Atom& atomA, 
                                                                   const MolDS_base_atoms::Atom& atomB) const;
   virtual double GetCoulombInt(MolDS_base::OrbitalType orbital1, 
                                MolDS_base::OrbitalType orbital2, 
                                const MolDS_base_atoms::Atom& atom) const; // Apendix in [BZ_1979]
   virtual double GetExchangeInt(MolDS_base::OrbitalType orbital1, 
                                 MolDS_base::OrbitalType orbital2, 
                                 const MolDS_base_atoms::Atom& atom) const; // Apendix in [BZ_1979]
   virtual void CalcTwoElecsTwoCores(double****** twoElecsTwoAtomCores, 
                                     double****** twoElecsAtomEpcCores,
                                     const MolDS_base::Molecule& molecule) const;
   virtual double GetMolecularIntegralElement(int moI, 
                                              int moJ, 
                                              int moK, 
                                              int moL, 
                                              const MolDS_base::Molecule& molecule, 
                                              double const* const* fockMatrix, 
                                              double const* const* gammaAB) const;
   virtual void CalcCISMatrix(double** matrixCIS) const;
   double GetCISDiagElement(double const* energiesMO,
                            double const* const* const* const* nishimotoMatagaMatrix,
                            const MolDS_base::Molecule& molecule,
                            double const* const* fockMatrix, 
                            int moI,
                            int moA) const;
   double GetCISOffDiagElement(double const* const* const* const* nishimotoMatagaMatrix,
                               const MolDS_base::Molecule& molecule,
                               double const* const* fockMatrix, 
                               int moI,
                               int moA,
                               int moJ,
                               int moB) const;
   bool RequiresExcitedStatesForce(const std::vector<int>& elecStates) const;
   virtual void CalcForce(const std::vector<int>& elecStates);
   int GetSlaterDeterminantIndex(int activeOccIndex, int activeVirIndex) const;
   int GetActiveOccIndex(const MolDS_base::Molecule& molecule, int matrixCISIndex) const;
   int GetActiveVirIndex(const MolDS_base::Molecule& molecule, int matrixCISIndex) const;
   void CheckMatrixForce(const std::vector<int>& elecStates);
   void CalcEtaMatrixForce(const std::vector<int>& elecStates);
   void CalcZMatrixForce(const std::vector<int>& elecStates);
   void CalcActiveSetVariablesQ(std::vector<MoIndexPair>* nonRedundantQIndeces, 
                                std::vector<MoIndexPair>* redundantQIndeces,
                                int numberActiveOcc,
                                int numberActiveVir) const;
   virtual double GetSmallQElement(int moI, 
                                   int moP, 
                                   double const* const* xiOcc, 
                                   double const* const* xiVir,
                                   double const* const* eta) const;
   double GetGammaNRElement(int moI, int moJ, int moK, int moL) const;
   double GetGammaRElement (int moI, int moJ, int moK, int moL) const;
   double GetNNRElement    (int moI, int moJ, int moK, int moL) const;
   double GetNRElement     (int moI, int moJ, int moK, int moL) const;
   double GetKNRElement    (int moI, int moJ, int moK, int moL) const;
   double GetKRElement     (int moI, int moJ, int moK, int moL) const;
   virtual double GetAuxiliaryKNRKRElement(int moI, int moJ, int moK, int moL) const;
   void CalcForceExcitedOverlapAOsPart(double* force, 
                                       int elecStateIndex,
                                       int indexAtomA,
                                       int indexAtomB,
                                       double const* const* const* diatomicOverlapAOs1stDerivs) const;
private:
   std::string errorMessageElecState;
   std::string errorMessageNishimotoMataga;
   std::string errorMessageDavidsonMaxIter;
   std::string errorMessageDavidsonMaxDim;
   std::string messageStartDirectCIS;
   std::string messageDoneDirectCIS;
   std::string messageStartDavidsonCIS;
   std::string messageDoneDavidsonCIS;
   std::string messageNumIterCIS;
   std::string messageResidualNorm;
   std::string messageDavidsonReachCISMatrix;
   std::string messageDavidsonGoToDirect;
   std::string messageExcitedStatesEnergies;
   std::string messageExcitedStatesEnergiesTitle;
   std::string messageExcitonEnergiesCIS;
   std::string messageExcitonEnergiesShortCIS;
   std::string messageExcitonEnergiesCISTitle;
   std::string messageTotalDipoleMomentsTitle;
   std::string messageTotalDipoleMoment;
   std::string messageElectronicDipoleMomentsTitle;
   std::string messageElectronicDipoleMoment;
   std::string messageTransitionDipoleMomentsTitle;
   std::string messageTransitionDipoleMoment;
   double**** nishimotoMatagaMatrix;
   int    matrixForceElecStatesNum;
   double nishimotoMatagaParamA;
   double nishimotoMatagaParamB;
   double overlapAOsCorrectionSigma;
   double overlapAOsCorrectionPi;
   int    zMatrixForceElecStatesNum;
   int    etaMatrixForceElecStatesNum;
   void DoCISDirect();
   void DoCISDavidson();
   void OutputCISDipole() const;
   void OutputCISTransitionDipole() const;
   void OutputCISMulliken() const;
   void OutputCISUnpairedPop() const;
   void CalcFreeExcitonEnergies(double** freeExcitonEnergiesCIS, 
                                const MolDS_base::Molecule& molecule, 
                                double const* energiesMO, 
                                double const* const* matrixCIS,
                                int matrixCISdimension) const;
   void CalcOrbitalElectronPopulationCIS(double**** orbitalElectronPopulationCIS, 
                                         double const* const* orbitalElectronPopulation, 
                                         const MolDS_base::Molecule& molecule, 
                                         double const* const* fockMatrix,
                                         double const* const* matrixCIS) const;
   void CalcAtomicElectronPopulationCIS(double*** atomicElectronPopulationCIS,
                                        double const* const* const* orbitalElectronPopulationCIS, 
                                        const MolDS_base::Molecule& molecule) const;
   void CalcAtomicUnpairedPopulationCIS(double*** atomicUnpairedPopulationCIS,
                                        double const* const* const* orbitalElectronPopulationCIS, 
                                        const MolDS_base::Molecule& molecule) const; 
   void CalcElectronicDipoleMomentsExcitedStates(double*** electronicTransitionDipoleMoments,
                                                 double const* const* fockMatrix,
                                                 double const* const* matrixCIS,
                                                 double const* const* const* cartesianMatrix,
                                                 const MolDS_base::Molecule& molecule, 
                                                 double const* const* orbitalElectronPopulation,
                                                 double const* const* overlapAOs) const;
   void CalcElectronicTransitionDipoleMoments(double*** electronicTransitionDipoleMoments,
                                              double const* const* fockMatrix,
                                              double const* const* matrixCIS,
                                              double const* const* const* cartesianMatrix,
                                              const MolDS_base::Molecule& molecule, 
                                              double const* const* orbitalElectronPopulation,
                                              double const* const* overlapAOs) const;
   double GetNishimotoMatagaTwoEleInt(const MolDS_base_atoms::Atom& atomA, 
                                      MolDS_base::OrbitalType orbitalA, 
                                      const MolDS_base_atoms::Atom& atomB, 
                                      MolDS_base::OrbitalType orbitalB) const; // ref. [MN_1957] and (5a) in [AEZ_1986]
   double GetNishimotoMatagaTwoEleInt(const MolDS_base_atoms::Atom& atomA, 
                                      MolDS_base::OrbitalType orbitalA, 
                                      const MolDS_base_atoms::Atom& atomB, 
                                      MolDS_base::OrbitalType orbitalB,
                                      const double rAB) const; // ref. [MN_1957] and (5a) in [AEZ_1986]
   double GetNishimotoMatagaTwoEleInt1stDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                   MolDS_base::OrbitalType orbitalA, 
                                                   const MolDS_base_atoms::Atom& atomB, 
                                                   MolDS_base::OrbitalType orbitalB,
                                                   MolDS_base::CartesianType axisA) const;// ref. [MN_1957] and (5a) in [AEZ_1986]
   double GetNishimotoMatagaTwoEleInt1stDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                   MolDS_base::OrbitalType orbitalA, 
                                                   const MolDS_base_atoms::Atom& atomB, 
                                                   MolDS_base::OrbitalType orbitalB,
                                                   const double rAB,
                                                   MolDS_base::CartesianType axisA) const;// ref. [MN_1957] and (5a) in [AEZ_1986]
   void CalcNishimotoMatagaMatrix(double**** nishimotoMatagaMatrix, 
                                  const MolDS_base::Molecule& molecule) const;
   void CalcRitzVector(double* ritzVector, 
                       double const* const* expansionVectors, 
                       double const* const* interactionMatrix, 
                       int interactionMatrixDimension, 
                       int ritzVectorIndex) const;
   void CalcResidualVectorAndNorm(double* residualVector, 
                                  double* norm, 
                                  double const* ritzVector, 
                                  double const* interactionEigenEnergies, 
                                  int residualVectorIndex) const;
   void SortCISEigenVectorCoefficients(std::vector<CISEigenVectorCoefficient>* cisEigenVectorCoefficients,
                                       double* cisEigenVector) const;
   void SortSingleExcitationSlaterDeterminants(std::vector<MoEnergyGap>* moEnergyGaps) const;
   void UpdateExpansionVectors(double** expansionVectors, 
                               int* notConvergedStates, 
                               double const* interactionEigenEnergies, 
                               double const* residualVector,
                               int interactionMatrixDimension, 
                               int residualVectorIndex) const;
   void CalcInteractionMatrix(double** interactionMatrix, 
                              double const* const* expansionVectors, 
                              int interactionMatrixDimension) const;
   void FreeDavidsonCISTemporaryMtrices(double*** expansionVectors, 
                                        double** residualVector, 
                                        double** ritzVector) const;
   void FreeDavidsonRoopCISTemporaryMtrices(double*** interactionMatrix, 
                                            int interactionMatrixDimension, 
                                            double** interactionEigenEnergies) const;
   void CalcDiatomicTwoElecsTwoCores1stDerivatives(double*** matrix, 
                                                   int indexAtomA, 
                                                   int indexAtomB) const;
   void MallocTempMatricesCalcForce(double**** diatomicOverlapAOs1stDerivs, 
                                    double**** diatomicTwoElecsTwoCores1stDerivs,
                                    double***  tmpDiaOverlapAOsInDiaFrame,       
                                    double***  tmpDiaOverlapAOs1stDerivInDiaFrame,
                                    double***  tmpRotMat,
                                    double***  tmpRotMat1stDeriv,
                                    double**** tmpRotMat1stDerivs,
                                    double***  tmpRotatedDiatomicOverlap,
                                    double**   tmpRotatedDiatomicOverlapVec,
                                    double***  tmpMatrixBC,
                                    double**   tmpVectorBC) const;         
   void FreeTempMatricesCalcForce(double**** diatomicOverlapAOs1stDerivs, 
                                  double**** diatomicTwoElecsTwoCores1stDerivs,
                                  double***  tmpDiaOverlapAOsInDiaFrame,       
                                  double***  tmpDiaOverlapAOs1stDerivInDiaFrame,
                                  double***  tmpRotMat,
                                  double***  tmpRotMat1stDeriv,
                                  double**** tmpRotMat1stDerivs,
                                  double***  tmpRotatedDiatomicOverlap,
                                  double**   tmpRotatedDiatomicOverlapVec,
                                  double***  tmpMatrixBC,
                                  double**   tmpVectorBC) const;         
   void CalcForceExcitedStaticPart(double* force, 
                                   int elecStateIndex,
                                   int indexAtomA,
                                   int indexAtomB,
                                   double const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
   void CalcForceExcitedElecCoreAttractionPart(double* force, 
                                               int elecStateIndex,
                                               int indexAtomA,
                                               int indexAtomB,
                                               double const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
   void CalcForceExcitedTwoElecPart(double* force, 
                                    int elecStateIndex,
                                    int indexAtomA,
                                    int indexAtomB,
                                    double const* const* const* diatomicTwoElecsTwoCores1stDerivs) const;
   void CheckZMatrixForce(const std::vector<int>& elecStates);
   void CheckEtaMatrixForce(const std::vector<int>& elecStates);
   double GetZMatrixForceElement(double const* y,
                                 double const* q,
                                 double const* const* transposedFockMatrix,
                                 const std::vector<MoIndexPair>& nonRedundantQIndeces,
                                 const std::vector<MoIndexPair>& redundantQIndeces,
                                 int mu, 
                                 int nu) const;
   void MallocTempMatrixForZMatrix(double** delta,
                                   double** q,
                                   double*** gammaNRMinusKNR, 
                                   double*** kRDag,
                                   double** y,
                                   double*** transposedFockMatrix,
                                   double*** xiOcc,
                                   double*** xiVir,
                                   int sizeQNR,
                                   int sizeQR) const;
   void FreeTempMatrixForZMatrix(double** delta,
                                 double** q,
                                 double*** gammaNRMinusKNR, 
                                 double*** kRDag,
                                 double** y,
                                 double*** transposedFockMatrix,
                                 double*** xiOcc,
                                 double*** xiVir,
                                 int sizeQNR,
                                 int sizeQR) const;
   void CalcDeltaVector(double* delta, int exciteState) const;
   void CalcQVector(double* q, 
                    double const* delta, 
                    double const* const* xiOcc,
                    double const* const* xiVir,
                    double const* const* eta,
                    const std::vector<MoIndexPair>& nonRedundantQIndeces,
                    const std::vector<MoIndexPair>& redundantQIndeces) const;
   void CalcXiMatrices(double** xiOcc, 
                       double** xiVir, 
                       int exciteState,
                       double const* const* transposedFockMatrix) const;
   void CalcAuxiliaryVector(double* y,
                            double const* q,
                            double const* const* kRDagerGammaRInv,
                            const std::vector<MoIndexPair>& nonRedundantQIndeces,
                            const std::vector<MoIndexPair>& redundantQIndeces) const;
   double GetKRDagerElement(int moI, int moJ, int moK, int moL) const;
   void CalcGammaNRMinusKNRMatrix(double** gammaNRMinusKNR, 
                                  const std::vector<MoIndexPair>& nonRedundantQIndeces) const;
   void CalcKRDagerGammaRInvMatrix(double** kRDagerGammaRInv, 
                                   const std::vector<MoIndexPair>& nonRedundantQIndeces,
                                   const std::vector<MoIndexPair>& redundantQIndeces) const;
};

}
#endif



