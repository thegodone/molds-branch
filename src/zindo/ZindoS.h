//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
// Copyright (C) 2012-2013 Michihiro Okuyama
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
 *  Main Refference for Zindo is [RZ_1973]
 */
class ZindoS : public MolDS_cndo::Cndo2{
public:
   ZindoS();
   virtual ~ZindoS();
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
   std::string messageStartCIS;
   std::string messageDoneCIS;
   std::string messageDavidsonConverge;
   std::string messageStartCalcCISMatrix;
   std::string messageOmpElapsedTimeCalcCISMarix;
   std::string messageOmpElapsedTimeCIS;
   std::string messageDoneCalcCISMatrix;
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual void CalcCISProperties();
   virtual double GetElectronicTransitionDipoleMoment(int to, int from, MolDS_base::CartesianType axis,
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
                                     double const* const* const* const* const* const* twoElecTwoCore,
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
                                        double const* const* const* const* const* const* twoElecTwoCore,
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
   virtual double GetMolecularIntegralElement(int moI, 
                                              int moJ, 
                                              int moK, 
                                              int moL, 
                                              const MolDS_base::Molecule& molecule, 
                                              double const* const* fockMatrix, 
                                              double const* const* gammaAB) const;
   virtual void CalcCISMatrix(double** matrixCIS) const;
   virtual void CalcForce(const std::vector<int>& elecStates);
   int GetSlaterDeterminantIndex(int activeOccIndex, int activeVirIndex) const;
   int GetActiveOccIndex(const MolDS_base::Molecule& molecule, int matrixCISIndex) const;
   int GetActiveVirIndex(const MolDS_base::Molecule& molecule, int matrixCISIndex) const;
   void CheckMatrixForce(const std::vector<int>& elecStates);
private:
   std::string errorMessageCalcForceNotGroundState;
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
   int    matrixForceElecStatesNum;
   double nishimotoMatagaParamA;
   double nishimotoMatagaParamB;
   double overlapAOsCorrectionSigma;
   double overlapAOsCorrectionPi;
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
   void CalcElectronicDipoleMomentsExcitedState(double*** electronicTransitionDipoleMoments,
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
   double GetNishimotoMatagaTwoEleInt1stDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                   MolDS_base::OrbitalType orbitalA, 
                                                   const MolDS_base_atoms::Atom& atomB, 
                                                   MolDS_base::OrbitalType orbitalB,
                                                   MolDS_base::CartesianType axisA) const;// ref. [MN_1957] and (5a) in [AEZ_1986]
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
};

}
#endif



