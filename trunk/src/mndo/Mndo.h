//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
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
#ifndef INCLUDED_MNDO
#define INCLUDED_MNDO
namespace MolDS_mndo{

/***
 *  Main Refferences for MNDO are [DT_1977, DT_1977-2, DT_1977-3]
 */
class Mndo : public MolDS_zindo::ZindoS{
public:
   Mndo();
   virtual ~Mndo();
   virtual void SetMolecule(MolDS_base::Molecule* molecule);
   virtual void OutputSCFResults() const;
protected:
   std::string errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles;
   std::string errorMessageGetSemiEmpiricalMultipoleInteraction1stDeriBadMultipoles;
   std::string errorMessageGetSemiEmpiricalMultipoleInteraction2ndDeriBadMultipoles;
   std::string errorMessageGetNddoRepulsionIntegral;
   std::string errorMessageGetNddoRepulsionIntegral1stDerivative;
   std::string errorMessageGetNddoRepulsionIntegral2ndDerivative;
   std::string errorMessageCalcDiatomicTwoElecTwoCoreNullMatrix;
   std::string errorMessageCalcTwoElecTwoCoreNullMatrix;
   std::string errorMessageCalcDiatomicTwoElecTwoCoreSameAtoms;
   std::string errorMessageCalcDiatomicTwoElecTwoCore1stDerivativesSameAtoms;
   std::string errorMessageCalcDiatomicTwoElecTwoCore2ndDerivativesSameAtoms;
   std::string errorMessageCalcDiatomicTwoElecTwoCore1stDerivativesNullMatrix;
   std::string errorMessageCalcDiatomicTwoElecTwoCore2ndDerivativesNullMatrix;
   std::string errorMessageCalcZMatrixForceEtaNull;
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual void CalcSCFProperties();
   virtual void CalcNormalModes(double** normalModes, double* normalForceConstants, const MolDS_base::Molecule& molecule) const;
   virtual void CalcForce(const std::vector<int>& elecStates);
   virtual double GetDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB) const;
   virtual double GetDiatomCoreRepulsion1stDerivative(int indexAtomA,
                                                      int indexAtomB, 
                                                      MolDS_base::CartesianType axisA) const;
   virtual double GetDiatomCoreRepulsion2ndDerivative(int indexAtomA,
                                                      int indexAtomB, 
                                                      MolDS_base::CartesianType axisA1,
                                                      MolDS_base::CartesianType axisA2) const;
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
   virtual void CalcDiatomicOverlapInDiatomicFrame(double** diatomicOverlap, 
                                                   const MolDS_base_atoms::Atom& atomA, 
                                                   const MolDS_base_atoms::Atom& atomB) const;
   virtual void CalcDiatomicOverlap1stDerivativeInDiatomicFrame(double** diatomicOverlapDeri, 
                                                                const MolDS_base_atoms::Atom& atomA, 
                                                                const MolDS_base_atoms::Atom& atomB) const;
   virtual void CalcDiatomicOverlap2ndDerivativeInDiatomicFrame(double** diatomicOverlap2ndDeri, 
                                                                const MolDS_base_atoms::Atom& atomA, 
                                                                const MolDS_base_atoms::Atom& atomB) const;
   virtual double GetCoulombInt(MolDS_base::OrbitalType orbital1, 
                                MolDS_base::OrbitalType orbital2, 
                                const MolDS_base_atoms::Atom& atom) const; 
   virtual double GetExchangeInt(MolDS_base::OrbitalType orbital1, 
                                 MolDS_base::OrbitalType orbital2, 
                                 const MolDS_base_atoms::Atom& atom) const; 
   virtual void CalcTwoElecTwoCore(double****** twoElecTwoCore, 
                                   const MolDS_base::Molecule& molecule) const;
   virtual double GetMolecularIntegralElement(int moI, 
                                              int moJ, 
                                              int moK, 
                                              int moL, 
                                              const MolDS_base::Molecule& molecule, 
                                              double const* const* fockMatrix, 
                                              double const* const* gammaAB) const;
   virtual void CalcCISMatrix(double** matrixCIS) const;
private:
   std::string errorMessageMultipoleA;
   std::string errorMessageMultipoleB;
   std::string messageHeatsFormation;
   std::string messageHeatsFormationTitle;
   struct MoIndexPair{int moI; int moJ; bool isMoICIMO; bool isMoJCIMO;};
   double    heatsFormation;
   int       zMatrixForceElecStatesNum;
   int       etaMatrixForceElecStatesNum;
   double*** zMatrixForce;
   double*** etaMatrixForce;
   double GetAuxiliaryDiatomCoreRepulsionEnergy(const MolDS_base_atoms::Atom& atomA,
                                                const MolDS_base_atoms::Atom& atomB,
                                                double distanceAB) const;
   double GetAuxiliaryDiatomCoreRepulsionEnergy1stDerivative(const MolDS_base_atoms::Atom& atomA,
                                                             const MolDS_base_atoms::Atom& atomB,
                                                             double distanceAB,
                                                             MolDS_base::CartesianType axisA) const;
   double GetAuxiliaryDiatomCoreRepulsionEnergy2ndDerivative(const MolDS_base_atoms::Atom& atomA,
                                                             const MolDS_base_atoms::Atom& atomB,
                                                             double distanceAB,
                                                             MolDS_base::CartesianType axisA1,
                                                             MolDS_base::CartesianType axisA2) const;
   double GetGammaNRElement(int moI, int moJ, int moK, int moL) const;
   double GetGammaRElement(int moI, int moJ, int moK, int moL) const;
   double GetNNRElement(int moI, int moJ, int moK, int moL) const;
   double GetNRElement(int moI, int moJ, int moK, int moL) const;
   double GetKNRElement(int moI, int moJ, int moK, int moL) const;
   double GetKRElement(int moI, int moJ, int moK, int moL) const;
   double GetKRDagerElement(int moI, int moJ, int moK, int moL) const;
   double GetAuxiliaryKNRKRElement(int moI, int moJ, int moK, int moL) const;
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
   double GetSmallQElement(int moI, 
                           int moP, 
                           double const* const* xiOcc, 
                           double const* const* xiVir,
                           double const* const* eta) const;
   void CalcQVector(double* q, 
                    double const* delta, 
                    double const* const* xiOcc,
                    double const* const* xiVir,
                    double const* const* eta,
                    const std::vector<MoIndexPair>& nonRedundantQIndeces,
                    const std::vector<MoIndexPair>& redundantQIndeces) const;
   void TransposeFockMatrixMatrix(double** transposedFockMatrix) const;
   void CalcGammaNRMinusKNRMatrix(double** gammaNRMinusKNR, 
                                  const std::vector<MoIndexPair>& nonRedundantQIndeces) const;
   void CalcKRDagerGammaRInvMatrix(double** kRDagerGammaRInv, 
                                   const std::vector<MoIndexPair>& nonRedundantQIndeces,
                                   const std::vector<MoIndexPair>& redundantQIndeces) const;
   void CalcAuxiliaryVector(double* y,
                            double const* q,
                            double const* const* kRDagerGammaRInv,
                            const std::vector<MoIndexPair>& nonRedundantQIndeces,
                            const std::vector<MoIndexPair>& redundantQIndeces) const;
   void CalcXiMatrices(double** xiOcc, 
                       double** xiVir, 
                       int exciteState,
                       double const* const* transposedFockMatrix) const;
   double GetZMatrixForceElement(double const* y,
                                 double const* q,
                                 double const* const* transposedFockMatrix,
                                 const std::vector<MoIndexPair>& nonRedundantQIndeces,
                                 const std::vector<MoIndexPair>& redundantQIndeces,
                                 int mu, 
                                 int nu) const;
   void CheckZMatrixForce(const std::vector<int>& elecStates);
   void CheckEtaMatrixForce(const std::vector<int>& elecStates);
   void CalcZMatrixForce(const std::vector<int>& elecStates);
   void CalcEtaMatrixForce(const std::vector<int>& elecStates);
   bool RequiresExcitedStatesForce(const std::vector<int>& elecStates) const;
   double GetCISCoefficientMOEnergy(int k, int l, int r, int numberActiveVir) const;
   double GetCISCoefficientTwoElecIntegral(int k, int l, int p, int q, int r, int s, int numberActiveVir) const;
   void CalcActiveSetVariablesQ(std::vector<MoIndexPair>* nonRedundantQIndeces, 
                                std::vector<MoIndexPair>* redundantQIndeces,
                                int numberActiveOcc,
                                int numberActiveVir) const;
   void CalcHessianSCF(double** hessianSCF, bool isMassWeighted) const;
   double GetHessianElementSameAtomsSCF(int indexAtomA, 
                                        MolDS_base::CartesianType axisA1,
                                        MolDS_base::CartesianType axisA2,
                                        double const* const*               orbitalElectronPopulation,
                                        double const* const* const* const* orbitalElectronPopulation1stDerivs,
                                        double const* const* const* const*                      diatomicOverlap1stDerivs,
                                        double const* const* const* const* const*               diatomicOverlap2ndDerivs,
                                        double const* const* const* const* const* const*        diatomicTwoElecTwoCore1stDerivs,
                                        double const* const* const* const* const* const* const* diatomicTwoElecTwoCore2ndDerivs) const;
   double GetHessianElementDifferentAtomsSCF(int indexAtomA, 
                                             int indexAtomB,
                                             MolDS_base::CartesianType axisA,
                                             MolDS_base::CartesianType axisB,
                                             double const* const*               orbitalElectronPopulation,
                                             double const* const* const* const* orbitalElectronPopulation1stDerivs,
                                             double const* const* const* const*                      diatomicOverlap1stDerivs,
                                             double const* const* const* const* const*               diatomicOverlap2ndDerivs,
                                             double const* const* const* const* const* const*        diatomicTwoElecTwoCore1stDerivs,
                                             double const* const* const* const* const* const* const* diatomicTwoElecTwoCore2ndDerivs) const;
   void MallocTempMatricesEachThreadCalcHessianSCF(double*****    diatomicOverlap1stDerivs,
                                                   double******   diatomicOverlap2ndDerivs,
                                                   double*******  diatomicTwoElecTwoCore1stDerivs,
                                                   double******** diatomicTwoElecTwoCore2ndDerivs) const;
   void FreeTempMatricesEachThreadCalcHessianSCF(double*****    diatomicOverlap1stDerivs,
                                                 double******   diatomicOverlap2ndDerivs,
                                                 double*******  diatomicTwoElecTwoCore1stDeriv,
                                                 double******** diatomicTwoElecTwoCore2ndDerivs) const;
   double GetAuxiliaryHessianElement1(int mu, 
                                      int nu, 
                                      int indexAtomA,
                                      int indexAtomC,
                                      MolDS_base::CartesianType axisA1,
                                      MolDS_base::CartesianType axisA2,
                                      double const* const* orbitalElectronPopulation,
                                      double const* const* const* const* const* const* diatomicTwoElecTwoCore2ndDerivs) const;
   double GetAuxiliaryHessianElement2(int mu, 
                                      int nu, 
                                      int indexAtomA,
                                      int indexAtomB,
                                      int indexAtomC,
                                      MolDS_base::CartesianType axisA,
                                      MolDS_base::CartesianType axisB,
                                      double const* const* const* const* orbitalElectronPopulation1stDerivs,
                                      double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const;
   double GetAuxiliaryHessianElement3(int lambda, 
                                      int sigma, 
                                      int indexAtomA,
                                      int indexAtomC,
                                      MolDS_base::CartesianType axisA1,
                                      MolDS_base::CartesianType axisA2,
                                      double const* const* orbitalElectronPopulation,
                                      double const* const* const* const* const* const* diatomicTwoElecTwoCore2ndDerivs) const;
   double GetAuxiliaryHessianElement4(int lambda, 
                                      int sigma, 
                                      int indexAtomA,
                                      int indexAtomB,
                                      int indexAtomC,
                                      MolDS_base::CartesianType axisA,
                                      MolDS_base::CartesianType axisB,
                                      double const* const* const* const* orbitalElectronPopulation1stDerivs,
                                      double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const;
   double GetAuxiliaryHessianElement5(int mu, 
                                      int lambda, 
                                      int indexAtomA,
                                      int indexAtomC,
                                      MolDS_base::CartesianType axisA1,
                                      MolDS_base::CartesianType axisA2,
                                      double const* const* orbitalElectronPopulation,
                                      double const* const* const* const* diatomicOverlap2ndDerivs) const;
   double GetAuxiliaryHessianElement6(int mu, 
                                      int lambda, 
                                      int indexAtomA,
                                      int indexAtomB,
                                      int indexAtomC,
                                      MolDS_base::CartesianType axisA,
                                      MolDS_base::CartesianType axisB,
                                      double const* const* const* const* orbitalElectronPopulation1stDerivs,
                                      double const* const* const* diatomicOverlap1stDerivs) const;
   double GetAuxiliaryHessianElement7(int mu, 
                                      int nu, 
                                      int lambda, 
                                      int sigma, 
                                      int indexAtomA,
                                      int indexAtomC,
                                      MolDS_base::CartesianType axisA1,
                                      MolDS_base::CartesianType axisA2,
                                      double const* const* orbitalElectronPopulation,
                                      double const* const* const* const* const* const* diatomicTwoElecTwoCore2ndDerivs) const;
   double GetAuxiliaryHessianElement8(int mu, 
                                      int nu, 
                                      int lambda, 
                                      int sigma, 
                                      int indexAtomA,
                                      int indexAtomB,
                                      int indexAtomC,
                                      MolDS_base::CartesianType axisA,
                                      MolDS_base::CartesianType axisB,
                                      double const* const* orbitalElectronPopulation,
                                      double const* const* const* const* orbitalElectronPopulation1stDerivs,
                                      double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const;
   void CalcOrbitalElectronPopulation1stDerivatives(double**** orbitalElectronPopulation1stDerivatives) const;
   void SolveCPHF(double** solutionsCPHF,
                  const std::vector<MoIndexPair>& nonRedundantQIndeces,
                  const std::vector<MoIndexPair>& redundantQIndeces) const;
   void CalcStaticFirstOrderFocks(double** staticFirstOrderFocks,
                                  const std::vector<MoIndexPair>& nonRedundantQIndeces,
                                  const std::vector<MoIndexPair>& redundantQIndeces) const;
   void CalcStaticFirstOrderFock(double* staticFirstOrderFock,
                                 const std::vector<MoIndexPair>& nonRedundantQIndeces,
                                 const std::vector<MoIndexPair>& redundantQIndeces,
                                 int indexAtomA,
                                 MolDS_base::CartesianType axisA) const;
   void MallocTempMatricesStaticFirstOrderFock(double****** diatomicTwoElecTwoCore1stDeriv,
                                               double****   diatomicOverlap1stDeriv) const;
   void FreeTempMatricesStaticFirstOrderFock(double****** diatomicTwoElecTwoCore1stDeriv,
                                             double****   diatomicOverlap1stDeriv) const;
   void CalcMatrixCPHF(double** matrixCPHF, 
                       const std::vector<MoIndexPair>& nonRedundantQIndeces,
                       const std::vector<MoIndexPair>& redundantQIndeces) const;
   void MallocTempMatricesSolveCPHF(double*** matrixCPHF,
                                    int dimensionCPHF) const;
   void FreeTempMatricesSolveCPHF(double*** matrixCPHF,
                                  int dimensionCPHF) const;
   void CalcHeatsFormation(double* heatsFormation, 
                           const MolDS_base::Molecule& molecule) const;
   double GetElectronCoreAttraction(int indexAtomA, 
                                    int indexAtomB, 
                                    int mu, 
                                    int nu, 
                                    double const* const* const* const* const* const* twoElecTwoCore) const;
   double GetElectronCoreAttraction1stDerivative(int indexAtomA, 
                                                 int indexAtomB, 
                                                 int mu, 
                                                 int nu, 
                                                 double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivatives,
                                                 MolDS_base::CartesianType axisA) const;
   void CalcDiatomicTwoElecTwoCore(double**** matrix, int indexAtomA, int indexAtomB) const;
   void CalcDiatomicTwoElecTwoCore1stDerivatives(double***** matrix, 
                                                 int indexAtomA, 
                                                 int indexAtomB) const;
   void CalcDiatomicTwoElecTwoCore2ndDerivatives(double****** matrix, 
                                                 int indexAtomA, 
                                                 int indexAtomB) const;
   void MallocDiatomicTwoElecTwoCore1stDeriTemps(double***   rotatingMatrix,
                                                 double****  rotMat1stDerivatives,
                                                 double***** diatomicTwoElecTwoCore) const;
   void MallocDiatomicTwoElecTwoCore2ndDeriTemps(double***    rotatingMatrix,
                                                 double****   rotMat1stDerivatives,
                                                 double*****  rotMat2ndDerivatives,
                                                 double*****  diatomicTwoElecTwoCore,
                                                 double****** diatomicTwoElecTwoCore1stDerivatives) const;
   void FreeDiatomicTwoElecTwoCore1stDeriTemps(double***   rotatingMatrix,
                                               double****  rotMat1stDerivatives,
                                               double***** diatomicTwoElecTwoCore) const;
   void FreeDiatomicTwoElecTwoCore2ndDeriTemps(double***    rotatingMatrix,
                                               double****   rotMat1stDerivatives,
                                               double*****  rotMat2ndDerivatives,
                                               double*****  diatomicTwoElecTwoCore,
                                               double****** diatomicTwoElecTwoCore1stDerivatives) const;
   void RotateDiatomicTwoElecTwoCoreToSpaceFrame(double**** matrix, 
                                                   double const* const* rotatingMatrix) const;
   void RotateDiatomicTwoElecTwoCore1stDerivativesToSpaceFrame(
        double***** matrix, 
        double const* const* const* const* diatomicTwoElecTwoCore,
        double const* const* rotatingMatrix,
        double const* const* const* rotMat1stDerivatives) const;
   void RotateDiatomicTwoElecTwoCore2ndDerivativesToSpaceFrame(
        double****** matrix, 
        double const* const* const* const*        diatomicTwoElecTwoCore,
        double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivatives,
        double const* const* rotatingMatrix,
        double const* const* const* rotMat1stDerivatives,
        double const* const* const* const* rotMat2ndDerivatives) const;
   double GetNddoRepulsionIntegral(const MolDS_base_atoms::Atom& atomA, 
                                   MolDS_base::OrbitalType mu, 
                                   MolDS_base::OrbitalType nu,
                                   const MolDS_base_atoms::Atom& atomB, 
                                   MolDS_base::OrbitalType lambda, 
                                   MolDS_base::OrbitalType sigma) const;
   double GetNddoRepulsionIntegral1stDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                MolDS_base::OrbitalType mu, 
                                                MolDS_base::OrbitalType nu,
                                                const MolDS_base_atoms::Atom& atomB, 
                                                MolDS_base::OrbitalType lambda, 
                                                MolDS_base::OrbitalType sigma,
                                                MolDS_base::CartesianType axisA) const;
   double GetNddoRepulsionIntegral2ndDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                MolDS_base::OrbitalType mu, 
                                                MolDS_base::OrbitalType nu,
                                                const MolDS_base_atoms::Atom& atomB, 
                                                MolDS_base::OrbitalType lambda, 
                                                MolDS_base::OrbitalType sigma,
                                                MolDS_base::CartesianType axisA1,
                                                MolDS_base::CartesianType axisA2) const;
   double GetSemiEmpiricalMultipoleInteraction(const MolDS_base_atoms::Atom& atomA,
                                               const MolDS_base_atoms::Atom& atomB,
                                               MolDS_base::MultipoleType multipoleA,
                                               MolDS_base::MultipoleType multipoleB,
                                               double Rab) const;
   double GetSemiEmpiricalMultipoleInteraction1stDerivative(const MolDS_base_atoms::Atom& atomA,
                                                            const MolDS_base_atoms::Atom& atomB,
                                                            MolDS_base::MultipoleType multipoleA,
                                                            MolDS_base::MultipoleType multipoleB,
                                                            double Rab) const;
   double GetSemiEmpiricalMultipoleInteraction2ndDerivative(const MolDS_base_atoms::Atom& atomA,
                                                            const MolDS_base_atoms::Atom& atomB,
                                                            MolDS_base::MultipoleType multipoleA,
                                                            MolDS_base::MultipoleType multipoleB,
                                                            double Rab) const;
   void MallocTempMatricesCalcForce(double****   diatomicOverlap1stDerivs, 
                                    double****** diatomiTwoElecTwoCore1stDerivs) const;
   void FreeTempMatricesCalcForce(double****   diatomicOverlap1stDerivs, 
                                  double****** diatomiTwoElecTwoCore1stDerivs) const;
   void CalcForceSCFElecCoreAttractionPart(double* force, 
                                           int indexAtomA,
                                           int indexAtomB,
                                           double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const;
   void CalcForceSCFOverlapPart(double* force, 
                                int indexAtomA,
                                int indexAtomB,
                                double const* const* const* diatomicOverlap1stDerivs) const;
   void CalcForceSCFTwoElecPart(double* force, 
                                int indexAtomA,
                                int indexAtomB,
                                double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const;
   void CalcForceExcitedStaticPart(double* force, 
                                   int elecStateIndex,
                                   int indexAtomA,
                                   int indexAtomB,
                                   double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const;
   void CalcForceExcitedElecCoreAttractionPart(double* force, 
                                               int elecStateIndex,
                                               int indexAtomA,
                                               int indexAtomB,
                                               double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const;
   void CalcForceExcitedOverlapPart(double* force, 
                                    int elecStateIndex,
                                    int indexAtomA,
                                    int indexAtomB,
                                    double const* const* const* diatomicOverlap1stDerivs) const;
   void CalcForceExcitedTwoElecPart(double* force, 
                                    int elecStateIndex,
                                    int indexAtomA,
                                    int indexAtomB,
                                    double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const;

};

}
#endif



