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
   std::string errorMessageGetSemiEmpiricalMultipoleInteractionFirstDeriBadMultipoles;
   std::string errorMessageGetSemiEmpiricalMultipoleInteractionSecondDeriBadMultipoles;
   std::string errorMessageGetNddoRepulsionIntegral;
   std::string errorMessageGetNddoRepulsionIntegralFirstDerivative;
   std::string errorMessageGetNddoRepulsionIntegralSecondDerivative;
   std::string errorMessageCalcTwoElecTwoCoreDiatomicNullMatrix;
   std::string errorMessageCalcTwoElecTwoCoreNullMatrix;
   std::string errorMessageCalcTwoElecTwoCoreDiatomicSameAtoms;
   std::string errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesSameAtoms;
   std::string errorMessageCalcTwoElecTwoCoreDiatomicSecondDerivativesSameAtoms;
   std::string errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesNullMatrix;
   std::string errorMessageCalcTwoElecTwoCoreDiatomicSecondDerivativesNullMatrix;
   std::string errorMessageCalcZMatrixForceEtaNull;
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual void CalcSCFProperties();
   virtual void CalcForce(const std::vector<int>& elecStates);
   virtual double GetDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB) const;
   virtual double GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
                                                        int atomBIndex, 
                                                        MolDS_base::CartesianType axisA) const;
   virtual double GetFockDiagElement(const MolDS_base_atoms::Atom& atomA, 
                                     int atomAIndex, 
                                     int mu, 
                                     const MolDS_base::Molecule& molecule, 
                                     double const* const* gammaAB,
                                     double const* const* orbitalElectronPopulation, 
                                     double const* atomicElectronPopulation,
                                     double const* const* const* const* const* const* twoElecTwoCore,
                                     bool isGuess) const;
   virtual double GetFockOffDiagElement(const MolDS_base_atoms::Atom& atomA, 
                                        const MolDS_base_atoms::Atom& atomB, 
                                        int atomAIndex, 
                                        int atomBIndex, 
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
   virtual void CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(double** diatomicOverlapDeri, 
                                                                  const MolDS_base_atoms::Atom& atomA, 
                                                                  const MolDS_base_atoms::Atom& atomB) const;
   virtual void CalcDiatomicOverlapSecondDerivativeInDiatomicFrame(double** diatomicOverlapSecondDeri, 
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
   double*** zMatrixForce;
   double*** etaMatrixForce;
   int zMatrixForceElecStatesNum;
   int etaMatrixForceElecStatesNum;
   double heatsFormation;
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
   double GetCISCoefficientMOEnergy(int k, 
                                    int l, 
                                    int r, 
                                    int numberActiveVir) const;
   double GetCISCoefficientTwoElecIntegral(int k, 
                                           int l, 
                                           int p, 
                                           int q, 
                                           int r, 
                                           int s, 
                                           int numberActiveVir) const;
   void CalcActiveSetVariablesQ(std::vector<MoIndexPair>* nonRedundantQIndeces, 
                                std::vector<MoIndexPair>* redundantQIndeces,
                                int numberActiveOcc,
                                int numberActiveVir) const;
   void SolveCPHF(double* solution, 
                  double const* const* transposedFockMatrix,
                  int atomAIndex, 
                  MolDS_base::CartesianType axis) const;
   void CalcStaticFirstOrderFock(double* staticFirstOrderFock,
                                 const std::vector<MoIndexPair>& nonRedundantQIndeces,
                                 const std::vector<MoIndexPair>& redundantQIndeces,
                                 double const* const* transposedFockMatrix,
                                 int atomAIndex,
                                 MolDS_base::CartesianType axis) const;
   void CalcMatrixCPHF(double** matrixCPHF, 
                       const std::vector<MoIndexPair>& nonRedundantQIndeces,
                       const std::vector<MoIndexPair>& redundantQIndeces) const;
   void MallocTempMatricesSolveCPHF(double** staticFirstOrderFock,
                                    double*** matrixCPHF,
                                    int dimensionCPHF) const;
   void FreeTempMatricesSolveCPHF(double** staticFirstOrderFock,
                                  double*** matrixCPHF,
                                  int dimensionCPHF) const;
   void CalcHeatsFormation(double* heatsFormation, 
                           const MolDS_base::Molecule& molecule) const;
   double GetElectronCoreAttraction(int atomAIndex, 
                                    int atomBIndex, 
                                    int mu, 
                                    int nu, 
                                    double const* const* const* const* const* const* twoElecTwoCore) const;
   double GetElectronCoreAttractionFirstDerivative(int atomAIndex, 
                                                   int atomBIndex, 
                                                   int mu, 
                                                   int nu, 
                                                   double const* const* const* const* const* twoElecTwoCoreFirstDerivative,
                                                   MolDS_base::CartesianType axisA) const;
   void CalcTwoElecTwoCoreDiatomic(double**** matrix, int atomAIndex, int atomBIndex) const;
   void CalcTwoElecTwoCoreDiatomicFirstDerivatives(double***** matrix, 
                                                   int atomAIndex, 
                                                   int atomBIndex) const;
   void CalcTwoElecTwoCoreDiatomicSecondDerivatives(double****** matrix, 
                                                    int atomAIndex, 
                                                    int atomBIndex) const;
   void MallocTwoElecTwoCoreDiatomicFirstDeriTemps(double*** rotatingMatrix,
                                                   double**** rotMatFirstDerivatives,
                                                   double***** twoElecTwoCoreDiatomic) const;
   void MallocTwoElecTwoCoreDiatomicSecondDeriTemps(double*** rotatingMatrix,
                                                    double**** rotMatFirstDerivatives,
                                                    double***** rotMatSecondDerivatives,
                                                    double***** twoElecTwoCoreDiatomic,
                                                    double****** twoElecTwoCoreDiatomicFirstDerivatives) const;
   void FreeTwoElecTwoCoreDiatomicFirstDeriTemps(double*** rotatingMatrix,
                                                 double**** rotMatFirstDerivatives,
                                                 double***** twoElecTwoCoreDiatomic) const;
   void FreeTwoElecTwoCoreDiatomicSecondDeriTemps(double*** rotatingMatrix,
                                                  double**** rotMatFirstDerivatives,
                                                  double***** rotMatSecondDerivatives,
                                                  double***** twoElecTwoCoreDiatomic,
                                                  double****** twoElecTwoCoreDiatomicFirstDerivatives) const;
   void RotateTwoElecTwoCoreDiatomicToSpaceFramegc(double**** matrix, 
                                                   double const* const* rotatingMatrix) const;
   void RotateTwoElecTwoCoreDiatomicFirstDerivativesToSpaceFramegc(
        double***** matrix, 
        double const* const* const* const* twoElecTwoCoreDiatomic,
        double const* const* rotatingMatrix,
        double const* const* const* rotMatFirstDerivatives) const;
   void RotateTwoElecTwoCoreDiatomicSecondDerivativesToSpaceFramegc(
        double****** matrix, 
        double const* const* const* const* twoElecTwoCoreDiatomic,
        double const* const* const* const* const* twoElecTwoCoreDiatomicFirstDerivatives,
        double const* const* rotatingMatrix,
        double const* const* const* rotMatFirstDerivatives,
        double const* const* const* const* rotMatSecondDerivatives) const;
   double GetNddoRepulsionIntegral(const MolDS_base_atoms::Atom& atomA, 
                                   MolDS_base::OrbitalType mu, 
                                   MolDS_base::OrbitalType nu,
                                   const MolDS_base_atoms::Atom& atomB, 
                                   MolDS_base::OrbitalType lambda, 
                                   MolDS_base::OrbitalType sigma) const;
   double GetNddoRepulsionIntegralFirstDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                  MolDS_base::OrbitalType mu, 
                                                  MolDS_base::OrbitalType nu,
                                                  const MolDS_base_atoms::Atom& atomB, 
                                                  MolDS_base::OrbitalType lambda, 
                                                  MolDS_base::OrbitalType sigma,
                                                  MolDS_base::CartesianType axisA) const;
   double GetNddoRepulsionIntegralSecondDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                   MolDS_base::OrbitalType mu, 
                                                   MolDS_base::OrbitalType nu,
                                                   const MolDS_base_atoms::Atom& atomB, 
                                                   MolDS_base::OrbitalType lambda, 
                                                   MolDS_base::OrbitalType sigma,
                                                   MolDS_base::CartesianType axisA1,
                                                   MolDS_base::CartesianType axisA2) const;
   double GetSemiEmpiricalMultipoleInteraction(MolDS_base::MultipoleType multipoleA,
                                               MolDS_base::MultipoleType multipoleB,
                                               double rhoA,
                                               double rhoB,
                                               double DA,
                                               double DB,
                                               double Rab) const;
   double GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                                               MolDS_base::MultipoleType multipoleA,
                                               MolDS_base::MultipoleType multipoleB,
                                               double rhoA,
                                               double rhoB,
                                               double DA,
                                               double DB,
                                               double Rab) const;
   double GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                                               MolDS_base::MultipoleType multipoleA,
                                               MolDS_base::MultipoleType multipoleB,
                                               double rhoA,
                                               double rhoB,
                                               double DA,
                                               double DB,
                                               double Rab) const;
   void FreeCalcForceTempMatrices(double**** overlapDer, 
                                  double****** twoElecTwoCoreFirstDeriv) const;
   void CalcForceSCFElecCoreAttractionPart(double* force, 
                                          int atomAIndex,
                                          int atomBIndex,
                                          double const* const* const* const* const* twoElecTwoCoreFirstDeriv) const;
   void CalcForceSCFOverlapPart(double* force, 
                                int atomAIndex,
                                int atomBIndex,
                                double const* const* const* overlapDer) const;
   void CalcForceSCFTwoElecPart(double* force, 
                                int atomAIndex,
                                int atomBIndex,
                                double const* const* const* const* const* twoElecTwoCoreFirstDeriv) const;
   void CalcForceExcitedStaticPart(double* force, 
                                   int elecStateIndex,
                                   int atomAIndex,
                                   int atomBIndex,
                                   double const* const* const* const* const* twoElecTwoCoreFirstDeriv) const;
   void CalcForceExcitedElecCoreAttractionPart(double* force, 
                                               int elecStateIndex,
                                               int atomAIndex,
                                               int atomBIndex,
                                               double const* const* const* const* const* twoElecTwoCoreFirstDeriv) const;
   void CalcForceExcitedOverlapPart(double* force, 
                                    int elecStateIndex,
                                    int atomAIndex,
                                    int atomBIndex,
                                    double const* const* const* overlapDer) const;
   void CalcForceExcitedTwoElecPart(double* force, 
                                    int elecStateIndex,
                                    int atomAIndex,
                                    int atomBIndex,
                                    double const* const* const* const* const* twoElecTwoCoreFirstDeriv) const;

};

}
#endif



