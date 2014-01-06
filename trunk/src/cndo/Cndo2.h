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
#ifndef INCLUDED_CNDO
#define INCLUDED_CNDO
namespace MolDS_cndo{

/***
 *  References for Cndo2 are [PB_1970], [PSS_1965], and [PS_1965].
 */
class Cndo2 : public MolDS_base::ElectronicStructure{
public:
   Cndo2();
   virtual ~Cndo2();
   virtual void SetMolecule(MolDS_base::Molecule* molecule);
   void DoSCF(bool requiresGuess=true);
   virtual void OutputSCFResults() const;
   double const* const* GetFockMatrix() const{return this->fockMatrix;};
   double const*        GetEnergiesMO() const{return this->energiesMO;};
   virtual void DoCIS();
   virtual void OutputCISResults() const;
   double const* const* GetMatrixCIS() const{return this->matrixCIS;};
   double const*        GetExcitedEnergies() const{return this->excitedEnergies;};
   double const* const*        GetForce(int elecState);
   double const* const* const* GetForce(const std::vector<int>& elecStates);
   double GetElectronicEnergy(int elecState) const;
   double GetCoreRepulsionEnergy() const;
   double GetVdWCorrectionEnergy() const;
   void CalcOverlapAOsWithAnotherConfiguration(double** overlapAOs,
                                               const MolDS_base::Molecule& lhsMolecule) const;
   void CalcOverlapMOsWithAnotherElectronicStructure(double** overlapMOs,
                                                     double const* const* overlapAOs,
                                                     const MolDS_base::ElectronicStructure& lhsElectronicStructure) const;
   virtual void CalcOverlapSingletSDsWithAnotherElectronicStructure(double** overlapSingletSDs, 
                                                                    double const* const* overlapMOs) const;
   virtual void CalcOverlapESsWithAnotherElectronicStructure(double** overlapESs, 
                                                             double const* const* overlapSingletSDs,
                                                             const MolDS_base::ElectronicStructure& lhsElectronicStructure) const;
   MolDS_base::TheoryType GetTheoryType() const;
protected:
   std::string errorMessageAtomA;
   std::string errorMessageAtomB;
   std::string errorMessageAtomType;
   std::string errorMessageOrbitalType;
   std::string errorMessageCartesianType;
   std::string errorMessageSCFNotConverged;
   std::string errorMessageMoleculeNotSet;
   std::string errorMessageOddTotalValenceElectrions;
   std::string errorMessageNotEnebleAtomType;
   std::string errorMessageCoulombInt;
   std::string errorMessageExchangeInt;
   std::string errorMessageMolecularIntegralElement;
   std::string errorMessageGetDiatomCoreRepulsion2ndDerivativeNotImplemented;
   std::string errorMessageGetGaussianCartesianMatrixBadOrbital;
   std::string errorMessageGetGaussianOverlapAOsBadOrbital;
   std::string errorMessageGetGaussianOverlapAOs1stDerivativeOrbitalD;
   std::string errorMessageCISNotImplemented;
   std::string errorMessageCalcForceNotImplemented;
   std::string errorMessageGetElectronicEnergyNULLCISEnergy;
   std::string errorMessageGetElectronicEnergyEnergyNotCalculated;
   std::string errorMessageGetElectronicEnergyNumberCISStates;
   std::string errorMessageGetElectronicEnergySetElecState;
   std::string errorMessageCalcElectronicTransitionDipoleMomentBadState;
   std::string errorMessageCalcFrequenciesNormalModesBadTheory;
   std::string errorMessageFromState;
   std::string errorMessageToState;
   std::string errorMessageNonExcitedStates;
   std::string messageSCFMetConvergence;
   std::string messageStartSCF;
   std::string messageDoneSCF;
   std::string messageOmpElapsedTimeSCF;
   std::string messageMullikenAtoms;
   std::string messageMullikenAtomsTitle;
   std::string messageUnpairedAtoms;
   std::string messageUnpairedAtomsTitle;
   std::string messageUnitSec; 
   std::vector<MolDS_base::AtomType> enableAtomTypes;
   MolDS_base::Molecule* molecule;
   MolDS_base::TheoryType theory;
   double       coreRepulsionEnergy;
   double       coreEpcCoulombEnergy;
   double       vdWCorrectionEnergy;
   int          matrixCISdimension;
   double**     fockMatrix;
   double*      energiesMO;
   double**     orbitalElectronPopulation; //P_{\mu\nu} of (2.50) in J. A. Pople book.
   double***    orbitalElectronPopulationCIS; 
   double*      atomicElectronPopulation; //P_{AB} of (3.21) in J. A. Pople book.
   double**     atomicElectronPopulationCIS; 
   double**     atomicUnpairedPopulationCIS; 
   double**     overlapAOs; // overlap integral between AOs
   double****** twoElecsTwoAtomCores;
   double****** twoElecsAtomEpcCores;
   double***    cartesianMatrix; // cartesian matrix represented by AOs
   double***    electronicTransitionDipoleMoments; // Diagnonal terms are electronic dipole moments of each eigenstates (i.e. electronicDipole[0][0][XAxis] is the x-component of the electronic dipole moment of the ground state. electronicDipole[10][10][XAxis] is the x-component of the electronic dipole moment of the 10-th excited state). Off-diagonal terms are transition dipole moments between eigenstates (i.e. electronicDipole[10][0][XAxis] is the x-component of the transition dipole moment from the ground state to 10-th excited state.).
   double*      coreDipoleMoment; // dipole moment of configuration.
   double*      normalForceConstants; // force constants of normal modes
   double**     normalModes; // in mass-weighted coordinates
   double**     matrixCIS;
   double*      excitedEnergies;
   double*      freeExcitonEnergiesCIS;
   double***    matrixForce;
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual void CalcSCFProperties();
   virtual void CalcCISProperties();
   virtual void CalcNormalModes(double** normalModes, double* normalForceConstants, const MolDS_base::Molecule& molecule) const;
   virtual void CalcElectronicTransitionDipoleMoment(double* transitionDipoleMoment,
                                                     int to, int from,
                                                     double const* const* fockMatrix,
                                                     double const* const* matrixCIS,
                                                     double const* const* const* cartesianMatrix,
                                                     const MolDS_base::Molecule& molecule, 
                                                     double const* const* orbitalElectronPopulation,
                                                     double const* const* overlapAOs,
                                                     double const* groundStateDipole) const;
   double GetBondingAdjustParameterK(MolDS_base::ShellType shellA, 
                                     MolDS_base::ShellType shellB) const;
   virtual double GetAtomCoreEpcCoulombEnergy (const MolDS_base_atoms::Atom& atom, 
                                               const MolDS_base_atoms::Atom& epc) const;
   virtual double GetDiatomCoreRepulsionEnergy(const MolDS_base_atoms::Atom& atomA,
                                               const MolDS_base_atoms::Atom& atomB) const;
   virtual double GetDiatomCoreRepulsion1stDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                      const MolDS_base_atoms::Atom& atomB, 
                                                      MolDS_base::CartesianType axisA) const;
   virtual double GetDiatomCoreRepulsion2ndDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                      const MolDS_base_atoms::Atom& atomB, 
                                                      MolDS_base::CartesianType axisA1,
                                                      MolDS_base::CartesianType axisA2) const;
   virtual double GetDiatomVdWCorrectionEnergy(const MolDS_base_atoms::Atom& atomA, 
                                               const MolDS_base_atoms::Atom& atomB) const;
   virtual double GetDiatomVdWCorrection1stDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                      const MolDS_base_atoms::Atom& atomB, 
                                                      MolDS_base::CartesianType axisA) const;
   virtual double GetDiatomVdWCorrection2ndDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                      const MolDS_base_atoms::Atom& atomB, 
                                                      MolDS_base::CartesianType axisA1,
                                                      MolDS_base::CartesianType axisA2) const;
   double GetReducedOverlapAOs                      (int na, int nb, double alpha, double beta) const;
   double GetReducedOverlapAOs                      (int na, int la, int m, int nb, int lb, double alpha, double beta) const;
   double GetReducedOverlapAOs1stDerivativeAlpha    (int na, int la, int m, int nb, int lb, double alpha, double beta) const;
   double GetReducedOverlapAOs1stDerivativeBeta     (int na, int la, int m, int nb, int lb, double alpha, double beta) const;
   double GetReducedOverlapAOs2ndDerivativeAlpha    (int na, int la, int m, int nb, int lb, double alpha, double beta) const;
   double GetReducedOverlapAOs2ndDerivativeBeta     (int na, int la, int m, int nb, int lb, double alpha, double beta) const;
   double GetReducedOverlapAOs2ndDerivativeAlphaBeta(int na, int la, int m, int nb, int lb, double alpha, double beta) const;
   double GetOverlapAOsElement1stDerivativeByGTOExpansion(const MolDS_base_atoms::Atom& atomA, 
                                                          int valenceIndexA, 
                                                          const MolDS_base_atoms::Atom& atomB, 
                                                          int valenceIndexB,
                                                          MolDS_base::STOnGType stonG, 
                                                          MolDS_base::CartesianType axisA) const; // See [DY_1977].
   void CalcRotatingMatrix(double** rotatingMatrix, 
                           const MolDS_base_atoms::Atom& atomA, 
                           const MolDS_base_atoms::Atom& atomB) const;
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
                                        int mu, 
                                        int nu, 
                                        const MolDS_base::Molecule& molecule, 
                                        double const* const* gammaAB, 
                                        double const* const* overlapAOs,
                                        double const* const* orbitalElectronPopulation, 
                                        double const* const* const* const* const* const* twoElecsTwoAtomCores, 
                                        bool isGuess) const;
   void TransposeFockMatrixMatrix(double** transposedFockMatrix) const;
   virtual void CalcDiatomicOverlapAOsInDiatomicFrame(double** diatomicOverlapAOs, 
                                                      const MolDS_base_atoms::Atom& atomA, 
                                                      const MolDS_base_atoms::Atom& atomB) const;
   virtual void CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(double** diatomicOverlapAOsDeri, 
                                                                   const MolDS_base_atoms::Atom& atomA, 
                                                                   const MolDS_base_atoms::Atom& atomB) const;
   virtual void CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(double** diatomicOverlapAOs2ndDeri, 
                                                                   const MolDS_base_atoms::Atom& atomA, 
                                                                   const MolDS_base_atoms::Atom& atomB) const;
   void CalcDiatomicOverlapAOs1stDerivatives(double*** diatomicOverlapAOs1stDerivs, 
                                             double**  tmpDiaOverlapAOsInDiaFrame,         
                                             double**  tmpDiaOverlapAOs1stDerivInDiaFrame,
                                             double**  tmpRotMat, 
                                             double**  tmpRotMat1stDeriv,
                                             double*** tmpRotMat1stDerivs,
                                             double**  tmpRotatedDiatomicOverlap,
                                             double*   tmpRotatedDiatomicOverlapVec,
                                             double**  tmpMatrixBC,
                                             double*   tmpVectorBC,
                                             const MolDS_base_atoms::Atom& atomA, 
                                             const MolDS_base_atoms::Atom& atomB) const;
   void CalcDiatomicOverlapAOs1stDerivatives(double*** diatomicOverlapAOs1stDerivs, 
                                             double**  tmpDiaOverlapAOsInDiaFrame,         
                                             double**  tmpDiaOverlapAOs1stDerivInDiaFrame,
                                             double**  tmpRotMat, 
                                             double**  tmpRotMat1stDeriv,
                                             double*** tmpRotMat1stDerivs,
                                             double**  tmpRotatedDiatomicOverlap,
                                             double*   tmpRotatedDiatomicOverlapVec,
                                             double**  tmpMatrixBC,
                                             double*   tmpVectorBC,
                                             int indexAtomA, 
                                             int indexAtomB) const;
   void CalcDiatomicOverlapAOs2ndDerivatives(double**** overlapAOs2ndDeri, 
                                             double**   tmpDiaOverlapAOsInDiaFrame,
                                             double**   tmpDiaOverlapAOs1stDerivInDiaFrame,
                                             double**   tmpDiaOverlapAOs2ndDerivInDiaFrame,
                                             double***  tmpDiaOverlapAOs1stDerivs,
                                             double**** tmpDiaOverlapAOs2ndDerivs,
                                             double**   tmpRotMat,
                                             double***  tmpRotMat1stDerivs,
                                             double**** tmpRotMat2ndDerivs,
                                             const MolDS_base_atoms::Atom& atomA, 
                                             const MolDS_base_atoms::Atom& atomB) const;
   void CalcDiatomicOverlapAOs2ndDerivatives(double**** overlapAOs2ndDeri, 
                                             double**   tmpDiaOverlapAOsInDiaFrame,
                                             double**   tmpDiaOverlapAOs1stDerivInDiaFrame,
                                             double**   tmpDiaOverlapAOs2ndDerivInDiaFrame,
                                             double***  tmpDiaOverlapAOs1stDerivs,
                                             double**** tmpDiaOverlapAOs2ndDerivs,
                                             double**   tmpRotMat,
                                             double***  tmpRotMat1stDerivs,
                                             double**** tmpRotMat2ndDerivs,
                                             int indexAtomA, 
                                             int indexAtomB) const;
   double Get2ndDerivativeElementFromDistanceDerivatives(double firstDistanceDeri,
                                                         double secondDistanceDeri,
                                                         MolDS_base::CartesianType axisA1,
                                                         MolDS_base::CartesianType axisA2,
                                                         double* cartesian,
                                                         double rAB) const;
   virtual double GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                              const MolDS_base::Molecule& molecule, 
                                              double const* const* fockMatrix, 
                                              double const* const* gammaAB) const;
   virtual void CalcTwoElecsTwoCores(double****** twoElecsTwoAtomCores, 
                                     double****** twoElecsAtomEpcCores,
                                     const MolDS_base::Molecule& molecule) const;
   virtual void CalcForce(const std::vector<int>& elecStates);
   void CalcRotatingMatrix1stDerivatives(double*** rotMat1stDerivatives, 
                                         const MolDS_base_atoms::Atom& atomA,
                                         const MolDS_base_atoms::Atom& atomB) const;
   void CalcRotatingMatrix2ndDerivatives(double**** rotMat2ndDerivatives, 
                                         const MolDS_base_atoms::Atom& atomA,
                                         const MolDS_base_atoms::Atom& atomB) const;
   struct MoEnergyGap{
      double energyGap;
      int occIndex;
      int virIndex;
      int slaterIndex;
   };
   struct LessMoEnergyGap { 
      bool operator()(const MoEnergyGap& left, const MoEnergyGap& right) 
      const { return left.energyGap < right.energyGap; } 
   };
   struct CISEigenVectorCoefficient{
      double coefficient;
      int occIndex;
      int virIndex;
      int slaterIndex;
   };
   struct MoreCISEigenVectorCoefficient { 
      bool operator()(const CISEigenVectorCoefficient& left, const CISEigenVectorCoefficient& right) 
      const { return fabs(left.coefficient) > fabs(right.coefficient); } 
   };
private:
   std::string errorMessageCalDiaOverlapAOsDiaFrameNullMatrix;
   std::string errorMessageCalcRotatingMatrixNullRotMatrix;
   std::string errorMessageRotDiaOverlapAOsToSpaceFrameNullDiaMatrix;
   std::string errorMessageRotDiaOverlapAOsToSpaceFrameNullRotMatrix;
   std::string errorMessageRotDiaOverlapAOsToSpaceFrameNullTmpDiaMatrix;
   std::string errorMessageRotDiaOverlapAOsToSpaceFrameNullTmpOldDiaMatrix;
   std::string errorMessageRotDiaOverlapAOsToSpaceFrameNullTmpMatrixBC;
   std::string errorMessageRotDiaOverlapAOsToSpaceFrameNullTmpVectorBC;
   std::string errorMessageSetOverlapAOsElementNullDiaMatrix;
   std::string errorMessageCalcOverlapAOsDifferentConfigurationsDiffAOs;
   std::string errorMessageCalcOverlapAOsDifferentConfigurationsDiffAtoms;
   std::string errorMessageCalcOverlapAOsDifferentConfigurationsOverlapAOsNULL;
   std::string errorMessageLhs;
   std::string errorMessageRhs;
   std::string messageIterSCFTitle;
   std::string messageIterSCF;
   std::string messageDiisApplied;
   std::string messageDampingApplied;
   std::string messageEnergyMO;
   std::string messageEnergyMOTitle;
   std::string messageElecEnergy;
   std::string messageNoteElecEnergy;
   std::string messageNoteElecEnergyVdW;
   std::string messageNoteElecEnergyEpcVdW;
   std::string messageNoteElecEnergyEpc;
   std::string messageElecEnergyTitle;
   std::string messageOcc;
   std::string messageUnOcc;
   std::string messageCoreRepulsionTitle;
   std::string messageCoreRepulsion;
   std::string messageCoreEpcCoulombTitle;
   std::string messageCoreEpcCoulomb;
   std::string messageVdWCorrectionTitle;
   std::string messageVdWCorrection;
   std::string messageElectronicDipoleMomentTitle;
   std::string messageElectronicDipoleMoment;
   std::string messageCoreDipoleMomentTitle;
   std::string messageCoreDipoleMoment;
   std::string messageTotalDipoleMomentTitle;
   std::string messageTotalDipoleMoment;
   std::string messageMullikenAtomsSCF;
   std::string messageNormalModesTitle;
   std::string messageNormalModesUnitsMassWeighted;
   std::string messageNormalModesUnitsNonMassWeighted;
   std::string messageNormalModesMassWeighted;
   std::string messageNormalModesNonMassWeighted;
   std::string messageNormalModesImaginaryFrequencies;
   double elecSCFEnergy;
   double bondingAdjustParameterK[2]; //see (3.79) in J. A. Pople book
   double** gammaAB;
   class ReducedOverlapAOsParameters : private MolDS_base::Uncopyable{
   public:
      // use Y[na][nb][la][lb][m][i][j] 
      // as Y_{ij\lammda} in (B.20) in Pople book for given na, nb, la, lb, m, i, and j.
      static const double Y[MolDS_base::ShellType_end+1]
                           [MolDS_base::ShellType_end+1]
                           [MolDS_base::ShellType_end]
                           [MolDS_base::ShellType_end]
                           [MolDS_base::ShellType_end]
                           [2*MolDS_base::ShellType_end+1]
                           [2*MolDS_base::ShellType_end+1];
      // use Z[na][nb][k] as Z_{k} in (B.30) in Pople book for given na, nb, and k. 
      static const double Z[2*MolDS_base::ShellType_end]
                           [2*MolDS_base::ShellType_end]
                           [4*MolDS_base::ShellType_end-1];
   private:
      ReducedOverlapAOsParameters();
      ~ReducedOverlapAOsParameters();
   };
   void OutputMOEnergies() const;
   void OutputSCFEnergies() const;
   void OutputSCFDipole() const;
   void OutputSCFMulliken() const;
   void OutputNormalModes(double const* const* normalModes, 
                          double const* normalForceConstants, 
                          const MolDS_base::Molecule& molecule) const;
   void CalcCoreRepulsionEnergy();
   void CalcVdWCorrectionEnergy();
   double GetVdwDampingValue(double vdWDistance, double distance) const;
   double GetVdwDampingValue1stDerivative(double vdWDistance, double distance) const;
   double GetVdwDampingValue2ndDerivative(double vdWDistance, double distance) const;
   void CalcElectronicDipoleMomentGroundState(double*** electronicTransitionDipoleMoments,
                                              double const* const* const* cartesianMatrix,
                                              const MolDS_base::Molecule& molecule, 
                                              double const* const* orbitalElectronPopulation,
                                              double const* const* overlapAOs) const;
   bool SatisfyConvergenceCriterion(double const* const* oldOrbitalElectronPopulation, 
                                    double const* const* orbitalElectronPopulation,
                                    int     numberAOs, 
                                    double* rmsDensity, 
                                    int     times,
                                    double  diisErro,
                                    bool    hasAppliedDIIS,
                                    bool    hasAppliedDamping) const;
   void UpdateOldOrbitalElectronPopulation(double** oldOrbitalElectronPopulation, 
                                           double const* const* orbitalElectronPopulation,
                                           int numberAOs) const;
   void CalcOrbitalElectronPopulation(double** orbitalElectronPopulation, 
                                      const MolDS_base::Molecule& molecule, 
                                      double const* const* fockMatrix) const;
   void CalcAtomicElectronPopulation(double* atomicElectronPopulation,
                                     double const* const* orbitalElectronPopulation, 
                                     const MolDS_base::Molecule& molecule) const;
   void CalcCoreDipoleMoment(double* coreDipoleMoment,
                       const MolDS_base::Molecule& molecule) const;
   void CalcCartesianMatrixByGTOExpansion(double*** cartesianMatrix,
                                          const MolDS_base::Molecule& molecule, 
                                          MolDS_base::STOnGType stonG) const; 
   void CalcCartesianMatrixElementsByGTOExpansion(double& xComponent,
                                                  double& yComponent,
                                                  double& zComponent,
                                                  const MolDS_base_atoms::Atom& atomA, 
                                                  int valenceIndexA, 
                                                  const MolDS_base_atoms::Atom& atomB, 
                                                  int valenceIndexB,
                                                  MolDS_base::STOnGType stonG) const;
   double GetGaussianCartesianMatrix(MolDS_base::AtomType atomTypeA, 
                                     MolDS_base::OrbitalType valenceOrbitalA, 
                                     double gaussianExponentA, 
                                     double const* xyzA,
                                     MolDS_base::AtomType atomTypeB, 
                                     MolDS_base::OrbitalType valenceOrbitalB, 
                                     double gaussianExponentB, 
                                     double const* xyzB,
                                     double rAB,
                                     MolDS_base::CartesianType axis) const;
   double GetGaussianCartesianMatrix(MolDS_base::AtomType atomTypeA, 
                                     MolDS_base::OrbitalType valenceOrbitalA, 
                                     double gaussianExponentA, 
                                     double const* xyzA,
                                     MolDS_base::AtomType atomTypeB, 
                                     MolDS_base::OrbitalType valenceOrbitalB, 
                                     double gaussianExponentB, 
                                     double const* xyzB,
                                     double rAB,
                                     double ovelapSASB,
                                     MolDS_base::CartesianType axis) const;
   void CalcOverlapAOs(double** overlapAOs, const MolDS_base::Molecule& molecule) const;
   void CalcOverlapAOsByGTOExpansion(double** overlapAOs, 
                                     const MolDS_base::Molecule& molecule, 
                                     MolDS_base::STOnGType stonG) const; //See [DY_1977]
   double GetOverlapAOsElementByGTOExpansion(const MolDS_base_atoms::Atom& atomA, 
                                             int valenceIndexA, 
                                             const MolDS_base_atoms::Atom& atomB, 
                                             int valenceIndexB,
                                             MolDS_base::STOnGType stonG) const; // see [DY_1977]
   double GetGaussianOverlapAOsSASB(double gaussianExponentA, 
                                    double gaussianExponentB, 
                                    double rAB) const; // see [DY_1977]
   double GetGaussianOverlapAOs(MolDS_base::AtomType atomTypeA, 
                                MolDS_base::OrbitalType valenceOrbitalA, 
                                double gaussianExponentA, 
                                MolDS_base::AtomType atomTypeB, 
                                MolDS_base::OrbitalType valenceOrbitalB, 
                                double gaussianExponentB, 
                                double dx, 
                                double dy, 
                                double dz, 
                                double rAB) const; // see [DY_1977]
   double GetGaussianOverlapAOs(MolDS_base::AtomType atomTypeA, 
                                MolDS_base::OrbitalType valenceOrbitalA, 
                                double gaussianExponentA, 
                                MolDS_base::AtomType atomTypeB, 
                                MolDS_base::OrbitalType valenceOrbitalB, 
                                double gaussianExponentB, 
                                double dx, 
                                double dy, 
                                double dz, 
                                double rAB,
                                double ovelapSASB) const; // see [DY_1977]
   double GetGaussianOverlapAOs1stDerivative(MolDS_base::AtomType atomTypeA, 
                                             MolDS_base::OrbitalType valenceOrbitalA, 
                                             double gaussianExponentA, 
                                             MolDS_base::AtomType atomTypeB, 
                                             MolDS_base::OrbitalType valenceOrbitalB, 
                                             double gaussianExponentB, 
                                             double dx, 
                                             double dy, 
                                             double dz, 
                                             double rAB, 
                                             MolDS_base::CartesianType axisA) const;// see [DY_1977]
   void CalcFockMatrix(double** fockMatrix, 
                       const MolDS_base::Molecule& molecule, 
                       double const* const* overlapAOs, 
                       double const* const* gammaAB,
                       double const* const* orbitalElectronPopulation, 
                       double const* atomicElectronPopulation,
                       double const* const* const* const* const* const* twoElecsTwoAtomCores,
                       bool isGuess) const;
   void RotateDiatmicOverlapAOsToSpaceFrame(double**             diatomicOverlapAOs, 
                                            double const* const* rotatingMatrix,
                                            double*              tmpDiatomicOverlapAOs,
                                            double**             tmpOldDiatomicOverlapAOs,
                                            double**             tmpMatrixBC,
                                            double*              tmpVectorBC) const;
   void SetOverlapAOsElement(double** overlapAOs, 
                             double const* const* diatomicOverlapAOs, 
                             const MolDS_base_atoms::Atom& atomA, 
                             const MolDS_base_atoms::Atom& atomB,
                             bool isSymmetricOverlapAOs = true) const;
   double GetAuxiliaryA(int k, double rho) const;
   double GetAuxiliaryB(int k, double rho) const;
   double GetAuxiliaryD(int la, int lb, int m) const;
   double GetAuxiliaryA1stDerivative(int k, double rho) const;
   double GetAuxiliaryA2ndDerivative(int k, double rho) const;
   double GetAuxiliaryB1stDerivative(int k, double rho) const;
   double GetAuxiliaryB2ndDerivative(int k, double rho) const;
   void DoDamp(double rmsDensity, 
               bool&  hasAppliedDamping,
               double** orbitalElectronPopulation, 
               double const* const* oldOrbitalElectronPopulation, 
               const MolDS_base::Molecule& molecule) const;
   void DoDIIS(double** orbitalElectronPopulation,
               double const* const* oldOrbitalElectronPopulation,
               double*** diisStoredDensityMatrix,
               double*** diisStoredErrorVect,
               double**  diisErrorProducts,
               double**  tmpDiisErrorProducts,
               double*   diisErrorCoefficients,
               double&   diisError,
               bool&     hasAppliedDIIS,
               int       diisNumErrorVect,
               const MolDS_base::Molecule& molecule, 
               int step) const;
   void CheckEnableAtomType(const MolDS_base::Molecule& molecule) const;
   void CheckNumberValenceElectrons(const MolDS_base::Molecule& molecule) const;
   void FreeDiatomicOverlapAOsAndRotatingMatrix(double*** diatomicOverlapAOs, 
                                                double*** rotatingMatrix) const;
   void CalcElecSCFEnergy(double* elecSCFEnergy, 
                          const MolDS_base::Molecule& molecule, 
                          double const* energiesMO, 
                          double const* const* fockMatrix, 
                          double const* const* gammaAB, 
                          double coreRepulsionEnergy,
                          double coreEpcCoulombEnergy,
                          double vdWCorrectionEnergy) const;
   void FreeElecEnergyMatrices(double*** fMatrix, 
                               double*** hMatrix, 
                               double*** dammyOrbitalElectronPopulation, 
                               double**  dammyAtomicElectronPopulation ) const;
   void FreeSCFTemporaryMatrices(double***  oldOrbitalElectronPopulation,
                                 double**** diisStoredDensityMatrix,
                                 double**** diisStoredErrorVect,
                                 double***  diisErrorProducts,
                                 double***  tmpDiisErrorProducts,
                                 double**   diisErrorCoefficients) const;
   void MallocSCFTemporaryMatrices(double***  oldOrbitalElectronPopulation,
                                   double**** diisStoredDensityMatrix,
                                   double**** diisStoredErrorVect,
                                   double***  diisErrorProducts,
                                   double***  tmpDiisErrorProducts,
                                   double**   diisErrorCoefficients);
};

}
#endif



