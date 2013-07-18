//************************************************************************//
// Copyright (C) 2011-2013 Mikiya Fujii                                   //
// Copyright (C) 2012-2013 Michihiro Okuyama                              //
// Copyright (C) 2013-2013 Katsuhiko Nishimra                             //
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
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<vector>
#include<stdexcept>
#include<omp.h>
#include<boost/format.hpp>
#include"../base/Uncopyable.h"
#include"../mpi/MpiProcess.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../wrappers/Blas.h"
#include"../wrappers/Lapack.h"
#include"../base/Enums.h"
#include"../base/MathUtilities.h"
#include"../base/MallocerFreer.h"
#include"../base/EularAngle.h"
#include"../base/Parameters.h"
#include"../base/atoms/Atom.h"
#include"../base/atoms/Hatom.h"
#include"../base/atoms/Liatom.h"
#include"../base/atoms/Catom.h"
#include"../base/atoms/Natom.h"
#include"../base/atoms/Oatom.h"
#include"../base/atoms/Satom.h"
#include"../base/Molecule.h"
#include"../base/GTOExpansionSTO.h"
#include"../base/RealSphericalHarmonicsIndex.h"
#include"../base/loggers/MOLogger.h"
#include"../base/ElectronicStructure.h"
#include"Cndo2.h"
#include"ReducedOverlapAOsParameters.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_cndo{

/***
 *  References for Cndo2 are [PB_1970], [PSS_1965], and [PS_1965].
 */
Cndo2::Cndo2(){
   //protected variables
   this->molecule = NULL;
   this->theory = CNDO2;
   this->coreRepulsionEnergy = 0.0;
   this->vdWCorrectionEnergy = 0.0;
   this->matrixCISdimension = 0;
   this->fockMatrix = NULL;
   this->energiesMO = NULL;
   this->orbitalElectronPopulation    = NULL;
   this->orbitalElectronPopulationCIS = NULL;
   this->atomicElectronPopulation     = NULL;
   this->atomicElectronPopulationCIS  = NULL;
   this->atomicUnpairedPopulationCIS  = NULL;
   this->overlapAOs = NULL;
   this->twoElecTwoCore = NULL;
   this->cartesianMatrix = NULL;
   this->electronicTransitionDipoleMoments = NULL;
   this->coreDipoleMoment = NULL;
   this->normalForceConstants = NULL;
   this->normalModes = NULL;
   this->matrixCIS = NULL;
   this->excitedEnergies = NULL;
   this->freeExcitonEnergiesCIS = NULL;
   this->matrixForce = NULL;

   //protected methods
   this->SetMessages();
   this->SetEnableAtomTypes();

   //private variables
   this->elecSCFEnergy = 0.0;
   this->bondingAdjustParameterK[0] = 1.000; //see (3.79) in J. A. Pople book
   this->bondingAdjustParameterK[1] = 0.750; //see (3.79) in J. A. Pople book
   this->gammaAB = NULL;
   //this->OutputLog("Cndo created\n");
}

Cndo2::~Cndo2(){
   MallocerFreer::GetInstance()->Free<double>(&this->fockMatrix, 
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Free<double>(&this->energiesMO, 
                                              this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Free<double>(&this->orbitalElectronPopulation, 
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Free<double>(&this->atomicElectronPopulation, 
                                              this->molecule->GetNumberAtoms());
   MallocerFreer::GetInstance()->Free<double>(&this->overlapAOs, 
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Free<double>(&this->cartesianMatrix, 
                                              CartesianType_end,
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs());
   int electronicTransitionDipoleMomentsDim = 1;
   if(Parameters::GetInstance()->RequiresCIS()){
      electronicTransitionDipoleMomentsDim += Parameters::GetInstance()->GetNumberExcitedStatesCIS();
   }
   MallocerFreer::GetInstance()->Free<double>(&this->electronicTransitionDipoleMoments, 
                                              electronicTransitionDipoleMomentsDim,
                                              electronicTransitionDipoleMomentsDim,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(&this->coreDipoleMoment, 
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(&this->gammaAB, 
                                              this->molecule->GetNumberAtoms(),
                                              this->molecule->GetNumberAtoms());
   //this->OutputLog("cndo deleted\n");
}

void Cndo2::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in cndo::Cndo2::DoSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in cndo::Cndo2::DoSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in cndo::Cndo2::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in cndo::Cndo2::ChecEnableAtomType: Non available atom is contained.\n";
   this->errorMessageAtomA = "Atom A is:\n";
   this->errorMessageAtomB = "Atom B is:\n";
   this->errorMessageAtomType = "\tatom type = ";
   this->errorMessageOrbitalType = "\torbital type = ";
   this->errorMessageCartesianType = "\tcartesian type = ";
   this->errorMessageMolecularIntegralElement
      = "Error in cndo::Cndo2::GetMolecularIntegralElement: Non available orbital is contained.\n";
   this->errorMessageGetDiatomCoreRepulsion2ndDerivativeNotImplemented
      = "Error in cndo::Cndo2::GetDiatomCoreRepulsion2ndDerivative: Second derivative is not implemented for CNDO2.\n";
   this->errorMessageGetGaussianCartesianMatrixBadOrbital 
      = "Error in cndo::Cndo2::GetGaussianCartesianMatrix: Untreatable orbital is contained in atom A or B.\n";
   this->errorMessageGetGaussianOverlapAOsBadOrbital 
      = "Error in cndo::Cndo2::GetGaussianOverlapAOs: Untreatable orbital is contained in atom A or B.\n";
   this->errorMessageGetGaussianOverlapAOs1stDerivativeOrbitalD 
      = "Error in cndo::Cndo2::GetGaussianOverlapAOs1stDerivative: d-orbital is not treatable. The d-orbital is contained in atom A or B.\n";
   this->errorMessageCISNotImplemented 
      = "Error in cndo::Cndo2: CIS is not implemented for CNDO2.\n";
   this->errorMessageCalcForceNotImplemented
      = "Error in cndo::Cndo2::CalcForce: Force is not available in CNDO2.\n";
   this->errorMessageGetElectronicEnergyNumberCISStates 
      = "\tNumber of calculated CIS states (excluding ground state) = ";
   this->errorMessageGetElectronicEnergySetElecState
      = "\tSet Electronic state = ";
   this->errorMessageGetElectronicEnergyEnergyNotCalculated
      = "Error in cndo::Cndo2::GetElectronicEnergy: Set electronic state is not calculated by CIS.\n";
   this->errorMessageGetElectronicEnergyNULLCISEnergy 
      = "Error in cndo::Cndo2::GetElectronicEnergy: excitedEnergies is NULL\n";
   this->errorMessageCalDiaOverlapAOsDiaFrameNullMatrix 
      = "Error in cndo::Cndo2::CalcDiatomicOverlapAOsInDiatomicFrame: diatomicOverlapAOs is NULL.\n";
   this->errorMessageCalcRotatingMatrixNullRotMatrix 
      = "Error in cndo::Cndo2::CalcRotatingMatrix: rotatingMatrix is NULL.\n";
   this->errorMessageRotDiaOverlapAOsToSpaceFrameNullDiaMatrix 
      = "Error in cndo::Cndo2::RotateDiatmicOverlapAOsToSpaceFrame diatomicOverlapAOs is NULL.\n";
   this->errorMessageRotDiaOverlapAOsToSpaceFrameNullRotMatrix 
      = "Error in cndo::Cndo2::RotateDiatmicOverlapAOsToSpaceFrame: rotatingMatrix is NULL.\n";
   this->errorMessageSetOverlapAOsElementNullDiaMatrix 
      = "Error in cndo::Cndo2::SetOverlapAOsElement: diatomicOverlapAOs is NULL.\n";
   this->errorMessageCalcElectronicTransitionDipoleMomentBadState
      = "Error in cndo::Cndo2::CalcElectronicTransitionDipoleMoment: Bad eigen state is set. In SCF module, the transition dipole moment of only between ground states can be calculated. Note taht state=0 means the ground state and other state = i means the i-th excited state in below.\n";
   this->errorMessageCalcFrequenciesNormalModesBadTheory
      = "Error in cndo::Cndo2::CalcFrequenciesNormalModesBadTheory: CNDO2 is not supported for frequency (normal mode) analysis.\n";
   this->errorMessageCalcOverlapAOsDifferentConfigurationsDiffAOs
      = "Error in cndo::Cndo2::CalcOverlapAOsDifferentConfigurations: Total number of AOs in lhs and rhs are different.\n";
   this->errorMessageCalcOverlapAOsDifferentConfigurationsDiffAtoms
      = "Error in cndo::Cndo2::CalcOverlapAOsDifferentConfigurations: Number Atoms in lhs and rhs are different.\n";
   this->errorMessageCalcOverlapAOsDifferentConfigurationsOverlapAOsNULL
      = "Error in cndo::Cndo2::CalcOverlapAOsDifferentConfigurations: ovelrapAOs is NULL.\n";
   this->errorMessageNonExcitedStates 
      = "Error in cndo::CNDO2::Excited states can not be calculated with CNDO2.\n";
   this->errorMessageLhs = "lhs: ";
   this->errorMessageRhs = "rhs: ";
   this->errorMessageFromState = "\tfrom state = ";
   this->errorMessageToState = "\tto state = ";
   this->messageSCFMetConvergence = "\n\n\n\t\tCNDO/2-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: CNDO/2-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: CNDO/2-SCF  **********\n\n\n";
   this->messageOmpElapsedTimeSCF = "\tElapsed time(omp) for the SCF = ";
   this->messageIterSCFTitle = "\t\t\t|  RMS density  | DIIS error | DIIS on/off | damping on/off |\n";
   this->messageIterSCF =      "\tSCF iter ";
   this->messageDiisApplied =     "on";
   this->messageDampingApplied =  "on";
   this->messageEnergyMO =     "\tEnergy of MO:";
   this->messageEnergyMOTitle = "\t\t\t| i-th | occ/unocc |  e[a.u.]  |  e[eV]  | \n";
   this->messageOcc   = "occ";
   this->messageUnOcc = "unocc";
   this->messageMullikenAtomsSCF   = "\tMulliken charge(SCF):";
   this->messageMullikenAtoms      = "\tMulliken charge:";
   this->messageMullikenAtomsTitle = "\t\t\t\t| k-th eigenstate | i-th atom | atom type | core charge[a.u.] | Mulliken charge[a.u.]| \n";
   this->messageUnpairedAtoms      = "\tUnpaired electron population:";
   this->messageUnpairedAtomsTitle = "\t\t\t\t| k-th eigenstate | i-th atom | atom type | Unpaired electron population[a.u.]| \n";
   this->messageElecEnergy = "\tElectronic energy(SCF):";
   this->messageNoteElecEnergy = "\tNote that this electronic energy includes core-repulsions.\n\n";
   this->messageNoteElecEnergyVdW = "\tNote that this electronic energy includes core-repulsions and vdW correction.\n\n";
   this->messageElecEnergyTitle = "\t\t\t\t|   [a.u.]   |   [eV]   |\n";
   this->messageUnitSec = "[s].";
   this->messageCoreRepulsionTitle = "\t\t\t\t|   [a.u.]   |   [eV]   |\n";
   this->messageCoreRepulsion = "\tCore repulsion energy:";
   this->messageVdWCorrectionTitle = "\t\t\t\t\t\t|   [a.u.]   |   [eV]   |\n";
   this->messageVdWCorrection = "\tEmpirical van der Waals correction:";
   this->messageElectronicDipoleMomentTitle = "\t\t\t\t\t|  x[a.u.]  |  y[a.u.]  |  z[a.u.]  |  magnitude[a.u.]  |\t\t|  x[debye]  |  y[debye]  |  z[debye]  |  magnitude[debye]  |\n";
   this->messageElectronicDipoleMoment = "\tElectronic Dipole moment(SCF):";
   this->messageCoreDipoleMomentTitle = "\t\t\t\t\t|  x[a.u.]  |  y[a.u.]  |  z[a.u.]  |  magnitude[a.u.]  |\t\t|  x[debye]  |  y[debye]  |  z[debye]  |  magnitude[debye]  |\n";
   this->messageCoreDipoleMoment = "\tCore Dipole moment:";
   this->messageTotalDipoleMomentTitle = "\t\t\t\t\t|   x[a.u.]   |   y[a.u.]   |   z[a.u.]   |  magnitude[a.u.]  |\t\t|  x[debye]  |  y[debye]  |  z[debye]  |  magnitude[debye]  |\n";
   this->messageTotalDipoleMoment = "\tTotal Dipole moment(SCF):";
   this->messageNormalModesTitle = "\t\t\t\t       |    normal frequencies          |   normalized normal mode ...\n";
   this->messageNormalModesUnitsNonMassWeighted = "\t\t\t\t| i-th |    [a.u.]    |    [cm-1]       |   [angst.] in non-mass-weighted coordinates ...\n";
   this->messageNormalModesUnitsMassWeighted = "\t\t\t\t| i-th |    [a.u.]    |    [cm-1]       |   [a.u.] in mass-weighted coordinates ...\n";
   this->messageNormalModesMassWeighted = "Normal mode(mw):";
   this->messageNormalModesNonMassWeighted = "Normal mode(nmw):";
   this->messageNormalModesImaginaryFrequencies = "\t\t\t\t\t'i' following the frequency means the imaginary frequency.\n";
}

void Cndo2::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(Li);
   //this->enableAtomTypes.push_back(Be);
   //this->enableAtomTypes.push_back(B);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   //this->enableAtomTypes.push_back(F);
   //this->enableAtomTypes.push_back(Na);
   //this->enableAtomTypes.push_back(Mg);
   //this->enableAtomTypes.push_back(Al);
   //this->enableAtomTypes.push_back(Si);
   //this->enableAtomTypes.push_back(P);
   this->enableAtomTypes.push_back(S);
   //this->enableAtomTypes.push_back(Cl);
}

TheoryType Cndo2::GetTheoryType() const{
   return this->theory;
}

void Cndo2::SetMolecule(Molecule* molecule){
   this->molecule = molecule;
   this->CheckNumberValenceElectrons(*molecule);
   this->CheckEnableAtomType(*molecule);

   // malloc
   MallocerFreer::GetInstance()->Malloc<double>(&this->fockMatrix,
                                                this->molecule->GetTotalNumberAOs(), 
                                                this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Malloc<double>(&this->energiesMO,
                                                this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Malloc<double>(&this->orbitalElectronPopulation,
                                                this->molecule->GetTotalNumberAOs(), 
                                                this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Malloc<double>(&this->atomicElectronPopulation,
                                                this->molecule->GetNumberAtoms());
   MallocerFreer::GetInstance()->Malloc<double>(&this->overlapAOs, 
                                                this->molecule->GetTotalNumberAOs(), 
                                                this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Malloc<double>(&this->cartesianMatrix, 
                                                CartesianType_end,
                                                this->molecule->GetTotalNumberAOs(), 
                                                this->molecule->GetTotalNumberAOs());
   int electronicTransitionDipoleMomentsDim = 1;
   if(Parameters::GetInstance()->RequiresCIS()){
      electronicTransitionDipoleMomentsDim += Parameters::GetInstance()->GetNumberExcitedStatesCIS();
   }
   MallocerFreer::GetInstance()->Malloc<double>(&this->electronicTransitionDipoleMoments, 
                                                electronicTransitionDipoleMomentsDim,
                                                electronicTransitionDipoleMomentsDim,
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(&this->coreDipoleMoment, 
                                                CartesianType_end);
   if(this->theory == CNDO2 || this->theory == INDO){
      MallocerFreer::GetInstance()->Malloc<double>(&this->gammaAB,
                                                   this->molecule->GetNumberAtoms(), 
                                                   this->molecule->GetNumberAtoms());
   }
}

void Cndo2::CheckNumberValenceElectrons(const Molecule& molecule) const{
   if(molecule.GetTotalNumberValenceElectrons() % 2 == 1){
      stringstream ss;
      ss << this->errorMessageOddTotalValenceElectrions << molecule.GetTotalNumberValenceElectrons() << "\n";
      throw MolDSException(ss.str());
   }
}

void Cndo2::CheckEnableAtomType(const Molecule& molecule) const{
   for(int i=0; i<molecule.GetNumberAtoms(); i++){
      AtomType atomType = molecule.GetAtom(i)->GetAtomType();
      bool enable = false;
      for(int j=0; j<this->enableAtomTypes.size(); j++){
         if(atomType == this->enableAtomTypes[j]){
            enable = true;
            break;
         }
      }
      if(!enable){
         stringstream ss;
         ss << this->errorMessageNotEnebleAtomType;
         ss << this->errorMessageAtomType << AtomTypeStr(atomType) << endl;
         throw MolDSException(ss.str());
      }
   }
}

void Cndo2::CalcCoreRepulsionEnergy(){
   double energy = 0.0;
   for(int i=0; i<this->molecule->GetNumberAtoms(); i++){
      for(int j=i+1; j<this->molecule->GetNumberAtoms(); j++){
         energy += this->GetDiatomCoreRepulsionEnergy(i, j);
      }
   }
   this->coreRepulsionEnergy = energy;
}

double Cndo2::GetDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   double distance = this->molecule->GetDistanceAtoms(indexAtomA, indexAtomB);
   return atomA.GetCoreCharge()*atomB.GetCoreCharge()/distance; 
}

// First derivative of diatomic core repulsion energy.
// This derivative is related to the coordinate of atomA.
double Cndo2::GetDiatomCoreRepulsion1stDerivative(int indexAtomA, int indexAtomB, 
                                                  CartesianType axisA) const{
   double value=0.0;
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   double distance = this->molecule->GetDistanceAtoms(indexAtomA, indexAtomB);
   value = atomA.GetCoreCharge()*atomB.GetCoreCharge();
   value *= (atomA.GetXyz()[axisA] - atomB.GetXyz()[axisA])/distance;
   value *= -1.0/pow(distance,2.0);
   return value;
}

// Second derivative of diatomic core repulsion energy.
// Both derivatives are related to the coordinate of atomA.
double Cndo2::GetDiatomCoreRepulsion2ndDerivative(int indexAtomA,
                                                  int indexAtomB, 
                                                  CartesianType axisA1,
                                                  CartesianType axisA2) const{
   stringstream ss;
   ss << this->errorMessageGetDiatomCoreRepulsion2ndDerivativeNotImplemented;
   throw MolDSException(ss.str());
}

// See (2) in [G_2004] ((11) in [G_2006])
void Cndo2::CalcVdWCorrectionEnergy(){
   double value = 0.0;
   for(int i=0; i<this->molecule->GetNumberAtoms(); i++){
      for(int j=i+1; j<this->molecule->GetNumberAtoms(); j++){
         value += this->GetDiatomVdWCorrectionEnergy(i, j);
      }
   }
   this->vdWCorrectionEnergy = value;
}

// See damping function in (2) in [G_2004] ((11) in [G_2006])
double Cndo2::GetVdwDampingValue(double vdWDistance, double distance) const{
   double dampingFactor = Parameters::GetInstance()->GetVdWDampingFactorSCF();
   return 1.0/(1.0+exp(-1.0*dampingFactor*(distance/vdWDistance - 1.0)));
}

// See damping function in (2) in [G_2004] ((11) in [G_2006])
double Cndo2::GetVdwDampingValue1stDerivative(double vdWDistance, double distance) const{
   double dampingFactor = Parameters::GetInstance()->GetVdWDampingFactorSCF();
   return (dampingFactor/vdWDistance)
         *exp(-1.0*dampingFactor*(distance/vdWDistance - 1.0))
         *pow(1.0+exp(-1.0*dampingFactor*(distance/vdWDistance - 1.0)),-2.0);
}

// See damping function in (2) in [G_2004] ((11) in [G_2006])
double Cndo2::GetVdwDampingValue2ndDerivative(double vdWDistance, double distance) const{
   double dampingFactor = Parameters::GetInstance()->GetVdWDampingFactorSCF();
   double exponent = -1.0*dampingFactor*(distance/vdWDistance - 1.0);
   double pre = dampingFactor/vdWDistance;
   double dominator = 1.0+exp(exponent);
   return 2.0*pow(dominator,-3.0)*pre*pre*exp(2.0*exponent) 
         -    pow(dominator,-2.0)*pre*pre*exp(    exponent);
}

// See (2) in [G_2004] ((11) in [G_2006])
double Cndo2::GetDiatomVdWCorrectionEnergy(int indexAtomA, int indexAtomB) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   double distance = this->molecule->GetDistanceAtoms(indexAtomA, indexAtomB);
   double vdWDistance = atomA.GetVdWRadii() + atomB.GetVdWRadii();
   double vdWCoefficients = 2.0*atomA.GetVdWCoefficient()*atomB.GetVdWCoefficient()
                           /(atomA.GetVdWCoefficient()+atomB.GetVdWCoefficient());
   double damping = this->GetVdwDampingValue(vdWDistance, distance);
   double scalingFactor = Parameters::GetInstance()->GetVdWScalingFactorSCF();
   return -1.0*scalingFactor*vdWCoefficients*pow(distance,-6.0)*damping;
}

// First derivative of the vdW correction related to the coordinate of atom A.
// See (2) in [G_2004] ((11) in [G_2006]).
double Cndo2::GetDiatomVdWCorrection1stDerivative(int indexAtomA, int indexAtomB, 
                                                  CartesianType axisA) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   double distance = this->molecule->GetDistanceAtoms(indexAtomA, indexAtomB);
   double vdWDistance = atomA.GetVdWRadii() + atomB.GetVdWRadii();
   double vdWCoefficients = 2.0*atomA.GetVdWCoefficient()*atomB.GetVdWCoefficient()
                           /(atomA.GetVdWCoefficient()+atomB.GetVdWCoefficient());
   double dampingFactor = Parameters::GetInstance()->GetVdWDampingFactorSCF();
   double damping = this->GetVdwDampingValue(vdWDistance, distance);
   double damping1stDerivative = this->GetVdwDampingValue1stDerivative(vdWDistance, distance);
   double value=0.0;
   value += 6.0*pow(distance,-7.0)*damping - pow(distance,-6.0)*damping1stDerivative;
   value *= vdWCoefficients;
   value *= Parameters::GetInstance()->GetVdWScalingFactorSCF();
   value *= (atomA.GetXyz()[axisA] - atomB.GetXyz()[axisA])/distance;
   return value;
}

// Second derivative of the vdW correction.
// Both derivative sare related to the coordinate of atom A.
// See (2) in [G_2004] ((11) in [G_2006]).
double Cndo2::GetDiatomVdWCorrection2ndDerivative(int indexAtomA, 
                                                     int indexAtomB, 
                                                     CartesianType axisA1,
                                                     CartesianType axisA2) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   double distance = this->molecule->GetDistanceAtoms(indexAtomA, indexAtomB);
   double dCartesian1 = atomA.GetXyz()[axisA1] - atomB.GetXyz()[axisA1];
   double dCartesian2 = atomA.GetXyz()[axisA2] - atomB.GetXyz()[axisA2];
   double vdWDistance = atomA.GetVdWRadii() + atomB.GetVdWRadii();
   double vdWScalingFacotor = Parameters::GetInstance()->GetVdWScalingFactorSCF();
   double vdWCoefficients = 2.0*atomA.GetVdWCoefficient()*atomB.GetVdWCoefficient()
                           /(atomA.GetVdWCoefficient()+atomB.GetVdWCoefficient());
   double dampingFactor = Parameters::GetInstance()->GetVdWDampingFactorSCF();
   double damping = this->GetVdwDampingValue(vdWDistance, distance);
   double damping1stDerivative = this->GetVdwDampingValue1stDerivative(vdWDistance, distance);
   double damping2ndDerivative = this->GetVdwDampingValue2ndDerivative(vdWDistance, distance);

   double temp1 = -6.0*pow(distance,-7.0)*damping 
                  +    pow(distance,-6.0)*damping1stDerivative;
   double temp2 = 42.0*pow(distance,-8.0)*damping 
                 -12.0*pow(distance,-7.0)*damping1stDerivative
                 +     pow(distance,-6.0)*damping2ndDerivative;

   double pre1=0.0;
   double pre2=0.0;
   if(axisA1 != axisA2){
      pre1 = -dCartesian1*dCartesian2/pow(distance,3.0);
      pre2 =  dCartesian1*dCartesian2/pow(distance,2.0);
   }
   else{
      pre1 = 1.0/distance - dCartesian1*dCartesian1/pow(distance,3.0);
      pre2 = pow(dCartesian1/distance,2.0);
   }

   double value= pre1*temp1 + pre2*temp2;
   value *= -1.0*vdWScalingFacotor*vdWCoefficients;
   return value;
}

/*******
 *
 * Call Cndo2::SetMolecule(Molecule* molecule) at least once, 
 * before this function is called.
 *
 *****/
void Cndo2::DoSCF(bool requiresGuess){
   this->OutputLog(this->messageStartSCF);
   double ompStartTime = omp_get_wtime();
   this->OutputLog(this->messageIterSCFTitle);
#ifdef MOLDS_DBG
   if(this->molecule == NULL){
      throw MolDSException(this->errorMessageMoleculeNotSet);
   }
#endif

   // temporary matrices for scf
   double**  oldOrbitalElectronPopulation = NULL;
   double*** diisStoredDensityMatrix = NULL;
   double*** diisStoredErrorVect = NULL;
   double**  diisErrorProducts = NULL;
   double*   diisErrorCoefficients = NULL;

   try{
      this->MallocSCFTemporaryMatrices(&oldOrbitalElectronPopulation,
                                       &diisStoredDensityMatrix,
                                       &diisStoredErrorVect,
                                       &diisErrorProducts,
                                       &diisErrorCoefficients);
      // calculate electron integral
      this->CalcGammaAB(this->gammaAB, *this->molecule);
      this->CalcOverlapAOs(this->overlapAOs, *this->molecule);
      this->CalcCartesianMatrixByGTOExpansion(this->cartesianMatrix, *this->molecule, STO6G);
      this->CalcTwoElecTwoCore(this->twoElecTwoCore, *this->molecule);

      // SCF
      double rmsDensity=0.0;
      int maxIterationsSCF = Parameters::GetInstance()->GetMaxIterationsSCF();
      bool isGuess=true;
      bool hasAppliedDIIS=false;
      bool hasAppliedDamping=false;
      double diisError=0.0;
      for(int iterationStep=0; iterationStep<maxIterationsSCF; iterationStep++){
         this->CalcAtomicElectronPopulation(this->atomicElectronPopulation, 
                                            this->orbitalElectronPopulation, 
                                            *this->molecule);
         this->UpdateOldOrbitalElectronPopulation(oldOrbitalElectronPopulation, 
                                                  this->orbitalElectronPopulation, 
                                                  this->molecule->GetTotalNumberAOs());
         isGuess = (iterationStep==0 && requiresGuess);
         this->CalcFockMatrix(this->fockMatrix, 
                              *this->molecule, 
                              this->overlapAOs, 
                              this->gammaAB,
                              this->orbitalElectronPopulation, 
                              this->atomicElectronPopulation,
                              this->twoElecTwoCore,
                              isGuess);

         // diagonalization of the Fock matrix
         bool calcEigenVectors = true;
         MolDS_wrappers::Lapack::GetInstance()->Dsyevd(this->fockMatrix, 
                                                       this->energiesMO, 
                                                       this->molecule->GetTotalNumberAOs(), 
                                                       calcEigenVectors);

         this->CalcOrbitalElectronPopulation(this->orbitalElectronPopulation, 
                                             *this->molecule, 
                                             this->fockMatrix);

         // check convergence
         bool hasConverged = this->SatisfyConvergenceCriterion(oldOrbitalElectronPopulation, 
                                                               this->orbitalElectronPopulation,
                                                               this->molecule->GetTotalNumberAOs(), 
                                                               &rmsDensity, 
                                                               iterationStep,
                                                               diisError,
                                                               hasAppliedDIIS,
                                                               hasAppliedDamping);
         if(hasConverged){
            this->OutputLog(this->messageSCFMetConvergence);
            this->CalcSCFProperties();
            this->OutputSCFResults();
            break;
         }
         else{
            if(!isGuess){ 
               this->DoDamp(rmsDensity, 
                            hasAppliedDamping,
                            this->orbitalElectronPopulation, 
                            oldOrbitalElectronPopulation, 
                            *this->molecule);
               this->DoDIIS(this->orbitalElectronPopulation,
                            oldOrbitalElectronPopulation,
                            diisStoredDensityMatrix,
                            diisStoredErrorVect,
                            diisErrorProducts,
                            diisErrorCoefficients,
                            diisError,
                            hasAppliedDIIS,
                            Parameters::GetInstance()->GetDiisNumErrorVectSCF(),
                            *this->molecule,
                            iterationStep);
            }
         }

         // SCF fails
         if(iterationStep==maxIterationsSCF-1){
            stringstream ss;
            ss << this->errorMessageSCFNotConverged << maxIterationsSCF << "\n";
            throw MolDSException(ss.str());
         }
      }
   }
   catch(MolDSException ex){
      this->FreeSCFTemporaryMatrices(&oldOrbitalElectronPopulation,
                                     &diisStoredDensityMatrix,
                                     &diisStoredErrorVect,
                                     &diisErrorProducts,
                                     &diisErrorCoefficients);

      throw ex;
   }
   this->FreeSCFTemporaryMatrices(&oldOrbitalElectronPopulation,
                                  &diisStoredDensityMatrix,
                                  &diisStoredErrorVect,
                                  &diisErrorProducts,
                                  &diisErrorCoefficients);

   double ompEndTime = omp_get_wtime();
   this->OutputLog(boost::format("%s%lf%s\n%s") % this->messageOmpElapsedTimeSCF.c_str()
                                                % (ompEndTime - ompStartTime)
                                                % this->messageUnitSec.c_str()
                                                % this->messageDoneSCF.c_str());

}

void Cndo2::CalcSCFProperties(){
   this->CalcAtomicElectronPopulation(this->atomicElectronPopulation, 
                                      this->orbitalElectronPopulation, 
                                      *this->molecule);
   this->CalcCoreRepulsionEnergy();
   if(Parameters::GetInstance()->RequiresVdWSCF()){
      this->CalcVdWCorrectionEnergy();
   }
   this->CalcElecSCFEnergy(&this->elecSCFEnergy, 
                          *this->molecule, 
                          this->energiesMO, 
                          this->fockMatrix, 
                          this->gammaAB,
                          this->coreRepulsionEnergy,
                          this->vdWCorrectionEnergy);
   this->CalcCoreDipoleMoment(this->coreDipoleMoment, *this->molecule);
   this->CalcElectronicDipoleMomentGroundState(this->electronicTransitionDipoleMoments, 
                                               this->cartesianMatrix,
                                               *this->molecule, 
                                               this->orbitalElectronPopulation,
                                               this->overlapAOs);
   const int groundState = 0;
   if(Parameters::GetInstance()->RequiresFrequencies() && 
      Parameters::GetInstance()->GetElectronicStateIndexFrequencies() == groundState){
      this->CalcNormalModes(this->normalModes, this->normalForceConstants, *this->molecule);
   }
}

void Cndo2::CalcNormalModes(double** normalModes, double* normalForceConstants, const Molecule& molecule) const{
   stringstream ss;
   ss << this->errorMessageCalcFrequenciesNormalModesBadTheory;
   throw MolDSException(ss.str());
}

double Cndo2::GetBondingAdjustParameterK(ShellType shellA, ShellType shellB) const{
   double value=1.0;
   if(shellA >= m || shellB >= m){
      return this->bondingAdjustParameterK[1];
   }
   return value;
}

void Cndo2::DoCIS(){
   stringstream ss;
   ss << this->errorMessageCISNotImplemented;
   throw MolDSException(ss.str());
}

void Cndo2::OutputCISResults() const{
   stringstream ss;
   ss << this->errorMessageCISNotImplemented;
   throw MolDSException(ss.str());
}

void Cndo2::CalcCISProperties(){
   stringstream ss;
   ss << this->errorMessageCISNotImplemented;
   throw MolDSException(ss.str());
}

// elecState=0 means ground state
double Cndo2::GetElectronicEnergy(int elecState) const{
   int groundState = 0;
   if(elecState==groundState){
      return this->elecSCFEnergy;
   }
   else{
#ifdef MOLDS_DBG
      if(this->excitedEnergies == NULL){
         throw MolDSException(this->errorMessageGetElectronicEnergyNULLCISEnergy);
      }
#endif
      int numberExcitedStates = Parameters::GetInstance()->GetNumberExcitedStatesCIS();
      if(numberExcitedStates < elecState){
         stringstream ss;
         ss << this->errorMessageGetElectronicEnergyEnergyNotCalculated;
         ss << this->errorMessageGetElectronicEnergySetElecState << elecState << endl;
         ss << errorMessageGetElectronicEnergyNumberCISStates << numberExcitedStates << endl;
         throw MolDSException(ss.str());
      }
      return this->elecSCFEnergy + this->excitedEnergies[elecState-1];
   }
}

double Cndo2::GetCoreRepulsionEnergy() const{
   return this->coreRepulsionEnergy;
}

double Cndo2::GetVdWCorrectionEnergy() const{
   return this->vdWCorrectionEnergy;
}

double*** Cndo2::GetForce(const vector<int>& elecStates){
   this->CalcForce(elecStates);
   return this->matrixForce;
}

double** Cndo2::GetForce(int elecState){
   vector<int> elecStates;
   elecStates.push_back(elecState);
   this->CalcForce(elecStates);
   return this->matrixForce[0];
}

void Cndo2::CalcTwoElecTwoCore(double****** twoElecTwoCore, 
                               const Molecule& molecule) const{
   // do nothing for CNDO, INDO, and ZINDO/S.
   // two electron two core integrals are not needed for CNDO, INDO, and ZINDO/S.
}

void Cndo2::CalcForce(const vector<int>& elecStates){
   stringstream ss;
   ss << this->errorMessageCalcForceNotImplemented;
   throw MolDSException(ss.str());
}

void Cndo2::FreeSCFTemporaryMatrices(double*** oldOrbitalElectronPopulation,
                                     double**** diisStoredDensityMatrix,
                                     double**** diisStoredErrorVect,
                                     double*** diisErrorProducts,
                                     double** diisErrorCoefficients) const{

   int diisNumErrorVect = Parameters::GetInstance()->GetDiisNumErrorVectSCF();
   MallocerFreer::GetInstance()->Free<double>(oldOrbitalElectronPopulation, 
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Free<double>(diisStoredDensityMatrix, 
                                              diisNumErrorVect, 
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Free<double>(diisStoredErrorVect, 
                                              diisNumErrorVect, 
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Free<double>(diisErrorProducts, 
                                              diisNumErrorVect+1,
                                              diisNumErrorVect+1);
   MallocerFreer::GetInstance()->Free<double>(diisErrorCoefficients,
                                              diisNumErrorVect+1);
}

void Cndo2::MallocSCFTemporaryMatrices(double*** oldOrbitalElectronPopulation,
                                       double**** diisStoredDensityMatrix,
                                       double**** diisStoredErrorVect,
                                       double*** diisErrorProducts,
                                       double** diisErrorCoefficients){

   int diisNumErrorVect = Parameters::GetInstance()->GetDiisNumErrorVectSCF();
   MallocerFreer::GetInstance()->Malloc<double>(oldOrbitalElectronPopulation, 
                                                this->molecule->GetTotalNumberAOs(), 
                                                this->molecule->GetTotalNumberAOs());
   if(0<diisNumErrorVect){
      MallocerFreer::GetInstance()->Malloc<double>(diisStoredDensityMatrix,
                                                   diisNumErrorVect, 
                                                   this->molecule->GetTotalNumberAOs(), 
                                                   this->molecule->GetTotalNumberAOs());
      MallocerFreer::GetInstance()->Malloc<double>(diisStoredErrorVect,
                                                   diisNumErrorVect, 
                                                   this->molecule->GetTotalNumberAOs(), 
                                                   this->molecule->GetTotalNumberAOs());
      MallocerFreer::GetInstance()->Malloc<double>(diisErrorProducts, diisNumErrorVect+1, diisNumErrorVect+1);
      MallocerFreer::GetInstance()->Malloc<double>(diisErrorCoefficients, diisNumErrorVect+1);
   }
}

/***
 *
 *  see ref. [P_1980] for diis methods.
 *
 */
void Cndo2::DoDIIS(double** orbitalElectronPopulation,
                   double const* const* oldOrbitalElectronPopulation,
                   double*** diisStoredDensityMatrix,
                   double*** diisStoredErrorVect,
                   double**  diisErrorProducts,
                   double*   diisErrorCoefficients,
                   double&   diisError,
                   bool&     hasAppliedDIIS,
                   int       diisNumErrorVect,
                   const     Molecule& molecule,
                   int       step) const{
   int totalNumberAOs = molecule.GetTotalNumberAOs();
   double diisStartError = Parameters::GetInstance()->GetDiisStartErrorSCF();
   double diisEndError = Parameters::GetInstance()->GetDiisEndErrorSCF();

   if( 0 < diisNumErrorVect){
      for(int m=0; m<diisNumErrorVect-1; m++){
         stringstream ompErrors;
#pragma omp parallel for schedule(auto)
         for(int j=0; j<totalNumberAOs; j++){
            try{
               for(int k=0; k<totalNumberAOs; k++){
                  diisStoredDensityMatrix[m][j][k] = diisStoredDensityMatrix[m+1][j][k];
                  diisStoredErrorVect[m][j][k] = diisStoredErrorVect[m+1][j][k];
               }
            }
            catch(MolDSException ex){
#pragma omp critical
               ex.Serialize(ompErrors);
            }
         }
         // Exception throwing for omp-region
         if(!ompErrors.str().empty()){
            throw MolDSException::Deserialize(ompErrors);
         }
      }
      {
         stringstream ompErrors;
#pragma omp parallel for schedule(auto)
         for(int j=0; j<totalNumberAOs; j++){
            try{
               for(int k=0; k<totalNumberAOs; k++){
                  diisStoredDensityMatrix[diisNumErrorVect-1][j][k] = orbitalElectronPopulation[j][k];
                  diisStoredErrorVect[diisNumErrorVect-1][j][k] = orbitalElectronPopulation[j][k] 
                                                                 -oldOrbitalElectronPopulation[j][k];
                     
               }
            }
            catch(MolDSException ex){
#pragma omp critical
               ex.Serialize(ompErrors);
            }
         }
         // Exception throwing for omp-region
         if(!ompErrors.str().empty()){
            throw MolDSException::Deserialize(ompErrors);
         }
      }
      for(int mi=0; mi<diisNumErrorVect-1; mi++){
         for(int mj=0; mj<diisNumErrorVect-1; mj++){
            diisErrorProducts[mi][mj] = diisErrorProducts[mi+1][mj+1];
         }
      }
               
      for(int mi=0; mi<diisNumErrorVect; mi++){
         double tempErrorProduct = 0.0;
         stringstream ompErrors;
#pragma omp parallel for schedule(auto) reduction(+:tempErrorProduct)
         for(int j=0; j<totalNumberAOs; j++){
            try{
               for(int k=0; k<totalNumberAOs; k++){
                  tempErrorProduct += diisStoredErrorVect[mi][j][k]
                                     *diisStoredErrorVect[diisNumErrorVect-1][j][k];
               }
            }
            catch(MolDSException ex){
#pragma omp critical
               ex.Serialize(ompErrors);
            }
         }
         // Exception throwing for omp-region
         if(!ompErrors.str().empty()){
            throw MolDSException::Deserialize(ompErrors);
         }
         diisErrorProducts[mi][diisNumErrorVect-1] = tempErrorProduct;
         diisErrorProducts[diisNumErrorVect-1][mi] = tempErrorProduct;
         diisErrorProducts[mi][diisNumErrorVect] = -1.0;
         diisErrorProducts[diisNumErrorVect][mi] = -1.0;
         diisErrorCoefficients[mi] = 0.0;
      }
      diisErrorProducts[diisNumErrorVect][diisNumErrorVect] = 0.0;
      diisErrorCoefficients[diisNumErrorVect] = -1.0;

      diisError = 0.0;
      for(int j=0; j<totalNumberAOs; j++){
         for(int k=0; k<totalNumberAOs; k++){
            diisError = max(diisError, fabs(diisStoredErrorVect[diisNumErrorVect-1][j][k]));
         }
      }

      hasAppliedDIIS = false;
      if(diisNumErrorVect <= step && diisEndError<diisError && diisError<diisStartError){
         hasAppliedDIIS = true;
         try{
            MolDS_wrappers::Lapack::GetInstance()->Dsysv(diisErrorProducts, 
                                                         diisErrorCoefficients, 
                                                         diisNumErrorVect+1);
         }catch(MolDSException ex){
            if(ex.HasKey(LapackInfo) && ex.GetKeyValue<int>(LapackInfo) > 0){
               // DIIS matrix is now singular, so not taking DIIS step.
               hasAppliedDIIS = false;
               return;
            }
            else{
               throw ex;
            }
         }
         for(int j=0; j<totalNumberAOs; j++){
            for(int k=0; k<totalNumberAOs; k++){
               orbitalElectronPopulation[j][k] = 0.0;
               for(int m=0; m<diisNumErrorVect; m++){
                  orbitalElectronPopulation[j][k] += diisErrorCoefficients[m]*diisStoredDensityMatrix[m][j][k];
               }
            }
         }
      }
   }
}

void Cndo2::DoDamp(double rmsDensity, 
                   bool&  hasAppliedDamping,
                   double** orbitalElectronPopulation, 
                   double const* const* oldOrbitalElectronPopulation, 
                   const Molecule& molecule) const{
   double dampingThresh = Parameters::GetInstance()->GetDampingThreshSCF();
   double dampingWeight = Parameters::GetInstance()->GetDampingWeightSCF();
   hasAppliedDamping = false;
   if(0.0 < dampingWeight && dampingThresh < rmsDensity){
      hasAppliedDamping = true;
      stringstream ompErrors;
#pragma omp parallel for schedule(auto)
      for(int j=0; j<molecule.GetTotalNumberAOs(); j++){
         try{
            for(int k=0; k<molecule.GetTotalNumberAOs(); k++){
               orbitalElectronPopulation[j][k] *= (1.0 - dampingWeight);
               orbitalElectronPopulation[j][k] += dampingWeight*oldOrbitalElectronPopulation[j][k];
            }
         }
         catch(MolDSException ex){
#pragma omp critical
            ex.Serialize(ompErrors);
         }
      }
      // Exception throwing for omp-region
      if(!ompErrors.str().empty()){
         throw MolDSException::Deserialize(ompErrors);
      }
   } 

}

void Cndo2::OutputMOEnergies() const{
   double eV2AU = Parameters::GetInstance()->GetEV2AU();
   this->OutputLog(this->messageEnergyMOTitle);
   for(int mo=0; mo<this->molecule->GetTotalNumberAOs(); mo++){
      string occUnOcc = this->messageUnOcc;
      if(mo < this->molecule->GetTotalNumberValenceElectrons()/2){
         occUnOcc = this->messageOcc;
      }
      this->OutputLog(boost::format("%s\t%d\t%s\t%e\t%e\n") % this->messageEnergyMO
                                                            % mo
                                                            % occUnOcc
                                                            % this->energiesMO[mo] 
                                                            % (this->energiesMO[mo]/eV2AU) );
   }
   this->OutputLog("\n");
}

void Cndo2::OutputSCFEnergies() const{
   double eV2AU = Parameters::GetInstance()->GetEV2AU();

   // electronic energy
   this->OutputLog(this->messageElecEnergyTitle);
   this->OutputLog(boost::format("%s\t%e\t%e\n") % this->messageElecEnergy
                                                 % this->elecSCFEnergy
                                                 % (this->elecSCFEnergy/eV2AU));
   if(Parameters::GetInstance()->RequiresVdWSCF()){
      this->OutputLog(this->messageNoteElecEnergyVdW);
   }
   else{
      this->OutputLog(this->messageNoteElecEnergy);
   }

   // output core repulsion energy
   this->OutputLog(this->messageCoreRepulsionTitle);
   this->OutputLog(boost::format("%s\t%e\t%e\n\n") % this->messageCoreRepulsion
                                                   % this->coreRepulsionEnergy 
                                                   % (this->coreRepulsionEnergy/eV2AU));

   // output van der Waals correction 
   if(Parameters::GetInstance()->RequiresVdWSCF()){
      this->OutputLog(this->messageVdWCorrectionTitle);
      this->OutputLog(boost::format("%s\t%e\t%e\n\n") % this->messageVdWCorrection
                                                      % this->vdWCorrectionEnergy 
                                                      % (this->vdWCorrectionEnergy/eV2AU));
   }
}

void Cndo2::OutputSCFDipole() const{
   int groundState=0;
   double debye2AU = Parameters::GetInstance()->GetDebye2AU();
   double magnitude = 0.0;
   double temp = 0.0;

   // output total dipole moment 
   temp = 0.0;
   temp += pow(this->electronicTransitionDipoleMoments[groundState][groundState][XAxis]+this->coreDipoleMoment[XAxis],2.0);
   temp += pow(this->electronicTransitionDipoleMoments[groundState][groundState][YAxis]+this->coreDipoleMoment[YAxis],2.0);
   temp += pow(this->electronicTransitionDipoleMoments[groundState][groundState][ZAxis]+this->coreDipoleMoment[ZAxis],2.0);
   magnitude = sqrt(temp);
   this->OutputLog(this->messageTotalDipoleMomentTitle);
   this->OutputLog(boost::format("%s\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n\n") 
      % this->messageTotalDipoleMoment
      % (this->electronicTransitionDipoleMoments[groundState][groundState][XAxis]+this->coreDipoleMoment[XAxis])
      % (this->electronicTransitionDipoleMoments[groundState][groundState][YAxis]+this->coreDipoleMoment[YAxis])
      % (this->electronicTransitionDipoleMoments[groundState][groundState][ZAxis]+this->coreDipoleMoment[ZAxis])
      % magnitude
      % ((this->electronicTransitionDipoleMoments[groundState][groundState][XAxis]+this->coreDipoleMoment[XAxis])/debye2AU)
      % ((this->electronicTransitionDipoleMoments[groundState][groundState][YAxis]+this->coreDipoleMoment[YAxis])/debye2AU)
      % ((this->electronicTransitionDipoleMoments[groundState][groundState][ZAxis]+this->coreDipoleMoment[ZAxis])/debye2AU)
      % (magnitude/debye2AU));

   // output electronic dipole moment 
   temp = 0.0;
   temp += pow(this->electronicTransitionDipoleMoments[groundState][groundState][XAxis],2.0);
   temp += pow(this->electronicTransitionDipoleMoments[groundState][groundState][YAxis],2.0);
   temp += pow(this->electronicTransitionDipoleMoments[groundState][groundState][ZAxis],2.0);
   magnitude = sqrt(temp);
   this->OutputLog(this->messageElectronicDipoleMomentTitle);
   this->OutputLog(boost::format("%s\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n\n") 
      % this->messageElectronicDipoleMoment
      % this->electronicTransitionDipoleMoments[groundState][groundState][XAxis]
      % this->electronicTransitionDipoleMoments[groundState][groundState][YAxis]
      % this->electronicTransitionDipoleMoments[groundState][groundState][ZAxis]
      % magnitude
      % (this->electronicTransitionDipoleMoments[groundState][groundState][XAxis]/debye2AU)
      % (this->electronicTransitionDipoleMoments[groundState][groundState][YAxis]/debye2AU)
      % (this->electronicTransitionDipoleMoments[groundState][groundState][ZAxis]/debye2AU)
      % (magnitude/debye2AU));

   // output core dipole moment 
   temp = 0.0;
   temp += pow(this->coreDipoleMoment[XAxis],2.0);
   temp += pow(this->coreDipoleMoment[YAxis],2.0);
   temp += pow(this->coreDipoleMoment[ZAxis],2.0);
   magnitude = sqrt(temp);
   this->OutputLog(this->messageCoreDipoleMomentTitle);
   this->OutputLog(boost::format("%s\t\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n\n") 
      % this->messageCoreDipoleMoment
      % this->coreDipoleMoment[XAxis]
      % this->coreDipoleMoment[YAxis]
      % this->coreDipoleMoment[ZAxis]
      % magnitude
      % (this->coreDipoleMoment[XAxis]/debye2AU)
      % (this->coreDipoleMoment[YAxis]/debye2AU)
      % (this->coreDipoleMoment[ZAxis]/debye2AU)
      % (magnitude/debye2AU));
}

void Cndo2::OutputSCFMulliken() const{
   int groundState = 0;
   this->OutputLog(this->messageMullikenAtomsTitle);
   for(int a=0; a<this->molecule->GetNumberAtoms(); a++){
      Atom* atom = this->molecule->GetAtom(a);
      this->OutputLog(boost::format("%s\t%d\t%d\t%s\t%e\t%e\n") % this->messageMullikenAtomsSCF
                                                                % groundState
                                                                % a
                                                                % AtomTypeStr(atom->GetAtomType())
                                                                % atom->GetCoreCharge()
                                                                % (atom->GetCoreCharge()-atomicElectronPopulation[a]));
   }
   this->OutputLog("\n");
}

void Cndo2::OutputNormalModes(double const* const* normalModes, 
                              double const* normalForceConstants, 
                              const Molecule& molecule) const{

   int hessianDim = CartesianType_end*molecule.GetNumberAtoms();
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double kayser2AU = Parameters::GetInstance()->GetKayser2AU();

   // output in mass-weighted coordinates
   this->OutputLog(this->messageNormalModesTitle);
   this->OutputLog(this->messageNormalModesUnitsMassWeighted);
   for(int i=0; i<hessianDim; i++){
      // normal frequencies
      if(normalForceConstants[i]>0)
         this->OutputLog(boost::format("\t%s\t%d\t%e\t%e\t") % this->messageNormalModesMassWeighted
                                                             % i
                                                             % sqrt(normalForceConstants[i])
                                                             % (sqrt(normalForceConstants[i])/kayser2AU));
      else
         this->OutputLog(boost::format("\t%s\t%d\t%ei\t%ei\t") % this->messageNormalModesMassWeighted
                                                               % i
                                                               % sqrt(fabs(normalForceConstants[i]))
                                                               % (sqrt(fabs(normalForceConstants[i]))/kayser2AU));
      // normal modes
      for(int a=0; a<molecule.GetNumberAtoms(); a++){
         const double sqrtCoreMass = sqrt(molecule.GetAtom(a)->GetCoreMass());
         for(int j=XAxis; j<CartesianType_end; j++){
            int hessianIndex = CartesianType_end*a+j;
            this->OutputLog(boost::format("\t%e") % normalModes[i][hessianIndex]);
         }
      }
      this->OutputLog("\n");
   }
   this->OutputLog(this->messageNormalModesImaginaryFrequencies);
   this->OutputLog("\n");

   // output in non-mass-weighted coordinates
   this->OutputLog(this->messageNormalModesTitle);
   this->OutputLog(this->messageNormalModesUnitsNonMassWeighted);
   for(int i=0; i<hessianDim; i++){
      // normal frequencies
      if(normalForceConstants[i]>0)
         this->OutputLog(boost::format("\t%s\t%d\t%e\t%e\t") % this->messageNormalModesNonMassWeighted
                                                             % i
                                                             % sqrt(normalForceConstants[i])
                                                             % (sqrt(normalForceConstants[i])/kayser2AU));
      else
         this->OutputLog(boost::format("\t%s\t%d\t%ei\t%ei\t") % this->messageNormalModesNonMassWeighted
                                                               % i
                                                               % sqrt(fabs(normalForceConstants[i]))
                                                               % (sqrt(fabs(normalForceConstants[i]))/kayser2AU));

      double normSquare=0.0;
      for(int a=0; a<molecule.GetNumberAtoms(); a++){
         const double sqrtCoreMass = sqrt(molecule.GetAtom(a)->GetCoreMass());
         for(int j=XAxis; j<CartesianType_end; j++){
            int hessianIndex = CartesianType_end*a+j;
            normSquare += pow(normalModes[i][hessianIndex]/(sqrtCoreMass*ang2AU),2.0);
         }
      }
      double norm = sqrt(normSquare);

      // normal modes
      for(int a=0; a<molecule.GetNumberAtoms(); a++){
         const double sqrtCoreMass = sqrt(molecule.GetAtom(a)->GetCoreMass());
         for(int j=XAxis; j<CartesianType_end; j++){
            int hessianIndex = CartesianType_end*a+j;
            this->OutputLog(boost::format("\t%e") % (normalModes[i][hessianIndex]/(sqrtCoreMass*ang2AU*norm)));
         }
      }
      this->OutputLog("\n");
   }
   this->OutputLog(this->messageNormalModesImaginaryFrequencies);
   this->OutputLog("\n");
}

void Cndo2::OutputSCFResults() const{
   this->OutputMOEnergies();
   this->OutputSCFEnergies();
   this->OutputSCFDipole();
   this->OutputSCFMulliken();
   // ToDo: output eigen-vectors of the Hartree Fock matrix

   // Normal modes and frequencies  
   const int groundState = 0;
   if(Parameters::GetInstance()->RequiresFrequencies() && 
      Parameters::GetInstance()->GetElectronicStateIndexFrequencies() == groundState){
      this->OutputNormalModes(this->normalModes, this->normalForceConstants, *this->molecule);
   }

   // output MOs
   if(Parameters::GetInstance()->RequiresMOPlot()){
      MolDS_base_loggers::MOLogger* moLogger = new MolDS_base_loggers::MOLogger(*this->molecule, 
                                                                                this->fockMatrix, 
                                                                                this->theory);
      moLogger->DrawMO(*(Parameters::GetInstance()->GetIndecesMOPlot()));
      delete moLogger;
   }
}

void Cndo2::CalcElecSCFEnergy(double* elecSCFEnergy, 
                             const Molecule& molecule, 
                             double const* energiesMO, 
                             double const* const* fockMatrix, 
                             double const* const* gammaAB, 
                             double coreRepulsionEnergy,
                             double vdWCorrectionEnergy) const{
   double electronicEnergy = 0.0;
   // use density matrix for electronic energy
   int totalNumberAOs = this->molecule->GetTotalNumberAOs();
   double** fMatrix = NULL;
   double** hMatrix = NULL;
   double** dammyOrbitalElectronPopulation = NULL;
   double* dammyAtomicElectronPopulation = NULL;

   try{
      MallocerFreer::GetInstance()->Malloc<double>(&fMatrix, totalNumberAOs, totalNumberAOs);
      MallocerFreer::GetInstance()->Malloc<double>(&hMatrix, totalNumberAOs,totalNumberAOs);
      MallocerFreer::GetInstance()->Malloc<double>(&dammyOrbitalElectronPopulation,
                                                   totalNumberAOs, 
                                                   totalNumberAOs);
      MallocerFreer::GetInstance()->Malloc<double>(&dammyAtomicElectronPopulation,
                                                   molecule.GetNumberAtoms());
      bool isGuess = false;
      this->CalcFockMatrix(fMatrix, 
                           molecule, 
                           this->overlapAOs, 
                           this->gammaAB,
                           this->orbitalElectronPopulation, 
                           this->atomicElectronPopulation,
                           this->twoElecTwoCore,
                           isGuess);
      this->CalcFockMatrix(hMatrix, 
                           molecule, 
                           this->overlapAOs, 
                           this->gammaAB,
                           dammyOrbitalElectronPopulation, 
                           dammyAtomicElectronPopulation,
                           this->twoElecTwoCore,
                           isGuess);

      for(int i=0; i<totalNumberAOs; i++){
         for(int j=i+1; j<totalNumberAOs; j++){
            fMatrix[j][i] = fMatrix[i][j];
            hMatrix[j][i] = hMatrix[i][j];
         }
      }

      for(int i=0; i<totalNumberAOs; i++){
         for(int j=0; j<totalNumberAOs; j++){
            electronicEnergy += this->orbitalElectronPopulation[j][i]*
                                 (fMatrix[i][j] + hMatrix[i][j]);
         }
      }
      electronicEnergy *= 0.5;
   }
   catch(MolDSException ex){
      this->FreeElecEnergyMatrices(&fMatrix, 
                                   &hMatrix, 
                                   &dammyOrbitalElectronPopulation, 
                                   &dammyAtomicElectronPopulation );
      throw ex;
   }
   this->FreeElecEnergyMatrices(&fMatrix, 
                                &hMatrix, 
                                &dammyOrbitalElectronPopulation, 
                                &dammyAtomicElectronPopulation );

   // use two electrons integrals for electronic energy
   /*
   for(int mo=0; mo<molecule.GetTotalNumberValenceElectrons()/2; mo++){
      electronicEnergy += 2.0*energiesMO[mo];
   }

   for(int moA=0; moA<molecule.GetTotalNumberValenceElectrons()/2; moA++){
      for(int moB=0; moB<molecule.GetTotalNumberValenceElectrons()/2; moB++){

         electronicEnergy -= 2.0*this->GetMolecularIntegralElement(moA, moA, moB, moB, 
                                                              molecule, fockMatrix, gammaAB);
         electronicEnergy += 1.0*this->GetMolecularIntegralElement(moA, moB, moB, moA, 
                                                              molecule, fockMatrix, gammaAB);
      }
   }
   */

   *elecSCFEnergy = electronicEnergy + coreRepulsionEnergy + vdWCorrectionEnergy;
}

void Cndo2::FreeElecEnergyMatrices(double*** fMatrix, 
                                   double*** hMatrix, 
                                   double*** dammyOrbitalElectronPopulation, 
                                   double**  dammyAtomicElectronPopulation ) const{
   int totalNumberAOs = this->molecule->GetTotalNumberAOs();
   MallocerFreer::GetInstance()->Free<double>(fMatrix, totalNumberAOs, totalNumberAOs);
   MallocerFreer::GetInstance()->Free<double>(hMatrix, totalNumberAOs, totalNumberAOs);
   MallocerFreer::GetInstance()->Free<double>(dammyOrbitalElectronPopulation, 
                                              totalNumberAOs,
                                              totalNumberAOs);
   MallocerFreer::GetInstance()->Free<double>(dammyAtomicElectronPopulation,
                                              this->molecule->GetNumberAtoms());
}

// The order of moI, moJ, moK, moL is consistent with Eq. (9) in [RZ_1973]
double Cndo2::GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                          const Molecule& molecule, 
                                          double const* const* fockMatrix, 
                                          double const* const* gammaAB) const{
   double value = 0.0;
   for(int A=0; A<molecule.GetNumberAtoms(); A++){
      const Atom& atomA = *molecule.GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();

      for(int B=0; B<molecule.GetNumberAtoms(); B++){
         const Atom& atomB = *molecule.GetAtom(B);
         int firstAOIndexB = atomB.GetFirstAOIndex();
         int lastAOIndexB  = atomB.GetLastAOIndex();
         double gamma = gammaAB[A][B];

         for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
            for(int nu=firstAOIndexB; nu<=lastAOIndexB; nu++){

               value += gamma*fockMatrix[moI][mu]*fockMatrix[moJ][mu]*fockMatrix[moK][nu]*fockMatrix[moL][nu];
            }
         }

      }
   }
   return value;
}

void Cndo2::UpdateOldOrbitalElectronPopulation(double** oldOrbitalElectronPopulation, 
                                               double const* const* orbitalElectronPopulation,
                                               int numberAOs) const{
   for(int i=0; i<numberAOs; i++){
      for(int j=0; j<numberAOs; j++){
         oldOrbitalElectronPopulation[i][j] = orbitalElectronPopulation[i][j];
      }
   }
}

bool Cndo2::SatisfyConvergenceCriterion(double const* const * oldOrbitalElectronPopulation,
                                        double const* const * orbitalElectronPopulation,
                                        int     numberAOs,
                                        double* rmsDensity,
                                        int     times,
                                        double  diisError,
                                        bool    hasAppliedDIIS,
                                        bool    hasAppliedDamping) const{
   bool satisfy = false;
   double change = 0.0;
   stringstream ompErrors;
#pragma omp parallel for schedule(auto) reduction(+:change)
   for(int i=0; i<numberAOs; i++){
      try{
         for(int j=0; j<numberAOs; j++){
            change += pow(oldOrbitalElectronPopulation[i][j] - orbitalElectronPopulation[i][j], 2.0);
         }
      }
      catch(MolDSException ex){
#pragma omp critical
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
   change /= numberAOs*numberAOs;
   *rmsDensity = sqrt(change);
 
   string diisOnOff    = hasAppliedDIIS    ? this->messageDiisApplied    : "";
   string dampingOnOff = hasAppliedDamping ? this->messageDampingApplied : "";
   this->OutputLog(boost::format("%s%d\t%e\t%e\t%s\t\t%s\n") % this->messageIterSCF.c_str()
                                                       % times
                                                       % *rmsDensity
                                                       % diisError
                                                       % diisOnOff
                                                       % dampingOnOff);

   if(*rmsDensity < Parameters::GetInstance()->GetThresholdSCF()){
      satisfy = true;
   }

   return satisfy; 
}

/*********
 *
 *
 * Upper right part of the Fock matrix is only caluculated.
 *
 *
 * ******/
void Cndo2::CalcFockMatrix(double** fockMatrix, 
                           const Molecule& molecule, 
                           double const* const* overlapAOs, 
                           double const* const* gammaAB,
                           double const* const* orbitalElectronPopulation, 
                           double const* atomicElectronPopulation,
                           double const* const* const* const* const* const* twoElecTwoCore, 
                           bool isGuess) const{
   MallocerFreer::GetInstance()->Initialize<double>(fockMatrix, 
                                                    molecule.GetTotalNumberAOs(), 
                                                    molecule.GetTotalNumberAOs());
   int totalNumberAtoms=molecule.GetNumberAtoms();
   stringstream ompErrors;
#pragma omp parallel for schedule(auto) 
   for(int A=0; A<totalNumberAtoms; A++){
      try{
        const Atom& atomA = *molecule.GetAtom(A);
         int firstAOIndexA = atomA.GetFirstAOIndex();
         int lastAOIndexA  = atomA.GetLastAOIndex();
         for(int B=A; B<totalNumberAtoms; B++){
            const Atom& atomB = *molecule.GetAtom(B);
            int firstAOIndexB = atomB.GetFirstAOIndex();
            int lastAOIndexB  = atomB.GetLastAOIndex();
            for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
               for(int nu=firstAOIndexB; nu<=lastAOIndexB; nu++){
                  if(mu == nu){
                     // diagonal part
                     fockMatrix[mu][mu] = this->GetFockDiagElement(atomA, 
                                                                   A, 
                                                                   mu, 
                                                                   molecule, 
                                                                   gammaAB,
                                                                   orbitalElectronPopulation, 
                                                                   atomicElectronPopulation,
                                                                   twoElecTwoCore,
                                                                   isGuess);
                  }
                  else if(mu < nu){
                     // upper right part
                     fockMatrix[mu][nu] = this->GetFockOffDiagElement(atomA, 
                                                                      atomB,
                                                                      A, 
                                                                      B, 
                                                                      mu, 
                                                                      nu, 
                                                                      molecule, 
                                                                      gammaAB,
                                                                      overlapAOs,
                                                                      orbitalElectronPopulation, 
                                                                      twoElecTwoCore,
                                                                      isGuess);
                  }
                  else{
                     // lower left part (not calculated)
                  }

               }
            }
         }
      }
      catch(MolDSException ex){
#pragma omp critical
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
   /*  
   this->OutputLog("fock matrix\n");
   for(int o=0; o<this->molecule.GetTotalNumberAOs(); o++){
      for(int p=0; p<this->molecule.GetTotalNumberAOs(); p++){
         this->OutputLog(boost::format("%lf\t") % fockMatrix[o][p]);
      }
      this->OutputLog("\n");
   }
   this->OutputLog("\n\n");
   */
}

double Cndo2::GetFockDiagElement(const Atom& atomA, 
                                 int indexAtomA, 
                                 int mu, 
                                 const Molecule& molecule, 
                                 double const* const* gammaAB,
                                 double const* const* orbitalElectronPopulation, 
                                 double const* atomicElectronPopulation,
                                 double const* const* const* const* const* const* twoElecTwoCore, 
                                 bool isGuess) const{
   double value;
   int firstAOIndexA = atomA.GetFirstAOIndex();
   value = atomA.GetCoreIntegral(atomA.GetValence(mu-firstAOIndexA), 
                                 gammaAB[indexAtomA][indexAtomA], 
                                 isGuess, this->theory);
   if(!isGuess){
      double temp = atomicElectronPopulation[indexAtomA] 
                   -0.5*orbitalElectronPopulation[mu][mu];
      value += temp*gammaAB[indexAtomA][indexAtomA];

      temp = 0.0;
      for(int BB=0; BB<molecule.GetNumberAtoms(); BB++){
         if(BB != indexAtomA){
            const Atom& atomBB = *molecule.GetAtom(BB);
            temp += ( atomicElectronPopulation[BB] - atomBB.GetCoreCharge()  )
                     *gammaAB[indexAtomA][BB];
         }
      }
      value += temp;
   }

   return value;
}

double Cndo2::GetFockOffDiagElement(const Atom& atomA, 
                                    const Atom& atomB, 
                                    int indexAtomA, 
                                    int indexAtomB, 
                                    int mu, 
                                    int nu, 
                                    const Molecule& molecule, 
                                    double const* const* gammaAB, 
                                    double const* const* overlapAOs,
                                    double const* const* orbitalElectronPopulation, 
                                    double const* const* const* const* const* const* twoElecTwoCore, 
                                    bool isGuess) const{
   double value;
   double K = this->GetBondingAdjustParameterK(atomA.GetValenceShellType(), atomB.GetValenceShellType());
   double bondParameter = 0.5*K*(atomA.GetBondingParameter() + atomB.GetBondingParameter()); 
   value =  bondParameter*overlapAOs[mu][nu];
   if(!isGuess){
      value -= 0.5*orbitalElectronPopulation[mu][nu]*gammaAB[indexAtomA][indexAtomB];
   }
   return value;
}

void Cndo2::TransposeFockMatrixMatrix(double** transposedFockMatrix) const{
   for(int i=0; i<this->molecule->GetTotalNumberAOs(); i++){
      for(int j=0; j<this->molecule->GetTotalNumberAOs(); j++){
         transposedFockMatrix[j][i] = this->fockMatrix[i][j];
      }
   }
}

void Cndo2::CalcOrbitalElectronPopulation(double** orbitalElectronPopulation, 
                                          const Molecule& molecule, 
                                          double const* const* fockMatrix) const{
   const int totalNumberAOs = molecule.GetTotalNumberAOs();
   const int numberTotalValenceElectrons = molecule.GetTotalNumberValenceElectrons();

   MallocerFreer::GetInstance()->Initialize<double>(orbitalElectronPopulation, totalNumberAOs, totalNumberAOs);

   bool isMatrixAColumnMajor = false;
   bool isMatrixATransposed = true;
   bool isLowerTriangularPartMatrixCUsed = false;
   double alpha = 2.0, beta = 0.0;
   MolDS_wrappers::Blas::GetInstance()->Dsyrk(totalNumberAOs, numberTotalValenceElectrons/2,
                                              isMatrixAColumnMajor,
                                              isMatrixATransposed,
                                              isLowerTriangularPartMatrixCUsed,
                                              alpha, fockMatrix,
                                              beta, orbitalElectronPopulation);

   /* 
   this->OutputLog("orbital population\n");
   for(int mu=0; mu<totalNumberAOs; mu++){
      for(int nu=0; nu<totalNumberAOs; nu++){
         this->OutputLog(boost::format("%lf\t") % orbitalElectronPopulation[mu][nu]);
      }
      this->OutputLog("\n");
   }
   this->OutputLog("\n");
   */
}

void Cndo2::CalcAtomicElectronPopulation(double* atomicElectronPopulation,
                                         double const* const* orbitalElectronPopulation, 
                                         const Molecule& molecule) const{
   int totalNumberAtoms = molecule.GetNumberAtoms();
   MallocerFreer::GetInstance()->Initialize<double>(atomicElectronPopulation, totalNumberAtoms);

   int firstAOIndex = 0;
   int numberAOs = 0;
   for(int A=0; A<totalNumberAtoms; A++){
      firstAOIndex = molecule.GetAtom(A)->GetFirstAOIndex();
      numberAOs = molecule.GetAtom(A)->GetValenceSize();
      for(int i=firstAOIndex; i<firstAOIndex+numberAOs; i++){
         atomicElectronPopulation[A] += orbitalElectronPopulation[i][i];
      }
      //this->OutputLog(boost::format("P_AA[%d]=%lf\n") % A % atomicElectronPopulation[A]);
   }
}

// calculate gammaAB matrix. (B.56) and (B.62) in J. A. Pople book.
void Cndo2::CalcGammaAB(double** gammaAB, const Molecule& molecule) const{
   int totalAtomNumber = molecule.GetNumberAtoms();
   stringstream ompErrors;
#pragma omp parallel for schedule(auto) 
   for(int A=0; A<totalAtomNumber; A++){
      try{
         const Atom& atomA = *molecule.GetAtom(A);
         int na = atomA.GetValenceShellType() + 1;
         double orbitalExponentA = atomA.GetOrbitalExponent(
                                         atomA.GetValenceShellType(), s, this->theory);
         for(int B=A; B<totalAtomNumber; B++){
            const Atom& atomB = *molecule.GetAtom(B);
            int nb = atomB.GetValenceShellType() + 1;
            double orbitalExponentB = atomB.GetOrbitalExponent(
                                            atomB.GetValenceShellType(), s, this->theory);

            double value = 0.0;
            double R = molecule.GetDistanceAtoms(A, B);
            double temp = 0.0;
            if(R>0.0){
               // (B.56)
               value = pow(0.5*R, 2.0*na);
               value *= this->GetReducedOverlapAOs(2*na-1, 0, 2.0*orbitalExponentA*R, 0);

               for(int l=1; l<=2*nb; l++){
                  temp = 0.0;
                  temp = l;
                  temp *= pow(2.0*orbitalExponentB, 2*nb-l);
                  temp /= Factorial(2*nb-l)*2.0*nb;
                  temp *= pow(0.5*R, 2.0*nb-l+2.0*na);
                  temp *= this->GetReducedOverlapAOs(2*na-1, 
                                                     2*nb-l, 
                                                     2.0*orbitalExponentA*R, 
                                                     2.0*orbitalExponentB*R);
                  value -= temp;
               }

               value *= pow(2.0*orbitalExponentA, 2.0*na+1.0);
               value /= Factorial(2*na);
            }
            else{
               // (B.62)
               value =  Factorial(2*na-1);
               value /= pow(2.0*orbitalExponentA, 2.0*na);

               for(int l=1; l<=2*nb; l++){
                  temp = l;
                  temp *= pow(2.0*orbitalExponentB, 2*nb-l);
                  temp *= Factorial(2*na+2*nb-l-1);
                  temp /= Factorial(2*nb-l);
                  temp /= 2.0*nb;
                  temp /= pow( 2.0*orbitalExponentA + 2.0*orbitalExponentB, 2.0*(na+nb)-l );
                  value -= temp;
               }
               value *= pow(2.0*orbitalExponentA, 2.0*na+1);
               value /= Factorial(2*na);
            }
            gammaAB[A][B] = value;
         }
      }
      catch(MolDSException ex){
#pragma omp critical
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }

#pragma omp parallel for schedule(auto)
   for(int A=0; A<totalAtomNumber; A++){
      for(int B=0; B<A; B++){
         gammaAB[A][B] = gammaAB[B][A];
      }
   }
   
   /* 
   this->OutputLog("gamma matrix\n");
   for(int A=0; A<totalAtomNumber; A++){
      for(int B=0; B<totalAtomNumber; B++){
         this->OutputLog(boost::format("gammaAB[%d][%d]=%lf\n") % A % B % gammaAB[A][B]);
      }
      this->OutputLog("\n");
   }
   this->OutputLog("\n");
   */

}

void Cndo2::CalcCoreDipoleMoment(double* coreDipoleMoment,
                                 const Molecule& molecule) const{

   for(int i=0; i<CartesianType_end; i++){
      coreDipoleMoment[i] = 0.0;
      for(int A=0; A<molecule.GetNumberAtoms(); A++){
         coreDipoleMoment[i] += molecule.GetAtom(A)->GetCoreCharge()
                               *(molecule.GetAtom(A)->GetXyz()[i] - molecule.GetXyzCOC()[i]);
      }
   }
}

void Cndo2::CalcElectronicDipoleMomentGroundState(double*** electronicTransitionDipoleMoments,
                                                  double const* const* const* cartesianMatrix,
                                                  const Molecule& molecule,
                                                  double const* const* orbitalElectronPopulation,
                                                  double const* const* overlapAOs) const{
   int groundState = 0;
   this->CalcElectronicTransitionDipoleMoment(electronicTransitionDipoleMoments[groundState][groundState],
                                              groundState,
                                              groundState,
                                              NULL,
                                              NULL,
                                              cartesianMatrix,
                                              molecule,
                                              orbitalElectronPopulation,
                                              overlapAOs,
                                              NULL);
}

void Cndo2::CalcElectronicTransitionDipoleMoment(double* transitionDipoleMoment,
                                                 int to, int from,
                                                 double const* const* fockMatrix,
                                                 double const* const* matrixCIS,
                                                 double const* const* const* cartesianMatrix,
                                                 const MolDS_base::Molecule& molecule, 
                                                 double const* const* orbitalElectronPopulation,
                                                 double const* const* overlapAOs,
                                                 double const* groundStateDipole) const{
   int groundState = 0;
   if(from == groundState && to == groundState){
      double const* centerOfDipole = molecule.GetXyzCOC();
      int totalAONumber = molecule.GetTotalNumberAOs();
      transitionDipoleMoment[XAxis] = 0.0;
      transitionDipoleMoment[YAxis] = 0.0;
      transitionDipoleMoment[ZAxis] = 0.0;
      transitionDipoleMoment[XAxis] -= MolDS_wrappers::Blas::GetInstance()->Ddot(totalAONumber*totalAONumber,
                                                                                 &orbitalElectronPopulation[0][0],
                                                                                 &cartesianMatrix[XAxis][0][0]);
      transitionDipoleMoment[YAxis] -= MolDS_wrappers::Blas::GetInstance()->Ddot(totalAONumber*totalAONumber,
                                                                                 &orbitalElectronPopulation[0][0],
                                                                                 &cartesianMatrix[YAxis][0][0]);
      transitionDipoleMoment[ZAxis] -= MolDS_wrappers::Blas::GetInstance()->Ddot(totalAONumber*totalAONumber,
                                                                                 &orbitalElectronPopulation[0][0],
                                                                                 &cartesianMatrix[ZAxis][0][0]);
      // set orign of dipole
      double temp = MolDS_wrappers::Blas::GetInstance()->Ddot(totalAONumber*totalAONumber,
                                                              &orbitalElectronPopulation[0][0],
                                                              &overlapAOs[0][0]);
      transitionDipoleMoment[XAxis] += centerOfDipole[XAxis]*temp;
      transitionDipoleMoment[YAxis] += centerOfDipole[YAxis]*temp;
      transitionDipoleMoment[ZAxis] += centerOfDipole[ZAxis]*temp;
   }
   else{
      stringstream ss;
      ss << this->errorMessageCalcElectronicTransitionDipoleMomentBadState;
      ss << this->errorMessageFromState << from << endl;
      ss << this->errorMessageToState << to << endl;
      throw MolDSException(ss.str());
   }
}

// calculate Cartesian matrix between atomic orbitals. 
// The analytic Cartesian matrix is calculated with Gaussian expansion technique written in [DY_1977]
void Cndo2::CalcCartesianMatrixByGTOExpansion(double*** cartesianMatrix, 
                                              const Molecule& molecule, 
                                              STOnGType stonG) const{
   int totalAONumber = molecule.GetTotalNumberAOs();
   int totalAtomNumber = molecule.GetNumberAtoms();

   stringstream ompErrors;
#pragma omp parallel for schedule(auto) 
   for(int A=0; A<totalAtomNumber; A++){
      try{
         const Atom& atomA = *molecule.GetAtom(A);
         int firstAOIndexAtomA = atomA.GetFirstAOIndex();
         for(int B=0; B<totalAtomNumber; B++){
            const Atom& atomB = *molecule.GetAtom(B);
            int firstAOIndexAtomB = atomB.GetFirstAOIndex();
            for(int a=0; a<atomA.GetValenceSize(); a++){
               for(int b=0; b<atomB.GetValenceSize(); b++){
                  int mu = firstAOIndexAtomA + a;      
                  int nu = firstAOIndexAtomB + b;      
                  this->CalcCartesianMatrixElementsByGTOExpansion(cartesianMatrix[XAxis][mu][nu], 
                                                                  cartesianMatrix[YAxis][mu][nu],
                                                                  cartesianMatrix[ZAxis][mu][nu],
                                                                  atomA, a, atomB, b, stonG);
               }
            }
            
         }
      }
      catch(MolDSException ex){
#pragma omp critical
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
}

// Calculate elements of Cartesian matrix between atomic orbitals. 
// The analytic Cartesian matrix is calculated with Gaussian expansion technique written in [DY_1977]
void Cndo2::CalcCartesianMatrixElementsByGTOExpansion(double& xComponent,
                                                      double& yComponent,
                                                      double& zComponent,
                                                      const Atom& atomA, int valenceIndexA, 
                                                      const Atom& atomB, int valenceIndexB,
                                                      STOnGType stonG) const{
   xComponent=0.0;
   yComponent=0.0;
   zComponent=0.0;
   ShellType shellTypeA = atomA.GetValenceShellType();
   ShellType shellTypeB = atomB.GetValenceShellType();
   OrbitalType valenceOrbitalA = atomA.GetValence(valenceIndexA);
   OrbitalType valenceOrbitalB = atomB.GetValence(valenceIndexB);
   double orbitalExponentA = atomA.GetOrbitalExponent(atomA.GetValenceShellType(), 
                                                      valenceOrbitalA, 
                                                      this->theory);
   double orbitalExponentB = atomB.GetOrbitalExponent(atomB.GetValenceShellType(), 
                                                      valenceOrbitalB, 
                                                      this->theory);
   double gaussianExponentA = 0.0;
   double gaussianExponentB = 0.0;
   double overlapSASB = 0.0;
   double rAB = sqrt( pow(atomA.GetXyz()[XAxis]-atomB.GetXyz()[XAxis], 2.0) 
                     +pow(atomA.GetXyz()[YAxis]-atomB.GetXyz()[YAxis], 2.0) 
                     +pow(atomA.GetXyz()[ZAxis]-atomB.GetXyz()[ZAxis], 2.0) );
   double temp  = 0.0;
   double tempX = 0.0;
   double tempY = 0.0;
   double tempZ = 0.0;
   for(int i=0; i<=stonG; i++){
      for(int j=0; j<=stonG; j++){
         temp = GTOExpansionSTO::GetInstance()->GetCoefficient(stonG, 
                                                               shellTypeA, 
                                                               valenceOrbitalA, 
                                                               i); 
         temp *= GTOExpansionSTO::GetInstance()->GetCoefficient(stonG, 
                                                                shellTypeB, 
                                                                valenceOrbitalB, 
                                                                j); 
         gaussianExponentA = pow(orbitalExponentA, 2.0) *
                             GTOExpansionSTO::GetInstance()->GetExponent(stonG, 
                                                                         shellTypeA, 
                                                                         valenceOrbitalA, 
                                                                         i);
         gaussianExponentB = pow(orbitalExponentB, 2.0) *
                             GTOExpansionSTO::GetInstance()->GetExponent(stonG, 
                                                                         shellTypeB, 
                                                                         valenceOrbitalB, 
                                                                         j);
         overlapSASB = this->GetGaussianOverlapAOsSASB(gaussianExponentA, gaussianExponentB, rAB);
         tempX = this->GetGaussianCartesianMatrix(atomA.GetAtomType(), valenceOrbitalA, gaussianExponentA, atomA.GetXyz(),
                                                  atomB.GetAtomType(), valenceOrbitalB, gaussianExponentB, atomB.GetXyz(), 
                                                  rAB, overlapSASB,
                                                  XAxis);
         tempY = this->GetGaussianCartesianMatrix(atomA.GetAtomType(), valenceOrbitalA, gaussianExponentA, atomA.GetXyz(),
                                                  atomB.GetAtomType(), valenceOrbitalB, gaussianExponentB, atomB.GetXyz(), 
                                                  rAB, overlapSASB, 
                                                  YAxis);
         tempZ = this->GetGaussianCartesianMatrix(atomA.GetAtomType(), valenceOrbitalA, gaussianExponentA, atomA.GetXyz(),
                                                  atomB.GetAtomType(), valenceOrbitalB, gaussianExponentB, atomB.GetXyz(), 
                                                  rAB, overlapSASB,
                                                  ZAxis);
         xComponent += temp*tempX;
         yComponent += temp*tempY;
         zComponent += temp*tempZ;
      }
   }
}

// calculate gaussian Caretesian integrals. 
double Cndo2::GetGaussianCartesianMatrix(AtomType atomTypeA, 
                                         OrbitalType valenceOrbitalA, 
                                         double gaussianExponentA, 
                                         double const* xyzA,
                                         AtomType atomTypeB, 
                                         OrbitalType valenceOrbitalB, 
                                         double gaussianExponentB,
                                         double const* xyzB,
                                         double rAB,
                                         CartesianType axis) const{
   double overlapSASB = this->GetGaussianOverlapAOsSASB(gaussianExponentA, gaussianExponentB, rAB);
   return this->GetGaussianCartesianMatrix(atomTypeA, valenceOrbitalA, gaussianExponentA, xyzA,
                                           atomTypeB, valenceOrbitalB, gaussianExponentB, xyzB,
                                           rAB, overlapSASB, axis);
}

// calculate gaussian Caretesian integrals. 
double Cndo2::GetGaussianCartesianMatrix(AtomType atomTypeA, 
                                         OrbitalType valenceOrbitalA, 
                                         double gaussianExponentA, 
                                         double const* xyzA,
                                         AtomType atomTypeB, 
                                         OrbitalType valenceOrbitalB, 
                                         double gaussianExponentB,
                                         double const* xyzB,
                                         double rAB,
                                         double overlapSASB,
                                         CartesianType axis) const{

   double value = 0.0;
   double beta = gaussianExponentA + gaussianExponentB;
   double dxyz[CartesianType_end] = {xyzA[XAxis] - xyzB[XAxis], 
                                     xyzA[YAxis] - xyzB[YAxis], 
                                     xyzA[ZAxis] - xyzB[ZAxis]};
   if(valenceOrbitalA == s && valenceOrbitalB == s){
      value = gaussianExponentA*xyzA[axis] + gaussianExponentB*xyzB[axis];
      value /= beta;
      value *= overlapSASB;
      return value;
   }
   else if( (valenceOrbitalA == s && axis == XAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == s && axis == XAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == s && axis == XAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == s && axis == XAxis && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == s && axis == XAxis && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == s && axis == XAxis && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == s && axis == XAxis && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == s && axis == XAxis && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == s && axis == YAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == s && axis == YAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == s && axis == YAxis && valenceOrbitalB == pz) ||
            (valenceOrbitalA == s && axis == YAxis && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == s && axis == YAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == s && axis == YAxis && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == s && axis == YAxis && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == s && axis == YAxis && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == s && axis == ZAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == s && axis == ZAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == s && axis == ZAxis && valenceOrbitalB == pz) ||
            (valenceOrbitalA == s && axis == ZAxis && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == s && axis == ZAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == s && axis == ZAxis && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == s && axis == ZAxis && valenceOrbitalB == dxxyy) ||
            (valenceOrbitalA == s && axis == ZAxis && valenceOrbitalB == dzz) ){
      OrbitalType pOrbital;
      if(axis == XAxis){
         pOrbital = px;
      }
      else if(axis == YAxis){
         pOrbital = py;
      }
      else if(axis == ZAxis){
         pOrbital = pz;
      }
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      double overlapAOs2 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       pOrbital, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = overlapAOs2/(2.0*sqrt(gaussianExponentA))+xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == px    && axis == XAxis && valenceOrbitalB == s) || 
            (valenceOrbitalA == py    && axis == XAxis && valenceOrbitalB == s) || 
            (valenceOrbitalA == pz    && axis == XAxis && valenceOrbitalB == s) || 
            (valenceOrbitalA == dxy   && axis == XAxis && valenceOrbitalB == s) || 
            (valenceOrbitalA == dyz   && axis == XAxis && valenceOrbitalB == s) || 
            (valenceOrbitalA == dzx   && axis == XAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == dxxyy && axis == XAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == dzz   && axis == XAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == px    && axis == YAxis && valenceOrbitalB == s) || 
            (valenceOrbitalA == py    && axis == YAxis && valenceOrbitalB == s) || 
            (valenceOrbitalA == pz    && axis == YAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == dxy   && axis == YAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == dyz   && axis == YAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == dzx   && axis == YAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == dxxyy && axis == YAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == dzz   && axis == YAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == px    && axis == ZAxis && valenceOrbitalB == s) || 
            (valenceOrbitalA == py    && axis == ZAxis && valenceOrbitalB == s) || 
            (valenceOrbitalA == pz    && axis == ZAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == dxy   && axis == ZAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == dyz   && axis == ZAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == dzx   && axis == ZAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == dxxyy && axis == ZAxis && valenceOrbitalB == s) ||
            (valenceOrbitalA == dzz   && axis == ZAxis && valenceOrbitalB == s) ){
      OrbitalType pOrbital;
      if(axis == XAxis){
         pOrbital = px;
      }
      else if(axis == YAxis){
         pOrbital = py;
      }
      else if(axis == ZAxis){
         pOrbital = pz;
      }
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      double overlapAOs2 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       pOrbital, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = overlapAOs2/(2.0*sqrt(gaussianExponentB))+xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == XAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == py && axis == YAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == pz && axis == ZAxis && valenceOrbitalB == valenceOrbitalA) ){
      double temp1 = gaussianExponentA*xyzA[axis] + gaussianExponentB*xyzB[axis];
      double temp2 = gaussianExponentA*xyzA[axis] - gaussianExponentA*xyzB[axis];
      double temp3 = gaussianExponentB*xyzA[axis] - gaussianExponentB*xyzB[axis];
      value = 0.5*(temp1+temp2-temp3);
      value -= temp1*temp2*temp3*pow(beta,-1.0);
      value *= 4.0*sqrt(gaussianExponentA*gaussianExponentB)*pow(beta,-2.0);
      value *= overlapSASB;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == YAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == px && axis == ZAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == py && axis == XAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == py && axis == ZAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == pz && axis == XAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == pz && axis == YAxis && valenceOrbitalB == valenceOrbitalA) ){
      CartesianType piDirection;
      if(valenceOrbitalA == px){
         piDirection = XAxis;
      }
      else if(valenceOrbitalA == py){
         piDirection = YAxis;
      }
      else if(valenceOrbitalA == pz){
         piDirection = ZAxis;
      }
      double temp1 = gaussianExponentA*xyzA[piDirection] - gaussianExponentA*xyzB[piDirection];
      double temp2 = gaussianExponentB*xyzA[piDirection] - gaussianExponentB*xyzB[piDirection];
      value = 0.5 - temp1*temp2*pow(beta,-1.0);
      value *= gaussianExponentA*xyzA[axis] + gaussianExponentB*xyzB[axis];
      value *= 4.0*sqrt(gaussianExponentA*gaussianExponentB)*pow(beta,-2.0);
      value *= overlapSASB;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == YAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == px && axis == ZAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == py && axis == XAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == py && axis == ZAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == pz && axis == XAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == pz && axis == YAxis && valenceOrbitalB == py) ){
      CartesianType piDirectionA;
      if(valenceOrbitalA == px){
         piDirectionA = XAxis;
      }
      else if(valenceOrbitalA == py){
         piDirectionA = YAxis;
      }
      else if(valenceOrbitalA == pz){
         piDirectionA = ZAxis;
      }
      double temp1 = gaussianExponentA*xyzA[axis] + gaussianExponentB*xyzB[axis];
      double temp2 = gaussianExponentA*xyzA[axis] - gaussianExponentA*xyzB[axis];
      value = 0.5 + temp1*temp2*pow(beta,-1.0);
      value *= gaussianExponentB*xyzA[piDirectionA] - gaussianExponentB*xyzB[piDirectionA];
      value *= -4.0*sqrt(gaussianExponentA*gaussianExponentB)*pow(beta,-2.0);
      value *= overlapSASB;
      return value;
   }
   else if( (valenceOrbitalA == py && axis == YAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == py && axis == YAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == px && axis == XAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == px && axis == XAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == pz && axis == ZAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == pz && axis == ZAxis && valenceOrbitalB == py) ){
      CartesianType piDirectionB;
      if(valenceOrbitalB == px){
         piDirectionB = XAxis;
      }
      else if(valenceOrbitalB == py){
         piDirectionB = YAxis;
      }
      else if(valenceOrbitalB == pz){
         piDirectionB = ZAxis;
      }
      double temp1 = gaussianExponentA*xyzA[axis] + gaussianExponentB*xyzB[axis];
      double temp2 = gaussianExponentB*xyzA[axis] - gaussianExponentB*xyzB[axis];
      value = 0.5 - temp1*temp2*pow(beta,-1.0);
      value *= gaussianExponentA*xyzA[piDirectionB] - gaussianExponentA*xyzB[piDirectionB];
      value *= 4.0*sqrt(gaussianExponentA*gaussianExponentB)*pow(beta,-2.0);
      value *= overlapSASB;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == YAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == py && axis == ZAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == pz && axis == XAxis && valenceOrbitalB == py) ||
            (valenceOrbitalA == px && axis == ZAxis && valenceOrbitalB == py) ||
            (valenceOrbitalA == pz && axis == YAxis && valenceOrbitalB == px) ||
            (valenceOrbitalA == py && axis == XAxis && valenceOrbitalB == pz) ){
      CartesianType piDirectionA;
      CartesianType piDirectionB;
      if(valenceOrbitalA == px){
         piDirectionA = XAxis;
      }
      else if(valenceOrbitalA == py){
         piDirectionA = YAxis;
      }
      else if(valenceOrbitalA == pz){
         piDirectionA = ZAxis;
      }
      if(valenceOrbitalB == px){
         piDirectionB = XAxis;
      }
      else if(valenceOrbitalB == py){
         piDirectionB = YAxis;
      }
      else if(valenceOrbitalB == pz){
         piDirectionB = ZAxis;
      }
      double temp1 = gaussianExponentB*xyzA[piDirectionA] - gaussianExponentB*xyzB[piDirectionA];
      double temp2 = gaussianExponentA*xyzA[axis]         + gaussianExponentB*xyzB[axis];
      double temp3 = gaussianExponentA*xyzA[piDirectionB] - gaussianExponentA*xyzB[piDirectionB];
      value = -4.0*sqrt(gaussianExponentA*gaussianExponentB)*pow(beta,-3.0);
      value *= temp1*temp2*temp3;
      value *= overlapSASB;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == XAxis && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == py && axis == YAxis && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == py && axis == YAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == pz && axis == ZAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == pz && axis == ZAxis && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == px && axis == XAxis && valenceOrbitalB == dzx) ){
      CartesianType anotherAxis;
      if(valenceOrbitalB == dxy){
         if(axis == XAxis){
            anotherAxis = YAxis;
         }
         else{
            anotherAxis = XAxis;
         }
      }
      else if(valenceOrbitalB == dyz){
         if(axis == YAxis){
            anotherAxis = ZAxis;
         }
         else{
            anotherAxis = YAxis;
         }
      }
      else if(valenceOrbitalB == dzx){
         if(axis == ZAxis){
            anotherAxis = XAxis;
         }
         else{
            anotherAxis = ZAxis;
         }
      }
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*gaussianExponentA*dxyz[axis]
             -gaussianExponentB*dxyz[axis]
             +pow(gaussianExponentB,2.0)*gaussianExponentA*pow(dxyz[axis],3.0)/beta;
      value *= 8.0*pow(gaussianExponentA, 1.5)*gaussianExponentB*pow(beta,-3.0)*dxyz[anotherAxis]*overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy && axis == XAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == dxy && axis == YAxis && valenceOrbitalB == py) ||
            (valenceOrbitalA == dyz && axis == YAxis && valenceOrbitalB == py) ||
            (valenceOrbitalA == dyz && axis == ZAxis && valenceOrbitalB == pz) ||
            (valenceOrbitalA == dzx && axis == ZAxis && valenceOrbitalB == pz) ||
            (valenceOrbitalA == dzx && axis == XAxis && valenceOrbitalB == px) ){
      CartesianType anotherAxis;
      if(valenceOrbitalA == dxy){
         if(axis == XAxis){
            anotherAxis = YAxis;
         }
         else{
            anotherAxis = XAxis;
         }
      }
      else if(valenceOrbitalA == dyz){
         if(axis == YAxis){
            anotherAxis = ZAxis;
         }
         else{
            anotherAxis = YAxis;
         }
      }
      else if(valenceOrbitalA == dzx){
         if(axis == ZAxis){
            anotherAxis = XAxis;
         }
         else{
            anotherAxis = ZAxis;
         }
      }
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*gaussianExponentB*dxyz[axis]
             -gaussianExponentA*dxyz[axis]
             +pow(gaussianExponentA,2.0)*gaussianExponentB*pow(dxyz[axis],3.0)/beta;
      value *= 8.0*pow(gaussianExponentB, 1.5)*gaussianExponentA*pow(beta,-3.0)*dxyz[anotherAxis]*overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == pz && axis == ZAxis && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == px && axis == XAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == py && axis == YAxis && valenceOrbitalB == dzx) ){
      CartesianType anotherAxis1;
      CartesianType anotherAxis2;
      if(axis == XAxis){
         anotherAxis1 = YAxis;
         anotherAxis2 = ZAxis;
      }
      else if(axis == YAxis){
         anotherAxis1 = XAxis;
         anotherAxis2 = ZAxis;
      }
      else if(axis == ZAxis){
         anotherAxis1 = XAxis;
         anotherAxis2 = YAxis;
      }
      double overlapAOs1=0.0; 
      overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                valenceOrbitalA, 
                                                gaussianExponentA, 
                                                atomTypeB, 
                                                valenceOrbitalB, 
                                                gaussianExponentB, 
                                                dxyz[XAxis], 
                                                dxyz[YAxis], 
                                                dxyz[ZAxis], 
                                                rAB,
                                                overlapSASB);
      value = 0.5+pow(gaussianExponentB*dxyz[axis], 2.0)/beta;
      value *= 8.0*pow(gaussianExponentA, 2.5)*gaussianExponentB*pow(beta, -3.0)*overlapSASB;
      value *= dxyz[anotherAxis1]*dxyz[anotherAxis2];
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy && axis == ZAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dyz && axis == XAxis && valenceOrbitalB == px) ||
            (valenceOrbitalA == dzx && axis == YAxis && valenceOrbitalB == py) ){
      CartesianType anotherAxis1;
      CartesianType anotherAxis2;
      if(axis == XAxis){
         anotherAxis1 = YAxis;
         anotherAxis2 = ZAxis;
      }
      else if(axis == YAxis){
         anotherAxis1 = XAxis;
         anotherAxis2 = ZAxis;
      }
      else if(axis == ZAxis){
         anotherAxis1 = XAxis;
         anotherAxis2 = YAxis;
      }
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB, 
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5+pow(gaussianExponentA*dxyz[axis], 2.0)/beta;
      value *= 8.0*pow(gaussianExponentB, 2.5)*gaussianExponentA*pow(beta, -3.0)*overlapSASB;
      value *= dxyz[anotherAxis1]*dxyz[anotherAxis2];
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == YAxis && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == px && axis == ZAxis && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == py && axis == XAxis && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == py && axis == ZAxis && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == pz && axis == XAxis && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == pz && axis == YAxis && valenceOrbitalB == dxy) || 
            (valenceOrbitalA == px && axis == YAxis && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == px && axis == ZAxis && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == py && axis == XAxis && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == py && axis == ZAxis && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == pz && axis == XAxis && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == pz && axis == YAxis && valenceOrbitalB == dyz) || 
            (valenceOrbitalA == px && axis == YAxis && valenceOrbitalB == dzx) || 
            (valenceOrbitalA == px && axis == ZAxis && valenceOrbitalB == dzx) || 
            (valenceOrbitalA == py && axis == XAxis && valenceOrbitalB == dzx) || 
            (valenceOrbitalA == py && axis == ZAxis && valenceOrbitalB == dzx) || 
            (valenceOrbitalA == pz && axis == XAxis && valenceOrbitalB == dzx) || 
            (valenceOrbitalA == pz && axis == YAxis && valenceOrbitalB == dzx) || 
            (valenceOrbitalA == px && axis == YAxis && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == px && axis == ZAxis && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == py && axis == XAxis && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == py && axis == ZAxis && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == pz && axis == XAxis && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == pz && axis == YAxis && valenceOrbitalB == dxxyy) || 
            (valenceOrbitalA == px && axis == YAxis && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == px && axis == ZAxis && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == py && axis == XAxis && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == py && axis == ZAxis && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == pz && axis == XAxis && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == pz && axis == YAxis && valenceOrbitalB == dzz) ){
      OrbitalType dOrbital;
      if( (valenceOrbitalA == py && axis == XAxis) || 
          (valenceOrbitalA == px && axis == YAxis) ){
         dOrbital = dxy;
      }
      else if( (valenceOrbitalA == py && axis == ZAxis) || 
               (valenceOrbitalA == pz && axis == YAxis) ){
         dOrbital = dyz;
      }
      else if( (valenceOrbitalA == px && axis == ZAxis) || 
               (valenceOrbitalA == pz && axis == XAxis) ){
         dOrbital = dzx;
      }
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      double overlapAOs2 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       dOrbital, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = overlapAOs2/(2.0*sqrt(gaussianExponentA))+xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy   && axis == YAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == dxy   && axis == ZAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == dxy   && axis == XAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == dxy   && axis == ZAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == dxy   && axis == XAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dxy   && axis == YAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dyz   && axis == YAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == dyz   && axis == ZAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == dyz   && axis == XAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == dyz   && axis == ZAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == dyz   && axis == XAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dyz   && axis == YAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dzx   && axis == YAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == dzx   && axis == ZAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == dzx   && axis == XAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == dzx   && axis == ZAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == dzx   && axis == XAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dzx   && axis == YAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dxxyy && axis == YAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == dxxyy && axis == ZAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == dxxyy && axis == XAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == dxxyy && axis == ZAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == dxxyy && axis == XAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dxxyy && axis == YAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dzz   && axis == YAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == dzz   && axis == ZAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == dzz   && axis == XAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == dzz   && axis == ZAxis && valenceOrbitalB == py) || 
            (valenceOrbitalA == dzz   && axis == XAxis && valenceOrbitalB == pz) || 
            (valenceOrbitalA == dzz   && axis == YAxis && valenceOrbitalB == pz) ){
      OrbitalType dOrbital;
      if( (valenceOrbitalB == py && axis == XAxis) || 
          (valenceOrbitalB == px && axis == YAxis) ){
         dOrbital = dxy;
      }
      else if( (valenceOrbitalB == py && axis == ZAxis) || 
               (valenceOrbitalB == pz && axis == YAxis) ){
         dOrbital = dyz;
      }
      else if( (valenceOrbitalB == px && axis == ZAxis) || 
               (valenceOrbitalB == pz && axis == XAxis) ){
         dOrbital = dzx;
      }
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      double overlapAOs2 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       dOrbital, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = overlapAOs2/(2.0*sqrt(gaussianExponentB))+xyzB[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == px && axis == XAxis && valenceOrbitalB == dxxyy){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5-2.0*gaussianExponentA*gaussianExponentB*pow(dxyz[XAxis],2.0)/beta;
      value += 0.5*pow(gaussianExponentA,2.0)*(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/beta;
      value += pow(gaussianExponentA*gaussianExponentB*dxyz[XAxis]/beta,2.0)
              *(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0));
      value *= 4.0*pow(gaussianExponentA,0.5)*gaussianExponentB*pow(beta,-2.0);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == XAxis && valenceOrbitalB == px){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5-2.0*gaussianExponentA*gaussianExponentB*pow(dxyz[XAxis],2.0)/beta;
      value += 0.5*pow(gaussianExponentB,2.0)*(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/beta;
      value += pow(gaussianExponentA*gaussianExponentB*dxyz[XAxis]/beta,2.0)
              *(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0));
      value *= 4.0*pow(gaussianExponentB,0.5)*gaussianExponentA*pow(beta,-2.0);
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == py && axis == YAxis && valenceOrbitalB == dxxyy){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5-2.0*gaussianExponentA*gaussianExponentB*pow(dxyz[YAxis],2.0)/beta;
      value += 0.5*pow(gaussianExponentA,2.0)*(pow(dxyz[YAxis],2.0)-pow(dxyz[XAxis],2.0))/beta;
      value += pow(gaussianExponentA*gaussianExponentB*dxyz[YAxis]/beta,2.0)
              *(pow(dxyz[YAxis],2.0)-pow(dxyz[XAxis],2.0));
      value *= -4.0*pow(gaussianExponentA,0.5)*gaussianExponentB*pow(beta,-2.0);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == YAxis && valenceOrbitalB == py){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5-2.0*gaussianExponentA*gaussianExponentB*pow(dxyz[YAxis],2.0)/beta;
      value += 0.5*pow(gaussianExponentB,2.0)*(pow(dxyz[YAxis],2.0)-pow(dxyz[XAxis],2.0))/beta;
      value += pow(gaussianExponentA*gaussianExponentB*dxyz[YAxis]/beta,2.0)
              *(pow(dxyz[YAxis],2.0)-pow(dxyz[XAxis],2.0));
      value *= -4.0*pow(gaussianExponentB,0.5)*gaussianExponentA*pow(beta,-2.0);
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == pz && axis == ZAxis && valenceOrbitalB == dxxyy){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0));
      value += pow(gaussianExponentB*dxyz[ZAxis],2.0)
              *(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))
              /beta;
      value *= 4.0*pow(gaussianExponentA,2.5)*gaussianExponentB*pow(beta,-3.0);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == ZAxis && valenceOrbitalB == pz){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0));
      value += pow(gaussianExponentA*dxyz[ZAxis],2.0)
              *(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))
              /beta;
      value *= 4.0*pow(gaussianExponentB,2.5)*gaussianExponentA*pow(beta,-3.0);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == px && axis == XAxis && valenceOrbitalB == dzz) || 
            (valenceOrbitalA == py && axis == YAxis && valenceOrbitalB == dzz) ){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5
             +2.0*gaussianExponentA*gaussianExponentB*pow(dxyz[axis],2.0)/beta
             +0.5*pow(gaussianExponentA,2.0)
                 *(2.0*pow(dxyz[ZAxis],2.0)-pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/beta
             +pow(gaussianExponentA*gaussianExponentB*dxyz[axis]/beta,2.0)
              *(2.0*pow(dxyz[ZAxis],2.0)-pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0));
      value *= 4.0*pow(gaussianExponentA,0.5)*gaussianExponentB*pow(beta,-2.0)/sqrt(3.0);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dzz && axis == XAxis && valenceOrbitalB == px) || 
            (valenceOrbitalA == dzz && axis == YAxis && valenceOrbitalB == py) ){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5
             +2.0*gaussianExponentA*gaussianExponentB*pow(dxyz[axis],2.0)/beta
             +0.5*pow(gaussianExponentB,2.0)
                 *(2.0*pow(dxyz[ZAxis],2.0)-pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/beta
             +pow(gaussianExponentA*gaussianExponentB*dxyz[axis]/beta,2.0)
              *(2.0*pow(dxyz[ZAxis],2.0)-pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0));
      value *= 4.0*pow(gaussianExponentB,0.5)*gaussianExponentA*pow(beta,-2.0)/sqrt(3.0);
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == pz && axis == ZAxis && valenceOrbitalB == dzz){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 1.0
             -4.0*gaussianExponentA*gaussianExponentB*pow(dxyz[axis],2.0)/beta
             +0.5*pow(gaussianExponentA,2.0)
                 *(2.0*pow(dxyz[ZAxis],2.0)-pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/beta
             +pow(gaussianExponentA*gaussianExponentB*dxyz[axis]/beta,2.0)
              *(2.0*pow(dxyz[ZAxis],2.0)-pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0));
      value *= 4.0*pow(gaussianExponentA,0.5)*gaussianExponentB*pow(beta,-2.0)/sqrt(3.0);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dzz && axis == ZAxis && valenceOrbitalB == pz){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 1.0
             -4.0*gaussianExponentA*gaussianExponentB*pow(dxyz[axis],2.0)/beta
             +0.5*pow(gaussianExponentB,2.0)
                 *(2.0*pow(dxyz[ZAxis],2.0)-pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/beta
             +pow(gaussianExponentA*gaussianExponentB*dxyz[axis]/beta,2.0)
              *(2.0*pow(dxyz[ZAxis],2.0)-pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0));
      value *= 4.0*pow(gaussianExponentB,0.5)*gaussianExponentA*pow(beta,-2.0)/sqrt(3.0);
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy && axis == XAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dxy && axis == YAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dyz && axis == YAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dyz && axis == ZAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dzx && axis == ZAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dzx && axis == XAxis && valenceOrbitalB == valenceOrbitalA) ){
      CartesianType anotherAxis;
      if(valenceOrbitalB == dxy){
         if(axis == XAxis){
            anotherAxis = YAxis;
         }
         else{
            anotherAxis = XAxis;
         }
      }
      else if(valenceOrbitalB == dyz){
         if(axis == YAxis){
            anotherAxis = ZAxis;
         }
         else{
            anotherAxis = YAxis;
         }
      }
      else if(valenceOrbitalB == dzx){
         if(axis == ZAxis){
            anotherAxis = XAxis;
         }
         else{
            anotherAxis = ZAxis;
         }
      }
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*(gaussianExponentA-gaussianExponentB)*dxyz[axis]
             +(gaussianExponentA-gaussianExponentB)*(gaussianExponentA*gaussianExponentB)
              *dxyz[axis]*pow(dxyz[anotherAxis],2.0)/beta;
      value *= 8.0*(gaussianExponentA*gaussianExponentB)*pow(beta,-3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy   && axis == ZAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dyz   && axis == XAxis && valenceOrbitalB == valenceOrbitalA) || 
            (valenceOrbitalA == dzx   && axis == YAxis && valenceOrbitalB == valenceOrbitalA) ||
            (valenceOrbitalA == dxxyy && axis == ZAxis && valenceOrbitalB == valenceOrbitalA) ){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxxyy && axis == XAxis && valenceOrbitalB == valenceOrbitalA) ||
            (valenceOrbitalA == dxxyy && axis == YAxis && valenceOrbitalB == valenceOrbitalA) ){
      CartesianType anotherAxis;
      if(axis == XAxis){
         anotherAxis = YAxis;
      }
      else{
         anotherAxis = XAxis;
      }
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = (gaussianExponentA-gaussianExponentB)*dxyz[axis]
             -(gaussianExponentA-gaussianExponentB)*(gaussianExponentA*gaussianExponentB)
              *(pow(dxyz[axis],2.0) - pow(dxyz[anotherAxis],2.0))*dxyz[axis]/beta;
      value *= 4.0*(gaussianExponentA*gaussianExponentB)*pow(beta,-3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dzz && axis == XAxis && valenceOrbitalB == valenceOrbitalA) ||
            (valenceOrbitalA == dzz && axis == YAxis && valenceOrbitalB == valenceOrbitalA) ){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = (gaussianExponentA-gaussianExponentB)*dxyz[axis]
             -(gaussianExponentA-gaussianExponentB)*(gaussianExponentA*gaussianExponentB)
              *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))*dxyz[axis]/beta;
      value *= 4.0*(gaussianExponentA*gaussianExponentB)*pow(beta,-3.0)/3.0;
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( valenceOrbitalA == dzz && axis == ZAxis && valenceOrbitalB == valenceOrbitalA) {
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 2.0*(gaussianExponentA-gaussianExponentB)*dxyz[axis]
             -(gaussianExponentA-gaussianExponentB)*(gaussianExponentA*gaussianExponentB)
              *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))*dxyz[axis]/beta;
      value *= 8.0*(gaussianExponentA*gaussianExponentB)*pow(beta,-3.0)/3.0;
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy && axis == YAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dyz && axis == ZAxis && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == dzx && axis == XAxis && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == dyz && axis == YAxis && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == dzx && axis == ZAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dxy && axis == XAxis && valenceOrbitalB == dzx) ){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = -8.0*pow(gaussianExponentA*gaussianExponentB,2.0)*dxyz[XAxis]*dxyz[YAxis]*dxyz[ZAxis]
             *(gaussianExponentA-gaussianExponentB)*pow(beta,-4.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy && axis == XAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dyz && axis == YAxis && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == dzx && axis == ZAxis && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == dyz && axis == ZAxis && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == dzx && axis == XAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dxy && axis == YAxis && valenceOrbitalB == dzx) ){
      CartesianType anotherAxis1;
      CartesianType anotherAxis2;
      if(valenceOrbitalA == dxy && valenceOrbitalB == dyz){
         anotherAxis1 = YAxis;
         anotherAxis2 = ZAxis;
      }
      else if(valenceOrbitalA == dyz && valenceOrbitalB == dzx){
         anotherAxis1 = ZAxis;
         anotherAxis2 = XAxis;
      }
      else if(valenceOrbitalA == dzx && valenceOrbitalB == dxy){
         anotherAxis1 = XAxis;
         anotherAxis2 = YAxis;
      }
      else if(valenceOrbitalA == dyz && valenceOrbitalB == dxy){
         anotherAxis1 = YAxis;
         anotherAxis2 = XAxis;
      }
      else if(valenceOrbitalA == dzx && valenceOrbitalB == dyz){
         anotherAxis1 = ZAxis;
         anotherAxis2 = YAxis;
      }
      else if(valenceOrbitalA == dxy && valenceOrbitalB == dzx){
         anotherAxis1 = XAxis;
         anotherAxis2 = ZAxis;
      }
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5-gaussianExponentA*gaussianExponentB*pow(dxyz[anotherAxis1],2.0)/beta;
      value *= 8.0*pow(gaussianExponentA,2.0)*gaussianExponentB*pow(beta,-3.0)*dxyz[anotherAxis2];
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dyz && axis == XAxis && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == dzx && axis == YAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dxy && axis == ZAxis && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == dxy && axis == ZAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dyz && axis == XAxis && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == dzx && axis == YAxis && valenceOrbitalB == dxy) ){
      CartesianType anotherAxis1;
      CartesianType anotherAxis2;
      if(valenceOrbitalA == dyz && valenceOrbitalB == dxy){
         anotherAxis1 = YAxis;
         anotherAxis2 = ZAxis;
      }
      else if(valenceOrbitalA == dzx && valenceOrbitalB == dyz){
         anotherAxis1 = ZAxis;
         anotherAxis2 = XAxis;
      }
      else if(valenceOrbitalA == dxy && valenceOrbitalB == dzx){
         anotherAxis1 = XAxis;
         anotherAxis2 = YAxis;
      }
      else if(valenceOrbitalA == dxy && valenceOrbitalB == dyz){
         anotherAxis1 = YAxis;
         anotherAxis2 = XAxis;
      }
      else if(valenceOrbitalA == dyz && valenceOrbitalB == dzx){
         anotherAxis1 = ZAxis;
         anotherAxis2 = YAxis;
      }
      else if(valenceOrbitalA == dzx && valenceOrbitalB == dxy){
         anotherAxis1 = XAxis;
         anotherAxis2 = ZAxis;
      }

      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5-gaussianExponentA*gaussianExponentB*pow(dxyz[anotherAxis1],2.0)/beta;
      value *= -8.0*pow(gaussianExponentB,2.0)*gaussianExponentA*pow(beta,-3.0)*dxyz[anotherAxis2];
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == XAxis && valenceOrbitalB == dxy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*beta
             -pow(gaussianExponentA,2.0)*gaussianExponentB*pow(dxyz[XAxis],2.0)/beta
             +pow(gaussianExponentB,2.0)*gaussianExponentA
             *(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/(2.0*beta);
      value *= 8.0*gaussianExponentA*gaussianExponentB*pow(beta,-3.0)*dxyz[YAxis];
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxy && axis == XAxis && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*beta
             -pow(gaussianExponentB,2.0)*gaussianExponentA*pow(dxyz[XAxis],2.0)/beta
             +pow(gaussianExponentA,2.0)*gaussianExponentB
             *(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/(2.0*beta);
      value *= -8.0*gaussianExponentA*gaussianExponentB*pow(beta,-3.0)*dxyz[YAxis];
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   } 
   else if(valenceOrbitalA == dxxyy && axis == YAxis && valenceOrbitalB == dxy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5*beta
             +pow(gaussianExponentA,2.0)*gaussianExponentB*pow(dxyz[YAxis],2.0)/beta
             +pow(gaussianExponentB,2.0)*gaussianExponentA
             *(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/(2.0*beta);
      value *= 8.0*gaussianExponentA*gaussianExponentB*pow(beta,-3.0)*dxyz[XAxis];
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxy && axis == YAxis && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5*beta
             +pow(gaussianExponentB,2.0)*gaussianExponentA*pow(dxyz[YAxis],2.0)/beta
             +pow(gaussianExponentA,2.0)*gaussianExponentB
             *(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/(2.0*beta);
      value *= -8.0*gaussianExponentA*gaussianExponentB*pow(beta,-3.0)*dxyz[XAxis];
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxxyy && axis == XAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dxxyy && axis == YAxis && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == dxxyy && axis == ZAxis && valenceOrbitalB == dxy) ){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 8.0*pow(gaussianExponentA,4.0)*pow(gaussianExponentB,3.0)
             *dxyz[XAxis]*dxyz[YAxis]*dxyz[ZAxis]
             *(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/pow(beta,5.0);
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dyz && axis == XAxis && valenceOrbitalB == dxxyy) ||
            (valenceOrbitalA == dzx && axis == YAxis && valenceOrbitalB == dxxyy) ||
            (valenceOrbitalA == dxy && axis == ZAxis && valenceOrbitalB == dxxyy) ){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = -8.0*pow(gaussianExponentA,3.0)*pow(gaussianExponentB,4.0)
             *dxyz[XAxis]*dxyz[YAxis]*dxyz[ZAxis]
             *(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/pow(beta,5.0);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == YAxis && valenceOrbitalB == dyz){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5*gaussianExponentA
             +pow(gaussianExponentA,2.0)*gaussianExponentB*pow(dxyz[YAxis],2.0)/beta
             +pow(gaussianExponentB,2.0)*gaussianExponentA
             *(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/(2.0*beta);
      value *= 8.0*gaussianExponentA*gaussianExponentB*dxyz[ZAxis]*pow(beta,-3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dyz && axis == YAxis && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5*gaussianExponentB
             +pow(gaussianExponentB,2.0)*gaussianExponentA*pow(dxyz[YAxis],2.0)/beta
             +pow(gaussianExponentA,2.0)*gaussianExponentB
             *(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/(2.0*beta);
      value *= -8.0*gaussianExponentA*gaussianExponentB*dxyz[ZAxis]*pow(beta,-3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == XAxis && valenceOrbitalB == dzx){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5*gaussianExponentA
             +pow(gaussianExponentA,2.0)*gaussianExponentB*pow(dxyz[XAxis],2.0)/beta
             +pow(gaussianExponentB,2.0)*gaussianExponentA
             *(pow(dxyz[YAxis],2.0)-pow(dxyz[XAxis],2.0))/(2.0*beta);
      value *= -8.0*gaussianExponentA*gaussianExponentB*dxyz[ZAxis]*pow(beta,-3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dzx && axis == XAxis && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = -0.5*gaussianExponentB
             +pow(gaussianExponentB,2.0)*gaussianExponentA*pow(dxyz[YAxis],2.0)/beta
             +pow(gaussianExponentA,2.0)*gaussianExponentB
             *(pow(dxyz[YAxis],2.0)-pow(dxyz[XAxis],2.0))/(2.0*beta);
      value *= 8.0*gaussianExponentA*gaussianExponentB*dxyz[ZAxis]*pow(beta,-3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == ZAxis && valenceOrbitalB == dyz){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentA*gaussianExponentB*(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/beta+1.0;
      value *= 4.0*gaussianExponentA*pow(gaussianExponentB,2.0)*dxyz[YAxis]*pow(beta,-3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dyz && axis == ZAxis && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentA*gaussianExponentB*(pow(dxyz[XAxis],2.0)-pow(dxyz[YAxis],2.0))/beta+1.0;
      value *= -4.0*gaussianExponentB*pow(gaussianExponentA,2.0)*dxyz[YAxis]*pow(beta,-3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == ZAxis && valenceOrbitalB == dzx){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentA*gaussianExponentB*(pow(dxyz[YAxis],2.0)-pow(dxyz[XAxis],2.0))/beta+1.0;
      value *= -4.0*gaussianExponentA*pow(gaussianExponentB,2.0)*dxyz[XAxis]*pow(beta,-3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dzx && axis == ZAxis && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentA*gaussianExponentB*(pow(dxyz[YAxis],2.0)-pow(dxyz[XAxis],2.0))/beta+1.0;
      value *= 4.0*gaussianExponentB*pow(gaussianExponentA,2.0)*dxyz[XAxis]*pow(beta,-3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dzz && axis == XAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dzz && axis == YAxis && valenceOrbitalB == dzx) ||
            (valenceOrbitalA == dzz && axis == ZAxis && valenceOrbitalB == dxy) ){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0);
      value *= dxyz[XAxis]*dxyz[YAxis]*dxyz[ZAxis];
      value *= 8.0*pow(gaussianExponentA,4.0)*pow(gaussianExponentB,3.0);
      value /= sqrt(3.0)*pow(beta,5.0);
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dyz && axis == XAxis && valenceOrbitalB == dzz) ||
            (valenceOrbitalA == dzx && axis == YAxis && valenceOrbitalB == dzz) ||
            (valenceOrbitalA == dxy && axis == ZAxis && valenceOrbitalB == dzz) ){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0);
      value *= dxyz[XAxis]*dxyz[YAxis]*dxyz[ZAxis];
      value *= -8.0*pow(gaussianExponentB,4.0)*pow(gaussianExponentA,3.0);
      value /= sqrt(3.0)*pow(beta,5.0);
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dzz && axis == XAxis && valenceOrbitalB == dxy) ||
            (valenceOrbitalA == dzz && axis == YAxis && valenceOrbitalB == dxy) ){
      CartesianType anotherAxis;
      if(axis == XAxis){
         anotherAxis = YAxis;
      }
      else if(axis == YAxis){
         anotherAxis = XAxis;
      }
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*(gaussianExponentB-gaussianExponentA)
             +3.0*pow(gaussianExponentA,2.0)*gaussianExponentB*pow(dxyz[axis],2.0)*pow(beta,-2.0)
             +gaussianExponentA*pow(gaussianExponentB,2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))/(2.0*beta)
             +pow(gaussianExponentA,3.0)*pow(gaussianExponentB,2.0)*pow(dxyz[axis],2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))*pow(beta,-2.0);
      value *= 8.0*gaussianExponentA*gaussianExponentB*dxyz[anotherAxis]/(sqrt(3.0)*pow(beta,3.0));
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dxy && axis == XAxis && valenceOrbitalB == dzz) ||
            (valenceOrbitalA == dxy && axis == YAxis && valenceOrbitalB == dzz) ){
      CartesianType anotherAxis;
      if(axis == XAxis){
         anotherAxis = YAxis;
      }
      else if(axis == YAxis){
         anotherAxis = XAxis;
      }
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = 0.5*(gaussianExponentA-gaussianExponentB)
             +3.0*pow(gaussianExponentB,2.0)*gaussianExponentA*pow(dxyz[axis],2.0)*pow(beta,-2.0)
             +gaussianExponentB*pow(gaussianExponentA,2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))/(2.0*beta)
             +pow(gaussianExponentB,3.0)*pow(gaussianExponentA,2.0)*pow(dxyz[axis],2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))*pow(beta,-2.0);
      value *= -8.0*gaussianExponentA*gaussianExponentB*dxyz[anotherAxis]/(sqrt(3.0)*pow(beta,3.0));
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dzz && axis == YAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dzz && axis == XAxis && valenceOrbitalB == dzx) ){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentB-0.5*gaussianExponentA
             +gaussianExponentA*pow(gaussianExponentB,2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))/(2.0*beta)
             +pow(gaussianExponentA,3.0)*pow(gaussianExponentB*dxyz[axis],2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))
             *pow(beta,-2.0);
      value *= 8.0*gaussianExponentA*gaussianExponentB*dxyz[ZAxis]
              /(sqrt(3.0)*pow(beta,3.0));
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dyz && axis == YAxis && valenceOrbitalB == dzz) ||
            (valenceOrbitalA == dzx && axis == XAxis && valenceOrbitalB == dzz) ){
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentA-0.5*gaussianExponentB
             +gaussianExponentB*pow(gaussianExponentA,2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))/(2.0*beta)
             +pow(gaussianExponentB,3.0)*pow(gaussianExponentA*dxyz[axis],2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))
             *pow(beta,-2.0);
      value *= -8.0*gaussianExponentA*gaussianExponentB*dxyz[ZAxis]
              /(sqrt(3.0)*pow(beta,3.0));
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dzz && axis == ZAxis && valenceOrbitalB == dyz) ||
            (valenceOrbitalA == dzz && axis == ZAxis && valenceOrbitalB == dzx) ){
      CartesianType anotherAxis;
      if(valenceOrbitalB == dyz){
         anotherAxis = YAxis;
      }
      else if(valenceOrbitalB == dzx){
         anotherAxis = XAxis;
      }
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentA-0.5*gaussianExponentB
             -3.0*pow(gaussianExponentA,2.0)*gaussianExponentB*pow(dxyz[axis]/beta,2.0)
             +gaussianExponentA*pow(gaussianExponentB,2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))/(2.0*beta)
             +pow(gaussianExponentA,3.0)*pow(gaussianExponentB*dxyz[axis],2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))
             *pow(beta,-2.0);
      value *= 8.0*gaussianExponentA*gaussianExponentB*dxyz[anotherAxis]
              /(sqrt(3.0)*pow(beta,3.0));
      value *= overlapSASB;
      value += xyzB[axis]*overlapAOs1;
      return value;
   }
   else if( (valenceOrbitalA == dyz && axis == ZAxis && valenceOrbitalB == dzz) ||
            (valenceOrbitalA == dzx && axis == ZAxis && valenceOrbitalB == dzz) ){
      CartesianType anotherAxis;
      if(valenceOrbitalA == dyz){
         anotherAxis = YAxis;
      }
      else if(valenceOrbitalA == dzx){
         anotherAxis = XAxis;
      }
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentB-0.5*gaussianExponentA
             -3.0*pow(gaussianExponentB,2.0)*gaussianExponentA*pow(dxyz[axis]/beta,2.0)
             +gaussianExponentB*pow(gaussianExponentA,2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))/(2.0*beta)
             +pow(gaussianExponentB,3.0)*pow(gaussianExponentA*dxyz[axis],2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))
             *pow(beta,-2.0);
      value *= -8.0*gaussianExponentA*gaussianExponentB*dxyz[anotherAxis]
              /(sqrt(3.0)*pow(beta,3.0));
      value *= overlapSASB;
      value += xyzA[axis]*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dzz && axis == XAxis && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentB - gaussianExponentA
             +gaussianExponentA*pow(gaussianExponentB,2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))/beta
             +pow(gaussianExponentA,2.0)*gaussianExponentB
             *(pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))/beta;
      value *= 4.0*gaussianExponentA*gaussianExponentB*dxyz[XAxis]*pow(beta,-3.0)/sqrt(3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == XAxis && valenceOrbitalB == dzz){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentA - gaussianExponentB
             +gaussianExponentB*pow(gaussianExponentA,2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))/beta
             +pow(gaussianExponentB,2.0)*gaussianExponentA
             *(pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))/beta;
      value *= -4.0*gaussianExponentA*gaussianExponentB*dxyz[XAxis]*pow(beta,-3.0)/sqrt(3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dzz && axis == YAxis && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentB - gaussianExponentA
             +gaussianExponentA*pow(gaussianExponentB,2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))/beta
             +pow(gaussianExponentA,2.0)*gaussianExponentB
             *(pow(dxyz[YAxis],2.0) - pow(dxyz[XAxis],2.0))/beta;
      value *= -4.0*gaussianExponentA*gaussianExponentB*dxyz[YAxis]*pow(beta,-3.0)/sqrt(3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == YAxis && valenceOrbitalB == dzz){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = gaussianExponentA - gaussianExponentB
             +gaussianExponentB*pow(gaussianExponentA,2.0)
             *(2.0*pow(dxyz[ZAxis],2.0) - pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0))/beta
             +pow(gaussianExponentB,2.0)*gaussianExponentA
             *(pow(dxyz[YAxis],2.0) - pow(dxyz[XAxis],2.0))/beta;
      value *= 4.0*gaussianExponentA*gaussianExponentB*dxyz[YAxis]*pow(beta,-3.0)/sqrt(3.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dzz && axis == ZAxis && valenceOrbitalB == dxxyy){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0);
      value *= -8.0*pow(gaussianExponentA,3.0)*pow(gaussianExponentB,2.0)*dxyz[ZAxis];
      value /= sqrt(3.0)*pow(beta,4.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else if(valenceOrbitalA == dxxyy && axis == ZAxis && valenceOrbitalB == dzz){
      double axisAverage = (gaussianExponentA*xyzA[axis]+gaussianExponentB*xyzB[axis])/beta;
      double overlapAOs1 = this->GetGaussianOverlapAOs(atomTypeA,
                                                       valenceOrbitalA, 
                                                       gaussianExponentA, 
                                                       atomTypeB, 
                                                       valenceOrbitalB, 
                                                       gaussianExponentB,
                                                       dxyz[XAxis], 
                                                       dxyz[YAxis], 
                                                       dxyz[ZAxis], 
                                                       rAB,
                                                       overlapSASB);
      value = pow(dxyz[XAxis],2.0) - pow(dxyz[YAxis],2.0);
      value *= 8.0*pow(gaussianExponentA,2.0)*pow(gaussianExponentB,3.0)*dxyz[ZAxis];
      value /= sqrt(3.0)*pow(beta,4.0);
      value *= overlapSASB;
      value += axisAverage*overlapAOs1;
      return value;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetGaussianCartesianMatrixBadOrbital;
      ss << this->errorMessageAtomA;
      ss << this->errorMessageAtomType << AtomTypeStr(atomTypeA) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(valenceOrbitalA) << endl;
      ss << this->errorMessageAtomB;
      ss << this->errorMessageAtomType << AtomTypeStr(atomTypeB) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(valenceOrbitalB) << endl;
      ss << this->errorMessageCartesianType << CartesianTypeStr(axis) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}

void Cndo2::FreeDiatomicOverlapAOsAndRotatingMatrix(double*** diatomicOverlapAOs, 
                                                    double*** rotatingMatrix) const{
   // free
   MallocerFreer::GetInstance()->Free<double>(diatomicOverlapAOs, OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Free<double>(rotatingMatrix,  OrbitalType_end, OrbitalType_end);
}

// calculate OverlapAOs matrix between different configurations, S^{AO}_{\mu\nu}.
// \mu and \nu are AOs belonging to left and right hand side configurations, respectively.
// S^{AO}_{\mu\nu} = 0 when \mu and \nu belong different atom.
// Note that rhs-moledule is this->molecule (current configuration) 
// and lhs-molecule is another molecule (another configuration).
void Cndo2::CalcOverlapAOsWithAnotherConfiguration(double** overlapAOs, 
                                                   const Molecule& lhsMolecule) const{
   const Molecule* rhsMolecule = this->molecule;
   if(lhsMolecule.GetTotalNumberAOs() != rhsMolecule->GetTotalNumberAOs()){
      stringstream ss;
      ss << this->errorMessageCalcOverlapAOsDifferentConfigurationsDiffAOs;
      ss << this->errorMessageLhs << lhsMolecule.GetTotalNumberAOs() << endl;
      ss << this->errorMessageRhs << rhsMolecule->GetTotalNumberAOs() << endl;
      throw MolDSException(ss.str());
   }
   if(lhsMolecule.GetNumberAtoms() != rhsMolecule->GetNumberAtoms()){
      stringstream ss;
      ss << this->errorMessageCalcOverlapAOsDifferentConfigurationsDiffAtoms;
      ss << this->errorMessageLhs << lhsMolecule.GetNumberAtoms() << endl;
      ss << this->errorMessageRhs << rhsMolecule->GetNumberAtoms() << endl;
      throw MolDSException(ss.str());
   }
#ifdef MOLDS_DBG
   if(overlapAOs == NULL){
      throw MolDSException(this->errorMessageCalcOverlapAOsDifferentConfigurationsOverlapAOsNULL);
   }
#endif
   
   int totalAONumber = lhsMolecule.GetTotalNumberAOs();
   int totalAtomNumber = lhsMolecule.GetNumberAtoms();
   MallocerFreer::GetInstance()->Initialize<double>(overlapAOs, totalAONumber, totalAONumber);

   stringstream ompErrors;
#pragma omp parallel 
   {
      double** diatomicOverlapAOs = NULL;
      double** rotatingMatrix = NULL;
      try{
         // malloc
         MallocerFreer::GetInstance()->Malloc<double>(&diatomicOverlapAOs,
                                                      OrbitalType_end, 
                                                      OrbitalType_end);
         MallocerFreer::GetInstance()->Malloc<double>(&rotatingMatrix,
                                                      OrbitalType_end, 
                                                      OrbitalType_end);
         // calculation overlapAOs matrix
         for(int mu=0; mu<totalAONumber; mu++){
            overlapAOs[mu][mu] = 1.0;
         }
         bool isSymmetricOverlapAOs = false;
#pragma omp for schedule(auto)
         for(int A=0; A<totalAtomNumber; A++){
            const Atom& lhsAtom = *lhsMolecule.GetAtom(A);
            const Atom& rhsAtom = *rhsMolecule->GetAtom(A);
            int firstAOIndexLhsAtom = lhsAtom.GetFirstAOIndex();
            int firstAOIndexRhsAtom = rhsAtom.GetFirstAOIndex();
            for(int i=0; i<lhsAtom.GetValenceSize(); i++){
               for(int j=0; j<rhsAtom.GetValenceSize(); j++){
                  int mu = firstAOIndexLhsAtom + i;      
                  int nu = firstAOIndexRhsAtom + j;      
                  double value = this->GetOverlapAOsElementByGTOExpansion(lhsAtom, i, rhsAtom, j, STO6G);
                  overlapAOs[mu][nu] = value;
               }
            }
         }
      }
      catch(MolDSException ex){
#pragma omp critical
         ex.Serialize(ompErrors);
      }
      this->FreeDiatomicOverlapAOsAndRotatingMatrix(&diatomicOverlapAOs, &rotatingMatrix);
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
}

// calculate OverlapMOs matrix between different electronic-structure, S^{MO}_{ij}.
// i and j are MOs belonging to left and right hand side electronic-structures, respectively.
// Note that rhs-electronic-structure is this electronic-structure  
// and lhs-electronic-structure is another electronic-structure.
void Cndo2::CalcOverlapMOsWithAnotherElectronicStructure(double** overlapMOs, 
                                                         double const* const* overlapAOs,
                                                         const ElectronicStructure& lhsElectronicStructure) const{
   const ElectronicStructure* rhsElectronicStructure = this;
   double const* const* rhsFockMatrix = this->fockMatrix;
   double const* const* lhsFockMatrix = lhsElectronicStructure.GetFockMatrix();
   int totalAONumber = this->molecule->GetTotalNumberAOs();
   int usedMONumber = this->molecule->GetTotalNumberValenceElectrons()/2
                     +Parameters::GetInstance()->GetActiveVirCIS();
   MallocerFreer::GetInstance()->Initialize<double>(overlapMOs, totalAONumber, totalAONumber);
   double** tmpMatrix=NULL;
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&tmpMatrix,totalAONumber,totalAONumber);
      bool isColumnMajorOverlapAOs = false;
      bool isColumnMajorRhsFock = true;
      double alpha=1.0;
      double beta=0.0;
      MolDS_wrappers::Blas::GetInstance()->Dgemm(isColumnMajorOverlapAOs,
                                                 isColumnMajorRhsFock,
                                                 totalAONumber,usedMONumber,totalAONumber,
                                                 alpha,
                                                 overlapAOs,
                                                 rhsFockMatrix,
                                                 beta,
                                                 tmpMatrix);
      MolDS_wrappers::Blas::GetInstance()->Dgemm(usedMONumber,totalAONumber,totalAONumber,
                                                 lhsFockMatrix,
                                                 tmpMatrix,
                                                 overlapMOs);
                                                 
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&tmpMatrix,totalAONumber,totalAONumber);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&tmpMatrix,totalAONumber,totalAONumber);
}

// calculate OverlapSingletSDs matrix between different electronic-structure, S^{SSD}_{ij}.
// i and j are singlet SDs belonging to left and right hand side electronic-structures, respectively.
// The index i=0 means the Hartree-Fock state.
// This overlapsingletSDs are calculated from overlapMOs.
// Note that rhs-electronic-structure is this electronic-structure  
// and lhs-electronic-structure is another electronic-structure.
void Cndo2::CalcOverlapSingletSDsWithAnotherElectronicStructure(double** overlapSingletSDs, 
                                                                double const* const* overlapMOs) const{
      stringstream ss;
      ss << this->errorMessageNonExcitedStates;
      throw MolDSException(ss.str());
}

// calculate overlapESs (ES means eigenstate) matrix between different electronic-structure, S^{ES}_{ij}.
// i and j are singlet SDs belonging to left and right hand side electronic-structures, respectively.
// The index i=0 means the ground state.
// This overlapESs is calculated from the overlapsingletSDs.
// Note that rhs-electronic-structure is this electronic-structure  
// and lhs-electronic-structure is another electronic-structure.
void Cndo2::CalcOverlapESsWithAnotherElectronicStructure(double** overlapESs, 
                                                         double const* const* overlapSingletSDs,
                                                         const MolDS_base::ElectronicStructure& lhsElectronicStructure) const{
      stringstream ss;
      ss << this->errorMessageNonExcitedStates;
      throw MolDSException(ss.str());
}

// calculate OverlapAOs matrix. E.g. S_{\mu\nu} in (3.74) in J. A. Pople book.
void Cndo2::CalcOverlapAOs(double** overlapAOs, const Molecule& molecule) const{
   int totalAONumber = molecule.GetTotalNumberAOs();
   int totalAtomNumber = molecule.GetNumberAtoms();

   stringstream ompErrors;
#pragma omp parallel 
   {
      double** diatomicOverlapAOs = NULL;
      double** rotatingMatrix = NULL;
      try{
         // malloc
         MallocerFreer::GetInstance()->Malloc<double>(&diatomicOverlapAOs,
                                                      OrbitalType_end, 
                                                      OrbitalType_end);
         MallocerFreer::GetInstance()->Malloc<double>(&rotatingMatrix,
                                                      OrbitalType_end, 
                                                      OrbitalType_end);
         // calculation overlapAOs matrix
         for(int mu=0; mu<totalAONumber; mu++){
            overlapAOs[mu][mu] = 1.0;
         }

#pragma omp for schedule(auto)
         for(int A=0; A<totalAtomNumber; A++){
            const Atom& atomA = *molecule.GetAtom(A);
            for(int B=A+1; B<totalAtomNumber; B++){
               const Atom& atomB = *molecule.GetAtom(B);
               this->CalcDiatomicOverlapAOsInDiatomicFrame(diatomicOverlapAOs, atomA, atomB);
               this->CalcRotatingMatrix(rotatingMatrix, atomA, atomB);
               this->RotateDiatmicOverlapAOsToSpaceFrame(diatomicOverlapAOs, rotatingMatrix);
               this->SetOverlapAOsElement(overlapAOs, diatomicOverlapAOs, atomA, atomB);
            }
         }
      }
      catch(MolDSException ex){
#pragma omp critical
         ex.Serialize(ompErrors);
      }
      this->FreeDiatomicOverlapAOsAndRotatingMatrix(&diatomicOverlapAOs, &rotatingMatrix);
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
   /* 
   this->OutputLog("overlapAOs matrix\n"); 
   for(int o=0; o<molecule.GetTotalNumberAOs(); o++){
      for(int p=0; p<molecule.GetTotalNumberAOs(); p++){
         this->OutputLog(boost::format("%lf\t") % overlapAOs[o][p]);
      }
      this->OutputLog("\n");
   }
   this->OutputLog("\n");
   */
}

// First derivative of diatomic overlapAOs integrals between AOs in space fixed flame.
// The OverlapAOs matrix is S_{\mu\nu} in (3.74) in J. A. Pople book.
// Note that this method can not treat d-obitals 
// because CalcRotatingMatrix1stDerivatives can not treat d-orbitals.
void Cndo2::CalcDiatomicOverlapAOs1stDerivatives(double*** diatomicOverlapAOs1stDerivs, 
                                                 const Atom& atomA, 
                                                 const Atom& atomB) const{
   double cartesian[CartesianType_end] = {atomA.GetXyz()[XAxis] - atomB.GetXyz()[XAxis], 
                                          atomA.GetXyz()[YAxis] - atomB.GetXyz()[YAxis],
                                          atomA.GetXyz()[ZAxis] - atomB.GetXyz()[ZAxis]};
   double R = sqrt( pow(cartesian[XAxis],2.0) + 
                    pow(cartesian[YAxis],2.0) + 
                    pow(cartesian[ZAxis],2.0) );
   
   double**  diaOverlapAOsInDiaFrame = NULL;  // diatomic overlapAOs in diatomic frame
   double**  diaOverlapAOs1stDerivInDiaFrame = NULL; // first derivative of the diaOverlapAOs. This derivative is related to the distance between two atoms.
   double**  rotMat = NULL; // rotating Matrix from the diatomic frame to space fixed frame.
   double*** rotMat1stDerivs = NULL; // first derivatives of the rotMat.
   double**  rotatedDiatomicOverlap = NULL;
   double**  tmpRotMat1stDeriv = NULL;
   double**  tmpMatrix         = NULL;

   try{
      this->MallocDiatomicOverlapAOs1stDeriTemps(&diaOverlapAOsInDiaFrame,
                                                 &diaOverlapAOs1stDerivInDiaFrame,
                                                 &rotMat,
                                                 &rotMat1stDerivs,
                                                 &rotatedDiatomicOverlap,
                                                 &tmpRotMat1stDeriv,
                                                 &tmpMatrix);
      this->CalcDiatomicOverlapAOsInDiatomicFrame(diaOverlapAOsInDiaFrame, atomA, atomB);
      this->CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(diaOverlapAOs1stDerivInDiaFrame, atomA, atomB);
      this->CalcRotatingMatrix(rotMat, atomA, atomB);
      this->CalcRotatingMatrix1stDerivatives(rotMat1stDerivs, atomA, atomB);

      // rotate (fast algorithm, see also slow algorithm shown later)
      int incrementOne = 1;
      bool isColumnMajorRotatingMatrix = false;
      bool isColumnMajorDiaOverlapAOs  = false;
      double alpha = 0.0;
      double beta  = 0.0;
      for(int c=0; c<CartesianType_end; c++){
         MolDS_wrappers::Blas::GetInstance()->Dcopy(OrbitalType_end*OrbitalType_end, 
                                                    &rotMat1stDerivs[0][0][c], CartesianType_end,
                                                    &tmpRotMat1stDeriv[0][0],  incrementOne);
         alpha = cartesian[c]/R;
         beta  = 0.0;
         MolDS_wrappers::Blas::GetInstance()->Dgemmm(isColumnMajorRotatingMatrix,
                                                     isColumnMajorDiaOverlapAOs,
                                                     !isColumnMajorRotatingMatrix,
                                                     OrbitalType_end, OrbitalType_end, OrbitalType_end, OrbitalType_end,
                                                     alpha,
                                                     rotMat,
                                                     diaOverlapAOs1stDerivInDiaFrame,
                                                     rotMat,
                                                     beta,
                                                     rotatedDiatomicOverlap,
                                                     tmpMatrix);
         alpha = 1.0;
         beta  = 1.0;
         MolDS_wrappers::Blas::GetInstance()->Dgemmm(isColumnMajorRotatingMatrix,
                                                     isColumnMajorDiaOverlapAOs,
                                                     !isColumnMajorRotatingMatrix,
                                                     OrbitalType_end, OrbitalType_end, OrbitalType_end, OrbitalType_end,
                                                     alpha,
                                                     tmpRotMat1stDeriv,
                                                     diaOverlapAOsInDiaFrame,
                                                     rotMat,
                                                     beta,
                                                     rotatedDiatomicOverlap,
                                                     tmpMatrix);
         MolDS_wrappers::Blas::GetInstance()->Dgemmm(isColumnMajorRotatingMatrix,
                                                     isColumnMajorDiaOverlapAOs,
                                                     !isColumnMajorRotatingMatrix,
                                                     OrbitalType_end, OrbitalType_end, OrbitalType_end, OrbitalType_end,
                                                     alpha,
                                                     rotMat,
                                                     diaOverlapAOsInDiaFrame,
                                                     tmpRotMat1stDeriv,
                                                     beta,
                                                     rotatedDiatomicOverlap,
                                                     tmpMatrix);
         MolDS_wrappers::Blas::GetInstance()->Dcopy(OrbitalType_end*OrbitalType_end, 
                                                    &rotatedDiatomicOverlap[0][0],                      incrementOne,
                                                    &diatomicOverlapAOs1stDerivs[0][0][c], CartesianType_end);
      }

      /* 
      // rotate (slow)
      for(int i=0; i<OrbitalType_end; i++){
         for(int j=0; j<OrbitalType_end; j++){
            for(int c=0; c<CartesianType_end; c++){
               diatomicOverlapAOs1stDerivs[i][j][c] = 0.0;

               double temp1 = 0.0;
               double temp2 = 0.0;
               double temp3 = 0.0;
               for(int k=0; k<OrbitalType_end; k++){
                  for(int l=0; l<OrbitalType_end; l++){
                     temp1 += rotMat[i][k] 
                             *rotMat[j][l]
                             *(cartesian[c]/R)
                             *diaOverlapAOs1stDerivInDiaFrame[k][l];
                     temp2 += rotMat1stDerivs[i][k][c] 
                             *rotMat[j][l]
                             *diaOverlapAOsInDiaFrame[k][l];
                     temp3 += rotMat[i][k] 
                             *rotMat1stDerivs[j][l][c]
                             *diaOverlapAOsInDiaFrame[k][l];
                  }
               }
               diatomicOverlapAOs1stDerivs[i][j][c] = temp1 + temp2 + temp3;
            }
         }
      }
      */

   }
   catch(MolDSException ex){
      this->FreeDiatomicOverlapAOs1stDeriTemps(&diaOverlapAOsInDiaFrame,
                                               &diaOverlapAOs1stDerivInDiaFrame,
                                               &rotMat,
                                               &rotMat1stDerivs,
                                               &rotatedDiatomicOverlap,
                                               &tmpRotMat1stDeriv,
                                               &tmpMatrix);
      throw ex;
   }
   // free
   this->FreeDiatomicOverlapAOs1stDeriTemps(&diaOverlapAOsInDiaFrame,
                                            &diaOverlapAOs1stDerivInDiaFrame,
                                            &rotMat,
                                            &rotMat1stDerivs,
                                            &rotatedDiatomicOverlap,
                                            &tmpRotMat1stDeriv,
                                            &tmpMatrix);
}

void Cndo2::CalcDiatomicOverlapAOs1stDerivatives(double*** diatomicOverlapAOs1stDerivs, 
                                                 int indexAtomA, 
                                                 int indexAtomB) const{
   this->CalcDiatomicOverlapAOs1stDerivatives(diatomicOverlapAOs1stDerivs,
                                              *this->molecule->GetAtom(indexAtomA),
                                              *this->molecule->GetAtom(indexAtomB));
}

// Second derivative of diatomic overlapAOs integrals between AOs in space fixed flame.
// The OverlapAOs matrix is S_{\mu\nu} in (3.74) in J. A. Pople book.
// Note that this method can not treat d-obitals 
// because CalcRotatingMatrix1stDerivatives can not treat d-orbitals.
void Cndo2::CalcDiatomicOverlapAOs2ndDerivatives(double**** diatomicOverlapAOs2ndDerivs, 
                                                 const Atom& atomA, 
                                                 const Atom& atomB) const{
   double cartesian[CartesianType_end] = {atomA.GetXyz()[XAxis] - atomB.GetXyz()[XAxis], 
                                          atomA.GetXyz()[YAxis] - atomB.GetXyz()[YAxis],
                                          atomA.GetXyz()[ZAxis] - atomB.GetXyz()[ZAxis]};
   double R = sqrt( pow(cartesian[XAxis],2.0) + 
                    pow(cartesian[YAxis],2.0) + 
                    pow(cartesian[ZAxis],2.0) );
   
   double**   diaOverlapAOsInDiaFrame = NULL;  // diatomic overlapAOs in diatomic frame
   double**   diaOverlapAOs1stDerivInDiaFrame = NULL; // first derivative of the diaOverlapAOs. This derivative is related to the distance between two atoms.
   double**   diaOverlapAOs2ndDerivInDiaFrame = NULL; // second derivative of the diaOverlapAOs. This derivative is related to the distance between two atoms.
   double**   rotMat = NULL; // rotating Matrix from the diatomic frame to space fixed frame.
   double***  rotMat1stDerivatives = NULL; //first derivatives of the rotMat
   double**** rotMat2ndDerivatives = NULL; //second derivatives of the rotMat
   double***  tempDiaOverlapAOs1stDerivs = NULL; // first derivatives of the diaOverlapAOs. This derivatives are related to the all Cartesian coordinates.
   double**** tempDiaOverlapAOs2ndDerivs = NULL; //sedond derivatives of the diaOverlapAOs. This derivatives are related to the all Cartesian coordinates.

   try{
      this->MallocDiatomicOverlapAOs2ndDeriTemps(&diaOverlapAOsInDiaFrame,
                                                 &diaOverlapAOs1stDerivInDiaFrame,
                                                 &diaOverlapAOs2ndDerivInDiaFrame,
                                                 &rotMat,
                                                 &rotMat1stDerivatives,
                                                 &rotMat2ndDerivatives,
                                                 &tempDiaOverlapAOs1stDerivs,
                                                 &tempDiaOverlapAOs2ndDerivs);
      this->CalcDiatomicOverlapAOsInDiatomicFrame(diaOverlapAOsInDiaFrame, atomA, atomB);
      this->CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(diaOverlapAOs1stDerivInDiaFrame, atomA, atomB);
      this->CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(diaOverlapAOs2ndDerivInDiaFrame, atomA, atomB);
      this->CalcRotatingMatrix(rotMat, atomA, atomB);
      this->CalcRotatingMatrix1stDerivatives(rotMat1stDerivatives, atomA, atomB);
      this->CalcRotatingMatrix2ndDerivatives(rotMat2ndDerivatives, atomA, atomB);

      // calculate each element of first derivatives
      for(int i=0; i<OrbitalType_end; i++){
         for(int j=0; j<OrbitalType_end; j++){
            for(int dimA1=0; dimA1<CartesianType_end; dimA1++){
               tempDiaOverlapAOs1stDerivs[i][j][dimA1] = (cartesian[dimA1]/R)*diaOverlapAOs1stDerivInDiaFrame[i][j];
            }
         }
      }

      // calculate each element of second derivatives
      for(int i=0; i<OrbitalType_end; i++){
         for(int j=0; j<OrbitalType_end; j++){
            for(int dimA1=XAxis; dimA1<CartesianType_end; dimA1++){
               for(int dimA2=XAxis; dimA2<CartesianType_end; dimA2++){
                  tempDiaOverlapAOs2ndDerivs[i][j][dimA1][dimA2] 
                     = this->Get2ndDerivativeElementFromDistanceDerivatives(diaOverlapAOs1stDerivInDiaFrame[i][j],
                                                                            diaOverlapAOs2ndDerivInDiaFrame[i][j],
                                                                            static_cast<CartesianType>(dimA1),
                                                                            static_cast<CartesianType>(dimA2),
                                                                            cartesian,
                                                                            R);
               }
            }
         }
      }

      // rotate
      for(int i=0; i<OrbitalType_end; i++){
         for(int j=0; j<OrbitalType_end; j++){
            for(int dimA1=XAxis; dimA1<CartesianType_end; dimA1++){
               for(int dimA2=XAxis; dimA2<CartesianType_end; dimA2++){
                  diatomicOverlapAOs2ndDerivs[i][j][dimA1][dimA2] = 0.0;
               
                  double temp1 = 0.0, temp2=0.0, temp3 = 0.0;
                  double temp4 = 0.0, temp5=0.0, temp6 = 0.0;
                  double temp7 = 0.0, temp8=0.0, temp9 = 0.0;
                  for(int k=0; k<OrbitalType_end; k++){
                     for(int l=0; l<OrbitalType_end; l++){
              
                        temp1 += rotMat2ndDerivatives   [i][k][dimA1][dimA2]
                                *rotMat                 [j][l]
                                *diaOverlapAOsInDiaFrame[k][l];
                        temp2 += rotMat                 [i][k]
                                *rotMat2ndDerivatives   [j][l][dimA1][dimA2]
                                *diaOverlapAOsInDiaFrame[k][l];
                        temp3 += rotMat                 [i][k]
                                *rotMat                 [j][l]
                                *tempDiaOverlapAOs2ndDerivs[k][l][dimA1][dimA2];
                        temp4 += rotMat1stDerivatives   [i][k][dimA1] 
                                *rotMat1stDerivatives   [j][l][dimA2]
                                *diaOverlapAOsInDiaFrame[k][l];
                        temp5 += rotMat1stDerivatives   [i][k][dimA1] 
                                *rotMat                 [j][l]
                                *tempDiaOverlapAOs1stDerivs[k][l][dimA2];
                        temp6 += rotMat1stDerivatives   [i][k][dimA2] 
                                *rotMat1stDerivatives   [j][l][dimA1]
                                *diaOverlapAOsInDiaFrame[k][l];
                        temp7 += rotMat                 [i][k] 
                                *rotMat1stDerivatives   [j][l][dimA1]
                                *tempDiaOverlapAOs1stDerivs[k][l][dimA2];
                        temp8 += rotMat1stDerivatives   [i][k][dimA2] 
                                *rotMat                 [j][l]
                                *tempDiaOverlapAOs1stDerivs[k][l][dimA1];
                        temp9 += rotMat                 [i][k] 
                                *rotMat1stDerivatives   [j][l][dimA2]
                                *tempDiaOverlapAOs1stDerivs[k][l][dimA1];
                     }
                  }

                  diatomicOverlapAOs2ndDerivs[i][j][dimA1][dimA2] = temp1+temp2+temp3 
                                                                   +temp4+temp5+temp6 
                                                                   +temp7+temp8+temp9;
               }
            }
         }
      }
      
   }
   catch(MolDSException ex){
      this->FreeDiatomicOverlapAOs2ndDeriTemps(&diaOverlapAOsInDiaFrame,
                                               &diaOverlapAOs1stDerivInDiaFrame,
                                               &diaOverlapAOs2ndDerivInDiaFrame,
                                               &rotMat,
                                               &rotMat1stDerivatives,
                                               &rotMat2ndDerivatives,
                                               &tempDiaOverlapAOs1stDerivs,
                                               &tempDiaOverlapAOs2ndDerivs);
      throw ex;
   }
   // free
   this->FreeDiatomicOverlapAOs2ndDeriTemps(&diaOverlapAOsInDiaFrame,
                                            &diaOverlapAOs1stDerivInDiaFrame,
                                            &diaOverlapAOs2ndDerivInDiaFrame,
                                            &rotMat,
                                            &rotMat1stDerivatives,
                                            &rotMat2ndDerivatives,
                                            &tempDiaOverlapAOs1stDerivs,
                                            &tempDiaOverlapAOs2ndDerivs);
   /*
   for(int i=0; i<OrbitalType_end; i++){
      for(int j=0; j<OrbitalType_end; j++){
         for(int dimA1=0; dimA1<CartesianType_end; dimA1++){
            for(int dimA2=0; dimA2<CartesianType_end; dimA2++){
               printf("i=%d j=%d dimA1=%d dimA2=%d: %e\n",i,j,dimA1,dimA2,overlapAOs2ndDeri[i][j][dimA1][dimA2]);
            }
         }
      }
   }
   */
}

void Cndo2::CalcDiatomicOverlapAOs2ndDerivatives(double**** diatomicOverlapAOs2ndDerivs, 
                                                 int indexAtomA, 
                                                 int indexAtomB) const{
   this->CalcDiatomicOverlapAOs2ndDerivatives(diatomicOverlapAOs2ndDerivs,
                                              *this->molecule->GetAtom(indexAtomA),
                                              *this->molecule->GetAtom(indexAtomB));
}

double Cndo2::Get2ndDerivativeElementFromDistanceDerivatives(double firstDistanceDeri,
                                                             double secondDistanceDeri,
                                                             CartesianType axisA1,
                                                             CartesianType axisA2,
                                                             double* cartesian,
                                                             double rAB) const{
   double value=0.0;               
   if(axisA1 != axisA2){
      value = -1.0*pow(rAB, -3.0)*firstDistanceDeri;
      value += pow(rAB, -2.0)*secondDistanceDeri;
      value *= cartesian[axisA1]*cartesian[axisA2];
   }
   else{
      value = (pow(rAB,2.0) - pow(cartesian[axisA1],2.0))*pow(rAB, -3.0)*firstDistanceDeri;
      value += pow(cartesian[axisA1]/rAB, 2.0)*secondDistanceDeri;
   }
   return value;
}

void Cndo2::MallocDiatomicOverlapAOs1stDeriTemps(double***  diaOverlapAOsInDiaFrame, 
                                                 double***  diaOverlapAOs1stDerivInDiaFrame,
                                                 double***  rotMat,
                                                 double**** rotMat1stDerivs,
                                                 double***  rotatedDiatomicOverlap,
                                                 double***  tmpRotMat1stDeriv,
                                                 double***  tmpMatrix) const{
   MallocerFreer::GetInstance()->Malloc<double>(diaOverlapAOsInDiaFrame,         OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diaOverlapAOs1stDerivInDiaFrame, OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Malloc<double>(rotMat,                          OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Malloc<double>(rotMat1stDerivs,                 OrbitalType_end, OrbitalType_end, CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(rotatedDiatomicOverlap,          OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Malloc<double>(tmpMatrix,                       OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Malloc<double>(tmpRotMat1stDeriv,               OrbitalType_end, OrbitalType_end);
}

void Cndo2::MallocDiatomicOverlapAOs2ndDeriTemps(double*** diaOverlapAOsInDiaFrame, 
                                                 double*** diaOverlapAOs1stDerivInDiaFrame,
                                                 double*** diaOverlapAOs2ndDerivInDiaFrame,
                                                 double***   rotMat,
                                                 double****  rotMat1stDerivs,
                                                 double***** rotMat2ndDerivs,
                                                 double****  tempDiaOverlapAOs1stDerivs,
                                                 double***** tempDiaOverlapAOs2ndDerivs) const{
   MallocerFreer::GetInstance()->Malloc<double>(diaOverlapAOsInDiaFrame,         OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diaOverlapAOs1stDerivInDiaFrame, OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diaOverlapAOs2ndDerivInDiaFrame, OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Malloc<double>(rotMat,                       OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Malloc<double>(rotMat1stDerivs,              OrbitalType_end, OrbitalType_end, CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(rotMat2ndDerivs,              OrbitalType_end, OrbitalType_end, CartesianType_end, CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(tempDiaOverlapAOs1stDerivs,      OrbitalType_end, OrbitalType_end, CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(tempDiaOverlapAOs2ndDerivs,      OrbitalType_end, OrbitalType_end, CartesianType_end, CartesianType_end);
}

void Cndo2::FreeDiatomicOverlapAOs1stDeriTemps(double*** diaOverlapAOsInDiaFrame, 
                                               double*** diaOverlapAOs1stDerivInDiaFrame,
                                               double*** rotMat,
                                               double**** rotMat1stDerivs,
                                               double***  rotatedDiatomicOverlap,
                                               double***  tmpRotMat1stDeriv,
                                               double***  tmpMatrix) const{
   MallocerFreer::GetInstance()->Free<double>(diaOverlapAOsInDiaFrame,         OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Free<double>(diaOverlapAOs1stDerivInDiaFrame, OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Free<double>(rotMat,                          OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Free<double>(rotMat1stDerivs,                 OrbitalType_end, OrbitalType_end, CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(rotatedDiatomicOverlap,          OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Free<double>(tmpMatrix,                       OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Free<double>(tmpRotMat1stDeriv,               OrbitalType_end, OrbitalType_end);
}

void Cndo2::FreeDiatomicOverlapAOs2ndDeriTemps(double*** diaOverlapAOsInDiaFrame, 
                                               double*** diaOverlapAOs1stDerivInDiaFrame,
                                               double*** diaOverlapAOs2ndDerivInDiaFrame,
                                               double***   rotMat,
                                               double****  rotMat1stDerivs,
                                               double***** rotMat2ndDerivs,
                                               double****  tempDiaOverlapAOs1stDerivs,
                                               double***** tempDiaOverlapAOs2ndDerivs) const{
   MallocerFreer::GetInstance()->Free<double>(diaOverlapAOsInDiaFrame,         OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Free<double>(diaOverlapAOs1stDerivInDiaFrame, OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Free<double>(diaOverlapAOs2ndDerivInDiaFrame, OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Free<double>(rotMat,                       OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Free<double>(rotMat1stDerivs,              OrbitalType_end, OrbitalType_end, CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(rotMat2ndDerivs,              OrbitalType_end, OrbitalType_end, CartesianType_end, CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(tempDiaOverlapAOs1stDerivs,      OrbitalType_end, OrbitalType_end, CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(tempDiaOverlapAOs2ndDerivs,      OrbitalType_end, OrbitalType_end, CartesianType_end, CartesianType_end);
}

// calculate OverlapAOs matrix. E.g. S_{\mu\nu} in (3.74) in J. A. Pople book by GTO expansion.
// See Eqs. (28) - (32) in [DY_1977]
void Cndo2::CalcOverlapAOsByGTOExpansion(double** overlapAOs, 
                                         const Molecule& molecule, 
                                         STOnGType stonG) const{
   int totalAONumber = molecule.GetTotalNumberAOs();
   int totalAtomNumber = molecule.GetNumberAtoms();

   // calculation overlapAOs matrix
   for(int mu=0; mu<totalAONumber; mu++){
      overlapAOs[mu][mu] = 1.0;
   }

   stringstream ompErrors;
#pragma omp parallel for schedule(auto) 
   for(int A=0; A<totalAtomNumber; A++){
      try{
         const Atom& atomA = *molecule.GetAtom(A);
         int firstAOIndexAtomA = atomA.GetFirstAOIndex();
         for(int B=A+1; B<totalAtomNumber; B++){
            const Atom& atomB = *molecule.GetAtom(B);
            int firstAOIndexAtomB = atomB.GetFirstAOIndex();
            for(int a=0; a<atomA.GetValenceSize(); a++){
               for(int b=0; b<atomB.GetValenceSize(); b++){
                  int mu = firstAOIndexAtomA + a;      
                  int nu = firstAOIndexAtomB + b;      
                  double value = this->GetOverlapAOsElementByGTOExpansion(atomA, a, atomB, b, stonG);
                  overlapAOs[mu][nu] = value;
                  overlapAOs[nu][mu] = value;
               }
            }
         }
      }
      catch(MolDSException ex){
#pragma omp critical
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
   /* 
   this->OutputLog("overlapAOs matrix by STOnG\n"); 
   for(int o=0; o<molecule.GetTotalNumberAOs(); o++){
      for(int p=0; p<molecule.GetTotalNumberAOs(); p++){
         this->OutputLog(boost::format("%lf\t") % overlapAOs[o][p]);
      }
      this->OutputLog("\n");
   }
   this->OutputLog("\n");
   */   
}

// calculate elements of overlapAOs matrix. 
// E.g. S_{\mu\nu} in (3.74) in J. A. Pople book by GTO expansion.
// See Eqs. (28) - (32) in [DY_1977]
double Cndo2::GetOverlapAOsElementByGTOExpansion(const Atom& atomA, int valenceIndexA, 
                                                 const Atom& atomB, int valenceIndexB,
                                                 STOnGType stonG) const{
   double value = 0.0;
   double dx = atomA.GetXyz()[XAxis] - atomB.GetXyz()[XAxis];
   double dy = atomA.GetXyz()[YAxis] - atomB.GetXyz()[YAxis];
   double dz = atomA.GetXyz()[ZAxis] - atomB.GetXyz()[ZAxis];
   double rAB = sqrt( pow(dx, 2.0) + pow(dy, 2.0) + pow(dz,2.0) );
   ShellType shellTypeA = atomA.GetValenceShellType();
   ShellType shellTypeB = atomB.GetValenceShellType();
   OrbitalType valenceOrbitalA = atomA.GetValence(valenceIndexA);
   OrbitalType valenceOrbitalB = atomB.GetValence(valenceIndexB);
   double orbitalExponentA = atomA.GetOrbitalExponent(atomA.GetValenceShellType(), 
                                                      valenceOrbitalA, 
                                                      this->theory);
   double orbitalExponentB = atomB.GetOrbitalExponent(atomB.GetValenceShellType(), 
                                                      valenceOrbitalB, 
                                                      this->theory);
   double gaussianExponentA = 0.0;
   double gaussianExponentB = 0.0;

   double temp = 0.0;
   for(int i=0; i<=stonG; i++){
      for(int j=0; j<=stonG; j++){
         temp = GTOExpansionSTO::GetInstance()->GetCoefficient(stonG, 
                                                               shellTypeA, 
                                                               valenceOrbitalA, 
                                                               i); 
         temp *= GTOExpansionSTO::GetInstance()->GetCoefficient(stonG, 
                                                                shellTypeB, 
                                                                valenceOrbitalB, 
                                                                j); 
         gaussianExponentA = pow(orbitalExponentA, 2.0) *
                             GTOExpansionSTO::GetInstance()->GetExponent(stonG, 
                                                                         shellTypeA, 
                                                                         valenceOrbitalA, 
                                                                         i);
         gaussianExponentB = pow(orbitalExponentB, 2.0) *
                             GTOExpansionSTO::GetInstance()->GetExponent(stonG, 
                                                                         shellTypeB, 
                                                                         valenceOrbitalB, 
                                                                         j);
         temp *= this->GetGaussianOverlapAOs(atomA.GetAtomType(), 
                                             valenceOrbitalA, 
                                             gaussianExponentA, 
                                             atomB.GetAtomType(), 
                                             valenceOrbitalB, 
                                             gaussianExponentB,
                                             dx, dy, dz, rAB);
         value += temp;
      }
   }
   return value;
}

// Calculate gaussian overlapAOs integrals of Sa and Sb.
// That is, calculate (S_A|S_B). See Eq. (28) in [DY_1977].
double Cndo2::GetGaussianOverlapAOsSASB(double gaussianExponentA, 
                                        double gaussianExponentB,
                                        double rAB) const{
   double value;
   double temp1 = 0.0;
   double temp2 = 0.0;
   temp1 = 2.0*pow(gaussianExponentA*gaussianExponentB, 0.5)
            /(gaussianExponentA+gaussianExponentB);
   temp2 = -1.0* gaussianExponentA*gaussianExponentB
            /(gaussianExponentA+gaussianExponentB);
   value = pow(temp1, 1.5)*exp(temp2*rAB*rAB);
   return value;
}


// calculate gaussian overlapAOs integrals. 
// See Eqs. (28) - (32) in [DY_1977].
// Although d-orbital is not calucluated in [DY_1977],
// the way to calculate overlapAOs related to d-orbital is 
// same to the one written in [DY_1977].
double Cndo2::GetGaussianOverlapAOs(AtomType atomTypeA, 
                                    OrbitalType valenceOrbitalA, 
                                    double gaussianExponentA, 
                                    AtomType atomTypeB, 
                                    OrbitalType valenceOrbitalB, 
                                    double gaussianExponentB,
                                    double dx, double dy, double dz, 
                                    double rAB) const{
   double overlapSASB = this->GetGaussianOverlapAOsSASB(gaussianExponentA, gaussianExponentB, rAB);
   return this->GetGaussianOverlapAOs(atomTypeA, valenceOrbitalA, gaussianExponentA,
                                      atomTypeB, valenceOrbitalB, gaussianExponentB,
                                      dx, dy, dz, rAB, overlapSASB);
}
// calculate gaussian overlapAOs integrals. 
// See Eqs. (28) - (32) in [DY_1977].
// Although d-orbital is not calucluated in [DY_1977],
// the way to calculate overlapAOs related to d-orbital is 
// same to the one written in [DY_1977].
double Cndo2::GetGaussianOverlapAOs(AtomType atomTypeA, 
                                    OrbitalType valenceOrbitalA, 
                                    double gaussianExponentA, 
                                    AtomType atomTypeB, 
                                    OrbitalType valenceOrbitalB, 
                                    double gaussianExponentB,
                                    double dx, double dy, double dz, 
                                    double rAB,
                                    double overlapSASB) const{

   double value = 0.0;
   if(valenceOrbitalA == s && valenceOrbitalB == s){
      value = 1.0;
   }

   else if(valenceOrbitalA == s && valenceOrbitalB == px){
      value = 2.0*gaussianExponentA*pow(gaussianExponentB, 0.5)*dx;
      value /= (gaussianExponentA+gaussianExponentB);
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == py){
      value = 2.0*gaussianExponentA*pow(gaussianExponentB, 0.5)*dy;
      value /= (gaussianExponentA+gaussianExponentB);
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == pz){
      value = 2.0*gaussianExponentA*pow(gaussianExponentB, 0.5)*dz;
      value /= (gaussianExponentA+gaussianExponentB);
   }

   else if(valenceOrbitalA == px && valenceOrbitalB == s){
      value = -2.0*pow(gaussianExponentA, 0.5)*gaussianExponentB*dx;
      value /= (gaussianExponentA+gaussianExponentB);
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == s){
      value = -2.0*pow(gaussianExponentA, 0.5)*gaussianExponentB*dy;
      value /= (gaussianExponentA+gaussianExponentB);
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == s){
      value = -2.0*pow(gaussianExponentA, 0.5)*gaussianExponentB*dz;
      value /= (gaussianExponentA+gaussianExponentB);
   }

   else if(valenceOrbitalA == px && valenceOrbitalB == px){
      double temp = 0.0;
      temp = -1.0*pow(dx,2.0)*gaussianExponentA*gaussianExponentB;
      temp /= (gaussianExponentA+gaussianExponentB);
      temp += 0.5;
      value = 4.0*pow(gaussianExponentA*gaussianExponentB, 0.5);
      value /= (gaussianExponentA+gaussianExponentB);
      value *= temp;
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == py){
      value = -4.0*pow(gaussianExponentA*gaussianExponentB, 1.5);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= dx*dy;
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == pz){
      value = -4.0*pow(gaussianExponentA*gaussianExponentB, 1.5);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= dx*dz;
   }

   else if(valenceOrbitalA == py && valenceOrbitalB == px){
      value = -4.0*pow(gaussianExponentA*gaussianExponentB, 1.5);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= dy*dx;
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == py){
      double temp = 0.0;
      temp = -1.0*pow(dy,2.0)*gaussianExponentA*gaussianExponentB;
      temp /= (gaussianExponentA+gaussianExponentB);
      temp += 0.5;
      value = 4.0*pow(gaussianExponentA*gaussianExponentB, 0.5);
      value /= (gaussianExponentA+gaussianExponentB);
      value *= temp;
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == pz){
      value = -4.0*pow(gaussianExponentA*gaussianExponentB, 1.5);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= dy*dz;
   }

   else if(valenceOrbitalA == pz && valenceOrbitalB == px){
      value = -4.0*pow(gaussianExponentA*gaussianExponentB, 1.5);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= dz*dx;
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == py){
      value = -4.0*pow(gaussianExponentA*gaussianExponentB, 1.5);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= dz*dy;
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == pz){
      double temp = 0.0;
      temp = -1.0*pow(dz,2.0)*gaussianExponentA*gaussianExponentB;
      temp /= (gaussianExponentA+gaussianExponentB);
      temp += 0.5;
      value = 4.0*pow(gaussianExponentA*gaussianExponentB, 0.5);
      value /= (gaussianExponentA+gaussianExponentB);
      value *= temp;
   }

   else if(valenceOrbitalA == dxy && valenceOrbitalB == s){
      value = 4.0*gaussianExponentA;
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= gaussianExponentB*dx;
      value *= gaussianExponentB*dy;
   }
   else if(valenceOrbitalA == dyz && valenceOrbitalB == s){
      value = 4.0*gaussianExponentA;
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= gaussianExponentB*dy;
      value *= gaussianExponentB*dz;
   }
   else if(valenceOrbitalA == dzx && valenceOrbitalB == s){
      value = 4.0*gaussianExponentA;
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= gaussianExponentB*dz;
      value *= gaussianExponentB*dx;
   }

   else if(valenceOrbitalA == s && valenceOrbitalB == dxy){
      value = 4.0*gaussianExponentB;
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= gaussianExponentA*dx;
      value *= gaussianExponentA*dy;
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == dyz){
      value = 4.0*gaussianExponentB;
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= gaussianExponentA*dy;
      value *= gaussianExponentA*dz;
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == dzx){
      value = 4.0*gaussianExponentB;
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= gaussianExponentA*dz;
      value *= gaussianExponentA*dx;
   }

   else if(valenceOrbitalA == dxy && valenceOrbitalB == px){
      double temp1 = -0.5*gaussianExponentB*dy;
      double temp2 = gaussianExponentB*dx*gaussianExponentA*dx*gaussianExponentB*dy;
      temp2 /= gaussianExponentA+gaussianExponentB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentA*pow(gaussianExponentB, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == dxy && valenceOrbitalB == py){
      double temp1 = -0.5*gaussianExponentB*dx;
      double temp2 = gaussianExponentB*dy*gaussianExponentA*dy*gaussianExponentB*dx;
      temp2 /= gaussianExponentA+gaussianExponentB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentA*pow(gaussianExponentB, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == dyz && valenceOrbitalB == py){
      double temp1 = -0.5*gaussianExponentB*dz;
      double temp2 = gaussianExponentB*dy*gaussianExponentA*dy*gaussianExponentB*dz;
      temp2 /= gaussianExponentA+gaussianExponentB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentA*pow(gaussianExponentB, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == dyz && valenceOrbitalB == pz){
      double temp1 = -0.5*gaussianExponentB*dy;
      double temp2 = gaussianExponentB*dz*gaussianExponentA*dz*gaussianExponentB*dy;
      temp2 /= gaussianExponentA+gaussianExponentB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentA*pow(gaussianExponentB, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == dzx && valenceOrbitalB == pz){
      double temp1 = -0.5*gaussianExponentB*dx;
      double temp2 = gaussianExponentB*dz*gaussianExponentA*dz*gaussianExponentB*dx;
      temp2 /= gaussianExponentA+gaussianExponentB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentA*pow(gaussianExponentB, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == dzx && valenceOrbitalB == px){
      double temp1 = -0.5*gaussianExponentB*dz;
      double temp2 = gaussianExponentB*dx*gaussianExponentA*dx*gaussianExponentB*dz;
      temp2 /= gaussianExponentA+gaussianExponentB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentA*pow(gaussianExponentB, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }

   else if(valenceOrbitalA == px && valenceOrbitalB == dxy){
      double temp1 = 0.5*gaussianExponentA*dy;
      double temp2 = -1.0*gaussianExponentA*dx*gaussianExponentB*dx*gaussianExponentA*dy;
      temp2 /= gaussianExponentA+gaussianExponentB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentB*pow(gaussianExponentA, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == dxy){
      double temp1 = 0.5*gaussianExponentA*dx;
      double temp2 = -1.0*gaussianExponentA*dy*gaussianExponentB*dy*gaussianExponentA*dx;
      temp2 /= gaussianExponentA+gaussianExponentB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentB*pow(gaussianExponentA, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == dyz){
      double temp1 = 0.5*gaussianExponentA*dz;
      double temp2 = -1.0*gaussianExponentA*dy*gaussianExponentB*dy*gaussianExponentA*dz;
      temp2 /= gaussianExponentA+gaussianExponentB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentB*pow(gaussianExponentA, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == dyz){
      double temp1 = 0.5*gaussianExponentA*dy;
      double temp2 = -1.0*gaussianExponentA*dz*gaussianExponentB*dz*gaussianExponentA*dy;
      temp2 /= gaussianExponentA+gaussianExponentB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentB*pow(gaussianExponentA, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == dzx){
      double temp1 = 0.5*gaussianExponentA*dx;
      double temp2 = -1.0*gaussianExponentA*dz*gaussianExponentB*dz*gaussianExponentA*dx;
      temp2 /= gaussianExponentA+gaussianExponentB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentB*pow(gaussianExponentA, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == dzx){
      double temp1 = 0.5*gaussianExponentA*dz;
      double temp2 = -1.0*gaussianExponentA*dx*gaussianExponentB*dx*gaussianExponentA*dz;
      temp2 /= gaussianExponentA+gaussianExponentB;
      value = temp1 + temp2;
      value *= 8.0*gaussianExponentB*pow(gaussianExponentA, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }

   else if(valenceOrbitalA == dxy && valenceOrbitalB == pz){
      value = 8.0*gaussianExponentA*pow(gaussianExponentB, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 3.0);
      value *= gaussianExponentB*dx*gaussianExponentB*dy*gaussianExponentA*dz;
   }
   else if(valenceOrbitalA == dyz && valenceOrbitalB == px){
      value = 8.0*gaussianExponentA*pow(gaussianExponentB, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 3.0);
      value *= gaussianExponentB*dy*gaussianExponentB*dz*gaussianExponentA*dx;
   }
   else if(valenceOrbitalA == dzx && valenceOrbitalB == py){
      value = 8.0*gaussianExponentA*pow(gaussianExponentB, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 3.0);
      value *= gaussianExponentB*dz*gaussianExponentB*dx*gaussianExponentA*dy;
   }

   else if(valenceOrbitalA == pz && valenceOrbitalB == dxy){
      value = -8.0*gaussianExponentB*pow(gaussianExponentA, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 3.0);
      value *= gaussianExponentA*dx*gaussianExponentA*dy*gaussianExponentB*dz;
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == dyz){
      value = -8.0*gaussianExponentB*pow(gaussianExponentA, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 3.0);
      value *= gaussianExponentA*dy*gaussianExponentA*dz*gaussianExponentB*dx;
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == dzx){
      value = -8.0*gaussianExponentB*pow(gaussianExponentA, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 3.0);
      value *= gaussianExponentA*dz*gaussianExponentA*dx*gaussianExponentB*dy;
   }

   else if(valenceOrbitalA == dxxyy && valenceOrbitalB == s){
      value = 2.0*gaussianExponentA;
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= pow(gaussianExponentB*dx, 2.0) - pow(gaussianExponentB*dy, 2.0);
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == dxxyy){
      value = 2.0*gaussianExponentB;
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= pow(gaussianExponentA*dx, 2.0) - pow(gaussianExponentA*dy, 2.0);
   }
   else if(valenceOrbitalA == dxxyy && valenceOrbitalB == px){
      value = gaussianExponentB*dx;
      value -= pow(gaussianExponentB*dx, 2.0)*gaussianExponentA*dx/(gaussianExponentA+gaussianExponentB);
      value += pow(gaussianExponentB*dy, 2.0)*gaussianExponentA*dx/(gaussianExponentA+gaussianExponentB);
      value *= -4.0*gaussianExponentA*pow(gaussianExponentB, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == dxxyy){
      value = gaussianExponentA*dx;
      value -= pow(gaussianExponentA*dx, 2.0)*gaussianExponentB*dx/(gaussianExponentA+gaussianExponentB);
      value += pow(gaussianExponentA*dy, 2.0)*gaussianExponentB*dx/(gaussianExponentA+gaussianExponentB);
      value *= 4.0*gaussianExponentB*pow(gaussianExponentA, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == dxxyy && valenceOrbitalB == py){
      value = gaussianExponentB*dy;
      value += pow(gaussianExponentB*dx, 2.0)*gaussianExponentA*dy/(gaussianExponentA+gaussianExponentB);
      value -= pow(gaussianExponentB*dy, 2.0)*gaussianExponentA*dy/(gaussianExponentA+gaussianExponentB);
      value *= 4.0*gaussianExponentA*pow(gaussianExponentB, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == dxxyy){
      value = gaussianExponentA*dy;
      value += pow(gaussianExponentA*dx, 2.0)*gaussianExponentB*dy/(gaussianExponentA+gaussianExponentB);
      value -= pow(gaussianExponentA*dy, 2.0)*gaussianExponentB*dy/(gaussianExponentA+gaussianExponentB);
      value *= -4.0*gaussianExponentB*pow(gaussianExponentA, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 2.0);
   }
   else if(valenceOrbitalA == dxxyy && valenceOrbitalB == pz){
      value = pow(gaussianExponentB*dx, 2.0) - pow(gaussianExponentB*dy, 2.0);
      value *= gaussianExponentA*dz;
      value *= 4.0*gaussianExponentA*pow(gaussianExponentB, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 3.0);
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == dxxyy){
      value = pow(gaussianExponentA*dx, 2.0) - pow(gaussianExponentA*dy, 2.0);
      value *= gaussianExponentB*dz;
      value *= -4.0*gaussianExponentB*pow(gaussianExponentA, 0.5);
      value /= pow(gaussianExponentA+gaussianExponentB, 3.0);
   }

   else if(valenceOrbitalA == dzz && valenceOrbitalB == s){
      double temp = 0.0;
      temp = 2.0*pow(gaussianExponentB*dz, 2.0) 
            -    pow(gaussianExponentB*dx, 2.0) 
            -    pow(gaussianExponentB*dy, 2.0);
      value = 2.0*gaussianExponentA/sqrt(3.0);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= temp;
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == dzz){
      double temp = 0.0;
      temp = 2.0*pow(gaussianExponentA*dz, 2.0) 
            -    pow(gaussianExponentA*dx, 2.0) 
            -    pow(gaussianExponentA*dy, 2.0);
      value = 2.0*gaussianExponentB/sqrt(3.0);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
      value *= temp;
   }
   else if(valenceOrbitalA == dzz && valenceOrbitalB == px){
      double temp = 0.0;
      temp = gaussianExponentB*dx;
      temp += 2.0*pow(gaussianExponentB*dz, 2.0)*gaussianExponentA*dx/(gaussianExponentA+gaussianExponentB);
      temp -=     pow(gaussianExponentB*dx, 2.0)*gaussianExponentA*dx/(gaussianExponentA+gaussianExponentB);
      temp -=     pow(gaussianExponentB*dy, 2.0)*gaussianExponentA*dx/(gaussianExponentA+gaussianExponentB);
      value = temp;
      value *= 4.0*gaussianExponentA*pow(gaussianExponentB, 0.5)/sqrt(3.0);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == dzz){
      double temp = 0.0;
      temp = gaussianExponentA*dx;
      temp += 2.0*pow(gaussianExponentA*dz, 2.0)*gaussianExponentB*dx/(gaussianExponentA+gaussianExponentB);
      temp -=     pow(gaussianExponentA*dx, 2.0)*gaussianExponentB*dx/(gaussianExponentA+gaussianExponentB);
      temp -=     pow(gaussianExponentA*dy, 2.0)*gaussianExponentB*dx/(gaussianExponentA+gaussianExponentB);
      value = temp;
      value *= -4.0*gaussianExponentB*pow(gaussianExponentA, 0.5)/sqrt(3.0);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
   }
   else if(valenceOrbitalA == dzz && valenceOrbitalB == py){
      double temp = 0.0;
      temp = gaussianExponentB*dy;
      temp += 2.0*pow(gaussianExponentB*dz, 2.0)*gaussianExponentA*dy/(gaussianExponentA+gaussianExponentB);
      temp -=     pow(gaussianExponentB*dx, 2.0)*gaussianExponentA*dy/(gaussianExponentA+gaussianExponentB);
      temp -=     pow(gaussianExponentB*dy, 2.0)*gaussianExponentA*dy/(gaussianExponentA+gaussianExponentB);
      value = temp;
      value *= 4.0*gaussianExponentA*pow(gaussianExponentB, 0.5)/sqrt(3.0);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == dzz){
      double temp = 0.0;
      temp = gaussianExponentA*dy;
      temp += 2.0*pow(gaussianExponentA*dz, 2.0)*gaussianExponentB*dy/(gaussianExponentA+gaussianExponentB);
      temp -=     pow(gaussianExponentA*dx, 2.0)*gaussianExponentB*dy/(gaussianExponentA+gaussianExponentB);
      temp -=     pow(gaussianExponentA*dy, 2.0)*gaussianExponentB*dy/(gaussianExponentA+gaussianExponentB);
      value = temp;
      value *= -4.0*gaussianExponentB*pow(gaussianExponentA, 0.5)/sqrt(3.0);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
   }
   else if(valenceOrbitalA == dzz && valenceOrbitalB == pz){
      double temp = 0.0;
      temp = -2.0*gaussianExponentB*dz;
      temp += 2.0*pow(gaussianExponentB*dz, 2.0)*gaussianExponentA*dz/(gaussianExponentA+gaussianExponentB);
      temp -=     pow(gaussianExponentB*dx, 2.0)*gaussianExponentA*dz/(gaussianExponentA+gaussianExponentB);
      temp -=     pow(gaussianExponentB*dy, 2.0)*gaussianExponentA*dz/(gaussianExponentA+gaussianExponentB);
      value = temp;
      value *= 4.0*gaussianExponentA*pow(gaussianExponentB, 0.5)/sqrt(3.0);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == dzz){
      double temp = 0.0;
      temp = -2.0*gaussianExponentA*dz;
      temp += 2.0*pow(gaussianExponentA*dz, 2.0)*gaussianExponentB*dz/(gaussianExponentA+gaussianExponentB);
      temp -=     pow(gaussianExponentA*dx, 2.0)*gaussianExponentB*dz/(gaussianExponentA+gaussianExponentB);
      temp -=     pow(gaussianExponentA*dy, 2.0)*gaussianExponentB*dz/(gaussianExponentA+gaussianExponentB);
      value = temp;
      value *= -4.0*gaussianExponentB*pow(gaussianExponentA, 0.5)/sqrt(3.0);
      value *= pow(gaussianExponentA+gaussianExponentB, -2.0);
   }

   else if(valenceOrbitalA == dxy && valenceOrbitalB == dxy){
      double temp = 0.25;
      temp -= 0.5*gaussianExponentB*dx*gaussianExponentA*dx/(gaussianExponentA+gaussianExponentB);
      temp -= 0.5*gaussianExponentB*dy*gaussianExponentA*dy/(gaussianExponentA+gaussianExponentB);
      temp += gaussianExponentB*dx*gaussianExponentA*dx
             *gaussianExponentB*dy*gaussianExponentA*dy
             *pow(gaussianExponentA+gaussianExponentB, -2.0);
      value = 16.0*temp*gaussianExponentA*gaussianExponentB
             *pow(gaussianExponentA+gaussianExponentB, -2.0);
   }
   else if(valenceOrbitalA == dyz && valenceOrbitalB == dyz){
      double temp = 0.25;
      temp -= 0.5*gaussianExponentB*dy*gaussianExponentA*dy/(gaussianExponentA+gaussianExponentB);
      temp -= 0.5*gaussianExponentB*dz*gaussianExponentA*dz/(gaussianExponentA+gaussianExponentB);
      temp += gaussianExponentB*dy*gaussianExponentA*dy
             *gaussianExponentB*dz*gaussianExponentA*dz
             *pow(gaussianExponentA+gaussianExponentB, -2.0);
      value = 16.0*temp*gaussianExponentA*gaussianExponentB
             *pow(gaussianExponentA+gaussianExponentB, -2.0);
   }
   else if(valenceOrbitalA == dzx && valenceOrbitalB == dzx){
      double temp = 0.25;
      temp -= 0.5*gaussianExponentB*dz*gaussianExponentA*dz/(gaussianExponentA+gaussianExponentB);
      temp -= 0.5*gaussianExponentB*dx*gaussianExponentA*dx/(gaussianExponentA+gaussianExponentB);
      temp += gaussianExponentB*dz*gaussianExponentA*dz
             *gaussianExponentB*dx*gaussianExponentA*dx
             *pow(gaussianExponentA+gaussianExponentB, -2.0);
      value = 16.0*temp*gaussianExponentA*gaussianExponentB
             *pow(gaussianExponentA+gaussianExponentB, -2.0);
   }
   else if(valenceOrbitalA == dxxyy && valenceOrbitalB == dxxyy){
      double temp1 = 1.0;
      temp1 -= 2.0*gaussianExponentB*dx*gaussianExponentA*dx/(gaussianExponentA+gaussianExponentB);
      temp1 -= 2.0*gaussianExponentB*dy*gaussianExponentA*dy/(gaussianExponentA+gaussianExponentB);
      double temp2 = gaussianExponentA*gaussianExponentB*(pow(dx,2.0)-pow(dy,2.0))
             /(gaussianExponentA+gaussianExponentB);
      temp1 += pow(temp2,2.0);
      value = 4.0*temp1*gaussianExponentA*gaussianExponentB
             *pow(gaussianExponentA+gaussianExponentB, -2.0);
   }
   else if(valenceOrbitalA == dzz && valenceOrbitalB == dzz){
      double temp = 3.0;
      temp -= gaussianExponentA*gaussianExponentB
             /(gaussianExponentA+gaussianExponentB)
             *(8.0*pow(dz,2.0)+2.0*pow(dx,2.0)+2.0*pow(dy,2.0));
      temp += pow(gaussianExponentA*gaussianExponentB,2.0)
             *pow(gaussianExponentA+gaussianExponentB, -2.0)
             *(4.0*pow(dz,4.0)
                  +pow(dx,4.0)
                  +pow(dy,4.0)
              -4.0*pow(dx*dz,2.0)
              -4.0*pow(dy*dz,2.0)
              +2.0*pow(dx*dy,2.0));
      value = 4.0*temp*gaussianExponentA*gaussianExponentB
             *pow(gaussianExponentA+gaussianExponentB, -2.0)
             /3.0;
   }

   else if((valenceOrbitalA == dxy && valenceOrbitalB == dyz) ||
           (valenceOrbitalA == dyz && valenceOrbitalB == dxy)){
      double temp = 0.5;
      temp -= gaussianExponentA*gaussianExponentB*pow(dy,2.0)/(gaussianExponentA+gaussianExponentB);
      value = -16.0*pow(gaussianExponentA*gaussianExponentB,2.0)
             *pow(gaussianExponentA+gaussianExponentB,-3.0)
             *dx*dz*temp;
   }
   else if((valenceOrbitalA == dyz && valenceOrbitalB == dzx) ||
           (valenceOrbitalA == dzx && valenceOrbitalB == dyz)){
      double temp = 0.5;
      temp -= gaussianExponentA*gaussianExponentB*pow(dz,2.0)/(gaussianExponentA+gaussianExponentB);
      value = -16.0*pow(gaussianExponentA*gaussianExponentB,2.0)
             *pow(gaussianExponentA+gaussianExponentB,-3.0)
             *dy*dx*temp;
   }
   else if((valenceOrbitalA == dzx && valenceOrbitalB == dxy) ||
           (valenceOrbitalA == dxy && valenceOrbitalB == dzx)){
      double temp = 0.5;
      temp -= gaussianExponentA*gaussianExponentB*pow(dx,2.0)/(gaussianExponentA+gaussianExponentB);
      value = -16.0*pow(gaussianExponentA*gaussianExponentB,2.0)
             *pow(gaussianExponentA+gaussianExponentB,-3.0)
             *dz*dy*temp;
   }

   else if((valenceOrbitalA == dxxyy && valenceOrbitalB == dxy) ||
           (valenceOrbitalA == dxy && valenceOrbitalB == dxxyy)){
      double temp = 2.0*gaussianExponentA*gaussianExponentB;
      value = pow(temp,3.0)*(dy*pow(dx,3.0)-dx*pow(dy,3.0))
             *pow(gaussianExponentA+gaussianExponentB,-4.0);
   }
   else if((valenceOrbitalA == dxxyy && valenceOrbitalB == dyz) ||
           (valenceOrbitalA == dyz && valenceOrbitalB == dxxyy)){
      double temp = 2.0*gaussianExponentA*gaussianExponentB;
      value = pow(temp,3.0)*(dy*dz*(gaussianExponentA+gaussianExponentB)
                            /(gaussianExponentA*gaussianExponentB)
                             +(pow(dx,2.0)*dy*dz - pow(dy,3.0)*dz))
             *pow(gaussianExponentA+gaussianExponentB,-4.0);
   }
   else if((valenceOrbitalA == dxxyy && valenceOrbitalB == dzx) ||
           (valenceOrbitalA == dzx && valenceOrbitalB == dxxyy)){
      double temp = 2.0*gaussianExponentA*gaussianExponentB;
      value = -1.0*pow(temp,3.0)*(dx*dz*(gaussianExponentA+gaussianExponentB)
                            /(gaussianExponentA*gaussianExponentB)
                             +(pow(dy,2.0)*dx*dz - pow(dx,3.0)*dz))
             *pow(gaussianExponentA+gaussianExponentB,-4.0);
   }

   else if((valenceOrbitalA == dzz && valenceOrbitalB == dxy) ||
           (valenceOrbitalA == dxy && valenceOrbitalB == dzz)){
      double temp = 2.0*dx*dy*pow(dz,2.0) - pow(dx,3.0)*dy - dx*pow(dy,3.0);
      temp *= gaussianExponentA*gaussianExponentB/(gaussianExponentA+gaussianExponentB);
      temp += 2.0*dx*dy;
      value = 8.0*pow(gaussianExponentA*gaussianExponentB,2.0)*temp;
      value /= sqrt(3.0)*pow(gaussianExponentA+gaussianExponentB,3.0);
   }
   else if((valenceOrbitalA == dzz && valenceOrbitalB == dyz) ||
           (valenceOrbitalA == dyz && valenceOrbitalB == dzz)){
      double temp1 = -1.0*dy*dz;
      double temp2 = 2.0*dy*pow(dz,3.0) - pow(dy,3.0)*dz - pow(dx,2.0)*dy*dz;
      temp2 *= gaussianExponentA*gaussianExponentB/(gaussianExponentA+gaussianExponentB);
      temp1 += temp2;
      value = 8.0*pow(gaussianExponentA*gaussianExponentB,2.0)*temp1;
      value /= sqrt(3.0)*pow(gaussianExponentA+gaussianExponentB,3.0);
   }
   else if((valenceOrbitalA == dzz && valenceOrbitalB == dzx) ||
           (valenceOrbitalA == dzx && valenceOrbitalB == dzz)){
      double temp1 = -1.0*dx*dz;
      double temp2 = 2.0*dx*pow(dz,3.0) - pow(dx,3.0)*dz - pow(dy,2.0)*dx*dz;
      temp2 *= gaussianExponentA*gaussianExponentB/(gaussianExponentA+gaussianExponentB);
      temp1 += temp2;
      value = 8.0*pow(gaussianExponentA*gaussianExponentB,2.0)*temp1;
      value /= sqrt(3.0)*pow(gaussianExponentA+gaussianExponentB,3.0);
   }
   else if((valenceOrbitalA == dxxyy && valenceOrbitalB == dzz) ||
           (valenceOrbitalA == dzz && valenceOrbitalB == dxxyy)){
      double temp = 2.0*pow(dz,2.0)-pow(dx,2.0)-pow(dy,2.0);
      temp *= gaussianExponentA*gaussianExponentB/(gaussianExponentA+gaussianExponentB);
      temp += 2.0;
      value = 4.0*pow(gaussianExponentA*gaussianExponentB,2.0);
      value /= sqrt(3.0)*pow(gaussianExponentA+gaussianExponentB,3.0);
      value *= (pow(dx,2.0)-pow(dy,2.0))*temp;
   }

   else{
      stringstream ss;
      ss << this->errorMessageGetGaussianOverlapAOsBadOrbital;
      ss << this->errorMessageAtomA;
      ss << this->errorMessageAtomType << AtomTypeStr(atomTypeA) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(valenceOrbitalA) << endl;
      ss << this->errorMessageAtomB;
      ss << this->errorMessageAtomType << AtomTypeStr(atomTypeB) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(valenceOrbitalB) << endl;
      throw MolDSException(ss.str());
      value = 0.0;
   }
   value *= overlapSASB;
   return value;
}

// calculate elements of analytic first derivative of the overlapAOs matrix. 
// The derivative is carried out related to the coordinate of atom A.
// See Eqs. (34) - (44) in [DY_1977]
double Cndo2::GetOverlapAOsElement1stDerivativeByGTOExpansion(const Atom& atomA, 
                                                              int valenceIndexA, 
                                                              const Atom& atomB, 
                                                              int valenceIndexB,
                                                              STOnGType stonG, 
                                                              CartesianType axisA) const{

   double value = 0.0;
   double dx = atomA.GetXyz()[XAxis] - atomB.GetXyz()[XAxis];
   double dy = atomA.GetXyz()[YAxis] - atomB.GetXyz()[YAxis];
   double dz = atomA.GetXyz()[ZAxis] - atomB.GetXyz()[ZAxis];
   double rAB = sqrt( pow(dx, 2.0) + pow(dy, 2.0) + pow(dz,2.0) );
   ShellType shellTypeA = atomA.GetValenceShellType();
   ShellType shellTypeB = atomB.GetValenceShellType();
   OrbitalType valenceOrbitalA = atomA.GetValence(valenceIndexA);
   OrbitalType valenceOrbitalB = atomB.GetValence(valenceIndexB);
   double orbitalExponentA = atomA.GetOrbitalExponent(atomA.GetValenceShellType(), 
                                                      valenceOrbitalA, 
                                                      this->theory);
   double orbitalExponentB = atomB.GetOrbitalExponent(atomB.GetValenceShellType(), 
                                                      valenceOrbitalB, 
                                                      this->theory);
   double gaussianExponentA = 0.0;
   double gaussianExponentB = 0.0;

   double temp = 0.0;
   for(int i=0; i<=stonG; i++){
      for(int j=0; j<=stonG; j++){
         temp = GTOExpansionSTO::GetInstance()->GetCoefficient(stonG, 
                                                               shellTypeA, 
                                                               valenceOrbitalA, 
                                                               i); 
         temp *= GTOExpansionSTO::GetInstance()->GetCoefficient(stonG, 
                                                                shellTypeB, 
                                                                valenceOrbitalB, 
                                                                j); 
         gaussianExponentA = pow(orbitalExponentA, 2.0) 
                            *GTOExpansionSTO::GetInstance()->GetExponent(stonG, 
                                                                         shellTypeA, 
                                                                         valenceOrbitalA, 
                                                                         i);
         gaussianExponentB = pow(orbitalExponentB, 2.0)
                            *GTOExpansionSTO::GetInstance()->GetExponent(stonG, 
                                                                         shellTypeB, 
                                                                         valenceOrbitalB, 
                                                                         j);
         temp *= this->GetGaussianOverlapAOs1stDerivative(atomA.GetAtomType(), 
                                                          valenceOrbitalA, 
                                                          gaussianExponentA, 
                                                          atomB.GetAtomType(), 
                                                          valenceOrbitalB, 
                                                          gaussianExponentB,
                                                          dx, 
                                                          dy, 
                                                          dz, 
                                                          rAB, 
                                                          axisA);
         value += temp;
      }
   }
   return value;
}

// calculate first derivative of gaussian overlapAOs integrals. 
// See Eqs. (35) - (44) in [DY_1977]
double Cndo2::GetGaussianOverlapAOs1stDerivative(AtomType atomTypeA, 
                                                 OrbitalType valenceOrbitalA, 
                                                 double gaussianExponentA, 
                                                 AtomType atomTypeB, 
                                                 OrbitalType valenceOrbitalB, 
                                                 double gaussianExponentB,
                                                 double dx, double dy, double dz, double rAB, 
                                                 CartesianType axisA) const{
   double value = 0.0;

   if(valenceOrbitalA == s && valenceOrbitalB == s){
      double temp = -2.0*gaussianExponentA*gaussianExponentB
                     /(gaussianExponentA+gaussianExponentB);
      value = temp;
      if(axisA == XAxis){
         value *= dx;
      }
      else if(axisA == YAxis){
         value *= dy;
      }
      else if(axisA == ZAxis){
         value *= dz;
      }
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == px){
      double temp1 = 4.0*pow(gaussianExponentA,2.0)*pow(gaussianExponentB, 1.5)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      if(axisA == XAxis){
         double temp2 = 2.0*gaussianExponentA*pow(gaussianExponentB, 0.5)
                        /(gaussianExponentA+gaussianExponentB);
         value = temp2-temp1*dx*dx;
      }
      else if(axisA == YAxis){
         value = -1.0*temp1*dx*dy;
      }
      else if(axisA == ZAxis){
         value = -1.0*temp1*dx*dz;
      }
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == py){
      double temp1 = 4.0*pow(gaussianExponentA,2.0)*pow(gaussianExponentB, 1.5)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      if(axisA == XAxis){
         value = -1.0*temp1*dx*dy;
      }
      else if(axisA == YAxis){
         double temp2 = 2.0*gaussianExponentA*pow(gaussianExponentB, 0.5)
                        /(gaussianExponentA+gaussianExponentB);
         value = temp2-temp1*dy*dy;
      }
      else if(axisA == ZAxis){
         value = -1.0*temp1*dy*dz;
      }
   }
   else if(valenceOrbitalA == s && valenceOrbitalB == pz){
      double temp1 = 4.0*pow(gaussianExponentA,2.0)*pow(gaussianExponentB, 1.5)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      if(axisA == XAxis){
         value = -1.0*temp1*dx*dz;
      }
      else if(axisA == YAxis){
         value = -1.0*temp1*dy*dz;
      }
      else if(axisA == ZAxis){
         double temp2 = 2.0*gaussianExponentA*pow(gaussianExponentB, 0.5)
                        /(gaussianExponentA+gaussianExponentB);
         value = temp2-temp1*dz*dz;
      }
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == s){
      double temp1 = 4.0*pow(gaussianExponentA,1.5)*pow(gaussianExponentB, 2.0)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      if(axisA == XAxis){
         double temp2 = 2.0*pow(gaussianExponentA,0.5)*gaussianExponentB
                        /(gaussianExponentA+gaussianExponentB);
         value = -1.0*temp2+temp1*dx*dx;
      }
      else if(axisA == YAxis){
         value = temp1*dx*dy;
      }
      else if(axisA == ZAxis){
         value = temp1*dx*dz;
      }
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == s){
      double temp1 = 4.0*pow(gaussianExponentA,1.5)*pow(gaussianExponentB, 2.0)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      if(axisA == XAxis){
         value = temp1*dx*dy;
      }
      else if(axisA == YAxis){
         double temp2 = 2.0*pow(gaussianExponentA,0.5)*gaussianExponentB
                        /(gaussianExponentA+gaussianExponentB);
         value = -1.0*temp2+temp1*dy*dy;
      }
      else if(axisA == ZAxis){
         value = temp1*dy*dz;
      }
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == s){
      double temp1 = 4.0*pow(gaussianExponentA,1.5)*pow(gaussianExponentB, 2.0)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      if(axisA == XAxis){
         value = temp1*dx*dz;
      }
      else if(axisA == YAxis){
         value = temp1*dy*dz;
      }
      else if(axisA == ZAxis){
         double temp2 = 2.0*pow(gaussianExponentA,0.5)*gaussianExponentB
                        /(gaussianExponentA+gaussianExponentB);
         value = -1.0*temp2+temp1*dz*dz;
      }
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == py){
      double temp1 = 8.0*pow(gaussianExponentA*gaussianExponentB,2.5)
                     /pow(gaussianExponentA+gaussianExponentB,3.0);
      double temp2 = 4.0*pow(gaussianExponentA*gaussianExponentB,1.5)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      if(axisA == XAxis){
         value = -1.0*temp2*dy+temp1*dx*dx*dy;
      }
      else if(axisA == YAxis){
         value = -1.0*temp2*dx+temp1*dx*dy*dy;
      }
      else if(axisA == ZAxis){
         value = temp1*dx*dy*dz;
      }
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == px){
      double temp1 = 8.0*pow(gaussianExponentA*gaussianExponentB,2.5)
                     /pow(gaussianExponentA+gaussianExponentB,3.0);
      double temp2 = 4.0*pow(gaussianExponentA*gaussianExponentB,1.5)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      if(axisA == XAxis){
         value = -1.0*temp2*dy+temp1*dy*dx*dx;
      }
      else if(axisA == YAxis){
         value = -1.0*temp2*dx+temp1*dy*dy*dx;
      }
      else if(axisA == ZAxis){
         value = temp1*dx*dy*dz;
      }
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == pz){
      double temp1 = 8.0*pow(gaussianExponentA*gaussianExponentB,2.5)
                     /pow(gaussianExponentA+gaussianExponentB,3.0);
      double temp2 = 4.0*pow(gaussianExponentA*gaussianExponentB,1.5)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      if(axisA == XAxis){
         value = -1.0*temp2*dz+temp1*dx*dx*dz;
      }
      else if(axisA == YAxis){
         value = temp1*dx*dy*dz;
      }
      else if(axisA == ZAxis){
         value = -1.0*temp2*dx+temp1*dx*dz*dz;
      }
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == px){
      double temp1 = 8.0*pow(gaussianExponentA*gaussianExponentB,2.5)
                     /pow(gaussianExponentA+gaussianExponentB,3.0);
      double temp2 = 4.0*pow(gaussianExponentA*gaussianExponentB,1.5)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      if(axisA == XAxis){
         value = -1.0*temp2*dz+temp1*dz*dx*dx;
      }
      else if(axisA == YAxis){
         value = temp1*dx*dy*dz;
      }
      else if(axisA == ZAxis){
         value = -1.0*temp2*dx+temp1*dz*dz*dx;
      }
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == pz){
      double temp1 = 8.0*pow(gaussianExponentA*gaussianExponentB,2.5)
                     /pow(gaussianExponentA+gaussianExponentB,3.0);
      double temp2 = 4.0*pow(gaussianExponentA*gaussianExponentB,1.5)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      if(axisA == XAxis){
         value = temp1*dx*dy*dz;
      }
      else if(axisA == YAxis){
         value = -1.0*temp2*dz+temp1*dy*dy*dz;
      }
      else if(axisA == ZAxis){
         value = -1.0*temp2*dy+temp1*dy*dz*dz;
      }
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == py){
      double temp1 = 8.0*pow(gaussianExponentA*gaussianExponentB,2.5)
                     /pow(gaussianExponentA+gaussianExponentB,3.0);
      double temp2 = 4.0*pow(gaussianExponentA*gaussianExponentB,1.5)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      if(axisA == XAxis){
         value = temp1*dx*dy*dz;
      }
      else if(axisA == YAxis){
         value = -1.0*temp2*dz+temp1*dz*dy*dy;
      }
      else if(axisA == ZAxis){
         value = -1.0*temp2*dy+temp1*dz*dz*dy;
      }
   }
   else if(valenceOrbitalA == px && valenceOrbitalB == px){
      double temp1 = 8.0*pow(gaussianExponentA*gaussianExponentB,1.5)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      double temp2 = gaussianExponentA*gaussianExponentB
                     /(gaussianExponentA+gaussianExponentB); 
      if(axisA == XAxis){
         value = -1.0*temp1*dx*(1.5-temp2*dx*dx);
      }
      else if(axisA == YAxis){
         value = -1.0*temp1*dy*(0.5-temp2*dx*dx);
      }
      else if(axisA == ZAxis){
         value = -1.0*temp1*dz*(0.5-temp2*dx*dx);
      }
   }
   else if(valenceOrbitalA == py && valenceOrbitalB == py){
      double temp1 = 8.0*pow(gaussianExponentA*gaussianExponentB,1.5)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      double temp2 = gaussianExponentA*gaussianExponentB
                     /(gaussianExponentA+gaussianExponentB); 
      if(axisA == XAxis){
         value = -1.0*temp1*dx*(0.5-temp2*dy*dy);
      }
      else if(axisA == YAxis){
         value = -1.0*temp1*dy*(1.5-temp2*dy*dy);
      }
      else if(axisA == ZAxis){
         value = -1.0*temp1*dz*(0.5-temp2*dy*dy);
      }
   }
   else if(valenceOrbitalA == pz && valenceOrbitalB == pz){
      double temp1 = 8.0*pow(gaussianExponentA*gaussianExponentB,1.5)
                     /pow(gaussianExponentA+gaussianExponentB,2.0);
      double temp2 = gaussianExponentA*gaussianExponentB
                     /(gaussianExponentA+gaussianExponentB); 
      if(axisA == XAxis){
         value = -1.0*temp1*dx*(0.5-temp2*dz*dz);
      }
      else if(axisA == YAxis){
         value = -1.0*temp1*dy*(0.5-temp2*dz*dz);
      }
      else if(axisA == ZAxis){
         value = -1.0*temp1*dz*(1.5-temp2*dz*dz);
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetGaussianOverlapAOs1stDerivativeOrbitalD;
      ss << this->errorMessageAtomA;
      ss << this->errorMessageAtomType << AtomTypeStr(atomTypeA) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(valenceOrbitalA) << endl;
      ss << this->errorMessageAtomB;
      ss << this->errorMessageAtomType << AtomTypeStr(atomTypeB) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(valenceOrbitalB) << endl;
      throw MolDSException(ss.str());
   }

   double overlapSASB = this->GetGaussianOverlapAOsSASB(gaussianExponentA,
                                                        gaussianExponentB, rAB);
   value *= overlapSASB;
   return value;
}

// see J. Mol. Struc. (Theochem), 419, 19 (1997) (ref. [BFB_1997])
// we set gamma=0 always.
void Cndo2::CalcRotatingMatrix(double** rotatingMatrix, 
                               const Atom& atomA, 
                               const Atom& atomB) const{
#ifdef MOLDS_DBG
   if(rotatingMatrix==NULL){
      throw MolDSException(this->errorMessageCalcRotatingMatrixNullRotMatrix);
   }
#endif
   MallocerFreer::GetInstance()->Initialize<double>(rotatingMatrix,  OrbitalType_end, OrbitalType_end);

   double x = atomB.GetXyz()[0] - atomA.GetXyz()[0];
   double y = atomB.GetXyz()[1] - atomA.GetXyz()[1];
   double z = atomB.GetXyz()[2] - atomA.GetXyz()[2];

   EularAngle eularAngle(x, y, z);
   double alpha = eularAngle.GetAlpha();
   double beta  = eularAngle.GetBeta();

   // rotating matrix for s-function
   rotatingMatrix[s][s] = 1.0;

   // rotating matrix for p-function
   // dMatrix is (53) with gamma=0 in J. Mol. Strct. 419, 19(1997) (ref. [BFB_1997])
   rotatingMatrix[py][py] = cos(alpha);
   rotatingMatrix[py][pz] = sin(alpha)*sin(beta);
   rotatingMatrix[py][px] = sin(alpha)*cos(beta);

   rotatingMatrix[pz][py] = 0.0;
   rotatingMatrix[pz][pz] = cos(beta);
   rotatingMatrix[pz][px] = -1.0*sin(beta);

   rotatingMatrix[px][py] = -1.0*sin(alpha);
   rotatingMatrix[px][pz] = cos(alpha)*sin(beta);
   rotatingMatrix[px][px] = cos(alpha)*cos(beta);

   // rotating matrix for d-function
   // dMatrix is (37) in J. Mol. Strct. 419, 19(1997) (ref. [BFB_1997])
   double** dMatrix = NULL;
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&dMatrix, OrbitalType_end, OrbitalType_end);
      dMatrix[dzz][dzz] = 0.5*(3.0*pow(cos(beta),2.0) - 1.0);
      dMatrix[dxxyy][dxxyy] = pow(cos(0.5*beta),4.0);
      dMatrix[dzx][dzx] = (2.0*cos(beta)-1.0)*pow(cos(0.5*beta),2.0);
      dMatrix[dxxyy][dzx] = -2.0*sin(0.5*beta)*pow(cos(0.5*beta),3.0);
      dMatrix[dxxyy][dzz] = sqrt(6.0)*pow(sin(0.5*beta),2.0)*pow(cos(0.5*beta),2.0);
      dMatrix[dxxyy][dyz] = -2.0*pow(sin(0.5*beta),3.0)*pow(cos(0.5*beta),1.0);
      dMatrix[dxxyy][dxy] = pow(sin(0.5*beta),4.0);
      dMatrix[dzx][dzz] = -sqrt(6.0)*cos(beta)*cos(0.5*beta)*sin(0.5*beta);
      dMatrix[dzx][dyz] = (2.0*cos(beta)+1.0)*pow(sin(0.5*beta),2.0);

      rotatingMatrix[dxy][dxy] = cos(2.0*alpha)*            (dMatrix[dxxyy][dxxyy] - dMatrix[dxxyy][dxy]);
      rotatingMatrix[dxy][dyz] = cos(2.0*alpha)*            (-1.0*dMatrix[dxxyy][dzx] - dMatrix[dxxyy][dyz]);
      rotatingMatrix[dxy][dzz] = sqrt(2.0)*sin(2.0*alpha)*  dMatrix[dxxyy][dzz];
      rotatingMatrix[dxy][dzx] = sin(2.0*alpha)*            (-1.0*dMatrix[dxxyy][dzx] + dMatrix[dxxyy][dyz]);
      rotatingMatrix[dxy][dxxyy] = sin(2.0*alpha)*          (dMatrix[dxxyy][dxxyy] + dMatrix[dxxyy][dxy]);

      rotatingMatrix[dyz][dxy] = cos(alpha)*                (dMatrix[dxxyy][dzx] + dMatrix[dxxyy][dyz]);
      rotatingMatrix[dyz][dyz] = cos(alpha)*                (dMatrix[dzx][dzx] + dMatrix[dzx][dyz]);
      rotatingMatrix[dyz][dzz] = -1.0*sqrt(2.0)*sin(alpha)* dMatrix[dzx][dzz];
      rotatingMatrix[dyz][dzx] = sin(alpha)*                (dMatrix[dzx][dzx] - dMatrix[dzx][dyz]);
      rotatingMatrix[dyz][dxxyy] = sin(alpha)*              (dMatrix[dxxyy][dzx] - dMatrix[dxxyy][dyz]);

      rotatingMatrix[dzz][dxy] = 0.0;
      rotatingMatrix[dzz][dyz] = 0.0;
      rotatingMatrix[dzz][dzz] = dMatrix[dzz][dzz];
      rotatingMatrix[dzz][dzx] = sqrt(2.0)*dMatrix[dzx][dzz];
      rotatingMatrix[dzz][dxxyy] = sqrt(2.0)*dMatrix[dxxyy][dzz];

      rotatingMatrix[dzx][dxy] = -1.0*sin(alpha)*           (dMatrix[dxxyy][dzx] + dMatrix[dxxyy][dyz]);
      rotatingMatrix[dzx][dyz] = -1.0*sin(alpha)*           (dMatrix[dzx][dzx] + dMatrix[dzx][dyz]);
      rotatingMatrix[dzx][dzz] = -1.0*sqrt(2.0)*cos(alpha)* dMatrix[dzx][dzz];
      rotatingMatrix[dzx][dzx] = cos(alpha)*                (dMatrix[dzx][dzx] - dMatrix[dzx][dyz]);
      rotatingMatrix[dzx][dxxyy] = cos(alpha)*              (dMatrix[dxxyy][dzx] - dMatrix[dxxyy][dyz]);

      rotatingMatrix[dxxyy][dxy] = -1.0*sin(2.0*alpha)*     (dMatrix[dxxyy][dxxyy] - dMatrix[dxxyy][dxy]);
      rotatingMatrix[dxxyy][dyz] = -1.0*sin(2.0*alpha)*     (-1.0*dMatrix[dxxyy][dzx] - dMatrix[dxxyy][dyz]);
      rotatingMatrix[dxxyy][dzz] = sqrt(2.0)*cos(2.0*alpha)*dMatrix[dxxyy][dzz];
      rotatingMatrix[dxxyy][dzx] = cos(2.0*alpha)*          (-1.0*dMatrix[dxxyy][dzx] + dMatrix[dxxyy][dyz]);
      rotatingMatrix[dxxyy][dxxyy] = cos(2.0*alpha)*        (dMatrix[dxxyy][dxxyy] + dMatrix[dxxyy][dxy]);
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&dMatrix, OrbitalType_end, OrbitalType_end);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&dMatrix, OrbitalType_end, OrbitalType_end);
}

// First derivative of rotating matirx. 
// This derivative is related to a coordinate of atom A.
// This method can not calculate d-orbital yet.
// For rotating matirxi, see J. Mol. Struc. (Theochem), 419, 19 (1997) (ref. [BFB_1997])
// we set gamma=0 always.
void Cndo2::CalcRotatingMatrix1stDerivatives(double*** rotMat1stDerivatives, 
                                             const Atom& atomA, 
                                             const Atom& atomB) const{

   MallocerFreer::GetInstance()->Initialize<double>(
                                 rotMat1stDerivatives,  
                                 OrbitalType_end, 
                                 OrbitalType_end,
                                 CartesianType_end);

   double x = atomB.GetXyz()[0] - atomA.GetXyz()[0];
   double y = atomB.GetXyz()[1] - atomA.GetXyz()[1];
   double z = atomB.GetXyz()[2] - atomA.GetXyz()[2];
   double r = sqrt( pow(x,2.0) + pow(y,2.0) );
   double R = sqrt( pow(x,2.0) + pow(y,2.0) + pow(z,2.0) );

   if(r==0e0){
      return;
   }

   // for s-function
   rotMat1stDerivatives[s][s][XAxis] = 0.0;
   rotMat1stDerivatives[s][s][YAxis] = 0.0;
   rotMat1stDerivatives[s][s][ZAxis] = 0.0;

   // for p-function
   rotMat1stDerivatives[py][py][XAxis] = -1.0/r + pow(x,2.0)/pow(r,3.0);
   rotMat1stDerivatives[py][pz][XAxis] = x*y/pow(R,3.0);
   rotMat1stDerivatives[py][px][XAxis] = (1.0/(pow(r,3.0)*R) + 1.0/(pow(R,3.0)*r))*x*y*z;

   rotMat1stDerivatives[pz][py][XAxis] = 0.0;
   rotMat1stDerivatives[pz][pz][XAxis] = x*z/pow(R,3.0);
   rotMat1stDerivatives[pz][px][XAxis] = x/(r*R) - x*r/pow(R,3.0);

   rotMat1stDerivatives[px][py][XAxis] = -1.0*x*y/pow(r,3.0);
   rotMat1stDerivatives[px][pz][XAxis] = -1.0/R + x*x/pow(R,3.0); 
   rotMat1stDerivatives[px][px][XAxis] = -1.0*z/(r*R) + 
                                  (1.0/(pow(r,3.0)*R) + 1.0/(pow(R,3.0)*r))*x*x*z;

   rotMat1stDerivatives[py][py][YAxis] = x*y/pow(r,3.0);
   rotMat1stDerivatives[py][pz][YAxis] = -1.0/R + y*y/pow(R,3.0);
   rotMat1stDerivatives[py][px][YAxis] = -1.0*z/(r*R) +
                                  (1.0/(pow(r,3.0)*R) + 1.0/(pow(R,3.0)*r))*y*y*z;

   rotMat1stDerivatives[pz][py][YAxis] = 0.0;
   rotMat1stDerivatives[pz][pz][YAxis] = y*z/pow(R,3.0); 
   rotMat1stDerivatives[pz][px][YAxis] = y/(r*R) - y*r/pow(R,3.0);

   rotMat1stDerivatives[px][py][YAxis] = 1.0/r - y*y/pow(r,3.0);
   rotMat1stDerivatives[px][pz][YAxis] = x*y/pow(R,3.0);
   rotMat1stDerivatives[px][px][YAxis] = (1.0/(pow(r,3.0)*R) + 1.0/(pow(R,3.0)*r))*x*y*z;

   rotMat1stDerivatives[py][py][ZAxis] = 0.0;
   rotMat1stDerivatives[py][pz][ZAxis] = y*z/pow(R,3.0);
   rotMat1stDerivatives[py][px][ZAxis] = -1.0*y/(r*R) + y*z*z/(r*pow(R,3.0));

   rotMat1stDerivatives[pz][py][ZAxis] = 0.0;
   rotMat1stDerivatives[pz][pz][ZAxis] = -1.0/R + z*z/pow(R,3.0); 
   rotMat1stDerivatives[pz][px][ZAxis] = -1.0*z*r/pow(R,3.0);

   rotMat1stDerivatives[px][py][ZAxis] = 0.0;
   rotMat1stDerivatives[px][pz][ZAxis] = x*z/pow(R,3.0);
   rotMat1stDerivatives[px][px][ZAxis] = -1.0*x/(r*R) + x*z*z/(r*pow(R,3.0));

   // for d-function
   // ToDo: First derivative of rotating matrix for d-orbital...

}

// Second derivative of rotating matirx. 
// Both derivatives are related to a coordinate of atom A.
// This method can not calculate d-orbital yet.
// For rotating matirxi, see J. Mol. Struc. (Theochem), 419, 19 (1997) (ref. [BFB_1997])
// we set gamma=0 always.
void Cndo2::CalcRotatingMatrix2ndDerivatives(double**** rotMat2ndDerivatives, 
                                                const Atom& atomA, 
                                                const Atom& atomB) const{

   MallocerFreer::GetInstance()->Initialize<double>(
                                 rotMat2ndDerivatives,  
                                 OrbitalType_end, 
                                 OrbitalType_end,
                                 CartesianType_end,
                                 CartesianType_end);

   double x = atomB.GetXyz()[0] - atomA.GetXyz()[0];
   double y = atomB.GetXyz()[1] - atomA.GetXyz()[1];
   double z = atomB.GetXyz()[2] - atomA.GetXyz()[2];
   double r = sqrt( pow(x,2.0) + pow(y,2.0) );
   double R = sqrt( pow(x,2.0) + pow(y,2.0) + pow(z,2.0) );

   if(r==0e0){
      return;
   }

   double temp1 = 1.0/(pow(r,3.0)*R) + 1.0/(r*pow(R,3.0));
   double temp2 = 2.0*pow(r*R,-3.0) + 3.0/(pow(r,5.0)*R) + 3.0/(r*pow(R,5.0));
   double temp3 = pow(r*R,-3.0) + 3.0/(r*pow(R,5.0));

   // for s-function
   rotMat2ndDerivatives[s][s][XAxis][XAxis] = 0.0;
   rotMat2ndDerivatives[s][s][XAxis][YAxis] = 0.0;
   rotMat2ndDerivatives[s][s][XAxis][ZAxis] = 0.0;
   rotMat2ndDerivatives[s][s][YAxis][XAxis] = 0.0;
   rotMat2ndDerivatives[s][s][YAxis][YAxis] = 0.0;
   rotMat2ndDerivatives[s][s][YAxis][ZAxis] = 0.0;
   rotMat2ndDerivatives[s][s][ZAxis][XAxis] = 0.0;
   rotMat2ndDerivatives[s][s][ZAxis][YAxis] = 0.0;
   rotMat2ndDerivatives[s][s][ZAxis][ZAxis] = 0.0;

   // for p-function, xx-derivatives
   rotMat2ndDerivatives[py][py][XAxis][XAxis] = -3.0*x*pow(r,-3.0) + 3.0*pow(x,3.0)*pow(r,-5.0);
   rotMat2ndDerivatives[py][pz][XAxis][XAxis] = -1.0*y*pow(R,-3.0) + 3.0*pow(x,2.0)*y*pow(R,-5.0);
   rotMat2ndDerivatives[py][px][XAxis][XAxis] = -1.0*temp1*y*z+temp2*pow(x,2.0)*y*z;
                                              
   rotMat2ndDerivatives[pz][py][XAxis][XAxis] = 0.0;
   rotMat2ndDerivatives[pz][pz][XAxis][XAxis] = -1.0*z*pow(R,-3.0) + 3.0*pow(x,2.0)*z*pow(R,-5.0);
   rotMat2ndDerivatives[pz][px][XAxis][XAxis] = -1.0*pow(r*R,-1.0) + temp1*pow(x,2.0)
                                                   +r*pow(R,-3.0) - 3.0*pow(x,2.0)*r*pow(R,-5.0) + pow(x,2.0)*pow(r,-1.0)*pow(R,-3.0);
                                              
   rotMat2ndDerivatives[px][py][XAxis][XAxis] = y*pow(r,-3.0) - 3.0*pow(x,2.0)*y*pow(r,-5.0);
   rotMat2ndDerivatives[px][pz][XAxis][XAxis] = -3.0*x*pow(R,-3.0) + 3.0*pow(x,3.0)*pow(R,-5.0);
   rotMat2ndDerivatives[px][px][XAxis][XAxis] = -3.0*temp1*x*z+temp2*pow(x,3.0)*z;

   // for p-function, xy-derivatives
   rotMat2ndDerivatives[py][py][XAxis][YAxis] = -1.0*y*pow(r,-3.0) + 3.0*pow(x,2.0)*y*pow(r,-5.0);
   rotMat2ndDerivatives[py][pz][XAxis][YAxis] = -1.0*x*pow(R,-3.0) + 3.0*x*pow(y,2.0)*pow(R,-5.0);  
   rotMat2ndDerivatives[py][px][XAxis][YAxis] = -1.0*temp1*x*z+temp2*x*pow(y,2.0)*z;
                                              
   rotMat2ndDerivatives[pz][py][XAxis][YAxis] = 0.0;
   rotMat2ndDerivatives[pz][pz][XAxis][YAxis] = 3.0*x*y*z*pow(R,-5.0);
   rotMat2ndDerivatives[pz][px][XAxis][YAxis] = temp1*x*y + x*y*pow(r,-1.0)*pow(R,-3.0) - 3.0*x*y*r*pow(R,-5.0);
                                              
   rotMat2ndDerivatives[px][py][XAxis][YAxis] = x*pow(r,-3.0) - 3.0*x*pow(y,2.0)*pow(r,-5.0);
   rotMat2ndDerivatives[px][pz][XAxis][YAxis] = rotMat2ndDerivatives[py][pz][XAxis][XAxis];
   rotMat2ndDerivatives[px][px][XAxis][YAxis] = rotMat2ndDerivatives[py][px][XAxis][XAxis];

   // for p-function, yx-derivatives
   for(int i=py; i<=px; i++){
      for(int j=py; j<=px; j++){
         rotMat2ndDerivatives[i][j][YAxis][XAxis] = rotMat2ndDerivatives[i][j][XAxis][YAxis];
      }
   }
   // for p-function, xz-derivatives
   rotMat2ndDerivatives[py][py][XAxis][ZAxis] = 0.0;
   rotMat2ndDerivatives[py][pz][XAxis][ZAxis] = rotMat2ndDerivatives[pz][pz][XAxis][YAxis];
   rotMat2ndDerivatives[py][px][XAxis][ZAxis] = -1.0*temp1*x*y +temp3*x*y*pow(z,2.0);
                                              
   rotMat2ndDerivatives[pz][py][XAxis][ZAxis] = 0.0;
   rotMat2ndDerivatives[pz][pz][XAxis][ZAxis] = -1.0*x*pow(R,-3.0) + 3.0*x*pow(z,2.0)*pow(R,-5.0); 
   rotMat2ndDerivatives[pz][px][XAxis][ZAxis] = x*z*pow(r,-1.0)*pow(R,-3.0) - 3.0*x*z*r*pow(R,-5.0);
                                              
   rotMat2ndDerivatives[px][py][XAxis][ZAxis] = 0.0;
   rotMat2ndDerivatives[px][pz][XAxis][ZAxis] = rotMat2ndDerivatives[pz][pz][XAxis][XAxis];
   rotMat2ndDerivatives[px][px][XAxis][ZAxis] = pow(r*R,-1.0) - pow(z,2.0)*pow(r,-1.0)*pow(R,-3.0)
                                                   -1.0*temp1*pow(x,2.0)+temp3*pow(x*z,2.0);


   // for p-function, zx-derivatives
   for(int i=py; i<=px; i++){
      for(int j=py; j<=px; j++){
         rotMat2ndDerivatives[i][j][ZAxis][XAxis] = rotMat2ndDerivatives[i][j][XAxis][ZAxis];
      }
   }

   // for p-function, yy-derivatives
   rotMat2ndDerivatives[py][py][YAxis][YAxis] = -1.0*x*pow(r,-3.0) + 3.0*x*pow(y,2.0)*pow(r,-5.0); 
   rotMat2ndDerivatives[py][pz][YAxis][YAxis] = -3.0*y*pow(R,-3.0) + 3.0*pow(y,3.0)*pow(R,-5.0);
   rotMat2ndDerivatives[py][px][YAxis][YAxis] = -3.0*temp1*y*z+temp2*pow(y,3.0)*z;
                                              
   rotMat2ndDerivatives[pz][py][YAxis][YAxis] = 0.0;
   rotMat2ndDerivatives[pz][pz][YAxis][YAxis] = -1.0*z*pow(R,-3.0) + 3.0*pow(y,2.0)*z*pow(R,-5.0);
   rotMat2ndDerivatives[pz][px][YAxis][YAxis] = -1.0*pow(r*R,-1.0) + temp1*pow(y,2.0)
                                                   +r*pow(R,-3.0) - 3.0*pow(y,2.0)*r*pow(R,-5.0) + pow(y,2.0)*pow(r,-1.0)*pow(R,-3.0);
                                              
   rotMat2ndDerivatives[px][py][YAxis][YAxis] = 3.0*y*pow(r,-3.0) - 3.0*pow(y,3.0)*pow(r,-5.0);
   rotMat2ndDerivatives[px][pz][YAxis][YAxis] = rotMat2ndDerivatives[py][pz][XAxis][YAxis];
   rotMat2ndDerivatives[px][px][YAxis][YAxis] = -1.0*temp1*x*z+temp2*x*pow(y,2.0)*z;
               
   // for p-function, yz-derivatives
   rotMat2ndDerivatives[py][py][YAxis][ZAxis] = 0.0;
   rotMat2ndDerivatives[py][pz][YAxis][ZAxis] = rotMat2ndDerivatives[pz][pz][YAxis][YAxis];
   rotMat2ndDerivatives[py][px][YAxis][ZAxis] = pow(r*R,-1.0) - pow(z,2.0)*pow(r,-1.0)*pow(R,-3.0)
                                                   -1.0*temp1*pow(y,2.0)+temp3*pow(y*z,2.0);
                                              
   rotMat2ndDerivatives[pz][py][YAxis][ZAxis] = 0.0;
   rotMat2ndDerivatives[pz][pz][YAxis][ZAxis] = -1.0*y*pow(R,-3.0) + 3.0*y*pow(z,2.0)*pow(R,-5.0);
   rotMat2ndDerivatives[pz][px][YAxis][ZAxis] = y*z*pow(r,-1.0)*pow(R,-3.0) - 3.0*y*z*r*pow(R,-5.0);
                                              
   rotMat2ndDerivatives[px][py][YAxis][ZAxis] = 0.0;
   rotMat2ndDerivatives[px][pz][YAxis][ZAxis] = rotMat2ndDerivatives[pz][pz][XAxis][YAxis];
   rotMat2ndDerivatives[px][px][YAxis][ZAxis] = -1.0*temp1*x*y+temp3*x*y*pow(z,2.0);
                                          
               
   // for p-function, zy-derivatives
   for(int i=py; i<=px; i++){
      for(int j=py; j<=px; j++){
         rotMat2ndDerivatives[i][j][ZAxis][YAxis] = rotMat2ndDerivatives[i][j][YAxis][ZAxis];
      }
   }

   // for p-function, zz-derivatives
   rotMat2ndDerivatives[py][py][ZAxis][ZAxis] = 0.0;
   rotMat2ndDerivatives[py][pz][ZAxis][ZAxis] = rotMat2ndDerivatives[pz][pz][YAxis][ZAxis];
   rotMat2ndDerivatives[py][px][ZAxis][ZAxis] = -3.0*y*z*pow(r,-1.0)*pow(R,-3.0) + 3.0*y*pow(z,3.0)*pow(r,-1.0)*pow(R,-5.0);
                                              
   rotMat2ndDerivatives[pz][py][ZAxis][ZAxis] = 0.0;
   rotMat2ndDerivatives[pz][pz][ZAxis][ZAxis] = -3.0*z*pow(R,-3.0) + 3.0*pow(z,3.0)*pow(R,-5.0); 
   rotMat2ndDerivatives[pz][px][ZAxis][ZAxis] = -3.0*pow(z,2.0)*r*pow(R,-5.0) + r*pow(R,-3.0);
                                              
   rotMat2ndDerivatives[px][py][ZAxis][ZAxis] = 0.0;
   rotMat2ndDerivatives[px][pz][ZAxis][ZAxis] = rotMat2ndDerivatives[pz][pz][XAxis][ZAxis];
   rotMat2ndDerivatives[px][px][ZAxis][ZAxis] = -3.0*x*z*pow(r,-1.0)*pow(R,-3.0) + 3.0*x*pow(z,3.0)*pow(r,-1.0)*pow(R,-5.0);

   // for d-function
   // ToDo: Second derivative of rotating matrix for d-orbital...
}

// see (B.40) in J. A. Pople book.
void Cndo2::CalcDiatomicOverlapAOsInDiatomicFrame(double** diatomicOverlapAOs, 
                                                  const Atom& atomA, 
                                                  const Atom& atomB) const{
#ifdef MOLDS_DBG
   if(diatomicOverlapAOs==NULL){
      throw MolDSException(this->errorMessageCalDiaOverlapAOsDiaFrameNullMatrix);
   }
#endif
   int na = atomA.GetValenceShellType() + 1;
   int nb = atomB.GetValenceShellType() + 1;
   int m = 0;
   double alpha = 0.0;
   double beta = 0.0;
   double pre = 0.0;
   double reducedOverlapAOs = 0.0;
   double orbitalExponentA = 0.0;
   double orbitalExponentB = 0.0;
   double rAB = 0.0; // Inter nuclear distance between aton A and B.

   MallocerFreer::GetInstance()->Initialize<double>(diatomicOverlapAOs, OrbitalType_end, OrbitalType_end);
   rAB = this->molecule->GetDistanceAtoms(atomA, atomB);

   for(int a=0; a<atomA.GetValenceSize(); a++){
      OrbitalType valenceOrbitalA = atomA.GetValence(a);
      RealSphericalHarmonicsIndex realShpericalHarmonicsA(valenceOrbitalA);
      orbitalExponentA = atomA.GetOrbitalExponent(
                               atomA.GetValenceShellType(), 
                               valenceOrbitalA, 
                               this->theory);

      for(int b=0; b<atomB.GetValenceSize(); b++){
         OrbitalType valenceOrbitalB = atomB.GetValence(b);
         RealSphericalHarmonicsIndex realShpericalHarmonicsB(valenceOrbitalB);
         orbitalExponentB = atomB.GetOrbitalExponent(
                                  atomB.GetValenceShellType(), 
                                  valenceOrbitalB, 
                                  this->theory);

         if(realShpericalHarmonicsA.GetM() == realShpericalHarmonicsB.GetM()){
            m = abs(realShpericalHarmonicsA.GetM());
            alpha = orbitalExponentA * rAB;
            beta =  orbitalExponentB * rAB;

            reducedOverlapAOs = this->GetReducedOverlapAOs(na, realShpericalHarmonicsA.GetL(), m,
                                                           nb, realShpericalHarmonicsB.GetL(), alpha, beta);


            pre =  pow(2.0*orbitalExponentA, na+0.5);
            pre *= pow(2.0*orbitalExponentB, nb+0.5);
            double factorials = Factorial(2*na)*Factorial(2*nb);
            pre /= sqrt(factorials);
            pre *= pow(rAB/2.0, na+nb+1.0);

            diatomicOverlapAOs[valenceOrbitalA][valenceOrbitalB] = pre*reducedOverlapAOs;
         }
         
      }
   }

   /*
   for(int i=0;i<OrbitalType_end;i++){
      for(int j=0;j<OrbitalType_end;j++){
         this->OutputLog(boost::format("diatomicOverlapAOs[%d][%d]=%lf\n") % i % j % diatomicOverlapAOs[i][j]);
      }
   }
   */
}

// First derivative of (B.40) in J. A. Pople book.
void Cndo2::CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(double** diatomicOverlapAOsDeri, 
                                                               const Atom& atomA, 
                                                               const Atom& atomB) const{

   int na = atomA.GetValenceShellType() + 1;
   int nb = atomB.GetValenceShellType() + 1;
   int m = 0;
   double alpha = 0.0;
   double beta = 0.0;
   double pre = 0.0;
   double reducedOverlapAOs = 0.0;
   double reducedOverlapAOs1stDerivAlpha = 0.0;
   double reducedOverlapAOs1stDerivBeta = 0.0;
   double orbitalExponentA = 0.0;
   double orbitalExponentB = 0.0;
   double temp1=0.0;
   double temp2=0.0;

   MallocerFreer::GetInstance()->Initialize<double>(diatomicOverlapAOsDeri, 
                                                    OrbitalType_end, 
                                                    OrbitalType_end);
   double R = this->molecule->GetDistanceAtoms(atomA, atomB);

   for(int a=0; a<atomA.GetValenceSize(); a++){
      OrbitalType valenceOrbitalA = atomA.GetValence(a);
      RealSphericalHarmonicsIndex realShpericalHarmonicsA(valenceOrbitalA);
      orbitalExponentA = atomA.GetOrbitalExponent(
                               atomA.GetValenceShellType(), 
                               valenceOrbitalA,
                               this->theory);

      for(int b=0; b<atomB.GetValenceSize(); b++){
         OrbitalType valenceOrbitalB = atomB.GetValence(b);
         RealSphericalHarmonicsIndex realShpericalHarmonicsB(valenceOrbitalB);
         orbitalExponentB = atomB.GetOrbitalExponent(
                                  atomB.GetValenceShellType(), 
                                  valenceOrbitalB,
                                  this->theory);

         if(realShpericalHarmonicsA.GetM() == realShpericalHarmonicsB.GetM()){
            m = abs(realShpericalHarmonicsA.GetM());
            alpha = orbitalExponentA * R;
            beta =  orbitalExponentB * R;

            reducedOverlapAOs = this->GetReducedOverlapAOs(na, realShpericalHarmonicsA.GetL(), m,
                                                           nb, realShpericalHarmonicsB.GetL(), alpha, beta);
            reducedOverlapAOs1stDerivAlpha = this->GetReducedOverlapAOs1stDerivativeAlpha(
                                                   na, 
                                                   realShpericalHarmonicsA.GetL(), 
                                                   m,
                                                   nb, 
                                                   realShpericalHarmonicsB.GetL(), 
                                                   alpha, 
                                                   beta);
            reducedOverlapAOs1stDerivBeta  = this->GetReducedOverlapAOs1stDerivativeBeta(
                                                   na, 
                                                   realShpericalHarmonicsA.GetL(), 
                                                   m,
                                                   nb, 
                                                   realShpericalHarmonicsB.GetL(), 
                                                   alpha, 
                                                   beta);

            temp1 = static_cast<double>(na+nb+1)*pow(R,na+nb)*reducedOverlapAOs;
            temp2 = pow(R,na+nb+1)*(orbitalExponentA*reducedOverlapAOs1stDerivAlpha
                                   +orbitalExponentB*reducedOverlapAOs1stDerivBeta);

            pre =  pow(2.0*orbitalExponentA, na+0.5);
            pre *= pow(2.0*orbitalExponentB, nb+0.5);
            double factorials = Factorial(2*na)*Factorial(2*nb);
            pre /= sqrt(factorials);
            pre /= pow(2.0, na+nb+1.0);

            diatomicOverlapAOsDeri[valenceOrbitalA][valenceOrbitalB] = pre*(temp1+temp2);
         }
         
      }
   }

   /*
   for(int i=0;i<OrbitalType_end;i++){
      for(int j=0;j<OrbitalType_end;j++){
         this->OutputLog(boost::format("diatomicOverlapAOs[%d][%d]=%lf\n") % i % j % diatomicOverlapAOs[i][j]);
      }
   }
   */

}

// Second derivative of (B.40) in J. A. Pople book.
void Cndo2::CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(double** diatomicOverlapAOs2ndDeri, 
                                                               const Atom& atomA, 
                                                               const Atom& atomB) const{

   int na = atomA.GetValenceShellType() + 1;
   int nb = atomB.GetValenceShellType() + 1;
   int m = 0;
   double alpha = 0.0;
   double beta = 0.0;
   double pre = 0.0;
   double reducedOverlapAOs = 0.0;
   double reducedOverlapAOs1stDerivAlpha = 0.0;
   double reducedOverlapAOs1stDerivBeta = 0.0;
   double reducedOverlapAOs2ndDerivAlpha = 0.0;
   double reducedOverlapAOs2ndDerivBeta = 0.0;
   double reducedOverlapAOs2ndDerivAlphaBeta = 0.0;
   double orbitalExponentA = 0.0;
   double orbitalExponentB = 0.0;
   double temp1=0.0;
   double temp2=0.0;
   double temp3=0.0;

   MallocerFreer::GetInstance()->Initialize<double>(diatomicOverlapAOs2ndDeri, 
                                                    OrbitalType_end, 
                                                    OrbitalType_end);
   double R = this->molecule->GetDistanceAtoms(atomA, atomB);

   for(int a=0; a<atomA.GetValenceSize(); a++){
      OrbitalType valenceOrbitalA = atomA.GetValence(a);
      RealSphericalHarmonicsIndex realShpericalHarmonicsA(valenceOrbitalA);
      orbitalExponentA = atomA.GetOrbitalExponent(atomA.GetValenceShellType(),
                                                  valenceOrbitalA,
                                                  this->theory);

      for(int b=0; b<atomB.GetValenceSize(); b++){
         OrbitalType valenceOrbitalB = atomB.GetValence(b);
         RealSphericalHarmonicsIndex realShpericalHarmonicsB(valenceOrbitalB);
         orbitalExponentB = atomB.GetOrbitalExponent(atomB.GetValenceShellType(),
                                                     valenceOrbitalB,
                                                     this->theory);

         if(realShpericalHarmonicsA.GetM() == realShpericalHarmonicsB.GetM()){
            m = abs(realShpericalHarmonicsA.GetM());
            alpha = orbitalExponentA * R;
            beta =  orbitalExponentB * R;

            reducedOverlapAOs = this->GetReducedOverlapAOs(na,
                                                           realShpericalHarmonicsA.GetL(),
                                                           m,
                                                           nb,
                                                           realShpericalHarmonicsB.GetL(),
                                                           alpha,
                                                           beta);
            reducedOverlapAOs1stDerivAlpha
               = this->GetReducedOverlapAOs1stDerivativeAlpha(na,
                                                              realShpericalHarmonicsA.GetL(),
                                                              m,
                                                              nb,
                                                              realShpericalHarmonicsB.GetL(),
                                                              alpha,
                                                              beta);
            reducedOverlapAOs1stDerivBeta
               = this->GetReducedOverlapAOs1stDerivativeBeta(na,
                                                             realShpericalHarmonicsA.GetL(),
                                                             m,
                                                             nb,
                                                             realShpericalHarmonicsB.GetL(),
                                                             alpha,
                                                             beta);
            reducedOverlapAOs2ndDerivAlpha
               = this->GetReducedOverlapAOs2ndDerivativeAlpha(na,
                                                              realShpericalHarmonicsA.GetL(),
                                                              m,
                                                              nb,
                                                              realShpericalHarmonicsB.GetL(),
                                                              alpha,
                                                              beta);
            reducedOverlapAOs2ndDerivBeta
               = this->GetReducedOverlapAOs2ndDerivativeBeta(na,
                                                             realShpericalHarmonicsA.GetL(),
                                                             m,
                                                             nb,
                                                             realShpericalHarmonicsB.GetL(),
                                                             alpha,
                                                             beta);
            reducedOverlapAOs2ndDerivAlphaBeta
               = this->GetReducedOverlapAOs2ndDerivativeAlphaBeta(na,
                                                                  realShpericalHarmonicsA.GetL(),
                                                                  m,
                                                                  nb,
                                                                  realShpericalHarmonicsB.GetL(),
                                                                  alpha,
                                                                  beta);

            temp1 = static_cast<double>(na+nb+1)
                   *static_cast<double>(na+nb)
                   *pow(R,na+nb-1)*reducedOverlapAOs;
            temp2 = 2.0*static_cast<double>(na+nb+1)*pow(R,na+nb)
                       *(orbitalExponentA*reducedOverlapAOs1stDerivAlpha
                        +orbitalExponentB*reducedOverlapAOs1stDerivBeta);
            temp3 = pow(R,na+nb+1)
                   *(pow(orbitalExponentA,2.0)*reducedOverlapAOs2ndDerivAlpha
                    +pow(orbitalExponentB,2.0)*reducedOverlapAOs2ndDerivBeta
                    +2.0*orbitalExponentA*orbitalExponentB*reducedOverlapAOs2ndDerivAlphaBeta);

            pre =  pow(2.0*orbitalExponentA, na+0.5);
            pre *= pow(2.0*orbitalExponentB, nb+0.5);
            double factorials = Factorial(2*na)*Factorial(2*nb);
            pre /= sqrt(factorials);
            pre /= pow(2.0, na+nb+1.0);

            diatomicOverlapAOs2ndDeri[valenceOrbitalA][valenceOrbitalB] = pre*(temp1+temp2+temp3);
         }
         
      }
   }

   /*
   for(int i=0;i<OrbitalType_end;i++){
      for(int j=0;j<OrbitalType_end;j++){
         this->OutputLog(boost::format("diatomicOverlapAOs[%d][%d]=%lf\n") % i % j % diatomicOverlapAOs[i][j]);
      }
   }
   */


}

// see (B.63) in Pople book.
void Cndo2::RotateDiatmicOverlapAOsToSpaceFrame(double** diatomicOverlapAOs, 
                                                double const* const* rotatingMatrix) const{
#ifdef MOLDS_DBG
   if(diatomicOverlapAOs==NULL){
      throw MolDSException(this->errorMessageRotDiaOverlapAOsToSpaceFrameNullDiaMatrix);
   }
   if(rotatingMatrix==NULL){
      throw MolDSException(this->errorMessageRotDiaOverlapAOsToSpaceFrameNullRotMatrix);
   }
#endif
   double** oldDiatomicOverlapAOs = NULL;
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&oldDiatomicOverlapAOs, OrbitalType_end, OrbitalType_end);
      for(int i=0; i<OrbitalType_end; i++){
         for(int j=0; j<OrbitalType_end; j++){
            oldDiatomicOverlapAOs[i][j] = diatomicOverlapAOs[i][j];
         }
      }
      // rotate
      bool isColumnMajorRotatingMatrix = false;
      bool isColumnMajorOldDiatomicOverlap        = false;
      double alpha = 1.0;
      double beta  = 0.0;
      MolDS_wrappers::Blas::GetInstance()->Dgemmm(isColumnMajorRotatingMatrix,
                                                  isColumnMajorOldDiatomicOverlap,
                                                  !isColumnMajorRotatingMatrix,
                                                  OrbitalType_end, OrbitalType_end, OrbitalType_end, OrbitalType_end, 
                                                  alpha,
                                                  rotatingMatrix,
                                                  oldDiatomicOverlapAOs,
                                                  rotatingMatrix,
                                                  beta, 
                                                  diatomicOverlapAOs);
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&oldDiatomicOverlapAOs, OrbitalType_end, OrbitalType_end);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&oldDiatomicOverlapAOs, OrbitalType_end, OrbitalType_end);

   /*
   for(int i=0;i<OrbitalType_end;i++){
      for(int j=0;j<OrbitalType_end;j++){
         this->OutputLog(boost::format("rotated diatomicOverlapAOs[%d][%d]=%lf\n") % i % j % diatomicOverlapAOs[i][j]);
         this->OutputLog(boost::format("rotating[%d][%d]=%lf\n") % i % j % rotating[i][j]);
      }
   }
   */

}

void Cndo2::SetOverlapAOsElement(double** overlapAOs, 
                                 double const* const* diatomicOverlapAOs, 
                                 const Atom& atomA, 
                                 const Atom& atomB,
                                 bool isSymmetricOverlapAOs) const{
#ifdef MOLDS_DBG
   if(diatomicOverlapAOs==NULL){
      throw MolDSException(this->errorMessageSetOverlapAOsElementNullDiaMatrix);
   }
#endif

   int firstAOIndexAtomA = atomA.GetFirstAOIndex();
   int firstAOIndexAtomB = atomB.GetFirstAOIndex();
   OrbitalType orbitalA;
   OrbitalType orbitalB;
   int mu=0;
   int nu=0;

   for(int i=0; i<atomA.GetValenceSize(); i++){
      orbitalA = atomA.GetValence(i);
      for(int j=0; j<atomB.GetValenceSize(); j++){
         orbitalB = atomB.GetValence(j);
         mu = firstAOIndexAtomA + i;      
         nu = firstAOIndexAtomB + j;      
         overlapAOs[mu][nu] = diatomicOverlapAOs[orbitalA][orbitalB];
         if(isSymmetricOverlapAOs){
            overlapAOs[nu][mu] = diatomicOverlapAOs[orbitalA][orbitalB];
         }
      }
   }

}

// see (B.24) in J. A. Pople book.
double Cndo2::GetReducedOverlapAOs(int na, int la, int m, int nb, int lb, double alpha, double beta) const{
   double value = 0.0;
   double temp = 0.0;
   int I = 2*ShellType_end+1;
   int J = 2*ShellType_end+1;

   for(int i=0; i<I; i++){
      for(int j=0; j<J; j++){
         temp = Cndo2::ReducedOverlapAOsParameters::Y[na][nb][la][lb][m][i][j];
         if(0e0<fabs(temp)){
            temp *= this->GetAuxiliaryA(i, 0.5*(alpha+beta));
            temp *= this->GetAuxiliaryB(j, 0.5*(alpha-beta));
            value += temp;
         }
      }
   }
   value *= this->GetAuxiliaryD(la, lb, m);

   return value;
}

// see (B.30) in J. A. Pople book.
double Cndo2::GetReducedOverlapAOs(int na, int nb, double alpha, double beta) const{
   double value = 0.0;
   double temp = 0.0;

   for(int k=0; k<=na+nb; k++){
      temp = Cndo2::ReducedOverlapAOsParameters::Z[na][nb][k];
      if(0e0<fabs(temp)){
         temp *= this->GetAuxiliaryA(k, 0.5*(alpha+beta));
         temp *= this->GetAuxiliaryB(na+nb-k, 0.5*(alpha-beta));
         value += temp;
      }
   }
   value *= 0.5;
   return value;
}

// First derivative of (B.24) in J. A. Pople book.
// This derivative is carried out by alpha.
double Cndo2::GetReducedOverlapAOs1stDerivativeAlpha(int na, 
                                                     int la, 
                                                     int m, 
                                                     int nb, 
                                                     int lb, 
                                                     double alpha, 
                                                     double beta) const{
   double value = 0.0;
   double temp1 = 0.0;
   double temp2 = 0.0;
   int I = 2*ShellType_end+1;
   int J = 2*ShellType_end+1;

   for(int i=0; i<I; i++){
      for(int j=0; j<J; j++){
         if(0e0<fabs(Cndo2::ReducedOverlapAOsParameters::Y[na][nb][la][lb][m][i][j])){
            temp1 = this->GetAuxiliaryA1stDerivative(i, 0.5*(alpha+beta))
                   *this->GetAuxiliaryB(j, 0.5*(alpha-beta));
            temp2 = this->GetAuxiliaryA(i, 0.5*(alpha+beta))
                   *this->GetAuxiliaryB1stDerivative(j, 0.5*(alpha-beta));
            value += Cndo2::ReducedOverlapAOsParameters::Y[na][nb][la][lb][m][i][j]*(temp1 + temp2);
         }
      }
   }
   value *= 0.5*this->GetAuxiliaryD(la, lb, m);

   return value;
}

// First derivative of (B.24) in J. A. Pople book.
// This derivative is carried out by Beta.
double Cndo2::GetReducedOverlapAOs1stDerivativeBeta(int na, 
                                                    int la, 
                                                    int m, 
                                                    int nb, 
                                                    int lb, 
                                                    double alpha, 
                                                    double beta) const{
   double value = 0.0;
   double temp1 = 0.0;
   double temp2 = 0.0;
   int I = 2*ShellType_end+1;
   int J = 2*ShellType_end+1;

   for(int i=0; i<I; i++){
      for(int j=0; j<J; j++){
         if(0e0<fabs(Cndo2::ReducedOverlapAOsParameters::Y[na][nb][la][lb][m][i][j])){
            temp1 = this->GetAuxiliaryA1stDerivative(i, 0.5*(alpha+beta))
                   *this->GetAuxiliaryB(j, 0.5*(alpha-beta));
            temp2 = this->GetAuxiliaryA(i, 0.5*(alpha+beta))
                   *this->GetAuxiliaryB1stDerivative(j, 0.5*(alpha-beta));
            value += Cndo2::ReducedOverlapAOsParameters::Y[na][nb][la][lb][m][i][j]*(temp1 - temp2);
         }
      }
   }
   value *= 0.5*this->GetAuxiliaryD(la, lb, m);

   return value;
}

// Second derivative of (B.24) in J. A. Pople book.
// This derivative is carried out by alpha twice.
double Cndo2::GetReducedOverlapAOs2ndDerivativeAlpha(int na, 
                                                     int la, 
                                                     int m, 
                                                     int nb, 
                                                     int lb, 
                                                     double alpha, 
                                                     double beta) const{
   double value = 0.0;
   double temp1 = 0.0;
   double temp2 = 0.0;
   double temp3 = 0.0;
   int I = 2*ShellType_end+1;
   int J = 2*ShellType_end+1;

   for(int i=0; i<I; i++){
      for(int j=0; j<J; j++){
         if(0e0<fabs(Cndo2::ReducedOverlapAOsParameters::Y[na][nb][la][lb][m][i][j])){
            temp1 = this->GetAuxiliaryA2ndDerivative(i, 0.5*(alpha+beta))
                   *this->GetAuxiliaryB(j, 0.5*(alpha-beta));
            temp2 = this->GetAuxiliaryA(i, 0.5*(alpha+beta))
                   *this->GetAuxiliaryB2ndDerivative(j, 0.5*(alpha-beta));
            temp3 = this->GetAuxiliaryA1stDerivative(i, 0.5*(alpha+beta))
                   *this->GetAuxiliaryB1stDerivative(j, 0.5*(alpha-beta));
            value += Cndo2::ReducedOverlapAOsParameters::Y[na][nb][la][lb][m][i][j]*(temp1 + temp2 + 2.0*temp3);
         }
      }
   }
   value *= 0.25*this->GetAuxiliaryD(la, lb, m);

   return value;
}

// Second derivative of (B.24) in J. A. Pople book.
// This derivative is carried out by beta twice.
double Cndo2::GetReducedOverlapAOs2ndDerivativeBeta(int na, 
                                                    int la, 
                                                    int m, 
                                                    int nb, 
                                                    int lb, 
                                                    double alpha, 
                                                    double beta) const{
   double value = 0.0;
   double temp1 = 0.0;
   double temp2 = 0.0;
   double temp3 = 0.0;
   int I = 2*ShellType_end+1;
   int J = 2*ShellType_end+1;

   for(int i=0; i<I; i++){
      for(int j=0; j<J; j++){
         if(0e0<fabs(Cndo2::ReducedOverlapAOsParameters::Y[na][nb][la][lb][m][i][j])){
            temp1 = this->GetAuxiliaryA2ndDerivative(i, 0.5*(alpha+beta))
                   *this->GetAuxiliaryB(j, 0.5*(alpha-beta));
            temp2 = this->GetAuxiliaryA(i, 0.5*(alpha+beta))
                   *this->GetAuxiliaryB2ndDerivative(j, 0.5*(alpha-beta));
            temp3 = this->GetAuxiliaryA1stDerivative(i, 0.5*(alpha+beta))
                   *this->GetAuxiliaryB1stDerivative(j, 0.5*(alpha-beta));
            value += Cndo2::ReducedOverlapAOsParameters::Y[na][nb][la][lb][m][i][j]*(temp1 + temp2 - 2.0*temp3);
         }
      }
   }
   value *= 0.25*this->GetAuxiliaryD(la, lb, m);

   return value;
}

// Second derivative of (B.24) in J. A. Pople book.
// This derivative is carried out by alpha and beta.
double Cndo2::GetReducedOverlapAOs2ndDerivativeAlphaBeta(int na, 
                                                         int la, 
                                                         int m, 
                                                         int nb, 
                                                         int lb, 
                                                         double alpha, 
                                                         double beta) const{
   double value = 0.0;
   double temp1 = 0.0;
   double temp2 = 0.0;
   int I = 2*ShellType_end+1;
   int J = 2*ShellType_end+1;

   for(int i=0; i<I; i++){
      for(int j=0; j<J; j++){
         if(0e0<fabs(Cndo2::ReducedOverlapAOsParameters::Y[na][nb][la][lb][m][i][j])){
            temp1 = this->GetAuxiliaryA2ndDerivative(i, 0.5*(alpha+beta))
                   *this->GetAuxiliaryB(j, 0.5*(alpha-beta));
            temp2 = this->GetAuxiliaryA(i, 0.5*(alpha+beta))
                   *this->GetAuxiliaryB2ndDerivative(j, 0.5*(alpha-beta));
            value += Cndo2::ReducedOverlapAOsParameters::Y[na][nb][la][lb][m][i][j]*(temp1 - temp2);
         }
      }
   }
   value *= 0.25*this->GetAuxiliaryD(la, lb, m);

   return value;
}

// see (B.22) in J. A. Pople book.
double Cndo2::GetAuxiliaryA(int k, double rho) const{
   double value = 0.0;
   double temp = 0.0;

   value = exp(-1.0*rho)*static_cast<double>(Factorial(k));
   for(int mu=1; mu<=k+1; mu++){
      temp += pow(rho,-1.0*mu)/static_cast<double>(Factorial(k-mu+1));
   }
   value *= temp;

   return value;
}

// First derivative of (B.22) in J. A. Pople book.
double Cndo2::GetAuxiliaryA1stDerivative(int k, double rho) const{
   return -1.0*this->GetAuxiliaryA(k+1, rho);
}

// Second derivative of (B.22) in J. A. Pople book.
double Cndo2::GetAuxiliaryA2ndDerivative(int k, double rho) const{
   return this->GetAuxiliaryA(k+2, rho);
}

// see (B.23) in J. A. Pople book.
double Cndo2::GetAuxiliaryB(int k, double rho) const{
   double value = 0.0;
   double pre1 = 0.0;
   double pre2 = 0.0;
   double temp1 = 0.0;
   double temp2 = 0.0;

   if(fabs(rho)>0){
      pre1 = -1.0*exp(-1.0*rho);
      pre2 = -1.0*exp(rho);
      
      for(int mu=1; mu<=k+1; mu++){
         temp1 += pow(rho,-1.0*mu)  *static_cast<double>(Factorial(k)/Factorial(k-mu+1)) ;
         temp2 += pow(rho,-1.0*mu)  *static_cast<double>(Factorial(k)/Factorial(k-mu+1)) *pow(-1.0,k-mu);
      }
      value = pre1*temp1 + pre2*temp2;
   }
   else{
      if(k%2 == 0){
         value = 2.0/(1.0+static_cast<double>(k));
      }
      else{
         value = 0;
      }
   }

   return value;
}

// First derivative of (B.23) in J. A. Pople book.
double Cndo2::GetAuxiliaryB1stDerivative(int k, double rho) const{
   return -1.0*this->GetAuxiliaryB(k+1, rho);
}

// Second derivative of (B.23) in J. A. Pople book.
double Cndo2::GetAuxiliaryB2ndDerivative(int k, double rho) const{
   return this->GetAuxiliaryB(k+2, rho);
}

// see (B.16) in J. A. Pople book.
double Cndo2::GetAuxiliaryD(int la, int lb, int m) const{
   string errorMessageAuxiliaryDNegativeM = "Error in cndo::Cndo2::GetAuxiliaryD: m<0\n";
   double value = 0.0;

   if(m<0){
      stringstream ss;
      ss << errorMessageAuxiliaryDNegativeM;
      throw MolDSException(ss.str());
   }

   double pre = pow(Factorial(m+1)/8.0, 2.0);
   double termA = ( (2.0*la+1.0)*Factorial(la-m) ) / ( 2.0*Factorial(la+m) );
   double termB = ( (2.0*lb+1.0)*Factorial(lb-m) ) / ( 2.0*Factorial(lb+m) );
   value = pre*sqrt(termA)*sqrt(termB);
   //this->OutputLog(boost::format("pre=%lf, termA=%lf, termB=%lf\n") % pre % termA % termB);
   
   return value;
}
}



