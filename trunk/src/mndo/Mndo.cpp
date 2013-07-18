//************************************************************************//
// Copyright (C) 2011-2013 Mikiya Fujii                                   //
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
#include"../base/ElectronicStructure.h"
#include"../cndo/Cndo2.h"
#include"../zindo/ZindoS.h"
#include"Mndo.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_mndo{

/***
 *  Main References for MNDO are [DT_1977, DT_1977-2, DT_1977-3]
 */
Mndo::Mndo() : MolDS_zindo::ZindoS(){
   // protedted variables and methods
   this->theory = MNDO;
   this->SetMessages();
   this->SetEnableAtomTypes();
   // private variables
   this->heatsFormation = 0.0;
   //this->OutputLog("Mndo created\n");
}

Mndo::~Mndo(){
   MallocerFreer::GetInstance()->Free<double>(
                                 &this->twoElecTwoCore, 
                                 this->molecule->GetNumberAtoms(),
                                 this->molecule->GetNumberAtoms(),
                                 dxy,
                                 dxy,
                                 dxy,
                                 dxy);
   MallocerFreer::GetInstance()->Free<double>(&this->normalForceConstants,
                                              CartesianType_end*molecule->GetNumberAtoms());
   MallocerFreer::GetInstance()->Free<double>(&this->normalModes,
                                              CartesianType_end*molecule->GetNumberAtoms(),
                                              CartesianType_end*molecule->GetNumberAtoms());
}

void Mndo::SetMolecule(Molecule* molecule){
   ZindoS::SetMolecule(molecule);
   MallocerFreer::GetInstance()->Malloc<double>(&this->twoElecTwoCore,
                                                molecule->GetNumberAtoms(),
                                                molecule->GetNumberAtoms(),
                                                dxy, dxy, dxy, dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&this->normalForceConstants,
                                                CartesianType_end*molecule->GetNumberAtoms());
   MallocerFreer::GetInstance()->Malloc<double>(&this->normalModes,
                                                CartesianType_end*molecule->GetNumberAtoms(),
                                                CartesianType_end*molecule->GetNumberAtoms());
}
void Mndo::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in mndo::Mndo::DoSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in mndo::Mndo::DoSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in mndo::Mndo::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in mndo::Mndo::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_mndo::Mndo::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_mndo::Mndo::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageCalcCISMatrix
      = "Error in mndo::Mndo::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in mndo::Mndo::DoCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles
      = "Error in mndo:: Mndo::GetSemiEmpiricalMultipoleInteraction: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteraction1stDeriBadMultipoles
      = "Error in mndo:: Mndo::GetSemiEmpiricalMultipoleInteraction1stDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteraction2ndDeriBadMultipoles
      = "Error in mndo:: Mndo::GetSemiEmpiricalMultipoleInteraction2ndDerivative: Bad multipole combintaion is set\n";
   this->errorMessageMultipoleA = "Multipole A is: ";
   this->errorMessageMultipoleB = "Multipole B is: ";
   this->errorMessageGetNddoRepulsionIntegral 
      = "Error in mndo::Mndo::GetNddoRepulsionIntegral: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegral1stDerivative 
      = "Error in mndo::Mndo::GetNddoRepulsionIntegral1stDerivative: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegral2ndDerivative 
      = "Error in mndo::Mndo::GetNddoRepulsionIntegral2ndDerivative: Bad orbital is set.\n";
   this->errorMessageCalcTwoElecTwoCoreNullMatrix 
      = "Error in mndo::Mndo::CalcTwoElecTwoCore: The two elec two core matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCoreSameAtoms
      = "Error in mndo::Mndo::CalcDiatomicTwoElecTwoCore: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCore1stDerivativesSameAtoms
      = "Error in mndo::Mndo::CalcDiatomicTwoElecTwoCore1stDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCore2ndDerivativesSameAtoms
      = "Error in mndo::Mndo::CalcDiatomicTwoElecTwoCore2ndDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCoreNullMatrix 
      = "Error in mndo::Mndo::CalcDiatomicTwoElecTwoCore: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCore1stDerivativesNullMatrix
      = "Error in mndo::Mndo::CalcDiatomicTwoElecTwoCore1stDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCore2ndDerivativesNullMatrix
      = "Error in mndo::Mndo::CalcDiatomicTwoElecTwoCore2ndDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageGetElectronicEnergyEnergyNotCalculated
      = "Error in mndo::Mndo::GetElectronicEnergy: Set electronic state is not calculated by CIS.\n";
   this->errorMessageGetElectronicEnergyNULLCISEnergy 
      = "Error in mndo::Mndo::GetElectronicEnergy: excitedEnergies is NULL\n";
   this->errorMessageCalcZMatrixForceEtaNull 
      = "Error in mndo::Mndo::CalcZMatrixForce: Nndo::etaMatrixForce is NULL. Call Mndo::CalcEtaMatrixForce before calling Mndo::CalcZMatrixForce.\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tMNDO-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: MNDO-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: MNDO-SCF  **********\n\n\n";
   this->messageHeatsFormation = "\tHeats of formation:";
   this->messageHeatsFormationTitle = "\t\t\t\t|  [a.u.]  |  [Kcal/mol]  | \n";
   this->messageStartCIS = "**********  START: MNDO-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: MNDO-CIS  **********\n\n\n";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for MNDO-CIS met convergence criterion(^^b\n\n\n";
}

void Mndo::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}

double Mndo::GetAuxiliaryDiatomCoreRepulsionEnergy(const Atom& atomA, 
                                                   const Atom& atomB,
                                                   double distanceAB) const{
   double value=0.0;
   double alphaA = atomA.GetNddoAlpha(this->theory);
   double alphaB = atomB.GetNddoAlpha(this->theory);
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   if(atomA.GetAtomType() == H && (atomB.GetAtomType() == N || atomB.GetAtomType() == O) ){
      value = 1.0 + (distanceAB/ang2AU)*exp(-alphaB*distanceAB) + exp(-alphaA*distanceAB);
   }
   else if(atomB.GetAtomType() == H && (atomA.GetAtomType() == N || atomA.GetAtomType() == O) ){
      value = 1.0 + (distanceAB/ang2AU)*exp(-alphaA*distanceAB) + exp(-alphaB*distanceAB);
   }
   else{
      value = 1.0 + exp(-alphaA*distanceAB) + exp(-alphaB*distanceAB);
   }
   return value;
}

// First derivative of Mndo::GetAuxiliaryDiatomCoreRepulsionEnergy.
// This deivative is related to the Cartesian coordinate of atomA.
double Mndo::GetAuxiliaryDiatomCoreRepulsionEnergy1stDerivative(const Atom& atomA, 
                                                                const Atom& atomB,
                                                                double distanceAB,
                                                                CartesianType axisA) const{
   double value=0.0;
   double alphaA = atomA.GetNddoAlpha(this->theory);
   double alphaB = atomB.GetNddoAlpha(this->theory);
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double dCartesian = (atomA.GetXyz()[axisA] - atomB.GetXyz()[axisA]);
   if(atomA.GetAtomType() == H && (atomB.GetAtomType() == N || 
                                   atomB.GetAtomType() == O)  ){
      value = ((1.0/ang2AU)-alphaB*(distanceAB/ang2AU))*exp(-alphaB*distanceAB) 
             -alphaA*exp(-alphaA*distanceAB);
   }
   else if(atomB.GetAtomType() == H && (atomA.GetAtomType() == N || 
                                        atomA.GetAtomType() == O)  ){
      value = ((1.0/ang2AU)-alphaA*(distanceAB/ang2AU))*exp(-alphaA*distanceAB)
             -alphaB*exp(-alphaB*distanceAB);
   }
   else{
      value = -alphaA*exp(-alphaA*distanceAB) 
              -alphaB*exp(-alphaB*distanceAB);
   }
   value *= dCartesian/distanceAB;
   return value;
}

// Second derivative of Mndo::GetAuxiliaryDiatomCoreRepulsionEnergy.
// Both deivatives are related to the Cartesian coordinate of atomA.
double Mndo::GetAuxiliaryDiatomCoreRepulsionEnergy2ndDerivative(const Atom& atomA, 
                                                                   const Atom& atomB,
                                                                   double distanceAB,
                                                                   CartesianType axisA1,
                                                                   CartesianType axisA2) const{
   double value=0.0;
   double dCartesian1 = (atomA.GetXyz()[axisA1] - atomB.GetXyz()[axisA1]);
   double dCartesian2 = (atomA.GetXyz()[axisA2] - atomB.GetXyz()[axisA2]);
   double pre1=0.0;
   double pre2=0.0;
   if(axisA1 == axisA2){
      pre1 = 1.0/distanceAB - pow(dCartesian1,2.0)/pow(distanceAB,3.0);
      pre2 = pow(dCartesian1/distanceAB,2.0);
   }
   else{
      pre1 = -dCartesian1*dCartesian2/pow(distanceAB,3.0);
      pre2 =  dCartesian1*dCartesian2/pow(distanceAB,2.0);
   }

   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double alphaA = atomA.GetNddoAlpha(this->theory);
   double alphaB = atomB.GetNddoAlpha(this->theory);
   double fact1=0.0;
   double fact2=0.0;
   if(atomA.GetAtomType() == H && (atomB.GetAtomType() == N || 
                                   atomB.GetAtomType() == O)  ){
      fact1 = -alphaA*exp(-alphaA*distanceAB)
             +((1.0/ang2AU) - alphaB*(distanceAB/ang2AU))*exp(-alphaB*distanceAB);
      fact2 = alphaA*alphaA*exp(-alphaA*distanceAB)
             +(-2.0*alphaB/ang2AU + (distanceAB/ang2AU)*alphaB*alphaB)*exp(-alphaB*distanceAB);
   }
   else if(atomB.GetAtomType() == H && (atomA.GetAtomType() == N || 
                                        atomA.GetAtomType() == O)  ){
      fact1 = -alphaB*exp(-alphaB*distanceAB)
             +((1.0/ang2AU) - alphaA*(distanceAB/ang2AU))*exp(-alphaA*distanceAB);
      fact2 = alphaB*alphaB*exp(-alphaB*distanceAB)
             +(-2.0*alphaA/ang2AU + (distanceAB/ang2AU)*alphaA*alphaA)*exp(-alphaA*distanceAB);
   }
   else{
      fact1 = -alphaA*exp(-alphaA*distanceAB) - alphaB*exp(-alphaB*distanceAB);
      fact2 = alphaA*alphaA*exp(-alphaA*distanceAB) + alphaB*alphaB*exp(-alphaB*distanceAB);
   }
   value = pre1*fact1 + pre2*fact2;
   return value;
}

double Mndo::GetDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   double temp = this->GetAuxiliaryDiatomCoreRepulsionEnergy(atomA,
                                                             atomB,
                                                             this->molecule->GetDistanceAtoms(atomA, atomB));
   double twoElecInt = this->GetNddoRepulsionIntegral(atomA, s, s, atomB, s, s);
   return  atomA.GetCoreCharge()*atomB.GetCoreCharge()*twoElecInt*temp; 
}

// First derivative of diatomic core repulsion energy.
// This derivative is related to the coordinate of atomA.
double Mndo::GetDiatomCoreRepulsion1stDerivative(int indexAtomA,
                                                 int indexAtomB, 
                                                 CartesianType axisA) const{
   double value =0.0;
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   double distanceAB = this->molecule->GetDistanceAtoms(indexAtomA, indexAtomB);
   double twoElecInt = this->GetNddoRepulsionIntegral(atomA, s, s, atomB, s, s);
   double twoElecInt1stDeriv = this->GetNddoRepulsionIntegral1stDerivative(
                                     atomA, s, s, atomB, s, s, axisA);
   double temp = this->GetAuxiliaryDiatomCoreRepulsionEnergy(atomA, atomB, distanceAB);
   double tempDeriv = this->GetAuxiliaryDiatomCoreRepulsionEnergy1stDerivative(atomA, atomB, distanceAB, axisA);
   value = atomA.GetCoreCharge()*atomB.GetCoreCharge()
          *(twoElecInt1stDeriv*temp + twoElecInt*tempDeriv); 
   return value;
}

// Second derivative of diatomic core repulsion energy.
// Both derivatives are related to the coordinate of atomA.
double Mndo::GetDiatomCoreRepulsion2ndDerivative(int indexAtomA,
                                                    int indexAtomB, 
                                                    CartesianType axisA1,
                                                    CartesianType axisA2) const{
   double value =0.0;
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   double distanceAB = this->molecule->GetDistanceAtoms(indexAtomA, indexAtomB);
   double twoElecInt = this->GetNddoRepulsionIntegral(atomA, s, s, 
                                                      atomB, s, s);
   double twoElecInt1stDeriv1 = this->GetNddoRepulsionIntegral1stDerivative(atomA, s, s, 
                                                                            atomB, s, s, 
                                                                            axisA1);
   double twoElecInt1stDeriv2 = this->GetNddoRepulsionIntegral1stDerivative(atomA, s, s, 
                                                                            atomB, s, s, 
                                                                            axisA2);
   double twoElecInt2ndDeriv = this->GetNddoRepulsionIntegral2ndDerivative(atomA, s, s, 
                                                                           atomB, s, s, 
                                                                           axisA1, 
                                                                           axisA2);

   double temp = this->GetAuxiliaryDiatomCoreRepulsionEnergy(atomA, 
                                                             atomB, 
                                                             distanceAB);
   double temp1stDeriv1 = this->GetAuxiliaryDiatomCoreRepulsionEnergy1stDerivative(atomA, 
                                                                                   atomB, 
                                                                                   distanceAB, 
                                                                                   axisA1);
   double temp1stDeriv2 = this->GetAuxiliaryDiatomCoreRepulsionEnergy1stDerivative(atomA, 
                                                                                   atomB, 
                                                                                   distanceAB, 
                                                                                   axisA2);
   double temp2ndDeriv = this->GetAuxiliaryDiatomCoreRepulsionEnergy2ndDerivative(atomA, 
                                                                                  atomB, 
                                                                                  distanceAB, 
                                                                                  axisA1, 
                                                                                  axisA2);

   value = atomA.GetCoreCharge()*atomB.GetCoreCharge();
   value *= twoElecInt*temp2ndDeriv 
           +twoElecInt1stDeriv1*temp1stDeriv2 
           +twoElecInt1stDeriv2*temp1stDeriv1
           +twoElecInt2ndDeriv*temp;
   return value;
}

void Mndo::CalcHeatsFormation(double* heatsFormation, 
                              const Molecule& molecule) const{
   int groundState = 0;
   *heatsFormation = this->GetElectronicEnergy(groundState);
   for(int A=0; A<molecule.GetNumberAtoms(); A++){
      const Atom& atom = *molecule.GetAtom(A);
      *heatsFormation -= atom.GetMndoElecEnergyAtom();
      *heatsFormation += atom.GetMndoHeatsFormAtom();
   }
}

void Mndo::CalcSCFProperties(){
   MolDS_cndo::Cndo2::CalcSCFProperties();
   this->CalcHeatsFormation(&this->heatsFormation, *this->molecule);
 
}

void Mndo::CalcNormalModes(double** normalModes, double* normalForceConstants, const Molecule& molecule) const{
   bool isMassWeighted = true;
   this->CalcHessianSCF(normalModes, isMassWeighted);
   bool calcEigenVectors = true;
   int hessianDim = CartesianType_end*molecule.GetNumberAtoms();
   MolDS_wrappers::Lapack::GetInstance()->Dsyevd(normalModes,
                                                 normalForceConstants,
                                                 hessianDim,
                                                 calcEigenVectors);
}

void Mndo::OutputSCFResults() const{
   MolDS_cndo::Cndo2::OutputSCFResults();
   // output heats of formation
   this->OutputLog(this->messageHeatsFormationTitle);
   this->OutputLog(boost::format("%s\t%e\t%e\n\n") % this->messageHeatsFormation
                                                   % this->heatsFormation 
                                                   % (this->heatsFormation/Parameters::GetInstance()->
                                                                                       GetKcalMolin2AU()));
}

double Mndo::GetFockDiagElement(const Atom& atomA, 
                                int indexAtomA, 
                                int mu, 
                                const Molecule& molecule, 
                                double const* const* gammaAB,
                                double const* const* orbitalElectronPopulation, 
                                double const* atomicElectronPopulation,
                                double const* const* const* const* const* const* twoElecTwoCore, 
                                bool isGuess) const{
   double value=0.0;
   int firstAOIndexA = atomA.GetFirstAOIndex();
   mu -= firstAOIndexA;
   value = atomA.GetCoreIntegral(atomA.GetValence(mu), isGuess, this->theory);
   if(!isGuess){
      double temp = 0.0;
      OrbitalType orbitalMu = atomA.GetValence(mu);
      for(int nu=0; nu<atomA.GetValenceSize(); nu++){
         OrbitalType orbitalNu = atomA.GetValence(nu);
         double coulomb  = this->GetCoulombInt(orbitalMu, orbitalNu, atomA);
         double exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
         temp += orbitalElectronPopulation[nu+firstAOIndexA]
                                          [nu+firstAOIndexA]
                *(coulomb - 0.5*exchange);
      }
      value += temp;

      temp = 0.0;
      int totalNumberAtoms=molecule.GetNumberAtoms();
      for(int B=0; B<totalNumberAtoms; B++){
         if(B != indexAtomA){
            const Atom& atomB = *molecule.GetAtom(B);
            int firstAOIndexB = atomB.GetFirstAOIndex();
            int valenceSizeB = atomB.GetValenceSize();
            for(int lambda=0; lambda<valenceSizeB; lambda++){
               for(int sigma=0; sigma<valenceSizeB; sigma++){
                  temp += orbitalElectronPopulation[lambda+firstAOIndexB]
                                                   [sigma+firstAOIndexB]
                         *twoElecTwoCore[indexAtomA][B][mu][mu][lambda][sigma];
               }
               /*
               temp += MolDS_wrappers::Blas::GetInstance()->Ddot(valenceSizeB, 
                                                                 &orbitalElectronPopulation[lambda+firstAOIndexB][firstAOIndexB],
                                                                 &twoElecTwoCore[indexAtomA][B][mu][mu][lambda][0]);
               */
            }
            temp += this->GetElectronCoreAttraction(indexAtomA, 
                                                    B, 
                                                    mu, 
                                                    mu, 
                                                    twoElecTwoCore);
         }
      }
      value += temp;
   }
   return value;
}

double Mndo::GetFockOffDiagElement(const Atom& atomA, 
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
   double value = 0.0;
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   mu -= firstAOIndexA;
   nu -= firstAOIndexB;
   OrbitalType orbitalMu = atomA.GetValence(mu);
   OrbitalType orbitalNu = atomB.GetValence(nu);
   double bondParameter = 0.5*(atomA.GetBondingParameter(this->theory, orbitalMu) 
                              +atomB.GetBondingParameter(this->theory, orbitalNu)); 
   if(isGuess){
      value = bondParameter*overlapAOs[mu+firstAOIndexA][nu+firstAOIndexB];
   }
   else{
      double coulomb = 0.0;
      double exchange = 0.0;
      double temp = 0.0;
      if(indexAtomA == indexAtomB){
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalNu, atomA); 
         exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA); 
         temp = (1.5*exchange - 0.5*coulomb)
               *orbitalElectronPopulation[mu+firstAOIndexA][nu+firstAOIndexB];
         int totalNumberAtoms = molecule.GetNumberAtoms();
         for(int BB=0; BB<totalNumberAtoms; BB++){
            if(BB != indexAtomA){
               const Atom& atomBB = *molecule.GetAtom(BB);
               int firstAOIndexBB = atomBB.GetFirstAOIndex();
               int valenceSizeBB = atomBB.GetValenceSize();
               for(int lambda=0; lambda<valenceSizeBB; lambda++){
                  for(int sigma=0; sigma<valenceSizeBB; sigma++){
                     temp += orbitalElectronPopulation[lambda+firstAOIndexBB]
                                                      [sigma+firstAOIndexBB]
                            *twoElecTwoCore[indexAtomA][BB][mu][nu][lambda][sigma];
                  }
                  /*
                  temp += MolDS_wrappers::Blas::GetInstance()->Ddot(valenceSizeBB, 
                                                                    &orbitalElectronPopulation[lambda+firstAOIndexBB][firstAOIndexBB],
                                                                    &twoElecTwoCore[indexAtomA][BB][mu][nu][lambda][0]);
                  */
               }
               temp += this->GetElectronCoreAttraction(indexAtomA, 
                                                       BB, 
                                                       mu, 
                                                       nu, 
                                                       twoElecTwoCore);
            }
         }
      }
      else{
         temp = bondParameter*overlapAOs[mu+firstAOIndexA][nu+firstAOIndexB];
         for(int sigma=0; sigma<atomA.GetValenceSize(); sigma++){
            int valenceSizeB = atomB.GetValenceSize();
            for(int lambda=0; lambda<valenceSizeB; lambda++){
               //temp -= 0.5*orbitalElectronPopulation[lambda+firstAOIndexB]
               //                                     [sigma+firstAOIndexA]
               //       *twoElecTwoCore[indexAtomA][indexAtomB][mu][sigma][nu][lambda];
               temp -= 0.5*orbitalElectronPopulation[sigma+firstAOIndexA]
                                                    [lambda+firstAOIndexB]
                      *twoElecTwoCore[indexAtomA][indexAtomB][mu][sigma][nu][lambda];
            }
            /*
            temp -= 0.5*MolDS_wrappers::Blas::GetInstance()->Ddot(valenceSizeB, 
                                                              &orbitalElectronPopulation[sigma+firstAOIndexA][firstAOIndexB],
                                                              &twoElecTwoCore[indexAtomA][indexAtomB][mu][sigma][nu][0]);
            */
         }
      }
      value += temp;
   }
   return value;
}

// NDDO Coulomb Interaction
double Mndo::GetCoulombInt(OrbitalType orbital1, OrbitalType orbital2, const Atom& atom) const{
   double value=0.0;
   if( orbital1 == s && orbital2 == s){ 
      value = atom.GetNddoGss(this->theory);
   }   
   else if( orbital1 == s && ( orbital2 == px || orbital2 == py || orbital2 == pz )){ 
      value = atom.GetNddoGsp(this->theory);
   }   
   else if( orbital2 == s && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = this->GetCoulombInt(orbital2, orbital1, atom);
   }   
   else if( (orbital1 == orbital2) && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = atom.GetNddoGpp(this->theory);
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom.GetNddoGpp2(this->theory);
   }   
   else{
      stringstream ss;
      ss << this->errorMessageCoulombInt;
      ss << this->errorMessageAtomType << AtomTypeStr(atom.GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital1) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital2) << "\n";
      throw MolDSException(ss.str());
   }   
   return value;
}

// NDDO Exchange Interaction
double Mndo::GetExchangeInt(OrbitalType orbital1, OrbitalType orbital2, const Atom& atom) const{
   double value=0.0;
   if( orbital1 == orbital2){
      value = this->GetCoulombInt(orbital1, orbital2, atom);
   }   
   else if( orbital1 == s && (orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom.GetNddoHsp(this->theory);
   }   
   else if( orbital2 == s && (orbital1 == px || orbital1 == py || orbital1 == pz ) ){
      value = this->GetExchangeInt(orbital2, orbital1, atom);
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom.GetNddoHpp(this->theory);
   }
   else{
      stringstream ss;
      ss << this->errorMessageExchangeInt;
      ss << this->errorMessageAtomType << AtomTypeStr(atom.GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital1) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital2) << "\n";
      throw MolDSException(ss.str());
   }   
   return value;
}

// electron in atom A (mu and nu) and core (atom B) attraction. 
// see Eq. (16) in [DT_1977-2] with f_2 = 0.
double Mndo::GetElectronCoreAttraction(int indexAtomA, 
                                       int indexAtomB, 
                                       int mu, 
                                       int nu, 
                                       double const* const* const* const* const* const* twoElecTwoCore) const{
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   return -1.0*atomB.GetCoreCharge()*twoElecTwoCore[indexAtomA][indexAtomB][mu][nu][s][s];
}

// First derivative of electron in atom A (mu and nu) and core (atom B) attraction. 
// This derivative is related to the coordinate of atomA.
// Note that diatomicTwoElecTwoCore1stDerivative is dioatomic one.
// see Eq. (16) in [DT_1977-2] with f_2 = 0.
double Mndo::GetElectronCoreAttraction1stDerivative(int indexAtomA, 
                                                    int indexAtomB, 
                                                    int mu, 
                                                    int nu, 
                                                    double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivatives,
                                                    CartesianType axisA) const{
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   double value = -1.0*atomB.GetCoreCharge()
                  *diatomicTwoElecTwoCore1stDerivatives[mu][nu][s][s][axisA];
   return value;
}

void Mndo::CalcDiatomicOverlapAOsInDiatomicFrame(double** diatomicOverlapAOs, 
                                                 const Atom& atomA, 
                                                 const Atom& atomB) const{
   MolDS_cndo::Cndo2::CalcDiatomicOverlapAOsInDiatomicFrame(diatomicOverlapAOs, atomA, atomB);
}

// First derivative of (B.40) in J. A. Pople book without bond corrections.
void Mndo::CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(double** diatomicOverlapAOsDeri, 
                                                              const Atom& atomA, 
                                                              const Atom& atomB) const{
   MolDS_cndo::Cndo2::CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(
                      diatomicOverlapAOsDeri,atomA, atomB);
}

// Second derivative of (B.40) in J. A. Pople book without bond corrections.
void Mndo::CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(double** diatomicOverlapAOs2ndDeri, 
                                                              const Atom& atomA, 
                                                              const Atom& atomB) const{
   MolDS_cndo::Cndo2::CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(
                      diatomicOverlapAOs2ndDeri,atomA, atomB);
}

// The order of mol, moJ, moK, moL is consistent with Eq. (9) in [RZ_1973]
double Mndo::GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                         const Molecule& molecule, 
                                         double const* const* fockMatrix, 
                                         double const* const* gammaAB) const{
   double value = 0.0;
   for(int A=0; A<molecule.GetNumberAtoms(); A++){
      const Atom& atomA = *molecule.GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();

      for(int B=A; B<molecule.GetNumberAtoms(); B++){
         const Atom& atomB = *molecule.GetAtom(B);
         int firstAOIndexB = atomB.GetFirstAOIndex();
         int lastAOIndexB  = atomB.GetLastAOIndex();

         double gamma = 0.0;
         if(A!=B){
            for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
               for(int nu=mu; nu<=lastAOIndexA; nu++){
                  for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                     for(int sigma=lambda; sigma<=lastAOIndexB; sigma++){
                        OrbitalType orbitalSigma = atomB.GetValence(sigma-firstAOIndexB);
                        gamma = this->twoElecTwoCore[A]
                                                    [B]
                                                    [mu-firstAOIndexA]
                                                    [nu-firstAOIndexA]
                                                    [lambda-firstAOIndexB]
                                                    [sigma-firstAOIndexB];

                        value += gamma*fockMatrix[moI][mu]
                                      *fockMatrix[moJ][nu]
                                      *fockMatrix[moK][lambda]
                                      *fockMatrix[moL][sigma];
                        value += gamma*fockMatrix[moI][lambda]
                                      *fockMatrix[moJ][sigma]
                                      *fockMatrix[moK][mu]
                                      *fockMatrix[moL][nu];
                        if(lambda != sigma){
                           value += gamma*fockMatrix[moI][mu]
                                         *fockMatrix[moJ][nu]
                                         *fockMatrix[moK][sigma]
                                         *fockMatrix[moL][lambda];
                           value += gamma*fockMatrix[moI][sigma]
                                         *fockMatrix[moJ][lambda]
                                         *fockMatrix[moK][mu]
                                         *fockMatrix[moL][nu];
                        }
                        if(mu != nu){
                           value += gamma*fockMatrix[moI][nu]
                                         *fockMatrix[moJ][mu]
                                         *fockMatrix[moK][lambda]
                                         *fockMatrix[moL][sigma];
                           value += gamma*fockMatrix[moI][lambda]
                                         *fockMatrix[moJ][sigma]
                                         *fockMatrix[moK][nu]
                                         *fockMatrix[moL][mu];
                        }
                        if(mu != nu && lambda != sigma){
                           value += gamma*fockMatrix[moI][nu]
                                         *fockMatrix[moJ][mu]
                                         *fockMatrix[moK][sigma]
                                         *fockMatrix[moL][lambda];
                           value += gamma*fockMatrix[moI][sigma]
                                         *fockMatrix[moJ][lambda]
                                         *fockMatrix[moK][nu]
                                         *fockMatrix[moL][mu];
                        }
                     }
                  }
               }
            }
         }
         else{
            for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
               for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
                  for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                     for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
                        if(mu==nu && lambda==sigma){
                           OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                           OrbitalType orbitalLambda = atomB.GetValence(lambda-firstAOIndexB);
                           gamma = this->GetCoulombInt(orbitalMu, orbitalLambda, atomA);
                        }
                        else if((mu==lambda && nu==sigma) || (nu==lambda && mu==sigma) ){
                           OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                           OrbitalType orbitalNu = atomA.GetValence(nu-firstAOIndexA);
                           gamma = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
                        }
                        else{
                           gamma = 0.0;
                        }
                        value += gamma*fockMatrix[moI][mu]
                                      *fockMatrix[moJ][nu]
                                      *fockMatrix[moK][lambda]
                                      *fockMatrix[moL][sigma];
                     }  
                  }
               }
            }
         }
      }
   }
   return value;
}

// right-upper part is only calculated by this method.
void Mndo::CalcCISMatrix(double** matrixCIS) const{
   this->OutputLog(this->messageStartCalcCISMatrix);
   double ompStartTime = omp_get_wtime();
   int mpiRank = MolDS_mpi::MpiProcess::GetInstance()->GetRank();
   int mpiSize = MolDS_mpi::MpiProcess::GetInstance()->GetSize();

   for(int k=0; k<this->matrixCISdimension; k++){
      if(k%mpiSize != mpiRank){continue;}

      // single excitation from I-th (occupied)MO to A-th (virtual)MO
      int moI = this->GetActiveOccIndex(*this->molecule, k);
      int moA = this->GetActiveVirIndex(*this->molecule, k);
      stringstream ompErrors;
#pragma omp parallel for schedule(auto)
      for(int l=k; l<this->matrixCISdimension; l++){
         try{
            // single excitation from J-th (occupied)MO to B-th (virtual)MO
            int moJ = this->GetActiveOccIndex(*this->molecule, l);
            int moB = this->GetActiveVirIndex(*this->molecule, l);
            double value=0.0;
             
            // Fast algorith, but this is not easy to read. 
            // Slow algorithm is alos written below.
            for(int A=0; A<molecule->GetNumberAtoms(); A++){
               const Atom& atomA = *molecule->GetAtom(A);
               int firstAOIndexA = atomA.GetFirstAOIndex();
               int lastAOIndexA  = atomA.GetLastAOIndex();

               for(int B=A; B<molecule->GetNumberAtoms(); B++){
                  const Atom& atomB = *molecule->GetAtom(B);
                  int firstAOIndexB = atomB.GetFirstAOIndex();
                  int lastAOIndexB  = atomB.GetLastAOIndex();

                  double gamma = 0.0;
                  if(A!=B){
                     for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
                        for(int nu=mu; nu<=lastAOIndexA; nu++){
                           double tmpMuNu01 = 2.0*fockMatrix[moA][mu]
                                                 *fockMatrix[moI][nu];
                           double tmpMuNu02 = 2.0*fockMatrix[moJ][mu]
                                                 *fockMatrix[moB][nu];
                           double tmpMuNu03 = fockMatrix[moA][mu]
                                             *fockMatrix[moB][nu];
                           double tmpMuNu04 = fockMatrix[moI][mu]
                                             *fockMatrix[moJ][nu];
                           double tmpMuNu09 = 2.0*fockMatrix[moI][mu]
                                                 *fockMatrix[moA][nu];
                           double tmpMuNu10 = 2.0*fockMatrix[moB][mu]
                                                 *fockMatrix[moJ][nu];
                           double tmpMuNu11 = fockMatrix[moB][mu]
                                             *fockMatrix[moA][nu];
                           double tmpMuNu12 = fockMatrix[moJ][mu]
                                             *fockMatrix[moI][nu];
                           for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                              double tmpMuNuLamda01 = tmpMuNu01*fockMatrix[moJ][lambda];
                              double tmpMuNuLamda02 = tmpMuNu02*fockMatrix[moA][lambda];
                              double tmpMuNuLamda03 = tmpMuNu03*fockMatrix[moI][lambda];
                              double tmpMuNuLamda04 = tmpMuNu04*fockMatrix[moA][lambda];
                              double tmpMuNuLamda05 = tmpMuNu01*fockMatrix[moB][lambda];
                              double tmpMuNuLamda06 = tmpMuNu02*fockMatrix[moI][lambda];
                              double tmpMuNuLamda07 = tmpMuNu03*fockMatrix[moJ][lambda];
                              double tmpMuNuLamda08 = tmpMuNu04*fockMatrix[moB][lambda];
                              double tmpMuNuLamda09 = tmpMuNu09*fockMatrix[moJ][lambda];
                              double tmpMuNuLamda10 = tmpMuNu10*fockMatrix[moA][lambda];
                              double tmpMuNuLamda11 = tmpMuNu11*fockMatrix[moI][lambda];
                              double tmpMuNuLamda12 = tmpMuNu12*fockMatrix[moA][lambda];
                              double tmpMuNuLamda13 = tmpMuNu09*fockMatrix[moB][lambda];
                              double tmpMuNuLamda14 = tmpMuNu10*fockMatrix[moI][lambda];
                              double tmpMuNuLamda15 = tmpMuNu11*fockMatrix[moJ][lambda];
                              double tmpMuNuLamda16 = tmpMuNu12*fockMatrix[moB][lambda];
                              for(int sigma=lambda; sigma<=lastAOIndexB; sigma++){
                                 OrbitalType orbitalSigma = atomB.GetValence(sigma-firstAOIndexB);
                                 gamma = this->twoElecTwoCore[A]
                                                             [B]
                                                             [mu-firstAOIndexA]
                                                             [nu-firstAOIndexA]
                                                             [lambda-firstAOIndexB]
                                                             [sigma-firstAOIndexB];
   
                                 value += gamma*tmpMuNuLamda01*fockMatrix[moB][sigma];
                                 value += gamma*tmpMuNuLamda02*fockMatrix[moI][sigma];
                                 value -= gamma*tmpMuNuLamda03*fockMatrix[moJ][sigma];
                                 value -= gamma*tmpMuNuLamda04*fockMatrix[moB][sigma];
                                 if(lambda != sigma){
                                    value += gamma*tmpMuNuLamda05*fockMatrix[moJ][sigma];
                                    value += gamma*tmpMuNuLamda06*fockMatrix[moA][sigma];
                                    value -= gamma*tmpMuNuLamda07*fockMatrix[moI][sigma];
                                    value -= gamma*tmpMuNuLamda08*fockMatrix[moA][sigma];
                                 }
                                 if(mu != nu){
                                    value += gamma*tmpMuNuLamda09*fockMatrix[moB][sigma];
                                    value += gamma*tmpMuNuLamda10*fockMatrix[moI][sigma];
                                    value -= gamma*tmpMuNuLamda11*fockMatrix[moJ][sigma];
                                    value -= gamma*tmpMuNuLamda12*fockMatrix[moB][sigma];
                                 }
                                 if(mu != nu && lambda != sigma){
                                    value += gamma*tmpMuNuLamda13*fockMatrix[moJ][sigma];
                                    value += gamma*tmpMuNuLamda14*fockMatrix[moA][sigma];
                                    value -= gamma*tmpMuNuLamda15*fockMatrix[moI][sigma];
                                    value -= gamma*tmpMuNuLamda16*fockMatrix[moA][sigma];
                                 }
                              }
                           }
                        }
                     }
                  }
                  else{
                     for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
                        for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
                           double tmpMuNu01 = 2.0*fockMatrix[moA][mu]
                                                 *fockMatrix[moI][nu];
                           double tmpMuNu02 = fockMatrix[moA][mu]
                                             *fockMatrix[moB][nu];
                           for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                              double tmpMuNuLamda01 = tmpMuNu01*fockMatrix[moJ][lambda];
                              double tmpMuNuLamda02 = tmpMuNu02*fockMatrix[moI][lambda];
                              for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
                                 if(mu==nu && lambda==sigma){
                                    OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                                    OrbitalType orbitalLambda = atomB.GetValence(lambda-firstAOIndexB);
                                    gamma = this->GetCoulombInt(orbitalMu, orbitalLambda, atomA);
                                 }
                                 else if((mu==lambda && nu==sigma) || (nu==lambda && mu==sigma) ){
                                    OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                                    OrbitalType orbitalNu = atomA.GetValence(nu-firstAOIndexA);
                                    gamma = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
                                 }
                                 else{
                                    gamma = 0.0;
                                 }
                                 value += gamma*tmpMuNuLamda01*fockMatrix[moB][sigma];
                                 value -= gamma*tmpMuNuLamda02*fockMatrix[moJ][sigma];
                              }  
                           }
                        }
                     }
                  }
               }
            }
            // End of the fast algorith.
            
            /* 
            // Slow algorith, but this is easy to read. Fast altorithm is also written above.
            value = 2.0*this->GetMolecularIntegralElement(moA, moI, moJ, moB, 
                                                          *this->molecule, this->fockMatrix, NULL)
                       -this->GetMolecularIntegralElement(moA, moB, moI, moJ, 
                                                          *this->molecule, this->fockMatrix, NULL);
            // End of the slow algorith.
            */
            // Diagonal term
            if(k==l){
               value += this->energiesMO[moA] - this->energiesMO[moI];
            }
            matrixCIS[k][l] = value;
         } 
         catch(MolDSException ex){
#pragma omp critical
            ex.Serialize(ompErrors);
         }
      }// end of l-loop
      // Exception throwing for omp-region
      if(!ompErrors.str().empty()){
         throw MolDSException::Deserialize(ompErrors);
      }
   } // end of k-loop

   // communication to collect all matrix data on rank 0
   if(mpiRank == 0){
      // receive the matrix data from other ranks
      for(int k=0; k<this->matrixCISdimension; k++){
         if(k%mpiSize == 0){continue;}
         int source = k%mpiSize;
         int tag = k;
         MolDS_mpi::MpiProcess::GetInstance()->Recv(source, tag, matrixCIS[k], this->matrixCISdimension);
      }
   }
   else{
      // send the matrix data to rank-0
      for(int k=0; k<this->matrixCISdimension; k++){
         if(k%mpiSize != mpiRank){continue;}
         int dest = 0;
         int tag = k;
         MolDS_mpi::MpiProcess::GetInstance()->Send(dest, tag, matrixCIS[k], this->matrixCISdimension);
      }
   }
   // broadcast all matrix data to all rank
   int root=0;
   MolDS_mpi::MpiProcess::GetInstance()->Broadcast(&matrixCIS[0][0], this->matrixCISdimension*this->matrixCISdimension, root);


   double ompEndTime = omp_get_wtime();
   this->OutputLog(boost::format("%s%lf%s\n%s") % this->messageOmpElapsedTimeCalcCISMarix.c_str()
                                                % (ompEndTime - ompStartTime)
                                                % this->messageUnitSec.c_str()
                                                % this->messageDoneCalcCISMatrix.c_str());
}

// \epsilon_{r}^{kl} in (1) in [PT_1997].
// k and l are index of CIS matrix.
double Mndo::GetCISCoefficientMOEnergy(int k, int l, int r, int numberActiveVir) const{
   double value=0.0;
   if(k==l){
      int moI = this->GetActiveOccIndex(*this->molecule, k);
      int moA = this->GetActiveVirIndex(*this->molecule, k);
      if(r==moI){
         // r is index of occupied MO.
         value = -1.0;
      }
      else if(r==moA){
         // r is index of virtual MO.
         value = 1.0;
      }
   }
   return value;
}

// \f_{pqrs}^{lm} in (1) in [PT_1997].
// k and l are index of CIS matrix.
double Mndo::GetCISCoefficientTwoElecIntegral(int k, 
                                              int l, 
                                              int p, 
                                              int q, 
                                              int r, 
                                              int s,
                                              int numberActiveVir) const{
   double value=0.0;
   // single excitation from I-th (occupied)MO to A-th (virtual)MO
   int moI = this->GetActiveOccIndex(*this->molecule, k);
   int moA = this->GetActiveVirIndex(*this->molecule, k);
   // single excitation from J-th (occupied)MO to B-th (virtual)MO
   int moJ = this->GetActiveOccIndex(*this->molecule, l);
   int moB = this->GetActiveVirIndex(*this->molecule, l);
   if(p==moI && q==moA && r==moJ && s==moB ){
      value = 2.0;
   }
   else if(p==moI && q==moJ && r==moA && s==moB ){
      value = -1.0;
   }
   return value;
}

void Mndo::MallocTempMatricesEachThreadCalcHessianSCF(double***** diatomicOverlapAOs1stDerivs,
                                                      double****** diatomicOverlapAOs2ndDerivs,
                                                      double******* diatomicTwoElecTwoCore1stDerivs,
                                                      double******** diatomicTwoElecTwoCore2ndDerivs) const{
   MallocerFreer::GetInstance()->Malloc<double>(diatomicOverlapAOs1stDerivs,
                                                this->molecule->GetNumberAtoms(),
                                                OrbitalType_end,
                                                OrbitalType_end,
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicOverlapAOs2ndDerivs,
                                                this->molecule->GetNumberAtoms(),
                                                OrbitalType_end,
                                                OrbitalType_end,
                                                CartesianType_end,
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicTwoElecTwoCore1stDerivs,
                                                this->molecule->GetNumberAtoms(),
                                                dxy,
                                                dxy,
                                                dxy,
                                                dxy,
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicTwoElecTwoCore2ndDerivs,
                                                this->molecule->GetNumberAtoms(),
                                                dxy,
                                                dxy,
                                                dxy,
                                                dxy,
                                                CartesianType_end,
                                                CartesianType_end);
}

void Mndo::FreeTempMatricesEachThreadCalcHessianSCF(double*****    diatomicOverlapAOs1stDerivs,
                                                    double******   diatomicOverlapAOs2ndDerivs,
                                                    double*******  diatomicTwoElecTwoCore1stDerivs,
                                                    double******** diatomicTwoElecTwoCore2ndDerivs) const{
   MallocerFreer::GetInstance()->Free<double>(diatomicOverlapAOs1stDerivs,
                                              this->molecule->GetNumberAtoms(),
                                              OrbitalType_end,
                                              OrbitalType_end,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicOverlapAOs2ndDerivs,
                                              this->molecule->GetNumberAtoms(),
                                              OrbitalType_end,
                                              OrbitalType_end,
                                              CartesianType_end,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicTwoElecTwoCore1stDerivs,
                                              this->molecule->GetNumberAtoms(),
                                              dxy,
                                              dxy,
                                              dxy,
                                              dxy,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicTwoElecTwoCore2ndDerivs,
                                              this->molecule->GetNumberAtoms(),
                                              dxy,
                                              dxy,
                                              dxy,
                                              dxy,
                                              CartesianType_end,
                                              CartesianType_end);
}

// mu and nu is included in atomA' AO. 
// s is included in atomC's AO.
// Both derivatives are related to the Cartesian coordinates (axisA1 and axisA2) of atomA.
double Mndo::GetAuxiliaryHessianElement1(int mu, 
                                         int nu, 
                                         int indexAtomA,
                                         int indexAtomC,
                                         CartesianType axisA1,
                                         CartesianType axisA2,
                                         double const* const* orbitalElectronPopulation,
                                         double const* const* const* const* const* const* diatomicTwoElecTwoCore2ndDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomC = *this->molecule->GetAtom(indexAtomC);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   double value = orbitalElectronPopulation[mu]
                                           [nu]
                 *diatomicTwoElecTwoCore2ndDerivs[mu-firstAOIndexA]
                                                    [nu-firstAOIndexA]
                                                    [s]
                                                    [s]
                                                    [axisA1]
                                                    [axisA2];
   return value*atomC.GetCoreCharge();
}

// mu and nu is included in atomA' AO. 
// s is included in atomC's AO.
// Derivtive of orbitalElectronPopulation is reralted to the Cartesian coordinate (axisB) of atomB.
double Mndo::GetAuxiliaryHessianElement2(int mu, 
                                         int nu, 
                                         int indexAtomA,
                                         int indexAtomB,
                                         int indexAtomC,
                                         CartesianType axisA,
                                         CartesianType axisB,
                                         double const* const* const* const* orbitalElectronPopulation1stDerivs,
                                         double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomC = *this->molecule->GetAtom(indexAtomC);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   double value = orbitalElectronPopulation1stDerivs[mu]
                                                    [nu]
                                                    [indexAtomB]
                                                    [axisB]
                 *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA]
                                                 [nu-firstAOIndexA]
                                                 [s]
                                                 [s]
                                                 [axisA];
   return value*atomC.GetCoreCharge();
}

// lambda and sigma is included in atomC' AO. 
// s is included in atomA's AO.
// Both derivatives are related to the Cartesian coordinates (axisA1 and axisA2) of atomA.
double Mndo::GetAuxiliaryHessianElement3(int lambda, 
                                         int sigma, 
                                         int indexAtomA,
                                         int indexAtomC,
                                         CartesianType axisA1,
                                         CartesianType axisA2,
                                         double const* const* orbitalElectronPopulation,
                                         double const* const* const* const* const* const* diatomicTwoElecTwoCore2ndDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomC = *this->molecule->GetAtom(indexAtomC);
   int firstAOIndexC = atomC.GetFirstAOIndex();
   double value = orbitalElectronPopulation[lambda]
                                           [sigma]
                 *diatomicTwoElecTwoCore2ndDerivs[s]
                                                    [s]
                                                    [lambda-firstAOIndexC]
                                                    [sigma-firstAOIndexC]
                                                    [axisA1]
                                                    [axisA2];
   return value*atomA.GetCoreCharge();
}

// lambda and sigma is included in atomC' AO. 
// s is included in atomA's AO.
// Derivtive of orbitalElectronPopulation is reralted to the Cartesian coordinate (axisB) of atomB.
double Mndo::GetAuxiliaryHessianElement4(int lambda, 
                                         int sigma, 
                                         int indexAtomA,
                                         int indexAtomB,
                                         int indexAtomC,
                                         CartesianType axisA,
                                         CartesianType axisB,
                                         double const* const* const* const* orbitalElectronPopulation1stDerivs,
                                         double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomC = *this->molecule->GetAtom(indexAtomC);
   int firstAOIndexC = atomC.GetFirstAOIndex();
   double value = orbitalElectronPopulation1stDerivs[lambda]
                                                    [sigma]
                                                    [indexAtomB]
                                                    [axisB]
                 *diatomicTwoElecTwoCore1stDerivs[s]
                                                 [s]
                                                 [lambda-firstAOIndexC]
                                                 [sigma-firstAOIndexC]
                                                 [axisA];
   return value*atomA.GetCoreCharge();
}

// mu is included in atomA's AO.
// lambda is included in atomC's AO.
// Both derivatives are related to the Cartesian coordinates (axisA1 and axisA2) of atomA.
double Mndo::GetAuxiliaryHessianElement5(int mu, 
                                         int lambda, 
                                         int indexAtomA,
                                         int indexAtomC,
                                         CartesianType axisA1,
                                         CartesianType axisA2,
                                         double const* const* orbitalElectronPopulation,
                                         double const* const* const* const* diatomicOverlapAOs2ndDerivs) const{
   const Atom& atomA = *molecule->GetAtom(indexAtomA);
   const Atom& atomC = *molecule->GetAtom(indexAtomC);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexC = atomC.GetFirstAOIndex();
   double bondParameterA = atomA.GetBondingParameter(this->theory, atomA.GetValence(mu-firstAOIndexA));
   double bondParameterC = atomC.GetBondingParameter(this->theory, atomC.GetValence(lambda-firstAOIndexC));
   double sumBondParameters = bondParameterA+bondParameterC;
   double value = orbitalElectronPopulation[mu][lambda]
                 *sumBondParameters
                 *diatomicOverlapAOs2ndDerivs[mu-firstAOIndexA]
                                             [lambda-firstAOIndexC]
                                             [axisA1]
                                             [axisA2];
   return value;
}

// mu is included in atomA's AO.
// lambda is included in atomC's AO.
// Derivtive of orbitalElectronPopulation is reralted to the Cartesian coordinate (axisB) of atomB.
double Mndo::GetAuxiliaryHessianElement6(int mu, 
                                         int lambda, 
                                         int indexAtomA,
                                         int indexAtomB,
                                         int indexAtomC,
                                         CartesianType axisA,
                                         CartesianType axisB,
                                         double const* const* const* const* orbitalElectronPopulation1stDerivs,
                                         double const* const* const* diatomicOverlapAOs1stDerivs) const{
   const Atom& atomA = *molecule->GetAtom(indexAtomA);
   const Atom& atomC = *molecule->GetAtom(indexAtomC);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexC = atomC.GetFirstAOIndex();
   double bondParameterA = atomA.GetBondingParameter(this->theory, atomA.GetValence(mu-firstAOIndexA));
   double bondParameterC = atomC.GetBondingParameter(this->theory, atomC.GetValence(lambda-firstAOIndexC));
   double sumBondParameters = bondParameterA+bondParameterC;
   double value = orbitalElectronPopulation1stDerivs[mu]
                                                    [lambda]
                                                    [indexAtomB]
                                                    [axisB]
                 *sumBondParameters
                 *diatomicOverlapAOs1stDerivs[mu-firstAOIndexA]
                                             [lambda-firstAOIndexC]
                                             [axisA];
   return value;
}

// mu and nu are included in atomA's AO.
// lambda and sigma are included in atomC's AO.
// Both derivatives are related to the Cartesian coordinates (axisA1 and axisA2) of atomA.
double Mndo::GetAuxiliaryHessianElement7(int mu, 
                                         int nu, 
                                         int lambda, 
                                         int sigma, 
                                         int indexAtomA,
                                         int indexAtomC,
                                         CartesianType axisA1,
                                         CartesianType axisA2,
                                         double const* const* orbitalElectronPopulation,
                                         double const* const* const* const* const* const* diatomicTwoElecTwoCore2ndDerivs) const{
   const Atom& atomA = *molecule->GetAtom(indexAtomA);
   const Atom& atomC = *molecule->GetAtom(indexAtomC);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexC = atomC.GetFirstAOIndex();
   double temp1 = orbitalElectronPopulation[mu][nu]*orbitalElectronPopulation[lambda][sigma];
   double temp2 = orbitalElectronPopulation[mu][lambda]*orbitalElectronPopulation[nu][sigma];
   double value = (temp1 - 0.5*temp2)
                 *diatomicTwoElecTwoCore2ndDerivs[mu-firstAOIndexA]
                                                    [nu-firstAOIndexA]
                                                    [lambda-firstAOIndexC]
                                                    [sigma-firstAOIndexC]
                                                    [axisA1]
                                                    [axisA2];
   return value;
}

// mu and nu are included in atomA's AO.
// lambda and sigma are included in atomC's AO.
// Derivtive of orbitalElectronPopulation is reralted to the Cartesian coordinate (axisB) of atomB.
double Mndo::GetAuxiliaryHessianElement8(int mu, 
                                         int nu, 
                                         int lambda, 
                                         int sigma, 
                                         int indexAtomA,
                                         int indexAtomB,
                                         int indexAtomC,
                                         CartesianType axisA,
                                         CartesianType axisB,
                                         double const* const* orbitalElectronPopulation,
                                         double const* const* const* const* orbitalElectronPopulation1stDerivs,
                                         double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const{
   const Atom& atomA = *molecule->GetAtom(indexAtomA);
   const Atom& atomC = *molecule->GetAtom(indexAtomC);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexC = atomC.GetFirstAOIndex();
   double temp1 = orbitalElectronPopulation1stDerivs[mu][nu]       [indexAtomB][axisB]
                 *orbitalElectronPopulation         [lambda][sigma];
   double temp2 = orbitalElectronPopulation         [mu][nu]
                 *orbitalElectronPopulation1stDerivs[lambda][sigma][indexAtomB][axisB];
   double temp3 = orbitalElectronPopulation1stDerivs[mu][lambda]   [indexAtomB][axisB]
                 *orbitalElectronPopulation         [nu][sigma];
   double temp4 = orbitalElectronPopulation         [mu][lambda]
                 *orbitalElectronPopulation1stDerivs[nu][sigma]    [indexAtomB][axisB];
   double value = ((temp1 + temp2) - 0.5*(temp3 + temp4))
                 *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA]
                                                 [nu-firstAOIndexA]
                                                 [lambda-firstAOIndexC]
                                                 [sigma-firstAOIndexC]
                                                 [axisA];
   return value;
}

// Return hessian element. 
// The Second derivative are related to axisA1 and axisA2.
// These axisA1 and axisA2 are the Cartesian coordinates of atomA labeled with indexAtomA.
double Mndo::GetHessianElementSameAtomsSCF(int indexAtomA, 
                                           CartesianType axisA1,
                                           CartesianType axisA2,
                                           double const* const*               orbitalElectronPopulation,
                                           double const* const* const* const* orbitalElectronPopulation1stDerivs,
                                           double const* const* const* const*        diatomicOverlapAOs1stDerivs,
                                           double const* const* const* const* const* diatomicOverlapAOs2ndDerivs,
                                           double const* const* const* const* const* const*        diatomicTwoElecTwoCore1stDerivs,
                                           double const* const* const* const* const* const* const* diatomicTwoElecTwoCore2ndDerivs) const{
   double value=0.0;
   int indexAtomB = indexAtomA;
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int lastAOIndexA  = atomA.GetLastAOIndex();
   for(int indexAtomC=0; indexAtomC<this->molecule->GetNumberAtoms(); indexAtomC++){
      if(indexAtomA != indexAtomC){
         const Atom& atomC = *this->molecule->GetAtom(indexAtomC);
         int firstAOIndexC = atomC.GetFirstAOIndex();
         int numberAOsC = atomC.GetValenceSize();

         // second derivative of electronic part
         for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
            for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
               value -= this->GetAuxiliaryHessianElement1(mu, 
                                                          nu, 
                                                          indexAtomA,
                                                          indexAtomC,
                                                          static_cast<CartesianType>(axisA1), 
                                                          static_cast<CartesianType>(axisA2), 
                                                          orbitalElectronPopulation,
                                                          diatomicTwoElecTwoCore2ndDerivs[indexAtomC]);
               value -= this->GetAuxiliaryHessianElement2(mu, 
                                                          nu, 
                                                          indexAtomA,
                                                          indexAtomB,
                                                          indexAtomC,
                                                          static_cast<CartesianType>(axisA1), 
                                                          static_cast<CartesianType>(axisA2), 
                                                          orbitalElectronPopulation1stDerivs,
                                                          diatomicTwoElecTwoCore1stDerivs[indexAtomC]);
            }
         }
         for(int lambda=firstAOIndexC; lambda<firstAOIndexC+numberAOsC; lambda++){
            for(int sigma=firstAOIndexC; sigma<firstAOIndexC+numberAOsC; sigma++){
               value -= this->GetAuxiliaryHessianElement3(lambda, 
                                                          sigma, 
                                                          indexAtomA,
                                                          indexAtomC,
                                                          static_cast<CartesianType>(axisA1), 
                                                          static_cast<CartesianType>(axisA2), 
                                                          orbitalElectronPopulation,
                                                          diatomicTwoElecTwoCore2ndDerivs[indexAtomC]);
               value -= this->GetAuxiliaryHessianElement4(lambda, 
                                                          sigma, 
                                                          indexAtomA,
                                                          indexAtomB,
                                                          indexAtomC,
                                                          static_cast<CartesianType>(axisA1), 
                                                          static_cast<CartesianType>(axisA2), 
                                                          orbitalElectronPopulation1stDerivs,
                                                          diatomicTwoElecTwoCore1stDerivs[indexAtomC]);
            }
         }
         for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
            for(int lambda=firstAOIndexC; lambda<firstAOIndexC+numberAOsC; lambda++){
               value += this->GetAuxiliaryHessianElement5(mu, 
                                                          lambda, 
                                                          indexAtomA,
                                                          indexAtomC,
                                                          static_cast<CartesianType>(axisA1), 
                                                          static_cast<CartesianType>(axisA2), 
                                                          orbitalElectronPopulation,
                                                          diatomicOverlapAOs2ndDerivs[indexAtomC]);
               value += this->GetAuxiliaryHessianElement6(mu, 
                                                          lambda, 
                                                          indexAtomA,
                                                          indexAtomB,
                                                          indexAtomC,
                                                          static_cast<CartesianType>(axisA1), 
                                                          static_cast<CartesianType>(axisA2), 
                                                          orbitalElectronPopulation1stDerivs,
                                                          diatomicOverlapAOs1stDerivs[indexAtomC]);
            }
         }
         for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
            for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
               for(int lambda=firstAOIndexC; lambda<firstAOIndexC+numberAOsC; lambda++){
                  for(int sigma=firstAOIndexC; sigma<firstAOIndexC+numberAOsC; sigma++){
                     value += this->GetAuxiliaryHessianElement7(mu, 
                                                                nu, 
                                                                lambda, 
                                                                sigma, 
                                                                indexAtomA,
                                                                indexAtomC,
                                                                static_cast<CartesianType>(axisA1), 
                                                                static_cast<CartesianType>(axisA2), 
                                                                orbitalElectronPopulation,
                                                                diatomicTwoElecTwoCore2ndDerivs[indexAtomC]);
                     value += this->GetAuxiliaryHessianElement8(mu, 
                                                                nu, 
                                                                lambda, 
                                                                sigma, 
                                                                indexAtomA,
                                                                indexAtomB,
                                                                indexAtomC,
                                                                static_cast<CartesianType>(axisA1), 
                                                                static_cast<CartesianType>(axisA2), 
                                                                orbitalElectronPopulation,
                                                                orbitalElectronPopulation1stDerivs,
                                                                diatomicTwoElecTwoCore1stDerivs[indexAtomC]);
                  }
               }
            }
         }

         // second derivatives of the nuclear repulsions
         value += this->GetDiatomCoreRepulsion2ndDerivative(indexAtomA, 
                                                            indexAtomC, 
                                                            static_cast<CartesianType>(axisA1), 
                                                            static_cast<CartesianType>(axisA2));
         // second derivatives of the van der waals corrections
         if(Parameters::GetInstance()->RequiresVdWSCF()){
            value += this->GetDiatomVdWCorrection2ndDerivative(indexAtomA, 
                                                               indexAtomC, 
                                                               static_cast<CartesianType>(axisA1), 
                                                               static_cast<CartesianType>(axisA2));
         }
      }
   }

   return value;
}

// Return hessian element. 
// The Second derivative are related to axisA and axisB.
// These axisA and axisB are the Cartesian coordinates of atomA and atomB, respectively.
double Mndo::GetHessianElementDifferentAtomsSCF(int indexAtomA, 
                                                int indexAtomB,
                                                CartesianType axisA,
                                                CartesianType axisB,
                                                double const* const*               orbitalElectronPopulation,
                                                double const* const* const* const* orbitalElectronPopulation1stDerivs,
                                                double const* const* const* const*        diatomicOverlapAOs1stDerivs,
                                                double const* const* const* const* const* diatomicOverlapAOs2ndDerivs,
                                                double const* const* const* const* const* const*        diatomicTwoElecTwoCore1stDerivs,
                                                double const* const* const* const* const* const* const* diatomicTwoElecTwoCore2ndDerivs) const{
   double value=0.0;
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int lastAOIndexA  = atomA.GetLastAOIndex();
   int lastAOIndexB  = atomB.GetLastAOIndex();

   // second derivative of electronic part
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
         value += this->GetAuxiliaryHessianElement1(mu, 
                                                    nu, 
                                                    indexAtomA,
                                                    indexAtomB,
                                                    static_cast<CartesianType>(axisA), 
                                                    static_cast<CartesianType>(axisB), 
                                                    orbitalElectronPopulation,
                                                    diatomicTwoElecTwoCore2ndDerivs[indexAtomB]);
      }
   }
   for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
      for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
         value += this->GetAuxiliaryHessianElement3(lambda, 
                                                    sigma, 
                                                    indexAtomA,
                                                    indexAtomB,
                                                    static_cast<CartesianType>(axisA), 
                                                    static_cast<CartesianType>(axisB), 
                                                    orbitalElectronPopulation,
                                                    diatomicTwoElecTwoCore2ndDerivs[indexAtomB]);
      }
   }
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
         value -= this->GetAuxiliaryHessianElement5(mu, 
                                                    lambda, 
                                                    indexAtomA,
                                                    indexAtomB,
                                                    static_cast<CartesianType>(axisA), 
                                                    static_cast<CartesianType>(axisB), 
                                                    orbitalElectronPopulation,
                                                    diatomicOverlapAOs2ndDerivs[indexAtomB]);
      }
   }
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
         for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
            for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
               value -= this->GetAuxiliaryHessianElement7(mu, 
                                                          nu, 
                                                          lambda, 
                                                          sigma, 
                                                          indexAtomA,
                                                          indexAtomB,
                                                          static_cast<CartesianType>(axisA), 
                                                          static_cast<CartesianType>(axisB), 
                                                          orbitalElectronPopulation,
                                                          diatomicTwoElecTwoCore2ndDerivs[indexAtomB]);
            }
         }
      }
   }

   for(int indexAtomC=0; indexAtomC<this->molecule->GetNumberAtoms(); indexAtomC++){
      if(indexAtomA != indexAtomC){
         const Atom& atomC = *this->molecule->GetAtom(indexAtomC);
         int firstAOIndexC = atomC.GetFirstAOIndex();
         int numberAOsC = atomC.GetValenceSize();

         // second derivative of electronic part
         for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
            for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
               value -= this->GetAuxiliaryHessianElement2(mu, 
                                                          nu, 
                                                          indexAtomA,
                                                          indexAtomB,
                                                          indexAtomC,
                                                          static_cast<CartesianType>(axisA), 
                                                          static_cast<CartesianType>(axisB), 
                                                          orbitalElectronPopulation1stDerivs,
                                                          diatomicTwoElecTwoCore1stDerivs[indexAtomC]);
            }
         }
         for(int lambda=firstAOIndexC; lambda<firstAOIndexC+numberAOsC; lambda++){
            for(int sigma=firstAOIndexC; sigma<firstAOIndexC+numberAOsC; sigma++){
               value -= this->GetAuxiliaryHessianElement4(lambda, 
                                                          sigma, 
                                                          indexAtomA,
                                                          indexAtomB,
                                                          indexAtomC,
                                                          static_cast<CartesianType>(axisA), 
                                                          static_cast<CartesianType>(axisB), 
                                                          orbitalElectronPopulation1stDerivs,
                                                          diatomicTwoElecTwoCore1stDerivs[indexAtomC]);
            }
         }
         for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
            for(int lambda=firstAOIndexC; lambda<firstAOIndexC+numberAOsC; lambda++){
               value += this->GetAuxiliaryHessianElement6(mu, 
                                                          lambda, 
                                                          indexAtomA,
                                                          indexAtomB,
                                                          indexAtomC,
                                                          static_cast<CartesianType>(axisA), 
                                                          static_cast<CartesianType>(axisB), 
                                                          orbitalElectronPopulation1stDerivs,
                                                          diatomicOverlapAOs1stDerivs[indexAtomC]);
            }
         }
         for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
            for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
               for(int lambda=firstAOIndexC; lambda<firstAOIndexC+numberAOsC; lambda++){
                  for(int sigma=firstAOIndexC; sigma<firstAOIndexC+numberAOsC; sigma++){
                     value += this->GetAuxiliaryHessianElement8(mu, 
                                                                nu, 
                                                                lambda, 
                                                                sigma, 
                                                                indexAtomA,
                                                                indexAtomB,
                                                                indexAtomC,
                                                                static_cast<CartesianType>(axisA), 
                                                                static_cast<CartesianType>(axisB), 
                                                                orbitalElectronPopulation,
                                                                orbitalElectronPopulation1stDerivs,
                                                                diatomicTwoElecTwoCore1stDerivs[indexAtomC]);
                  }
               }
            }
         }
      }
   }

   // second derivatives of the nuclear repulsions
   value -= this->GetDiatomCoreRepulsion2ndDerivative(indexAtomA, 
                                                      indexAtomB, 
                                                      static_cast<CartesianType>(axisA), 
                                                      static_cast<CartesianType>(axisB));
   // second derivatives of the van der waals corrections
   if(Parameters::GetInstance()->RequiresVdWSCF()){
      value -= this->GetDiatomVdWCorrection2ndDerivative(indexAtomA, 
                                                         indexAtomB, 
                                                         static_cast<CartesianType>(axisA), 
                                                         static_cast<CartesianType>(axisB));
   }

   return value;
}

void Mndo::CalcHessianSCF(double** hessianSCF, bool isMassWeighted) const{
   int totalNumberAOs = this->molecule->GetTotalNumberAOs();
   double**** orbitalElectronPopulation1stDerivs = NULL;

   try{
      MallocerFreer::GetInstance()->Malloc<double>(&orbitalElectronPopulation1stDerivs, 
                                                   totalNumberAOs,
                                                   totalNumberAOs,
                                                   this->molecule->GetNumberAtoms(),
                                                   CartesianType_end);
      this->CalcOrbitalElectronPopulation1stDerivatives(orbitalElectronPopulation1stDerivs);

      stringstream ompErrors;
#pragma omp parallel
      {
         double****    diatomicOverlapAOs1stDerivs = NULL;
         double*****   diatomicOverlapAOs2ndDerivs = NULL;
         double******  diatomicTwoElecTwoCore1stDerivs = NULL;
         double******* diatomicTwoElecTwoCore2ndDerivs = NULL;
         this->MallocTempMatricesEachThreadCalcHessianSCF(&diatomicOverlapAOs1stDerivs,
                                                          &diatomicOverlapAOs2ndDerivs,
                                                          &diatomicTwoElecTwoCore1stDerivs, 
                                                          &diatomicTwoElecTwoCore2ndDerivs);
         try{
#pragma omp for schedule(auto)                                                 
            for(int indexAtomA=0; indexAtomA<this->molecule->GetNumberAtoms(); indexAtomA++){
               const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
               int firstAOIndexA = atomA.GetFirstAOIndex();
               int lastAOIndexA  = atomA.GetLastAOIndex();
               for(int axisA = XAxis; axisA<CartesianType_end; axisA++){
            
                  // calculation of derivatives of the overlapAOss and two electron integrals
                  for(int indexAtomB=0; indexAtomB<this->molecule->GetNumberAtoms(); indexAtomB++){
                     if(indexAtomA != indexAtomB){
                        this->CalcDiatomicOverlapAOs1stDerivatives(diatomicOverlapAOs1stDerivs[indexAtomB], 
                                                                   indexAtomA, 
                                                                   indexAtomB);
                        this->CalcDiatomicOverlapAOs2ndDerivatives(diatomicOverlapAOs2ndDerivs[indexAtomB], 
                                                                   indexAtomA, 
                                                                   indexAtomB);
                        this->CalcDiatomicTwoElecTwoCore1stDerivatives(diatomicTwoElecTwoCore1stDerivs[indexAtomB], 
                                                                       indexAtomA, 
                                                                       indexAtomB);
                        this->CalcDiatomicTwoElecTwoCore2ndDerivatives(diatomicTwoElecTwoCore2ndDerivs[indexAtomB], 
                                                                       indexAtomA, 
                                                                       indexAtomB);
                     }
                  }
            
                  // calculation of each hessian element
                  int k = indexAtomA*CartesianType_end + axisA; // hessian index, i.e. hessian[k][l]
                  for(int indexAtomB=indexAtomA; indexAtomB<this->molecule->GetNumberAtoms(); indexAtomB++){
                     // hessian element (atomA != atomB)
                     if(indexAtomA!=indexAtomB){
                        const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
                        for(int axisB = XAxis; axisB<CartesianType_end; axisB++){
                           int l = indexAtomB*CartesianType_end + axisB; // hessian index, i.e. hessian[k][l]
                           hessianSCF[k][l] = this->GetHessianElementDifferentAtomsSCF(indexAtomA, 
                                                                                       indexAtomB,
                                                                                       static_cast<CartesianType>(axisA), 
                                                                                       static_cast<CartesianType>(axisB), 
                                                                                       orbitalElectronPopulation,
                                                                                       orbitalElectronPopulation1stDerivs,
                                                                                       diatomicOverlapAOs1stDerivs,
                                                                                       diatomicOverlapAOs2ndDerivs,
                                                                                       diatomicTwoElecTwoCore1stDerivs,
                                                                                       diatomicTwoElecTwoCore2ndDerivs);
                           if(isMassWeighted){
                              hessianSCF[k][l] /= sqrt(atomA.GetCoreMass()*atomB.GetCoreMass());
                           }
                        }
                     }
                     // hessian element (atomA == atomB)
                     else{
                        for(int axisA2 = axisA; axisA2<CartesianType_end; axisA2++){
                           int l = indexAtomA*CartesianType_end + axisA2; // hessian index, i.e. hessian[k][l]
                           hessianSCF[k][l] = this->GetHessianElementSameAtomsSCF(indexAtomA, 
                                                                                  static_cast<CartesianType>(axisA), 
                                                                                  static_cast<CartesianType>(axisA2), 
                                                                                  orbitalElectronPopulation,
                                                                                  orbitalElectronPopulation1stDerivs,
                                                                                  diatomicOverlapAOs1stDerivs,
                                                                                  diatomicOverlapAOs2ndDerivs,
                                                                                  diatomicTwoElecTwoCore1stDerivs,
                                                                                  diatomicTwoElecTwoCore2ndDerivs);
                           if(isMassWeighted){
                              hessianSCF[k][l] /= atomA.GetCoreMass();
                           }
                        }
                     }
                  }
            
               }
            }
         }
         catch(MolDSException ex){
#pragma omp critical
            ex.Serialize(ompErrors);
         }
         this->FreeTempMatricesEachThreadCalcHessianSCF(&diatomicOverlapAOs1stDerivs,
                                                        &diatomicOverlapAOs2ndDerivs,
                                                        &diatomicTwoElecTwoCore1stDerivs, 
                                                        &diatomicTwoElecTwoCore2ndDerivs);
      }// end of omp-region
      // Exception throwing for omp-region
      if(!ompErrors.str().empty()){
         throw MolDSException::Deserialize(ompErrors);
      }

   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&orbitalElectronPopulation1stDerivs, 
                                                 totalNumberAOs,
                                                 totalNumberAOs,
                                                 this->molecule->GetNumberAtoms(),
                                                 CartesianType_end);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&orbitalElectronPopulation1stDerivs, 
                                              totalNumberAOs,
                                              totalNumberAOs,
                                              this->molecule->GetNumberAtoms(),
                                              CartesianType_end);
   int hessianDim = this->molecule->GetNumberAtoms()*CartesianType_end;
   for(int k=0; k<hessianDim; k++){
      for(int l=k; l<hessianDim; l++){
         hessianSCF[l][k] = hessianSCF[k][l];
      }
   }
   /*
   int hessianDim = this->molecule->GetNumberAtoms()*CartesianType_end;
   for(int i=0; i<hessianDim; i++){
      for(int j=0; j<hessianDim; j++){
         printf("hess elem: %d %d %e\n",i,j,hessianSCF[i][j]);
         printf("%e ",hessianSCF[i][j]);
      }
      cout << endl;
   }
   cout << endl << endl;
   */
}

void Mndo::CalcOrbitalElectronPopulation1stDerivatives(double**** orbitalElectronPopulation1stDerivs) const{
   int totalNumberAOs = this->molecule->GetTotalNumberAOs();
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int numberVir = this->molecule->GetTotalNumberAOs() - numberOcc;
   vector<MoIndexPair> nonRedundantQIndeces;
   vector<MoIndexPair> redundantQIndeces;
   this->CalcActiveSetVariablesQ(&nonRedundantQIndeces, &redundantQIndeces, numberOcc, numberVir);
   int dimensionCPHF = nonRedundantQIndeces.size() + redundantQIndeces.size();
   int numberCPHFs = this->molecule->GetNumberAtoms()*CartesianType_end;
   double** solutionsCPHF = NULL; // solutions of CPHF
   double** transposedFockMatrix = NULL; // transposed Fock matrix
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&solutionsCPHF, numberCPHFs, dimensionCPHF);
      MallocerFreer::GetInstance()->Malloc<double>(&transposedFockMatrix, totalNumberAOs, totalNumberAOs);
      this->SolveCPHF(solutionsCPHF, nonRedundantQIndeces, redundantQIndeces);
      this->TransposeFockMatrixMatrix(transposedFockMatrix);
      stringstream ompErrors;
#pragma omp parallel for schedule(auto)
      for(int mu=0; mu<totalNumberAOs; mu++){
         try{
            for(int nu=0; nu<totalNumberAOs; nu++){
               for(int indexAtomA=0; indexAtomA<this->molecule->GetNumberAtoms(); indexAtomA++){
                  for(int axis=XAxis; axis<CartesianType_end; axis++){
         
                     int moI, moJ;
                     double nI, nJ;
                     int indexSolutionCPHF = indexAtomA*CartesianType_end+axis;
                     orbitalElectronPopulation1stDerivs[mu][nu][indexAtomA][axis] = 0.0;
                     for(int k=0; k<nonRedundantQIndeces.size(); k++){
                        moI = nonRedundantQIndeces[k].moI;
                        moJ = nonRedundantQIndeces[k].moJ;
                        nI = moI<numberOcc ? 2.0 : 0.0;
                        nJ = moJ<numberOcc ? 2.0 : 0.0;
                        orbitalElectronPopulation1stDerivs[mu][nu][indexAtomA][axis]
                           += (nJ-nI)*
                              (transposedFockMatrix[mu][moJ]*transposedFockMatrix[nu][moI]+
                               transposedFockMatrix[mu][moI]*transposedFockMatrix[nu][moJ])*
                              solutionsCPHF[indexSolutionCPHF][k];
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
      // check the CPHF's solutions 
      for(int indexAtomA=0; indexAtomA<this->molecule->GetNumberAtoms(); indexAtomA++){
         for(int axis=XAxis; axis<CartesianType_end; axis++){
            double temp=0.0;
            printf("hoge-cphf: atom:%d axis:%s start\n ",indexAtomA,CartesianTypeStr(axis));
            for(int mu=0; mu<totalNumberAOs; mu++){
               temp += orbitalElectronPopulation1stDerivs[mu][mu][indexAtomA][axis];
               printf("%e\n",orbitalElectronPopulation1stDerivs[mu][mu][indexAtomA][axis]);
            }
            printf("hoge-cphf: atom:%d axis:%s %e\n\n",indexAtomA,CartesianTypeStr(axis),temp);
         }
      }
      */
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&solutionsCPHF, numberCPHFs, dimensionCPHF);
      MallocerFreer::GetInstance()->Free<double>(&transposedFockMatrix, totalNumberAOs, totalNumberAOs);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&solutionsCPHF, numberCPHFs, dimensionCPHF);
   MallocerFreer::GetInstance()->Free<double>(&transposedFockMatrix, totalNumberAOs, totalNumberAOs);
}

// Solve CPHF (34) in [PT_1996].
// Derivative coordinates is "axis" of atomA.
// The solution of the CPHF is set to solution.
// solutionsCPHF[i][j] is the j-th element of i-th CPHF solution.
void Mndo::SolveCPHF(double** solutionsCPHF,
                     const vector<MoIndexPair>& nonRedundantQIndeces,
                     const vector<MoIndexPair>& redundantQIndeces) const{
   int dimensionCPHF = nonRedundantQIndeces.size() + redundantQIndeces.size();
   int numberCPHFs = this->molecule->GetNumberAtoms()*CartesianType_end;
   double** matrixCPHF = NULL; // (Gmamma - K matrix)N, see (40) - (46) to slove (34) in [PT_1996].
   try{
      this->MallocTempMatricesSolveCPHF(&matrixCPHF, dimensionCPHF);
      this->CalcMatrixCPHF(matrixCPHF, nonRedundantQIndeces, redundantQIndeces);
      // Static first order focks are temporary stored in solutionsCPHF.
      // This focks in solutionsCPHF are overwritten with solutions of the CPHF by Lapack.
      this->CalcStaticFirstOrderFocks(solutionsCPHF, nonRedundantQIndeces,redundantQIndeces);
      MolDS_wrappers::Lapack::GetInstance()->Dgetrs(matrixCPHF, solutionsCPHF, dimensionCPHF,numberCPHFs);
   }
   catch(MolDSException ex){
      this->FreeTempMatricesSolveCPHF(&matrixCPHF, dimensionCPHF);
      throw ex;
   }
   this->FreeTempMatricesSolveCPHF(&matrixCPHF, dimensionCPHF);
}

// clac right side hands of CPHF, (34) in [PT_1996]
void Mndo::CalcStaticFirstOrderFocks(double** staticFirstOrderFocks,
                                     const vector<MoIndexPair>& nonRedundantQIndeces,
                                     const vector<MoIndexPair>& redundantQIndeces) const{
   stringstream ompErrors;
#pragma omp parallel for schedule(auto)
   for(int indexAtomA=0; indexAtomA<this->molecule->GetNumberAtoms(); indexAtomA++){
      try{
         for(int axisA=XAxis; axisA<CartesianType_end; axisA++){
            int k=indexAtomA*CartesianType_end + axisA;
            this->CalcStaticFirstOrderFock(staticFirstOrderFocks[k], 
                                           nonRedundantQIndeces,
                                           redundantQIndeces,
                                           indexAtomA,
                                           static_cast<CartesianType>(axisA));
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

// clac right side hand of CPHF, (34) in [PT_1996]
// Derivative coordinates is "axisA" of atomA.
void Mndo::CalcStaticFirstOrderFock(double* staticFirstOrderFock,
                                    const vector<MoIndexPair>& nonRedundantQIndeces,
                                    const vector<MoIndexPair>& redundantQIndeces,
                                    int indexAtomA,
                                    CartesianType axisA) const{
   MallocerFreer::GetInstance()->Initialize<double>(staticFirstOrderFock,
                                                    nonRedundantQIndeces.size()+redundantQIndeces.size());
   double***** diatomicTwoElecTwoCore1stDerivs = NULL;
   double***   diatomicOverlapAOs1stDerivs = NULL;
   
   try{
      this->MallocTempMatricesStaticFirstOrderFock(&diatomicTwoElecTwoCore1stDerivs, &diatomicOverlapAOs1stDerivs);
      const Atom& atomA = *molecule->GetAtom(indexAtomA);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();
      int coreChargeA   = atomA.GetCoreCharge();
      for(int indexAtomB=0; indexAtomB<this->molecule->GetNumberAtoms(); indexAtomB++){
         if(indexAtomA != indexAtomB){
            const Atom& atomB = *molecule->GetAtom(indexAtomB);
            int firstAOIndexB = atomB.GetFirstAOIndex();
            int lastAOIndexB  = atomB.GetLastAOIndex();
            int coreChargeB   = atomB.GetCoreCharge();

            // calc. first derivative of two elec two core interaction
            this->CalcDiatomicTwoElecTwoCore1stDerivatives(diatomicTwoElecTwoCore1stDerivs, indexAtomA, indexAtomB);
            // calc. first derivative of overlapAOs.
            this->CalcDiatomicOverlapAOs1stDerivatives(diatomicOverlapAOs1stDerivs, atomA, atomB);

            for(int i=0; i<nonRedundantQIndeces.size()+redundantQIndeces.size();i++){
               int moI=0, moJ=0;
               if(i<nonRedundantQIndeces.size()){
                  moI = nonRedundantQIndeces[i].moI;
                  moJ = nonRedundantQIndeces[i].moJ;
               }
               else{
                  moI = redundantQIndeces[i-nonRedundantQIndeces.size()].moI;
                  moJ = redundantQIndeces[i-nonRedundantQIndeces.size()].moJ;
               }
               // calc. static first order Fock;
               for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
                  for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
                     for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                        for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){

                           double temp1 = this->fockMatrix[moI][mu]
                                         *this->fockMatrix[moJ][nu]
                                         *this->orbitalElectronPopulation[lambda][sigma]
                                         +this->fockMatrix[moI][lambda]
                                         *this->fockMatrix[moJ][sigma]
                                         *this->orbitalElectronPopulation[mu][nu]
                                         -0.5
                                         *this->fockMatrix[moI][mu]
                                         *this->fockMatrix[moJ][lambda]
                                         *this->orbitalElectronPopulation[nu][sigma]
                                         -0.5
                                         *this->fockMatrix[moI][lambda]
                                         *this->fockMatrix[moJ][mu]
                                         *this->orbitalElectronPopulation[nu][sigma];
                           staticFirstOrderFock[i] += temp1*diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA]
                                                                                           [nu-firstAOIndexA]
                                                                                           [lambda-firstAOIndexB]
                                                                                           [sigma-firstAOIndexB]
                                                                                           [axisA];
                        } //sigma-loop
                     } // lambda-loop

                     double temp2 = this->fockMatrix[moI][mu]
                                   *this->fockMatrix[moJ][nu]
                                   *coreChargeB
                                   *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA]
                                                                   [nu-firstAOIndexA]
                                                                   [s]
                                                                   [s]
                                                                   [axisA];
                     staticFirstOrderFock[i] -= temp2;

                  } // nu-loop
               } // mu-loop

               for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                  for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
                     
                     double temp3 = this->fockMatrix[moI][lambda]
                                   *this->fockMatrix[moJ][sigma]
                                   *coreChargeA
                                   *diatomicTwoElecTwoCore1stDerivs[s]
                                                                   [s]
                                                                   [lambda-firstAOIndexB]
                                                                   [sigma-firstAOIndexB]
                                                                   [axisA];
                     staticFirstOrderFock[i] -= temp3;
            
                  } //sigma-loop
               } // lambda-loop
            
               for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
                  for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                     double bondParameter = 0.5*(atomA.GetBondingParameter(this->theory, 
                                                                           atomA.GetValence(mu-firstAOIndexA)) 
                                                +atomB.GetBondingParameter(this->theory, 
                                                                           atomB.GetValence(lambda-firstAOIndexB))); 
                     double temp4 = ( this->fockMatrix[moI][mu]
                                     *this->fockMatrix[moJ][lambda]
                                     +this->fockMatrix[moI][lambda]
                                     *this->fockMatrix[moJ][mu]
                                    )
                                    *bondParameter
                                    *diatomicOverlapAOs1stDerivs[mu-firstAOIndexA][lambda-firstAOIndexB][axisA];
                     staticFirstOrderFock[i] += temp4;
                        
                  } //lambda-loop
               } // mu-loop
            } // i-loop
         }
      }
   }
   catch(MolDSException ex){
      this->FreeTempMatricesStaticFirstOrderFock(&diatomicTwoElecTwoCore1stDerivs, &diatomicOverlapAOs1stDerivs);
      throw ex;
   }
   this->FreeTempMatricesStaticFirstOrderFock(&diatomicTwoElecTwoCore1stDerivs, &diatomicOverlapAOs1stDerivs);

   /*
   printf("staticFirstOrderFock(atomA:%d axis:%s)\n",indexAtomA,CartesianTypeStr(axisA));
   for(int i=0; i<nonRedundantQIndeces.size()+redundantQIndeces.size();i++){
      printf("i:%d %e\n",i,staticFirstOrderFock[i]);
   }
   */
}

void Mndo::MallocTempMatricesStaticFirstOrderFock(double****** diatomicTwoElecTwoCore1stDeriv,
                                                  double**** diatomicOverlapAOs1stDeriv)const{
   MallocerFreer::GetInstance()->Malloc<double>(diatomicTwoElecTwoCore1stDeriv,
                                                dxy,
                                                dxy,
                                                dxy,
                                                dxy,
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicOverlapAOs1stDeriv,
                                                OrbitalType_end, 
                                                OrbitalType_end, 
                                                CartesianType_end);
}

void Mndo::FreeTempMatricesStaticFirstOrderFock(double****** diatomicTwoElecTwoCore1stDeriv,
                                                double****   diatomicOverlapAOs1stDeriv)const{
   MallocerFreer::GetInstance()->Free<double>(diatomicTwoElecTwoCore1stDeriv,
                                              dxy,
                                              dxy,
                                              dxy,
                                              dxy,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicOverlapAOs1stDeriv,
                                              OrbitalType_end, 
                                              OrbitalType_end, 
                                              CartesianType_end);
}

// see (40) - (46) in [PT_1996].
// This method calculates "(\Gamma - K)N" to solve CPHF (34) in [PT_1966]
void Mndo::CalcMatrixCPHF(double** matrixCPHF, 
                          const vector<MoIndexPair>& nonRedundantQIndeces,
                          const vector<MoIndexPair>& redundantQIndeces) const{
   int dimensionCPHF = nonRedundantQIndeces.size() + redundantQIndeces.size();
   double* occupations = NULL;
   MallocerFreer::GetInstance()->Malloc<double>(&occupations, dimensionCPHF);
   stringstream ompErrors;
#pragma omp parallel 
   {
      try{
         // calc diagonal part of N
#pragma omp for schedule(auto)
         for(int i=0; i<dimensionCPHF; i++){
            if(i<nonRedundantQIndeces.size()){
               int moI = nonRedundantQIndeces[i].moI;
               int moJ = nonRedundantQIndeces[i].moJ;
               occupations[i] = this->GetNNRElement(moI, moJ, moI, moJ);
            }
            else{
               int moI = redundantQIndeces[i-nonRedundantQIndeces.size()].moI;
               int moJ = redundantQIndeces[i-nonRedundantQIndeces.size()].moJ;
               occupations[i] = this->GetNRElement(moI, moJ, moI, moJ);
            }
         }

         // calc (\Gamma - K)N
#pragma omp for schedule(auto)
         for(int i=0; i<nonRedundantQIndeces.size(); i++){
            int moI = nonRedundantQIndeces[i].moI;
            int moJ = nonRedundantQIndeces[i].moJ;
            for(int j=0; j<nonRedundantQIndeces.size(); j++){
               int moK = nonRedundantQIndeces[j].moI;
               int moL = nonRedundantQIndeces[j].moJ;
               matrixCPHF[i][j] = (this->GetGammaNRElement(moI, moJ, moK, moL)-this->GetKNRElement(moI, moJ, moK, moL))
                                 *occupations[j];
            }    
         }  
   
#pragma omp for schedule(auto)
         for(int i=nonRedundantQIndeces.size(); i<dimensionCPHF; i++){
            int moI = redundantQIndeces[i-nonRedundantQIndeces.size()].moI;
            int moJ = redundantQIndeces[i-nonRedundantQIndeces.size()].moJ;
            for(int j=0; j<nonRedundantQIndeces.size(); j++){
               int moK = nonRedundantQIndeces[j].moI;
               int moL = nonRedundantQIndeces[j].moJ;
               matrixCPHF[i][j] = -this->GetKRElement(moI, moJ, moK, moL)*occupations[j];
            }
         }
   
#pragma omp for schedule(auto)
         for(int i=nonRedundantQIndeces.size(); i<dimensionCPHF; i++){
            int moI = redundantQIndeces[i-nonRedundantQIndeces.size()].moI;
            int moJ = redundantQIndeces[i-nonRedundantQIndeces.size()].moJ;
            matrixCPHF[i][i] = this->GetGammaRElement(moI, moJ, moI, moJ)*occupations[i];
         }
      }
      catch(MolDSException ex){
#pragma omp critical
         ex.Serialize(ompErrors);
      }
   }
   MallocerFreer::GetInstance()->Free<double>(&occupations, dimensionCPHF);
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }

   /*
   printf("matrixCPHF\n");
   for(int i=0; i<nonRedundantQIndeces.size()+redundantQIndeces.size(); i++){
      for(int j=0; j<nonRedundantQIndeces.size()+redundantQIndeces.size(); j++){
         printf("i:%d j:%d %e\n",i,j,matrixCPHF[i][j]);
      }
   }
   */
}

void Mndo::MallocTempMatricesSolveCPHF(double*** matrixCPHF,
                                       int dimensionCPHF) const{
   MallocerFreer::GetInstance()->Malloc<double>(matrixCPHF, dimensionCPHF, dimensionCPHF);
}

void Mndo::FreeTempMatricesSolveCPHF(double*** matrixCPHF,
                                     int dimensionCPHF) const{
   MallocerFreer::GetInstance()->Free<double>(matrixCPHF, dimensionCPHF, dimensionCPHF);
}

void Mndo::CalcForceSCFElecCoreAttractionPart(double* force, 
                                             int indexAtomA, 
                                             int indexAtomB,
                                             double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int lastAOIndexA  = atomA.GetLastAOIndex();
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
         for(int i=0; i<CartesianType_end; i++){
            force[i] += -1.0
                       *this->orbitalElectronPopulation[mu][nu]
                       *this->GetElectronCoreAttraction1stDerivative(indexAtomA, 
                                                                     indexAtomB, 
                                                                     mu-firstAOIndexA, 
                                                                     nu-firstAOIndexA,
                                                                     diatomicTwoElecTwoCore1stDerivs,
                                                                     (CartesianType)i);
         }
      }
   }
}

void Mndo::CalcForceSCFOverlapAOsPart(double* force, 
                                     int indexAtomA, 
                                     int indexAtomB,
                                     double const* const* const* diatomicOverlapAOs1stDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int lastAOIndexA  = atomA.GetLastAOIndex();
   int lastAOIndexB  = atomB.GetLastAOIndex();
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      for(int nu=firstAOIndexB; nu<=lastAOIndexB; nu++){
         double bondParameter = atomA.GetBondingParameter(
                                      this->theory, 
                                      atomA.GetValence(mu-firstAOIndexA)) 
                               +atomB.GetBondingParameter(
                                      this->theory, 
                                      atomB.GetValence(nu-firstAOIndexB)); 
         bondParameter*=0.5;
         for(int i=0; i<CartesianType_end; i++){
            force[i] += -1.0
                       *this->orbitalElectronPopulation[mu][nu]
                       *bondParameter
                       *diatomicOverlapAOs1stDerivs[mu-firstAOIndexA][nu-firstAOIndexB][i];
         }
      }
   }
}

void Mndo::CalcForceSCFTwoElecPart(double* force, 
                                  int indexAtomA, 
                                  int indexAtomB,
                                  double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int lastAOIndexA  = atomA.GetLastAOIndex();
   int lastAOIndexB  = atomB.GetLastAOIndex();
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
         for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
            for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
               for(int i=0; i<CartesianType_end; i++){
                  force[i] -= 0.5
                             *this->orbitalElectronPopulation[mu][nu]
                             *this->orbitalElectronPopulation[lambda][sigma]
                             *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA]
                                                             [nu-firstAOIndexA]
                                                             [lambda-firstAOIndexB]
                                                             [sigma-firstAOIndexB]
                                                             [(CartesianType)i];
                  force[i] += 0.25
                             *this->orbitalElectronPopulation[mu][lambda]
                             *this->orbitalElectronPopulation[nu][sigma]
                             *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA]
                                                             [nu-firstAOIndexA]
                                                             [lambda-firstAOIndexB]
                                                             [sigma-firstAOIndexB]
                                                             [(CartesianType)i];
               }
            }
         }
      }
   }
}

void Mndo::CalcForceExcitedStaticPart(double* force, 
                                      int elecStateIndex,
                                      int indexAtomA, 
                                      int indexAtomB,
                                      double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int lastAOIndexA  = atomA.GetLastAOIndex();
   int lastAOIndexB  = atomB.GetLastAOIndex();
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
         for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
            for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
               for(int i=0; i<CartesianType_end; i++){
                  double temp= 2.0*this->etaMatrixForce[elecStateIndex][mu][nu]
                                  *this->etaMatrixForce[elecStateIndex][lambda][sigma]
                              -1.0*this->etaMatrixForce[elecStateIndex][mu][lambda]
                                  *this->etaMatrixForce[elecStateIndex][nu][sigma];
                  force[i] += temp
                             *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA]
                                                             [nu-firstAOIndexA]
                                                             [lambda-firstAOIndexB]
                                                             [sigma-firstAOIndexB]
                                                             [i];
               }
            }
         }
      }
   }
}

void Mndo::CalcForceExcitedElecCoreAttractionPart(double* force, 
                                                  int elecStateIndex,
                                                  int indexAtomA, 
                                                  int indexAtomB,
                                                  double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int lastAOIndexA  = atomA.GetLastAOIndex();
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
         for(int i=0; i<CartesianType_end; i++){
            force[i] += -1.0
                       *this->zMatrixForce[elecStateIndex][mu][nu]
                       *this->GetElectronCoreAttraction1stDerivative(indexAtomA, 
                                                                     indexAtomB, 
                                                                     mu-firstAOIndexA, 
                                                                     nu-firstAOIndexA,
                                                                     diatomicTwoElecTwoCore1stDerivs,
                                                                     (CartesianType)i);
         }
      }
   }
}

void Mndo::CalcForceExcitedTwoElecPart(double* force, 
                                       int elecStateIndex,
                                       int indexAtomA, 
                                       int indexAtomB,
                                       double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int lastAOIndexA  = atomA.GetLastAOIndex();
   int lastAOIndexB  = atomB.GetLastAOIndex();
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
         for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
            for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
               for(int i=0; i<CartesianType_end; i++){
                  force[i] -= this->zMatrixForce[elecStateIndex][mu][nu]
                             *this->orbitalElectronPopulation[lambda][sigma]
                             *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA]
                                                             [nu-firstAOIndexA]
                                                             [lambda-firstAOIndexB]
                                                             [sigma-firstAOIndexB]
                                                             [i];
                  force[i] += 0.50
                             *this->zMatrixForce[elecStateIndex][mu][lambda]
                             *this->orbitalElectronPopulation[nu][sigma]
                             *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA]
                                                             [nu-firstAOIndexA]
                                                             [lambda-firstAOIndexB]
                                                             [sigma-firstAOIndexB]
                                                             [(CartesianType)i];
               }
            }
         }
      }
   }
}

// electronicStateIndex is index of the electroinc eigen state.
// "electronicStateIndex = 0" means electronic ground state. 
void Mndo::CalcForce(const vector<int>& elecStates){
   this->CheckMatrixForce(elecStates);
   if(this->RequiresExcitedStatesForce(elecStates)){
      this->CalcEtaMatrixForce(elecStates);
      this->CalcZMatrixForce(elecStates);
   }
   stringstream ompErrors;
#pragma omp parallel
   {
      double***** diatomicTwoElecTwoCore1stDerivs = NULL;
      double***   diatomicOverlapAOs1stDerivs = NULL;
      try{
         this->MallocTempMatricesCalcForce(&diatomicOverlapAOs1stDerivs, &diatomicTwoElecTwoCore1stDerivs);
#pragma omp for schedule(auto)
         for(int a=0; a<this->molecule->GetNumberAtoms(); a++){
            const Atom& atomA = *molecule->GetAtom(a);
            int firstAOIndexA = atomA.GetFirstAOIndex();
            int lastAOIndexA  = atomA.GetLastAOIndex();
            for(int b=0; b<this->molecule->GetNumberAtoms(); b++){
               if(a == b){continue;}
               const Atom& atomB = *molecule->GetAtom(b);
               int firstAOIndexB = atomB.GetFirstAOIndex();
               int lastAOIndexB  = atomB.GetLastAOIndex();

               // calc. first derivative of overlapAOs.
               this->CalcDiatomicOverlapAOs1stDerivatives(diatomicOverlapAOs1stDerivs, atomA, atomB);
               // calc. first derivative of two elec two core interaction
               this->CalcDiatomicTwoElecTwoCore1stDerivatives(diatomicTwoElecTwoCore1stDerivs, a, b);

               // core repulsion part
               double coreRepulsion[CartesianType_end] = {0.0,0.0,0.0};
               for(int i=0; i<CartesianType_end; i++){
                  coreRepulsion[i] += this->GetDiatomCoreRepulsion1stDerivative(
                                            a, b, (CartesianType)i);
                  if(Parameters::GetInstance()->RequiresVdWSCF()){
                     coreRepulsion[i] += this->GetDiatomVdWCorrection1stDerivative(
                                               a, b, (CartesianType)i);
                  }
               }  
               // electron core attraction part (ground state)
               double forceElecCoreAttPart[CartesianType_end] = {0.0,0.0,0.0};
               this->CalcForceSCFElecCoreAttractionPart(forceElecCoreAttPart,
                                                        a,
                                                        b,
                                                        diatomicTwoElecTwoCore1stDerivs);
               // overlapAOs part (ground state)
               double forceOverlapAOsPart[CartesianType_end] = {0.0,0.0,0.0};
               this->CalcForceSCFOverlapAOsPart(forceOverlapAOsPart, 
                                                a,
                                                b,
                                                diatomicOverlapAOs1stDerivs);
               // two electron part (ground state)
               double forceTwoElecPart[CartesianType_end] = {0.0,0.0,0.0};
               this->CalcForceSCFTwoElecPart(forceTwoElecPart,
                                             a,
                                             b,
                                             diatomicTwoElecTwoCore1stDerivs);
               // sum up contributions from each part (ground state)
#pragma omp critical
               {
                  for(int n=0; n<elecStates.size(); n++){
                     for(int i=0; i<CartesianType_end; i++){
                        this->matrixForce[n][a][i] -= coreRepulsion[i];
                        this->matrixForce[n][a][i] += forceElecCoreAttPart[i];
                        this->matrixForce[n][a][i] += forceOverlapAOsPart[i];
                        this->matrixForce[n][a][i] += forceTwoElecPart[i];
                        this->matrixForce[n][b][i] -= forceElecCoreAttPart[i];
                        this->matrixForce[n][b][i] -= forceOverlapAOsPart[i];
                        this->matrixForce[n][b][i] -= forceTwoElecPart[i];
                     }
                  }
               }
               // excited state force
               for(int n=0; n<elecStates.size(); n++){
                  if(elecStates[n]<=0){continue;}
                  // static part
                  double forceExcitedStaticPart[CartesianType_end] = {0.0,0.0,0.0};
                  this->CalcForceExcitedStaticPart(forceExcitedStaticPart,
                                                   n,
                                                   a,
                                                   b,
                                                   diatomicTwoElecTwoCore1stDerivs);
                  // sum up contributions from static part (excited state)
#pragma omp critical
                  {
                     for(int i=0; i<CartesianType_end; i++){
                        this->matrixForce[n][b][i] += forceExcitedStaticPart[i];
                        this->matrixForce[n][a][i] -= forceExcitedStaticPart[i];
                     }
                  }

                  // response part
                  // electron core attraction part (excited states)
                  double forceExcitedElecCoreAttPart[CartesianType_end]={0.0,0.0,0.0};
                  this->CalcForceExcitedElecCoreAttractionPart(
                                             forceExcitedElecCoreAttPart,
                                             n,
                                             a,
                                             b,
                                             diatomicTwoElecTwoCore1stDerivs);
                  // overlapAOs part (excited states)
                  double forceExcitedOverlapAOsPart[CartesianType_end] = {0.0,0.0,0.0};
                  this->CalcForceExcitedOverlapAOsPart(forceExcitedOverlapAOsPart, 
                                                       n,
                                                       a,
                                                       b,
                                                       diatomicOverlapAOs1stDerivs);
                  // two electron part (excited states)
                  double forceExcitedTwoElecPart[CartesianType_end] = {0.0,0.0,0.0};
                  this->CalcForceExcitedTwoElecPart(forceExcitedTwoElecPart,
                                                       n,
                                                       a,
                                                       b,
                                                       diatomicTwoElecTwoCore1stDerivs);
                  // sum up contributions from response part (excited state)
#pragma omp critical
                  {
                     for(int i=0; i<CartesianType_end; i++){
                        this->matrixForce[n][a][i] += forceExcitedElecCoreAttPart[i];
                        this->matrixForce[n][a][i] += forceExcitedOverlapAOsPart[i];
                        this->matrixForce[n][a][i] += forceExcitedTwoElecPart[i];
                        this->matrixForce[n][b][i] -= forceExcitedElecCoreAttPart[i];
                        this->matrixForce[n][b][i] -= forceExcitedOverlapAOsPart[i];
                        this->matrixForce[n][b][i] -= forceExcitedTwoElecPart[i];
                     }
                  }
               } // end of excited state force
            }    // end of for(int b)
         }       // end of for(int a)
      }          // end of try
      catch(MolDSException ex){
#pragma omp critical
         ex.Serialize(ompErrors);
      }
      this->FreeTempMatricesCalcForce(&diatomicOverlapAOs1stDerivs, &diatomicTwoElecTwoCore1stDerivs);
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
}

void Mndo::MallocTempMatricesCalcForce(double**** diatomicOverlapAOs1stDerivs, double****** diatomicTwoElecTwoCore1stDerivs) const{
   MallocerFreer::GetInstance()->Malloc<double>(diatomicOverlapAOs1stDerivs, 
                                                OrbitalType_end,
                                                OrbitalType_end,
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicTwoElecTwoCore1stDerivs,
                                                dxy,
                                                dxy,
                                                dxy,
                                                dxy,
                                                CartesianType_end);
}

void Mndo::FreeTempMatricesCalcForce(double**** diatomicOverlapAOs1stDerivs, double****** diatomicTwoElecTwoCore1stDerivs) const{
   MallocerFreer::GetInstance()->Free<double>(diatomicOverlapAOs1stDerivs, 
                                              OrbitalType_end,
                                              OrbitalType_end,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicTwoElecTwoCore1stDerivs,
                                              dxy,
                                              dxy,
                                              dxy,
                                              dxy,
                                              CartesianType_end);
}

// see (18) in [PT_1997]
double Mndo::GetSmallQElement(int moI, 
                              int moP, 
                              double const* const* xiOcc, 
                              double const* const* xiVir, 
                              double const* const* eta) const{
   double value = 0.0;
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   bool isMoPOcc = moP<numberOcc ? true : false;
   
   for(int A=0; A<molecule->GetNumberAtoms(); A++){
      const Atom& atomA = *molecule->GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();

      for(int B=A; B<molecule->GetNumberAtoms(); B++){
         const Atom& atomB = *molecule->GetAtom(B);
         int firstAOIndexB = atomB.GetFirstAOIndex();
         int lastAOIndexB  = atomB.GetLastAOIndex();

         if(A!=B){
            for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
               for(int nu=mu; nu<=lastAOIndexA; nu++){
                  for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                     for(int sigma=lambda; sigma<=lastAOIndexB; sigma++){
                        double twoElecInt = 0.0;
                        twoElecInt = this->twoElecTwoCore[A]
                                                         [B]
                                                         [mu-firstAOIndexA]
                                                         [nu-firstAOIndexA]
                                                         [lambda-firstAOIndexB]
                                                         [sigma-firstAOIndexB];
                        double temp = 0.0;
                        if(isMoPOcc){
                           int p = numberOcc - (moP+1);
                           temp = 4.0*xiOcc[p][nu]*eta[lambda][sigma]
                                 -1.0*xiOcc[p][lambda]*eta[nu][sigma]
                                 -1.0*xiOcc[p][sigma]*eta[nu][lambda];
                           value += twoElecInt*this->fockMatrix[moI][mu]*temp;
                           temp = 4.0*xiOcc[p][sigma]*eta[mu][nu]
                                 -1.0*xiOcc[p][mu]*eta[sigma][nu]
                                 -1.0*xiOcc[p][nu]*eta[sigma][mu];
                           value += twoElecInt*this->fockMatrix[moI][lambda]*temp;
                        }
                        else{
                           int p = moP - numberOcc;
                           temp = 4.0*xiVir[p][nu]*eta[lambda][sigma]
                                 -1.0*xiVir[p][lambda]*eta[sigma][nu]
                                 -1.0*xiVir[p][sigma]*eta[lambda][nu];
                           value += twoElecInt*this->fockMatrix[moI][mu]*temp;
                           temp = 4.0*xiVir[p][sigma]*eta[mu][nu]
                                 -1.0*xiVir[p][mu]*eta[nu][sigma]
                                 -1.0*xiVir[p][nu]*eta[mu][sigma];
                           value += twoElecInt*this->fockMatrix[moI][lambda]*temp;
                        }
                         
                        if(lambda!=sigma){
                           if(isMoPOcc){
                              int p = numberOcc - (moP+1);
                              temp = 4.0*xiOcc[p][nu]*eta[sigma][lambda]
                                    -1.0*xiOcc[p][sigma]*eta[nu][lambda]
                                    -1.0*xiOcc[p][lambda]*eta[nu][sigma];
                              value += twoElecInt*this->fockMatrix[moI][mu]*temp;
                              temp = 4.0*xiOcc[p][lambda]*eta[mu][nu]
                                    -1.0*xiOcc[p][mu]*eta[lambda][nu]
                                    -1.0*xiOcc[p][nu]*eta[lambda][mu];
                              value += twoElecInt*this->fockMatrix[moI][sigma]*temp;
                           }
                           else{
                              int p = moP - numberOcc;
                              temp = 4.0*xiVir[p][nu]*eta[sigma][lambda]
                                    -1.0*xiVir[p][sigma]*eta[lambda][nu]
                                    -1.0*xiVir[p][lambda]*eta[sigma][nu];
                              value += twoElecInt*this->fockMatrix[moI][mu]*temp;
                              temp = 4.0*xiVir[p][lambda]*eta[mu][nu]
                                    -1.0*xiVir[p][mu]*eta[nu][lambda]
                                    -1.0*xiVir[p][nu]*eta[mu][lambda];
                              value += twoElecInt*this->fockMatrix[moI][sigma]*temp;
                           }
                        }
                        
                        if(mu!=nu){
                           if(isMoPOcc){
                              int p = numberOcc - (moP+1);
                              temp = 4.0*xiOcc[p][mu]*eta[lambda][sigma]
                                    -1.0*xiOcc[p][lambda]*eta[mu][sigma]
                                    -1.0*xiOcc[p][sigma]*eta[mu][lambda];
                              value += twoElecInt*this->fockMatrix[moI][nu]*temp;
                              temp = 4.0*xiOcc[p][sigma]*eta[nu][mu]
                                    -1.0*xiOcc[p][nu]*eta[sigma][mu]
                                    -1.0*xiOcc[p][mu]*eta[sigma][nu];
                              value += twoElecInt*this->fockMatrix[moI][lambda]*temp;
                           }
                           else{
                              int p = moP - numberOcc;
                              temp = 4.0*xiVir[p][mu]*eta[lambda][sigma]
                                    -1.0*xiVir[p][lambda]*eta[sigma][mu]
                                    -1.0*xiVir[p][sigma]*eta[lambda][mu];
                              value += twoElecInt*this->fockMatrix[moI][nu]*temp;
                              temp = 4.0*xiVir[p][sigma]*eta[nu][mu]
                                    -1.0*xiVir[p][nu]*eta[mu][sigma]
                                    -1.0*xiVir[p][mu]*eta[nu][sigma];
                              value += twoElecInt*this->fockMatrix[moI][lambda]*temp;
                           }
                        }

                        if(mu!=nu && lambda!=sigma){
                           if(isMoPOcc){
                              int p = numberOcc - (moP+1);
                              temp = 4.0*xiOcc[p][mu]*eta[sigma][lambda]
                                    -1.0*xiOcc[p][sigma]*eta[mu][lambda]
                                    -1.0*xiOcc[p][lambda]*eta[mu][sigma];
                              value += twoElecInt*this->fockMatrix[moI][nu]*temp;
                              temp = 4.0*xiOcc[p][lambda]*eta[nu][mu]
                                    -1.0*xiOcc[p][nu]*eta[lambda][mu]
                                    -1.0*xiOcc[p][mu]*eta[lambda][nu];
                              value += twoElecInt*this->fockMatrix[moI][sigma]*temp;
                           }
                           else{
                              int p = moP - numberOcc;
                              temp = 4.0*xiVir[p][mu]*eta[sigma][lambda]
                                    -1.0*xiVir[p][sigma]*eta[lambda][mu]
                                    -1.0*xiVir[p][lambda]*eta[sigma][mu];
                              value += twoElecInt*this->fockMatrix[moI][nu]*temp;
                              temp = 4.0*xiVir[p][lambda]*eta[nu][mu]
                                    -1.0*xiVir[p][nu]*eta[mu][lambda]
                                    -1.0*xiVir[p][mu]*eta[nu][lambda];
                              value += twoElecInt*this->fockMatrix[moI][sigma]*temp;
                           }
                        }
                        
                     }
                  }
               }
            }
         }
         else{
            for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
               for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
                  for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                     for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
                        double twoElecInt = 0.0;
                        if(mu==nu && lambda==sigma){
                           OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                           OrbitalType orbitalLambda = atomB.GetValence(lambda-firstAOIndexB);
                           twoElecInt = this->GetCoulombInt(orbitalMu, 
                                                            orbitalLambda, 
                                                            atomA);
                        }
                        else if((mu==lambda && nu==sigma) || (nu==lambda && mu==sigma) ){
                           OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                           OrbitalType orbitalNu = atomA.GetValence(nu-firstAOIndexA);
                           twoElecInt = this->GetExchangeInt(orbitalMu, 
                                                             orbitalNu, 
                                                             atomA);
                        }
                        else{
                           twoElecInt = 0.0;
                        }

                        double temp = 0.0;
                        if(isMoPOcc){
                           int p = numberOcc - (moP+1);
                           temp = 4.0*xiOcc[p][nu]*eta[lambda][sigma]
                                 -1.0*xiOcc[p][lambda]*eta[nu][sigma]
                                 -1.0*xiOcc[p][sigma]*eta[nu][lambda];
                        }
                        else{
                           int p = moP - numberOcc;
                           temp = 4.0*xiVir[p][nu]*eta[lambda][sigma]
                                 -1.0*xiVir[p][lambda]*eta[sigma][nu]
                                 -1.0*xiVir[p][sigma]*eta[lambda][nu];
                        }
                        value += twoElecInt*this->fockMatrix[moI][mu]*temp;
                     }  
                  }
               }
            }
         }
      }
   }
   return value;
}

// see common term in eqs. (45) and (46) in [PT_1996],
// that is, 4.0(ij|kl) - (ik|jl) - (il|jk).
double Mndo::GetAuxiliaryKNRKRElement(int moI, int moJ, int moK, int moL) const{
   double value = 0.0;

   // Fast algorith, but this is not easy to read. 
   // Slow algorithm is alos written below.
   for(int A=0; A<this->molecule->GetNumberAtoms(); A++){
      const Atom& atomA = *this->molecule->GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();

      for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
         int muOffSet = mu - firstAOIndexA;
         for(int nu=mu; nu<=lastAOIndexA; nu++){
            int nuOffSet = nu - firstAOIndexA;
            double tmpMN01 = 0.0, tmpMN02 = 0.0, tmpMN03 = 0.0, 
                   tmpMN04 = 0.0, tmpMN05 = 0.0, tmpMN06 = 0.0, 
                   tmpMN13 = 0.0, tmpMN14 = 0.0, tmpMN15 = 0.0, 
                   tmpMN16 = 0.0, tmpMN17 = 0.0, tmpMN18 = 0.0;
            tmpMN01 = 4.0
                     *this->fockMatrix[moI][mu]
                     *this->fockMatrix[moJ][nu];
            tmpMN02 = 4.0
                     *this->fockMatrix[moK][mu]
                     *this->fockMatrix[moL][nu];
            tmpMN03 = this->fockMatrix[moI][mu]
                     *this->fockMatrix[moK][nu];
            tmpMN04 = this->fockMatrix[moJ][mu]
                     *this->fockMatrix[moL][nu];
            tmpMN05 = this->fockMatrix[moI][mu]
                     *this->fockMatrix[moL][nu];
            tmpMN06 = this->fockMatrix[moJ][mu]
                     *this->fockMatrix[moK][nu];
            if(mu != nu){
               tmpMN13 = 4.0
                        *this->fockMatrix[moI][nu]
                        *this->fockMatrix[moJ][mu];
               tmpMN14 = 4.0
                        *this->fockMatrix[moK][nu]
                        *this->fockMatrix[moL][mu];
               tmpMN15 = this->fockMatrix[moI][nu]
                        *this->fockMatrix[moK][mu];
               tmpMN16 = this->fockMatrix[moJ][nu]
                        *this->fockMatrix[moL][mu];
               tmpMN17 = this->fockMatrix[moI][nu]
                        *this->fockMatrix[moL][mu];
               tmpMN18 = this->fockMatrix[moJ][nu]
                        *this->fockMatrix[moK][mu];
            }

            for(int B=A; B<this->molecule->GetNumberAtoms(); B++){
               const Atom& atomB = *this->molecule->GetAtom(B);
               int firstAOIndexB = atomB.GetFirstAOIndex();
               int lastAOIndexB  = atomB.GetLastAOIndex();

               for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                  int lambdaOffSet = lambda - firstAOIndexB;
                  double tmpMNL01 = 0.0, tmpMNL02 = 0.0, tmpMNL03 = 0.0, tmpMNL04 = 0.0, 
                         tmpMNL05 = 0.0, tmpMNL06 = 0.0, tmpMNL07 = 0.0, tmpMNL08 = 0.0, 
                         tmpMNL09 = 0.0, tmpMNL10 = 0.0, tmpMNL11 = 0.0, tmpMNL12 = 0.0, 
                         tmpMNL13 = 0.0, tmpMNL14 = 0.0, tmpMNL15 = 0.0, tmpMNL16 = 0.0, 
                         tmpMNL17 = 0.0, tmpMNL18 = 0.0, tmpMNL19 = 0.0, tmpMNL20 = 0.0,
                         tmpMNL21 = 0.0, tmpMNL22 = 0.0, tmpMNL23 = 0.0, tmpMNL24 = 0.0;
                  tmpMNL01 = tmpMN01*this->fockMatrix[moK][lambda];
                  tmpMNL02 = tmpMN02*this->fockMatrix[moI][lambda];
                  tmpMNL03 = tmpMN03*this->fockMatrix[moJ][lambda];
                  tmpMNL04 = tmpMN04*this->fockMatrix[moI][lambda];
                  tmpMNL05 = tmpMN05*this->fockMatrix[moJ][lambda];
                  tmpMNL06 = tmpMN06*this->fockMatrix[moI][lambda];
                  tmpMNL07 = tmpMN01*this->fockMatrix[moL][lambda];
                  tmpMNL08 = tmpMN02*this->fockMatrix[moJ][lambda];
                  tmpMNL09 = tmpMN03*this->fockMatrix[moL][lambda];
                  tmpMNL10 = tmpMN04*this->fockMatrix[moK][lambda];
                  tmpMNL11 = tmpMN05*this->fockMatrix[moK][lambda];
                  tmpMNL12 = tmpMN06*this->fockMatrix[moL][lambda];
                  tmpMNL01 -= tmpMNL03 + tmpMNL06;
                  tmpMNL04 += tmpMNL05;
                  tmpMNL08 -= tmpMNL10 + tmpMNL12;
                  tmpMNL09 += tmpMNL11;
                  if(mu != nu){
                     tmpMNL13 = tmpMN13*this->fockMatrix[moK][lambda];
                     tmpMNL14 = tmpMN14*this->fockMatrix[moI][lambda];
                     tmpMNL15 = tmpMN15*this->fockMatrix[moJ][lambda];
                     tmpMNL16 = tmpMN16*this->fockMatrix[moI][lambda];
                     tmpMNL17 = tmpMN17*this->fockMatrix[moJ][lambda];
                     tmpMNL18 = tmpMN18*this->fockMatrix[moI][lambda];
                     tmpMNL19 = tmpMN13*this->fockMatrix[moL][lambda];
                     tmpMNL20 = tmpMN14*this->fockMatrix[moJ][lambda];
                     tmpMNL21 = tmpMN15*this->fockMatrix[moL][lambda];
                     tmpMNL22 = tmpMN16*this->fockMatrix[moK][lambda];
                     tmpMNL23 = tmpMN17*this->fockMatrix[moK][lambda];
                     tmpMNL24 = tmpMN18*this->fockMatrix[moL][lambda];
                     tmpMNL13 -= tmpMNL15 + tmpMNL18;
                     tmpMNL16 += tmpMNL17;
                     tmpMNL20 -= tmpMNL22 + tmpMNL24;
                     tmpMNL21 += tmpMNL23;
                     tmpMNL01 += tmpMNL13;
                     tmpMNL02 += tmpMNL14;
                     tmpMNL04 += tmpMNL16;
                     tmpMNL07 += tmpMNL19;
                     tmpMNL08 += tmpMNL20;
                     tmpMNL09 += tmpMNL21;
                  }
                  for(int sigma=lambda; sigma<=lastAOIndexB; sigma++){
                     int sigmaOffSet = sigma - firstAOIndexB;
                     double tmpValue = 0.0;
                     tmpValue += tmpMNL01*this->fockMatrix[moL][sigma];
                     tmpValue += tmpMNL02*this->fockMatrix[moJ][sigma];
                     tmpValue -= tmpMNL04*this->fockMatrix[moK][sigma];
                     if(lambda != sigma){
                        tmpValue += tmpMNL07*this->fockMatrix[moK][sigma];
                        tmpValue += tmpMNL08*this->fockMatrix[moI][sigma];
                        tmpValue -= tmpMNL09*this->fockMatrix[moJ][sigma];
                     }
                     double gamma = 0.0;
                     if(A!=B){
                        gamma = this->twoElecTwoCore[A][B][muOffSet][nuOffSet][lambdaOffSet][sigmaOffSet];
                     }
                     else{
                        if(mu==nu && lambda==sigma){
                           OrbitalType orbitalMu     = atomA.GetValence(muOffSet);
                           OrbitalType orbitalLambda = atomA.GetValence(lambdaOffSet);
                           gamma = this->GetCoulombInt(orbitalMu, orbitalLambda, atomA);
                        }
                        else if((mu==lambda && nu==sigma) || (nu==lambda && mu==sigma) ){
                           OrbitalType orbitalMu = atomA.GetValence(muOffSet);
                           OrbitalType orbitalNu = atomA.GetValence(nuOffSet);
                           gamma = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
                        }
                        else{
                           gamma = 0.0;
                        }
                        gamma *= 0.5;
                     }
                     value += tmpValue*gamma;
                  }
               }
            }
         }
      }
   }
   // End of the fast algorith.

   /*
   // Algorithm using blas
   double** twoElec = NULL;
   double*  twiceMoIJ = NULL;
   double*  twiceMoIK = NULL;
   double*  twiceMoIL = NULL;
   double*  twiceMoKL = NULL;
   double*  twiceMoJL = NULL;
   double*  twiceMoJK = NULL;
   double*  tmpVector = NULL;
   int numAOs = this->molecule->GetTotalNumberAOs();
   MallocerFreer::GetInstance()->Malloc<double>(&twoElec,   this->molecule->GetNumberAtoms()*dxy*dxy, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&twiceMoIJ, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&twiceMoIK, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&twiceMoIL, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&twiceMoKL, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&twiceMoJL, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&twiceMoJK, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&tmpVector, this->molecule->GetNumberAtoms()*dxy*dxy);
   for(int A=0; A<this->molecule->GetNumberAtoms(); A++){
      const Atom& atomA = *this->molecule->GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();
      for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
         for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
            twiceMoIJ[A*dxy*dxy+(mu    -firstAOIndexA)*dxy+(nu   -firstAOIndexA)]=fockMatrix[moI][mu    ]*fockMatrix[moJ][nu   ];
            twiceMoIK[A*dxy*dxy+(mu    -firstAOIndexA)*dxy+(nu   -firstAOIndexA)]=fockMatrix[moI][mu    ]*fockMatrix[moK][nu   ];
            twiceMoIL[A*dxy*dxy+(mu    -firstAOIndexA)*dxy+(nu   -firstAOIndexA)]=fockMatrix[moI][mu    ]*fockMatrix[moL][nu   ];
         }
      }
   }

   for(int B=0; B<this->molecule->GetNumberAtoms(); B++){
      const Atom& atomB = *this->molecule->GetAtom(B);
      int firstAOIndexB = atomB.GetFirstAOIndex();
      int lastAOIndexB  = atomB.GetLastAOIndex();
      for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
         for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
            twiceMoKL[B*dxy*dxy+(lambda-firstAOIndexB)*dxy+(sigma-firstAOIndexB)]=fockMatrix[moK][lambda]*fockMatrix[moL][sigma];
            twiceMoJL[B*dxy*dxy+(lambda-firstAOIndexB)*dxy+(sigma-firstAOIndexB)]=fockMatrix[moJ][lambda]*fockMatrix[moL][sigma];
            twiceMoJK[B*dxy*dxy+(lambda-firstAOIndexB)*dxy+(sigma-firstAOIndexB)]=fockMatrix[moJ][lambda]*fockMatrix[moK][sigma];
         }
      }
   }

   for(int A=0; A<this->molecule->GetNumberAtoms(); A++){
      const Atom& atomA = *this->molecule->GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();
      for(int B=A; B<this->molecule->GetNumberAtoms(); B++){
         const Atom& atomB = *this->molecule->GetAtom(B);
         int firstAOIndexB = atomB.GetFirstAOIndex();
         int lastAOIndexB  = atomB.GetLastAOIndex();
         double gamma = 0.0;
         if(A!=B){
            for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
               for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
                  for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                     for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
                        twoElec[A*dxy*dxy+(mu-firstAOIndexA)*dxy+(nu-firstAOIndexA)]
                               [B*dxy*dxy+(lambda-firstAOIndexB)*dxy+(sigma-firstAOIndexB)] = 
                            this->twoElecTwoCore[A]
                                                [B]
                                                [mu-firstAOIndexA]
                                                [nu-firstAOIndexA]
                                                [lambda-firstAOIndexB]
                                                [sigma-firstAOIndexB];
                     }
                  }
               }
            }
         }
         else{
            for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
               for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
                  for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                     for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
                        if(mu==nu && lambda==sigma){
                           OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                           OrbitalType orbitalLambda = atomB.GetValence(lambda-firstAOIndexB);
                           gamma = this->GetCoulombInt(orbitalMu, orbitalLambda, atomA);
                        }
                        else if((mu==lambda && nu==sigma) || (nu==lambda && mu==sigma) ){
                           OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                           OrbitalType orbitalNu = atomA.GetValence(nu-firstAOIndexA);
                           gamma = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
                        }
                        else{
                           gamma = 0.0;
                        }
                        twoElec[A*dxy*dxy+(mu-firstAOIndexA)*dxy+(nu-firstAOIndexA)]
                               [B*dxy*dxy+(lambda-firstAOIndexB)*dxy+(sigma-firstAOIndexB)] = gamma;
                     }  
                  }
               }
            }
         }
      }
   }
   MolDS_wrappers::Blas::GetInstance()->Dsymv(this->molecule->GetNumberAtoms()*dxy*dxy, 
                                              twoElec, 
                                              twiceMoKL,
                                              tmpVector);
   value = 4.0*MolDS_wrappers::Blas::GetInstance()->Ddot(this->molecule->GetNumberAtoms()*dxy*dxy,twiceMoIJ, tmpVector);
   MolDS_wrappers::Blas::GetInstance()->Dsymv(this->molecule->GetNumberAtoms()*dxy*dxy, 
                                              twoElec, 
                                              twiceMoJL,
                                              tmpVector);
   value -= MolDS_wrappers::Blas::GetInstance()->Ddot(this->molecule->GetNumberAtoms()*dxy*dxy,twiceMoIK, tmpVector);
   MolDS_wrappers::Blas::GetInstance()->Dsymv(this->molecule->GetNumberAtoms()*dxy*dxy, 
                                              twoElec, 
                                              twiceMoJK,
                                              tmpVector);
   value -= MolDS_wrappers::Blas::GetInstance()->Ddot(this->molecule->GetNumberAtoms()*dxy*dxy,twiceMoIL, tmpVector);
   MallocerFreer::GetInstance()->Free<double>(&twoElec,   this->molecule->GetNumberAtoms()*dxy*dxy, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(&twiceMoIJ, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(&twiceMoIK, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(&twiceMoIL, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(&twiceMoKL, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(&twiceMoJL, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(&twiceMoJK, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(&tmpVector, this->molecule->GetNumberAtoms()*dxy*dxy);
   // End of algorithm using blas
   */

   /*
   // Second algorithm using blas.
   // This algorithm uses DGEMM.
   double** twoElec = NULL;
   double*  twiceMoIJ = NULL;
   double*  twiceMoIK = NULL;
   double*  twiceMoIL = NULL;
   double** twiceMoB  = NULL;
   double** tmpMatrix = NULL;
   int numAOs = this->molecule->GetTotalNumberAOs();
   MallocerFreer::GetInstance()->Malloc<double>(&twoElec,   this->molecule->GetNumberAtoms()*dxy*dxy, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&twiceMoIJ, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&twiceMoIK, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&twiceMoIL, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&twiceMoB, 3, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(&tmpMatrix,3, this->molecule->GetNumberAtoms()*dxy*dxy);
   for(int A=0; A<this->molecule->GetNumberAtoms(); A++){
      const Atom& atomA = *this->molecule->GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();
      for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
         for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
            twiceMoIJ[A*dxy*dxy+(mu    -firstAOIndexA)*dxy+(nu   -firstAOIndexA)]=fockMatrix[moI][mu    ]*fockMatrix[moJ][nu   ];
            twiceMoIK[A*dxy*dxy+(mu    -firstAOIndexA)*dxy+(nu   -firstAOIndexA)]=fockMatrix[moI][mu    ]*fockMatrix[moK][nu   ];
            twiceMoIL[A*dxy*dxy+(mu    -firstAOIndexA)*dxy+(nu   -firstAOIndexA)]=fockMatrix[moI][mu    ]*fockMatrix[moL][nu   ];
         }
      }
   }

   for(int B=0; B<this->molecule->GetNumberAtoms(); B++){
      const Atom& atomB = *this->molecule->GetAtom(B);
      int firstAOIndexB = atomB.GetFirstAOIndex();
      int lastAOIndexB  = atomB.GetLastAOIndex();
      for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
         for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
            twiceMoB[0][B*dxy*dxy+(lambda-firstAOIndexB)*dxy+(sigma-firstAOIndexB)]=fockMatrix[moK][lambda]*fockMatrix[moL][sigma];
            twiceMoB[1][B*dxy*dxy+(lambda-firstAOIndexB)*dxy+(sigma-firstAOIndexB)]=fockMatrix[moJ][lambda]*fockMatrix[moL][sigma];
            twiceMoB[2][B*dxy*dxy+(lambda-firstAOIndexB)*dxy+(sigma-firstAOIndexB)]=fockMatrix[moJ][lambda]*fockMatrix[moK][sigma];
         }
      }
   }

   for(int A=0; A<this->molecule->GetNumberAtoms(); A++){
      const Atom& atomA = *this->molecule->GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();
      for(int B=0; B<this->molecule->GetNumberAtoms(); B++){
         const Atom& atomB = *this->molecule->GetAtom(B);
         int firstAOIndexB = atomB.GetFirstAOIndex();
         int lastAOIndexB  = atomB.GetLastAOIndex();
         double gamma = 0.0;
         if(A!=B){
            for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
               for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
                  for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                     for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
                        twoElec[A*dxy*dxy+(mu-firstAOIndexA)*dxy+(nu-firstAOIndexA)]
                               [B*dxy*dxy+(lambda-firstAOIndexB)*dxy+(sigma-firstAOIndexB)] = 
                            this->twoElecTwoCore[A]
                                                [B]
                                                [mu-firstAOIndexA]
                                                [nu-firstAOIndexA]
                                                [lambda-firstAOIndexB]
                                                [sigma-firstAOIndexB];
                     }
                  }
               }
            }
         }
         else{
            for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
               for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
                  for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                     for(int sigma=firstAOIndexB; sigma<=lastAOIndexB; sigma++){
                        if(mu==nu && lambda==sigma){
                           OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                           OrbitalType orbitalLambda = atomB.GetValence(lambda-firstAOIndexB);
                           gamma = this->GetCoulombInt(orbitalMu, orbitalLambda, atomA);
                        }
                        else if((mu==lambda && nu==sigma) || (nu==lambda && mu==sigma) ){
                           OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                           OrbitalType orbitalNu = atomA.GetValence(nu-firstAOIndexA);
                           gamma = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
                        }
                        else{
                           gamma = 0.0;
                        }
                        twoElec[A*dxy*dxy+(mu-firstAOIndexA)*dxy+(nu-firstAOIndexA)]
                               [B*dxy*dxy+(lambda-firstAOIndexB)*dxy+(sigma-firstAOIndexB)] = gamma;
                     }  
                  }
               }
            }
         }
      }
   }

   MolDS_wrappers::Blas::GetInstance()->Dgemm(false, true, true,
                                              this->molecule->GetNumberAtoms()*dxy*dxy,
                                              3,
                                              this->molecule->GetNumberAtoms()*dxy*dxy,
                                              1.0,
                                              twoElec,
                                              twiceMoB,
                                              0.0,
                                              tmpMatrix);
   value = 4.0*MolDS_wrappers::Blas::GetInstance()->Ddot(this->molecule->GetNumberAtoms()*dxy*dxy,twiceMoIJ, &tmpMatrix[0][0]);
   value -=    MolDS_wrappers::Blas::GetInstance()->Ddot(this->molecule->GetNumberAtoms()*dxy*dxy,twiceMoIK, &tmpMatrix[1][0]);
   value -=    MolDS_wrappers::Blas::GetInstance()->Ddot(this->molecule->GetNumberAtoms()*dxy*dxy,twiceMoIL, &tmpMatrix[2][0]);
   MallocerFreer::GetInstance()->Free<double>(&twoElec,   this->molecule->GetNumberAtoms()*dxy*dxy, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(&twiceMoIJ, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(&twiceMoIK, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(&twiceMoIL, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(&twiceMoB, 3, this->molecule->GetNumberAtoms()*dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(&tmpMatrix,3, this->molecule->GetNumberAtoms()*dxy*dxy);
   // End of second algorithm using blas
   */

   /*
   // slow algorithm
   value = 4.0*this->GetMolecularIntegralElement(moI, moJ, moK, moL, 
                                                 *this->molecule, 
                                                 this->fockMatrix, NULL)
          -1.0*this->GetMolecularIntegralElement(moI, moK, moJ, moL, 
                                                 *this->molecule, 
                                                 this->fockMatrix, NULL)
          -1.0*this->GetMolecularIntegralElement(moI, moL, moJ, moK, 
                                                 *this->molecule, 
                                                 this->fockMatrix, NULL);
   */
   return value;
}

void Mndo::CalcTwoElecTwoCore(double****** twoElecTwoCore, 
                              const Molecule& molecule) const{
#ifdef MOLDS_DBG
   if(twoElecTwoCore == NULL){
      throw MolDSException(this->errorMessageCalcTwoElecTwoCoreNullMatrix);
   }
#endif
   MallocerFreer::GetInstance()->Initialize<double>(twoElecTwoCore, 
                                                    molecule.GetNumberAtoms(),
                                                    molecule.GetNumberAtoms(),
                                                    dxy, dxy, dxy, dxy);

   stringstream ompErrors;
#pragma omp parallel
   {
      double**** diatomicTwoElecTwoCore = NULL;
      try{
         MallocerFreer::GetInstance()->Malloc<double>(&diatomicTwoElecTwoCore, dxy, dxy, dxy, dxy);
         // note that terms with condition a==b are not needed to calculate. 
#pragma omp for schedule(auto)
         for(int a=0; a<molecule.GetNumberAtoms(); a++){
            for(int b=a+1; b<molecule.GetNumberAtoms(); b++){
               this->CalcDiatomicTwoElecTwoCore(diatomicTwoElecTwoCore, a, b);
               for(int mu=0; mu<dxy; mu++){
                  for(int nu=mu; nu<dxy; nu++){
                     for(int lambda=0; lambda<dxy; lambda++){
                        for(int sigma=lambda; sigma<dxy; sigma++){
                           double value = diatomicTwoElecTwoCore[mu][nu][lambda][sigma];
                           twoElecTwoCore[a][b][mu][nu][lambda][sigma] = value;
                           twoElecTwoCore[a][b][mu][nu][sigma][lambda] = value;
                           twoElecTwoCore[a][b][nu][mu][lambda][sigma] = value;
                           twoElecTwoCore[a][b][nu][mu][sigma][lambda] = value;
                           twoElecTwoCore[b][a][lambda][sigma][mu][nu] = value;
                           twoElecTwoCore[b][a][lambda][sigma][nu][mu] = value;
                           twoElecTwoCore[b][a][sigma][lambda][mu][nu] = value;
                           twoElecTwoCore[b][a][sigma][lambda][nu][mu] = value;
                        }
                     }
                  }
               }
            }
         }
      }
      catch(MolDSException ex){
#pragma omp critical
         ex.Serialize(ompErrors);
      }
      MallocerFreer::GetInstance()->Free<double>(&diatomicTwoElecTwoCore, dxy, dxy, dxy, dxy);
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
}

// Calculation of two electrons two cores integral (mu, nu | lambda, sigma) in space fixed frame, 
// taht is, Eq. (9) in ref. [DT_1977-2].
// mu and nu are included in atomA's AOs. 
// lambda and sigma are included in atomB's AOs.
// Note that atomA != atomB.
// Note taht d-orbital cannot be treated, 
// that is, matrix[dxy][dxy][dxy][dxy] cannot be treatable.
void Mndo::CalcDiatomicTwoElecTwoCore(double**** matrix, int indexAtomA, int indexAtomB) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   if(indexAtomA == indexAtomB){
      stringstream ss;
      ss << this->errorMessageCalcDiatomicTwoElecTwoCoreSameAtoms;
      ss << this->errorMessageAtomA << indexAtomA 
                                    << AtomTypeStr(atomA.GetAtomType()) << endl;
      ss << this->errorMessageAtomB << indexAtomB 
                                    << AtomTypeStr(atomB.GetAtomType()) << endl;
      throw MolDSException(ss.str());
   }

#ifdef MOLDS_DBG
   if(matrix == NULL){
      throw MolDSException(this->errorMessageCalcDiatomicTwoElecTwoCoreNullMatrix);
   }
#endif
   MallocerFreer::GetInstance()->Initialize<double>(matrix, dxy, dxy, dxy, dxy);

   // calclation in diatomic frame
   for(int mu=0; mu<atomA.GetValenceSize(); mu++){
      for(int nu=0; nu<atomA.GetValenceSize(); nu++){
         for(int lambda=0; lambda<atomB.GetValenceSize(); lambda++){
            for(int sigma=0; sigma<atomB.GetValenceSize(); sigma++){
               matrix[mu][nu][lambda][sigma] = this->GetNddoRepulsionIntegral(
                                                     atomA, 
                                                     atomA.GetValence(mu),
                                                     atomA.GetValence(nu),
                                                     atomB, 
                                                     atomB.GetValence(lambda),
                                                     atomB.GetValence(sigma));
                     
            }
         }
      }
   }

   // rotate matirix into the space frame
   double** rotatingMatrix = NULL;
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&rotatingMatrix,
                                                   OrbitalType_end, OrbitalType_end);
      this->CalcRotatingMatrix(rotatingMatrix, atomA, atomB);
      this->RotateDiatomicTwoElecTwoCoreToSpaceFrame(matrix, rotatingMatrix);
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&rotatingMatrix, OrbitalType_end, OrbitalType_end);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&rotatingMatrix, OrbitalType_end, OrbitalType_end);

   /* 
   this->OutputLog("(mu, nu | lambda, sigma) matrix\n");
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               this->OutputLog(boost::format("mu=%d nu=%d lambda=%d sigma=%d $e\n") % mu
                                                                                    % nu
                                                                                    % lambda
                                                                                    % sigma 
                                                                                    % matrix[mu][nu][lambda][sigma]);
            }
         }
      }
   }
   */
}

// Calculation of first derivatives of the two electrons two cores integral in space fixed frame,
// (mu, nu | lambda, sigma), taht is, Eq. (9) in ref. [DT_1977-2].
// mu and nu are included in atomA's AOs. 
// lambda and sigma are included in atomB's AOs.
// This derivative is related to the coordinates of atomA.
// Note that atomA != atomB.
// Note taht d-orbital cannot be treated, 
// that is, matrix[dxy][dxy][dxy][dxy][CartesianType_end] cannot be treatable.
void Mndo::CalcDiatomicTwoElecTwoCore1stDerivatives(double***** matrix, 
                                                    int indexAtomA, 
                                                    int indexAtomB) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   if(indexAtomA == indexAtomB){
      stringstream ss;
      ss << this->errorMessageCalcDiatomicTwoElecTwoCore1stDerivativesSameAtoms;
      ss << this->errorMessageAtomA << indexAtomA 
                                    << AtomTypeStr(atomA.GetAtomType()) << endl;
      ss << this->errorMessageAtomB << indexAtomB 
                                    << AtomTypeStr(atomB.GetAtomType()) << endl;
      throw MolDSException(ss.str());
   }

#ifdef MOLDS_DBG
   if(matrix == NULL){
      throw MolDSException(this->errorMessageCalcDiatomicTwoElecTwoCore1stDerivativesNullMatrix);
   }
#endif
   MallocerFreer::GetInstance()->Initialize<double>(matrix, 
                                                    dxy, 
                                                    dxy, 
                                                    dxy, 
                                                    dxy, 
                                                    CartesianType_end);

   double**   rotatingMatrix = NULL;
   double***  rotMat1stDerivatives = NULL;
   double**** diatomicTwoElecTwoCore = NULL;
   try{
      this->MallocDiatomicTwoElecTwoCore1stDeriTemps(&rotatingMatrix,
                                                     &rotMat1stDerivatives,
                                                     &diatomicTwoElecTwoCore);
      // calclation in diatomic frame
      for(int mu=0; mu<atomA.GetValenceSize(); mu++){
         for(int nu=mu; nu<atomA.GetValenceSize(); nu++){
            for(int lambda=0; lambda<atomB.GetValenceSize(); lambda++){
               for(int sigma=lambda; sigma<atomB.GetValenceSize(); sigma++){
                  for(int dimA=0; dimA<CartesianType_end; dimA++){
                     matrix[mu][nu][lambda][sigma][dimA] 
                        = this->GetNddoRepulsionIntegral1stDerivative(
                                atomA, 
                                atomA.GetValence(mu),
                                atomA.GetValence(nu),
                                atomB, 
                                atomB.GetValence(lambda),
                                atomB.GetValence(sigma),
                                static_cast<CartesianType>(dimA));
                     matrix[nu][mu][lambda][sigma][dimA] = matrix[mu][nu][lambda][sigma][dimA];
                     matrix[nu][mu][sigma][lambda][dimA] = matrix[mu][nu][lambda][sigma][dimA];
                     matrix[mu][nu][sigma][lambda][dimA] = matrix[mu][nu][lambda][sigma][dimA];
                  }  
                  diatomicTwoElecTwoCore[mu][nu][lambda][sigma] 
                     = this->GetNddoRepulsionIntegral(
                             atomA, 
                             atomA.GetValence(mu),
                             atomA.GetValence(nu),
                             atomB, 
                             atomB.GetValence(lambda),
                             atomB.GetValence(sigma));
                  diatomicTwoElecTwoCore[nu][mu][lambda][sigma] = diatomicTwoElecTwoCore[mu][nu][lambda][sigma];
                  diatomicTwoElecTwoCore[nu][mu][sigma][lambda] = diatomicTwoElecTwoCore[mu][nu][lambda][sigma];
                  diatomicTwoElecTwoCore[mu][nu][sigma][lambda] = diatomicTwoElecTwoCore[mu][nu][lambda][sigma];
               }
            }
         }
      }

      // rotate matirix into the space frame
      this->CalcRotatingMatrix(rotatingMatrix, atomA, atomB);
      this->CalcRotatingMatrix1stDerivatives(rotMat1stDerivatives, atomA, atomB);
      this->RotateDiatomicTwoElecTwoCore1stDerivativesToSpaceFrame(matrix, 
                                                                   diatomicTwoElecTwoCore,
                                                                   rotatingMatrix,
                                                                   rotMat1stDerivatives);
   }
   catch(MolDSException ex){
      this->FreeDiatomicTwoElecTwoCore1stDeriTemps(&rotatingMatrix,
                                                   &rotMat1stDerivatives,
                                                   &diatomicTwoElecTwoCore);
      throw ex;
   }
   this->FreeDiatomicTwoElecTwoCore1stDeriTemps(&rotatingMatrix,
                                                &rotMat1stDerivatives,
                                                &diatomicTwoElecTwoCore);
}

// Calculation of second derivatives of the two electrons two cores integral in space fixed frame,
// (mu, nu | lambda, sigma), taht is, Eq. (9) in ref. [DT_1977-2].
// mu and nu are included in atomA's AOs. 
// lambda and sigma are included in atomB's AOs.
// Both derivative is related to the coordinates of atomA.
// Note that atomA != atomB.
// Note taht d-orbital cannot be treated, 
// that is, matrix[dxy][dxy][dxy][dxy][CartesianType_end][CartesianType_end] cannot be treatable.
void Mndo::CalcDiatomicTwoElecTwoCore2ndDerivatives(double****** matrix, 
                                                       int indexAtomA, 
                                                       int indexAtomB) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   if(indexAtomA == indexAtomB){
      stringstream ss;
      ss << this->errorMessageCalcDiatomicTwoElecTwoCore2ndDerivativesSameAtoms;
      ss << this->errorMessageAtomA << indexAtomA 
                                    << AtomTypeStr(atomA.GetAtomType()) << endl;
      ss << this->errorMessageAtomB << indexAtomB 
                                    << AtomTypeStr(atomB.GetAtomType()) << endl;
      throw MolDSException(ss.str());
   }

#ifdef MOLDS_DBG
   if(matrix == NULL){
      throw MolDSException(this->errorMessageCalcDiatomicTwoElecTwoCore2ndDerivativesNullMatrix);
   }
#endif
   MallocerFreer::GetInstance()->Initialize<double>(matrix, 
                                                    dxy, 
                                                    dxy, 
                                                    dxy, 
                                                    dxy, 
                                                    CartesianType_end,
                                                    CartesianType_end);

   double** rotatingMatrix = NULL;
   double*** rotMat1stDerivatives = NULL;
   double**** rotMat2ndDerivatives = NULL;
   double**** diatomicTwoElecTwoCore = NULL;
   double***** diatomicTwoElecTwoCore1stDerivatives = NULL;
   try{
      this->MallocDiatomicTwoElecTwoCore2ndDeriTemps(&rotatingMatrix,
                                                     &rotMat1stDerivatives,
                                                     &rotMat2ndDerivatives,
                                                     &diatomicTwoElecTwoCore,
                                                     &diatomicTwoElecTwoCore1stDerivatives);
      // calclation in diatomic frame
      for(int mu=0; mu<atomA.GetValenceSize(); mu++){
         for(int nu=0; nu<atomA.GetValenceSize(); nu++){
            for(int lambda=0; lambda<atomB.GetValenceSize(); lambda++){
               for(int sigma=0; sigma<atomB.GetValenceSize(); sigma++){
                  for(int dimA1=0; dimA1<CartesianType_end; dimA1++){
                     for(int dimA2=0; dimA2<CartesianType_end; dimA2++){
                        matrix[mu][nu][lambda][sigma][dimA1][dimA2]
                           = this->GetNddoRepulsionIntegral2ndDerivative(
                                   atomA, 
                                   atomA.GetValence(mu),
                                   atomA.GetValence(nu),
                                   atomB, 
                                   atomB.GetValence(lambda),
                                   atomB.GetValence(sigma),
                                   static_cast<CartesianType>(dimA1),
                                   static_cast<CartesianType>(dimA2));
                     }
                     diatomicTwoElecTwoCore1stDerivatives[mu][nu][lambda][sigma][dimA1] 
                        = this->GetNddoRepulsionIntegral1stDerivative(
                                atomA, 
                                atomA.GetValence(mu),
                                atomA.GetValence(nu),
                                atomB, 
                                atomB.GetValence(lambda),
                                atomB.GetValence(sigma),
                                static_cast<CartesianType>(dimA1));
                  }  
                  diatomicTwoElecTwoCore[mu][nu][lambda][sigma] 
                     = this->GetNddoRepulsionIntegral(
                             atomA, 
                             atomA.GetValence(mu),
                             atomA.GetValence(nu),
                             atomB, 
                             atomB.GetValence(lambda),
                             atomB.GetValence(sigma));
               }
            }
         }
      }

      // rotate matirix into the space frame
      this->CalcRotatingMatrix(rotatingMatrix, atomA, atomB);
      this->CalcRotatingMatrix1stDerivatives(rotMat1stDerivatives, atomA, atomB);
      this->CalcRotatingMatrix2ndDerivatives(rotMat2ndDerivatives, atomA, atomB);
      this->RotateDiatomicTwoElecTwoCore2ndDerivativesToSpaceFrame(matrix, 
                                                                   diatomicTwoElecTwoCore,
                                                                   diatomicTwoElecTwoCore1stDerivatives,
                                                                   rotatingMatrix,
                                                                   rotMat1stDerivatives,
                                                                   rotMat2ndDerivatives);
   }
   catch(MolDSException ex){
      this->FreeDiatomicTwoElecTwoCore2ndDeriTemps(&rotatingMatrix,
                                                   &rotMat1stDerivatives,
                                                   &rotMat2ndDerivatives,
                                                   &diatomicTwoElecTwoCore,
                                                   &diatomicTwoElecTwoCore1stDerivatives);
      throw ex;
   }
   this->FreeDiatomicTwoElecTwoCore2ndDeriTemps(&rotatingMatrix,
                                                &rotMat1stDerivatives,
                                                &rotMat2ndDerivatives,
                                                &diatomicTwoElecTwoCore,
                                                &diatomicTwoElecTwoCore1stDerivatives);
}

void Mndo::MallocDiatomicTwoElecTwoCore1stDeriTemps(double*** rotatingMatrix,
                                                      double**** rotMat1stDerivatives,
                                                      double***** diatomicTwoElecTwoCore) const{
   MallocerFreer::GetInstance()->Malloc<double>(rotatingMatrix,
                                                OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Malloc<double>(rotMat1stDerivatives, 
                                                OrbitalType_end, 
                                                OrbitalType_end, 
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicTwoElecTwoCore, dxy, dxy, dxy, dxy);
}

void Mndo::MallocDiatomicTwoElecTwoCore2ndDeriTemps(double*** rotatingMatrix,
                                                    double**** rotMat1stDerivatives,
                                                    double***** rotMat2ndDerivatives,
                                                    double***** diatomicTwoElecTwoCore,
                                                    double****** diatomicTwoElecTwoCore1stDerivatives) const{
   this->MallocDiatomicTwoElecTwoCore1stDeriTemps(rotatingMatrix,
                                                  rotMat1stDerivatives,
                                                  diatomicTwoElecTwoCore);
   MallocerFreer::GetInstance()->Malloc<double>(rotMat2ndDerivatives, 
                                                OrbitalType_end, 
                                                OrbitalType_end, 
                                                CartesianType_end,
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicTwoElecTwoCore1stDerivatives, 
                                                dxy, 
                                                dxy, 
                                                dxy, 
                                                dxy,
                                                CartesianType_end);
}

void Mndo::FreeDiatomicTwoElecTwoCore1stDeriTemps(double*** rotatingMatrix,
                                                  double**** rotMat1stDerivatives,
                                                  double***** diatomicTwoElecTwoCore) const{
   MallocerFreer::GetInstance()->Free<double>(rotatingMatrix,
                                              OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Free<double>(rotMat1stDerivatives, 
                                              OrbitalType_end, 
                                              OrbitalType_end, 
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicTwoElecTwoCore, dxy, dxy, dxy, dxy);
}

void Mndo::FreeDiatomicTwoElecTwoCore2ndDeriTemps(double*** rotatingMatrix,
                                                  double**** rotMat1stDerivatives,
                                                  double***** rotMat2ndDerivatives,
                                                  double***** diatomicTwoElecTwoCore,
                                                  double****** diatomicTwoElecTwoCore1stDerivatives) const{
   this->FreeDiatomicTwoElecTwoCore1stDeriTemps(rotatingMatrix,
                                                rotMat1stDerivatives,
                                                diatomicTwoElecTwoCore);
   MallocerFreer::GetInstance()->Free<double>(rotMat2ndDerivatives, 
                                              OrbitalType_end, 
                                              OrbitalType_end, 
                                              CartesianType_end,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicTwoElecTwoCore1stDerivatives, 
                                              dxy, 
                                              dxy, 
                                              dxy, 
                                              dxy,
                                              CartesianType_end);
}

// Rotate 4-dimensional matrix from diatomic frame to space frame
// Note tha in this method d-orbitals can not be treatable.
void Mndo::RotateDiatomicTwoElecTwoCoreToSpaceFrame(double**** matrix, 
                                                    double const* const* rotatingMatrix) const{
   double oldMatrix[dxy][dxy][dxy][dxy];
   MolDS_wrappers::Blas::GetInstance()->Dcopy(dxy*dxy*dxy*dxy, &matrix[0][0][0][0], &oldMatrix[0][0][0][0]);

   // rotate (fast algorithm, see also slow algorithm shown later)
   double** twiceRotatingMatrix = NULL;
   double** ptrOldMatrix        = NULL;
   double** ptrMatrix           = NULL;
   try{
      MallocTempMatricesRotateDiatomicTwoElecTwoCore(&twiceRotatingMatrix,
                                                     &ptrOldMatrix,
                                                     &ptrMatrix);
      for(int mu=0; mu<dxy; mu++){
         for(int nu=0; nu<dxy; nu++){
            int i=mu*dxy+nu;
            for(int lambda=0; lambda<dxy; lambda++){
               for(int sigma=0; sigma<dxy; sigma++){
                  int j=lambda*dxy+sigma;
                  twiceRotatingMatrix[i][j] = rotatingMatrix[mu][lambda]*rotatingMatrix[nu][sigma];
               }
            }
            ptrOldMatrix[i] = &oldMatrix[mu][nu][0][0];
            ptrMatrix   [i] = &matrix   [mu][nu][0][0];
         }
      }
      bool isColumnMajorTwiceRotatingMatrix = false;
      bool isColumnMajorPtrOldMatrix        = false;
      double alpha = 1.0;
      double beta  = 0.0;
      MolDS_wrappers::Blas::GetInstance()->Dgemmm(isColumnMajorTwiceRotatingMatrix,
                                                  isColumnMajorPtrOldMatrix,
                                                  !isColumnMajorTwiceRotatingMatrix,
                                                  dxy*dxy, dxy*dxy, dxy*dxy, dxy*dxy,
                                                  alpha,
                                                  twiceRotatingMatrix,
                                                  ptrOldMatrix,
                                                  twiceRotatingMatrix,
                                                  beta, 
                                                  ptrMatrix);
   }
   catch(MolDSException ex){
      FreeTempMatricesRotateDiatomicTwoElecTwoCore(&twiceRotatingMatrix,
                                                   &ptrOldMatrix,
                                                   &ptrMatrix);
   }
   FreeTempMatricesRotateDiatomicTwoElecTwoCore(&twiceRotatingMatrix,
                                                &ptrOldMatrix,
                                                &ptrMatrix);
   /*
   // rotate (slow algorithm)
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               matrix[mu][nu][lambda][sigma] = 0.0;
               for(int i=0; i<dxy; i++){
                  for(int j=0; j<dxy; j++){
                     for(int k=0; k<dxy; k++){
                        for(int l=0; l<dxy; l++){
                           matrix[mu][nu][lambda][sigma] += oldMatrix[i][j][k][l] 
                                                            *rotatingMatrix[mu][i] 
                                                            *rotatingMatrix[nu][j] 
                                                            *rotatingMatrix[lambda][k] 
                                                            *rotatingMatrix[sigma][l];
                        }
                     }
                  }
               }
            }
         }
      }
   }
   */
}

void Mndo::MallocTempMatricesRotateDiatomicTwoElecTwoCore(double*** twiceRotatingMatrix,
                                                          double*** ptrOldMatrix,
                                                          double*** ptrMatrix) const{
   MallocerFreer::GetInstance()->Malloc<double>(twiceRotatingMatrix, dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double*>(ptrOldMatrix,       dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double*>(ptrMatrix,          dxy*dxy);
}

void Mndo::FreeTempMatricesRotateDiatomicTwoElecTwoCore(double*** twiceRotatingMatrix,
                                                          double*** ptrOldMatrix,
                                                          double*** ptrMatrix) const{
   MallocerFreer::GetInstance()->Free<double>(twiceRotatingMatrix, dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Free<double*>(ptrOldMatrix,       dxy*dxy);
   MallocerFreer::GetInstance()->Free<double*>(ptrMatrix,          dxy*dxy);
}

// Rotate 5-dimensional matrix from diatomic frame to space frame
// Note tha in this method d-orbitals can not be treatable.
void Mndo::RotateDiatomicTwoElecTwoCore1stDerivativesToSpaceFrame(
           double***** matrix, 
           double const* const*const* const* diatomicTwoElecTwoCore,
           double const* const* rotatingMatrix,
           double const* const* const* rotMat1stDerivatives) const{

   // rotate (fast algorithm, see also slow algorithm shown later)
   int incrementOne = 1;
   bool isColumnMajorTwiceRotatingMatrix = false;
   bool isColumnMajorOldMatrix           = false;
   double alpha;
   double beta;
   double** twiceRotatingMatrix       = NULL;
   double** twiceRotatingMatrixDerivA = NULL;
   double** twiceRotatingMatrixDerivB = NULL;
   double** oldMatrix                 = NULL;
   double** rotatedMatrix             = NULL;
   double** tmpMatrix                 = NULL;
   double** ptrDiatomic               = NULL;
   try{
      this->MallocTempMatricesRotateDiatomicTwoElecTwoCore1stDerivs(&twiceRotatingMatrix,
                                                                    &twiceRotatingMatrixDerivA,
                                                                    &twiceRotatingMatrixDerivB,
                                                                    &oldMatrix,
                                                                    &rotatedMatrix,
                                                                    &tmpMatrix,                
                                                                    &ptrDiatomic);
      for(int mu=0; mu<dxy; mu++){
         for(int nu=0; nu<dxy; nu++){
            int i=mu*dxy+nu;
            for(int lambda=0; lambda<dxy; lambda++){
               for(int sigma=0; sigma<dxy; sigma++){
                  int j=lambda*dxy+sigma;
                  twiceRotatingMatrix[i][j] = rotatingMatrix[mu][lambda]
                                             *rotatingMatrix[nu][sigma ];
               }
            }
            ptrDiatomic[i] = const_cast<double*>(&diatomicTwoElecTwoCore[mu][nu][0][0]);
         }
      }
      for(int axis=0; axis<CartesianType_end; axis++){
         for(int mu=0; mu<dxy; mu++){
            for(int nu=0; nu<dxy; nu++){
               int i=mu*dxy+nu;
               for(int lambda=0; lambda<dxy; lambda++){
                  for(int sigma=0; sigma<dxy; sigma++){
                     int j=lambda*dxy+sigma;
                     twiceRotatingMatrixDerivA[i][j] = rotMat1stDerivatives[mu][lambda][axis]
                                                      *rotatingMatrix      [nu][sigma ];
                     twiceRotatingMatrixDerivB[i][j] = rotatingMatrix      [mu][lambda] 
                                                      *rotMat1stDerivatives[nu][sigma ][axis];
                     oldMatrix[i][j] = matrix[mu][nu][lambda][sigma][axis];
                  }
               }
            }
         }
         alpha = 1.0;
         beta  = 0.0;
         MolDS_wrappers::Blas::GetInstance()->Dgemmm(isColumnMajorTwiceRotatingMatrix,
                                                     isColumnMajorOldMatrix,
                                                     !isColumnMajorTwiceRotatingMatrix,
                                                     dxy*dxy, dxy*dxy, dxy*dxy, dxy*dxy,
                                                     alpha,
                                                     twiceRotatingMatrix,
                                                     oldMatrix,
                                                     twiceRotatingMatrix,
                                                     beta, 
                                                     rotatedMatrix,
                                                     tmpMatrix);
         alpha = 1.0;
         beta  = 1.0;
         MolDS_wrappers::Blas::GetInstance()->Dgemmm(isColumnMajorTwiceRotatingMatrix,
                                                     isColumnMajorOldMatrix,
                                                     !isColumnMajorTwiceRotatingMatrix,
                                                     dxy*dxy, dxy*dxy, dxy*dxy, dxy*dxy,
                                                     alpha,
                                                     twiceRotatingMatrixDerivA,
                                                     ptrDiatomic,
                                                     twiceRotatingMatrix,
                                                     beta, 
                                                     rotatedMatrix,
                                                     tmpMatrix);
         MolDS_wrappers::Blas::GetInstance()->Dgemmm(isColumnMajorTwiceRotatingMatrix,
                                                     isColumnMajorOldMatrix,
                                                     !isColumnMajorTwiceRotatingMatrix,
                                                     dxy*dxy, dxy*dxy, dxy*dxy, dxy*dxy,
                                                     alpha,
                                                     twiceRotatingMatrixDerivB,
                                                     ptrDiatomic,
                                                     twiceRotatingMatrix,
                                                     beta, 
                                                     rotatedMatrix,
                                                     tmpMatrix);
         MolDS_wrappers::Blas::GetInstance()->Dgemmm(isColumnMajorTwiceRotatingMatrix,
                                                     isColumnMajorOldMatrix,
                                                     !isColumnMajorTwiceRotatingMatrix,
                                                     dxy*dxy, dxy*dxy, dxy*dxy, dxy*dxy,
                                                     alpha,
                                                     twiceRotatingMatrix,
                                                     ptrDiatomic,
                                                     twiceRotatingMatrixDerivA,
                                                     beta, 
                                                     rotatedMatrix,
                                                     tmpMatrix);
         MolDS_wrappers::Blas::GetInstance()->Dgemmm(isColumnMajorTwiceRotatingMatrix,
                                                     isColumnMajorOldMatrix,
                                                     !isColumnMajorTwiceRotatingMatrix,
                                                     dxy*dxy, dxy*dxy, dxy*dxy, dxy*dxy,
                                                     alpha,
                                                     twiceRotatingMatrix,
                                                     ptrDiatomic,
                                                     twiceRotatingMatrixDerivB,
                                                     beta, 
                                                     rotatedMatrix,
                                                     tmpMatrix);

         MolDS_wrappers::Blas::GetInstance()->Dcopy(dxy*dxy*dxy*dxy, 
                                                    &rotatedMatrix[0][0]     , incrementOne,
                                                    &matrix[0][0][0][0][axis], CartesianType_end);
      }
   }
   catch(MolDSException ex){
      this->FreeTempMatricesRotateDiatomicTwoElecTwoCore1stDerivs(&twiceRotatingMatrix,
                                                                  &twiceRotatingMatrixDerivA,
                                                                  &twiceRotatingMatrixDerivB,
                                                                  &oldMatrix,
                                                                  &rotatedMatrix,
                                                                  &tmpMatrix,                
                                                                  &ptrDiatomic);
      throw ex;
   }
   this->FreeTempMatricesRotateDiatomicTwoElecTwoCore1stDerivs(&twiceRotatingMatrix,
                                                               &twiceRotatingMatrixDerivA,
                                                               &twiceRotatingMatrixDerivB,
                                                               &oldMatrix,
                                                               &rotatedMatrix,
                                                               &tmpMatrix,                
                                                               &ptrDiatomic);

   /*
   // rotate (slow algorithm)
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               for(int c=0; c<CartesianType_end; c++){
                  matrix[mu][nu][lambda][sigma][c] = 0.0;
                  for(int i=0; i<dxy; i++){
                     for(int j=0; j<dxy; j++){
                        for(int k=0; k<dxy; k++){
                           for(int l=0; l<dxy; l++){
                              matrix[mu][nu][lambda][sigma][c] 
                                 += oldMatrix[i][j][k][l][c]
                                   *rotatingMatrix[mu][i] 
                                   *rotatingMatrix[nu][j] 
                                   *rotatingMatrix[lambda][k] 
                                   *rotatingMatrix[sigma][l];
                              matrix[mu][nu][lambda][sigma][c] 
                                 += diatomicTwoElecTwoCore[i][j][k][l]
                                   *rotMat1stDerivatives[mu][i][c]
                                   *rotatingMatrix[nu][j] 
                                   *rotatingMatrix[lambda][k] 
                                   *rotatingMatrix[sigma][l];
                              matrix[mu][nu][lambda][sigma][c] 
                                 += diatomicTwoElecTwoCore[i][j][k][l]
                                   *rotatingMatrix[mu][i] 
                                   *rotMat1stDerivatives[nu][j][c]
                                   *rotatingMatrix[lambda][k] 
                                   *rotatingMatrix[sigma][l];
                              matrix[mu][nu][lambda][sigma][c] 
                                 += diatomicTwoElecTwoCore[i][j][k][l]
                                   *rotatingMatrix[mu][i] 
                                   *rotatingMatrix[nu][j] 
                                   *rotMat1stDerivatives[lambda][k][c]
                                   *rotatingMatrix[sigma][l];
                              matrix[mu][nu][lambda][sigma][c] 
                                 += diatomicTwoElecTwoCore[i][j][k][l]
                                   *rotatingMatrix[mu][i] 
                                   *rotatingMatrix[nu][j] 
                                   *rotatingMatrix[lambda][k] 
                                   *rotMat1stDerivatives[sigma][l][c];
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   */
}

void Mndo::MallocTempMatricesRotateDiatomicTwoElecTwoCore1stDerivs(double*** twiceRotatingMatrix,
                                                                   double*** twiceRotatingMatrixDerivA,
                                                                   double*** twiceRotatingMatrixDerivB,
                                                                   double*** oldMatrix,
                                                                   double*** rotatedMatrix,
                                                                   double*** tmpMatrix,                
                                                                   double*** ptrDiatomic) const{
   MallocerFreer::GetInstance()->Malloc<double>(twiceRotatingMatrix,       dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(twiceRotatingMatrixDerivA, dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(twiceRotatingMatrixDerivB, dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(oldMatrix,                 dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(rotatedMatrix,             dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double>(tmpMatrix,                 dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Malloc<double*>(ptrDiatomic,              dxy*dxy);
}

void Mndo::FreeTempMatricesRotateDiatomicTwoElecTwoCore1stDerivs(double*** twiceRotatingMatrix,
                                                                 double*** twiceRotatingMatrixDerivA,
                                                                 double*** twiceRotatingMatrixDerivB,
                                                                 double*** oldMatrix,
                                                                 double*** rotatedMatrix,
                                                                 double*** tmpMatrix,                
                                                                 double*** ptrDiatomic) const{
   MallocerFreer::GetInstance()->Free<double>(twiceRotatingMatrix,       dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(twiceRotatingMatrixDerivA, dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(twiceRotatingMatrixDerivB, dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(oldMatrix,                 dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(rotatedMatrix,                 dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Free<double>(tmpMatrix,                 dxy*dxy, dxy*dxy);
   MallocerFreer::GetInstance()->Free<double*>(ptrDiatomic,              dxy*dxy);
}

// Rotate 6-dimensional matrix from diatomic frame to space frame
// Note tha in this method d-orbitals can not be treatable.
void Mndo::RotateDiatomicTwoElecTwoCore2ndDerivativesToSpaceFrame(
           double****** matrix, 
           double const* const* const* const* diatomicTwoElecTwoCore,
           double const* const* const* const* const* diatomicTwoElecTwoCore1stDerivatives,
           double const* const* rotatingMatrix,
           double const* const* const* rotMat1stDerivatives,
           double const* const* const* const* rotMat2ndDerivatives) const{
   double oldMatrix[dxy][dxy][dxy][dxy][CartesianType_end][CartesianType_end];
   for(int mu=s; mu<dxy; mu++){
      for(int nu=s; nu<dxy; nu++){
         for(int lambda=s; lambda<dxy; lambda++){
            for(int sigma=s; sigma<dxy; sigma++){
               for(int dimA1=XAxis; dimA1<CartesianType_end; dimA1++){
                  for(int dimA2=XAxis; dimA2<CartesianType_end; dimA2++){
                     oldMatrix[mu][nu][lambda][sigma][dimA1][dimA2] = matrix[mu][nu][lambda][sigma][dimA1][dimA2];
                  }
               }
            }
         }
      }
   }

   // rotate (fast algorithm, see also slow algorithm shown later)
   int numberTerms = 25;
   double* tempIJK = NULL;
   double* tempIJ = NULL;
   double* tempI = NULL;
   MallocerFreer::GetInstance()->Malloc<double>(&tempIJK, numberTerms);
   MallocerFreer::GetInstance()->Malloc<double>(&tempIJ, numberTerms);
   MallocerFreer::GetInstance()->Malloc<double>(&tempI, numberTerms);
   try{ 
      for(int mu=s; mu<dxy; mu++){
         for(int nu=mu; nu<dxy; nu++){
            for(int lambda=s; lambda<dxy; lambda++){
               for(int sigma=lambda; sigma<dxy; sigma++){
                  for(int dimA1=XAxis; dimA1<CartesianType_end; dimA1++){
                     for(int dimA2=dimA1; dimA2<CartesianType_end; dimA2++){

                        double value=0.0;
                        for(int i=s; i<dxy; i++){
                           MallocerFreer::GetInstance()->Initialize<double>(tempI, numberTerms);
                           for(int j=s; j<dxy; j++){
                              MallocerFreer::GetInstance()->Initialize<double>(tempIJ, numberTerms);
                              for(int k=s; k<dxy; k++){
                                 MallocerFreer::GetInstance()->Initialize<double>(tempIJK, numberTerms);
                                 for(int l=s; l<dxy; l++){
                                    
                                    tempIJK[0]  += oldMatrix             [i][j][k][l][dimA1][dimA2]*rotatingMatrix      [sigma][l];
                                    tempIJK[1]  += diatomicTwoElecTwoCore[i][j][k][l]              *rotatingMatrix      [sigma][l];
                                    tempIJK[4]  += diatomicTwoElecTwoCore[i][j][k][l]              *rotMat2ndDerivatives[sigma][l][dimA1][dimA2];

                                    tempIJK[5]  += diatomicTwoElecTwoCore1stDerivatives[i][j][k][l][dimA1]*rotatingMatrix      [sigma][l];
                                    tempIJK[8]  += diatomicTwoElecTwoCore1stDerivatives[i][j][k][l][dimA1]*rotMat1stDerivatives[sigma][l][dimA2];

                                    tempIJK[9]  += diatomicTwoElecTwoCore1stDerivatives[i][j][k][l][dimA2]*rotatingMatrix      [sigma][l];
                                    tempIJK[12] += diatomicTwoElecTwoCore              [i][j][k][l]       *rotMat1stDerivatives[sigma][l][dimA2];

                                    tempIJK[21] += diatomicTwoElecTwoCore1stDerivatives[i][j][k][l][dimA2]*rotMat1stDerivatives[sigma][l][dimA1];
                                    tempIJK[22] += diatomicTwoElecTwoCore              [i][j][k][l]       *rotMat1stDerivatives[sigma][l][dimA1];
                                 }
                                 tempIJ[0]  += tempIJK[0] *rotatingMatrix      [lambda][k];
                                 tempIJ[1]  += tempIJK[1] *rotatingMatrix      [lambda][k];
                                 tempIJ[3]  += tempIJK[1] *rotMat2ndDerivatives[lambda][k][dimA1][dimA2];
                                 tempIJ[4]  += tempIJK[4] *rotatingMatrix      [lambda][k];
                                           
                                 tempIJ[5]  += tempIJK[5] *rotatingMatrix      [lambda][k];
                                 tempIJ[7]  += tempIJK[5] *rotMat1stDerivatives[lambda][k][dimA2];
                                 tempIJ[8]  += tempIJK[8] *rotatingMatrix      [lambda][k];

                                 tempIJ[9]  += tempIJK[9] *rotatingMatrix      [lambda][k];
                                 tempIJ[11] += tempIJK[1] *rotMat1stDerivatives[lambda][k][dimA2];
                                 tempIJ[12] += tempIJK[12]*rotatingMatrix      [lambda][k];

                                 tempIJ[17] += tempIJK[9] *rotMat1stDerivatives[lambda][k][dimA1];
                                 tempIJ[18] += tempIJK[1 ]*rotMat1stDerivatives[lambda][k][dimA1];
                                 tempIJ[20] += tempIJK[12]*rotMat1stDerivatives[lambda][k][dimA1];

                                 tempIJ[21] += tempIJK[21]*rotatingMatrix      [lambda][k];
                                 tempIJ[22] += tempIJK[22]*rotatingMatrix      [lambda][k];
                                 tempIJ[24] += tempIJK[22]*rotMat1stDerivatives[lambda][k][dimA2];
                              }
                              tempI[0]  += tempIJ[0] *rotatingMatrix      [nu][j];
                              tempI[1]  += tempIJ[1] *rotatingMatrix      [nu][j];
                              tempI[2]  += tempIJ[1] *rotMat2ndDerivatives[nu][j][dimA1][dimA2];
                              tempI[3]  += tempIJ[3] *rotatingMatrix      [nu][j];
                              tempI[4]  += tempIJ[4] *rotatingMatrix      [nu][j];
                                       
                              tempI[5]  += tempIJ[5] *rotatingMatrix      [nu][j];
                              tempI[6]  += tempIJ[5] *rotMat1stDerivatives[nu][j][dimA2];
                              tempI[7]  += tempIJ[7] *rotatingMatrix      [nu][j];
                              tempI[8]  += tempIJ[8] *rotatingMatrix      [nu][j];
                                       
                              tempI[9]  += tempIJ[9] *rotatingMatrix      [nu][j];
                              tempI[10] += tempIJ[1] *rotMat1stDerivatives[nu][j][dimA2];
                              tempI[11] += tempIJ[11]*rotatingMatrix      [nu][j];
                              tempI[12] += tempIJ[12]*rotatingMatrix      [nu][j];

                              tempI[13] += tempIJ[9] *rotMat1stDerivatives[nu][j][dimA1];
                              tempI[14] += tempIJ[1] *rotMat1stDerivatives[nu][j][dimA1];
                              tempI[15] += tempIJ[11]*rotMat1stDerivatives[nu][j][dimA1];
                              tempI[16] += tempIJ[12]*rotMat1stDerivatives[nu][j][dimA1];

                              tempI[17] += tempIJ[17]*rotatingMatrix      [nu][j];
                              tempI[18] += tempIJ[18]*rotatingMatrix      [nu][j];
                              tempI[19] += tempIJ[18]*rotMat1stDerivatives[nu][j][dimA2];
                              tempI[20] += tempIJ[20]*rotatingMatrix      [nu][j];

                              tempI[21] += tempIJ[21]*rotatingMatrix      [nu][j];
                              tempI[22] += tempIJ[22]*rotatingMatrix      [nu][j];
                              tempI[23] += tempIJ[22]*rotMat1stDerivatives[nu][j][dimA2];
                              tempI[24] += tempIJ[24]*rotatingMatrix      [nu][j];
                           }
                           value += (tempI[0] + tempI[2] + tempI[3] + tempI[4] + 
                                     tempI[6] + tempI[7] + tempI[8] + tempI[13] + 
                                     tempI[15] + tempI[16] + tempI[17] + tempI[19] + 
                                     tempI[20] + tempI[21] + tempI[23] + tempI[24])*rotatingMatrix[mu][i];
                           value += (tempI[9] + tempI[10] + tempI[11] + tempI[12]) *rotMat1stDerivatives[mu][i][dimA1];
                           value += (tempI[5] + tempI[14] + tempI[18] + tempI[22]) *rotMat1stDerivatives[mu][i][dimA2];
                           value += tempI[1] *rotMat2ndDerivatives[mu][i][dimA1][dimA2];
                        }
                        matrix[mu][nu][lambda][sigma][dimA1][dimA2] = value;
                        matrix[mu][nu][sigma][lambda][dimA1][dimA2] = value;
                        matrix[nu][mu][lambda][sigma][dimA1][dimA2] = value;
                        matrix[nu][mu][sigma][lambda][dimA1][dimA2] = value;

                        matrix[mu][nu][lambda][sigma][dimA2][dimA1] = value;
                        matrix[mu][nu][sigma][lambda][dimA2][dimA1] = value;
                        matrix[nu][mu][lambda][sigma][dimA2][dimA1] = value;
                        matrix[nu][mu][sigma][lambda][dimA2][dimA1] = value;
                     }
                  }
               }
            }
         }
      }
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&tempIJK, numberTerms);
      MallocerFreer::GetInstance()->Free<double>(&tempIJ, numberTerms);
      MallocerFreer::GetInstance()->Free<double>(&tempI, numberTerms);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&tempIJK, numberTerms);
   MallocerFreer::GetInstance()->Free<double>(&tempIJ, numberTerms);
   MallocerFreer::GetInstance()->Free<double>(&tempI, numberTerms);

   /*
   // rotate (slow algorithm shown later)
   for(int mu=s; mu<dxy; mu++){
      for(int nu=s; nu<dxy; nu++){
         for(int lambda=s; lambda<dxy; lambda++){
            for(int sigma=s; sigma<dxy; sigma++){
               for(int dimA1=XAxis; dimA1<CartesianType_end; dimA1++){
                  for(int dimA2=XAxis; dimA2<CartesianType_end; dimA2++){


                     matrix[mu][nu][lambda][sigma][dimA1][dimA2] = 0.0;
                     double value=0.0;
                     for(int i=s; i<dxy; i++){
                        for(int j=s; j<dxy; j++){
                           for(int k=s; k<dxy; k++){
                              for(int l=s; l<dxy; l++){
                                 // term0
                                 value += oldMatrix[i][j][k][l][dimA1][dimA2]
                                         *rotatingMatrix         [mu    ][i] 
                                         *rotatingMatrix         [nu    ][j] 
                                         *rotatingMatrix         [lambda][k] 
                                         *rotatingMatrix         [sigma ][l];
                                 // term1
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMat2ndDerivatives[mu    ][i][dimA1][dimA2]
                                         *rotatingMatrix         [nu    ][j] 
                                         *rotatingMatrix         [lambda][k] 
                                         *rotatingMatrix         [sigma ][l];
                                 // term2
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix         [mu    ][i] 
                                         *rotMat2ndDerivatives[nu    ][j][dimA1][dimA2]
                                         *rotatingMatrix         [lambda][k] 
                                         *rotatingMatrix         [sigma ][l];
                                 // term3
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix         [mu    ][i] 
                                         *rotatingMatrix         [nu    ][j] 
                                         *rotMat2ndDerivatives[lambda][k][dimA1][dimA2]
                                         *rotatingMatrix         [sigma ][l];
                                 // term4
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix         [mu    ][i] 
                                         *rotatingMatrix         [nu    ][j] 
                                         *rotatingMatrix         [lambda][k] 
                                         *rotMat2ndDerivatives[sigma ][l][dimA1][dimA2];
                                 
                                 // term5
                                 value += diatomicTwoElecTwoCore1stDerivatives[i][j][k][l][dimA1]
                                         *rotMat1stDerivatives[mu    ][i][dimA2]
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotatingMatrix        [lambda][k] 
                                         *rotatingMatrix        [sigma ][l];
                                 // term6
                                 value += diatomicTwoElecTwoCore1stDerivatives[i][j][k][l][dimA1]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotMat1stDerivatives[nu    ][j][dimA2]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotatingMatrix        [sigma ][l];
                                 // term7
                                 value += diatomicTwoElecTwoCore1stDerivatives[i][j][k][l][dimA1]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotMat1stDerivatives[lambda][k][dimA2]
                                         *rotatingMatrix        [sigma ][l];
                                 // term8
                                 value += diatomicTwoElecTwoCore1stDerivatives[i][j][k][l][dimA1]
                                         *rotatingMatrix        [mu    ][i]
                                         *rotatingMatrix        [nu    ][j]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotMat1stDerivatives[sigma ][l][dimA2];

                                 // term9
                                 value += diatomicTwoElecTwoCore1stDerivatives[i][j][k][l][dimA2]
                                         *rotMat1stDerivatives[mu    ][i][dimA1]
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotatingMatrix        [lambda][k] 
                                         *rotatingMatrix        [sigma ][l];
                                 // term10
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMat1stDerivatives[mu    ][i][dimA1]
                                         *rotMat1stDerivatives[nu    ][j][dimA2]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotatingMatrix        [sigma ][l];
                                 // term11
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMat1stDerivatives[mu    ][i][dimA1]
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotMat1stDerivatives[lambda][k][dimA2]
                                         *rotatingMatrix        [sigma ][l];
                                 // term12
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMat1stDerivatives[mu    ][i][dimA1]
                                         *rotatingMatrix        [nu    ][j]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotMat1stDerivatives[sigma ][l][dimA2];

                                 // term13
                                 value += diatomicTwoElecTwoCore1stDerivatives[i][j][k][l][dimA2]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotMat1stDerivatives[nu    ][j][dimA1]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotatingMatrix        [sigma ][l];
                                 // term14
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMat1stDerivatives[mu    ][i][dimA2]
                                         *rotMat1stDerivatives[nu    ][j][dimA1]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotatingMatrix        [sigma ][l];
                                 // term15
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotMat1stDerivatives[nu    ][j][dimA1]
                                         *rotMat1stDerivatives[lambda][k][dimA2]
                                         *rotatingMatrix        [sigma ][l];
                                 // term16
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotMat1stDerivatives[nu    ][j][dimA1]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotMat1stDerivatives[sigma ][l][dimA2];

                                 // term17
                                 value += diatomicTwoElecTwoCore1stDerivatives[i][j][k][l][dimA2]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotMat1stDerivatives[lambda][k][dimA1]
                                         *rotatingMatrix        [sigma ][l];
                                 // term18
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMat1stDerivatives[mu    ][i][dimA2]
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotMat1stDerivatives[lambda][k][dimA1]
                                         *rotatingMatrix        [sigma ][l];
                                 // term19
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotMat1stDerivatives[nu    ][j][dimA2]
                                         *rotMat1stDerivatives[lambda][k][dimA1]
                                         *rotatingMatrix        [sigma ][l];
                                 // term20
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotMat1stDerivatives[lambda][k][dimA1]
                                         *rotMat1stDerivatives[sigma ][l][dimA2];
                                 
                                 // term21
                                 value += diatomicTwoElecTwoCore1stDerivatives[i][j][k][l][dimA2]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotatingMatrix        [lambda][k] 
                                         *rotMat1stDerivatives[sigma ][l][dimA1];
                                 // term22
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMat1stDerivatives[mu    ][i][dimA2]
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotatingMatrix        [lambda][k] 
                                         *rotMat1stDerivatives[sigma ][l][dimA1];
                                 // term23
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotMat1stDerivatives[nu    ][j][dimA2]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotMat1stDerivatives[sigma ][l][dimA1];
                                 // term24
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotMat1stDerivatives[lambda][k][dimA2]
                                         *rotMat1stDerivatives[sigma ][l][dimA1];
                              }
                           }
                        }
                     }
                     matrix[mu][nu][lambda][sigma][dimA1][dimA2] = value;

                  }
               }
            }
         }
      }
   }
   */ 
}

// See Apendix in [DT_1977]
// Orbital mu and nu belong atom A, 
// orbital lambda and sigma belong atomB.
double Mndo::GetNddoRepulsionIntegral(const Atom& atomA, 
                                      OrbitalType mu, 
                                      OrbitalType nu,
                                      const Atom& atomB, 
                                      OrbitalType lambda, 
                                      OrbitalType sigma) const{
   double value = 0.0;
   double DA=0.0;
   double DB=0.0;
   double rhoA = 0.0;
   double rhoB = 0.0;
   double rAB = this->molecule->GetDistanceAtoms(atomA, atomB);
   int lA = 0;
   int lB = 0;
   // (28) in [DT_1977]
   if(mu == s && nu == s && lambda == s && sigma == s){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
   }
   // (29) in [DT_1977]
   else if(mu == s && nu == s && lambda == px && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, Qxx, rAB);
      value = temp1 + temp2;
   }
   else if(mu == s && nu == s && lambda == py && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, Qyy, rAB);
      value = temp1 + temp2;
   }
   // (30) in [DT_1977]
   else if(mu == s && nu == s && lambda == pz && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, Qzz, rAB);
      value = temp1 + temp2;
   }
   // (31) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == s){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxx, sQ, rAB);
      value = temp1 + temp2;
   }
   else if(mu == py && nu == py && lambda == s && sigma == s){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qyy, sQ, rAB);
      value = temp1 + temp2;
   }
   // (32) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == s){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qzz, sQ, rAB);
      value = temp1 + temp2;
   }
   // (33) in [DT_1977]
   else if(mu == px && nu == px && lambda == px && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, Qxx, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxx, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxx, Qxx, rAB);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == py && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, Qyy, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qyy, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qyy, Qyy, rAB);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (34) in [DT_1977]
   else if(mu == px && nu == px && lambda == py && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, Qyy, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxx, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxx, Qyy, rAB);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == px && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, Qxx, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qyy, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qyy, Qxx, rAB);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (35) in [DT_1977]
   else if(mu == px && nu == px && lambda == pz && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, Qzz, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxx, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxx, Qzz, rAB);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, Qzz, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qyy, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qyy, Qzz, rAB);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (36) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == px && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, Qxx, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qzz, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qzz, Qxx, rAB);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == pz && nu == pz && lambda == py && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, Qyy, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qzz, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qzz, Qyy, rAB);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (37) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == pz && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, Qzz, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qzz, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qzz, Qzz, rAB);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (38) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == s){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, muz, sQ, rAB);
      value = temp1;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   // (39) in [DT_1977]
   else if(mu == s && nu == pz && lambda == px && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, muz, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, muz, Qxx, rAB);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == px && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == pz && lambda == py && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, muz, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, muz, Qyy, rAB);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == py && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   // (40) in [DT_1977]
   else if(mu == s && nu == pz && lambda == pz && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, muz, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, muz, Qzz, rAB);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   // (41) in [DT_1977]
   else if(mu == s && nu == s && lambda == s && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, muz, rAB);
      value = temp1;
   }
   else if(mu == s && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   // (42) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, muz, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxx, muz, rAB);
      value = temp1 + temp2;
   }
   else if(mu == px && nu == px && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == py && lambda == s && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, muz, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qyy, muz, rAB);
      value = temp1 + temp2;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   // (43) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, sQ, muz, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qzz, muz, rAB);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   // (44) in [DT_1977]
   else if(mu == s && nu == px && lambda == s && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, mux, mux, rAB);
      value = temp1;
   }
   else if(mu == px && nu == s && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == px && nu == s && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   else if(mu == s && nu == py && lambda == s && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, muy, muy, rAB);
      value = temp1;
   }
   else if(mu == py && nu == s && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == s && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (45) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, muz, muz, rAB);
      value = temp1;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (46) in [DT_1977]
   else if(mu == s && nu == px && lambda == px && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, mux, Qxz, rAB);
      value = temp1;
   }
   else if(mu == px && nu == s && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == px && nu == s && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   else if(mu == s && nu == py && lambda == py && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, muy, Qyz, rAB);
      value = temp1;
   }
   else if(mu == py && nu == s && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == s && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (47) in [DT_1977]
   else if(mu == px && nu == pz && lambda == s && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxz, mux, rAB);
      value = temp1;
   }
   else if(mu == pz && nu == px && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == px && nu == pz && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == pz && lambda == s && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qyz, muy, rAB);
      value = temp1;
   }
   else if(mu == pz && nu == py && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (48) in [DT_1977]
   else if(mu == px && nu == pz && lambda == px && sigma == pz){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxz, Qxz, rAB);
   }
   else if(mu == pz && nu == px && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == px && nu == pz && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == pz){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qyz, Qyz, rAB);
   }
   else if(mu == pz && nu == py && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == py && nu == pz && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == pz && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // (49) in [DT_1977] and p19 in [MOPAC_1990]
   else if(mu == px && nu == py && lambda == px && sigma == py){
      value = 0.5*(this->GetNddoRepulsionIntegral(atomA, mu, mu, atomB, mu, mu)
                  -this->GetNddoRepulsionIntegral(atomA, mu, mu, atomB, nu, nu));
   }
   else if(mu == py && nu == px && lambda == px && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == px && nu == py && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == px && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, sigma, lambda);
   }
   // d-orbitals
   else if(mu == dxy || mu == dyz || mu == dzz || mu == dzx || mu == dxxyy ||
           nu == dxy || nu == dyz || nu == dzz || nu == dzx || nu == dxxyy ||
           lambda == dxy || lambda == dyz || lambda == dzz || lambda  == dzx || lambda == dxxyy ||
           sigma == dxy || sigma == dyz || sigma == dzz || sigma  == dzx || sigma == dxxyy){

      stringstream ss;
      ss << this->errorMessageGetNddoRepulsionIntegral;
      ss << this->errorMessageAtomA << AtomTypeStr(atomA.GetAtomType()) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(mu) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(nu) << endl;
      ss << this->errorMessageAtomB << AtomTypeStr(atomB.GetAtomType()) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(lambda) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(sigma) << endl;
      throw MolDSException(ss.str());
   }
   else{
      value = 0.0;
   }
   return value;
}

// First derivative of NDDO repulsion integral.
// This derivation is related to the coordinate of atomA
// See Apendix in [DT_1977]
// Orbital mu and nu belong atom A, 
// orbital lambda and sigma belong atomB.
double Mndo::GetNddoRepulsionIntegral1stDerivative(
                                       const Atom& atomA, OrbitalType mu, OrbitalType nu,
                                       const Atom& atomB, OrbitalType lambda, OrbitalType sigma,
                                       CartesianType axisA) const{
   double value = 0.0;
   double DA=0.0;
   double DB=0.0;
   double rhoA = 0.0;
   double rhoB = 0.0;
   double rAB = this->molecule->GetDistanceAtoms(atomA, atomB);
   double drABDa = (atomA.GetXyz()[axisA] - atomB.GetXyz()[axisA])/rAB;
   int lA = 0;
   int lB = 0;
   // (28) in [DT_1977]
   if(mu == s && nu == s && lambda == s && sigma == s){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      value *= drABDa;
   }
   // (29) in [DT_1977]
   else if(mu == s && nu == s && lambda == px && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qxx, rAB);
      value = temp1 + temp2;
      value *= drABDa;
   }
   else if(mu == s && nu == s && lambda == py && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qyy, rAB);
      value = temp1 + temp2;
      value *= drABDa;
   }
   // (30) in [DT_1977]
   else if(mu == s && nu == s && lambda == pz && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qzz, rAB);
      value = temp1 + temp2;
      value *= drABDa;
   }
   // (31) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == s){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, sQ, rAB);
      value = temp1 + temp2;
      value *= drABDa;
   }
   else if(mu == py && nu == py && lambda == s && sigma == s){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, sQ, rAB);
      value = temp1 + temp2;
      value *= drABDa;
   }
   // (32) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == s){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, sQ, rAB);
      value = temp1 + temp2;
      value *= drABDa;
   }
   // (33) in [DT_1977]
   else if(mu == px && nu == px && lambda == px && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qxx, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, Qxx, rAB);
      value = temp1 + temp2 + temp3 + temp4;
      value *= drABDa;
   }
   else if(mu == py && nu == py && lambda == py && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qyy, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, Qyy, rAB);
      value = temp1 + temp2 + temp3 + temp4;
      value *= drABDa;
   }
   // (34) in [DT_1977]
   else if(mu == px && nu == px && lambda == py && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qyy, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, Qyy, rAB);
      value = temp1 + temp2 + temp3 + temp4;
      value *= drABDa;
   }
   else if(mu == py && nu == py && lambda == px && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qxx, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, Qxx, rAB);
      value = temp1 + temp2 + temp3 + temp4;
      value *= drABDa;
   }
   // (35) in [DT_1977]
   else if(mu == px && nu == px && lambda == pz && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qzz, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, Qzz, rAB);
      value = temp1 + temp2 + temp3 + temp4;
      value *= drABDa;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qzz, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, Qzz, rAB);
      value = temp1 + temp2 + temp3 + temp4;
      value *= drABDa;
   }
   // (36) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == px && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qxx, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, Qxx, rAB);
      value = temp1 + temp2 + temp3 + temp4;
      value *= drABDa;
   }
   else if(mu == pz && nu == pz && lambda == py && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qyy, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, Qyy, rAB);
      value = temp1 + temp2 + temp3 + temp4;
      value *= drABDa;
   }
   // (37) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == pz && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qzz, rAB);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, sQ, rAB);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, Qzz, rAB);
      value = temp1 + temp2 + temp3 + temp4;
      value *= drABDa;
   }
   // (38) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == s){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, sQ, rAB);
      value = temp1;
      value *= drABDa;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   // (39) in [DT_1977]
   else if(mu == s && nu == pz && lambda == px && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, Qxx, rAB);
      value = temp1 + temp2;
      value *= drABDa;
   }
   else if(mu == pz && nu == s && lambda == px && sigma == px){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == pz && lambda == py && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, Qyy, rAB);
      value = temp1 + temp2;
      value *= drABDa;
   }
   else if(mu == pz && nu == s && lambda == py && sigma == py){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   // (40) in [DT_1977]
   else if(mu == s && nu == pz && lambda == pz && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, sQ, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, Qzz, rAB);
      value = temp1 + temp2;
      value *= drABDa;
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == pz){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   // (41) in [DT_1977]
   else if(mu == s && nu == s && lambda == s && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, muz, rAB);
      value = temp1;
      value *= drABDa;
   }
   else if(mu == s && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   // (42) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, muz, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, muz, rAB);
      value = temp1 + temp2;
      value *= drABDa;
   }
   else if(mu == px && nu == px && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == py && lambda == s && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, muz, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, muz, rAB);
      value = temp1 + temp2;
      value *= drABDa;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   // (43) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, muz, rAB);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, muz, rAB);
      value = temp1 + temp2;
      value *= drABDa;
   }
   else if(mu == pz && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   // (44) in [DT_1977]
   else if(mu == s && nu == px && lambda == s && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, mux, mux, rAB);
      value = temp1;
      value *= drABDa;
   }
   else if(mu == px && nu == s && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == px && nu == s && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   else if(mu == s && nu == py && lambda == s && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muy, muy, rAB);
      value = temp1;
      value *= drABDa;
   }
   else if(mu == py && nu == s && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == s && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (45) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, muz, rAB);
      value = temp1;
      value *= drABDa;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == pz){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (46) in [DT_1977]
   else if(mu == s && nu == px && lambda == px && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, mux, Qxz, rAB);
      value = temp1;
      value *= drABDa;
   }
   else if(mu == px && nu == s && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == px && nu == s && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   else if(mu == s && nu == py && lambda == py && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muy, Qyz, rAB);
      value = temp1;
      value *= drABDa;
   }
   else if(mu == py && nu == s && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == s && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (47) in [DT_1977]
   else if(mu == px && nu == pz && lambda == s && sigma == px){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxz, mux, rAB);
      value = temp1;
      value *= drABDa;
   }
   else if(mu == pz && nu == px && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == px && nu == pz && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == pz && lambda == s && sigma == py){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyz, muy, rAB);
      value = temp1;
      value *= drABDa;
   }
   else if(mu == pz && nu == py && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (48) in [DT_1977]
   else if(mu == px && nu == pz && lambda == px && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxz, Qxz, rAB);
      value = temp1;
      value *= drABDa;
   }
   else if(mu == pz && nu == px && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == px && nu == pz && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == pz){
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyz, Qyz, rAB);
      value = temp1;
      value *= drABDa;
   }
   else if(mu == pz && nu == py && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == py && nu == pz && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (49) in [DT_1977] and p19 in [MOPAC_1990]
   else if(mu == px && nu == py && lambda == px && sigma == py){
      value = 0.5*(this->GetNddoRepulsionIntegral1stDerivative(
                         atomA, mu, mu, atomB, mu, mu, axisA)
                  -this->GetNddoRepulsionIntegral1stDerivative(
                         atomA, mu, mu, atomB, nu, nu, axisA));
   }
   else if(mu == py && nu == px && lambda == px && sigma == py){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == px && nu == py && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == px && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegral1stDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // d-orbitals
   else if(mu == dxy || mu == dyz || mu == dzz || mu == dzx || mu == dxxyy ||
           nu == dxy || nu == dyz || nu == dzz || nu == dzx || nu == dxxyy ||
           lambda == dxy || lambda == dyz || lambda == dzz || lambda  == dzx || lambda == dxxyy ||
           sigma == dxy || sigma == dyz || sigma == dzz || sigma  == dzx || sigma == dxxyy){

      stringstream ss;
      ss << this->errorMessageGetNddoRepulsionIntegral1stDerivative;
      ss << this->errorMessageAtomA << AtomTypeStr(atomA.GetAtomType()) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(mu) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(nu) << endl;
      ss << this->errorMessageAtomB << AtomTypeStr(atomB.GetAtomType()) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(lambda) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(sigma) << endl;
      throw MolDSException(ss.str());
   }
   else{
      value = 0.0;
   }
   return value;
}

// Second derivative of NDDO repulsion integral.
// Both derivation are related to the coordinate of atomA
// See Apendix in [DT_1977]
// Orbital mu and nu belong atom A, 
// orbital lambda and sigma belong atomB.
double Mndo::GetNddoRepulsionIntegral2ndDerivative(
                                       const Atom& atomA, OrbitalType mu, OrbitalType nu,
                                       const Atom& atomB, OrbitalType lambda, OrbitalType sigma,
                                       CartesianType axisA1,
                                       CartesianType axisA2) const{
   double value = 0.0;
   double DA=0.0;
   double DB=0.0;
   double rhoA = 0.0;
   double rhoB = 0.0;
   double rAB = this->molecule->GetDistanceAtoms(atomA, atomB);
   double cartesian[CartesianType_end] = {atomA.GetXyz()[XAxis] - atomB.GetXyz()[XAxis], 
                                          atomA.GetXyz()[YAxis] - atomB.GetXyz()[YAxis],
                                          atomA.GetXyz()[ZAxis] - atomB.GetXyz()[ZAxis]};
   double deriv1st=0.0; // first derivative of semi empirical multipole interaction
   double deriv2nd=0.0; // second derivative of semi empirical multipole interaction
   int lA = 0;
   int lB = 0;
   // (28) in [DT_1977]
   if(mu == s && nu == s && lambda == s && sigma == s){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      value = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                   deriv2nd,
                                                                   axisA1,
                                                                   axisA2,
                                                                   cartesian,
                                                                   rAB);
   }
   // (29) in [DT_1977]
   else if(mu == s && nu == s && lambda == px && sigma == px){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qxx, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, Qxx, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2;
   }
   else if(mu == s && nu == s && lambda == py && sigma == py){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qyy, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, Qyy, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2;
   }
   // (30) in [DT_1977]
   else if(mu == s && nu == s && lambda == pz && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qzz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, Qzz, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2;
   }
   // (31) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == s){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxx, sQ, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2;
   }
   else if(mu == py && nu == py && lambda == s && sigma == s){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qyy, sQ, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2;
   }
   // (32) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == s){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qzz, sQ, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2;
   }
   // (33) in [DT_1977]
   else if(mu == px && nu == px && lambda == px && sigma == px){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qxx, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, Qxx, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxx, sQ, rAB);
      double temp3 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, Qxx, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxx, Qxx, rAB);
      double temp4 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == py && sigma == py){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qyy, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, Qyy, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qyy, sQ, rAB);
      double temp3 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, Qyy, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qyy, Qyy, rAB);
      double temp4 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      value = temp1 + temp2 + temp3 + temp4;
   }
   // (34) in [DT_1977]
   else if(mu == px && nu == px && lambda == py && sigma == py){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qyy, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, Qyy, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxx, sQ, rAB);
      double temp3 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, Qyy, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxx, Qyy, rAB);
      double temp4 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == px && sigma == px){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qxx, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, Qxx, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qyy, sQ, rAB);
      double temp3 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, Qxx, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qyy, Qxx, rAB);
      double temp4 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      value = temp1 + temp2 + temp3 + temp4;
   }
   // (35) in [DT_1977]
   else if(mu == px && nu == px && lambda == pz && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qzz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, Qzz, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxx, sQ, rAB);
      double temp3 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, Qzz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxx, Qzz, rAB);
      double temp4 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qzz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, Qzz, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qyy, sQ, rAB);
      double temp3 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, Qzz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qyy, Qzz, rAB);
      double temp4 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      value = temp1 + temp2 + temp3 + temp4;
   }
   // (36) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == px && sigma == px){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qxx, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, Qxx, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qzz, sQ, rAB);
      double temp3 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, Qxx, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qzz, Qxx, rAB);
      double temp4 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == pz && nu == pz && lambda == py && sigma == py){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, Qyy, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, Qyy, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qzz, sQ, rAB);
      double temp3 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, Qyy, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qzz, Qyy, rAB);
      double temp4 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      value = temp1 + temp2 + temp3 + temp4;
   }
   // (37) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == pz && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, 
                                           sQ, Qzz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, Qzz, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, 
                                           Qzz, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qzz, sQ, rAB);
      double temp3 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, Qzz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qzz, Qzz, rAB);
      double temp4 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      value = temp1 + temp2 + temp3 + temp4;
   }
   // (38) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == s){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, muz, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   // (39) in [DT_1977]
   else if(mu == s && nu == pz && lambda == px && sigma == px){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, muz, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, Qxx, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, muz, Qxx, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == px && sigma == px){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == s && nu == pz && lambda == py && sigma == py){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, muz, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, Qyy, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, muz, Qyy, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == py && sigma == py){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   // (40) in [DT_1977]
   else if(mu == s && nu == pz && lambda == pz && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, sQ, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, muz, sQ, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, Qzz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, muz, Qzz, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == pz){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   // (41) in [DT_1977]
   else if(mu == s && nu == s && lambda == s && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, muz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, muz, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1;
   }
   else if(mu == s && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (42) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, muz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, muz, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, muz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxx, muz, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2;
   }
   else if(mu == px && nu == px && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == py && nu == py && lambda == s && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, muz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, muz, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyy, muz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qyy, muz, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (43) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, sQ, muz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, sQ, muz, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);

      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qzz, muz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qzz, muz, rAB);
      double temp2 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (44) in [DT_1977]
   else if(mu == s && nu == px && lambda == s && sigma == px){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, mux, mux, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, mux, mux, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1;
   }
   else if(mu == px && nu == s && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == s && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == px && nu == s && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == s && nu == py && lambda == s && sigma == py){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muy, muy, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, muy, muy, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1;
   }
   else if(mu == py && nu == s && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == s && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == py && nu == s && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (45) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, muz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, muz, muz, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                           rAB);
      value = temp1;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == pz){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == s && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (46) in [DT_1977]
   else if(mu == s && nu == px && lambda == px && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, mux, Qxz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, mux, Qxz, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1;
   }
   else if(mu == px && nu == s && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == s && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == px && nu == s && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == s && nu == py && lambda == py && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muy, Qyz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, muy, Qyz, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1;
   }
   else if(mu == py && nu == s && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == s && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == py && nu == s && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (47) in [DT_1977]
   else if(mu == px && nu == pz && lambda == s && sigma == px){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxz, mux, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxz, mux, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1;
   }
   else if(mu == pz && nu == px && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == px && nu == pz && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == pz && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == py && nu == pz && lambda == s && sigma == py){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyz, muy, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qyz, muy, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1;
   }
   else if(mu == pz && nu == py && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == pz && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (48) in [DT_1977]
   else if(mu == px && nu == pz && lambda == px && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxz, Qxz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxz, Qxz, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1;
   }
   else if(mu == pz && nu == px && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == px && nu == pz && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == pz && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == pz){
      deriv1st = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qyz, Qyz, rAB);
      deriv2nd = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qyz, Qyz, rAB);
      double temp1 = this->Get2ndDerivativeElementFromDistanceDerivatives(deriv1st,
                                                                          deriv2nd,
                                                                          axisA1,
                                                                          axisA2,
                                                                          cartesian,
                                                                          rAB);
      value = temp1;
   }
   else if(mu == pz && nu == py && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == py && nu == pz && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == pz && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (49) in [DT_1977] and p19 in [MOPAC_1990]
   else if(mu == px && nu == py && lambda == px && sigma == py){
      value = 0.5*(this->GetNddoRepulsionIntegral2ndDerivative(
                         atomA, mu, mu, atomB, mu, mu, axisA1, axisA2)
                  -this->GetNddoRepulsionIntegral2ndDerivative(
                         atomA, mu, mu, atomB, nu, nu, axisA1, axisA2));
   }
   else if(mu == py && nu == px && lambda == px && sigma == py){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == px && nu == py && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == py && nu == px && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegral2ndDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // d-orbitals
   else if(mu == dxy || mu == dyz || mu == dzz || mu == dzx || mu == dxxyy ||
           nu == dxy || nu == dyz || nu == dzz || nu == dzx || nu == dxxyy ||
           lambda == dxy || lambda == dyz || lambda == dzz || lambda  == dzx || lambda == dxxyy ||
           sigma == dxy || sigma == dyz || sigma == dzz || sigma  == dzx || sigma == dxxyy){

      stringstream ss;
      ss << this->errorMessageGetNddoRepulsionIntegral2ndDerivative;
      ss << this->errorMessageAtomA << AtomTypeStr(atomA.GetAtomType()) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(mu) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(nu) << endl;
      ss << this->errorMessageAtomB << AtomTypeStr(atomB.GetAtomType()) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(lambda) << endl;
      ss << "\t" << this->errorMessageOrbitalType << OrbitalTypeStr(sigma) << endl;
      throw MolDSException(ss.str());
   }
   else{
      value = 0.0;
   }
   return value;
}

// See Apendix in [DT_1977]
double Mndo::GetSemiEmpiricalMultipoleInteraction(const Atom& atomA,
                                                  const Atom& atomB,
                                                  MultipoleType multipoleA,
                                                  MultipoleType multipoleB,
                                                  double rAB) const{
   double value = 0.0;
   double DA = atomA.GetNddoDerivedParameterD(this->theory, multipoleA);
   double DB = atomB.GetNddoDerivedParameterD(this->theory, multipoleB);
   double rhoA = atomA.GetNddoDerivedParameterRho(this->theory, multipoleA);
   double rhoB = atomB.GetNddoDerivedParameterRho(this->theory, multipoleB);
   double a = rhoA + rhoB;

   // Eq. (52) in [DT_1977]
   if(multipoleA == sQ && multipoleB == sQ){
      value = pow(pow(rAB,2.0) + pow(a,2.0), -0.5);
   }
   // Eq. (53) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == muz){
      double temp1 = pow(rAB+DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB-DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/2.0 - pow(temp2,-0.5)/2.0;
   }
   else if(multipoleA == muz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,1.0);
   }
   // Eq. (54) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == Qxx){
      double temp1 = pow(rAB,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/2.0 - pow(temp2,-0.5)/2.0;
   }
   else if(multipoleA == Qxx && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == sQ && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, multipoleA, Qxx, rAB);
   }
   else if(multipoleA == Qyy && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,2.0);
   }
   // Eq. (55) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == Qzz){
      double temp1 = pow(rAB+2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB-2.0*DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 - pow(temp2,-0.5)/2.0 + pow(temp3,-0.5)/4.0;
   }
   else if(multipoleA == Qzz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,2.0);
   }
   // Eq. (56) in [DT_1977]
   else if(multipoleA == mux && multipoleB == mux){
      double temp1 = pow(rAB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/2.0 - pow(temp2,-0.5)/2.0;
   }
   else if(multipoleA == muy && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, mux, mux, rAB);
   }
   // Eq. (57) in [DT_1977]
   else if(multipoleA == muz && multipoleB == muz){
      double temp1 = pow(rAB+DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB+DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB-DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB-DA+DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 - pow(temp2,-0.5)/4.0 
             -pow(temp3,-0.5)/4.0 + pow(temp4,-0.5)/4.0;
   }
   // Eq. (58) in [DT_1977]
   else if(multipoleA == mux && multipoleB == Qxz){
      double temp1 = pow(rAB-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value =-pow(temp1,-0.5)/4.0 + pow(temp2,-0.5)/4.0 
             +pow(temp3,-0.5)/4.0 - pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == Qxz && multipoleB == mux){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muy && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, mux, Qxz, rAB);
   }
   else if(multipoleA == Qyz && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   // Eq. (59) in [DT_1977]
   else if(multipoleA == muz && multipoleB == Qxx){
      double temp1 = pow(rAB+DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB-DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB+DA,2.0) + pow(a,2.0);
      double temp4 = pow(rAB-DA,2.0) + pow(a,2.0);
      value =-pow(temp1,-0.5)/4.0 + pow(temp2,-0.5)/4.0 
             +pow(temp3,-0.5)/4.0 - pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == Qxx && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, muz, Qxx, rAB);
   }
   else if(multipoleA == Qyy && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   // Eq. (60) in [DT_1977]
   else if(multipoleA == muz && multipoleB == Qzz){
      double temp1 = pow(rAB+DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB-DA-2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB+DA+2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB-DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(rAB+DA,2.0) + pow(a,2.0);
      double temp6 = pow(rAB-DA,2.0) + pow(a,2.0);
      value =-pow(temp1,-0.5)/8.0 + pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/8.0 + pow(temp4,-0.5)/8.0
             +pow(temp5,-0.5)/4.0 - pow(temp6,-0.5)/4.0;
   }
   else if(multipoleA == Qzz && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   // Eq. (61) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qxx){
      double temp1 = pow(rAB,2.0) + 4.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB,2.0) + 4.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp4 = pow(rAB,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(rAB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/8.0 + pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/4.0 - pow(temp4,-0.5)/4.0
             +pow(temp5,-0.5)/4.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxx, Qxx, rAB);
   }
   // Eq. (62) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qyy){
      double temp1 = pow(rAB,2.0) + pow(2.0*DA,2.0) + pow(2.0*DB,2.0)+ pow(a,2.0);
      double temp2 = pow(rAB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(rAB,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 - pow(temp2,-0.5)/4.0 
             -pow(temp3,-0.5)/4.0 + pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,4.0);
   }
   // Eq. (63) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qzz){
      double temp1 = pow(rAB-2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp2 = pow(rAB+2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(rAB-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(rAB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(rAB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/8.0 + pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/8.0 - pow(temp4,-0.5)/8.0
             -pow(temp5,-0.5)/4.0 + pow(temp6,-0.5)/4.0;
   }
   else if(multipoleA == Qzz && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qyy && multipoleB == Qzz){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxx, multipoleB, rAB);
   }
   else if(multipoleA == Qzz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,4.0);
   }
   // Eq. (64) in [DT_1977]
   else if(multipoleA == Qzz && multipoleB == Qzz){
      double temp1 = pow(rAB+2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB+2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB-2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB-2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(rAB+2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(rAB-2.0*DA,2.0) + pow(a,2.0);
      double temp7 = pow(rAB+2.0*DB,2.0) + pow(a,2.0);
      double temp8 = pow(rAB-2.0*DB,2.0) + pow(a,2.0);
      double temp9 = pow(rAB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/16.0 + pow(temp2,-0.5)/16.0 
             +pow(temp3,-0.5)/16.0 + pow(temp4,-0.5)/16.0
             -pow(temp5,-0.5)/8.0 - pow(temp6,-0.5)/8.0
             -pow(temp7,-0.5)/8.0 - pow(temp8,-0.5)/8.0
             +pow(temp9,-0.5)/4.0;
   }
   // Eq. (65) in [DT_1977]
   else if(multipoleA == Qxz && multipoleB == Qxz){
      double temp1 = pow(rAB+DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB+DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB+DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB+DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp5 = pow(rAB-DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp6 = pow(rAB-DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp7 = pow(rAB-DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp8 = pow(rAB-DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/8.0 - pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/8.0 + pow(temp4,-0.5)/8.0
             -pow(temp5,-0.5)/8.0 + pow(temp6,-0.5)/8.0
             +pow(temp7,-0.5)/8.0 - pow(temp8,-0.5)/8.0;
   }
   else if(multipoleA == Qyz && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteraction(atomA, atomB, Qxz, Qxz, rAB);
   }
   // Eq. (66) in [DT_1977]
   else if(multipoleA == Qxy && multipoleB == Qxy){
      double temp1 = pow(rAB,2.0) + 2.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB,2.0) + 2.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB,2.0) + 2.0*pow(DA,2.0) + 2.0*pow(DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 + pow(temp2,-0.5)/4.0 
             -pow(temp3,-0.5)/2.0;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles;
      ss << this->errorMessageMultipoleA << MultipoleTypeStr(multipoleA) << endl;
      ss << this->errorMessageMultipoleB << MultipoleTypeStr(multipoleB) << endl;
      throw MolDSException(ss.str());
   }
   return value;
}

// First derivative of semiempirical multipole-multipole interactions.
// This derivativ is related to the nuclear distance rAB.
// See Apendix in [DT_1977]
double Mndo::GetSemiEmpiricalMultipoleInteraction1stDerivative(const Atom& atomA,
                                                               const Atom& atomB,
                                                               MultipoleType multipoleA,
                                                               MultipoleType multipoleB,
                                                               double rAB) const{
   double value = 0.0;
   double DA = atomA.GetNddoDerivedParameterD(this->theory, multipoleA);
   double DB = atomB.GetNddoDerivedParameterD(this->theory, multipoleB);
   double rhoA = atomA.GetNddoDerivedParameterRho(this->theory, multipoleA);
   double rhoB = atomB.GetNddoDerivedParameterRho(this->theory, multipoleB);
   double a = rhoA + rhoB;

   // Eq. (52) in [DT_1977]
   if(multipoleA == sQ && multipoleB == sQ){
      value = -1.0*rAB*pow(pow(rAB,2.0) + pow(a,2.0), -1.5);
   }
   // Eq. (53) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == muz){
      double temp1 = pow(rAB+DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB-DB,2.0) + pow(a,2.0);
      value = (rAB+DB)*pow(temp1,-1.5)/2.0 
             -(rAB-DB)*pow(temp2,-1.5)/2.0;
      value *= -1.0;
   }
   else if(multipoleA == muz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,1.0);
   }
   // Eq. (54) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == Qxx){
      double temp1 = pow(rAB,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB,2.0) + pow(a,2.0);
      value = rAB*pow(temp1,-1.5)/2.0 
             -rAB*pow(temp2,-1.5)/2.0;
      value *= -1.0;
   }
   else if(multipoleA == Qxx && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == sQ && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, multipoleA, Qxx, rAB);
   }
   else if(multipoleA == Qyy && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,2.0);
   }
   // Eq. (55) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == Qzz){
      double temp1 = pow(rAB+2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB-2.0*DB,2.0) + pow(a,2.0);
      value = (rAB+2.0*DB)*pow(temp1,-1.5)/4.0 
             -(rAB)*pow(temp2,-1.5)/2.0 
             +(rAB-2.0*DB)*pow(temp3,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qzz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,2.0);
   }
   // Eq. (56) in [DT_1977]
   else if(multipoleA == mux && multipoleB == mux){
      double temp1 = pow(rAB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = (rAB)*pow(temp1,-1.5)/2.0 
             -(rAB)*pow(temp2,-1.5)/2.0;
      value *= -1.0;
   }
   else if(multipoleA == muy && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, mux, mux, rAB);
   }
   // Eq. (57) in [DT_1977]
   else if(multipoleA == muz && multipoleB == muz){
      double temp1 = pow(rAB+DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB+DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB-DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB-DA+DB,2.0) + pow(a,2.0);
      value = (rAB+DA-DB)*pow(temp1,-1.5)/4.0 
             -(rAB+DA+DB)*pow(temp2,-1.5)/4.0 
             -(rAB-DA-DB)*pow(temp3,-1.5)/4.0 
             +(rAB-DA+DB)*pow(temp4,-1.5)/4.0;
      value *= -1.0;
   }
   // Eq. (58) in [DT_1977]
   else if(multipoleA == mux && multipoleB == Qxz){
      double temp1 = pow(rAB-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value =-(rAB-DB)*pow(temp1,-1.5)/4.0 
             +(rAB-DB)*pow(temp2,-1.5)/4.0 
             +(rAB+DB)*pow(temp3,-1.5)/4.0 
             -(rAB+DB)*pow(temp4,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qxz && multipoleB == mux){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muy && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, mux, Qxz, rAB);
   }
   else if(multipoleA == Qyz && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   // Eq. (59) in [DT_1977]
   else if(multipoleA == muz && multipoleB == Qxx){
      double temp1 = pow(rAB+DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB-DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB+DA,2.0) + pow(a,2.0);
      double temp4 = pow(rAB-DA,2.0) + pow(a,2.0);
      value =-(rAB+DA)*pow(temp1,-1.5)/4.0 
             +(rAB-DA)*pow(temp2,-1.5)/4.0 
             +(rAB+DA)*pow(temp3,-1.5)/4.0 
             -(rAB-DA)*pow(temp4,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qxx && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, muz, Qxx, rAB);
   }
   else if(multipoleA == Qyy && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   // Eq. (60) in [DT_1977]
   else if(multipoleA == muz && multipoleB == Qzz){
      double temp1 = pow(rAB+DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB-DA-2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB+DA+2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB-DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(rAB+DA,2.0) + pow(a,2.0);
      double temp6 = pow(rAB-DA,2.0) + pow(a,2.0);
      value =-(rAB+DA-2.0*DB)*pow(temp1,-1.5)/8.0 
             +(rAB-DA-2.0*DB)*pow(temp2,-1.5)/8.0 
             -(rAB+DA+2.0*DB)*pow(temp3,-1.5)/8.0 
             +(rAB-DA+2.0*DB)*pow(temp4,-1.5)/8.0
             +(rAB+DA       )*pow(temp5,-1.5)/4.0 
             -(rAB-DA       )*pow(temp6,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qzz && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   // Eq. (61) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qxx){
      double temp1 = pow(rAB,2.0) + 4.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB,2.0) + 4.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp4 = pow(rAB,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(rAB,2.0) + pow(a,2.0);
      value = (rAB)*pow(temp1,-1.5)/8.0 
             +(rAB)*pow(temp2,-1.5)/8.0 
             -(rAB)*pow(temp3,-1.5)/4.0 
             -(rAB)*pow(temp4,-1.5)/4.0
             +(rAB)*pow(temp5,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, Qxx, rAB);
   }
   // Eq. (62) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qyy){
      double temp1 = pow(rAB,2.0) + pow(2.0*DA,2.0) + pow(2.0*DB,2.0)+ pow(a,2.0);
      double temp2 = pow(rAB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(rAB,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB,2.0) + pow(a,2.0);
      value = (rAB)*pow(temp1,-1.5)/4.0 
             -(rAB)*pow(temp2,-1.5)/4.0 
             -(rAB)*pow(temp3,-1.5)/4.0 
             +(rAB)*pow(temp4,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,4.0);
   }
   // Eq. (63) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qzz){
      double temp1 = pow(rAB-2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp2 = pow(rAB+2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(rAB-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(rAB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(rAB,2.0) + pow(a,2.0);
      value = (rAB-2.0*DB)*pow(temp1,-1.5)/8.0 
             +(rAB+2.0*DB)*pow(temp2,-1.5)/8.0 
             -(rAB-2.0*DB)*pow(temp3,-1.5)/8.0 
             -(rAB+2.0*DB)*pow(temp4,-1.5)/8.0
             -(rAB       )*pow(temp5,-1.5)/4.0 
             +(rAB       )*pow(temp6,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qzz && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qyy && multipoleB == Qzz){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxx, multipoleB, rAB);
   }
   else if(multipoleA == Qzz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,4.0);
   }
   // Eq. (64) in [DT_1977]
   else if(multipoleA == Qzz && multipoleB == Qzz){
      double temp1 = pow(rAB+2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB+2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB-2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB-2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(rAB+2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(rAB-2.0*DA,2.0) + pow(a,2.0);
      double temp7 = pow(rAB+2.0*DB,2.0) + pow(a,2.0);
      double temp8 = pow(rAB-2.0*DB,2.0) + pow(a,2.0);
      double temp9 = pow(rAB,2.0) + pow(a,2.0);
      value = (rAB+2.0*DA-2.0*DB)*pow(temp1,-1.5)/16.0 
             +(rAB+2.0*DA+2.0*DB)*pow(temp2,-1.5)/16.0 
             +(rAB-2.0*DA-2.0*DB)*pow(temp3,-1.5)/16.0 
             +(rAB-2.0*DA+2.0*DB)*pow(temp4,-1.5)/16.0
             -(rAB+2.0*DA)*pow(temp5,-1.5)/8.0 
             -(rAB-2.0*DA)*pow(temp6,-1.5)/8.0
             -(rAB+2.0*DB)*pow(temp7,-1.5)/8.0 
             -(rAB-2.0*DB)*pow(temp8,-1.5)/8.0
             +(rAB)*pow(temp9,-1.5)/4.0;
      value *= -1.0;
   }
   // Eq. (65) in [DT_1977]
   else if(multipoleA == Qxz && multipoleB == Qxz){
      double temp1 = pow(rAB+DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB+DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB+DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(rAB+DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp5 = pow(rAB-DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp6 = pow(rAB-DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp7 = pow(rAB-DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp8 = pow(rAB-DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = (rAB+DA-DB)*pow(temp1,-1.5)/8.0 
             -(rAB+DA-DB)*pow(temp2,-1.5)/8.0 
             -(rAB+DA+DB)*pow(temp3,-1.5)/8.0 
             +(rAB+DA+DB)*pow(temp4,-1.5)/8.0
             -(rAB-DA-DB)*pow(temp5,-1.5)/8.0 
             +(rAB-DA-DB)*pow(temp6,-1.5)/8.0
             +(rAB-DA+DB)*pow(temp7,-1.5)/8.0 
             -(rAB-DA+DB)*pow(temp8,-1.5)/8.0;
      value *= -1.0;
   }
   else if(multipoleA == Qyz && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteraction1stDerivative(atomA, atomB, Qxz, Qxz, rAB);
   }
   // Eq. (66) in [DT_1977]
   else if(multipoleA == Qxy && multipoleB == Qxy){
      double temp1 = pow(rAB,2.0) + 2.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(rAB,2.0) + 2.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(rAB,2.0) + 2.0*pow(DA,2.0) + 2.0*pow(DB,2.0) + pow(a,2.0);
      value = (rAB)*pow(temp1,-1.5)/4.0 
             +(rAB)*pow(temp2,-1.5)/4.0 
             -(rAB)*pow(temp3,-1.5)/2.0;
      value *= -1.0;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetSemiEmpiricalMultipoleInteraction1stDeriBadMultipoles;
      ss << this->errorMessageMultipoleA << MultipoleTypeStr(multipoleA) << endl;
      ss << this->errorMessageMultipoleB << MultipoleTypeStr(multipoleB) << endl;
      throw MolDSException(ss.str());
   }
   return value;
}

// Second derivative of semiempirical multipole-multipole interactions.
// This derivativ is related to the nuclear distance rAB.
// See Apendix in [DT_1977]
double Mndo::GetSemiEmpiricalMultipoleInteraction2ndDerivative(const Atom& atomA,
                                                                  const Atom& atomB,
                                                                  MultipoleType multipoleA,
                                                                  MultipoleType multipoleB,
                                                                  double rAB) const{
   double value = 0.0;
   double DA = atomA.GetNddoDerivedParameterD(this->theory, multipoleA);
   double DB = atomB.GetNddoDerivedParameterD(this->theory, multipoleB);
   double rhoA = atomA.GetNddoDerivedParameterRho(this->theory, multipoleA);
   double rhoB = atomB.GetNddoDerivedParameterRho(this->theory, multipoleB);
   double a = rhoA + rhoB;

   // Eq. (52) in [DT_1977]
   if(multipoleA == sQ && multipoleB == sQ){
      double c1 = 1.0;
      double f1 = pow(rAB,2.0);
      double a1 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
   }
   // Eq. (53) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == muz){
      double c1 = 0.5;
      double c2 = -0.5;
      double f1 = pow(rAB+DB,2.0);
      double f2 = pow(rAB-DB,2.0);
      double a1 = pow(a,2.0);
      double a2 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
   }
   else if(multipoleA == muz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,1.0);
   }
   // Eq. (54) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == Qxx){
      double c1 = 0.5;
      double c2 = -0.5;
      double f1 = pow(rAB,2.0);
      double f2 = pow(rAB,2.0);
      double a1 = pow(2.0*DB,2.0) + pow(a,2.0);
      double a2 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
   }
   else if(multipoleA == Qxx && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == sQ && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, multipoleA, Qxx, rAB);
   }
   else if(multipoleA == Qyy && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,2.0);
   }
   // Eq. (55) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == Qzz){
      double c1 = 0.25;
      double c2 = -0.50;
      double c3 = 0.25;
      double f1 = pow(rAB+2.0*DB,2.0);
      double f2 = pow(rAB,2.0);
      double f3 = pow(rAB-2.0*DB,2.0);
      double a1 = pow(a,2.0);
      double a2 = pow(a,2.0);
      double a3 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
   }
   else if(multipoleA == Qzz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,2.0);
   }
   // Eq. (56) in [DT_1977]
   else if(multipoleA == mux && multipoleB == mux){
      double c1 = 0.50;
      double c2 = -0.50;
      double f1 = pow(rAB,2.0);
      double f2 = pow(rAB,2.0);
      double a1 = pow(DA-DB,2.0) + pow(a,2.0);
      double a2 = pow(DA+DB,2.0) + pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
   }
   else if(multipoleA == muy && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, mux, mux, rAB);
   }
   // Eq. (57) in [DT_1977]
   else if(multipoleA == muz && multipoleB == muz){
      double c1 =  0.25;
      double c2 = -0.25;
      double c3 = -0.25;
      double c4 =  0.25;
      double f1 = pow(rAB+DA-DB,2.0);
      double f2 = pow(rAB+DA+DB,2.0);
      double f3 = pow(rAB-DA-DB,2.0);
      double f4 = pow(rAB-DA+DB,2.0);
      double a1 = pow(a,2.0);
      double a2 = pow(a,2.0);
      double a3 = pow(a,2.0);
      double a4 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
      value += c4*(3.0*f4*pow(f4+a4,-2.5) - pow(f4+a4,-1.5));
   }
   // Eq. (58) in [DT_1977]
   else if(multipoleA == mux && multipoleB == Qxz){
      double c1 = -0.25;
      double c2 = 0.25;
      double c3 = 0.25;
      double c4 = -0.25;
      double f1 = pow(rAB-DB,2.0);
      double f2 = pow(rAB-DB,2.0);
      double f3 = pow(rAB+DB,2.0);
      double f4 = pow(rAB+DB,2.0);
      double a1 = pow(DA-DB,2.0) + pow(a,2.0);
      double a2 = pow(DA+DB,2.0) + pow(a,2.0);
      double a3 = pow(DA-DB,2.0) + pow(a,2.0);
      double a4 = pow(DA+DB,2.0) + pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
      value += c4*(3.0*f4*pow(f4+a4,-2.5) - pow(f4+a4,-1.5));
   }
   else if(multipoleA == Qxz && multipoleB == mux){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muy && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, mux, Qxz, rAB);
   }
   else if(multipoleA == Qyz && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   // Eq. (59) in [DT_1977]
   else if(multipoleA == muz && multipoleB == Qxx){
      double c1 = -0.25;
      double c2 = 0.25;
      double c3 = 0.25;
      double c4 = -0.25;
      double f1 = pow(rAB+DA,2.0);
      double f2 = pow(rAB-DA,2.0);
      double f3 = pow(rAB+DA,2.0);
      double f4 = pow(rAB-DA,2.0);
      double a1 = pow(2.0*DB,2.0) + pow(a,2.0);
      double a2 = pow(2.0*DB,2.0) + pow(a,2.0);
      double a3 = pow(a,2.0);
      double a4 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
      value += c4*(3.0*f4*pow(f4+a4,-2.5) - pow(f4+a4,-1.5));
   }
   else if(multipoleA == Qxx && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, muz, Qxx, rAB);
   }
   else if(multipoleA == Qyy && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   // Eq. (60) in [DT_1977]
   else if(multipoleA == muz && multipoleB == Qzz){
      double c1 = -0.125;
      double c2 =  0.125;
      double c3 = -0.125;
      double c4 =  0.125;
      double c5 =  0.25;
      double c6 = -0.25;
      double f1 = pow(rAB+DA-2.0*DB,2.0);
      double f2 = pow(rAB-DA-2.0*DB,2.0);
      double f3 = pow(rAB+DA+2.0*DB,2.0);
      double f4 = pow(rAB-DA+2.0*DB,2.0);
      double f5 = pow(rAB+DA,2.0);
      double f6 = pow(rAB-DA,2.0);
      double a1 = pow(a,2.0);
      double a2 = pow(a,2.0);
      double a3 = pow(a,2.0);
      double a4 = pow(a,2.0);
      double a5 = pow(a,2.0);
      double a6 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
      value += c4*(3.0*f4*pow(f4+a4,-2.5) - pow(f4+a4,-1.5));
      value += c5*(3.0*f5*pow(f5+a5,-2.5) - pow(f5+a5,-1.5));
      value += c6*(3.0*f6*pow(f6+a6,-2.5) - pow(f6+a6,-1.5));
   }
   else if(multipoleA == Qzz && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,3.0);
   }
   // Eq. (61) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qxx){
      double c1 =  0.125;
      double c2 =  0.125;
      double c3 = -0.25;
      double c4 = -0.25;
      double c5 =  0.25;
      double f1 = pow(rAB,2.0);
      double f2 = pow(rAB,2.0);
      double f3 = pow(rAB,2.0);
      double f4 = pow(rAB,2.0);
      double f5 = pow(rAB,2.0);
      double a1 = 4.0*pow(DA-DB,2.0) + pow(a,2.0);
      double a2 = 4.0*pow(DA+DB,2.0) + pow(a,2.0);
      double a3 = pow(2.0*DA,2.0)    + pow(a,2.0);
      double a4 = pow(2.0*DB,2.0)    + pow(a,2.0);
      double a5 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
      value += c4*(3.0*f4*pow(f4+a4,-2.5) - pow(f4+a4,-1.5));
      value += c5*(3.0*f5*pow(f5+a5,-2.5) - pow(f5+a5,-1.5));
   }
   else if(multipoleA == Qyy && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxx, Qxx, rAB);
   }
   // Eq. (62) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qyy){
      double c1 =  0.25;
      double c2 = -0.25;
      double c3 = -0.25;
      double c4 =  0.25;
      double f1 = pow(rAB,2.0);
      double f2 = pow(rAB,2.0);
      double f3 = pow(rAB,2.0);
      double f4 = pow(rAB,2.0);
      double a1 = pow(2.0*DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double a2 = pow(2.0*DA,2.0) +                   pow(a,2.0);
      double a3 = pow(2.0*DB,2.0) +                   pow(a,2.0);
      double a4 =                                     pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
      value += c4*(3.0*f4*pow(f4+a4,-2.5) - pow(f4+a4,-1.5));
   }
   else if(multipoleA == Qyy && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,4.0);
   }
   // Eq. (63) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qzz){
      double c1 =  0.125;
      double c2 =  0.125;
      double c3 = -0.125;
      double c4 = -0.125;
      double c5 = -0.25;
      double c6 =  0.25;
      double f1 = pow(rAB-2.0*DB,2.0);
      double f2 = pow(rAB+2.0*DB,2.0);
      double f3 = pow(rAB-2.0*DB,2.0);
      double f4 = pow(rAB+2.0*DB,2.0);
      double f5 = pow(rAB       ,2.0);
      double f6 = pow(rAB       ,2.0);
      double a1 = pow(2.0*DA,2.0) + pow(a,2.0);
      double a2 = pow(2.0*DA,2.0) + pow(a,2.0);
      double a3 =                   pow(a,2.0);
      double a4 =                   pow(a,2.0);
      double a5 = pow(2.0*DA,2.0) + pow(a,2.0);
      double a6 =                   pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
      value += c4*(3.0*f4*pow(f4+a4,-2.5) - pow(f4+a4,-1.5));
      value += c5*(3.0*f5*pow(f5+a5,-2.5) - pow(f5+a5,-1.5));
      value += c6*(3.0*f6*pow(f6+a6,-2.5) - pow(f6+a6,-1.5));
   }
   else if(multipoleA == Qzz && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qyy && multipoleB == Qzz){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxx, multipoleB, rAB);
   }
   else if(multipoleA == Qzz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomB, atomA, multipoleB, multipoleA, rAB);
      value *= pow(-1.0,4.0);
   }
   // Eq. (64) in [DT_1977]
   else if(multipoleA == Qzz && multipoleB == Qzz){
      double c1 =  0.0625;
      double c2 =  0.0625;
      double c3 =  0.0625;
      double c4 =  0.0625;
      double c5 = -0.125;
      double c6 = -0.125;
      double c7 = -0.125;
      double c8 = -0.125;
      double c9 =  0.25;
      double f1 = pow(rAB+2.0*DA-2.0*DB,2.0);
      double f2 = pow(rAB+2.0*DA+2.0*DB,2.0);
      double f3 = pow(rAB-2.0*DA-2.0*DB,2.0);
      double f4 = pow(rAB-2.0*DA+2.0*DB,2.0);
      double f5 = pow(rAB+2.0*DA       ,2.0);
      double f6 = pow(rAB-2.0*DA       ,2.0);
      double f7 = pow(rAB+2.0*DB       ,2.0);
      double f8 = pow(rAB-2.0*DB       ,2.0);
      double f9 = pow(rAB              ,2.0);
      double a1 = pow(a,2.0);
      double a2 = pow(a,2.0);
      double a3 = pow(a,2.0);
      double a4 = pow(a,2.0);
      double a5 = pow(a,2.0);
      double a6 = pow(a,2.0);
      double a7 = pow(a,2.0);
      double a8 = pow(a,2.0);
      double a9 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
      value += c4*(3.0*f4*pow(f4+a4,-2.5) - pow(f4+a4,-1.5));
      value += c5*(3.0*f5*pow(f5+a5,-2.5) - pow(f5+a5,-1.5));
      value += c6*(3.0*f6*pow(f6+a6,-2.5) - pow(f6+a6,-1.5));
      value += c7*(3.0*f7*pow(f7+a7,-2.5) - pow(f7+a7,-1.5));
      value += c8*(3.0*f8*pow(f8+a8,-2.5) - pow(f8+a8,-1.5));
      value += c9*(3.0*f9*pow(f9+a9,-2.5) - pow(f9+a9,-1.5));
   }
   // Eq. (65) in [DT_1977]
   else if(multipoleA == Qxz && multipoleB == Qxz){
      double c1 =  0.125;
      double c2 = -0.125;
      double c3 = -0.125;
      double c4 =  0.125;
      double c5 = -0.125;
      double c6 =  0.125;
      double c7 =  0.125;
      double c8 = -0.125;
      double f1 = pow(rAB+DA-DB,2.0);
      double f2 = pow(rAB+DA-DB,2.0);
      double f3 = pow(rAB+DA+DB,2.0);
      double f4 = pow(rAB+DA+DB,2.0);
      double f5 = pow(rAB-DA-DB,2.0);
      double f6 = pow(rAB-DA-DB,2.0);
      double f7 = pow(rAB-DA+DB,2.0);
      double f8 = pow(rAB-DA+DB,2.0);
      double a1 = pow(DA-DB,2.0) + pow(a,2.0);
      double a2 = pow(DA+DB,2.0) + pow(a,2.0);
      double a3 = pow(DA-DB,2.0) + pow(a,2.0);
      double a4 = pow(DA+DB,2.0) + pow(a,2.0);
      double a5 = pow(DA-DB,2.0) + pow(a,2.0);
      double a6 = pow(DA+DB,2.0) + pow(a,2.0);
      double a7 = pow(DA-DB,2.0) + pow(a,2.0);
      double a8 = pow(DA+DB,2.0) + pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
      value += c4*(3.0*f4*pow(f4+a4,-2.5) - pow(f4+a4,-1.5));
      value += c5*(3.0*f5*pow(f5+a5,-2.5) - pow(f5+a5,-1.5));
      value += c6*(3.0*f6*pow(f6+a6,-2.5) - pow(f6+a6,-1.5));
      value += c7*(3.0*f7*pow(f7+a7,-2.5) - pow(f7+a7,-1.5));
      value += c8*(3.0*f8*pow(f8+a8,-2.5) - pow(f8+a8,-1.5));
   }
   else if(multipoleA == Qyz && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteraction2ndDerivative(atomA, atomB, Qxz, Qxz, rAB);
   }
   // Eq. (66) in [DT_1977]
   else if(multipoleA == Qxy && multipoleB == Qxy){
      double c1 =  0.25;
      double c2 =  0.25;
      double c3 = -0.50;
      double f1 = pow(rAB,2.0);
      double f2 = pow(rAB,2.0);
      double f3 = pow(rAB,2.0);
      double a1 = 2.0*pow(DA-DB,2.0) + pow(a,2.0);
      double a2 = 2.0*pow(DA+DB,2.0) + pow(a,2.0);
      double a3 = 2.0*pow(DA,2.0) + 2.0*pow(DB,2.0) + pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetSemiEmpiricalMultipoleInteraction2ndDeriBadMultipoles;
      ss << this->errorMessageMultipoleA << MultipoleTypeStr(multipoleA) << endl;
      ss << this->errorMessageMultipoleB << MultipoleTypeStr(multipoleB) << endl;
      throw MolDSException(ss.str());
   }
   return value;
}

}
