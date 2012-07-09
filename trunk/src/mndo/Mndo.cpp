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
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/Uncopyable.h"
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
 *  Main Refferences for MNDO are [DT_1977, DT_1977-2, DT_1977-3]
 */
Mndo::Mndo() : MolDS_zindo::ZindoS(){
   this->theory = MNDO;
   this->SetMessages();
   this->SetEnableAtomTypes();
   this->heatsFormation = 0.0;
   this->zMatrixForceElecStatesNum = 0;
   this->etaMatrixForceElecStatesNum = 0;
   this->zMatrixForce = NULL;
   this->etaMatrixForce = NULL;
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
   MallocerFreer::GetInstance()->Free<double>(&this->zMatrixForce, 
                                              this->zMatrixForceElecStatesNum,
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Free<double>(&this->etaMatrixForce, 
                                              this->etaMatrixForceElecStatesNum,
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs());
}

void Mndo::SetMolecule(Molecule* molecule){
   Cndo2::SetMolecule(molecule);
   MallocerFreer::GetInstance()->Malloc<double>(&this->twoElecTwoCore,
                                                molecule->GetNumberAtoms(),
                                                molecule->GetNumberAtoms(),
                                                dxy, dxy, dxy, dxy);
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
   this->errorMessageGetSemiEmpiricalMultipoleInteractionFirstDeriBadMultipoles
      = "Error in mndo:: Mndo::GetSemiEmpiricalMultipoleInteractionFirstDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionSecondDeriBadMultipoles
      = "Error in mndo:: Mndo::GetSemiEmpiricalMultipoleInteractionSecondDerivative: Bad multipole combintaion is set\n";
   this->errorMessageMultipoleA = "Multipole A is: ";
   this->errorMessageMultipoleB = "Multipole B is: ";
   this->errorMessageGetNddoRepulsionIntegral 
      = "Error in mndo::Mndo::GetNddoRepulsionIntegral: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegralFirstDerivative 
      = "Error in mndo::Mndo::GetNddoRepulsionIntegralFirstDerivative: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegralSecondDerivative 
      = "Error in mndo::Mndo::GetNddoRepulsionIntegralSecondDerivative: Bad orbital is set.\n";
   this->errorMessageCalcTwoElecTwoCoreNullMatrix 
      = "Error in mndo::Mndo::CalcTwoElecTwoCore: The two elec two core matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCoreSameAtoms
      = "Error in mndo::Mndo::CalcDiatomicTwoElecTwoCore: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCoreFirstDerivativesSameAtoms
      = "Error in mndo::Mndo::CalcDiatomicTwoElecTwoCoreFirstDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCoreSecondDerivativesSameAtoms
      = "Error in mndo::Mndo::CalcDiatomicTwoElecTwoCoreSecondDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCoreNullMatrix 
      = "Error in mndo::Mndo::CalcDiatomicTwoElecTwoCore: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCoreFirstDerivativesNullMatrix
      = "Error in mndo::Mndo::CalcDiatomicTwoElecTwoCoreFirstDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCoreSecondDerivativesNullMatrix
      = "Error in mndo::Mndo::CalcDiatomicTwoElecTwoCoreSecondDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
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
double Mndo::GetAuxiliaryDiatomCoreRepulsionEnergyFirstDerivative(const Atom& atomA, 
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
double Mndo::GetAuxiliaryDiatomCoreRepulsionEnergySecondDerivative(const Atom& atomA, 
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
      pre2 = pow(dCartesian1/distanceAB,2.0);
   }

   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double alphaA = atomA.GetNddoAlpha(this->theory);
   double alphaB = atomB.GetNddoAlpha(this->theory);
   double fact1=0.0;
   double fact2=0.0;
   if(atomA.GetAtomType() == H && (atomB.GetAtomType() == N || 
                                   atomB.GetAtomType() == O)  ){
      fact1 = -alphaA*exp(-alphaA*distanceAB);
             +((1.0/ang2AU) - alphaB*(distanceAB/ang2AU))*exp(-alphaB*distanceAB);
      fact2 = alphaA*alphaA*exp(-alphaA*distanceAB);
             +(-2.0*alphaA/ang2AU + (distanceAB/ang2AU)*alphaA*alphaA)*exp(-alphaB*distanceAB);
   }
   else if(atomB.GetAtomType() == H && (atomA.GetAtomType() == N || 
                                        atomA.GetAtomType() == O)  ){
      fact1 = -alphaB*exp(-alphaB*distanceAB);
             +((1.0/ang2AU) - alphaA*(distanceAB/ang2AU))*exp(-alphaA*distanceAB);
      fact2 = alphaB*alphaB*exp(-alphaB*distanceAB);
             +(-2.0*alphaB/ang2AU + (distanceAB/ang2AU)*alphaB*alphaB)*exp(-alphaA*distanceAB);
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
double Mndo::GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
                                                   int atomBIndex, 
                                                   CartesianType axisA) const{
   double value =0.0;
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   double distanceAB = this->molecule->GetDistanceAtoms(atomAIndex, atomBIndex);
   double twoElecInt = this->GetNddoRepulsionIntegral(atomA, s, s, atomB, s, s);
   double twoElecIntFirstDeriv = this->GetNddoRepulsionIntegralFirstDerivative(
                                       atomA, s, s, atomB, s, s, axisA);
   double temp = this->GetAuxiliaryDiatomCoreRepulsionEnergy(atomA, atomB, distanceAB);
   double tempDeriv = this->GetAuxiliaryDiatomCoreRepulsionEnergyFirstDerivative(atomA, atomB, distanceAB, axisA);
   value = atomA.GetCoreCharge()*atomB.GetCoreCharge()
          *(twoElecIntFirstDeriv*temp + twoElecInt*tempDeriv); 
   return value;
}

// Second derivative of diatomic core repulsion energy.
// Both derivatives are related to the coordinate of atomA.
double Mndo::GetDiatomCoreRepulsionSecondDerivative(int atomAIndex,
                                                    int atomBIndex, 
                                                    CartesianType axisA1,
                                                    CartesianType axisA2) const{
   double value =0.0;
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   double distanceAB = this->molecule->GetDistanceAtoms(atomAIndex, atomBIndex);
   double twoElecInt = this->GetNddoRepulsionIntegral(atomA, s, s, 
                                                      atomB, s, s);
   double twoElecIntFirstDeriv1 = this->GetNddoRepulsionIntegralFirstDerivative(atomA, s, s, 
                                                                                atomB, s, s, 
                                                                                axisA1);
   double twoElecIntFirstDeriv2 = this->GetNddoRepulsionIntegralFirstDerivative(atomA, s, s, 
                                                                                atomB, s, s, 
                                                                                axisA2);
   double twoElecIntSecondDeriv = this->GetNddoRepulsionIntegralSecondDerivative(atomA, s, s, 
                                                                                 atomB, s, s, 
                                                                                 axisA1, 
                                                                                 axisA2);

   double temp = this->GetAuxiliaryDiatomCoreRepulsionEnergy(atomA, 
                                                             atomB, 
                                                             distanceAB);
   double tempFirstDeriv1 = this->GetAuxiliaryDiatomCoreRepulsionEnergyFirstDerivative(atomA, 
                                                                                       atomB, 
                                                                                       distanceAB, 
                                                                                       axisA1);
   double tempFirstDeriv2 = this->GetAuxiliaryDiatomCoreRepulsionEnergyFirstDerivative(atomA, 
                                                                                       atomB, 
                                                                                       distanceAB, 
                                                                                       axisA2);
   double tempSecondDeriv = this->GetAuxiliaryDiatomCoreRepulsionEnergySecondDerivative(atomA, 
                                                                                        atomB, 
                                                                                        distanceAB, 
                                                                                        axisA1, 
                                                                                        axisA2);

   value = atomA.GetCoreCharge()*atomB.GetCoreCharge();
   value *= twoElecInt*tempSecondDeriv 
           +twoElecIntFirstDeriv1*tempFirstDeriv2 
           +twoElecIntFirstDeriv2*tempFirstDeriv1
           +twoElecIntSecondDeriv*temp;
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

   /*
   // test code for hessian
   int hessianDim = this->molecule->GetNumberAtoms()*3;
   double** hessian = NULL;
   double* forceCons = NULL;
   bool isMassWeighted = true;
   MallocerFreer::GetInstance()->Malloc<double>(&hessian, hessianDim, hessianDim);
   MallocerFreer::GetInstance()->Malloc<double>(&forceCons, hessianDim);
   this->CalcHessianSCF(hessian, isMassWeighted);
   for(int i=0; i<hessianDim; i++){
      for(int j=0; j<hessianDim; j++){
         printf("hess elem: %d %d %e\n",i,j,hessian[i][j]);
      }
   }
   cout << endl << endl;
   bool calcEigenVectors = true;
   MolDS_wrappers::Lapack::GetInstance()->Dsyevd(hessian,
                                                 forceCons,
                                                 hessianDim,
                                                 calcEigenVectors);
   for(int i=0; i<hessianDim; i++){
      printf("force cons: %d %e\n",i,forceCons[i]);
   }
   cout << endl << endl;
   MallocerFreer::GetInstance()->Free<double>(&hessian, hessianDim, hessianDim);
   MallocerFreer::GetInstance()->Free<double>(&forceCons, hessianDim);
   */
}

void Mndo::OutputSCFResults() const{
   MolDS_cndo::Cndo2::OutputSCFResults();
   // output heats of formation
   this->OutputLog(this->messageHeatsFormationTitle);
   this->OutputLog((boost::format("%s\t%e\t%e\n\n") % this->messageHeatsFormation
                                                    % this->heatsFormation 
                                                    % (this->heatsFormation/Parameters::GetInstance()->
                                                                                        GetKcalMolin2AU())).str());
}

double Mndo::GetFockDiagElement(const Atom& atomA, 
                                int atomAIndex, 
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
      for(int B=0; B<molecule.GetNumberAtoms(); B++){
         if(B != atomAIndex){
            const Atom& atomB = *molecule.GetAtom(B);
            int firstAOIndexB = atomB.GetFirstAOIndex();
            for(int lambda=0; lambda<atomB.GetValenceSize(); lambda++){
               for(int sigma=0; sigma<atomB.GetValenceSize(); sigma++){
                  temp += orbitalElectronPopulation[lambda+firstAOIndexB]
                                                   [sigma+firstAOIndexB]
                         *twoElecTwoCore[atomAIndex][B][mu][mu][lambda][sigma];
               }
            }
            temp += this->GetElectronCoreAttraction(atomAIndex, 
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
                                   int atomAIndex, 
                                   int atomBIndex, 
                                   int mu, 
                                   int nu, 
                                   const Molecule& molecule, 
                                   double const* const* gammaAB, 
                                   double const* const* overlap,
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
      value = bondParameter*overlap[mu+firstAOIndexA][nu+firstAOIndexB];
   }
   else{
      double coulomb = 0.0;
      double exchange = 0.0;
      double temp = 0.0;
      if(atomAIndex == atomBIndex){
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalNu, atomA); 
         exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA); 
         temp = (1.5*exchange - 0.5*coulomb)
               *orbitalElectronPopulation[mu+firstAOIndexA][nu+firstAOIndexB];
         for(int BB=0; BB<molecule.GetNumberAtoms(); BB++){
            if(BB != atomAIndex){
               const Atom& atomBB = *molecule.GetAtom(BB);
               int firstAOIndexBB = atomBB.GetFirstAOIndex();
               for(int lambda=0; lambda<atomBB.GetValenceSize(); lambda++){
                  for(int sigma=0; sigma<atomBB.GetValenceSize(); sigma++){
                     temp += orbitalElectronPopulation[lambda+firstAOIndexBB]
                                                      [sigma+firstAOIndexBB]
                            *twoElecTwoCore[atomAIndex][BB][mu][nu][lambda][sigma];
                  }
               }
               temp += this->GetElectronCoreAttraction(atomAIndex, 
                                                       BB, 
                                                       mu, 
                                                       nu, 
                                                       twoElecTwoCore);
            }
         }
      }
      else{
         temp = bondParameter*overlap[mu+firstAOIndexA][nu+firstAOIndexB];
         for(int sigma=0; sigma<atomA.GetValenceSize(); sigma++){
            for(int lambda=0; lambda<atomB.GetValenceSize(); lambda++){
               temp -= 0.5*orbitalElectronPopulation[lambda+firstAOIndexB]
                                                    [sigma+firstAOIndexA]
                      *twoElecTwoCore[atomAIndex][atomBIndex][mu][sigma][nu][lambda];
            }
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
double Mndo::GetElectronCoreAttraction(int atomAIndex, 
                                       int atomBIndex, 
                                       int mu, 
                                       int nu, 
                                       double const* const* const* const* const* const* twoElecTwoCore) const{
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   return -1.0*atomB.GetCoreCharge()*twoElecTwoCore[atomAIndex][atomBIndex][mu][nu][s][s];
}

// First derivative of electron in atom A (mu and nu) and core (atom B) attraction. 
// This derivative is related to the coordinate of atomA.
// Note that diatomicTwoElecTwoCoreFirstDerivative is dioatomic one.
// see Eq. (16) in [DT_1977-2] with f_2 = 0.
double Mndo::GetElectronCoreAttractionFirstDerivative(int atomAIndex, 
                                                      int atomBIndex, 
                                                      int mu, 
                                                      int nu, 
                                                      double const* const* const* const* const* diatomicTwoElecTwoCoreFirstDerivatives,
                                                      CartesianType axisA) const{
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   double value = -1.0*atomB.GetCoreCharge()
                  *diatomicTwoElecTwoCoreFirstDerivatives[mu][nu][s][s][axisA];
   return value;
}

void Mndo::CalcDiatomicOverlapInDiatomicFrame(double** diatomicOverlap, 
                                              const Atom& atomA, 
                                              const Atom& atomB) const{
   MolDS_cndo::Cndo2::CalcDiatomicOverlapInDiatomicFrame(diatomicOverlap, atomA, atomB);
}

// First derivative of (B.40) in J. A. Pople book without bond corrections.
void Mndo::CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(double** diatomicOverlapDeri, 
                                                             const Atom& atomA, 
                                                             const Atom& atomB) const{
   MolDS_cndo::Cndo2::CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(
                      diatomicOverlapDeri,atomA, atomB);
}

// Second derivative of (B.40) in J. A. Pople book without bond corrections.
void Mndo::CalcDiatomicOverlapSecondDerivativeInDiatomicFrame(double** diatomicOverlapSecondDeri, 
                                                              const Atom& atomA, 
                                                              const Atom& atomB) const{
   MolDS_cndo::Cndo2::CalcDiatomicOverlapSecondDerivativeInDiatomicFrame(
                      diatomicOverlapSecondDeri,atomA, atomB);
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
      int numberAOsA = atomA.GetValenceSize();

      for(int B=A; B<molecule.GetNumberAtoms(); B++){
         const Atom& atomB = *molecule.GetAtom(B);
         int firstAOIndexB = atomB.GetFirstAOIndex();
         int numberAOsB = atomB.GetValenceSize();

         double gamma = 0.0;
         if(A!=B){
            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               for(int nu=mu; nu<firstAOIndexA+numberAOsA; nu++){
                  for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                     for(int sigma=lambda; sigma<firstAOIndexB+numberAOsB; sigma++){
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
            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
                  for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                     for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
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

   stringstream ompErrors;
   #pragma omp parallel for schedule(auto)
   for(int k=0; k<this->matrixCISdimension; k++){
      try{
         // single excitation from I-th (occupied)MO to A-th (virtual)MO
         int moI = this->GetActiveOccIndex(*this->molecule, k);
         int moA = this->GetActiveVirIndex(*this->molecule, k);

         for(int l=k; l<this->matrixCISdimension; l++){
            // single excitation from J-th (occupied)MO to B-th (virtual)MO
            int moJ = this->GetActiveOccIndex(*this->molecule, l);
            int moB = this->GetActiveVirIndex(*this->molecule, l);
            double value=0.0;
          
            // Fast algorith, but this is not easy to read. 
            // Slow algorithm is alos written below.
            for(int A=0; A<molecule->GetNumberAtoms(); A++){
               const Atom& atomA = *molecule->GetAtom(A);
               int firstAOIndexA = atomA.GetFirstAOIndex();
               int numberAOsA = atomA.GetValenceSize();

               for(int B=A; B<molecule->GetNumberAtoms(); B++){
                  const Atom& atomB = *molecule->GetAtom(B);
                  int firstAOIndexB = atomB.GetFirstAOIndex();
                  int numberAOsB = atomB.GetValenceSize();

                  double gamma = 0.0;
                  if(A!=B){
                     for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
                        for(int nu=mu; nu<firstAOIndexA+numberAOsA; nu++){
                           for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                              for(int sigma=lambda; sigma<firstAOIndexB+numberAOsB; sigma++){
                                 OrbitalType orbitalSigma = atomB.GetValence(sigma-firstAOIndexB);
                                 gamma = this->twoElecTwoCore[A]
                                                             [B]
                                                             [mu-firstAOIndexA]
                                                             [nu-firstAOIndexA]
                                                             [lambda-firstAOIndexB]
                                                             [sigma-firstAOIndexB];
   
                                 value += 2.0*gamma*fockMatrix[moA][mu]
                                                   *fockMatrix[moI][nu]
                                                   *fockMatrix[moJ][lambda]
                                                      *fockMatrix[moB][sigma];
                                 value += 2.0*gamma*fockMatrix[moA][lambda]
                                                   *fockMatrix[moI][sigma]
                                                   *fockMatrix[moJ][mu]
                                                   *fockMatrix[moB][nu];
                                 value -= gamma*fockMatrix[moA][mu]
                                               *fockMatrix[moB][nu]
                                               *fockMatrix[moI][lambda]
                                               *fockMatrix[moJ][sigma];
                                 value -= gamma*fockMatrix[moA][lambda]
                                               *fockMatrix[moB][sigma]
                                               *fockMatrix[moI][mu]
                                               *fockMatrix[moJ][nu];
                                 if(lambda != sigma){
                                    value += 2.0*gamma*fockMatrix[moA][mu]
                                                      *fockMatrix[moI][nu]
                                                      *fockMatrix[moJ][sigma]
                                                      *fockMatrix[moB][lambda];
                                    value += 2.0*gamma*fockMatrix[moA][sigma]
                                                      *fockMatrix[moI][lambda]
                                                      *fockMatrix[moJ][mu]
                                                      *fockMatrix[moB][nu];
                                    value -= gamma*fockMatrix[moA][mu]
                                                  *fockMatrix[moB][nu]
                                                  *fockMatrix[moI][sigma]
                                                  *fockMatrix[moJ][lambda];
                                    value -= gamma*fockMatrix[moA][sigma]
                                                  *fockMatrix[moB][lambda]
                                                  *fockMatrix[moI][mu]
                                                  *fockMatrix[moJ][nu];
                                 }
                                 if(mu != nu){
                                    value += 2.0*gamma*fockMatrix[moA][nu]
                                                      *fockMatrix[moI][mu]
                                                      *fockMatrix[moJ][lambda]
                                                      *fockMatrix[moB][sigma];
                                    value += 2.0*gamma*fockMatrix[moA][lambda]
                                                      *fockMatrix[moI][sigma]
                                                      *fockMatrix[moJ][nu]
                                                      *fockMatrix[moB][mu];
                                    value -= gamma*fockMatrix[moA][nu]
                                                  *fockMatrix[moB][mu]
                                                  *fockMatrix[moI][lambda]
                                                  *fockMatrix[moJ][sigma];
                                    value -= gamma*fockMatrix[moA][lambda]
                                                  *fockMatrix[moB][sigma]
                                                  *fockMatrix[moI][nu]
                                                  *fockMatrix[moJ][mu];
                                 }
                                 if(mu != nu && lambda != sigma){
                                    value += 2.0*gamma*fockMatrix[moA][nu]
                                                      *fockMatrix[moI][mu]
                                                      *fockMatrix[moJ][sigma]
                                                      *fockMatrix[moB][lambda];
                                    value += 2.0*gamma*fockMatrix[moA][sigma]
                                                      *fockMatrix[moI][lambda]
                                                      *fockMatrix[moJ][nu]
                                                      *fockMatrix[moB][mu];
                                    value -= gamma*fockMatrix[moA][nu]
                                                  *fockMatrix[moB][mu]
                                                  *fockMatrix[moI][sigma]
                                                  *fockMatrix[moJ][lambda];
                                    value -= gamma*fockMatrix[moA][sigma]
                                                  *fockMatrix[moB][lambda]
                                                  *fockMatrix[moI][nu]
                                                  *fockMatrix[moJ][mu];
                                 }
                              }
                           }
                        }
                     }
                  }
                  else{
                     for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
                        for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
                           for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                              for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
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
                                 value += 2.0*gamma*fockMatrix[moA][mu]
                                                   *fockMatrix[moI][nu]
                                                   *fockMatrix[moJ][lambda]
                                                   *fockMatrix[moB][sigma];
                                 value -= gamma*fockMatrix[moA][mu]
                                               *fockMatrix[moB][nu]
                                               *fockMatrix[moI][lambda]
                                               *fockMatrix[moJ][sigma];
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
                                                          this->molecule, this->fockMatrix, NULL)
                       -this->GetMolecularIntegralElement(moA, moB, moI, moJ, 
                                                          this->molecule, this->fockMatrix, NULL);
            // End of the slow algorith.
            */
            // Diagonal term
            if(k==l){
               value += this->energiesMO[moA] - this->energiesMO[moI];
            }
            matrixCIS[k][l] = value;
         }
      }
      catch(MolDSException ex){
         #pragma omp critical
         ompErrors << ex.what() << endl ;
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
   double ompEndTime = omp_get_wtime();
   this->OutputLog((boost::format("%s%lf%s\n%s") % this->messageOmpElapsedTimeCalcCISMarix.c_str()
                                                 % (ompEndTime - ompStartTime)
                                                 % this->messageUnitSec.c_str()
                                                 % this->messageDoneCalcCISMatrix.c_str()).str());
}

void Mndo::CheckZMatrixForce(const vector<int>& elecStates){
   // malloc or initialize Z matrix
   if(this->zMatrixForce == NULL){
      MallocerFreer::GetInstance()->Malloc<double>(&this->zMatrixForce, 
                                                   elecStates.size(),
                                                   this->molecule->GetTotalNumberAOs(), 
                                                   this->molecule->GetTotalNumberAOs());
      this->zMatrixForceElecStatesNum = elecStates.size();
   }
   else{
      MallocerFreer::GetInstance()->
      Initialize<double>(this->zMatrixForce,
                         elecStates.size(),
                         this->molecule->GetTotalNumberAOs(), 
                         this->molecule->GetTotalNumberAOs());
   }
}

void Mndo::CheckEtaMatrixForce(const vector<int>& elecStates){
   // malloc or initialize eta matrix
   if(this->etaMatrixForce == NULL){
      MallocerFreer::GetInstance()->Malloc<double>(&this->etaMatrixForce, 
                                                   elecStates.size(),
                                                   this->molecule->GetTotalNumberAOs(), 
                                                   this->molecule->GetTotalNumberAOs());
      this->etaMatrixForceElecStatesNum = elecStates.size();
   }
   else{
      MallocerFreer::GetInstance()->
      Initialize<double>(this->etaMatrixForce,
                         elecStates.size(),
                         this->molecule->GetTotalNumberAOs(), 
                         this->molecule->GetTotalNumberAOs());
   }
}

// see variable Q-vector in [PT_1996, PT_1997]
void Mndo::CalcActiveSetVariablesQ(vector<MoIndexPair>* nonRedundantQIndeces, 
                                   vector<MoIndexPair>* redundantQIndeces,
                                   int numberActiveOcc,
                                   int numberActiveVir) const{
   int numberAOs = this->molecule->GetTotalNumberAOs();
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   for(int moI=0; moI<numberOcc; moI++){
      bool isMoICIMO = numberOcc-numberActiveOcc<=moI ? true : false;
      for(int moJ=numberOcc; moJ<numberAOs; moJ++){
         bool isMoJCIMO = moJ<numberOcc+numberActiveVir ? true : false;
         MoIndexPair moIndexPair = {moI, moJ, isMoICIMO, isMoJCIMO};
         nonRedundantQIndeces->push_back(moIndexPair);
      }
   }
   for(int moI=numberOcc-numberActiveOcc; moI<numberOcc; moI++){
      bool isMoICIMO = true;
      for(int moJ=moI; moJ<numberOcc; moJ++){
         bool isMoJCIMO = true;
         MoIndexPair moIndexPair = {moI, moJ, isMoICIMO, isMoJCIMO};
         redundantQIndeces->push_back(moIndexPair);
      }
   }
   for(int moI=numberOcc; moI<numberOcc+numberActiveVir; moI++){
      bool isMoICIMO = true;
      for(int moJ=moI; moJ<numberOcc+numberActiveVir; moJ++){
         bool isMoJCIMO = true;
         MoIndexPair moIndexPair = {moI, moJ, isMoICIMO, isMoJCIMO};
         redundantQIndeces->push_back(moIndexPair);
      }
   }
}

void Mndo::MallocTempMatrixForZMatrix(double** delta,
                                      double** q,
                                      double*** gammaNRMinusKNR,
                                      double*** kRDag,
                                      double** y,
                                      double*** transposedFockMatrix,
                                      double*** xiOcc,
                                      double*** xiVir,
                                      int sizeQNR,
                                      int sizeQR) const{
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   int numberActiveMO = numberActiveOcc + numberActiveVir;
   int numberAOs = this->molecule->GetTotalNumberAOs();
   MallocerFreer::GetInstance()->Malloc<double>(delta, numberActiveMO);
   MallocerFreer::GetInstance()->Malloc<double>(q, sizeQNR+sizeQR);
   MallocerFreer::GetInstance()->Malloc<double>(gammaNRMinusKNR, sizeQNR, sizeQNR);
   MallocerFreer::GetInstance()->Malloc<double>(kRDag, sizeQNR, sizeQR);
   MallocerFreer::GetInstance()->Malloc<double>(y, sizeQNR);
   MallocerFreer::GetInstance()->Malloc<double>(transposedFockMatrix,
                                                numberAOs,
                                                numberAOs);
   MallocerFreer::GetInstance()->Malloc<double>(xiOcc, numberActiveOcc,numberAOs);
   MallocerFreer::GetInstance()->Malloc<double>(xiVir,numberActiveVir,numberAOs);
}

void Mndo::FreeTempMatrixForZMatrix(double** delta,
                                    double** q,
                                    double*** gammaNRMinusKNR,
                                    double*** kRDag,
                                    double** y,
                                    double*** transposedFockMatrix,
                                    double*** xiOcc,
                                    double*** xiVir,
                                    int sizeQNR,
                                    int sizeQR) const{
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   int numberActiveMO = numberActiveOcc + numberActiveVir;
   int numberAOs = this->molecule->GetTotalNumberAOs();
   MallocerFreer::GetInstance()->Free<double>(delta, numberActiveMO);
   MallocerFreer::GetInstance()->Free<double>(q, sizeQNR+sizeQR);
   MallocerFreer::GetInstance()->Free<double>(gammaNRMinusKNR, sizeQNR, sizeQNR);
   MallocerFreer::GetInstance()->Free<double>(kRDag, sizeQNR, sizeQR);
   MallocerFreer::GetInstance()->Free<double>(y, sizeQNR);
   MallocerFreer::GetInstance()->Free<double>(transposedFockMatrix, numberAOs, numberAOs);
   MallocerFreer::GetInstance()->Free<double>(xiOcc, numberActiveOcc, numberAOs);
   MallocerFreer::GetInstance()->Free<double>(xiVir, numberActiveVir, numberAOs);
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

// see (40) in [PT_1996]
double Mndo::GetGammaNRElement(int moI, int moJ, int moK, int moL) const{
   double value=0.0;
   if(moI==moK && moJ==moL){
      int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
      double nI = moI<numberOcc ? 2.0 : 0.0;
      double nJ = moJ<numberOcc ? 2.0 : 0.0;
      value = (this->energiesMO[moJ]-this->energiesMO[moI])/(nJ-nI);
   }
   return value;
}

// see (41) & (42) in [PT_1996]
double Mndo::GetGammaRElement(int moI, int moJ, int moK, int moL) const{
   double value=0.0;
   if(moI==moK && moJ==moL){
      value = moI==moJ ? 1.0 : this->energiesMO[moJ]-this->energiesMO[moI];
   }
   return value;
}

// see (43) in [PT_1996]
double Mndo::GetNNRElement(int moI, int moJ, int moK, int moL) const{
   double value=0.0;
   if(moI==moK && moJ==moL){
      int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
      double nI = moI<numberOcc ? 2.0 : 0.0;
      double nJ = moJ<numberOcc ? 2.0 : 0.0;
      value = (nJ-nI);
   }
   return value;
}

// see (44) in [PT_1996]
double Mndo::GetNRElement(int moI, int moJ, int moK, int moL) const{
   double value=0.0;
   if(moI==moK && moJ==moL){
      value = 1.0;
   }
   return value;
}

// see (44) in [PT_1996]
double Mndo::GetKNRElement(int moI, int moJ, int moK, int moL) const{
   double value=0.0;
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int nI = moI<numberOcc ? 2 : 0;
   int nJ = moJ<numberOcc ? 2 : 0;
   int nK = moK<numberOcc ? 2 : 0;
   int nL = moL<numberOcc ? 2 : 0;
   
   if(nI!=nJ && nK!=nL){
      value = this->GetAuxiliaryKNRKRElement(moI, moJ, moK, moL);
   }
   return 0.5*value;
}

// Dager of (45) in [PT_1996]. Note taht the (45) is real number.
double Mndo::GetKRDagerElement(int moI, int moJ, int moK, int moL) const{
   return this->GetKRElement(moK, moL, moI, moJ);
}

// see (45) in [PT_1996]
double Mndo::GetKRElement(int moI, int moJ, int moK, int moL) const{
   double value=0.0;
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int nI = moI<numberOcc ? 2 : 0;
   int nJ = moJ<numberOcc ? 2 : 0;
   int nK = moK<numberOcc ? 2 : 0;
   int nL = moL<numberOcc ? 2 : 0;

   if(nI==nJ && nK!=nL){
      value = this->GetAuxiliaryKNRKRElement(moI, moJ, moK, moL);
   }
   return 0.5*value;
}

double Mndo::GetAuxiliaryKNRKRElement(int moI, int moJ, int moK, int moL) const{
   double value = 0.0;

   // Fast algorith, but this is not easy to read. 
   // Slow algorithm is alos written below.
   for(int A=0; A<this->molecule->GetNumberAtoms(); A++){
      const Atom& atomA = *this->molecule->GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int numberAOsA = atomA.GetValenceSize();

      for(int B=A; B<this->molecule->GetNumberAtoms(); B++){
         const Atom& atomB = *this->molecule->GetAtom(B);
         int firstAOIndexB = atomB.GetFirstAOIndex();
         int numberAOsB = atomB.GetValenceSize();

         double gamma = 0.0;
         if(A!=B){
            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               for(int nu=mu; nu<firstAOIndexA+numberAOsA; nu++){
                  for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                     for(int sigma=lambda; sigma<firstAOIndexB+numberAOsB; sigma++){
                        OrbitalType orbitalSigma = atomB.GetValence(sigma-firstAOIndexB);
                        gamma = this->twoElecTwoCore[A]
                                                    [B]
                                                    [mu-firstAOIndexA]
                                                    [nu-firstAOIndexA]
                                                    [lambda-firstAOIndexB]
                                                    [sigma-firstAOIndexB];
                        // Order in moI, moJ, moK, then moL 
                        value += 4.0*gamma*this->fockMatrix[moI][mu]
                                          *this->fockMatrix[moJ][nu]
                                          *this->fockMatrix[moK][lambda]
                                          *this->fockMatrix[moL][sigma];
                        value += 4.0*gamma*this->fockMatrix[moI][lambda]
                                          *this->fockMatrix[moJ][sigma]
                                          *this->fockMatrix[moK][mu]
                                          *this->fockMatrix[moL][nu];
                        // Order in moI, moK, moJ, then moL 
                        value -= gamma*this->fockMatrix[moI][mu]
                                      *this->fockMatrix[moK][nu]
                                      *this->fockMatrix[moJ][lambda]
                                      *this->fockMatrix[moL][sigma];
                        value -= gamma*this->fockMatrix[moI][lambda]
                                      *this->fockMatrix[moK][sigma]
                                      *this->fockMatrix[moJ][mu]
                                      *this->fockMatrix[moL][nu];
                        // Order in moI, moL, moJ, then moK 
                        value -= gamma*this->fockMatrix[moI][mu]
                                      *this->fockMatrix[moL][nu]
                                      *this->fockMatrix[moJ][lambda]
                                      *this->fockMatrix[moK][sigma];
                        value -= gamma*this->fockMatrix[moI][lambda]
                                      *this->fockMatrix[moL][sigma]
                                      *this->fockMatrix[moJ][mu]
                                      *this->fockMatrix[moK][nu];
                        if(lambda != sigma){
                           // Order in moI, moJ, moK, then moL 
                           value += 4.0*gamma*this->fockMatrix[moI][mu]
                                             *this->fockMatrix[moJ][nu]
                                             *this->fockMatrix[moK][sigma]
                                             *this->fockMatrix[moL][lambda];
                           value += 4.0*gamma*this->fockMatrix[moI][sigma]
                                             *this->fockMatrix[moJ][lambda]
                                             *this->fockMatrix[moK][mu]
                                             *this->fockMatrix[moL][nu];
                           // Order in moI, moK, moJ, then moL 
                           value -= gamma*this->fockMatrix[moI][mu]
                                         *this->fockMatrix[moK][nu]
                                         *this->fockMatrix[moJ][sigma]
                                         *this->fockMatrix[moL][lambda];
                           value -= gamma*this->fockMatrix[moI][sigma]
                                         *this->fockMatrix[moK][lambda]
                                         *this->fockMatrix[moJ][mu]
                                         *this->fockMatrix[moL][nu];
                           // Order in moI, moL, moJ, then moK 
                           value -= gamma*this->fockMatrix[moI][mu]
                                         *this->fockMatrix[moL][nu]
                                         *this->fockMatrix[moJ][sigma]
                                         *this->fockMatrix[moK][lambda];
                           value -= gamma*this->fockMatrix[moI][sigma]
                                         *this->fockMatrix[moL][lambda]
                                         *this->fockMatrix[moJ][mu]
                                         *this->fockMatrix[moK][nu];
                        }
                        if(mu != nu){
                           // Order in moI, moJ, moK, then moL 
                           value += 4.0*gamma*this->fockMatrix[moI][nu]
                                             *this->fockMatrix[moJ][mu]
                                             *this->fockMatrix[moK][lambda]
                                             *this->fockMatrix[moL][sigma];
                           value += 4.0*gamma*this->fockMatrix[moI][lambda]
                                             *this->fockMatrix[moJ][sigma]
                                             *this->fockMatrix[moK][nu]
                                             *this->fockMatrix[moL][mu];
                           // Order in moI, moK, moJ, then moL 
                           value -= gamma*this->fockMatrix[moI][nu]
                                         *this->fockMatrix[moK][mu]
                                         *this->fockMatrix[moJ][lambda]
                                         *this->fockMatrix[moL][sigma];
                           value -= gamma*this->fockMatrix[moI][lambda]
                                         *this->fockMatrix[moK][sigma]
                                         *this->fockMatrix[moJ][nu]
                                         *this->fockMatrix[moL][mu];
                           // Order in moI, moL, moJ, then moK 
                           value -= gamma*this->fockMatrix[moI][nu]
                                         *this->fockMatrix[moL][mu]
                                         *this->fockMatrix[moJ][lambda]
                                         *this->fockMatrix[moK][sigma];
                           value -= gamma*this->fockMatrix[moI][lambda]
                                         *this->fockMatrix[moL][sigma]
                                         *this->fockMatrix[moJ][nu]
                                         *this->fockMatrix[moK][mu];
                        }
                        if(mu != nu && lambda != sigma){
                           // Order in moI, moJ, moK, then moL 
                           value += 4.0*gamma*this->fockMatrix[moI][nu]
                                             *this->fockMatrix[moJ][mu]
                                             *this->fockMatrix[moK][sigma]
                                             *this->fockMatrix[moL][lambda];
                           value += 4.0*gamma*this->fockMatrix[moI][sigma]
                                             *this->fockMatrix[moJ][lambda]
                                             *this->fockMatrix[moK][nu]
                                             *this->fockMatrix[moL][mu];
                           // Order in moI, moK, moJ, then moL 
                           value -= gamma*this->fockMatrix[moI][nu]
                                         *this->fockMatrix[moK][mu]
                                         *this->fockMatrix[moJ][sigma]
                                         *this->fockMatrix[moL][lambda];
                           value -= gamma*this->fockMatrix[moI][sigma]
                                         *this->fockMatrix[moK][lambda]
                                         *this->fockMatrix[moJ][nu]
                                         *this->fockMatrix[moL][mu];
                           // Order in moI, moL, moJ, then moK 
                           value -= gamma*this->fockMatrix[moI][nu]
                                         *this->fockMatrix[moL][mu]
                                         *this->fockMatrix[moJ][sigma]
                                         *this->fockMatrix[moK][lambda];
                           value -= gamma*this->fockMatrix[moI][sigma]
                                         *this->fockMatrix[moL][lambda]
                                         *this->fockMatrix[moJ][nu]
                                         *this->fockMatrix[moK][mu];
                        }
                     }
                  }
               }
            }
         }
         else{
            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
                  for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                     for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
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
                        // Order in moI, moJ, moK, then moL 
                        value += 4.0*gamma*this->fockMatrix[moI][mu]
                                          *this->fockMatrix[moJ][nu]
                                          *this->fockMatrix[moK][lambda]
                                          *this->fockMatrix[moL][sigma];
                        // Order in moI, moK, moJ, then moL 
                        value -= gamma*this->fockMatrix[moI][mu]
                                      *this->fockMatrix[moK][nu]
                                      *this->fockMatrix[moJ][lambda]
                                      *this->fockMatrix[moL][sigma];
                        // Order in moI, moL, moJ, then moK 
                        value -= gamma*this->fockMatrix[moI][mu]
                                      *this->fockMatrix[moL][nu]
                                      *this->fockMatrix[moJ][lambda]
                                      *this->fockMatrix[moK][sigma];
                     }  
                  }
               }
            }
         }
      }
   }
   // End of the fast algorith.
   
   /*
   // slow algorythm
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

// see (9) in [PT_1997]
void Mndo::CalcDeltaVector(double* delta, int exciteState) const{
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   int numberActiveMO = numberActiveOcc + numberActiveVir;
   MallocerFreer::GetInstance()->Initialize<double>(delta, numberActiveMO);
   stringstream ompErrors;
   #pragma omp parallel for schedule(auto)
   for(int r=0; r<numberActiveMO; r++){
      try{
         double value = 0.0;
         if(r<numberActiveOcc){
            // r is active occupied MO
            int rr=numberActiveOcc-(r+1);
            for(int a=0; a<numberActiveVir; a++){
               int slaterDeterminantIndex = this->GetSlaterDeterminantIndex(rr,a);
               value -= pow(this->matrixCIS[exciteState][slaterDeterminantIndex],2.0);
            }
         }
         else{
            // r is active virtual MO
            int rr=r-numberActiveOcc;
            for(int i=0; i<numberActiveOcc; i++){
               int slaterDeterminantIndex = this->GetSlaterDeterminantIndex(i,rr);
               value += pow(this->matrixCIS[exciteState][slaterDeterminantIndex],2.0);
            }
         }
         delta[r] = value;
      }
      catch(MolDSException ex){
         #pragma omp critical
         ompErrors << ex.what() << endl ;
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
}

// see (18) in [PT_1977]
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
      int numberAOsA = atomA.GetValenceSize();

      for(int B=A; B<molecule->GetNumberAtoms(); B++){
         const Atom& atomB = *molecule->GetAtom(B);
         int firstAOIndexB = atomB.GetFirstAOIndex();
         int numberAOsB = atomB.GetValenceSize();

         if(A!=B){
            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               for(int nu=mu; nu<firstAOIndexA+numberAOsA; nu++){
                  for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                     for(int sigma=lambda; sigma<firstAOIndexB+numberAOsB; sigma++){
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
            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
                  for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                     for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
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

// see (20) - (23) in [PT_1997]
void Mndo::CalcQVector(double* q, 
                       double const* delta, 
                       double const* const* xiOcc,
                       double const* const* xiVir,
                       double const* const* eta,
                       const vector<MoIndexPair>& nonRedundantQIndeces,
                       const vector<MoIndexPair>& redundantQIndeces) const{
   MallocerFreer::GetInstance()->Initialize<double>(
                                 q,
                                 nonRedundantQIndeces.size()+redundantQIndeces.size());

   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   stringstream ompErrors;
   #pragma omp parallel for schedule(auto)
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      try{
         int moI = nonRedundantQIndeces[i].moI;
         int moJ = nonRedundantQIndeces[i].moJ;
         bool isMoICIMO = nonRedundantQIndeces[i].isMoICIMO;
         bool isMoJCIMO = nonRedundantQIndeces[i].isMoJCIMO;
         if(!isMoICIMO && isMoJCIMO){
            q[i] = this->GetSmallQElement(moI, moJ, xiOcc, xiVir, eta);
         }
         else if(isMoICIMO && !isMoJCIMO){
            q[i] = -1.0*this->GetSmallQElement(moJ, moI, xiOcc, xiVir, eta);
         }
         else if(isMoICIMO && isMoJCIMO){
            q[i] = this->GetSmallQElement(moI, moJ, xiOcc, xiVir, eta)
                  -this->GetSmallQElement(moJ, moI, xiOcc, xiVir, eta);
         }
         else{
            q[i] = 0.0;
         }
      }
      catch(MolDSException ex){
         #pragma omp critical
         ompErrors << ex.what() << endl ;
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
   #pragma omp parallel for schedule(auto)
   for(int i=0; i<redundantQIndeces.size(); i++){
      try{
         int r = nonRedundantQIndeces.size() + i;
         int moI = redundantQIndeces[i].moI;
         int moJ = redundantQIndeces[i].moJ;
         if(moI == moJ){
            int rr = moI - (numberOcc-numberActiveOcc);
            q[r] = delta[rr];
         }
         else{
            q[r] = this->GetSmallQElement(moI, moJ, xiOcc, xiVir, eta)
                  -this->GetSmallQElement(moJ, moI, xiOcc, xiVir, eta);
         }
      }
      catch(MolDSException ex){
         #pragma omp critical
         ompErrors << ex.what() << endl ;
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
   /* 
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      this->OutputLog((boost::format("q[%d] = %e\n") % i % q[i]).str());
   }
   for(int i=0; i<redundantQIndeces.size(); i++){
      int r = nonRedundantQIndeces.size() + i;
      this->OutputLog((boost::format("q[%d] = %e\n") % r % q[r]).str());
   }
   */
}

// see (40) and (45) in [PT_1996].
// This method calculates "\Gamma_{NR} - K_{NR}" to solve (54) in [PT_1966]
// Note taht K_{NR} is not calculated.
void Mndo::CalcGammaNRMinusKNRMatrix(double** gammaNRMinusKNR, const vector<MoIndexPair>& nonRedundantQIndeces) const{
   stringstream ompErrors;
   #pragma omp parallel for schedule(auto)
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      try{
         int moI = nonRedundantQIndeces[i].moI;
         int moJ = nonRedundantQIndeces[i].moJ;
         for(int j=i; j<nonRedundantQIndeces.size(); j++){
            int moK = nonRedundantQIndeces[j].moI;
            int moL = nonRedundantQIndeces[j].moJ;
            gammaNRMinusKNR[i][j] = this->GetGammaNRElement(moI, moJ, moK, moL)
                                   -this->GetKNRElement(moI, moJ, moK, moL);
         }
      }
      catch(MolDSException ex){
         #pragma omp critical
         ompErrors << ex.what() << endl ;
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
}

// see (41), (42), and (46) in [PT_1996].
// This method calculates "K_{R}^{\dagger} * Gamma_{R}" matrix, see (41), (42), and (46) to solve (54) in [PT_1996]
// Note taht K_{R}^{\dager} is not calculated.
void Mndo::CalcKRDagerGammaRInvMatrix(double** kRDagerGammaRInv, 
                                      const vector<MoIndexPair>& nonRedundantQIndeces,
                                      const vector<MoIndexPair>& redundantQIndeces) const{
   stringstream ompErrors;
   #pragma omp parallel for schedule(auto)
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      try{
         int moI = nonRedundantQIndeces[i].moI;
         int moJ = nonRedundantQIndeces[i].moJ;
         for(int j=0; j<redundantQIndeces.size(); j++){
            int moK = redundantQIndeces[j].moI;
            int moL = redundantQIndeces[j].moJ;
            kRDagerGammaRInv[i][j] = this->GetKRDagerElement(moI, moJ, moK, moL)
                                    /this->GetGammaRElement(moK, moL, moK, moL);
         }
      }
      catch(MolDSException ex){
         #pragma omp critical
         ompErrors << ex.what() << endl ;
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
}

// right hand side of (54) in [PT_1996]      
void Mndo::CalcAuxiliaryVector(double* y, 
                               double const* q, 
                               double const* const* kRDagerGammaRInv, 
                               const vector<MoIndexPair>& nonRedundantQIndeces, 
                               const vector<MoIndexPair>& redundantQIndeces) const{
   MallocerFreer::GetInstance()->Initialize<double>(
                                 y,
                                 nonRedundantQIndeces.size());
   stringstream ompErrors;
   #pragma omp parallel for schedule(auto)
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      try{
         int moI = nonRedundantQIndeces[i].moI;
         int moJ = nonRedundantQIndeces[i].moJ;
         y[i] += q[i]/this->GetNNRElement(moI, moJ, moI, moJ);
         for(int j=0; j<redundantQIndeces.size(); j++){
            int k = nonRedundantQIndeces.size() + j; 
            y[i] += kRDagerGammaRInv[i][j]*q[k];
         }
      }
      catch(MolDSException ex){
         #pragma omp critical
         ompErrors << ex.what() << endl ;
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
}

void Mndo::TransposeFockMatrixMatrix(double** transposedFockMatrix) const{
   for(int i=0; i<this->molecule->GetTotalNumberAOs(); i++){
      for(int j=0; j<this->molecule->GetTotalNumberAOs(); j++){
         transposedFockMatrix[j][i] = this->fockMatrix[i][j];
      }
   }
}

// each element (mu, nu) of z matrix.
// see (57) in [PT_1996]
double Mndo::GetZMatrixForceElement(double const* y,
                                    double const* q,
                                    double const* const* transposedFockMatrix,
                                    const vector<MoIndexPair>& nonRedundantQIndeces,
                                    const vector<MoIndexPair>& redundantQIndeces,
                                    int mu,
                                    int nu) const{
   double value=0.0;
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      int moI = nonRedundantQIndeces[i].moI;
      int moJ = nonRedundantQIndeces[i].moJ;
      value += y[i]
              *transposedFockMatrix[mu][moI]
              *transposedFockMatrix[nu][moJ];
   }
   for(int i=0; i<redundantQIndeces.size(); i++){
      int j = nonRedundantQIndeces.size() + i;
      int moI = redundantQIndeces[i].moI;
      int moJ = redundantQIndeces[i].moJ;
      value += (q[j]/this->GetGammaRElement(moI, moJ, moI, moJ))
              *transposedFockMatrix[mu][moI]
              *transposedFockMatrix[nu][moJ];
   }
   return value;
}

void Mndo::CalcXiMatrices(double** xiOcc, 
                          double** xiVir, 
                          int exciteState, 
                          double const* const* transposedFockMatrix) const{
   int numberAOs = this->molecule->GetTotalNumberAOs();
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   MallocerFreer::GetInstance()->Initialize<double>(
                                 xiOcc, numberActiveOcc, numberAOs);
   MallocerFreer::GetInstance()->Initialize<double>(
                                 xiVir, numberActiveVir, numberAOs);
   stringstream ompErrors;
   // xiOcc
   #pragma omp parallel for schedule(auto)
   for(int p=0; p<numberActiveOcc; p++){
      try{
         for(int mu=0; mu<numberAOs; mu++){
            for(int a=0; a<numberActiveVir; a++){
               int moA = numberOcc + a;
               int slaterDeterminantIndex = this->GetSlaterDeterminantIndex(p,a);
               xiOcc[p][mu] += this->matrixCIS[exciteState][slaterDeterminantIndex]
                              *transposedFockMatrix[mu][moA];
            }
         }
      }
      catch(MolDSException ex){
         #pragma omp critical
         ompErrors << ex.what() << endl ;
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
   // xiVir
   #pragma omp parallel for schedule(auto)
   for(int p=0; p<numberActiveVir; p++){
      try{
         for(int mu=0; mu<numberAOs; mu++){
            for(int i=0; i<numberActiveOcc; i++){
               int moI = numberOcc - (i+1);
               int slaterDeterminantIndex = this->GetSlaterDeterminantIndex(i,p);
               xiVir[p][mu] += this->matrixCIS[exciteState][slaterDeterminantIndex]
                              *transposedFockMatrix[mu][moI];
            }
         }
      }
      catch(MolDSException ex){
         #pragma omp critical
         ompErrors << ex.what() << endl ;
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
}

void Mndo::MallocTempMatricesEachThreadCalcHessianSCF(double***** diatomicOverlapFirstDerivs,
                                                      double****** diatomicOverlapSecondDerivs,
                                                      double******* diatomicTwoElecTwoCoreFirstDerivs,
                                                      double******** diatomicTwoElecTwoCoreSecondDerivs) const{
   MallocerFreer::GetInstance()->Malloc<double>(diatomicOverlapFirstDerivs,
                                                this->molecule->GetNumberAtoms(),
                                                this->molecule->GetTotalNumberAOs(),
                                                this->molecule->GetTotalNumberAOs(),
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicOverlapSecondDerivs,
                                                this->molecule->GetNumberAtoms(),
                                                this->molecule->GetTotalNumberAOs(),
                                                this->molecule->GetTotalNumberAOs(),
                                                CartesianType_end,
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicTwoElecTwoCoreFirstDerivs,
                                                this->molecule->GetNumberAtoms(),
                                                dxy,
                                                dxy,
                                                dxy,
                                                dxy,
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicTwoElecTwoCoreSecondDerivs,
                                                this->molecule->GetNumberAtoms(),
                                                dxy,
                                                dxy,
                                                dxy,
                                                dxy,
                                                CartesianType_end,
                                                CartesianType_end);
}

void Mndo::FreeTempMatricesEachThreadCalcHessianSCF(double***** diatomicOverlapFirstDerivs,
                                                    double****** diatomicOverlapSecondDerivs,
                                                    double******* diatomicTwoElecTwoCoreFirstDerivs,
                                                    double******** diatomicTwoElecTwoCoreSecondDerivs) const{
   MallocerFreer::GetInstance()->Free<double>(diatomicOverlapFirstDerivs,
                                              this->molecule->GetNumberAtoms(),
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs(),
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicOverlapSecondDerivs,
                                              this->molecule->GetNumberAtoms(),
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs(),
                                              CartesianType_end,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicTwoElecTwoCoreFirstDerivs,
                                              this->molecule->GetNumberAtoms(),
                                              dxy,
                                              dxy,
                                              dxy,
                                              dxy,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicTwoElecTwoCoreSecondDerivs,
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
                                         int atomAIndex,
                                         int atomCIndex,
                                         CartesianType axisA1,
                                         CartesianType axisA2,
                                         double const* const* orbitalElectronPopulation,
                                         double const* const* const* const* const* const* diatomicTwoElecTwoCoreSecondDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomC = *this->molecule->GetAtom(atomCIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   double value = orbitalElectronPopulation[mu]
                                           [nu]
                 *diatomicTwoElecTwoCoreSecondDerivs[mu-firstAOIndexA]
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
                                         int atomAIndex,
                                         int atomBIndex,
                                         int atomCIndex,
                                         CartesianType axisA,
                                         CartesianType axisB,
                                         double const* const* const* const* orbitalElectronPopulationFirstDerivs,
                                         double const* const* const* const* const* diatomicTwoElecTwoCoreFirstDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomC = *this->molecule->GetAtom(atomCIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   double value = orbitalElectronPopulationFirstDerivs[mu]
                                                      [nu]
                                                      [atomBIndex]
                                                      [axisB]
                 *diatomicTwoElecTwoCoreFirstDerivs[mu-firstAOIndexA]
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
                                         int atomAIndex,
                                         int atomCIndex,
                                         CartesianType axisA1,
                                         CartesianType axisA2,
                                         double const* const* orbitalElectronPopulation,
                                         double const* const* const* const* const* const* diatomicTwoElecTwoCoreSecondDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomC = *this->molecule->GetAtom(atomCIndex);
   int firstAOIndexC = atomC.GetFirstAOIndex();
   double value = orbitalElectronPopulation[lambda]
                                           [sigma]
                 *diatomicTwoElecTwoCoreSecondDerivs[s]
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
                                         int atomAIndex,
                                         int atomBIndex,
                                         int atomCIndex,
                                         CartesianType axisA,
                                         CartesianType axisB,
                                         double const* const* const* const* orbitalElectronPopulationFirstDerivs,
                                         double const* const* const* const* const* diatomicTwoElecTwoCoreFirstDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomC = *this->molecule->GetAtom(atomCIndex);
   int firstAOIndexC = atomC.GetFirstAOIndex();
   double value = orbitalElectronPopulationFirstDerivs[lambda]
                                                      [sigma]
                                                      [atomBIndex]
                                                      [axisB]
                 *diatomicTwoElecTwoCoreFirstDerivs[s]
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
                                         int atomAIndex,
                                         int atomCIndex,
                                         CartesianType axisA1,
                                         CartesianType axisA2,
                                         double const* const* orbitalElectronPopulation,
                                         double const* const* const* const* diatomicOverlapSecondDerivs) const{
   const Atom& atomA = *molecule->GetAtom(atomAIndex);
   const Atom& atomC = *molecule->GetAtom(atomCIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexC = atomC.GetFirstAOIndex();
   double bondParameterA = atomA.GetBondingParameter(this->theory, atomA.GetValence(mu-firstAOIndexA));
   double bondParameterC = atomC.GetBondingParameter(this->theory, atomC.GetValence(lambda-firstAOIndexC));
   double sumBondParameters = bondParameterA+bondParameterC;
   double value = orbitalElectronPopulation[mu][lambda]
                 *sumBondParameters
                 *diatomicOverlapSecondDerivs[mu-firstAOIndexA]
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
                                         int atomAIndex,
                                         int atomBIndex,
                                         int atomCIndex,
                                         CartesianType axisA,
                                         CartesianType axisB,
                                         double const* const* const* const* orbitalElectronPopulationFirstDerivs,
                                         double const* const* const* diatomicOverlapFirstDerivs) const{
   const Atom& atomA = *molecule->GetAtom(atomAIndex);
   const Atom& atomC = *molecule->GetAtom(atomCIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexC = atomC.GetFirstAOIndex();
   double bondParameterA = atomA.GetBondingParameter(this->theory, atomA.GetValence(mu-firstAOIndexA));
   double bondParameterC = atomC.GetBondingParameter(this->theory, atomC.GetValence(lambda-firstAOIndexC));
   double sumBondParameters = bondParameterA+bondParameterC;
   double value = orbitalElectronPopulationFirstDerivs[mu]
                                                      [lambda]
                                                      [atomBIndex]
                                                      [axisB]
                 *sumBondParameters
                 *diatomicOverlapFirstDerivs[mu-firstAOIndexA]
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
                                         int atomAIndex,
                                         int atomCIndex,
                                         CartesianType axisA1,
                                         CartesianType axisA2,
                                         double const* const* orbitalElectronPopulation,
                                         double const* const* const* const* const* const* diatomicTwoElecTwoCoreSecondDerivs) const{
   const Atom& atomA = *molecule->GetAtom(atomAIndex);
   const Atom& atomC = *molecule->GetAtom(atomCIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexC = atomC.GetFirstAOIndex();
   double temp1 = orbitalElectronPopulation[mu][nu]*orbitalElectronPopulation[lambda][sigma];
   double temp2 = orbitalElectronPopulation[mu][lambda]*orbitalElectronPopulation[nu][sigma];
   double value = (temp1 - 0.5*temp2)
                 *diatomicTwoElecTwoCoreSecondDerivs[mu-firstAOIndexA]
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
                                         int atomAIndex,
                                         int atomBIndex,
                                         int atomCIndex,
                                         CartesianType axisA,
                                         CartesianType axisB,
                                         double const* const* orbitalElectronPopulation,
                                         double const* const* const* const* orbitalElectronPopulationFirstDerivs,
                                         double const* const* const* const* const* diatomicTwoElecTwoCoreFirstDerivs) const{
   const Atom& atomA = *molecule->GetAtom(atomAIndex);
   const Atom& atomC = *molecule->GetAtom(atomCIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexC = atomC.GetFirstAOIndex();
   double temp1 = orbitalElectronPopulationFirstDerivs[mu][nu]       [atomBIndex][axisB]
                 *orbitalElectronPopulation           [lambda][sigma];
   double temp2 = orbitalElectronPopulation           [mu][nu]
                 *orbitalElectronPopulationFirstDerivs[lambda][sigma][atomBIndex][axisB];
   double temp3 = orbitalElectronPopulationFirstDerivs[mu][lambda]   [atomBIndex][axisB]
                 *orbitalElectronPopulation           [nu][sigma];
   double temp4 = orbitalElectronPopulation           [mu][lambda]
                 *orbitalElectronPopulationFirstDerivs[nu][sigma]    [atomBIndex][axisB];
   double value = (0.5*(temp1 + temp2) - 0.25*(temp3 + temp4))
                 *diatomicTwoElecTwoCoreFirstDerivs[mu-firstAOIndexA]
                                                   [nu-firstAOIndexA]
                                                   [lambda-firstAOIndexC]
                                                   [sigma-firstAOIndexC]
                                                   [axisA];
   return value;
}

// Return hessian element. 
// The Second derivative are related to axisA1 and axisA2.
// These axisA1 and axisA2 are the Cartesian coordinates of atomA labeled with atomAIndex.
double Mndo::GetHessianElementSameAtomsSCF(int atomAIndex, 
                                           CartesianType axisA1,
                                           CartesianType axisA2,
                                           bool isMassWeighted,
                                           double const* const* orbitalElectronPopulation,
                                           double const* const* const* const* orbitalElectronPopulationFirstDerivs,
                                           double const* const* const* const* diatomicOverlapFirstDerivs,
                                           double const* const* const* const* const* diatomicOverlapSecondDerivs,
                                           double const* const* const* const* const* const* diatomicTwoElecTwoCoreFirstDerivs,
                                           double const* const* const* const* const* const* const* diatomicTwoElecTwoCoreSecondDerivs) const{
   double value=0.0;
   int atomBIndex = atomAIndex;
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int numberAOsA = atomA.GetValenceSize();
   for(int atomCIndex=0; atomCIndex<this->molecule->GetNumberAtoms(); atomCIndex++){
      if(atomAIndex != atomCIndex){
         const Atom& atomC = *this->molecule->GetAtom(atomCIndex);
         int firstAOIndexC = atomC.GetFirstAOIndex();
         int numberAOsC = atomC.GetValenceSize();

         // second derivative of electronic part
         for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
            for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
               value -= this->GetAuxiliaryHessianElement1(mu, 
                                                          nu, 
                                                          atomAIndex,
                                                          atomCIndex,
                                                          static_cast<CartesianType>(axisA1), 
                                                          static_cast<CartesianType>(axisA2), 
                                                          orbitalElectronPopulation,
                                                          diatomicTwoElecTwoCoreSecondDerivs[atomCIndex]);
               value -= this->GetAuxiliaryHessianElement2(mu, 
                                                          nu, 
                                                          atomAIndex,
                                                          atomBIndex,
                                                          atomCIndex,
                                                          static_cast<CartesianType>(axisA1), 
                                                          static_cast<CartesianType>(axisA2), 
                                                          orbitalElectronPopulationFirstDerivs,
                                                          diatomicTwoElecTwoCoreFirstDerivs[atomCIndex]);
            }
         }
         for(int lambda=firstAOIndexC; lambda<firstAOIndexC+numberAOsC; lambda++){
            for(int sigma=firstAOIndexC; sigma<firstAOIndexC+numberAOsC; sigma++){
               value -= this->GetAuxiliaryHessianElement3(lambda, 
                                                          sigma, 
                                                          atomAIndex,
                                                          atomCIndex,
                                                          static_cast<CartesianType>(axisA1), 
                                                          static_cast<CartesianType>(axisA2), 
                                                          orbitalElectronPopulation,
                                                          diatomicTwoElecTwoCoreSecondDerivs[atomCIndex]);
               value -= this->GetAuxiliaryHessianElement4(lambda, 
                                                          sigma, 
                                                          atomAIndex,
                                                          atomBIndex,
                                                          atomCIndex,
                                                          static_cast<CartesianType>(axisA1), 
                                                          static_cast<CartesianType>(axisA2), 
                                                          orbitalElectronPopulationFirstDerivs,
                                                          diatomicTwoElecTwoCoreFirstDerivs[atomCIndex]);
            }
         }
         for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
            for(int lambda=firstAOIndexC; lambda<firstAOIndexC+numberAOsC; lambda++){
               value += this->GetAuxiliaryHessianElement5(mu, 
                                                          lambda, 
                                                          atomAIndex,
                                                          atomCIndex,
                                                          static_cast<CartesianType>(axisA1), 
                                                          static_cast<CartesianType>(axisA2), 
                                                          orbitalElectronPopulation,
                                                          diatomicOverlapSecondDerivs[atomCIndex]);
               value += this->GetAuxiliaryHessianElement6(mu, 
                                                          lambda, 
                                                          atomAIndex,
                                                          atomBIndex,
                                                          atomCIndex,
                                                          static_cast<CartesianType>(axisA1), 
                                                          static_cast<CartesianType>(axisA2), 
                                                          orbitalElectronPopulationFirstDerivs,
                                                          diatomicOverlapFirstDerivs[atomCIndex]);
            }
         }
         for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
            for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
               for(int lambda=firstAOIndexC; lambda<firstAOIndexC+numberAOsC; lambda++){
                  for(int sigma=firstAOIndexC; sigma<firstAOIndexC+numberAOsC; sigma++){
                     value += this->GetAuxiliaryHessianElement7(mu, 
                                                                nu, 
                                                                lambda, 
                                                                sigma, 
                                                                atomAIndex,
                                                                atomCIndex,
                                                                static_cast<CartesianType>(axisA1), 
                                                                static_cast<CartesianType>(axisA2), 
                                                                orbitalElectronPopulation,
                                                                diatomicTwoElecTwoCoreSecondDerivs[atomCIndex]);
                     value += this->GetAuxiliaryHessianElement8(mu, 
                                                                nu, 
                                                                lambda, 
                                                                sigma, 
                                                                atomAIndex,
                                                                atomBIndex,
                                                                atomCIndex,
                                                                static_cast<CartesianType>(axisA1), 
                                                                static_cast<CartesianType>(axisA2), 
                                                                orbitalElectronPopulation,
                                                                orbitalElectronPopulationFirstDerivs,
                                                                diatomicTwoElecTwoCoreFirstDerivs[atomCIndex]);
                  }
               }
            }
         }

         // second derivatives of the nuclear repulsions
         value += this->GetDiatomCoreRepulsionSecondDerivative(atomAIndex, 
                                                               atomCIndex, 
                                                               static_cast<CartesianType>(axisA1), 
                                                               static_cast<CartesianType>(axisA2));
         // second derivatives of the van der waals corrections
         if(Parameters::GetInstance()->RequiresVdWSCF()){
            value += this->GetDiatomVdWCorrectionSecondDerivative(atomAIndex, 
                                                                  atomCIndex, 
                                                                  static_cast<CartesianType>(axisA1), 
                                                                  static_cast<CartesianType>(axisA2));
         }
      }
   }

   if(isMassWeighted){
      value /= atomA.GetCoreMass();
   }
   return value;
}

// Return hessian element. 
// The Second derivative are related to axisA and axisB.
// These axisA and axisB are the Cartesian coordinates of atomA and atomB, respectively.
double Mndo::GetHessianElementDifferentAtomsSCF(int atomAIndex, 
                                                int atomBIndex,
                                                CartesianType axisA,
                                                CartesianType axisB,
                                                bool isMassWeighted,
                                                double const* const* orbitalElectronPopulation,
                                                double const* const* const* const* orbitalElectronPopulationFirstDerivs,
                                                double const* const* const* const* diatomicOverlapFirstDerivs,
                                                double const* const* const* const* const* diatomicOverlapSecondDerivs,
                                                double const* const* const* const* const* const* diatomicTwoElecTwoCoreFirstDerivs,
                                                double const* const* const* const* const* const* const* diatomicTwoElecTwoCoreSecondDerivs) const{
   double value=0.0;
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int numberAOsA = atomA.GetValenceSize();
   int numberAOsB = atomB.GetValenceSize();

   // second derivative of electronic part
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
         value += this->GetAuxiliaryHessianElement1(mu, 
                                                    nu, 
                                                    atomAIndex,
                                                    atomBIndex,
                                                    static_cast<CartesianType>(axisA), 
                                                    static_cast<CartesianType>(axisB), 
                                                    orbitalElectronPopulation,
                                                    diatomicTwoElecTwoCoreSecondDerivs[atomBIndex]);
      }
   }
   for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
      for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
         value += this->GetAuxiliaryHessianElement3(lambda, 
                                                    sigma, 
                                                    atomAIndex,
                                                    atomBIndex,
                                                    static_cast<CartesianType>(axisA), 
                                                    static_cast<CartesianType>(axisB), 
                                                    orbitalElectronPopulation,
                                                    diatomicTwoElecTwoCoreSecondDerivs[atomBIndex]);
      }
   }
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
         value -= this->GetAuxiliaryHessianElement5(mu, 
                                                    lambda, 
                                                    atomAIndex,
                                                    atomBIndex,
                                                    static_cast<CartesianType>(axisA), 
                                                    static_cast<CartesianType>(axisB), 
                                                    orbitalElectronPopulation,
                                                    diatomicOverlapSecondDerivs[atomBIndex]);
      }
   }
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
         for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
            for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
               value += this->GetAuxiliaryHessianElement7(mu, 
                                                          nu, 
                                                          lambda, 
                                                          sigma, 
                                                          atomAIndex,
                                                          atomBIndex,
                                                          static_cast<CartesianType>(axisA), 
                                                          static_cast<CartesianType>(axisB), 
                                                          orbitalElectronPopulation,
                                                          diatomicTwoElecTwoCoreSecondDerivs[atomBIndex]);
            }
         }
      }
   }

   for(int atomCIndex=0; atomCIndex<this->molecule->GetNumberAtoms(); atomCIndex++){
      if(atomAIndex != atomCIndex){
         const Atom& atomC = *this->molecule->GetAtom(atomCIndex);
         int firstAOIndexC = atomC.GetFirstAOIndex();
         int numberAOsC = atomC.GetValenceSize();

         // second derivative of electronic part
         for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
            for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
               value -= this->GetAuxiliaryHessianElement2(mu, 
                                                          nu, 
                                                          atomAIndex,
                                                          atomBIndex,
                                                          atomCIndex,
                                                          static_cast<CartesianType>(axisA), 
                                                          static_cast<CartesianType>(axisB), 
                                                          orbitalElectronPopulationFirstDerivs,
                                                          diatomicTwoElecTwoCoreFirstDerivs[atomCIndex]);
            }
         }
         for(int lambda=firstAOIndexC; lambda<firstAOIndexC+numberAOsC; lambda++){
            for(int sigma=firstAOIndexC; sigma<firstAOIndexC+numberAOsC; sigma++){
               value -= this->GetAuxiliaryHessianElement4(lambda, 
                                                          sigma, 
                                                          atomAIndex,
                                                          atomBIndex,
                                                          atomCIndex,
                                                          static_cast<CartesianType>(axisA), 
                                                          static_cast<CartesianType>(axisB), 
                                                          orbitalElectronPopulationFirstDerivs,
                                                          diatomicTwoElecTwoCoreFirstDerivs[atomCIndex]);
            }
         }
         for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
            for(int lambda=firstAOIndexC; lambda<firstAOIndexC+numberAOsC; lambda++){
               value += this->GetAuxiliaryHessianElement6(mu, 
                                                          lambda, 
                                                          atomAIndex,
                                                          atomBIndex,
                                                          atomCIndex,
                                                          static_cast<CartesianType>(axisA), 
                                                          static_cast<CartesianType>(axisB), 
                                                          orbitalElectronPopulationFirstDerivs,
                                                          diatomicOverlapFirstDerivs[atomCIndex]);
            }
         }
         for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
            for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
               for(int lambda=firstAOIndexC; lambda<firstAOIndexC+numberAOsC; lambda++){
                  for(int sigma=firstAOIndexC; sigma<firstAOIndexC+numberAOsC; sigma++){
                     value += this->GetAuxiliaryHessianElement8(mu, 
                                                                nu, 
                                                                lambda, 
                                                                sigma, 
                                                                atomAIndex,
                                                                atomBIndex,
                                                                atomCIndex,
                                                                static_cast<CartesianType>(axisA), 
                                                                static_cast<CartesianType>(axisB), 
                                                                orbitalElectronPopulation,
                                                                orbitalElectronPopulationFirstDerivs,
                                                                diatomicTwoElecTwoCoreFirstDerivs[atomCIndex]);
                  }
               }
            }
         }
      }
   }

   // second derivatives of the nuclear repulsions
   value -= this->GetDiatomCoreRepulsionSecondDerivative(atomAIndex, 
                                                         atomBIndex, 
                                                         static_cast<CartesianType>(axisA), 
                                                         static_cast<CartesianType>(axisB));
   // second derivatives of the van der waals corrections
   if(Parameters::GetInstance()->RequiresVdWSCF()){
      value -= this->GetDiatomVdWCorrectionSecondDerivative(atomAIndex, 
                                                            atomBIndex, 
                                                            static_cast<CartesianType>(axisA), 
                                                            static_cast<CartesianType>(axisB));
   }

   if(isMassWeighted){
      value /= sqrt(atomA.GetCoreMass()*atomB.GetCoreMass());
   }
   return value;
}
void Mndo::CalcHessianSCF(double** hessianSCF, bool isMassWeighted) const{
   int totalNumberAOs = this->molecule->GetTotalNumberAOs();
   double**** orbitalElectronPopulationFirstDerivs = NULL;

   try{
      MallocerFreer::GetInstance()->Malloc<double>(&orbitalElectronPopulationFirstDerivs, 
                                                   totalNumberAOs,
                                                   totalNumberAOs,
                                                   this->molecule->GetNumberAtoms(),
                                                   CartesianType_end);
      this->CalcOrbitalElectronPopulationFirstDerivatives(orbitalElectronPopulationFirstDerivs);

//#pragma omp parallel
//{
      double**** diatomicOverlapFirstDerivs = NULL;
      double***** diatomicOverlapSecondDerivs = NULL;
      double****** diatomicTwoElecTwoCoreFirstDerivs = NULL;
      double******* diatomicTwoElecTwoCoreSecondDerivs = NULL;
      this->MallocTempMatricesEachThreadCalcHessianSCF(&diatomicOverlapFirstDerivs,
                                                       &diatomicOverlapSecondDerivs,
                                                       &diatomicTwoElecTwoCoreFirstDerivs, 
                                                       &diatomicTwoElecTwoCoreSecondDerivs);

//#pragma omp for schedule(auto)                                                 
      for(int atomAIndex=0; atomAIndex<this->molecule->GetNumberAtoms(); atomAIndex++){
         const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
         int firstAOIndexA = atomA.GetFirstAOIndex();
         int numberAOsA = atomA.GetValenceSize();
         for(int axisA = XAxis; axisA<CartesianType_end; axisA++){

            // calculation of derivatives of the overlaps and two electron integrals
            for(int atomBIndex=0; atomBIndex<this->molecule->GetNumberAtoms(); atomBIndex++){
               if(atomAIndex != atomBIndex){
                  this->CalcDiatomicOverlapFirstDerivatives(diatomicOverlapFirstDerivs[atomBIndex], 
                                                            atomAIndex, 
                                                            atomBIndex);
                  this->CalcDiatomicOverlapSecondDerivatives(diatomicOverlapSecondDerivs[atomBIndex], 
                                                             atomAIndex, 
                                                             atomBIndex);
                  this->CalcDiatomicTwoElecTwoCoreFirstDerivatives(diatomicTwoElecTwoCoreFirstDerivs[atomBIndex], 
                                                                   atomAIndex, 
                                                                   atomBIndex);
                  this->CalcDiatomicTwoElecTwoCoreSecondDerivatives(diatomicTwoElecTwoCoreSecondDerivs[atomBIndex], 
                                                                    atomAIndex, 
                                                                    atomBIndex);
               }
            }

            int k = atomAIndex*CartesianType_end + axisA; // hessian index, i.e. hessian[k][l]
            // hessian element (atomA == atomB)
            for(int axisA2 = axisA; axisA2<CartesianType_end; axisA2++){
               int l = atomAIndex*CartesianType_end + axisA2; // hessian index, i.e. hessian[k][l]
               hessianSCF[l][k] = 
               hessianSCF[k][l] = this->GetHessianElementSameAtomsSCF(atomAIndex, 
                                                                      static_cast<CartesianType>(axisA), 
                                                                      static_cast<CartesianType>(axisA2), 
                                                                      isMassWeighted,
                                                                      orbitalElectronPopulation,
                                                                      orbitalElectronPopulationFirstDerivs,
                                                                      diatomicOverlapFirstDerivs,
                                                                      diatomicOverlapSecondDerivs,
                                                                      diatomicTwoElecTwoCoreFirstDerivs,
                                                                      diatomicTwoElecTwoCoreSecondDerivs);
            }
            // hessian element atomA < atomB
            for(int atomBIndex=atomAIndex+1; atomBIndex<this->molecule->GetNumberAtoms(); atomBIndex++){
               const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
               for(int axisB = XAxis; axisB<CartesianType_end; axisB++){
                  int l = atomBIndex*CartesianType_end + axisB;
                  hessianSCF[l][k] = 
                  hessianSCF[k][l] = this->GetHessianElementDifferentAtomsSCF(atomAIndex, 
                                                                              atomBIndex,
                                                                              static_cast<CartesianType>(axisA), 
                                                                              static_cast<CartesianType>(axisB), 
                                                                              isMassWeighted,
                                                                              orbitalElectronPopulation,
                                                                              orbitalElectronPopulationFirstDerivs,
                                                                              diatomicOverlapFirstDerivs,
                                                                              diatomicOverlapSecondDerivs,
                                                                              diatomicTwoElecTwoCoreFirstDerivs,
                                                                              diatomicTwoElecTwoCoreSecondDerivs);
               }
            }

         }
      }
      this->FreeTempMatricesEachThreadCalcHessianSCF(&diatomicOverlapFirstDerivs,
                                                     &diatomicOverlapSecondDerivs,
                                                     &diatomicTwoElecTwoCoreFirstDerivs, 
                                                     &diatomicTwoElecTwoCoreSecondDerivs);
//}

   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&orbitalElectronPopulationFirstDerivs, 
                                                 totalNumberAOs,
                                                 totalNumberAOs,
                                                 this->molecule->GetNumberAtoms(),
                                                 CartesianType_end);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&orbitalElectronPopulationFirstDerivs, 
                                              totalNumberAOs,
                                              totalNumberAOs,
                                              this->molecule->GetNumberAtoms(),
                                              CartesianType_end);
}

void Mndo::CalcOrbitalElectronPopulationFirstDerivatives(double**** orbitalElectronPopulationFirstDerivs) const{
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
      for(int mu=0; mu<totalNumberAOs; mu++){
         for(int nu=0; nu<totalNumberAOs; nu++){
            for(int atomAIndex=0; atomAIndex<this->molecule->GetNumberAtoms(); atomAIndex++){
               for(int axis=XAxis; axis<CartesianType_end; axis++){

                  int moI, moJ;
                  double nI, nJ;
                  int indexSolutionCPHF = atomAIndex*CartesianType_end+axis;
                  orbitalElectronPopulationFirstDerivs[mu][nu][atomAIndex][axis] = 0.0;
                  for(int k=0; k<nonRedundantQIndeces.size(); k++){
                     moI = nonRedundantQIndeces[k].moI;
                     moJ = nonRedundantQIndeces[k].moJ;
                     nI = moI<numberOcc ? 2.0 : 0.0;
                     nJ = moJ<numberOcc ? 2.0 : 0.0;
                     orbitalElectronPopulationFirstDerivs[mu][nu][atomAIndex][axis]
                        += (nJ-nI)*
                           (transposedFockMatrix[mu][moJ]*transposedFockMatrix[nu][moI]+
                            transposedFockMatrix[mu][moI]*transposedFockMatrix[nu][moJ])*
                           solutionsCPHF[indexSolutionCPHF][k];
                  }

               }
            }
         }
      }
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
   for(int atomAIndex=0; atomAIndex<this->molecule->GetNumberAtoms(); atomAIndex++){
      for(int axisA=XAxis; axisA<CartesianType_end; axisA++){
         int k=atomAIndex*CartesianType_end + axisA;
         this->CalcStaticFirstOrderFock(staticFirstOrderFocks[k], 
                                        nonRedundantQIndeces,
                                        redundantQIndeces,
                                        atomAIndex,
                                        static_cast<CartesianType>(axisA));
      }
   }
}

// clac right side hand of CPHF, (34) in [PT_1996]
// Derivative coordinates is "axisA" of atomA.
void Mndo::CalcStaticFirstOrderFock(double* staticFirstOrderFock,
                                    const vector<MoIndexPair>& nonRedundantQIndeces,
                                    const vector<MoIndexPair>& redundantQIndeces,
                                    int atomAIndex,
                                    CartesianType axisA) const{
   MallocerFreer::GetInstance()->Initialize<double>(staticFirstOrderFock,
                                                    nonRedundantQIndeces.size()+redundantQIndeces.size());
   double***** diatomicTwoElecTwoCoreFirstDerivs = NULL;
   double*** diatomicOverlapFirstDerivs = NULL;
   
   try{
      this->MallocTempMatricesStaticFirstOrderFock(&diatomicTwoElecTwoCoreFirstDerivs, &diatomicOverlapFirstDerivs);
      const Atom& atomA = *molecule->GetAtom(atomAIndex);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int numberAOsA = atomA.GetValenceSize();
      int coreChargeA = atomA.GetCoreCharge();
      for(int atomBIndex=0; atomBIndex<this->molecule->GetNumberAtoms(); atomBIndex++){
         if(atomAIndex != atomBIndex){
            const Atom& atomB = *molecule->GetAtom(atomBIndex);
            int firstAOIndexB = atomB.GetFirstAOIndex();
            int numberAOsB = atomB.GetValenceSize();
            int coreChargeB = atomB.GetCoreCharge();

            // calc. first derivative of two elec two core interaction
            this->CalcDiatomicTwoElecTwoCoreFirstDerivatives(diatomicTwoElecTwoCoreFirstDerivs, atomAIndex, atomBIndex);
            // calc. first derivative of overlap.
            this->CalcDiatomicOverlapFirstDerivatives(diatomicOverlapFirstDerivs, atomA, atomB);

            // calc. static first order Fock;
            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
                  for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                     for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){

                        for(int i=0; i<nonRedundantQIndeces.size()+redundantQIndeces.size();i++){
                           int moI;
                           int moJ;
                           if(i<nonRedundantQIndeces.size()){
                              moI = nonRedundantQIndeces[i].moI;
                              moJ = nonRedundantQIndeces[i].moJ;
                           }
                           else{
                              moI = redundantQIndeces[i-nonRedundantQIndeces.size()].moI;
                              moJ = redundantQIndeces[i-nonRedundantQIndeces.size()].moJ;
                           }
                           double temp1 = this->fockMatrix[mu][moI]
                                         *this->fockMatrix[nu][moJ]
                                         *this->orbitalElectronPopulation[lambda][sigma]
                                         +this->fockMatrix[lambda][moI]
                                         *this->fockMatrix[sigma][moJ]
                                         *this->orbitalElectronPopulation[mu][nu]
                                         -0.5
                                         *this->fockMatrix[mu][moI]
                                         *this->fockMatrix[lambda][moJ]
                                         *this->orbitalElectronPopulation[nu][sigma]
                                         -0.5
                                         *this->fockMatrix[lambda][moI]
                                         *this->fockMatrix[mu][moJ]
                                         *this->orbitalElectronPopulation[nu][sigma];
                           staticFirstOrderFock[i] += temp1*diatomicTwoElecTwoCoreFirstDerivs[mu-firstAOIndexA]
                                                                                             [nu-firstAOIndexA]
                                                                                             [lambda-firstAOIndexB]
                                                                                             [sigma-firstAOIndexB]
                                                                                             [axisA];
                        } //i-loop

                     } //sigma-loop
                  } // lambda-loop

                  for(int i=0; i<nonRedundantQIndeces.size()+redundantQIndeces.size();i++){
                     int moI;
                     int moJ;
                     if(i<nonRedundantQIndeces.size()){
                        moI = nonRedundantQIndeces[i].moI;
                        moJ = nonRedundantQIndeces[i].moJ;
                     }
                     else{
                        moI = redundantQIndeces[i-nonRedundantQIndeces.size()].moI;
                        moJ = redundantQIndeces[i-nonRedundantQIndeces.size()].moJ;
                     }
                     double temp2 = this->fockMatrix[mu][moI]
                                   *this->fockMatrix[nu][moJ]
                                   *coreChargeB
                                   *diatomicTwoElecTwoCoreFirstDerivs[mu-firstAOIndexA]
                                                                     [nu-firstAOIndexA]
                                                                     [s]
                                                                     [s]
                                                                     [axisA];
                     staticFirstOrderFock[i] -= temp2;
                  }

               } // nu-loop
            } // mu-loop

            for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
               for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
                  
                  for(int i=0; i<nonRedundantQIndeces.size()+redundantQIndeces.size();i++){
                     int moI;
                     int moJ;
                     if(i<nonRedundantQIndeces.size()){
                        moI = nonRedundantQIndeces[i].moI;
                        moJ = nonRedundantQIndeces[i].moJ;
                     }
                     else{
                        moI = redundantQIndeces[i-nonRedundantQIndeces.size()].moI;
                        moJ = redundantQIndeces[i-nonRedundantQIndeces.size()].moJ;
                     }
                     double temp3 = this->fockMatrix[lambda][moI]
                                   *this->fockMatrix[sigma][moJ]
                                   *coreChargeA
                                   *diatomicTwoElecTwoCoreFirstDerivs[s]
                                                                     [s]
                                                                     [lambda-firstAOIndexB]
                                                                     [sigma-firstAOIndexB]
                                                                     [axisA];
                     staticFirstOrderFock[i] -= temp3;
                  } //i-loop

               } //sigma-loop
            } // lambda-loop

            for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
               for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
                  double bondParameter = 0.5*(atomA.GetBondingParameter(this->theory, 
                                                                        atomA.GetValence(mu-firstAOIndexA)) 
                                             +atomB.GetBondingParameter(this->theory, 
                                                                        atomB.GetValence(lambda-firstAOIndexB))); 
                  for(int i=0; i<nonRedundantQIndeces.size()+redundantQIndeces.size();i++){
                     int moI;
                     int moJ;
                     if(i<nonRedundantQIndeces.size()){
                        moI = nonRedundantQIndeces[i].moI;
                        moJ = nonRedundantQIndeces[i].moJ;
                     }
                     else{
                        moI = redundantQIndeces[i-nonRedundantQIndeces.size()].moI;
                        moJ = redundantQIndeces[i-nonRedundantQIndeces.size()].moJ;
                     }
                     double temp4 = ( this->fockMatrix[mu][moI]
                                     *this->fockMatrix[lambda][moJ]
                                     +this->fockMatrix[lambda][moI]
                                     *this->fockMatrix[mu][moJ]
                                    )
                                    *bondParameter
                                    *diatomicOverlapFirstDerivs[mu-firstAOIndexA][lambda-firstAOIndexB][axisA];
                     staticFirstOrderFock[i] += temp4;
                  }
               } //lambda-loop
            } // mu-loop
         }
      }

   }
   catch(MolDSException ex){
      this->FreeTempMatricesStaticFirstOrderFock(&diatomicTwoElecTwoCoreFirstDerivs, &diatomicOverlapFirstDerivs);
      throw ex;
   }
   this->FreeTempMatricesStaticFirstOrderFock(&diatomicTwoElecTwoCoreFirstDerivs, &diatomicOverlapFirstDerivs);

}

void Mndo::MallocTempMatricesStaticFirstOrderFock(double****** diatomicTwoElecTwoCoreFirstDeriv,
                                                  double**** diatomicOverlapFirstDeriv)const{
   MallocerFreer::GetInstance()->Malloc<double>(diatomicTwoElecTwoCoreFirstDeriv,
                                                dxy,
                                                dxy,
                                                dxy,
                                                dxy,
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicOverlapFirstDeriv,
                                                OrbitalType_end, 
                                                OrbitalType_end, 
                                                CartesianType_end);
}

void Mndo::FreeTempMatricesStaticFirstOrderFock(double****** diatomicTwoElecTwoCoreFirstDeriv,
                                                double**** diatomicOverlapFirstDeriv)const{
   MallocerFreer::GetInstance()->Free<double>(diatomicTwoElecTwoCoreFirstDeriv,
                                              dxy,
                                              dxy,
                                              dxy,
                                              dxy,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicOverlapFirstDeriv,
                                              OrbitalType_end, 
                                              OrbitalType_end, 
                                              CartesianType_end);
}

// see (40) - (46) in [PT_1996].
// This method calculates "(\Gamma - K)N" to solve CPHF (34) in [PT_1966]
void Mndo::CalcMatrixCPHF(double** matrixCPHF, 
                          const vector<MoIndexPair>& nonRedundantQIndeces,
                          const vector<MoIndexPair>& redundantQIndeces) const{
   double* occupations = NULL;
   MallocerFreer::GetInstance()->Malloc<double>(&occupations, nonRedundantQIndeces.size()+redundantQIndeces.size());
   try{   
      // calc diagonal part of N
      for(int i=0; i<nonRedundantQIndeces.size(); i++){
         int moI = nonRedundantQIndeces[i].moI;
         int moJ = nonRedundantQIndeces[i].moJ;
         occupations[i] = this->GetNNRElement(moI, moJ, moI, moJ);
      }
      for(int i=nonRedundantQIndeces.size(); i<nonRedundantQIndeces.size()+redundantQIndeces.size(); i++){
         int moI = redundantQIndeces[i-nonRedundantQIndeces.size()].moI;
         int moJ = redundantQIndeces[i-nonRedundantQIndeces.size()].moJ;
         occupations[i] = this->GetNRElement(moI, moJ, moI, moJ);
      }

      // calc (\Gamma - K)N
      for(int i=0; i<nonRedundantQIndeces.size(); i++){
         int moI = nonRedundantQIndeces[i].moI;
         int moJ = nonRedundantQIndeces[i].moJ;
         for(int j=i; j<nonRedundantQIndeces.size(); j++){
            int moK = nonRedundantQIndeces[j].moI;
            int moL = nonRedundantQIndeces[j].moJ;
            matrixCPHF[i][j] = (this->GetGammaNRElement(moI, moJ, moK, moL)-this->GetKNRElement(moI, moJ, moK, moL))
                              *occupations[j];

         }
         for(int j=0; j<i; j++){
            matrixCPHF[j][i] = matrixCPHF[i][j];
         }
      }

      for(int i=nonRedundantQIndeces.size(); i<nonRedundantQIndeces.size()+redundantQIndeces.size(); i++){
         int moI = redundantQIndeces[i-nonRedundantQIndeces.size()].moI;
         int moJ = redundantQIndeces[i-nonRedundantQIndeces.size()].moJ;
         for(int j=0; j<nonRedundantQIndeces.size(); j++){
            int moK = nonRedundantQIndeces[j].moI;
            int moL = nonRedundantQIndeces[j].moJ;
            matrixCPHF[i][j] = -1.0*this->GetKRElement(moI, moJ, moK, moL)*occupations[j];
         }
      }

      for(int i=nonRedundantQIndeces.size(); i<nonRedundantQIndeces.size()+redundantQIndeces.size(); i++){
         int moI = redundantQIndeces[i-nonRedundantQIndeces.size()].moI;
         int moJ = redundantQIndeces[i-nonRedundantQIndeces.size()].moJ;
         matrixCPHF[i][i] = this->GetGammaRElement(moI, moJ, moI, moJ)*occupations[i];
      }
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&occupations, nonRedundantQIndeces.size()+redundantQIndeces.size());
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&occupations, nonRedundantQIndeces.size()+redundantQIndeces.size());
}

void Mndo::MallocTempMatricesSolveCPHF(double*** matrixCPHF,
                                       int dimensionCPHF) const{
   MallocerFreer::GetInstance()->Malloc<double>(matrixCPHF, dimensionCPHF, dimensionCPHF);
}

void Mndo::FreeTempMatricesSolveCPHF(double*** matrixCPHF,
                                     int dimensionCPHF) const{
   MallocerFreer::GetInstance()->Free<double>(matrixCPHF, dimensionCPHF, dimensionCPHF);
}

// see [PT_1996, PT_1997]
void Mndo::CalcZMatrixForce(const vector<int>& elecStates){
   if(this->etaMatrixForce == NULL){
      stringstream ss;
      ss << this->errorMessageCalcZMatrixForceEtaNull;
      throw MolDSException(ss.str());
   }
   this->CheckZMatrixForce(elecStates); 

   // creat MO-index-pair for Q variables. 
   vector<MoIndexPair> nonRedundantQIndeces;
   vector<MoIndexPair> redundantQIndeces;
   this->CalcActiveSetVariablesQ(&nonRedundantQIndeces, 
                                 &redundantQIndeces,
                                 Parameters::GetInstance()->GetActiveOccCIS(),
                                 Parameters::GetInstance()->GetActiveVirCIS());

   // malloc temporary arraies
   double* delta = NULL; // Delta matrix, see (9) in [PT_1997]
   double* q = NULL; //// Q-vector in (19) in [PT_1997]
   double** gammaNRMinusKNR = NULL; // Gmamma_{NR} - K_{NR} matrix, see (40) and (45) to slove (54) in [PT_1996]
   double** kRDagerGammaRInv = NULL; // K_{R}^{\dagger} * Gamma_{R} matrix, see (41), (42), and (46) to solve (54) in [PT_1996]
   double* y = NULL; // y-vector in (54) in [PT_1996]
   double** transposedFockMatrix = NULL; // transposed Fock matrix
   double** xiOcc = NULL;
   double** xiVir = NULL;
   try{
      this->MallocTempMatrixForZMatrix(&delta,
                                       &q,
                                       &gammaNRMinusKNR,
                                       &kRDagerGammaRInv,
                                       &y,
                                       &transposedFockMatrix,
                                       &xiOcc,
                                       &xiVir,
                                       nonRedundantQIndeces.size(),
                                       redundantQIndeces.size());
      this->TransposeFockMatrixMatrix(transposedFockMatrix);
      this->CalcGammaNRMinusKNRMatrix(gammaNRMinusKNR, nonRedundantQIndeces);
      this->CalcKRDagerGammaRInvMatrix(kRDagerGammaRInv, nonRedundantQIndeces,redundantQIndeces);
      int groundState=0;
      for(int n=0; n<elecStates.size(); n++){
         if(groundState < elecStates[n]){
            int exciteState = elecStates[n]-1;
            this->CalcDeltaVector(delta, exciteState);
            this->CalcXiMatrices(xiOcc, xiVir, exciteState, transposedFockMatrix);
            this->CalcQVector(q, 
                              delta, 
                              xiOcc, 
                              xiVir,
                              this->etaMatrixForce[n],
                              nonRedundantQIndeces, 
                              redundantQIndeces);
            this->CalcAuxiliaryVector(y, q, kRDagerGammaRInv, nonRedundantQIndeces, redundantQIndeces);
            // solve (54) in [PT_1996]
            MolDS_wrappers::Lapack::GetInstance()->Dsysv(gammaNRMinusKNR, 
                                                         y, 
                                                         nonRedundantQIndeces.size());
            // calculate each element of Z matrix.
            stringstream ompErrors;
            #pragma omp parallel for schedule(auto)
            for(int mu=0; mu<this->molecule->GetTotalNumberAOs(); mu++){
               try{
                  for(int nu=0; nu<this->molecule->GetTotalNumberAOs(); nu++){
                     this->zMatrixForce[n][mu][nu] = this->GetZMatrixForceElement(
                                                           y,
                                                           q,
                                                           transposedFockMatrix,
                                                           nonRedundantQIndeces,
                                                           redundantQIndeces,
                                                           mu,
                                                           nu);
                  }
               }
               catch(MolDSException ex){
                  #pragma omp critical
                  ompErrors << ex.what() << endl ;
               }
            }
            // Exception throwing for omp-region
            if(!ompErrors.str().empty()){
               throw MolDSException(ompErrors.str());
            }

         }
      }
   }
   catch(MolDSException ex){
      this->FreeTempMatrixForZMatrix(&delta,
                                     &q,
                                     &gammaNRMinusKNR,
                                     &kRDagerGammaRInv,
                                     &y,
                                     &transposedFockMatrix,
                                     &xiOcc,
                                     &xiVir,
                                     nonRedundantQIndeces.size(),
                                     redundantQIndeces.size());
      throw ex;
   }
   this->FreeTempMatrixForZMatrix(&delta,
                                  &q,
                                  &gammaNRMinusKNR,
                                  &kRDagerGammaRInv,
                                  &y,
                                  &transposedFockMatrix,
                                  &xiOcc,
                                  &xiVir,
                                  nonRedundantQIndeces.size(),
                                  redundantQIndeces.size());
}

void Mndo::CalcEtaMatrixForce(const vector<int>& elecStates){
   this->CheckEtaMatrixForce(elecStates); 
   int numberAOs = this->molecule->GetTotalNumberAOs();
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   int groundState = 0;
   double** transposedFockMatrix = NULL; // transposed Fock matrix
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&transposedFockMatrix,
                                                   numberAOs,
                                                   numberAOs);
      this->TransposeFockMatrixMatrix(transposedFockMatrix);
      for(int n=0; n<elecStates.size(); n++){
         if(groundState < elecStates[n]){
            int exciteState = elecStates[n]-1;

            // calc each element
            stringstream ompErrors;
            #pragma omp parallel for schedule(auto)
            for(int mu=0; mu<numberAOs; mu++){
               try{
                  for(int nu=0; nu<numberAOs; nu++){
                     for(int i=0; i<numberActiveOcc; i++){
                        int moI = numberOcc-(i+1);
                        for(int a=0; a<numberActiveVir; a++){
                           int moA = numberOcc+a;
                           int slaterDeterminantIndex = this->GetSlaterDeterminantIndex(i,a);
                           this->etaMatrixForce[n][mu][nu] 
                                    += this->matrixCIS[exciteState][slaterDeterminantIndex]
                                      *transposedFockMatrix[mu][moI]
                                      *transposedFockMatrix[nu][moA];
                        }
                     }
                  }
               }
               catch(MolDSException ex){
                  #pragma omp critical
                  ompErrors << ex.what() << endl ;
               }
            }
            // Exception throwing for omp-region
            if(!ompErrors.str().empty()){
               throw MolDSException(ompErrors.str());
            }

         }
      }
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&transposedFockMatrix,numberAOs,numberAOs);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&transposedFockMatrix,numberAOs, numberAOs);
}

bool Mndo::RequiresExcitedStatesForce(const vector<int>& elecStates) const{
   bool requires = true;
   if(elecStates.size()==1 && elecStates[0]==0){
      requires = false;
   }
   return requires;
}

void Mndo::CalcForceSCFElecCoreAttractionPart(double* force, 
                                             int atomAIndex, 
                                             int atomBIndex,
                                             double const* const* const* const* const* diatomicTwoElecTwoCoreFirstDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int numberAOsA = atomA.GetValenceSize();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
         for(int i=0; i<CartesianType_end; i++){
            force[i] += -1.0
                       *this->orbitalElectronPopulation[mu][nu]
                       *this->GetElectronCoreAttractionFirstDerivative(
                                   atomAIndex, 
                                   atomBIndex, 
                                   mu-firstAOIndexA, 
                                   nu-firstAOIndexA,
                                   diatomicTwoElecTwoCoreFirstDerivs,
                                   (CartesianType)i);
         }
      }
   }
}

void Mndo::CalcForceSCFOverlapPart(double* force, 
                                  int atomAIndex, 
                                  int atomBIndex,
                                  double const* const* const* diatomicOverlapFirstDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int numberAOsA = atomA.GetValenceSize();
   int numberAOsB = atomB.GetValenceSize();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
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
                       *diatomicOverlapFirstDerivs[mu-firstAOIndexA][nu-firstAOIndexB][i];
         }
      }
   }
}

void Mndo::CalcForceSCFTwoElecPart(double* force, 
                                  int atomAIndex, 
                                  int atomBIndex,
                                  double const* const* const* const* const* diatomicTwoElecTwoCoreFirstDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int numberAOsA = atomA.GetValenceSize();
   int numberAOsB = atomB.GetValenceSize();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
         for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
            for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
               for(int i=0; i<CartesianType_end; i++){
                  force[i] -= 0.5
                             *this->orbitalElectronPopulation[mu][nu]
                             *this->orbitalElectronPopulation[lambda][sigma]
                             *diatomicTwoElecTwoCoreFirstDerivs[mu-firstAOIndexA]
                                                               [nu-firstAOIndexA]
                                                               [lambda-firstAOIndexB]
                                                               [sigma-firstAOIndexB]
                                                               [(CartesianType)i];
                  force[i] += 0.25
                             *this->orbitalElectronPopulation[mu][lambda]
                             *this->orbitalElectronPopulation[nu][sigma]
                             *diatomicTwoElecTwoCoreFirstDerivs[mu-firstAOIndexA]
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
                                      int atomAIndex, 
                                      int atomBIndex,
                                      double const* const* const* const* const* diatomicTwoElecTwoCoreFirstDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int numberAOsA = atomA.GetValenceSize();
   int numberAOsB = atomB.GetValenceSize();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
         for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
            for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
               for(int i=0; i<CartesianType_end; i++){
                  double temp= 2.0*this->etaMatrixForce[elecStateIndex][mu][nu]
                                  *this->etaMatrixForce[elecStateIndex][lambda][sigma]
                              -1.0*this->etaMatrixForce[elecStateIndex][mu][lambda]
                                  *this->etaMatrixForce[elecStateIndex][nu][sigma];
                  force[i] += temp
                             *diatomicTwoElecTwoCoreFirstDerivs[mu-firstAOIndexA]
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
                                                  int atomAIndex, 
                                                  int atomBIndex,
                                                  double const* const* const* const* const* diatomicTwoElecTwoCoreFirstDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int numberAOsA = atomA.GetValenceSize();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
         for(int i=0; i<CartesianType_end; i++){
            force[i] += -1.0
                       *this->zMatrixForce[elecStateIndex][mu][nu]
                       *this->GetElectronCoreAttractionFirstDerivative(
                                   atomAIndex, 
                                   atomBIndex, 
                                   mu-firstAOIndexA, 
                                   nu-firstAOIndexA,
                                   diatomicTwoElecTwoCoreFirstDerivs,
                                   (CartesianType)i);
         }
      }
   }
}

void Mndo::CalcForceExcitedOverlapPart(double* force, 
                                       int elecStateIndex,
                                       int atomAIndex, 
                                       int atomBIndex,
                                       double const* const* const* diatomicOverlapFirstDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int numberAOsA = atomA.GetValenceSize();
   int numberAOsB = atomB.GetValenceSize();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
         double bondParameter = atomA.GetBondingParameter(
                                      this->theory, 
                                      atomA.GetValence(mu-firstAOIndexA)) 
                               +atomB.GetBondingParameter(
                                      this->theory, 
                                      atomB.GetValence(nu-firstAOIndexB)); 
         bondParameter *= 0.5;
         for(int i=0; i<CartesianType_end; i++){
            force[i] += -1.0
                       *this->zMatrixForce[elecStateIndex][mu][nu]
                       *bondParameter
                       *diatomicOverlapFirstDerivs[mu-firstAOIndexA][nu-firstAOIndexB][i];
         }
      }
   }
}

void Mndo::CalcForceExcitedTwoElecPart(double* force, 
                                       int elecStateIndex,
                                       int atomAIndex, 
                                       int atomBIndex,
                                       double const* const* const* const* const* diatomicTwoElecTwoCoreFirstDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int numberAOsA = atomA.GetValenceSize();
   int numberAOsB = atomB.GetValenceSize();
   for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
      for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
         for(int lambda=firstAOIndexB; lambda<firstAOIndexB+numberAOsB; lambda++){
            for(int sigma=firstAOIndexB; sigma<firstAOIndexB+numberAOsB; sigma++){
               for(int i=0; i<CartesianType_end; i++){
                  force[i] -= this->zMatrixForce[elecStateIndex][mu][nu]
                             *this->orbitalElectronPopulation[lambda][sigma]
                             *diatomicTwoElecTwoCoreFirstDerivs[mu-firstAOIndexA]
                                                               [nu-firstAOIndexA]
                                                               [lambda-firstAOIndexB]
                                                               [sigma-firstAOIndexB]
                                                               [i];
                  force[i] += 0.50
                             *this->zMatrixForce[elecStateIndex][mu][lambda]
                             *this->orbitalElectronPopulation[nu][sigma]
                             *diatomicTwoElecTwoCoreFirstDerivs[mu-firstAOIndexA]
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
      double***** diatomicTwoElecTwoCoreFirstDerivs = NULL;
      double*** diatomicOverlapFirstDerivs = NULL;
      try{
         this->MallocTempMatricesCalcForce(&diatomicOverlapFirstDerivs, &diatomicTwoElecTwoCoreFirstDerivs);
         #pragma omp for schedule(auto)
         for(int a=0; a<this->molecule->GetNumberAtoms(); a++){
            const Atom& atomA = *molecule->GetAtom(a);
            int firstAOIndexA = atomA.GetFirstAOIndex();
            int numberAOsA = atomA.GetValenceSize();
            for(int b=0; b<this->molecule->GetNumberAtoms(); b++){
               if(a != b){
                  const Atom& atomB = *molecule->GetAtom(b);
                  int firstAOIndexB = atomB.GetFirstAOIndex();
                  int numberAOsB = atomB.GetValenceSize();

                  // calc. first derivative of overlap.
                  this->CalcDiatomicOverlapFirstDerivatives(diatomicOverlapFirstDerivs, atomA, atomB);
                  // calc. first derivative of two elec two core interaction
                  this->CalcDiatomicTwoElecTwoCoreFirstDerivatives(diatomicTwoElecTwoCoreFirstDerivs, 
                                                                   a, 
                                                                   b);

                  // core repulsion part
                  double coreRepulsion[CartesianType_end] = {0.0,0.0,0.0};
                  for(int i=0; i<CartesianType_end; i++){
                     coreRepulsion[i] += this->GetDiatomCoreRepulsionFirstDerivative(
                                               a, b, (CartesianType)i);
                     if(Parameters::GetInstance()->RequiresVdWSCF()){
                        coreRepulsion[i] += this->GetDiatomVdWCorrectionFirstDerivative(
                                                  a, b, (CartesianType)i);
                     }
                  }  
                  // electron core attraction part (ground state)
                  double forceElecCoreAttPart[CartesianType_end] = {0.0,0.0,0.0};
                  this->CalcForceSCFElecCoreAttractionPart(forceElecCoreAttPart,
                                                          a,
                                                          b,
                                                          diatomicTwoElecTwoCoreFirstDerivs);
                  // overlap part (ground state)
                  double forceOverlapPart[CartesianType_end] = {0.0,0.0,0.0};
                  this->CalcForceSCFOverlapPart(forceOverlapPart, 
                                               a,
                                               b,
                                               diatomicOverlapFirstDerivs);
                  // two electron part (ground state)
                  double forceTwoElecPart[CartesianType_end] = {0.0,0.0,0.0};
                  this->CalcForceSCFTwoElecPart(forceTwoElecPart,
                                               a,
                                               b,
                                               diatomicTwoElecTwoCoreFirstDerivs);
                  // sum up contributions from each part (ground state)
                  #pragma omp critical
                  {
                     for(int n=0; n<elecStates.size(); n++){
                        for(int i=0; i<CartesianType_end; i++){
                           this->matrixForce[n][a][i] -= coreRepulsion[i];
                           this->matrixForce[n][a][i] += forceElecCoreAttPart[i];
                           this->matrixForce[n][a][i] += forceOverlapPart[i];
                           this->matrixForce[n][a][i] += forceTwoElecPart[i];
                           this->matrixForce[n][b][i] -= forceElecCoreAttPart[i];
                           this->matrixForce[n][b][i] -= forceOverlapPart[i];
                           this->matrixForce[n][b][i] -= forceTwoElecPart[i];
                        }
                     }
                  }
                  // excited state potential
                  for(int n=0; n<elecStates.size(); n++){
                     if(0<elecStates[n]){
                        // static part
                        double forceExcitedStaticPart[CartesianType_end] = {0.0,0.0,0.0};
                        this->CalcForceExcitedStaticPart(forceExcitedStaticPart,
                                                         n,
                                                         a,
                                                         b,
                                                         diatomicTwoElecTwoCoreFirstDerivs);
                        // sum up contributions from static part (excited state)
                        #pragma omp critical
                        {
                           for(int i=0; i<CartesianType_end; i++){
                              this->matrixForce[n][b][i] += forceExcitedStaticPart[i];
                              this->matrixForce[n][a][i] -= forceExcitedStaticPart[i];
                           }
                        }

                        // response part
                        // electron core attraction part (excited state)
                        double forceExcitedElecCoreAttPart[CartesianType_end]={0.0,0.0,0.0};
                        this->CalcForceExcitedElecCoreAttractionPart(
                                                   forceExcitedElecCoreAttPart,
                                                   n,
                                                   a,
                                                   b,
                                                   diatomicTwoElecTwoCoreFirstDerivs);
                        // overlap part (excited state)
                        double forceExcitedOverlapPart[CartesianType_end] = {0.0,0.0,0.0};
                        this->CalcForceExcitedOverlapPart(forceExcitedOverlapPart, 
                                                          n,
                                                          a,
                                                          b,
                                                          diatomicOverlapFirstDerivs);
                        // two electron part (ground state)
                        double forceExcitedTwoElecPart[CartesianType_end] = {0.0,0.0,0.0};
                        this->CalcForceExcitedTwoElecPart(forceExcitedTwoElecPart,
                                                          n,
                                                          a,
                                                          b,
                                                          diatomicTwoElecTwoCoreFirstDerivs);
                        // sum up contributions from response part (excited state)
                        #pragma omp critical
                        {
                           for(int i=0; i<CartesianType_end; i++){
                              this->matrixForce[n][a][i] += forceExcitedElecCoreAttPart[i];
                              this->matrixForce[n][a][i] += forceExcitedOverlapPart[i];
                              this->matrixForce[n][a][i] += forceExcitedTwoElecPart[i];
                              this->matrixForce[n][b][i] -= forceExcitedElecCoreAttPart[i];
                              this->matrixForce[n][b][i] -= forceExcitedOverlapPart[i];
                              this->matrixForce[n][b][i] -= forceExcitedTwoElecPart[i];
                           }
                        }

                     }
                  }
               }
            }
         }
      }
      catch(MolDSException ex){
         #pragma omp critical
         ompErrors << ex.what() << endl ;
      }
      this->FreeTempMatricesCalcForce(&diatomicOverlapFirstDerivs, &diatomicTwoElecTwoCoreFirstDerivs);
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
}

void Mndo::MallocTempMatricesCalcForce(double**** diatomicOverlapFirstDerivs, double****** diatomicTwoElecTwoCoreFirstDerivs) const{
   MallocerFreer::GetInstance()->Malloc<double>(diatomicOverlapFirstDerivs, 
                                                OrbitalType_end,
                                                OrbitalType_end,
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicTwoElecTwoCoreFirstDerivs,
                                                dxy,
                                                dxy,
                                                dxy,
                                                dxy,
                                                CartesianType_end);
}

void Mndo::FreeTempMatricesCalcForce(double**** diatomicOverlapFirstDerivs, double****** diatomicTwoElecTwoCoreFirstDerivs) const{
   MallocerFreer::GetInstance()->Free<double>(diatomicOverlapFirstDerivs, 
                                              OrbitalType_end,
                                              OrbitalType_end,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicTwoElecTwoCoreFirstDerivs,
                                              dxy,
                                              dxy,
                                              dxy,
                                              dxy,
                                              CartesianType_end);
}

void Mndo::CalcTwoElecTwoCore(double****** twoElecTwoCore, 
                              const Molecule& molecule) const{
   if(twoElecTwoCore == NULL){
      stringstream ss;
      ss << this->errorMessageCalcTwoElecTwoCoreNullMatrix;
      throw MolDSException(ss.str());
   }
   else{
      MallocerFreer::GetInstance()->Initialize<double>(twoElecTwoCore, 
                                                       molecule.GetNumberAtoms(),
                                                       molecule.GetNumberAtoms(),
                                                       dxy, dxy, dxy, dxy);
   } 

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
         ompErrors << ex.what() << endl ;
      }
      MallocerFreer::GetInstance()->Free<double>(&diatomicTwoElecTwoCore, dxy, dxy, dxy, dxy);
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
}

// Calculation of two electrons two cores integral (mu, nu | lambda, sigma) in space fixed frame, 
// taht is, Eq. (9) in ref. [DT_1977-2].
// mu and nu are included in atomA's AOs. 
// lambda and sigma are included in atomB's AOs.
// Note that atomA != atomB.
// Note taht d-orbital cannot be treated, 
// that is, matrix[dxy][dxy][dxy][dxy] cannot be treatable.
void Mndo::CalcDiatomicTwoElecTwoCore(double**** matrix, int atomAIndex, int atomBIndex) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   if(atomAIndex == atomBIndex){
      stringstream ss;
      ss << this->errorMessageCalcDiatomicTwoElecTwoCoreSameAtoms;
      ss << this->errorMessageAtomA << atomAIndex 
                                    << AtomTypeStr(atomA.GetAtomType()) << endl;
      ss << this->errorMessageAtomB << atomBIndex 
                                    << AtomTypeStr(atomB.GetAtomType()) << endl;
      throw MolDSException(ss.str());
   }

   if(matrix == NULL){
      stringstream ss;
      ss << this->errorMessageCalcDiatomicTwoElecTwoCoreNullMatrix;
      throw MolDSException(ss.str());
   }
   else{
      MallocerFreer::GetInstance()->Initialize<double>(matrix, dxy, dxy, dxy, dxy);
   } 

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
      this->RotateDiatomicTwoElecTwoCoreToSpaceFramegc(matrix, rotatingMatrix);
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
               this->OutputLog((boost::format("mu=%d nu=%d lambda=%d sigma=%d $e\n") % mu
                                                                                     % nu
                                                                                     % lambda
                                                                                     % sigma 
                                                                                     % matrix[mu][nu][lambda][sigma]).str());
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
void Mndo::CalcDiatomicTwoElecTwoCoreFirstDerivatives(double***** matrix, 
                                                      int atomAIndex, 
                                                      int atomBIndex) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   if(atomAIndex == atomBIndex){
      stringstream ss;
      ss << this->errorMessageCalcDiatomicTwoElecTwoCoreFirstDerivativesSameAtoms;
      ss << this->errorMessageAtomA << atomAIndex 
                                    << AtomTypeStr(atomA.GetAtomType()) << endl;
      ss << this->errorMessageAtomB << atomBIndex 
                                    << AtomTypeStr(atomB.GetAtomType()) << endl;
      throw MolDSException(ss.str());
   }

   if(matrix == NULL){
      stringstream ss;
      ss << this->errorMessageCalcDiatomicTwoElecTwoCoreFirstDerivativesNullMatrix;
      throw MolDSException(ss.str());
   }
   else{
      MallocerFreer::GetInstance()->Initialize<double>(matrix, 
                                                       dxy, 
                                                       dxy, 
                                                       dxy, 
                                                       dxy, 
                                                       CartesianType_end);
   } 

   double** rotatingMatrix = NULL;
   double*** rotMatFirstDerivatives = NULL;
   double**** diatomicTwoElecTwoCore = NULL;
   try{
      this->MallocDiatomicTwoElecTwoCoreFirstDeriTemps(&rotatingMatrix,
                                                       &rotMatFirstDerivatives,
                                                       &diatomicTwoElecTwoCore);
      // calclation in diatomic frame
      for(int mu=0; mu<atomA.GetValenceSize(); mu++){
         for(int nu=0; nu<atomA.GetValenceSize(); nu++){
            for(int lambda=0; lambda<atomB.GetValenceSize(); lambda++){
               for(int sigma=0; sigma<atomB.GetValenceSize(); sigma++){
                  for(int dimA=0; dimA<CartesianType_end; dimA++){
                     matrix[mu][nu][lambda][sigma][dimA] 
                        = this->GetNddoRepulsionIntegralFirstDerivative(
                                atomA, 
                                atomA.GetValence(mu),
                                atomA.GetValence(nu),
                                atomB, 
                                atomB.GetValence(lambda),
                                atomB.GetValence(sigma),
                                static_cast<CartesianType>(dimA));
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
      this->CalcRotatingMatrixFirstDerivatives(rotMatFirstDerivatives, atomA, atomB);
      this->RotateDiatomicTwoElecTwoCoreFirstDerivativesToSpaceFramegc(matrix, 
                                                                       diatomicTwoElecTwoCore,
                                                                       rotatingMatrix,
                                                                       rotMatFirstDerivatives);
   }
   catch(MolDSException ex){
      this->FreeDiatomicTwoElecTwoCoreFirstDeriTemps(&rotatingMatrix,
                                                     &rotMatFirstDerivatives,
                                                     &diatomicTwoElecTwoCore);
      throw ex;
   }
   this->FreeDiatomicTwoElecTwoCoreFirstDeriTemps(&rotatingMatrix,
                                                  &rotMatFirstDerivatives,
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
void Mndo::CalcDiatomicTwoElecTwoCoreSecondDerivatives(double****** matrix, 
                                                       int atomAIndex, 
                                                       int atomBIndex) const{
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   if(atomAIndex == atomBIndex){
      stringstream ss;
      ss << this->errorMessageCalcDiatomicTwoElecTwoCoreSecondDerivativesSameAtoms;
      ss << this->errorMessageAtomA << atomAIndex 
                                    << AtomTypeStr(atomA.GetAtomType()) << endl;
      ss << this->errorMessageAtomB << atomBIndex 
                                    << AtomTypeStr(atomB.GetAtomType()) << endl;
      throw MolDSException(ss.str());
   }

   if(matrix == NULL){
      stringstream ss;
      ss << this->errorMessageCalcDiatomicTwoElecTwoCoreSecondDerivativesNullMatrix;
      throw MolDSException(ss.str());
   }
   else{
      MallocerFreer::GetInstance()->Initialize<double>(matrix, 
                                                       dxy, 
                                                       dxy, 
                                                       dxy, 
                                                       dxy, 
                                                       CartesianType_end,
                                                       CartesianType_end);
   } 

   double** rotatingMatrix = NULL;
   double*** rotMatFirstDerivatives = NULL;
   double**** rotMatSecondDerivatives = NULL;
   double**** diatomicTwoElecTwoCore = NULL;
   double***** diatomicTwoElecTwoCoreFirstDerivatives = NULL;
   try{
      this->MallocDiatomicTwoElecTwoCoreSecondDeriTemps(&rotatingMatrix,
                                                        &rotMatFirstDerivatives,
                                                        &rotMatSecondDerivatives,
                                                        &diatomicTwoElecTwoCore,
                                                        &diatomicTwoElecTwoCoreFirstDerivatives);
      // calclation in diatomic frame
      for(int mu=0; mu<atomA.GetValenceSize(); mu++){
         for(int nu=0; nu<atomA.GetValenceSize(); nu++){
            for(int lambda=0; lambda<atomB.GetValenceSize(); lambda++){
               for(int sigma=0; sigma<atomB.GetValenceSize(); sigma++){
                  for(int dimA1=0; dimA1<CartesianType_end; dimA1++){
                     for(int dimA2=0; dimA2<CartesianType_end; dimA2++){
                        matrix[mu][nu][lambda][sigma][dimA1][dimA2]
                           = this->GetNddoRepulsionIntegralSecondDerivative(
                                   atomA, 
                                   atomA.GetValence(mu),
                                   atomA.GetValence(nu),
                                   atomB, 
                                   atomB.GetValence(lambda),
                                   atomB.GetValence(sigma),
                                   static_cast<CartesianType>(dimA1),
                                   static_cast<CartesianType>(dimA2));
                     }
                     diatomicTwoElecTwoCoreFirstDerivatives[mu][nu][lambda][sigma][dimA1] 
                        = this->GetNddoRepulsionIntegralFirstDerivative(
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
      this->CalcRotatingMatrixFirstDerivatives(rotMatFirstDerivatives, atomA, atomB);
      this->CalcRotatingMatrixSecondDerivatives(rotMatSecondDerivatives, atomA, atomB);
      this->RotateDiatomicTwoElecTwoCoreSecondDerivativesToSpaceFramegc(matrix, 
                                                                        diatomicTwoElecTwoCore,
                                                                        diatomicTwoElecTwoCoreFirstDerivatives,
                                                                        rotatingMatrix,
                                                                        rotMatFirstDerivatives,
                                                                        rotMatSecondDerivatives);
   }
   catch(MolDSException ex){
      this->FreeDiatomicTwoElecTwoCoreSecondDeriTemps(&rotatingMatrix,
                                                      &rotMatFirstDerivatives,
                                                      &rotMatSecondDerivatives,
                                                      &diatomicTwoElecTwoCore,
                                                      &diatomicTwoElecTwoCoreFirstDerivatives);
      throw ex;
   }
   this->FreeDiatomicTwoElecTwoCoreSecondDeriTemps(&rotatingMatrix,
                                                   &rotMatFirstDerivatives,
                                                   &rotMatSecondDerivatives,
                                                   &diatomicTwoElecTwoCore,
                                                   &diatomicTwoElecTwoCoreFirstDerivatives);
}

void Mndo::MallocDiatomicTwoElecTwoCoreFirstDeriTemps(double*** rotatingMatrix,
                                                      double**** rotMatFirstDerivatives,
                                                      double***** diatomicTwoElecTwoCore) const{
   MallocerFreer::GetInstance()->Malloc<double>(rotatingMatrix,
                                                OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Malloc<double>(rotMatFirstDerivatives, 
                                                OrbitalType_end, 
                                                OrbitalType_end, 
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicTwoElecTwoCore, dxy, dxy, dxy, dxy);
}

void Mndo::MallocDiatomicTwoElecTwoCoreSecondDeriTemps(double*** rotatingMatrix,
                                                       double**** rotMatFirstDerivatives,
                                                       double***** rotMatSecondDerivatives,
                                                       double***** diatomicTwoElecTwoCore,
                                                       double****** diatomicTwoElecTwoCoreFirstDerivatives) const{
   this->MallocDiatomicTwoElecTwoCoreFirstDeriTemps(rotatingMatrix,
                                                    rotMatFirstDerivatives,
                                                    diatomicTwoElecTwoCore);
   MallocerFreer::GetInstance()->Malloc<double>(rotMatSecondDerivatives, 
                                                OrbitalType_end, 
                                                OrbitalType_end, 
                                                CartesianType_end,
                                                CartesianType_end);
   MallocerFreer::GetInstance()->Malloc<double>(diatomicTwoElecTwoCoreFirstDerivatives, 
                                                dxy, 
                                                dxy, 
                                                dxy, 
                                                dxy,
                                                CartesianType_end);
}

void Mndo::FreeDiatomicTwoElecTwoCoreFirstDeriTemps(double*** rotatingMatrix,
                                                    double**** rotMatFirstDerivatives,
                                                    double***** diatomicTwoElecTwoCore) const{
   MallocerFreer::GetInstance()->Free<double>(rotatingMatrix,
                                              OrbitalType_end, OrbitalType_end);
   MallocerFreer::GetInstance()->Free<double>(rotMatFirstDerivatives, 
                                              OrbitalType_end, 
                                              OrbitalType_end, 
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicTwoElecTwoCore, dxy, dxy, dxy, dxy);
}

void Mndo::FreeDiatomicTwoElecTwoCoreSecondDeriTemps(double*** rotatingMatrix,
                                                     double**** rotMatFirstDerivatives,
                                                     double***** rotMatSecondDerivatives,
                                                     double***** diatomicTwoElecTwoCore,
                                                     double****** diatomicTwoElecTwoCoreFirstDerivatives) const{
   this->FreeDiatomicTwoElecTwoCoreFirstDeriTemps(rotatingMatrix,
                                                  rotMatFirstDerivatives,
                                                  diatomicTwoElecTwoCore);
   MallocerFreer::GetInstance()->Free<double>(rotMatSecondDerivatives, 
                                              OrbitalType_end, 
                                              OrbitalType_end, 
                                              CartesianType_end,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(diatomicTwoElecTwoCoreFirstDerivatives, 
                                              dxy, 
                                              dxy, 
                                              dxy, 
                                              dxy,
                                              CartesianType_end);
}

// Rotate 4-dimensional matrix from diatomic frame to space frame
// Note tha in this method d-orbitals can not be treatable.
void Mndo::RotateDiatomicTwoElecTwoCoreToSpaceFramegc(double**** matrix, 
                                                      double const* const* rotatingMatrix) const{
   double oldMatrix[dxy][dxy][dxy][dxy];
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               oldMatrix[mu][nu][lambda][sigma] = matrix[mu][nu][lambda][sigma];
            }
         }
      }
   }
   
   // rotate (fast algorythm, see also slow algorythm shown later)
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               matrix[mu][nu][lambda][sigma] = 0.0;
               for(int i=0; i<dxy; i++){
                  double tempI = 0.0;
                  for(int j=0; j<dxy; j++){
                     double tempIJ = 0.0;
                     for(int k=0; k<dxy; k++){
                        double tempIJK = 0.0;
                        for(int l=0; l<dxy; l++){
                           tempIJK += oldMatrix[i][j][k][l]*rotatingMatrix[sigma][l];
                        }
                        tempIJ += tempIJK*rotatingMatrix[lambda][k];
                     }
                     tempI += tempIJ*rotatingMatrix[nu][j];
                  }
                  matrix[mu][nu][lambda][sigma] += tempI*rotatingMatrix[mu][i];
               }
            }
         }
      }
   }

   /*
   // rotate (slow algorythm)
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

// Rotate 5-dimensional matrix from diatomic frame to space frame
// Note tha in this method d-orbitals can not be treatable.
void Mndo::RotateDiatomicTwoElecTwoCoreFirstDerivativesToSpaceFramegc(
           double***** matrix, 
           double const* const* const* const* diatomicTwoElecTwoCore,
           double const* const* rotatingMatrix,
           double const* const* const* rotMatFirstDerivatives) const{
   double oldMatrix[dxy][dxy][dxy][dxy][CartesianType_end];
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               for(int c=0; c<CartesianType_end; c++){
                  oldMatrix[mu][nu][lambda][sigma][c] = matrix[mu][nu][lambda][sigma][c];
               }
            }
         }
      }
   }
   
   // rotate (fast algorythm, see also slow algorythm shown later)
   for(int mu=0; mu<dxy; mu++){
      for(int nu=0; nu<dxy; nu++){
         for(int lambda=0; lambda<dxy; lambda++){
            for(int sigma=0; sigma<dxy; sigma++){
               for(int c=0; c<CartesianType_end; c++){

                  matrix[mu][nu][lambda][sigma][c] = 0.0;
                  for(int i=0; i<dxy; i++){
                     double tempI_1 = 0.0;
                     double tempI_2 = 0.0;
                     double tempI_3 = 0.0;
                     double tempI_4 = 0.0;
                     double tempI_5 = 0.0;
                     for(int j=0; j<dxy; j++){
                        double tempIJ_1 = 0.0;
                        double tempIJ_2 = 0.0;
                        double tempIJ_3 = 0.0;
                        double tempIJ_4 = 0.0;
                        double tempIJ_5 = 0.0;
                        for(int k=0; k<dxy; k++){
                           double tempIJK_1 = 0.0;
                           double tempIJK_2 = 0.0;
                           double tempIJK_3 = 0.0;
                           double tempIJK_4 = 0.0;
                           double tempIJK_5 = 0.0;
                           for(int l=0; l<dxy; l++){
                              tempIJK_1 += oldMatrix[i][j][k][l][c]*rotatingMatrix[sigma][l];
                              tempIJK_2 += diatomicTwoElecTwoCore[i][j][k][l]*rotatingMatrix[sigma][l];
                              tempIJK_3 += diatomicTwoElecTwoCore[i][j][k][l]*rotatingMatrix[sigma][l];
                              tempIJK_4 += diatomicTwoElecTwoCore[i][j][k][l]*rotatingMatrix[sigma][l];
                              tempIJK_5 += diatomicTwoElecTwoCore[i][j][k][l]*rotMatFirstDerivatives[sigma][l][c];
                           }
                           tempIJ_1 += tempIJK_1*rotatingMatrix[lambda][k];
                           tempIJ_2 += tempIJK_2*rotatingMatrix[lambda][k];
                           tempIJ_3 += tempIJK_3*rotatingMatrix[lambda][k];
                           tempIJ_4 += tempIJK_4*rotMatFirstDerivatives[lambda][k][c];
                           tempIJ_5 += tempIJK_5*rotatingMatrix[lambda][k];
                        }
                        tempI_1 += tempIJ_1*rotatingMatrix[nu][j];
                        tempI_2 += tempIJ_2*rotatingMatrix[nu][j];
                        tempI_3 += tempIJ_3*rotMatFirstDerivatives[nu][j][c];
                        tempI_4 += tempIJ_4*rotatingMatrix[nu][j];
                        tempI_5 += tempIJ_5*rotatingMatrix[nu][j];
                     }
                     matrix[mu][nu][lambda][sigma][c] += tempI_1*rotatingMatrix[mu][i];
                     matrix[mu][nu][lambda][sigma][c] += tempI_2*rotMatFirstDerivatives[mu][i][c];
                     matrix[mu][nu][lambda][sigma][c] += tempI_3*rotatingMatrix[mu][i];
                     matrix[mu][nu][lambda][sigma][c] += tempI_4*rotatingMatrix[mu][i];
                     matrix[mu][nu][lambda][sigma][c] += tempI_5*rotatingMatrix[mu][i];
                  }

               }
            }
         }
      }
   }

   /*
   // rotate (slow algorythm)
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
                                   *rotMatFirstDerivatives[mu][i][c]
                                   *rotatingMatrix[nu][j] 
                                   *rotatingMatrix[lambda][k] 
                                   *rotatingMatrix[sigma][l];
                              matrix[mu][nu][lambda][sigma][c] 
                                 += diatomicTwoElecTwoCore[i][j][k][l]
                                   *rotatingMatrix[mu][i] 
                                   *rotMatFirstDerivatives[nu][j][c]
                                   *rotatingMatrix[lambda][k] 
                                   *rotatingMatrix[sigma][l];
                              matrix[mu][nu][lambda][sigma][c] 
                                 += diatomicTwoElecTwoCore[i][j][k][l]
                                   *rotatingMatrix[mu][i] 
                                   *rotatingMatrix[nu][j] 
                                   *rotMatFirstDerivatives[lambda][k][c]
                                   *rotatingMatrix[sigma][l];
                              matrix[mu][nu][lambda][sigma][c] 
                                 += diatomicTwoElecTwoCore[i][j][k][l]
                                   *rotatingMatrix[mu][i] 
                                   *rotatingMatrix[nu][j] 
                                   *rotatingMatrix[lambda][k] 
                                   *rotMatFirstDerivatives[sigma][l][c];
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

// Rotate 6-dimensional matrix from diatomic frame to space frame
// Note tha in this method d-orbitals can not be treatable.
void Mndo::RotateDiatomicTwoElecTwoCoreSecondDerivativesToSpaceFramegc(
           double****** matrix, 
           double const* const* const* const* diatomicTwoElecTwoCore,
           double const* const* const* const* const* diatomicTwoElecTwoCoreFirstDerivatives,
           double const* const* rotatingMatrix,
           double const* const* const* rotMatFirstDerivatives,
           double const* const* const* const* rotMatSecondDerivatives) const{
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
   
   // rotate (fast algorythm, see also slow algorythm shown later)
   int numberTerms = 25;
   double* tempIJK = NULL;
   double* tempIJ = NULL;
   double* tempI = NULL;
   MallocerFreer::GetInstance()->Malloc<double>(&tempIJK, numberTerms);
   MallocerFreer::GetInstance()->Malloc<double>(&tempIJ, numberTerms);
   MallocerFreer::GetInstance()->Malloc<double>(&tempI, numberTerms);
   try{ 
      for(int mu=s; mu<dxy; mu++){
         for(int nu=s; nu<dxy; nu++){
            for(int lambda=s; lambda<dxy; lambda++){
               for(int sigma=s; sigma<dxy; sigma++){
                  for(int dimA1=XAxis; dimA1<CartesianType_end; dimA1++){
                     for(int dimA2=XAxis; dimA2<CartesianType_end; dimA2++){

                        matrix[mu][nu][lambda][sigma][dimA1][dimA2] = 0.0;
                        double value=0.0;
                        for(int i=s; i<dxy; i++){
                           MallocerFreer::GetInstance()->Initialize<double>(tempI, numberTerms);
                           for(int j=s; j<dxy; j++){
                              MallocerFreer::GetInstance()->Initialize<double>(tempIJ, numberTerms);
                              for(int k=s; k<dxy; k++){
                                 MallocerFreer::GetInstance()->Initialize<double>(tempIJK, numberTerms);
                                 for(int l=s; l<dxy; l++){
                                    
                                    tempIJK[0]  += oldMatrix             [i][j][k][l][dimA1][dimA2]*rotatingMatrix         [sigma][l];
                                    tempIJK[1]  += diatomicTwoElecTwoCore[i][j][k][l]              *rotatingMatrix         [sigma][l];
                                    tempIJK[4]  += diatomicTwoElecTwoCore[i][j][k][l]              *rotMatSecondDerivatives[sigma][l][dimA1][dimA2];

                                    tempIJK[5]  += diatomicTwoElecTwoCoreFirstDerivatives[i][j][k][l][dimA1]*rotatingMatrix        [sigma][l];
                                    tempIJK[8]  += diatomicTwoElecTwoCoreFirstDerivatives[i][j][k][l][dimA1]*rotMatFirstDerivatives[sigma][l][dimA2];

                                    tempIJK[9]  += diatomicTwoElecTwoCoreFirstDerivatives[i][j][k][l][dimA2]*rotatingMatrix        [sigma][l];
                                    tempIJK[12] += diatomicTwoElecTwoCore                [i][j][k][l]       *rotMatFirstDerivatives[sigma][l][dimA2];

                                    tempIJK[21] += diatomicTwoElecTwoCoreFirstDerivatives[i][j][k][l][dimA2]*rotMatFirstDerivatives[sigma][l][dimA1];
                                    tempIJK[22] += diatomicTwoElecTwoCore                [i][j][k][l]       *rotMatFirstDerivatives[sigma][l][dimA1];
                                 }
                                 tempIJ[0]  += tempIJK[0] *rotatingMatrix         [lambda][k];
                                 tempIJ[1]  += tempIJK[1] *rotatingMatrix         [lambda][k];
                                 tempIJ[3]  += tempIJK[1] *rotMatSecondDerivatives[lambda][k][dimA1][dimA2];
                                 tempIJ[4]  += tempIJK[4] *rotatingMatrix         [lambda][k];
                                           
                                 tempIJ[5]  += tempIJK[5] *rotatingMatrix        [lambda][k];
                                 tempIJ[7]  += tempIJK[5] *rotMatFirstDerivatives[lambda][k][dimA2];
                                 tempIJ[8]  += tempIJK[8] *rotatingMatrix        [lambda][k];

                                 tempIJ[9]  += tempIJK[9] *rotatingMatrix        [lambda][k];
                                 tempIJ[11] += tempIJK[1] *rotMatFirstDerivatives[lambda][k][dimA2];
                                 tempIJ[12] += tempIJK[12]*rotatingMatrix        [lambda][k];

                                 tempIJ[17] += tempIJK[9] *rotMatFirstDerivatives[lambda][k][dimA1];
                                 tempIJ[18] += tempIJK[1 ]*rotMatFirstDerivatives[lambda][k][dimA1];
                                 tempIJ[20] += tempIJK[12]*rotMatFirstDerivatives[lambda][k][dimA1];

                                 tempIJ[21] += tempIJK[21]*rotatingMatrix        [lambda][k];
                                 tempIJ[22] += tempIJK[22]*rotatingMatrix        [lambda][k];
                                 tempIJ[24] += tempIJK[22]*rotMatFirstDerivatives[lambda][k][dimA2];
                              }
                              tempI[0]  += tempIJ[0] *rotatingMatrix         [nu][j];
                              tempI[1]  += tempIJ[1] *rotatingMatrix         [nu][j];
                              tempI[2]  += tempIJ[1] *rotMatSecondDerivatives[nu][j][dimA1][dimA2];
                              tempI[3]  += tempIJ[3] *rotatingMatrix         [nu][j];
                              tempI[4]  += tempIJ[4] *rotatingMatrix         [nu][j];
                                       
                              tempI[5]  += tempIJ[5] *rotatingMatrix        [nu][j];
                              tempI[6]  += tempIJ[5] *rotMatFirstDerivatives[nu][j][dimA2];
                              tempI[7]  += tempIJ[7] *rotatingMatrix        [nu][j];
                              tempI[8]  += tempIJ[8] *rotatingMatrix        [nu][j];
                                       
                              tempI[9]  += tempIJ[9] *rotatingMatrix        [nu][j];
                              tempI[10] += tempIJ[1] *rotMatFirstDerivatives[nu][j][dimA2];
                              tempI[11] += tempIJ[11]*rotatingMatrix        [nu][j];
                              tempI[12] += tempIJ[12]*rotatingMatrix        [nu][j];

                              tempI[13] += tempIJ[9] *rotMatFirstDerivatives[nu][j][dimA1];
                              tempI[14] += tempIJ[1] *rotMatFirstDerivatives[nu][j][dimA1];
                              tempI[15] += tempIJ[11]*rotMatFirstDerivatives[nu][j][dimA1];
                              tempI[16] += tempIJ[12]*rotMatFirstDerivatives[nu][j][dimA1];

                              tempI[17] += tempIJ[17]*rotatingMatrix        [nu][j];
                              tempI[18] += tempIJ[18]*rotatingMatrix        [nu][j];
                              tempI[19] += tempIJ[18]*rotMatFirstDerivatives[nu][j][dimA2];
                              tempI[20] += tempIJ[20]*rotatingMatrix        [nu][j];

                              tempI[21] += tempIJ[21]*rotatingMatrix        [nu][j];
                              tempI[22] += tempIJ[22]*rotatingMatrix        [nu][j];
                              tempI[23] += tempIJ[22]*rotMatFirstDerivatives[nu][j][dimA2];
                              tempI[24] += tempIJ[24]*rotatingMatrix        [nu][j];
                           }
                           value += (tempI[0] + tempI[2] + tempI[3] + tempI[4] + 
                                     tempI[6] + tempI[7] + tempI[8] + tempI[13] + 
                                     tempI[15] + tempI[16] + tempI[17] + tempI[19] + 
                                     tempI[20] + tempI[21] + tempI[23] + tempI[24])*rotatingMatrix[mu][i];
                           value += (tempI[9] + tempI[10] + tempI[11] + tempI[12]) *rotMatFirstDerivatives[mu][i][dimA1];
                           value += (tempI[5] + tempI[14] + tempI[18] + tempI[22]) *rotMatFirstDerivatives[mu][i][dimA2];
                           value += tempI[1] *rotMatSecondDerivatives[mu][i][dimA1][dimA2];
                        }
                        matrix[mu][nu][lambda][sigma][dimA1][dimA2] = value;


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
   // rotate (slow algorythm shown later)
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
                                         *rotMatSecondDerivatives[mu    ][i][dimA1][dimA2]
                                         *rotatingMatrix         [nu    ][j] 
                                         *rotatingMatrix         [lambda][k] 
                                         *rotatingMatrix         [sigma ][l];
                                 // term2
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix         [mu    ][i] 
                                         *rotMatSecondDerivatives[nu    ][j][dimA1][dimA2]
                                         *rotatingMatrix         [lambda][k] 
                                         *rotatingMatrix         [sigma ][l];
                                 // term3
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix         [mu    ][i] 
                                         *rotatingMatrix         [nu    ][j] 
                                         *rotMatSecondDerivatives[lambda][k][dimA1][dimA2]
                                         *rotatingMatrix         [sigma ][l];
                                 // term4
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix         [mu    ][i] 
                                         *rotatingMatrix         [nu    ][j] 
                                         *rotatingMatrix         [lambda][k] 
                                         *rotMatSecondDerivatives[sigma ][l][dimA1][dimA2];
                                 
                                 // term5
                                 value += diatomicTwoElecTwoCoreFirstDerivatives[i][j][k][l][dimA1]
                                         *rotMatFirstDerivatives[mu    ][i][dimA2]
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotatingMatrix        [lambda][k] 
                                         *rotatingMatrix        [sigma ][l];
                                 // term6
                                 value += diatomicTwoElecTwoCoreFirstDerivatives[i][j][k][l][dimA1]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotMatFirstDerivatives[nu    ][j][dimA2]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotatingMatrix        [sigma ][l];
                                 // term7
                                 value += diatomicTwoElecTwoCoreFirstDerivatives[i][j][k][l][dimA1]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotMatFirstDerivatives[lambda][k][dimA2]
                                         *rotatingMatrix        [sigma ][l];
                                 // term8
                                 value += diatomicTwoElecTwoCoreFirstDerivatives[i][j][k][l][dimA1]
                                         *rotatingMatrix        [mu    ][i]
                                         *rotatingMatrix        [nu    ][j]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotMatFirstDerivatives[sigma ][l][dimA2];

                                 // term9
                                 value += diatomicTwoElecTwoCoreFirstDerivatives[i][j][k][l][dimA2]
                                         *rotMatFirstDerivatives[mu    ][i][dimA1]
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotatingMatrix        [lambda][k] 
                                         *rotatingMatrix        [sigma ][l];
                                 // term10
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMatFirstDerivatives[mu    ][i][dimA1]
                                         *rotMatFirstDerivatives[nu    ][j][dimA2]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotatingMatrix        [sigma ][l];
                                 // term11
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMatFirstDerivatives[mu    ][i][dimA1]
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotMatFirstDerivatives[lambda][k][dimA2]
                                         *rotatingMatrix        [sigma ][l];
                                 // term12
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMatFirstDerivatives[mu    ][i][dimA1]
                                         *rotatingMatrix        [nu    ][j]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotMatFirstDerivatives[sigma ][l][dimA2];

                                 // term13
                                 value += diatomicTwoElecTwoCoreFirstDerivatives[i][j][k][l][dimA2]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotMatFirstDerivatives[nu    ][j][dimA1]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotatingMatrix        [sigma ][l];
                                 // term14
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMatFirstDerivatives[mu    ][i][dimA2]
                                         *rotMatFirstDerivatives[nu    ][j][dimA1]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotatingMatrix        [sigma ][l];
                                 // term15
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotMatFirstDerivatives[nu    ][j][dimA1]
                                         *rotMatFirstDerivatives[lambda][k][dimA2]
                                         *rotatingMatrix        [sigma ][l];
                                 // term16
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotMatFirstDerivatives[nu    ][j][dimA1]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotMatFirstDerivatives[sigma ][l][dimA2];

                                 // term17
                                 value += diatomicTwoElecTwoCoreFirstDerivatives[i][j][k][l][dimA2]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotMatFirstDerivatives[lambda][k][dimA1]
                                         *rotatingMatrix        [sigma ][l];
                                 // term18
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMatFirstDerivatives[mu    ][i][dimA2]
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotMatFirstDerivatives[lambda][k][dimA1]
                                         *rotatingMatrix        [sigma ][l];
                                 // term19
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotMatFirstDerivatives[nu    ][j][dimA2]
                                         *rotMatFirstDerivatives[lambda][k][dimA1]
                                         *rotatingMatrix        [sigma ][l];
                                 // term20
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotMatFirstDerivatives[lambda][k][dimA1]
                                         *rotMatFirstDerivatives[sigma ][l][dimA2];
                                 
                                 // term21
                                 value += diatomicTwoElecTwoCoreFirstDerivatives[i][j][k][l][dimA2]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotatingMatrix        [lambda][k] 
                                         *rotMatFirstDerivatives[sigma ][l][dimA1];
                                 // term22
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotMatFirstDerivatives[mu    ][i][dimA2]
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotatingMatrix        [lambda][k] 
                                         *rotMatFirstDerivatives[sigma ][l][dimA1];
                                 // term23
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotMatFirstDerivatives[nu    ][j][dimA2]
                                         *rotatingMatrix        [lambda][k] 
                                         *rotMatFirstDerivatives[sigma ][l][dimA1];
                                 // term24
                                 value += diatomicTwoElecTwoCore[i][j][k][l]
                                         *rotatingMatrix        [mu    ][i] 
                                         *rotatingMatrix        [nu    ][j] 
                                         *rotMatFirstDerivatives[lambda][k][dimA2]
                                         *rotMatFirstDerivatives[sigma ][l][dimA1];
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
   double Rab = this->molecule->GetDistanceAtoms(atomA, atomB);
   int lA = 0;
   int lB = 0;
   // (28) in [DT_1977]
   if(mu == s && nu == s && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      value = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
   }
   // (29) in [DT_1977]
   else if(mu == s && nu == s && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == s && nu == s && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   // (30) in [DT_1977]
   else if(mu == s && nu == s && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   // (31) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == py && nu == py && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   // (32) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   // (33) in [DT_1977]
   else if(mu == px && nu == px && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (34) in [DT_1977]
   else if(mu == px && nu == px && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (35) in [DT_1977]
   else if(mu == px && nu == px && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (36) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == pz && nu == pz && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (37) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   // (38) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   // (39) in [DT_1977]
   else if(mu == s && nu == pz && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(muz, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == px && sigma == px){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   else if(mu == s && nu == pz && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(muz, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == py && sigma == py){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   // (40) in [DT_1977]
   else if(mu == s && nu == pz && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(muz, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == pz){
      value = this->GetNddoRepulsionIntegral(atomA, nu, mu, atomB, lambda, sigma);
   }
   // (41) in [DT_1977]
   else if(mu == s && nu == s && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
   }
   else if(mu == s && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   // (42) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qxx, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == px && nu == px && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   else if(mu == py && nu == py && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qyy, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   // (43) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteraction(Qzz, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegral(atomA, mu, nu, atomB, sigma, lambda);
   }
   // (44) in [DT_1977]
   else if(mu == s && nu == px && lambda == s && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(mux, mux, rhoA, rhoB, DA, DB, Rab);
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
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muy, muy, rhoA, rhoB, DA, DB, Rab);
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
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muz, muz, rhoA, rhoB, DA, DB, Rab);
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
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(mux, Qxz, rhoA, rhoB, DA, DB, Rab);
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
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(muy, Qyz, rhoA, rhoB, DA, DB, Rab);
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
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(Qxz, mux, rhoA, rhoB, DA, DB, Rab);
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
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(Qyz, muy, rhoA, rhoB, DA, DB, Rab);
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
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(Qxz, Qxz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
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
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteraction(Qyz, Qyz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
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
double Mndo::GetNddoRepulsionIntegralFirstDerivative(
                                       const Atom& atomA, OrbitalType mu, OrbitalType nu,
                                       const Atom& atomB, OrbitalType lambda, OrbitalType sigma,
                                       CartesianType axisA) const{
   double value = 0.0;
   double DA=0.0;
   double DB=0.0;
   double rhoA = 0.0;
   double rhoB = 0.0;
   double Rab = this->molecule->GetDistanceAtoms(atomA, atomB);
   double dRabDa = (atomA.GetXyz()[axisA] - atomB.GetXyz()[axisA])/Rab;
   int lA = 0;
   int lB = 0;
   // (28) in [DT_1977]
   if(mu == s && nu == s && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      value *= dRabDa;
   }
   // (29) in [DT_1977]
   else if(mu == s && nu == s && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == s && nu == s && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   // (30) in [DT_1977]
   else if(mu == s && nu == s && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   // (31) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == py && nu == py && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   // (32) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   // (33) in [DT_1977]
   else if(mu == px && nu == px && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   else if(mu == py && nu == py && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   // (34) in [DT_1977]
   else if(mu == px && nu == px && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   else if(mu == py && nu == py && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   // (35) in [DT_1977]
   else if(mu == px && nu == px && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   // (36) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   else if(mu == pz && nu == pz && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   // (37) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp3 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp4 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2 + temp3 + temp4;
      value *= dRabDa;
   }
   // (38) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, sQ, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   // (39) in [DT_1977]
   else if(mu == s && nu == pz && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, Qxx, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == pz && nu == s && lambda == px && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == pz && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, Qyy, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == pz && nu == s && lambda == py && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   // (40) in [DT_1977]
   else if(mu == s && nu == pz && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, sQ, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, Qzz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == pz){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   // (41) in [DT_1977]
   else if(mu == s && nu == s && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == s && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   // (42) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxx, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == px && nu == px && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == py && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyy, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   // (43) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           sQ, muz, rhoA, rhoB, DA, DB, Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp2 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qzz, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1 + temp2;
      value *= dRabDa;
   }
   else if(mu == pz && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   // (44) in [DT_1977]
   else if(mu == s && nu == px && lambda == s && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           mux, mux, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == px && nu == s && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == px && nu == s && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   else if(mu == s && nu == py && lambda == s && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muy, muy, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == py && nu == s && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == s && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (45) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muz, muz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == pz){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (46) in [DT_1977]
   else if(mu == s && nu == px && lambda == px && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           mux, Qxz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == px && nu == s && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == px && nu == s && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   else if(mu == s && nu == py && lambda == py && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           muy, Qyz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == py && nu == s && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == s && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == s && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (47) in [DT_1977]
   else if(mu == px && nu == pz && lambda == s && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxz, mux, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == pz && nu == px && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == px && nu == pz && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == pz && lambda == s && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyz, muy, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == pz && nu == py && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (48) in [DT_1977]
   else if(mu == px && nu == pz && lambda == px && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qxz, Qxz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == pz && nu == px && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == px && nu == pz && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      double temp1 = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                           Qyz, Qyz, rhoA, rhoB, DA, DB, Rab);
      value = temp1;
      value *= dRabDa;
   }
   else if(mu == pz && nu == py && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == py && nu == pz && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == pz && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // (49) in [DT_1977] and p19 in [MOPAC_1990]
   else if(mu == px && nu == py && lambda == px && sigma == py){
      value = 0.5*(this->GetNddoRepulsionIntegralFirstDerivative(
                         atomA, mu, mu, atomB, mu, mu, axisA)
                  -this->GetNddoRepulsionIntegralFirstDerivative(
                         atomA, mu, mu, atomB, nu, nu, axisA));
   }
   else if(mu == py && nu == px && lambda == px && sigma == py){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA);
   }
   else if(mu == px && nu == py && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA);
   }
   else if(mu == py && nu == px && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegralFirstDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA);
   }
   // d-orbitals
   else if(mu == dxy || mu == dyz || mu == dzz || mu == dzx || mu == dxxyy ||
           nu == dxy || nu == dyz || nu == dzz || nu == dzx || nu == dxxyy ||
           lambda == dxy || lambda == dyz || lambda == dzz || lambda  == dzx || lambda == dxxyy ||
           sigma == dxy || sigma == dyz || sigma == dzz || sigma  == dzx || sigma == dxxyy){

      stringstream ss;
      ss << this->errorMessageGetNddoRepulsionIntegralFirstDerivative;
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
double Mndo::GetNddoRepulsionIntegralSecondDerivative(
                                       const Atom& atomA, OrbitalType mu, OrbitalType nu,
                                       const Atom& atomB, OrbitalType lambda, OrbitalType sigma,
                                       CartesianType axisA1,
                                       CartesianType axisA2) const{
   double value = 0.0;
   double DA=0.0;
   double DB=0.0;
   double rhoA = 0.0;
   double rhoB = 0.0;
   double Rab = this->molecule->GetDistanceAtoms(atomA, atomB);
   double cartesian[CartesianType_end] = {atomA.GetXyz()[XAxis] - atomB.GetXyz()[XAxis], 
                                          atomA.GetXyz()[YAxis] - atomB.GetXyz()[YAxis],
                                          atomA.GetXyz()[ZAxis] - atomB.GetXyz()[ZAxis]};
   double firstDeri=0.0;
   double secondDeri=0.0;
   int lA = 0;
   int lB = 0;
   // (28) in [DT_1977]
   if(mu == s && nu == s && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      value = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                      secondDeri,
                                                                      axisA1,
                                                                      axisA2,
                                                                      cartesian,
                                                                      Rab);
   }
   // (29) in [DT_1977]
   else if(mu == s && nu == s && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2;
   }
   else if(mu == s && nu == s && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2;
   }
   // (30) in [DT_1977]
   else if(mu == s && nu == s && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2;
   }
   // (31) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2;
   }
   else if(mu == py && nu == py && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2;
   }
   // (32) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2;
   }
   // (33) in [DT_1977]
   else if(mu == px && nu == px && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp3 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qxx, Qxx, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qxx, Qxx, rhoA, rhoB, DA, DB, Rab);
      double temp4 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp3 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qyy, Qyy, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qyy, Qyy, rhoA, rhoB, DA, DB, Rab);
      double temp4 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      value = temp1 + temp2 + temp3 + temp4;
   }
   // (34) in [DT_1977]
   else if(mu == px && nu == px && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp3 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qxx, Qyy, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qxx, Qyy, rhoA, rhoB, DA, DB, Rab);
      double temp4 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp3 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qyy, Qxx, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qyy, Qxx, rhoA, rhoB, DA, DB, Rab);
      double temp4 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      value = temp1 + temp2 + temp3 + temp4;
   }
   // (35) in [DT_1977]
   else if(mu == px && nu == px && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qxx, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp3 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qxx, Qzz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qxx, Qzz, rhoA, rhoB, DA, DB, Rab);
      double temp4 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qyy, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp3 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qyy, Qzz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qyy, Qzz, rhoA, rhoB, DA, DB, Rab);
      double temp4 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      value = temp1 + temp2 + temp3 + temp4;
   }
   // (36) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, Qxx, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp3 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qzz, Qxx, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qzz, Qxx, rhoA, rhoB, DA, DB, Rab);
      double temp4 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      value = temp1 + temp2 + temp3 + temp4;
   }
   else if(mu == pz && nu == pz && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, Qyy, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp3 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qzz, Qyy, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qzz, Qyy, rhoA, rhoB, DA, DB, Rab);
      double temp4 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      value = temp1 + temp2 + temp3 + temp4;
   }
   // (37) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, Qzz, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qzz, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp3 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qzz, Qzz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qzz, Qzz, rhoA, rhoB, DA, DB, Rab);
      double temp4 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      value = temp1 + temp2 + temp3 + temp4;
   }
   // (38) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == s){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        muz, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         muz, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   // (39) in [DT_1977]
   else if(mu == s && nu == pz && lambda == px && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        muz, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         muz, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        muz, Qxx, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         muz, Qxx, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == px && sigma == px){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == s && nu == pz && lambda == py && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        muz, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         muz, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        muz, Qyy, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         muz, Qyy, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == py && sigma == py){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   // (40) in [DT_1977]
   else if(mu == s && nu == pz && lambda == pz && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 0);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 0);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        muz, sQ, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         muz, sQ, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        muz, Qzz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         muz, Qzz, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == pz){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   // (41) in [DT_1977]
   else if(mu == s && nu == s && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, muz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, muz, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1;
   }
   else if(mu == s && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (42) in [DT_1977]
   else if(mu == px && nu == px && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, muz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, muz, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qxx, muz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qxx, muz, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2;
   }
   else if(mu == px && nu == px && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == py && nu == py && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, muz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, muz, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qyy, muz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qyy, muz, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2;
   }
   else if(mu == py && nu == py && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (43) in [DT_1977]
   else if(mu == pz && nu == pz && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 0);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 0);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        sQ, muz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         sQ, muz, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);

      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qzz, muz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qzz, muz, rhoA, rhoB, DA, DB, Rab);
      double temp2 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1 + temp2;
   }
   else if(mu == pz && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (44) in [DT_1977]
   else if(mu == s && nu == px && lambda == s && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        mux, mux, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         mux, mux, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1;
   }
   else if(mu == px && nu == s && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == s && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == px && nu == s && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == s && nu == py && lambda == s && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        muy, muy, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         muy, muy, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1;
   }
   else if(mu == py && nu == s && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == s && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == py && nu == s && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (45) in [DT_1977]
   else if(mu == s && nu == pz && lambda == s && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        muz, muz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         muz, muz, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1;
   }
   else if(mu == pz && nu == s && lambda == s && sigma == pz){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == s && nu == pz && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == pz && nu == s && lambda == pz && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (46) in [DT_1977]
   else if(mu == s && nu == px && lambda == px && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        mux, Qxz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         mux, Qxz, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1;
   }
   else if(mu == px && nu == s && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == s && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == px && nu == s && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == s && nu == py && lambda == py && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 1);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 1);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        muy, Qyz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         muy, Qyz, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1;
   }
   else if(mu == py && nu == s && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == s && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == py && nu == s && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (47) in [DT_1977]
   else if(mu == px && nu == pz && lambda == s && sigma == px){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qxz, mux, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qxz, mux, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1;
   }
   else if(mu == pz && nu == px && lambda == s && sigma == px){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == px && nu == pz && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == pz && nu == px && lambda == px && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == py && nu == pz && lambda == s && sigma == py){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 1);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 1);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qyz, muy, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qyz, muy, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1;
   }
   else if(mu == pz && nu == py && lambda == s && sigma == py){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == pz && nu == py && lambda == py && sigma == s){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (48) in [DT_1977]
   else if(mu == px && nu == pz && lambda == px && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qxz, Qxz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qxz, Qxz, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1;
   }
   else if(mu == pz && nu == px && lambda == px && sigma == pz){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == px && nu == pz && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == pz && nu == px && lambda == pz && sigma == px){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == py && nu == pz && lambda == py && sigma == pz){
      DA = atomA.GetNddoDerivedParameterD(this->theory, 2);
      DB = atomB.GetNddoDerivedParameterD(this->theory, 2);
      rhoA = atomA.GetNddoDerivedParameterRho(this->theory, 2);
      rhoB = atomB.GetNddoDerivedParameterRho(this->theory, 2);
      firstDeri = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                        Qyz, Qyz, rhoA, rhoB, DA, DB, Rab);
      secondDeri = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                         Qyz, Qyz, rhoA, rhoB, DA, DB, Rab);
      double temp1 = this->GetSecondDerivativeElementFromDistanceDerivatives(firstDeri,
                                                                             secondDeri,
                                                                             axisA1,
                                                                             axisA2,
                                                                             cartesian,
                                                                             Rab);
      value = temp1;
   }
   else if(mu == pz && nu == py && lambda == py && sigma == pz){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == py && nu == pz && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == pz && nu == py && lambda == pz && sigma == py){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // (49) in [DT_1977] and p19 in [MOPAC_1990]
   else if(mu == px && nu == py && lambda == px && sigma == py){
      value = 0.5*(this->GetNddoRepulsionIntegralSecondDerivative(
                         atomA, mu, mu, atomB, mu, mu, axisA1, axisA2)
                  -this->GetNddoRepulsionIntegralSecondDerivative(
                         atomA, mu, mu, atomB, nu, nu, axisA1, axisA2));
   }
   else if(mu == py && nu == px && lambda == px && sigma == py){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, lambda, sigma, axisA1, axisA2);
   }
   else if(mu == px && nu == py && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, mu, nu, atomB, sigma, lambda, axisA1, axisA2);
   }
   else if(mu == py && nu == px && lambda == py && sigma == px){
      value = this->GetNddoRepulsionIntegralSecondDerivative(
                    atomA, nu, mu, atomB, sigma, lambda, axisA1, axisA2);
   }
   // d-orbitals
   else if(mu == dxy || mu == dyz || mu == dzz || mu == dzx || mu == dxxyy ||
           nu == dxy || nu == dyz || nu == dzz || nu == dzx || nu == dxxyy ||
           lambda == dxy || lambda == dyz || lambda == dzz || lambda  == dzx || lambda == dxxyy ||
           sigma == dxy || sigma == dyz || sigma == dzz || sigma  == dzx || sigma == dxxyy){

      stringstream ss;
      ss << this->errorMessageGetNddoRepulsionIntegralSecondDerivative;
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
double Mndo::GetSemiEmpiricalMultipoleInteraction(MultipoleType multipoleA,
                                                  MultipoleType multipoleB,
                                                  double rhoA,
                                                  double rhoB,
                                                  double DA,
                                                  double DB,
                                                  double Rab) const{
   double value = 0.0;
   double a = rhoA + rhoB;

   // Eq. (52) in [DT_1977]
   if(multipoleA == sQ && multipoleB == sQ){
      value = pow(pow(Rab,2.0) + pow(a,2.0), -0.5);
   }
   // Eq. (53) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == muz){
      double temp1 = pow(Rab+DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/2.0 - pow(temp2,-0.5)/2.0;
   }
   else if(multipoleA == muz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,1.0);
   }
   // Eq. (54) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == Qxx){
      double temp1 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/2.0 - pow(temp2,-0.5)/2.0;
   }
   else if(multipoleA == Qxx && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == sQ && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleA, Qxx,
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyy && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   // Eq. (55) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == Qzz){
      double temp1 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 - pow(temp2,-0.5)/2.0 + pow(temp3,-0.5)/4.0;
   }
   else if(multipoleA == Qzz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   // Eq. (56) in [DT_1977]
   else if(multipoleA == mux && multipoleB == mux){
      double temp1 = pow(Rab,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/2.0 - pow(temp2,-0.5)/2.0;
   }
   else if(multipoleA == muy && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteraction(mux, mux, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   // Eq. (57) in [DT_1977]
   else if(multipoleA == muz && multipoleB == muz){
      double temp1 = pow(Rab+DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA+DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 - pow(temp2,-0.5)/4.0 
             -pow(temp3,-0.5)/4.0 + pow(temp4,-0.5)/4.0;
   }
   // Eq. (58) in [DT_1977]
   else if(multipoleA == mux && multipoleB == Qxz){
      double temp1 = pow(Rab-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value =-pow(temp1,-0.5)/4.0 + pow(temp2,-0.5)/4.0 
             +pow(temp3,-0.5)/4.0 - pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == Qxz && multipoleB == mux){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muy && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteraction(mux, Qxz, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyz && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   // Eq. (59) in [DT_1977]
   else if(multipoleA == muz && multipoleB == Qxx){
      double temp1 = pow(Rab+DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA,2.0) + pow(a,2.0);
      value =-pow(temp1,-0.5)/4.0 + pow(temp2,-0.5)/4.0 
             +pow(temp3,-0.5)/4.0 - pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == Qxx && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(muz, Qxx, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyy && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   // Eq. (60) in [DT_1977]
   else if(multipoleA == muz && multipoleB == Qzz){
      double temp1 = pow(Rab+DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DA-2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA+2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab+DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-DA,2.0) + pow(a,2.0);
      value =-pow(temp1,-0.5)/8.0 + pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/8.0 + pow(temp4,-0.5)/8.0
             +pow(temp5,-0.5)/4.0 - pow(temp6,-0.5)/4.0;
   }
   else if(multipoleA == Qzz && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   // Eq. (61) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qxx){
      double temp1 = pow(Rab,2.0) + 4.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + 4.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp4 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/8.0 + pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/4.0 - pow(temp4,-0.5)/4.0
             +pow(temp5,-0.5)/4.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(Qxx, Qxx, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   // Eq. (62) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qyy){
      double temp1 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(2.0*DB,2.0)+ pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/4.0 - pow(temp2,-0.5)/4.0 
             -pow(temp3,-0.5)/4.0 + pow(temp4,-0.5)/4.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   // Eq. (63) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qzz){
      double temp1 = pow(Rab-2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/8.0 + pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/8.0 - pow(temp4,-0.5)/8.0
             -pow(temp5,-0.5)/4.0 + pow(temp6,-0.5)/4.0;
   }
   else if(multipoleA == Qzz && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qyy && multipoleB == Qzz){
      value = this->GetSemiEmpiricalMultipoleInteraction(Qxx, multipoleB, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qzz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteraction(multipoleB, multipoleA,
                                                         rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   // Eq. (64) in [DT_1977]
   else if(multipoleA == Qzz && multipoleB == Qzz){
      double temp1 = pow(Rab+2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab+2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-2.0*DA,2.0) + pow(a,2.0);
      double temp7 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp8 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      double temp9 = pow(Rab,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/16.0 + pow(temp2,-0.5)/16.0 
             +pow(temp3,-0.5)/16.0 + pow(temp4,-0.5)/16.0
             -pow(temp5,-0.5)/8.0 - pow(temp6,-0.5)/8.0
             -pow(temp7,-0.5)/8.0 - pow(temp8,-0.5)/8.0
             +pow(temp9,-0.5)/4.0;
   }
   // Eq. (65) in [DT_1977]
   else if(multipoleA == Qxz && multipoleB == Qxz){
      double temp1 = pow(Rab+DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab-DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp7 = pow(Rab-DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp8 = pow(Rab-DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = pow(temp1,-0.5)/8.0 - pow(temp2,-0.5)/8.0 
             -pow(temp3,-0.5)/8.0 + pow(temp4,-0.5)/8.0
             -pow(temp5,-0.5)/8.0 + pow(temp6,-0.5)/8.0
             +pow(temp7,-0.5)/8.0 - pow(temp8,-0.5)/8.0;
   }
   else if(multipoleA == Qyz && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteraction(Qxz, Qxz, 
                                                         rhoA, rhoB, DA, DB, Rab);
   }
   // Eq. (66) in [DT_1977]
   else if(multipoleA == Qxy && multipoleB == Qxy){
      double temp1 = pow(Rab,2.0) + 2.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + 2.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + 2.0*pow(DA,2.0) + 2.0*pow(DB,2.0) + pow(a,2.0);
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
// This derivativ is related to the nuclear distance Rab.
// See Apendix in [DT_1977]
double Mndo::GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                                                  MultipoleType multipoleA,
                                                  MultipoleType multipoleB,
                                                  double rhoA,
                                                  double rhoB,
                                                  double DA,
                                                  double DB,
                                                  double Rab) const{
   double value = 0.0;
   double a = rhoA + rhoB;

   // Eq. (52) in [DT_1977]
   if(multipoleA == sQ && multipoleB == sQ){
      value = -1.0*Rab*pow(pow(Rab,2.0) + pow(a,2.0), -1.5);
   }
   // Eq. (53) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == muz){
      double temp1 = pow(Rab+DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DB,2.0) + pow(a,2.0);
      value = (Rab+DB)*pow(temp1,-1.5)/2.0 
             -(Rab-DB)*pow(temp2,-1.5)/2.0;
      value *= -1.0;
   }
   else if(multipoleA == muz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,1.0);
   }
   // Eq. (54) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == Qxx){
      double temp1 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(a,2.0);
      value = Rab*pow(temp1,-1.5)/2.0 
             -Rab*pow(temp2,-1.5)/2.0;
      value *= -1.0;
   }
   else if(multipoleA == Qxx && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == sQ && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleA, Qxx, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyy && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   // Eq. (55) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == Qzz){
      double temp1 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      value = (Rab+2.0*DB)*pow(temp1,-1.5)/4.0 
             -(Rab)*pow(temp2,-1.5)/2.0 
             +(Rab-2.0*DB)*pow(temp3,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qzz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   // Eq. (56) in [DT_1977]
   else if(multipoleA == mux && multipoleB == mux){
      double temp1 = pow(Rab,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = (Rab)*pow(temp1,-1.5)/2.0 
             -(Rab)*pow(temp2,-1.5)/2.0;
      value *= -1.0;
   }
   else if(multipoleA == muy && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    mux, mux, rhoA, rhoB, DA, DB, Rab);
   }
   // Eq. (57) in [DT_1977]
   else if(multipoleA == muz && multipoleB == muz){
      double temp1 = pow(Rab+DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA+DB,2.0) + pow(a,2.0);
      value = (Rab+DA-DB)*pow(temp1,-1.5)/4.0 
             -(Rab+DA+DB)*pow(temp2,-1.5)/4.0 
             -(Rab-DA-DB)*pow(temp3,-1.5)/4.0 
             +(Rab-DA+DB)*pow(temp4,-1.5)/4.0;
      value *= -1.0;
   }
   // Eq. (58) in [DT_1977]
   else if(multipoleA == mux && multipoleB == Qxz){
      double temp1 = pow(Rab-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value =-(Rab-DB)*pow(temp1,-1.5)/4.0 
             +(Rab-DB)*pow(temp2,-1.5)/4.0 
             +(Rab+DB)*pow(temp3,-1.5)/4.0 
             -(Rab+DB)*pow(temp4,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qxz && multipoleB == mux){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muy && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    mux, Qxz, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyz && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   // Eq. (59) in [DT_1977]
   else if(multipoleA == muz && multipoleB == Qxx){
      double temp1 = pow(Rab+DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DA,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA,2.0) + pow(a,2.0);
      value =-(Rab+DA)*pow(temp1,-1.5)/4.0 
             +(Rab-DA)*pow(temp2,-1.5)/4.0 
             +(Rab+DA)*pow(temp3,-1.5)/4.0 
             -(Rab-DA)*pow(temp4,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qxx && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    muz, Qxx, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyy && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   // Eq. (60) in [DT_1977]
   else if(multipoleA == muz && multipoleB == Qzz){
      double temp1 = pow(Rab+DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab-DA-2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA+2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab+DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-DA,2.0) + pow(a,2.0);
      value =-(Rab+DA-2.0*DB)*pow(temp1,-1.5)/8.0 
             +(Rab-DA-2.0*DB)*pow(temp2,-1.5)/8.0 
             -(Rab+DA+2.0*DB)*pow(temp3,-1.5)/8.0 
             +(Rab-DA+2.0*DB)*pow(temp4,-1.5)/8.0
             +(Rab+DA       )*pow(temp5,-1.5)/4.0 
             -(Rab-DA       )*pow(temp6,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qzz && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   // Eq. (61) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qxx){
      double temp1 = pow(Rab,2.0) + 4.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + 4.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp4 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab,2.0) + pow(a,2.0);
      value = (Rab)*pow(temp1,-1.5)/8.0 
             +(Rab)*pow(temp2,-1.5)/8.0 
             -(Rab)*pow(temp3,-1.5)/4.0 
             -(Rab)*pow(temp4,-1.5)/4.0
             +(Rab)*pow(temp5,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    Qxx, Qxx, rhoA, rhoB, DA, DB, Rab);
   }
   // Eq. (62) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qyy){
      double temp1 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(2.0*DB,2.0)+ pow(a,2.0);
      double temp2 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + pow(2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab,2.0) + pow(a,2.0);
      value = (Rab)*pow(temp1,-1.5)/4.0 
             -(Rab)*pow(temp2,-1.5)/4.0 
             -(Rab)*pow(temp3,-1.5)/4.0 
             +(Rab)*pow(temp4,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qyy && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   // Eq. (63) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qzz){
      double temp1 = pow(Rab-2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+2.0*DB,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab,2.0) + pow(2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab,2.0) + pow(a,2.0);
      value = (Rab-2.0*DB)*pow(temp1,-1.5)/8.0 
             +(Rab+2.0*DB)*pow(temp2,-1.5)/8.0 
             -(Rab-2.0*DB)*pow(temp3,-1.5)/8.0 
             -(Rab+2.0*DB)*pow(temp4,-1.5)/8.0
             -(Rab       )*pow(temp5,-1.5)/4.0 
             +(Rab       )*pow(temp6,-1.5)/4.0;
      value *= -1.0;
   }
   else if(multipoleA == Qzz && multipoleB == Qxx){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qyy && multipoleB == Qzz){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    Qxx, multipoleB, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qzz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   // Eq. (64) in [DT_1977]
   else if(multipoleA == Qzz && multipoleB == Qzz){
      double temp1 = pow(Rab+2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab-2.0*DA-2.0*DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab-2.0*DA+2.0*DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab+2.0*DA,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-2.0*DA,2.0) + pow(a,2.0);
      double temp7 = pow(Rab+2.0*DB,2.0) + pow(a,2.0);
      double temp8 = pow(Rab-2.0*DB,2.0) + pow(a,2.0);
      double temp9 = pow(Rab,2.0) + pow(a,2.0);
      value = (Rab+2.0*DA-2.0*DB)*pow(temp1,-1.5)/16.0 
             +(Rab+2.0*DA+2.0*DB)*pow(temp2,-1.5)/16.0 
             +(Rab-2.0*DA-2.0*DB)*pow(temp3,-1.5)/16.0 
             +(Rab-2.0*DA+2.0*DB)*pow(temp4,-1.5)/16.0
             -(Rab+2.0*DA)*pow(temp5,-1.5)/8.0 
             -(Rab-2.0*DA)*pow(temp6,-1.5)/8.0
             -(Rab+2.0*DB)*pow(temp7,-1.5)/8.0 
             -(Rab-2.0*DB)*pow(temp8,-1.5)/8.0
             +(Rab)*pow(temp9,-1.5)/4.0;
      value *= -1.0;
   }
   // Eq. (65) in [DT_1977]
   else if(multipoleA == Qxz && multipoleB == Qxz){
      double temp1 = pow(Rab+DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab+DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab+DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp4 = pow(Rab+DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp5 = pow(Rab-DA-DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp6 = pow(Rab-DA-DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      double temp7 = pow(Rab-DA+DB,2.0) + pow(DA-DB,2.0) + pow(a,2.0);
      double temp8 = pow(Rab-DA+DB,2.0) + pow(DA+DB,2.0) + pow(a,2.0);
      value = (Rab+DA-DB)*pow(temp1,-1.5)/8.0 
             -(Rab+DA-DB)*pow(temp2,-1.5)/8.0 
             -(Rab+DA+DB)*pow(temp3,-1.5)/8.0 
             +(Rab+DA+DB)*pow(temp4,-1.5)/8.0
             -(Rab-DA-DB)*pow(temp5,-1.5)/8.0 
             +(Rab-DA-DB)*pow(temp6,-1.5)/8.0
             +(Rab-DA+DB)*pow(temp7,-1.5)/8.0 
             -(Rab-DA+DB)*pow(temp8,-1.5)/8.0;
      value *= -1.0;
   }
   else if(multipoleA == Qyz && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteractionFirstDerivative(
                    Qxz, Qxz, rhoA, rhoB, DA, DB, Rab);
   }
   // Eq. (66) in [DT_1977]
   else if(multipoleA == Qxy && multipoleB == Qxy){
      double temp1 = pow(Rab,2.0) + 2.0*pow(DA-DB,2.0) + pow(a,2.0);
      double temp2 = pow(Rab,2.0) + 2.0*pow(DA+DB,2.0) + pow(a,2.0);
      double temp3 = pow(Rab,2.0) + 2.0*pow(DA,2.0) + 2.0*pow(DB,2.0) + pow(a,2.0);
      value = (Rab)*pow(temp1,-1.5)/4.0 
             +(Rab)*pow(temp2,-1.5)/4.0 
             -(Rab)*pow(temp3,-1.5)/2.0;
      value *= -1.0;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetSemiEmpiricalMultipoleInteractionFirstDeriBadMultipoles;
      ss << this->errorMessageMultipoleA << MultipoleTypeStr(multipoleA) << endl;
      ss << this->errorMessageMultipoleB << MultipoleTypeStr(multipoleB) << endl;
      throw MolDSException(ss.str());
   }
   return value;
}

// Second derivative of semiempirical multipole-multipole interactions.
// This derivativ is related to the nuclear distance Rab.
// See Apendix in [DT_1977]
double Mndo::GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                                                  MultipoleType multipoleA,
                                                  MultipoleType multipoleB,
                                                  double rhoA,
                                                  double rhoB,
                                                  double DA,
                                                  double DB,
                                                  double Rab) const{
   double value = 0.0;
   double a = rhoA + rhoB;

   // Eq. (52) in [DT_1977]
   if(multipoleA == sQ && multipoleB == sQ){
      double c1 = 1.0;
      double f1 = pow(Rab,2.0);
      double a1 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
   }
   // Eq. (53) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == muz){
      double c1 = 0.5;
      double c2 = -0.5;
      double f1 = pow(Rab+DB,2.0);
      double f2 = pow(Rab-DB,2.0);
      double a1 = pow(a,2.0);
      double a2 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
   }
   else if(multipoleA == muz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,1.0);
   }
   // Eq. (54) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == Qxx){
      double c1 = 0.5;
      double c2 = -0.5;
      double f1 = pow(Rab,2.0);
      double f2 = pow(Rab,2.0);
      double a1 = pow(2.0*DB,2.0) + pow(a,2.0);
      double a2 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
   }
   else if(multipoleA == Qxx && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   else if(multipoleA == sQ && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleA, Qxx, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyy && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   // Eq. (55) in [DT_1977]
   else if(multipoleA == sQ && multipoleB == Qzz){
      double c1 = 0.25;
      double c2 = -0.50;
      double c3 = 0.25;
      double f1 = pow(Rab+2.0*DB,2.0);
      double f2 = pow(Rab,2.0);
      double f3 = pow(Rab-2.0*DB,2.0);
      double a1 = pow(a,2.0);
      double a2 = pow(a,2.0);
      double a3 = pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
   }
   else if(multipoleA == Qzz && multipoleB == sQ){
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,2.0);
   }
   // Eq. (56) in [DT_1977]
   else if(multipoleA == mux && multipoleB == mux){
      double c1 = 0.50;
      double c2 = -0.50;
      double f1 = pow(Rab,2.0);
      double f2 = pow(Rab,2.0);
      double a1 = pow(DA-DB,2.0) + pow(a,2.0);
      double a2 = pow(DA+DB,2.0) + pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
   }
   else if(multipoleA == muy && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    mux, mux, rhoA, rhoB, DA, DB, Rab);
   }
   // Eq. (57) in [DT_1977]
   else if(multipoleA == muz && multipoleB == muz){
      double c1 =  0.25;
      double c2 = -0.25;
      double c3 = -0.25;
      double c4 =  0.25;
      double f1 = pow(Rab+DA-DB,2.0);
      double f2 = pow(Rab+DA+DB,2.0);
      double f3 = pow(Rab-DA-DB,2.0);
      double f4 = pow(Rab-DA+DB,2.0);
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
      double f1 = pow(Rab-DB,2.0);
      double f2 = pow(Rab-DB,2.0);
      double f3 = pow(Rab+DB,2.0);
      double f4 = pow(Rab+DB,2.0);
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
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muy && multipoleB == Qyz){
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    mux, Qxz, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyz && multipoleB == muy){
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   // Eq. (59) in [DT_1977]
   else if(multipoleA == muz && multipoleB == Qxx){
      double c1 = -0.25;
      double c2 = 0.25;
      double c3 = 0.25;
      double c4 = -0.25;
      double f1 = pow(Rab+DA,2.0);
      double f2 = pow(Rab-DA,2.0);
      double f3 = pow(Rab+DA,2.0);
      double f4 = pow(Rab-DA,2.0);
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
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   else if(multipoleA == muz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    muz, Qxx, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qyy && multipoleB == muz){
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
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
      double f1 = pow(Rab+DA-2.0*DB,2.0);
      double f2 = pow(Rab-DA-2.0*DB,2.0);
      double f3 = pow(Rab+DA+2.0*DB,2.0);
      double f4 = pow(Rab-DA+2.0*DB,2.0);
      double f5 = pow(Rab+DA,2.0);
      double f6 = pow(Rab-DA,2.0);
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
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,3.0);
   }
   // Eq. (61) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qxx){
      double c1 =  0.125;
      double c2 =  0.125;
      double c3 = -0.25;
      double c4 = -0.25;
      double c5 =  0.25;
      double f1 = pow(Rab,2.0);
      double f2 = pow(Rab,2.0);
      double f3 = pow(Rab,2.0);
      double f4 = pow(Rab,2.0);
      double f5 = pow(Rab,2.0);
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
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    Qxx, Qxx, rhoA, rhoB, DA, DB, Rab);
   }
   // Eq. (62) in [DT_1977]
   else if(multipoleA == Qxx && multipoleB == Qyy){
      double c1 =  0.25;
      double c2 = -0.25;
      double c3 = -0.25;
      double c4 =  0.25;
      double f1 = pow(Rab,2.0);
      double f2 = pow(Rab,2.0);
      double f3 = pow(Rab,2.0);
      double f4 = pow(Rab,2.0);
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
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
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
      double f1 = pow(Rab-2.0*DB,2.0);
      double f2 = pow(Rab+2.0*DB,2.0);
      double f3 = pow(Rab-2.0*DB,2.0);
      double f4 = pow(Rab+2.0*DB,2.0);
      double f5 = pow(Rab       ,2.0);
      double f6 = pow(Rab       ,2.0);
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
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
      value *= pow(-1.0,4.0);
   }
   else if(multipoleA == Qyy && multipoleB == Qzz){
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    Qxx, multipoleB, rhoA, rhoB, DA, DB, Rab);
   }
   else if(multipoleA == Qzz && multipoleB == Qyy){
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    multipoleB, multipoleA, rhoB, rhoA, DB, DA, Rab);
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
      double f1 = pow(Rab+2.0*DA-2.0*DB,2.0);
      double f2 = pow(Rab+2.0*DA+2.0*DB,2.0);
      double f3 = pow(Rab-2.0*DA-2.0*DB,2.0);
      double f4 = pow(Rab-2.0*DA+2.0*DB,2.0);
      double f5 = pow(Rab+2.0*DA       ,2.0);
      double f6 = pow(Rab-2.0*DA       ,2.0);
      double f7 = pow(Rab+2.0*DB       ,2.0);
      double f8 = pow(Rab-2.0*DB       ,2.0);
      double f9 = pow(Rab              ,2.0);
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
      double f1 = pow(Rab+DA-DB,2.0);
      double f2 = pow(Rab+DA-DB,2.0);
      double f3 = pow(Rab+DA+DB,2.0);
      double f4 = pow(Rab+DA+DB,2.0);
      double f5 = pow(Rab-DA-DB,2.0);
      double f6 = pow(Rab-DA-DB,2.0);
      double f7 = pow(Rab-DA+DB,2.0);
      double f8 = pow(Rab-DA+DB,2.0);
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
      value = this->GetSemiEmpiricalMultipoleInteractionSecondDerivative(
                    Qxz, Qxz, rhoA, rhoB, DA, DB, Rab);
   }
   // Eq. (66) in [DT_1977]
   else if(multipoleA == Qxy && multipoleB == Qxy){
      double c1 =  0.25;
      double c2 =  0.25;
      double c3 = -0.50;
      double f1 = pow(Rab,2.0);
      double f2 = pow(Rab,2.0);
      double f3 = pow(Rab,2.0);
      double a1 = 2.0*pow(DA-DB,2.0) + pow(a,2.0);
      double a2 = 2.0*pow(DA+DB,2.0) + pow(a,2.0);
      double a3 = 2.0*pow(DA,2.0) + 2.0*pow(DB,2.0) + pow(a,2.0);
      value  = c1*(3.0*f1*pow(f1+a1,-2.5) - pow(f1+a1,-1.5));
      value += c2*(3.0*f2*pow(f2+a2,-2.5) - pow(f2+a2,-1.5));
      value += c3*(3.0*f3*pow(f3+a3,-2.5) - pow(f3+a3,-1.5));
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetSemiEmpiricalMultipoleInteractionSecondDeriBadMultipoles;
      ss << this->errorMessageMultipoleA << MultipoleTypeStr(multipoleA) << endl;
      ss << this->errorMessageMultipoleB << MultipoleTypeStr(multipoleB) << endl;
      throw MolDSException(ss.str());
   }
   return value;
}

}


