//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   //
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
#include<boost/format.hpp>
#include"../base/Enums.h"
#include"../base/Uncopyable.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
#include"../base/EularAngle.h"
#include"../base/Parameters.h"
#include"../base/RealSphericalHarmonicsIndex.h"
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
#include"../mndo/Mndo.h"
#include"Am1.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_am1{

/***
 *  Main References for AM1 are [DZHS_1985, DY_1990]
 */
Am1::Am1() : MolDS_mndo::Mndo(){
   this->theory = AM1;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //this->OutputLog("Am1 created\n");
}

Am1::~Am1(){
   //this->OutputLog("Am1 deleted\n");
}

void Am1::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in am1::Am1::DoSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in am1::Am1::DoSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in am1::Am1::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in am1::Am1::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_am1::Am1::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_am1::Am1::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageCalcCISMatrix
      = "Error in am1::Am1::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in am1::Am1::DoCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles
      = "Error in am1::Am1::GetSemiEmpiricalMultipoleInteraction: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadAtomTypes
      = "Error in am1::Am1::GetSemiEmpiricalMultipoleInteraction: Bad atom types are set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteraction1stDeriBadMultipoles
      = "Error in am1::Am1::GetSemiEmpiricalMultipoleInteraction1stDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteraction2ndDeriBadMultipoles
      = "Error in am1::Am1::GetSemiEmpiricalMultipoleInteraction2ndDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetNddoRepulsionIntegral 
      = "Error in am1::Am1::GetNddoRepulsionIntegral: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegralBadAtomTypes
      = "Error in am1::Am1::GetNddoRepulsionIntegral: Bad atom types are set.\n";
   this->errorMessageGetNddoRepulsionIntegral1stDerivative 
      = "Error in am1::Am1::GetNddoRepulsionIntegral1stDerivative: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegral2ndDerivative 
      = "Error in am1::Am1::GetNddoRepulsionIntegral2ndDerivative: Bad orbital is set.\n";
   this->errorMessageCalcTwoElecsTwoAtomCoresNullMatrix 
      = "Error in am1::Am1::CalcTwoElecsTwoAtomCores: The two elec two atom core matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecsAtomEpcCoresNullMatrix 
      = "Error in am1::Am1::CalcTwoElecsAtomEpcCores: The two elec atom-epc core matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCoresSameAtoms
      = "Error in am1::Am1::CalcDiatomicTwoElecsTwoCores: Atom A and B is same atom (not EPC).\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCoresSameEpcs
      = "Error in am1::Am1::CalcDiatomicTwoElecsTwoCores: Atom A and B is same EPC.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores1stDerivativesSameAtoms
      = "Error in am1::Am1::CalcDiatomicTwoElecsTwoCores1stDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores2ndDerivativesSameAtoms
      = "Error in am1::Am1::CalcDiatomicTwoElecsTwoCores2ndDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCoresNullMatrix 
      = "Error in am1::Am1::CalcDiatomicTwoElecsTwoCores: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores1stDerivativesNullMatrix
      = "Error in am1::Am1::CalcDiatomicTwoElecsTwoCores1stDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores2ndDerivativesNullMatrix
      = "Error in am1::Am1::CalcDiatomicTwoElecsTwoCores2ndDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageGetElectronicEnergyEnergyNotCalculated
      = "Error in am1::Am1::GetElectronicEnergy: Set electronic state is not calculated by CIS.\n";
   this->errorMessageGetElectronicEnergyNULLCISEnergy 
      = "Error in am1::Am1::GetElectronicEnergy: excitedEnergies is NULL\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tAM1-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: AM1-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: AM1-SCF  **********\n\n\n";
   this->messageStartCIS = "**********  START: AM1-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: AM1-CIS  **********\n\n\n";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for AM1-CIS met convergence criterion(^^b\n\n\n";
}

void Am1::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
   this->enableAtomTypes.push_back(Zn);
}

double Am1::GetDiatomCoreRepulsionEnergy(const Atom& atomA, const Atom& atomB) const{
   // MNDO term
   double mndoTerm = Mndo::GetDiatomCoreRepulsionEnergy(atomA, atomB);

   // additional term, Eq. (4) in [S_1989].
   double distance   = this->molecule->GetDistanceAtoms(atomA, atomB);
   double ang2AU     = Parameters::GetInstance()->GetAngstrom2AU();
   double alphaA     = atomA.GetNddoAlpha(this->theory);
   double alphaB     = atomB.GetNddoAlpha(this->theory);
   double kA, lA, mA;
   double kB, lB, mB;
   double temp = 0.0;
   for(int i=0; i<4; i++){
      kA = atomA.GetNddoParameterK(this->theory, i);
      lA = atomA.GetNddoParameterL(this->theory, i);
      mA = atomA.GetNddoParameterM(this->theory, i);
      kB = atomB.GetNddoParameterK(this->theory, i);
      lB = atomB.GetNddoParameterL(this->theory, i);
      mB = atomB.GetNddoParameterM(this->theory, i);
      temp += this->GetAdditionalDiatomCoreRepulsionTerm(kA, lA, mA, distance);
      temp += this->GetAdditionalDiatomCoreRepulsionTerm(kB, lB, mB, distance);
   }
   double additionalTerm = atomA.GetCoreCharge()*atomB.GetCoreCharge()*temp*ang2AU/distance;

   return mndoTerm + additionalTerm;
}

// First derivative of diatomic core repulsion energy.
// This derivative is related to the coordinate of atomA.
double Am1::GetDiatomCoreRepulsion1stDerivative(const Atom& atomA, const Atom& atomB, 
                                                CartesianType axisA) const{
   // MNDO term
   double mndoTerms = Mndo::GetDiatomCoreRepulsion1stDerivative(atomA, atomB, axisA);

   // additional term, first derivative of eq. (4) in [S_1989]
   double ang2AU     = Parameters::GetInstance()->GetAngstrom2AU();
   double alphaA     = atomA.GetNddoAlpha(this->theory);
   double alphaB     = atomB.GetNddoAlpha(this->theory);
   double distance   = this->molecule->GetDistanceAtoms(atomA, atomB);
   double dCartesian = (atomA.GetXyz()[axisA] - atomB.GetXyz()[axisA]);
   double kA, lA, mA;
   double kB, lB, mB;
   double temp1 = 0.0;
   double temp2 = 0.0;
   for(int i=0; i<4; i++){
      kA = atomA.GetNddoParameterK(this->theory, i);
      lA = atomA.GetNddoParameterL(this->theory, i);
      mA = atomA.GetNddoParameterM(this->theory, i);
      kB = atomB.GetNddoParameterK(this->theory, i);
      lB = atomB.GetNddoParameterL(this->theory, i);
      mB = atomB.GetNddoParameterM(this->theory, i);
      temp1 += this->GetAdditionalDiatomCoreRepulsionTerm(kA, lA, mA, distance);
      temp1 += this->GetAdditionalDiatomCoreRepulsionTerm(kB, lB, mB, distance);
      temp2 += this->GetAdditionalDiatomCoreRepulsionTerm1stDerivative(kA, lA, mA, distance);
      temp2 += this->GetAdditionalDiatomCoreRepulsionTerm1stDerivative(kB, lB, mB, distance);
   }
   double additionalTerm = 0.0;
   additionalTerm  = -temp1/pow(distance,3.0)
                     +temp2/pow(distance,2.0);
   additionalTerm *= dCartesian*atomA.GetCoreCharge()*atomB.GetCoreCharge()*ang2AU;
   
   return mndoTerms + additionalTerm;
}

// Second derivative of diatomic core repulsion energy.
// Both derivatives are related to the coordinate of atomA.
double Am1::GetDiatomCoreRepulsion2ndDerivative(const Atom& atomA,
                                                const Atom& atomB, 
                                                CartesianType axisA1,
                                                CartesianType axisA2) const{
   // MNDO term
   double mndoTerm = Mndo::GetDiatomCoreRepulsion2ndDerivative(atomA, atomB, axisA1, axisA2);

   // additional term, first derivative of eq. (4) in [S_1989]
   double ang2AU     = Parameters::GetInstance()->GetAngstrom2AU();
   double alphaA     = atomA.GetNddoAlpha(this->theory);
   double alphaB     = atomB.GetNddoAlpha(this->theory);
   double distance   = this->molecule->GetDistanceAtoms(atomA, atomB);
   double kA, lA, mA;
   double kB, lB, mB;
   double temp1 = 0.0;
   double temp2 = 0.0;
   double temp3 = 0.0;
   for(int i=0; i<4; i++){
      kA = atomA.GetNddoParameterK(this->theory, i);
      lA = atomA.GetNddoParameterL(this->theory, i);
      mA = atomA.GetNddoParameterM(this->theory, i);
      kB = atomB.GetNddoParameterK(this->theory, i);
      lB = atomB.GetNddoParameterL(this->theory, i);
      mB = atomB.GetNddoParameterM(this->theory, i);
      temp1 += this->GetAdditionalDiatomCoreRepulsionTerm(kA, lA, mA, distance);
      temp1 += this->GetAdditionalDiatomCoreRepulsionTerm(kB, lB, mB, distance);
      temp2 += this->GetAdditionalDiatomCoreRepulsionTerm1stDerivative(kA, lA, mA, distance);
      temp2 += this->GetAdditionalDiatomCoreRepulsionTerm1stDerivative(kB, lB, mB, distance);
      temp3 += this->GetAdditionalDiatomCoreRepulsionTerm2ndDerivative(kA, lA, mA, distance);
      temp3 += this->GetAdditionalDiatomCoreRepulsionTerm2ndDerivative(kB, lB, mB, distance);
   }
   double additionalTerm = 0.0;
   if(axisA1 != axisA2){
      double dCartesian1 = (atomA.GetXyz()[axisA1] - atomB.GetXyz()[axisA1]);
      double dCartesian2 = (atomA.GetXyz()[axisA2] - atomB.GetXyz()[axisA2]);
      additionalTerm = 3.0*dCartesian1*dCartesian2*temp1/pow(distance,5.0)
                      -3.0*dCartesian1*dCartesian2*temp2/pow(distance,4.0)
                      +1.0*dCartesian1*dCartesian2*temp3/pow(distance,3.0);
   }
   else{
      double dCartesian = (atomA.GetXyz()[axisA1] - atomB.GetXyz()[axisA1]);
      additionalTerm =-(1.0/pow(distance,3.0) - 3.0*pow(dCartesian,2.0)/pow(distance,5.0))*temp1
                      +(1.0/pow(distance,2.0) - 3.0*pow(dCartesian,2.0)/pow(distance,4.0))*temp2
                      +(                            pow(dCartesian,2.0)/pow(distance,3.0))*temp3;
   }
   additionalTerm *= atomA.GetCoreCharge()*atomB.GetCoreCharge()*ang2AU;

   return mndoTerm + additionalTerm;
}

void Am1::CalcSCFProperties(){
   MolDS_cndo::Cndo2::CalcSCFProperties();
}

void Am1::OutputSCFResults() const{
   MolDS_cndo::Cndo2::OutputSCFResults();
}

}



