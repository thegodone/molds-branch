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
#include<boost/format.hpp>
#include"../base/PrintController.h"
#include"../base/Uncopyable.h"
#include"../base/Enums.h"
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
#include"../mndo/Mndo.h"
#include"Am1.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_am1{

/***
 *  Main Refferences for AM1 are [DZHS_1985, DY_1990]
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
      = "Error in am1:: Am1::GetSemiEmpiricalMultipoleInteraction: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionFirstDeriBadMultipoles
      = "Error in am1:: Am1::GetSemiEmpiricalMultipoleInteractionFirstDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionSecondDeriBadMultipoles
      = "Error in am1:: Am1::GetSemiEmpiricalMultipoleInteractionSecondDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetNddoRepulsionIntegral 
      = "Error in am1::Am1::GetNddoRepulsionIntegral: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegralFirstDerivative 
      = "Error in am1::Am1::GetNddoRepulsionIntegralFirstDerivative: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegralSecondDerivative 
      = "Error in am1::Am1::GetNddoRepulsionIntegralSecondDerivative: Bad orbital is set.\n";
   this->errorMessageCalcTwoElecTwoCoreNullMatrix 
      = "Error in am1::Am1::CalcTwoElecTwoCore: The two elec two core matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicSameAtoms
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomic: Atom A and B is same.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesSameAtoms
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomicFirstDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicSecondDerivativesSameAtoms
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomicSecondDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicNullMatrix 
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomic: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicFirstDerivativesNullMatrix
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomicFirstDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecTwoCoreDiatomicSecondDerivativesNullMatrix
      = "Error in am1::Am1::CalcTwoElecTwoCoreDiatomicSecondDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
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
}

double Am1::GetDiatomCoreRepulsionEnergy(int indexAtomA, int indexAtomB) const{
   double energy = Mndo::GetDiatomCoreRepulsionEnergy(indexAtomA, indexAtomB);
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   double distance = this->molecule->GetDistanceAtoms(indexAtomA, indexAtomB);
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double alphaA = atomA.GetNddoAlpha(this->theory);
   double alphaB = atomB.GetNddoAlpha(this->theory);
   double temp = 0.0;
   for(int i=0; i<4; i++){
      double kA = atomA.GetNddoParameterK(this->theory, i);
      double lA = atomA.GetNddoParameterL(this->theory, i);
      double mA = atomA.GetNddoParameterM(this->theory, i);
      double kB = atomB.GetNddoParameterK(this->theory, i);
      double lB = atomB.GetNddoParameterL(this->theory, i);
      double mB = atomB.GetNddoParameterM(this->theory, i);
      temp += kA*exp(-lA*pow(distance-mA,2.0));
      temp += kB*exp(-lB*pow(distance-mB,2.0));
   }
   energy += atomA.GetCoreCharge()*atomB.GetCoreCharge()*temp/(distance/ang2AU);
   return energy;
}

// First derivative of diatomic core repulsion energy.
// This derivative is related to the coordinate of atomA.
double Am1::GetDiatomCoreRepulsionFirstDerivative(int atomAIndex,
                                                  int atomBIndex, 
                                                  CartesianType axisA) const{
   double value = Mndo::GetDiatomCoreRepulsionFirstDerivative(atomAIndex,
                                                              atomBIndex,
                                                              axisA);
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   const Atom& atomA = *this->molecule->GetAtom(atomAIndex);
   const Atom& atomB = *this->molecule->GetAtom(atomBIndex);
   double alphaA = atomA.GetNddoAlpha(this->theory);
   double alphaB = atomB.GetNddoAlpha(this->theory);
   double Rab = this->molecule->GetDistanceAtoms(atomAIndex, atomBIndex);
   double dRabDa = (atomA.GetXyz()[axisA] - atomB.GetXyz()[axisA])/Rab;
   double temp1 = 0.0;
   double temp2 = 0.0;
   for(int i=0; i<4; i++){
      double kA = atomA.GetNddoParameterK(this->theory, i);
      double lA = atomA.GetNddoParameterL(this->theory, i);
      double mA = atomA.GetNddoParameterM(this->theory, i);
      double kB = atomB.GetNddoParameterK(this->theory, i);
      double lB = atomB.GetNddoParameterL(this->theory, i);
      double mB = atomB.GetNddoParameterM(this->theory, i);
      temp1 += kA*exp(-lA*pow(Rab-mA,2.0));
      temp1 += kB*exp(-lB*pow(Rab-mB,2.0));
      temp2 += -2.0*lA*(Rab-mA)*kA*exp(-lA*pow(Rab-mA,2.0));
      temp2 += -2.0*lB*(Rab-mB)*kB*exp(-lB*pow(Rab-mB,2.0));
   }
   value -= dRabDa
           *atomA.GetCoreCharge()
           *atomB.GetCoreCharge()
           *temp1
           /(pow(Rab,2.0)/ang2AU);
   value += dRabDa
           *atomA.GetCoreCharge()
           *atomB.GetCoreCharge()
           *temp2/(Rab/ang2AU);
   return value;
}

void Am1::CalcSCFProperties(){
   MolDS_cndo::Cndo2::CalcSCFProperties();
}

void Am1::OutputSCFResults() const{
   MolDS_cndo::Cndo2::OutputSCFResults();
}

}



