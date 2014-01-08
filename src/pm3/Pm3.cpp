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
#include"../am1/Am1.h"
#include"Pm3.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_pm3{

/***
 *  Main References for PM3 are [S_1989, S_1989-2, S_1991, S_2004, S_2007]
 */
Pm3::Pm3() : MolDS_am1::Am1(){
   this->theory = PM3;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //this->OutputLog("Pm3 created\n");
}

Pm3::~Pm3(){
   //this->OutputLog("Pm3 deleted\n");
}

void Pm3::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in pm3::Pm3::DoSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in pm3::Pm3::DoSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in pm3::Pm3::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in pm3::Pm3::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_pm3::Pm3::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_pm3::Pm3::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageCalcCISMatrix
      = "Error in pm3::Pm3::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in pm3::Pm3::DoCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles
      = "Error in pm3::Pm3::GetSemiEmpiricalMultipoleInteraction: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadAtomTypes
      = "Error in pm3::Pm3::GetSemiEmpiricalMultipoleInteraction: Bad atom types are set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteraction1stDeriBadMultipoles
      = "Error in pm3::Pm3::GetSemiEmpiricalMultipoleInteraction1stDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteraction2ndDeriBadMultipoles
      = "Error in pm3::Pm3::GetSemiEmpiricalMultipoleInteraction2ndDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetNddoRepulsionIntegral 
      = "Error in pm3::Pm3::GetNddoRepulsionIntegral: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegralBadAtomTypes
      = "Error in pm3::Pm3::GetNddoRepulsionIntegral: Bad atom types are set.\n";
   this->errorMessageGetNddoRepulsionIntegral1stDerivative 
      = "Error in pm3::Pm3::GetNddoRepulsionIntegral1stDerivative: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegral2ndDerivative 
      = "Error in pm3::Pm3::GetNddoRepulsionIntegral2ndDerivative: Bad orbital is set.\n";
   this->errorMessageCalcTwoElecsTwoAtomCoresNullMatrix 
      = "Error in pm3::Pm3::CalcTwoElecsTwoAtomCores: The two elec two atom core matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecsAtomEpcCoresNullMatrix 
      = "Error in pm3::Pm3::CalcTwoElecsAtomEpcCores: The two elec atom-epc core matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCoresSameAtoms
      = "Error in pm3::Pm3::CalcDiatomicTwoElecsTwoCores: Atom A and B is same atom (not EPC).\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCoresSameEpcs
      = "Error in pm3::Pm3::CalcDiatomicTwoElecsTwoCores: Atom A and B is same EPC.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores1stDerivativesSameAtoms
      = "Error in pm3::Pm3::CalcDiatomicTwoElecsTwoCores1stDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores2ndDerivativesSameAtoms
      = "Error in pm3::Pm3::CalcDiatomicTwoElecsTwoCores2ndDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCoresNullMatrix 
      = "Error in pm3::Pm3::CalcDiatomicTwoElecsTwoCores: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores1stDerivativesNullMatrix
      = "Error in pm3::Pm3::CalcDiatomicTwoElecsTwoCores1stDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores2ndDerivativesNullMatrix
      = "Error in pm3::Pm3::CalcDiatomicTwoElecsTwoCores2ndDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageGetElectronicEnergyEnergyNotCalculated
      = "Error in pm3::Pm3::GetElectronicEnergy: Set electronic state is not calculated by CIS.\n";
   this->errorMessageGetElectronicEnergyNULLCISEnergy 
      = "Error in pm3::Pm3::GetElectronicEnergy: excitedEnergies is NULL\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tPM3-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: PM3-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: PM3-SCF  **********\n\n\n";
   this->messageStartCIS = "**********  START: PM3-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: PM3-CIS  **********\n\n\n";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for PM3-CIS met convergence criterion(^^b\n\n\n";
}

void Pm3::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
   this->enableAtomTypes.push_back(Zn);
}

}



