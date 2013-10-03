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
#include<boost/format.hpp>
#include"../base/Uncopyable.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
#include"../base/Enums.h"
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
#include"Am1D.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_am1{

/***
 *  Main References for AM1-D are [MH_2007, MMHBV_2007]
 */
Am1D::Am1D() : MolDS_am1::Am1(){
   this->theory = AM1D;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //this->OutputLog("Am1D created\n");
}

Am1D::~Am1D(){
   //this->OutputLog("Am1d deleted\n");
}

void Am1D::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in am1::Am1D::DoSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in am1::Am1D::DoSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in am1::Am1D::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in am1::Am1D::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in am1::Am1D::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in am1::Am1D::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageCalcCISMatrix
      = "Error in am1::Am1D::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in am1::Am1D::DoCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles
      = "Error in am1::Am1D::GetSemiEmpiricalMultipoleInteraction: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteraction1stDeriBadMultipoles
      = "Error in am1::Am1D::GetSemiEmpiricalMultipoleInteraction1stDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteraction2ndDeriBadMultipoles
      = "Error in am1::Am1D::GetSemiEmpiricalMultipoleInteraction2ndDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetNddoRepulsionIntegral 
      = "Error in am1::Am1D::GetNddoRepulsionIntegral: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegral1stDerivative 
      = "Error in am1::Am1D::GetNddoRepulsionIntegral1stDerivative: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegral2ndDerivative 
      = "Error in am1::Am1D::GetNddoRepulsionIntegral2ndDerivative: Bad orbital is set.\n";
   this->errorMessageCalcTwoElecTwoCoreNullMatrix 
      = "Error in am1::Am1D::CalcTwoElecTwoCore: The two elec two core matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCoreSameAtoms
      = "Error in am1::Am1D::CalcDiatomicTwoElecTwoCore: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCore1stDerivativesSameAtoms
      = "Error in am1::Am1D::CalcDiatomicTwoElecTwoCore1stDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCore2ndDerivativesSameAtoms
      = "Error in am1::Am1D::CalcDiatomicTwoElecTwoCore2ndDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCoreNullMatrix 
      = "Error in am1::Am1D::CalcDiatomicTwoElecTwoCore: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCore1stDerivativesNullMatrix
      = "Error in am1::Am1D::CalcDiatomicTwoElecTwoCore1stDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecTwoCore2ndDerivativesNullMatrix
      = "Error in am1::Am1D::CalcDiatomicTwoElecTwoCore2ndDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageGetElectronicEnergyEnergyNotCalculated
      = "Error in am1::Am1D::GetElectronicEnergy: Set electronic state is not calculated by CIS.\n";
   this->errorMessageGetElectronicEnergyNULLCISEnergy 
      = "Error in am1::Am1D::GetElectronicEnergy: excitedEnergies is NULL\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tAM1-D-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: AM1-D-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: AM1-D-SCF  **********\n\n\n";
   this->messageStartCIS = "**********  START: AM1-D-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: AM1-D-CIS  **********\n\n\n";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for AM1-D-CIS met convergence criterion(^^b\n\n\n";
}

void Am1D::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}


}



