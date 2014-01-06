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
#include"Pm3D.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_pm3{

/***
 *  Main References for PM3-D are [MH_2007, MMHBV_2007]
 */
Pm3D::Pm3D() : MolDS_pm3::Pm3(){
   this->theory = PM3D;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //this->OutputLog("Pm3D created\n");
}

Pm3D::~Pm3D(){
   //this->OutputLog("Pm3d deleted\n");
}

void Pm3D::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in pm3::Pm3D::DoSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in pm3::Pm3D::DoSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in pm3::Pm3D::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in pm3::Pm3D::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in pm3::Pm3D::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in pm3::Pm3D::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageCalcCISMatrix
      = "Error in pm3::Pm3D::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in pm3::Pm3D::DoCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles
      = "Error in pm3::Pm3D::GetSemiEmpiricalMultipoleInteraction: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadAtomTypes
      = "Error in pm3::Pm3D::GetSemiEmpiricalMultipoleInteraction: Bad atom types are set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteraction1stDeriBadMultipoles
      = "Error in pm3::Pm3D::GetSemiEmpiricalMultipoleInteraction1stDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteraction2ndDeriBadMultipoles
      = "Error in pm3::Pm3D::GetSemiEmpiricalMultipoleInteraction2ndDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetNddoRepulsionIntegral 
      = "Error in pm3::Pm3D::GetNddoRepulsionIntegral: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegralBadAtomTypes
      = "Error in pm3::Pm3D::GetNddoRepulsionIntegral: Bad atom types are set.\n";
   this->errorMessageGetNddoRepulsionIntegral1stDerivative 
      = "Error in pm3::Pm3D::GetNddoRepulsionIntegral1stDerivative: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegral2ndDerivative 
      = "Error in pm3::Pm3D::GetNddoRepulsionIntegral2ndDerivative: Bad orbital is set.\n";
   this->errorMessageCalcTwoElecsTwoAtomCoresNullMatrix 
      = "Error in pm3::Pm3D::CalcTwoElecsTwoAtomCores: The two elec two atom core matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecsAtomEpcCoresNullMatrix 
      = "Error in pm3::Pm3D::CalcTwoElecsAtomEpcCores: The two elec atom-epc core matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCoresSameAtoms
      = "Error in pm3::Pm3D::CalcDiatomicTwoElecsTwoCores: Atom A and B is same atom (not EPC).\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCoresSameEpcs
      = "Error in pm3::Pm3D::CalcDiatomicTwoElecsTwoCores: Atom A and B is same EPC.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores1stDerivativesSameAtoms
      = "Error in pm3::Pm3D::CalcDiatomicTwoElecsTwoCores1stDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores2ndDerivativesSameAtoms
      = "Error in pm3::Pm3D::CalcDiatomicTwoElecsTwoCores2ndDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCoresNullMatrix 
      = "Error in pm3::Pm3D::CalcDiatomicTwoElecsTwoCores: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores1stDerivativesNullMatrix
      = "Error in pm3::Pm3D::CalcDiatomicTwoElecsTwoCores1stDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores2ndDerivativesNullMatrix
      = "Error in pm3::Pm3D::CalcDiatomicTwoElecsTwoCores2ndDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageGetElectronicEnergyEnergyNotCalculated
      = "Error in pm3::Pm3D::GetElectronicEnergy: Set electronic state is not calculated by CIS.\n";
   this->errorMessageGetElectronicEnergyNULLCISEnergy 
      = "Error in pm3::Pm3D::GetElectronicEnergy: excitedEnergies is NULL\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tPM3-D-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: PM3-D-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: PM3-D-SCF  **********\n\n\n";
   this->messageStartCIS = "**********  START: PM3-D-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: PM3-D-CIS  **********\n\n\n";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for PM3-D-CIS met convergence criterion(^^b\n\n\n";
}

void Pm3D::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}


}



