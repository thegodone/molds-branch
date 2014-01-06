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
#include"Pm3Pddg.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_pm3{

/***
 *  Main References for PM3/PDDG are [RCJ_2002, BGRJ_2003, and BGJ_2003]
 */
Pm3Pddg::Pm3Pddg() : MolDS_pm3::Pm3(){
   this->theory = PM3PDDG;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //this->OutputLog("Pm3Pddg created\n");
}

Pm3Pddg::~Pm3Pddg(){
   //this->OutputLog("Pm3Pddg deleted\n");
}

void Pm3Pddg::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in pm3::Pm3Pddg::DoSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in pm3::Pm3Pddg::DoSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in pm3::Pm3Pddg::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in pm3::Pm3Pddg::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_pm3::Pm3Pddg::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_pm3::Pm3Pddg::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageCalcCISMatrix
      = "Error in pm3::Pm3Pddg::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in pm3::Pm3Pddg::DoCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadMultipoles
      = "Error in pm3::Pm3Pddg::GetSemiEmpiricalMultipoleInteraction: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteractionBadAtomTypes
      = "Error in pm3::Pm3Pddg::GetSemiEmpiricalMultipoleInteraction: Bad atom types are set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteraction1stDeriBadMultipoles
      = "Error in pm3::Pm3Pddg::GetSemiEmpiricalMultipoleInteraction1stDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetSemiEmpiricalMultipoleInteraction2ndDeriBadMultipoles
      = "Error in pm3::Pm3Pddg::GetSemiEmpiricalMultipoleInteraction2ndDerivative: Bad multipole combintaion is set\n";
   this->errorMessageGetNddoRepulsionIntegral 
      = "Error in pm3::Pm3Pddg::GetNddoRepulsionIntegral: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegralBadAtomTypes
      = "Error in pm3::Pm3Pddg::GetNddoRepulsionIntegral: Bad atom types are set.\n";
   this->errorMessageGetNddoRepulsionIntegral1stDerivative 
      = "Error in pm3::Pm3Pddg::GetNddoRepulsionIntegral1stDerivative: Bad orbital is set.\n";
   this->errorMessageGetNddoRepulsionIntegral2ndDerivative 
      = "Error in pm3::Pm3Pddg::GetNddoRepulsionIntegral2ndDerivative: Bad orbital is set.\n";
   this->errorMessageCalcTwoElecsTwoAtomCoresNullMatrix 
      = "Error in pm3::Pm3Pddg::CalcTwoElecsTwoAtomCores: The two elec two atom core matrix is NULL.\n"; 
   this->errorMessageCalcTwoElecsAtomEpcCoresNullMatrix 
      = "Error in pm3::Pm3Pddg::CalcTwoElecsAtomEpcCores: The two elec atom-epc core matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCoresSameAtoms
      = "Error in pm3::Pm3Pddg::CalcDiatomicTwoElecsTwoCores: Atom A and B is same atom (not EPC).\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCoresSameEpcs
      = "Error in pm3::Pm3Pddg::CalcDiatomicTwoElecsTwoCores: Atom A and B is same EPC.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores1stDerivativesSameAtoms
      = "Error in pm3::Pm3Pddg::CalcDiatomicTwoElecsTwoCores1stDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores2ndDerivativesSameAtoms
      = "Error in pm3::Pm3Pddg::CalcDiatomicTwoElecsTwoCores2ndDerivatives: Atom A and B is same.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCoresNullMatrix 
      = "Error in pm3::Pm3Pddg::CalcDiatomicTwoElecsTwoCores: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores1stDerivativesNullMatrix
      = "Error in pm3::Pm3Pddg::CalcDiatomicTwoElecsTwoCores1stDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageCalcDiatomicTwoElecsTwoCores2ndDerivativesNullMatrix
      = "Error in pm3::Pm3Pddg::CalcDiatomicTwoElecsTwoCores2ndDerivatives: The two elec two core diatomic matrix is NULL.\n"; 
   this->errorMessageGetElectronicEnergyEnergyNotCalculated
      = "Error in pm3::Pm3Pddg::GetElectronicEnergy: Set electronic state is not calculated by CIS.\n";
   this->errorMessageGetElectronicEnergyNULLCISEnergy 
      = "Error in pm3::Pm3Pddg::GetElectronicEnergy: excitedEnergies is NULL\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tPM3/PDDG-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: PM3/PDDG-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: PM3/PDDG-SCF  **********\n\n\n";
   this->messageStartCIS = "**********  START: PM3/PDDG-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: PM3/PDDG-CIS  **********\n\n\n";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for PM3/PDDG-CIS met convergence criterion(^^b\n\n\n";
}

void Pm3Pddg::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}

double Pm3Pddg::GetDiatomCoreRepulsionEnergy(const Atom& atomA, const Atom& atomB) const{
   // PM3 term
   double pm3Term = Pm3::GetDiatomCoreRepulsionEnergy(atomA, atomB);

   // pddg additional term, eq. (4) in [RCJ_2002]
   int na = atomA.GetNumberValenceElectrons();
   int nb = atomB.GetNumberValenceElectrons();
   double distance = this->molecule->GetDistanceAtoms(atomA, atomB);
   double temp = 0.0;
   for(int i=0; i<2; i++){
      double pa = atomA.GetPm3PddgParameterPa(i);
      double da = atomA.GetPm3PddgParameterDa(i);
      for(int j=0; j<2; j++){
         double pb = atomB.GetPm3PddgParameterPa(j);
         double db = atomB.GetPm3PddgParameterDa(j);
         temp += this->GetPddgAdditonalDiatomCoreRepulsionTerm(na, pa, da, nb, pb, db, distance);
      }
   }
   double additionalTerm = temp/(na+nb);
   
   return pm3Term + additionalTerm;
}

// First derivative of diatomic core repulsion energy.
// This derivative is related to the coordinate of atomA.
double Pm3Pddg::GetDiatomCoreRepulsion1stDerivative(const Atom& atomA, const Atom& atomB, 
                                                    CartesianType axisA) const{
   // PM3 term
   double pm3Term = Pm3::GetDiatomCoreRepulsion1stDerivative(atomA, atomB, axisA);

   // pddg additional term, first derivative of eq. (4) in [RCJ_2002]
   int na = atomA.GetNumberValenceElectrons();
   int nb = atomB.GetNumberValenceElectrons();
   double distance = this->molecule->GetDistanceAtoms(atomA, atomB);
   double temp = 0.0;
   for(int i=0; i<2; i++){
      double pa = atomA.GetPm3PddgParameterPa(i);
      double da = atomA.GetPm3PddgParameterDa(i);
      for(int j=0; j<2; j++){
         double pb = atomB.GetPm3PddgParameterPa(j);
         double db = atomB.GetPm3PddgParameterDa(j);
         temp += this->GetPddgAdditonalDiatomCoreRepulsionTerm1stDerivative(na, pa, da, nb, pb, db, distance);
      }
   }
   double dCartesian = (atomA.GetXyz()[axisA] - atomB.GetXyz()[axisA]);
   double additionalTerm = temp*dCartesian/(distance*(na+nb));

   return pm3Term + additionalTerm;
}

// Second derivative of diatomic core repulsion energy.
// Both derivative are related to the coordinate of atomA.
double Pm3Pddg::GetDiatomCoreRepulsion2ndDerivative(const Atom& atomA,
                                                    const Atom& atomB, 
                                                    CartesianType axisA1,
                                                    CartesianType axisA2) const{
   // PM3 term
   double pm3Term = Pm3::GetDiatomCoreRepulsion2ndDerivative(atomA, atomB, axisA1, axisA2);

   // pddg additional term, first derivative of eq. (4) in [RCJ_2002]
   double distance = this->molecule->GetDistanceAtoms(atomA, atomB);
   double dCartesian1 = (atomA.GetXyz()[axisA1] - atomB.GetXyz()[axisA1]);
   double dCartesian2 = (atomA.GetXyz()[axisA2] - atomB.GetXyz()[axisA2]);
   int na = atomA.GetNumberValenceElectrons();
   int nb = atomB.GetNumberValenceElectrons();
   double pddgExponent = -10.0;
   double temp1stDeriv = 0.0;
   double temp2ndDeriv = 0.0;
   for(int i=0; i<2; i++){
      double pa = atomA.GetPm3PddgParameterPa(i);
      double da = atomA.GetPm3PddgParameterDa(i);
      for(int j=0; j<2; j++){
         double pb = atomB.GetPm3PddgParameterPa(j);
         double db = atomB.GetPm3PddgParameterDa(j);
         temp1stDeriv +=  this->GetPddgAdditonalDiatomCoreRepulsionTerm1stDerivative(na, pa, da, nb, pb, db, distance);
         temp2ndDeriv += this->GetPddgAdditonalDiatomCoreRepulsionTerm2ndDerivative(na, pa, da, nb, pb, db, distance);
      }
   }
   double pre1stDeriv = 0.0;
   double pre2ndDeriv = 0.0;
   if(axisA1 != axisA2){
      pre1stDeriv = -dCartesian1*dCartesian2/pow(distance,3.0);
      pre2ndDeriv = dCartesian1*dCartesian2/pow(distance,2.0);
   }
   else{
      pre1stDeriv = 1.0/distance - dCartesian1*dCartesian1/pow(distance,3.0);
      pre2ndDeriv = pow(dCartesian1/distance,2.0);
   }
   pre1stDeriv  /= static_cast<double>(na+nb);
   pre2ndDeriv /= static_cast<double>(na+nb);
   double additionalTerm = pre1stDeriv*temp1stDeriv + pre2ndDeriv*temp2ndDeriv;

   return pm3Term + additionalTerm;
}

// see eq. (4) in [RCJ_2002]
double Pm3Pddg::GetPddgAdditonalDiatomCoreRepulsionTerm(int na, double pa, double da,
                                                        int nb, double pb, double db,
                                                        double distance) const{
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double pddgExponent = -10.0/(ang2AU*ang2AU);
   return (static_cast<double>(na)*pa +static_cast<double>(nb)*pb)*exp(pddgExponent*pow((distance-da-db),2.0));
}

// see eq. (4) in [RCJ_2002]
double Pm3Pddg::GetPddgAdditonalDiatomCoreRepulsionTerm1stDerivative(int na, double pa, double da,
                                                                       int nb, double pb, double db,
                                                                       double distance) const{
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double pddgExponent = -10.0/(ang2AU*ang2AU);
   return (static_cast<double>(na)*pa +static_cast<double>(nb)*pb)*exp(pddgExponent*pow((distance-da-db),2.0))
         *(2.0*pddgExponent*(distance-da-db));
}

// see eq. (4) in [RCJ_2002]
double Pm3Pddg::GetPddgAdditonalDiatomCoreRepulsionTerm2ndDerivative(int na, double pa, double da,
                                                                        int nb, double pb, double db,
                                                                        double distance) const{
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double pddgExponent = -10.0/(ang2AU*ang2AU);
   return (static_cast<double>(na)*pa +static_cast<double>(nb)*pb)
         *(2.0*pddgExponent + pow(2.0*pddgExponent*(distance-da-db),2.0))
         *exp(pddgExponent*pow((distance-da-db),2.0));
}

}



