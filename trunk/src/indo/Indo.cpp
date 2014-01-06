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
#include<stdexcept>
#include<boost/format.hpp>
#include"../config.h"
#include"../base/Enums.h"
#include"../base/Uncopyable.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
#include"../base/EularAngle.h"
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
#include"Indo.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_indo{

Indo::Indo() : MolDS_cndo::Cndo2(){
   this->theory = INDO;
   this->SetMessages();
   this->SetEnableAtomTypes();
   //this->OutputLog("Indo created\n");
}

Indo::~Indo(){
   //this->OutputLog("Indo deleted\n");
}

void Indo::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in indo::Indo::DoSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in indo::Indo::DoSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in indo::Indo::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in indo::Indo::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in base_indo::Indo::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in base_indo::Indo::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageMolecularIntegralElement
      = "Error in indo::Indo::GetMolecularIntegralElement: Non available orbital is contained.\n";
   this->errorMessageGetDiatomCoreRepulsion2ndDerivativeNotImplemented
      = "Error in indo::Indo::GetDiatomCoreRepulsion2ndDerivative: Second derivative is not implemented for INDO.\n";
   this->errorMessageCISNotImplemented 
      = "Error in indo::Indo::DoCIS: CIS is not implemented for INDO.\n";
   this->errorMessageCalcForceNotImplemented
      = "Error in indo::Indo::CalcForce: Force is not available in INDO.\n";
   this->errorMessageGetElectronicEnergyEnergyNotCalculated
      = "Error in indo::Indo::GetElectronicEnergy: Set electronic state is not calculated by CIS.\n";
   this->errorMessageGetElectronicEnergyNULLCISEnergy 
      = "Error in indo::Indo::GetElectronicEnergy: excitedEnergies is NULL\n";
   this->errorMessageCalcFrequenciesNormalModesBadTheory
      = "Error in indo::Indo::CalcFrequenciesNormalModesBadTheory: INDO is not supported for frequency (normal mode) analysis.\n";
   this->errorMessageNonExcitedStates 
      = "Error in indo::Indo::Excited states can not be calculated with INDO.\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tINDO-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: INDO-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: INDO-SCF  **********\n\n\n";
}

void Indo::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(Li);
   //this->enableAtomTypes.push_back(Be);
   //this->enableAtomTypes.push_back(B);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   //this->enableAtomTypes.push_back(F);
}

double Indo::GetFockDiagElement(const Atom& atomA, 
                                int indexAtomA, 
                                int mu, 
                                const Molecule& molecule, 
                                double const* const* gammaAB,
                                double const* const* orbitalElectronPopulation, 
                                double const* atomicElectronPopulation,
                                double const* const* const* const* const* const* twoElecsTwoAtomCores,
                                bool isGuess) const{
   double value;
   int firstAOIndexA = atomA.GetFirstAOIndex();
   value = atomA.GetCoreIntegral(atomA.GetValence(mu-firstAOIndexA), 
                                 gammaAB[indexAtomA][indexAtomA], 
                                 isGuess, this->theory);

   if(!isGuess){
      double temp = 0.0;
      double coulomb = 0.0;
      double exchange = 0.0;
      int lammda = 0;
      OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
      for(int v=0; v<atomA.GetValenceSize(); v++){
         OrbitalType orbitalLam = atomA.GetValence(v);
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalLam, gammaAB[indexAtomA][indexAtomA], atomA);
         exchange = this->GetExchangeInt(orbitalMu, orbitalLam, gammaAB[indexAtomA][indexAtomA], atomA);
         lammda = firstAOIndexA + v;
         temp += orbitalElectronPopulation[lammda][lammda]*(coulomb - 0.5*exchange);
      }
      value += temp;
   
      temp = 0.0;
      for(int B=0; B<molecule.GetAtomVect().size(); B++){
         if(B != indexAtomA){
            const Atom& atomB = *molecule.GetAtomVect()[B];
            temp += ( atomicElectronPopulation[B] - atomB.GetCoreCharge()  )
                     *gammaAB[indexAtomA][B];
         }
      }
      value += temp;
   }

   return value;
}

double Indo::GetFockOffDiagElement(const Atom& atomA, 
                                   const Atom& atomB, 
                                   int indexAtomA, 
                                   int indexAtomB, 
                                   int mu, 
                                   int nu, 
                                   const Molecule& molecule, 
                                   double const* const* gammaAB, 
                                   double const* const* overlapAOs,
                                   double const* const* orbitalElectronPopulation, 
                                   double const* const* const* const* const* const* twoElecsTwoAtomCores,
                                   bool isGuess) const{
   double value;
   double K = this->GetBondingAdjustParameterK(atomA.GetValenceShellType(), atomB.GetValenceShellType());
   double bondParameter = 0.5*K*(atomA.GetBondingParameter() + atomB.GetBondingParameter()); 

   if(isGuess){
      value = bondParameter*overlapAOs[mu][nu];
   }
   else{
      double coulomb = 0.0;
      double exchange = 0.0;
      if(indexAtomA == indexAtomB){
         OrbitalType orbitalMu = atomA.GetValence(mu-atomA.GetFirstAOIndex());
         OrbitalType orbitalNu = atomA.GetValence(nu-atomA.GetFirstAOIndex());
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalNu, gammaAB[indexAtomA][indexAtomA], atomA); 
         exchange = this->GetExchangeInt(orbitalMu, orbitalNu, gammaAB[indexAtomA][indexAtomA], atomA); 
         value = (1.5*exchange - 0.5*coulomb)*orbitalElectronPopulation[mu][nu];
      }
      else{
         value = bondParameter*overlapAOs[mu][nu];
         value -= 0.5*orbitalElectronPopulation[mu][nu]*gammaAB[indexAtomA][indexAtomB];
      }
   }

   return value;
}

// The order of mol, moJ, moK, moL is consistent with Eq. (9) in [RZ_1973]
double Indo::GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                         const Molecule& molecule, 
                                         double const* const* fockMatrix, 
                                         double const* const* gammaAB) const{
   double value = 0.0;
   Atom* atomA;
   Atom* atomB;;
   int firstAOIndexA;
   int numberAOsA;
   double exchange;
   double coulomb;
   OrbitalType orbitalMu;
   OrbitalType orbitalNu;

   // CNDO terms
   value = Cndo2::GetMolecularIntegralElement(moI, moJ, moK, moL, molecule, fockMatrix, gammaAB);

   // Aditional terms for INDO, see Eq. (10) in [RZ_1973]
   for(int A=0; A<molecule.GetAtomVect().size(); A++){
      const Atom& atomA = *molecule.GetAtomVect()[A];
      firstAOIndexA = atomA.GetFirstAOIndex();
      numberAOsA = atomA.GetValenceSize();

      for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
         orbitalMu = atomA.GetValence(mu-firstAOIndexA);
         for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
            orbitalNu = atomA.GetValence(nu-firstAOIndexA);

            if(mu!=nu){
               exchange = this->GetExchangeInt(orbitalMu, orbitalNu, gammaAB[A][A], atomA);
               value += exchange
                       *fockMatrix[moI][mu]
                       *fockMatrix[moJ][nu]
                       *fockMatrix[moK][nu]
                       *fockMatrix[moL][mu];
               value += exchange
                       *fockMatrix[moI][mu]
                       *fockMatrix[moJ][nu]
                       *fockMatrix[moK][mu]
                       *fockMatrix[moL][nu];
            }

            coulomb = this->GetCoulombInt(orbitalMu, orbitalNu, gammaAB[A][A], atomA);
            value += (coulomb-gammaAB[A][A])
                     *fockMatrix[moI][mu]
                     *fockMatrix[moJ][mu]
                     *fockMatrix[moK][nu]
                     *fockMatrix[moL][nu];
         }
      }
   }

   return value;
}

// (3.87) - (3.91) in J. A. Pople book.
// Indo Coulomb Interaction
double Indo::GetCoulombInt(OrbitalType orbital1, OrbitalType orbital2, double gamma, const Atom& atom) const{

   double value=0.0;
   if( orbital1 == s && orbital2 == s){ 
      value = gamma;
   }   
   else if( orbital1 == s && ( orbital2 == px || orbital2 == py || orbital2 == pz )){ 
      value = gamma;
   }   
   else if( (orbital1 == px || orbital1 == py || orbital1 == pz ) && orbital2 == s){ 
      value = gamma;
   }   
   else if( (orbital1 == orbital2) && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = gamma + 4.0*atom.GetIndoF2()/25.0;
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = gamma - 2.0*atom.GetIndoF2()/25.0;
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

// (3.87) - (3.91) in J. A. Pople book.
// Indo Exchange Interaction
double Indo::GetExchangeInt(OrbitalType orbital1, OrbitalType orbital2, double gamma, const Atom& atom) const{

   double value=0.0;

   if( orbital1 == orbital2){
      value = this->GetCoulombInt(orbital1, orbital2, gamma, atom);
   }   
   else if( (orbital1 == s) && (orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom.GetIndoG1()/3.0;
   }   
   else if( (orbital1 == px || orbital1 == py || orbital1 == pz) && orbital2 == s  ){  
      value = atom.GetIndoG1()/3.0;
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = 3.0*atom.GetIndoF2()/25.0;
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


}



