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
#include<algorithm>
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
#include"../base/loggers/DensityLogger.h"
#include"../base/loggers/HoleDensityLogger.h"
#include"../base/loggers/ParticleDensityLogger.h"
#include"../cndo/Cndo2.h"
#include"ZindoS.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_zindo{

/***
 *  Main Refference for Zindo is [RZ_1973]
 */

ZindoS::ZindoS() : MolDS_cndo::Cndo2(){
   this->theory = ZINDOS;
   this->SetMessages();
   this->SetEnableAtomTypes();
   this->matrixForceElecStatesNum = 0;
   this->nishimotoMatagaParamA = 1.2;
   this->nishimotoMatagaParamB = 2.4;
   this->overlapCorrectionSigma = 1.267;
   this->overlapCorrectionPi = 0.585;
   //this->OutputLog("ZindoS created\n");
}

ZindoS::~ZindoS(){
   MallocerFreer::GetInstance()->Free<double>(&this->matrixCIS, 
                                              this->matrixCISdimension,
                                              this->matrixCISdimension);
   MallocerFreer::GetInstance()->Free<double>(&this->excitedEnergies, 
                                              this->matrixCISdimension);
   MallocerFreer::GetInstance()->Free<double>(&this->freeExcitonEnergiesCIS, 
                                              this->matrixCISdimension);
   MallocerFreer::GetInstance()->Free<double>(&this->matrixForce, 
                                              this->matrixForceElecStatesNum,
                                              this->molecule->GetNumberAtoms(),
                                              CartesianType_end);
   //this->OutputLog("ZindoS deleted\n");
}

void ZindoS::SetMessages(){
   this->errorMessageSCFNotConverged 
      = "Error in zindo::ZindoS::DoSCF: SCF did not met convergence criterion. maxIterationsSCF=";
   this->errorMessageMoleculeNotSet 
      = "Error in zindo::ZindoS::DoSCF: A molecule is not set.\n";
   this->errorMessageOddTotalValenceElectrions 
      = "Error in zindo::ZindoS::SetMolecule: Total number of valence electrons is odd. totalNumberValenceElectrons=";
   this->errorMessageNotEnebleAtomType  
      = "Error in zindo::ZindoS::CheckEnableAtomType: Non available atom is contained.\n";
   this->errorMessageCoulombInt = "Error in zindo::ZindoS::GetCoulombInt: Invalid orbitalType.\n";
   this->errorMessageExchangeInt = "Error in zindo::ZindoS::GetExchangeInt: Invalid orbitalType.\n";
   this->errorMessageNishimotoMataga = "Error in zindo::ZindoS::GetNishimotoMatagaTwoEleInt: Invalid orbitalType.\n";
   this->errorMessageMolecularIntegralElement
      = "Error in zindo::ZindoS::GetMolecularIntegralElement: Non available orbital is contained.\n";
   this->errorMessageCalcCISMatrix
      = "Error in zindo::ZindoS::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in zindo::ZindoS::DoCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageDavidsonMaxIter = "Davidson roop reaches max_iter=";
   this->errorMessageDavidsonMaxDim = "Dimension of the expansion vectors reaches max_dim=";
   this->errorMessageCalcForceNotGroundState 
      = "Error in zindo::ZindoS::CalcForce: Only ground state is enable in ZindoS.";
   this->errorMessageElecState = "Electronic State = ";
   this->errorMessageGetElectronicEnergyEnergyNotCalculated
      = "Error in zindo::ZindoS::GetElectronicEnergy: Set electronic state is not calculated by CIS.\n";
   this->errorMessageGetElectronicEnergyNULLCISEnergy 
      = "Error in zindo::ZindoS::GetElectronicEnergy: excitedEnergies is NULL\n";
   this->errorMessageGetElectronicTransitionDipoleMomentBadState
      = "Error in zindo::ZindoS::GetElectronicTransitionDipoleMoment: Bad eigen state is set to calculate the transition dipole moment. Note taht state=0 means the ground state and other state = i means the i-th excited state in below.\n";
   this->messageSCFMetConvergence = "\n\n\n\t\tZINDO/S-SCF met convergence criterion(^^b\n\n\n";
   this->messageStartSCF = "**********  START: ZINDO/S-SCF  **********\n";
   this->messageDoneSCF = "**********  DONE: ZINDO/S-SCF  **********\n\n\n";
   this->messageStartCIS = "**********  START: ZINDO/S-CIS  **********\n";
   this->messageDoneCIS = "**********  DONE: ZINDO/S-CIS  **********\n\n\n";
   this->messageOmpElapsedTimeCalcCISMarix = "\tElapsed time(omp) for the calc. of the CIS matrix = ";
   this->messageOmpElapsedTimeCIS = "\tElapsed time(omp) for the CIS = ";
   this->messageStartCalcCISMatrix = "----------- START: Calculation of the CIS matrix -----------\n";
   this->messageDoneCalcCISMatrix  = "----------- DONE: Calculation of the CIS matrix -----------\n\n";
   this->messageStartDirectCIS = "\t======  START: Direct-CIS  =====\n\n";
   this->messageDoneDirectCIS =  "\t======  DONE: Direct-CIS  =====\n\n\n";
   this->messageStartDavidsonCIS = "\t======  START: Davidson-CIS  =====\n";
   this->messageDoneDavidsonCIS =  "\t======  DONE: Davidson-CIS  =====\n\n\n";
   this->messageNumIterCIS = "\tDavidson iter=";
   this->messageResidualNorm = "-th excited: norm of the residual = ";
   this->messageDavidsonConverge = "\n\n\t\tDavidson for ZINDO/S-CIS met convergence criterion(^^b\n\n\n";
   this->messageDavidsonReachCISMatrix = "\n\t\tDimension of the expansion vectors reaches to the dimension of the CIS-matrix.\n";
   this->messageDavidsonGoToDirect = "\t\tHence, we go to the Direct-CIS.\n\n";
   this->messageExcitedStatesEnergies = "\tExcitation energies:";
   this->messageExcitedStatesEnergiesTitle = "\t\t\t\t|   i-th   |   e[a.u.]   |   e[eV]   | dominant eigenvector coefficients (occ. -> vir.) |\n";
   this->messageExcitonEnergiesCIS = "\tFree exciton (Ef) and exciton binding (Eb) energies:\n";
   this->messageExcitonEnergiesShortCIS = "\tEf and Eb:";
   this->messageExcitonEnergiesCISTitle = "\t\t\t|   i-th   |   Ef[a.u.]   |   Ef[eV]   |   Eb[a.u.]   |   Eb[eV]   |\n";
   this->messageTotalDipoleMomentsTitle = "\t\t\t\t| i-th eigenstate |  x[a.u.]  |  y[a.u.]  |  z[a.u.]  |  magnitude[a.u.]  |\t\t|  x[debye]  |  y[debye]  |  z[debye]  |  magnitude[debye]  |\n";
   this->messageTotalDipoleMoment = "Total dipole moment:";
   this->messageElectronicDipoleMomentsTitle = "\t\t\t\t\t| i-th eigenstate |  x[a.u.]  |  y[a.u.]  |  z[a.u.]  |  magnitude[a.u.]  |\t\t|  x[debye]  |  y[debye]  |  z[debye]  |  magnitude[debye]  |\n";
   this->messageElectronicDipoleMoment = "Electronic dipole moment:";
   this->messageTransitionDipoleMomentsTitle = "\t\t\t\t\t| from and to eigenstates |  x[a.u.]  |  y[a.u.]  |  z[a.u.]  |  magnitude[a.u.]  |\t\t|  x[debye]  |  y[debye]  |  z[debye]  |  magnitude[debye]  |\n";
   this->messageTransitionDipoleMoment = "Transition dipole moment:";

}

void ZindoS::SetEnableAtomTypes(){
   this->enableAtomTypes.clear();
   this->enableAtomTypes.push_back(H);
   this->enableAtomTypes.push_back(C);
   this->enableAtomTypes.push_back(N);
   this->enableAtomTypes.push_back(O);
   this->enableAtomTypes.push_back(S);
}

double ZindoS::GetFockDiagElement(const Atom& atomA, 
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
   value = atomA.GetCoreIntegral(atomA.GetValence(mu-firstAOIndexA), 
                                  isGuess, this->theory);
   if(!isGuess){
      double temp = 0.0;
      double coulomb = 0.0;
      double exchange = 0.0;
      int lammda = 0;
      int totalNumberAOs = molecule.GetTotalNumberAOs();
      double orbitalElectronPopulationDiagPart[totalNumberAOs];

      for(int i=0; i<totalNumberAOs; i++){
         orbitalElectronPopulationDiagPart[i] = orbitalElectronPopulation[i][i];
      }

      OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
      OrbitalType orbitalLam;
      int atomANumberValence = atomA.GetValenceSize();
      for(int v=0; v<atomANumberValence; v++){
         orbitalLam = atomA.GetValence(v);
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalLam, atomA);
         exchange = this->GetExchangeInt(orbitalMu, orbitalLam, atomA);
         lammda = v + firstAOIndexA;
         temp += orbitalElectronPopulationDiagPart[lammda]*(coulomb - 0.5*exchange);
      }
      value += temp;
   
      temp = 0.0;
      int totalNumberAtoms = molecule.GetNumberAtoms();
      for(int B=0; B<totalNumberAtoms; B++){
         if(B != atomAIndex){
            const Atom& atomB = *molecule.GetAtom(B);
            OrbitalType orbitalSigma;
            int sigma;
            int atomBNumberValence = atomB.GetValenceSize();
            for(int i=0; i<atomBNumberValence; i++){
               sigma = i + atomB.GetFirstAOIndex();
               orbitalSigma = atomB.GetValence(i);
               temp += orbitalElectronPopulationDiagPart[sigma]
                      *this->GetNishimotoMatagaTwoEleInt(atomA, 
                                                         orbitalMu, 
                                                         atomB, 
                                                         orbitalSigma);
            }
            temp -= atomB.GetCoreCharge() 
                   *this->GetNishimotoMatagaTwoEleInt(atomA, s, atomB, s);
         }
      }
      value += temp;
   }

   return value;
}

double ZindoS::GetFockOffDiagElement(const Atom& atomA, 
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
   OrbitalType orbitalMu = atomA.GetValence(mu-atomA.GetFirstAOIndex());
   OrbitalType orbitalNu = atomB.GetValence(nu-atomB.GetFirstAOIndex());
   double bondParameter = 0.5*(atomA.GetBondingParameter(this->theory, orbitalMu) 
                              +atomB.GetBondingParameter(this->theory, orbitalNu)); 

   if(isGuess){
      value = bondParameter*overlap[mu][nu];
   }
   else{
      double coulomb = 0.0;
      double exchange = 0.0;
      if(atomAIndex == atomBIndex){
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalNu, atomA); 
         exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA); 
         value = (1.5*exchange - 0.5*coulomb)*orbitalElectronPopulation[mu][nu];
      }
      else{
         value = bondParameter*overlap[mu][nu];
         value -= 0.5*orbitalElectronPopulation[mu][nu]
                  *this->GetNishimotoMatagaTwoEleInt(atomA, orbitalMu, atomB, orbitalNu);
      }
   }
   return value;
}

void ZindoS::CalcGammaAB(double** gammaAB, const Molecule& molecule) const{
   // Do nothing;
}

// Apendix in [BZ_1972]
// ZINDO Coulomb Interaction
double ZindoS::GetCoulombInt(OrbitalType orbital1, OrbitalType orbital2, const Atom& atom) const{

   double value=0.0;
   
   if( orbital1 == s && orbital2 == s){ 
      value = atom.GetZindoF0ssLower();
   }   
   else if( orbital1 == s && ( orbital2 == px || orbital2 == py || orbital2 == pz )){ 
      value = atom.GetZindoF0ssLower();
   }   
   else if( orbital2 == s && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = atom.GetZindoF0ssLower();
   }   
   else if( (orbital1 == orbital2) && ( orbital1 == px || orbital1 == py || orbital1 == pz )){ 
      value = atom.GetZindoF0ssLower()
             +atom.GetZindoF2ppLower()*4.0;
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom.GetZindoF0ssLower()
             -atom.GetZindoF2ppLower()*2.0;
   }   
   // ToDo: There are bugs for d-orbitals.
   /*
   else if( orbital1 == s && ( orbital2 == dxy || 
                               orbital2 == dyz || 
                               orbital2 == dzz || 
                               orbital2 == dzx || 
                               orbital2 == dxxyy )){ 
      value = atom.GetZindoF0sdLower();
   }   
   else if( orbital2 == s && ( orbital1 == dxy || 
                               orbital1 == dyz || 
                               orbital1 == dzz || 
                               orbital1 == dzx || 
                               orbital1 == dxxyy )){ 
      value = atom.GetZindoF0sdLower();
   }
   else if( orbital1 == dzz && (orbital2 == px || orbital2==py) ){
      value = atom.GetZindoF0sdLower()
             -atom.GetZindoF2pdLower()*2.0;
   }
   else if( orbital2 == dzz && (orbital1 == px || orbital1==py) ){
      value = atom.GetZindoF0sdLower()
             -atom.GetZindoF2pdLower()*2.0;
   }
   else if( (orbital1 == dzz && orbital2 == pz) ||
            (orbital2 == dzz && orbital1 == pz) ){
      value = atom.GetZindoF0sdLower()
             +atom.GetZindoF2pdLower()*4.0;
   }
   else if( (orbital1 == orbital2) && ( orbital1 == dxy || 
                                        orbital1 == dyz || 
                                        orbital1 == dzz || 
                                        orbital1 == dzx || 
                                        orbital1 == dxxyy )){ 
      value = atom.GetZindoF0ddLower()
             +atom.GetZindoF2ddLower()*4.0
             +atom.GetZindoF4ddLower()*36.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == px) ||
            (orbital2 == dxxyy && orbital1 == px) || 
            (orbital1 == dxxyy && orbital2 == py) ||
            (orbital2 == dxxyy && orbital1 == py) || 
            (orbital1 == dxy && orbital2 == px) ||
            (orbital2 == dxy && orbital1 == px) || 
            (orbital1 == dxy && orbital2 == py) ||
            (orbital2 == dxy && orbital1 == py) ||
            (orbital1 == dzx && orbital2 == px) ||
            (orbital2 == dzx && orbital1 == px) || 
            (orbital1 == dzx && orbital2 == pz) ||
            (orbital2 == dzx && orbital1 == pz) || 
            (orbital1 == dyz && orbital2 == py) ||
            (orbital2 == dyz && orbital1 == py) || 
            (orbital1 == dyz && orbital2 == pz) ||
            (orbital2 == dyz && orbital1 == pz) ){
      value = atom.GetZindoF0sdLower()
             +atom.GetZindoF2pdLower()*2.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == pz) ||
            (orbital2 == dxxyy && orbital1 == pz) || 
            (orbital1 == dxy && orbital2 == pz) ||
            (orbital2 == dxy && orbital1 == pz) ||
            (orbital1 == dzx && orbital2 == py) ||
            (orbital2 == dzx && orbital1 == py) ||
            (orbital1 == dyz && orbital2 == px) ||
            (orbital2 == dyz && orbital1 == px)  ){
      value = atom.GetZindoF0sdLower()
             -atom.GetZindoF2pdLower()*4.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == dzz) ||
            (orbital2 == dxxyy && orbital1 == dzz) || 
            (orbital1 == dxy && orbital2 == dzz) ||
            (orbital2 == dxy && orbital1 == dzz)  ){
      value = atom.GetZindoF0ddLower()
             -atom.GetZindoF2ddLower()*4.0
             +atom.GetZindoF4ddLower()*6.0;
   }
   else if( (orbital1 == dxy && orbital2 == dxxyy) ||
            (orbital2 == dxy && orbital1 == dxxyy)  ){
      value = atom.GetZindoF0ddLower()
             +atom.GetZindoF2ddLower()*4.0
             -atom.GetZindoF4ddLower()*34.0;
   }
   else if( (orbital1 == dzx && orbital2 == dzz) ||
            (orbital2 == dzx && orbital1 == dzz) || 
            (orbital1 == dyz && orbital2 == dzz) ||
            (orbital2 == dyz && orbital1 == dzz)  ){
      value = atom.GetZindoF0ddLower()
             +atom.GetZindoF2ddLower()*2.0
             -atom.GetZindoF4ddLower()*24.0;
   }
   else if( (orbital1 == dzx && orbital2 == dxxyy) ||
            (orbital2 == dzx && orbital1 == dxxyy) || 
            (orbital1 == dzx && orbital2 == dxy) || 
            (orbital2 == dzx && orbital1 == dxy) || 
            (orbital1 == dyz && orbital2 == dxxyy) || 
            (orbital2 == dyz && orbital1 == dxxyy) || 
            (orbital1 == dyz && orbital2 == dxy) || 
            (orbital2 == dyz && orbital1 == dxy) || 
            (orbital1 == dyz && orbital2 == dzx) || 
            (orbital2 == dyz && orbital1 == dzx) ){
      value = atom.GetZindoF0ddLower()
             -atom.GetZindoF2ddLower()*2.0
             -atom.GetZindoF4ddLower()*4.0;
   }
   */
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

// Apendix in [BZ_1972]
// ZINDO Exchange Interaction
double ZindoS::GetExchangeInt(OrbitalType orbital1, OrbitalType orbital2, const Atom& atom) const{

   double value=0.0;

   if( orbital1 == orbital2){
      value = this->GetCoulombInt(orbital1, orbital2, atom);
   }   
   else if( orbital1 == s && (orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom.GetZindoG1spLower();
   }   
   else if( orbital2 == s && (orbital1 == px || orbital1 == py || orbital1 == pz ) ){
      value = atom.GetZindoG1spLower();
   }   
   else if( (orbital1 != orbital2) 
         && ( orbital1 == px || orbital1 == py || orbital1 == pz )
         && ( orbital2 == px || orbital2 == py || orbital2 == pz ) ){
      value = atom.GetZindoF2ppLower()*3.0;
   }
   // ToDo: There are bugs for d-orbitals.
   /*
   else if( (orbital1 == s) && (orbital2 == dxy || 
                                orbital2 == dyz || 
                                orbital2 == dzz || 
                                orbital2 == dzx || 
                                orbital2 == dxxyy ) ){
      value = atom->GetZindoG2sdLower();
   }   
   else if( (orbital2 == s) && (orbital1 == dxy || 
                                orbital1 == dyz || 
                                orbital1 == dzz || 
                                orbital1 == dzx || 
                                orbital1 == dxxyy ) ){
      value = atom->GetZindoG2sdLower();
   }   
   else if( (orbital1 == px && orbital2 == dzz) ||
            (orbital2 == px && orbital1 == dzz) ||
            (orbital1 == py && orbital2 == dzz) ||
            (orbital2 == py && orbital1 == dzz) ){
      value = atom->GetZindoG1pdLower()
             +atom->GetZindoG3pdLower()*18.0;
   }
   else if( (orbital1 == px && orbital2 == dxxyy) ||
            (orbital2 == px && orbital1 == dxxyy) ||
            (orbital1 == px && orbital2 == dxy) ||
            (orbital2 == px && orbital1 == dxy) ||
            (orbital1 == px && orbital2 == dzx) ||
            (orbital2 == px && orbital1 == dzx) ||
            (orbital1 == py && orbital2 == dxxyy) ||
            (orbital2 == py && orbital1 == dxxyy) ||
            (orbital1 == py && orbital2 == dxy) ||
            (orbital2 == py && orbital1 == dxy) ||
            (orbital1 == py && orbital2 == dyz) ||
            (orbital2 == py && orbital1 == dyz) ||
            (orbital1 == pz && orbital2 == dzx) ||
            (orbital2 == pz && orbital1 == dzx) ||
            (orbital1 == pz && orbital2 == dyz) ||
            (orbital2 == pz && orbital1 == dyz) ){
      value = atom->GetZindoG1pdLower()*3.0
             +atom->GetZindoG3pdLower()*24.0;
   }
   else if( (orbital1 == px && orbital2 == dyz) ||
            (orbital2 == px && orbital1 == dyz) ||
            (orbital1 == py && orbital2 == dzx) ||
            (orbital2 == py && orbital1 == dzx) ||
            (orbital1 == pz && orbital2 == dxxyy) ||
            (orbital2 == pz && orbital1 == dxxyy) ||
            (orbital1 == pz && orbital2 == dxy) ||
            (orbital2 == pz && orbital1 == dxy) ){
      value = atom->GetZindoG3pdLower()*15.0;
   }
   else if( (orbital1 == pz && orbital2 == dzz) ||
            (orbital2 == pz && orbital1 == dzz) ){
      value = atom->GetZindoG1pdLower()*4.0
             +atom->GetZindoG3pdLower()*27.0;
   }
   else if( (orbital1 == dzz && orbital2 == dxxyy) ||
            (orbital2 == dzz && orbital1 == dxxyy) ||
            (orbital1 == dzz && orbital2 == dxy) ||
            (orbital2 == dzz && orbital1 == dxy) ){
      value = atom->GetZindoF2ddLower()*4.0
             +atom->GetZindoF4ddLower()*15.0;
   }
   else if( (orbital1 == dzz && orbital2 == dzx) ||
            (orbital2 == dzz && orbital1 == dzx) ||
            (orbital1 == dzz && orbital2 == dyz) ||
            (orbital2 == dzz && orbital1 == dyz) ){
      value = atom->GetZindoF2ddLower()
             +atom->GetZindoF4ddLower()*30.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == dxy) ||
            (orbital2 == dxxyy && orbital1 == dxy) ){
      value = atom->GetZindoF4ddLower()*35.0;
   }
   else if( (orbital1 == dxxyy && orbital2 == dzx) ||
            (orbital2 == dxxyy && orbital1 == dzx) ||
            (orbital1 == dxxyy && orbital2 == dyz) ||
            (orbital2 == dxxyy && orbital1 == dyz) ||
            (orbital1 == dxy && orbital2 == dzx) ||
            (orbital2 == dxy && orbital1 == dzx) ||
            (orbital1 == dxy && orbital2 == dyz) ||
            (orbital2 == dxy && orbital1 == dyz) ||
            (orbital1 == dzx && orbital2 == dyz) ||
            (orbital2 == dzx && orbital1 == dyz) ){
      value = atom->GetZindoF2ddLower()*3.0
             +atom->GetZindoF4ddLower()*20.0;
   }
   */
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

// ref. [MN_1957] and (5a) in [AEZ_1986]
double ZindoS::GetNishimotoMatagaTwoEleInt(const Atom& atomA, OrbitalType orbitalA, 
                                           const Atom& atomB, OrbitalType orbitalB) const{
   double r = this->molecule->GetDistanceAtoms(atomA, atomB);
   double gammaAA;
   if(orbitalA == s || 
      orbitalA == px ||
      orbitalA == py ||
      orbitalA == pz ){
      gammaAA = atomA.GetZindoF0ss();
   }
   /*
   else if(orbitalA == dxy ||
           orbitalA == dyz ||
           orbitalA == dzz ||
           orbitalA == dzx ||
           orbitalA == dxxyy ){
      gammaAA = atomA->GetZindoF0dd();
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageNishimotoMataga;
      ss << this->errorMessageAtomType << AtomTypeStr(atomA.GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalA) << "\n";
      throw MolDSException(ss.str());
   }   

   double gammaBB;
   if(orbitalB == s || 
      orbitalB == px ||
      orbitalB == py ||
      orbitalB == pz ){
      gammaBB = atomB.GetZindoF0ss();
   }
   /*
   else if(orbitalB == dxy ||
           orbitalB == dyz ||
           orbitalB == dzz ||
           orbitalB == dzx ||
           orbitalB == dxxyy ){
      gammaBB = atomB->GetZindoF0dd();
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageNishimotoMataga;
      ss << this->errorMessageAtomType << AtomTypeStr(atomB.GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalB) << "\n";
      throw MolDSException(ss.str());
   }  

   return this->nishimotoMatagaParamA/( r+this->nishimotoMatagaParamB/(gammaAA+gammaBB) );

}

// First derivative of Nishimoto-Mataga related to the coordinate of atom A.
// For Nishimoto-Mataga, See ZindoS::GetNishimotoMatagaTwoEleInt 
// or ref. [MN_1957] and (5a) in [AEZ_1986]
double ZindoS::GetNishimotoMatagaTwoEleIntFirstDerivative(const Atom& atomA, 
                                                          OrbitalType orbitalA, 
                                                          const Atom& atomB, 
                                                          OrbitalType orbitalB,
                                                          CartesianType axisA) const{
   double gammaAA;
   if(orbitalA == s || 
      orbitalA == px ||
      orbitalA == py ||
      orbitalA == pz ){
      gammaAA = atomA.GetZindoF0ss();
   }
   /*
   else if(orbitalA == dxy ||
           orbitalA == dyz ||
           orbitalA == dzz ||
           orbitalA == dzx ||
           orbitalA == dxxyy ){
      gammaAA = atomA.GetZindoF0dd();
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageNishimotoMataga;
      ss << this->errorMessageAtomType << AtomTypeStr(atomA.GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalA) << "\n";
      throw MolDSException(ss.str());
   }   

   double gammaBB;
   if(orbitalB == s || 
      orbitalB == px ||
      orbitalB == py ||
      orbitalB == pz ){
      gammaBB = atomB.GetZindoF0ss();
   }
   /*
   else if(orbitalB == dxy ||
           orbitalB == dyz ||
           orbitalB == dzz ||
           orbitalB == dzx ||
           orbitalB == dxxyy ){
      gammaBB = atomB.GetZindoF0dd();
   }
   */
   else{
      stringstream ss;
      ss << this->errorMessageNishimotoMataga;
      ss << this->errorMessageAtomType << AtomTypeStr(atomB.GetAtomType()) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalB) << "\n";
      throw MolDSException(ss.str());
   }  

   double r = this->molecule->GetDistanceAtoms(atomA, atomB);
   double dCartesian = atomA.GetXyz()[axisA] - atomB.GetXyz()[axisA];
   double value = -1.0*dCartesian/r;
   value *= this->nishimotoMatagaParamA;
   value *= pow( r+this->nishimotoMatagaParamB/(gammaAA+gammaBB) ,-2.0);
   return value;
}

void ZindoS::CalcDiatomicOverlapInDiatomicFrame(double** diatomicOverlap, 
                                                const Atom& atomA, 
                                                const Atom& atomB) const{
   MolDS_cndo::Cndo2::CalcDiatomicOverlapInDiatomicFrame(diatomicOverlap, atomA, atomB);

   // see (4f) in [AEZ_1986]
   diatomicOverlap[pz][pz] *= this->overlapCorrectionSigma;
   diatomicOverlap[py][py] *= this->overlapCorrectionPi;
   diatomicOverlap[px][px] *= this->overlapCorrectionPi;
   
   for(int i=0;i<OrbitalType_end;i++){
      for(int j=0;j<OrbitalType_end;j++){
         //this->OutputLog((boost::format("diatomicOverlap[%d][%d]=%lf\n") % i % j % diatomicOverlap[i][j]).str());
      }
   }
   
}

// First derivative of (B.40) in J. A. Pople book with bond correction.
void ZindoS::CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(double** diatomicOverlapDeri, 
                                                               const Atom& atomA, 
                                                               const Atom& atomB) const{

   MolDS_cndo::Cndo2::CalcDiatomicOverlapFirstDerivativeInDiatomicFrame(
                      diatomicOverlapDeri,atomA, atomB);

   // see (4f) in [AEZ_1986] like as overlap integlral
   diatomicOverlapDeri[pz][pz] *= this->overlapCorrectionSigma;
   diatomicOverlapDeri[py][py] *= this->overlapCorrectionPi;
   diatomicOverlapDeri[px][px] *= this->overlapCorrectionPi;
}

// Second derivative of (B.40) in J. A. Pople book with bond correction.
void ZindoS::CalcDiatomicOverlapSecondDerivativeInDiatomicFrame(double** diatomicOverlapSecondDeri, 
                                                                const Atom& atomA, 
                                                                const Atom& atomB) const{

   MolDS_cndo::Cndo2::CalcDiatomicOverlapSecondDerivativeInDiatomicFrame(
                      diatomicOverlapSecondDeri,atomA, atomB);

   // see (4f) in [AEZ_1986] like as overlap integlral
   diatomicOverlapSecondDeri[pz][pz] *= this->overlapCorrectionSigma;
   diatomicOverlapSecondDeri[py][py] *= this->overlapCorrectionPi;
   diatomicOverlapSecondDeri[px][px] *= this->overlapCorrectionPi;
}

// The order of mol, moJ, moK, moL is consistent with Eq. (9) in [RZ_1973]
double ZindoS::GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                           const Molecule& molecule, 
                                           double const* const* fockMatrix, 
                                           double const* const* gammaAB) const{
   double value = 0.0;
   int firstAOIndexA;
   int firstAOIndexB;
   int numberAOsA;
   int numberAOsB;
   double gamma;
   double exchange;
   double coulomb;
   OrbitalType orbitalMu;
   OrbitalType orbitalNu;

   for(int A=0; A<molecule.GetNumberAtoms(); A++){
      const Atom& atomA = *molecule.GetAtom(A);
      firstAOIndexA = atomA.GetFirstAOIndex();
      numberAOsA = atomA.GetValenceSize();

      for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
         orbitalMu = atomA.GetValence(mu-firstAOIndexA);

         // CNDO term
         for(int B=A; B<molecule.GetNumberAtoms(); B++){
            const Atom& atomB = *molecule.GetAtom(B);
            firstAOIndexB = atomB.GetFirstAOIndex();
            numberAOsB = atomB.GetValenceSize();

            for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
               orbitalNu = atomB.GetValence(nu-firstAOIndexB);

               if(A<B){
                  gamma = this->GetNishimotoMatagaTwoEleInt(atomA, 
                                                            orbitalMu, 
                                                            atomB, 
                                                            orbitalNu);
                  value += gamma
                          *fockMatrix[moI][mu]
                          *fockMatrix[moJ][mu]
                          *fockMatrix[moK][nu]
                          *fockMatrix[moL][nu];
                  value += gamma
                          *fockMatrix[moI][nu]
                          *fockMatrix[moJ][nu]
                          *fockMatrix[moK][mu]
                          *fockMatrix[moL][mu];
               }
               else{
                  gamma = atomA.GetZindoF0ss();
                  value += gamma
                          *fockMatrix[moI][mu]
                          *fockMatrix[moJ][mu]
                          *fockMatrix[moK][nu]
                          *fockMatrix[moL][nu];
               }  
            }
         }

         // Aditional term for INDO or ZIND/S, see Eq. (10) in [RZ_1973]
         for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
            orbitalNu = atomA.GetValence(nu-firstAOIndexA);

            if(mu!=nu){
               exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
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

            coulomb = this->GetCoulombInt(orbitalMu, orbitalNu, atomA);

            if( (orbitalMu == s || orbitalMu == px || orbitalMu == py || orbitalMu == pz) &&
                (orbitalNu == s || orbitalNu == px || orbitalNu == py || orbitalNu == pz) ){
                  gamma = atomA.GetZindoF0ss();
            }
            else{
               stringstream ss;
               ss << this->errorMessageMolecularIntegralElement;
               ss << this->errorMessageAtomType << AtomTypeStr(atomA.GetAtomType()) << "\n";
               ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalMu) << "\n";
               ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalNu) << "\n";
               throw MolDSException(ss.str());
            }   

            value += (coulomb-gamma)
                    *fockMatrix[moI][mu]
                    *fockMatrix[moJ][mu]
                    *fockMatrix[moK][nu]
                    *fockMatrix[moL][nu];
         }
      }
   }

   return value;
}

void ZindoS::DoCIS(){
   this->OutputLog(this->messageStartCIS);
   double ompStartTime = omp_get_wtime();

   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   // malloc or initialize  CIS matrix
   if(this->matrixCIS == NULL){
      this->matrixCISdimension = numberActiveOcc*numberActiveVir;
      MallocerFreer::GetInstance()->Malloc<double>(&this->matrixCIS, 
                                                   this->matrixCISdimension, 
                                                   this->matrixCISdimension);
   }
   else{
      MallocerFreer::GetInstance()->Initialize<double>(this->matrixCIS, 
                                                       this->matrixCISdimension, 
                                                       this->matrixCISdimension);
   }
   // malloc or initialize CIS eigen vector
   if(this->excitedEnergies == NULL){
      MallocerFreer::GetInstance()->Malloc<double>(&this->excitedEnergies,
                                                   this->matrixCISdimension);
   }
   else{
      MallocerFreer::GetInstance()->Initialize<double>(this->excitedEnergies, 
                                                       this->matrixCISdimension);
   }
   // calculate CIS matrix
   this->CalcCISMatrix(matrixCIS);
   // calculate excited energies
   if(Parameters::GetInstance()->IsDavidsonCIS()){
      this->DoCISDavidson();
   }
   else{
      this->DoCISDirect();
   }
   this->CalcCISProperties();
   this->OutputCISResults();

   double ompEndTime = omp_get_wtime();
   this->OutputLog((boost::format("%s%lf%s\n%s") % this->messageOmpElapsedTimeCIS.c_str()
                                                 % (ompEndTime - ompStartTime)
                                                 % this->messageUnitSec.c_str()
                                                 % this->messageDoneCIS.c_str() ).str());
}

void ZindoS::CalcCISProperties(){

   // dipole moments of excited states
   this->CalcElectronicDipoleMomentsExcitedState(this->electronicTransitionDipoleMoments,
                                                 this->fockMatrix,
                                                 this->matrixCIS,
                                                 this->cartesianMatrix,
                                                 *this->molecule, 
                                                 this->orbitalElectronPopulation,
                                                 this->overlap);
   
   // transition dipole moment
   this->CalcElectronicTransitionDipoleMoments(this->electronicTransitionDipoleMoments,
                                               this->fockMatrix,
                                               this->matrixCIS,
                                               this->cartesianMatrix,
                                               *this->molecule, 
                                               this->orbitalElectronPopulation,
                                               this->overlap);


   // free exciton energies
   this->CalcFreeExcitonEnergies(&this->freeExcitonEnergiesCIS, 
                                 *this->molecule, 
                                 this->energiesMO, 
                                 this->matrixCIS,
                                 this->matrixCISdimension);
}

void ZindoS::CalcElectronicDipoleMomentsExcitedState(double*** electronicTransitionDipoleMoments,
                                                     double const* const* fockMatrix,
                                                     double const* const* matrixCIS,
                                                     double const* const* const* cartesianMatrix,
                                                     const MolDS_base::Molecule& molecule, 
                                                     double const* const* orbitalElectronPopulation,
                                                     double const* const* overlap) const{
   int groundState = 0;
   // dipole moment of excited states
   for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
      int excitedState = k+1; // (k+1)-th excited state
      for(int axis=0; axis<CartesianType_end; axis++){
         electronicTransitionDipoleMoments[excitedState][excitedState][axis] = this->GetElectronicTransitionDipoleMoment(
                                                                                     excitedState,
                                                                                     excitedState,
                                                                                     static_cast<CartesianType>(axis),
                                                                                     fockMatrix,
                                                                                     matrixCIS,
                                                                                     cartesianMatrix,
                                                                                     molecule,
                                                                                     orbitalElectronPopulation,
                                                                                     overlap,
                                                                                     electronicTransitionDipoleMoments[groundState][groundState]);
      }
   }
}

void ZindoS::CalcElectronicTransitionDipoleMoments(double*** electronicTransitionDipoleMoments,
                                                   double const* const* fockMatrix,
                                                   double const* const* matrixCIS,
                                                   double const* const* const* cartesianMatrix,
                                                   const MolDS_base::Molecule& molecule, 
                                                   double const* const* orbitalElectronPopulation,
                                                   double const* const* overlap) const{
   int groundState = 0;
   // transition dipole moments from ground state to excited states
   for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
      int excitedState = k+1; // (k+1)-th excited state
      for(int axis=0; axis<CartesianType_end; axis++){
         this->electronicTransitionDipoleMoments[excitedState][groundState][axis] = this->GetElectronicTransitionDipoleMoment(
                                                                                          excitedState,
                                                                                          groundState,
                                                                                          static_cast<CartesianType>(axis),
                                                                                          this->fockMatrix,
                                                                                          this->matrixCIS,
                                                                                          this->cartesianMatrix,
                                                                                          *this->molecule,
                                                                                          this->orbitalElectronPopulation,
                                                                                          this->overlap,
                                                                                          this->electronicTransitionDipoleMoments[groundState][groundState]);
      }
   }

   if(Parameters::GetInstance()->RequiresAllTransitionDipoleMomentsCIS()){
      // transition dipole moments from excited states to ground state
      for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
         int excitedState = k+1; // (k+1)-th excited state
         for(int axis=0; axis<CartesianType_end; axis++){
            this->electronicTransitionDipoleMoments[groundState][excitedState][axis] = this->GetElectronicTransitionDipoleMoment(
                                                                                             groundState,
                                                                                             excitedState,
                                                                                             static_cast<CartesianType>(axis),
                                                                                             this->fockMatrix,
                                                                                             this->matrixCIS,
                                                                                             this->cartesianMatrix,
                                                                                             *this->molecule,
                                                                                             this->orbitalElectronPopulation,
                                                                                             this->overlap,
                                                                                             this->electronicTransitionDipoleMoments[groundState][groundState]);
         }
      }

      // transition dipole moments between excited states
      for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
         int departureExcitedState = k+1; // (k+1)-th excited state
         for(int l=0; l<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); l++){
            if(k != l){
               int destinationExcitedState = l+1; // (l+1)-th excited state
               for(int axis=0; axis<CartesianType_end; axis++){
                  this->electronicTransitionDipoleMoments[destinationExcitedState][departureExcitedState][axis] 
                     = this->GetElectronicTransitionDipoleMoment(destinationExcitedState,
                                                                 departureExcitedState,
                                                                 static_cast<CartesianType>(axis),
                                                                 this->fockMatrix,
                                                                 this->matrixCIS,
                                                                 this->cartesianMatrix,
                                                                 *this->molecule,
                                                                 this->orbitalElectronPopulation,
                                                                 this->overlap,
                                                                 this->electronicTransitionDipoleMoments[groundState][groundState]);
               }
            }
         }
      }
   }

}

double ZindoS::GetElectronicTransitionDipoleMoment(int to, int from, CartesianType axis,
                                                   double const* const* fockMatrix,
                                                   double const* const* matrixCIS,
                                                   double const* const* const* cartesianMatrix,
                                                   const MolDS_base::Molecule& molecule, 
                                                   double const* const* orbitalElectronPopulation,
                                                   double const* const* overlap,
                                                   double const* groundStateDipole) const{
   double value = 0.0;
   int groundState = 0;

   if(Parameters::GetInstance()->GetNumberExcitedStatesCIS() < from ||
      Parameters::GetInstance()->GetNumberExcitedStatesCIS() < to ){
      stringstream ss;
      ss << this->errorMessageGetElectronicTransitionDipoleMomentBadState;
      ss << this->errorMessageFromState << from << endl;
      ss << this->errorMessageToState << to << endl;
      ss << this->errorMessageCartesianType << CartesianTypeStr(axis) << endl;
      throw MolDSException(ss.str());
   }

   if(from != groundState && to != groundState && from != to){
      // transition dipole moment between excited states
      for(int l=0; l<this->matrixCISdimension; l++){
         double temp = 0.0;
         // single excitation from I-th (occupied)MO to A-th (virtual)MO
         int moI = this->GetActiveOccIndex(molecule, l);
         int moA = this->GetActiveVirIndex(molecule, l);
         for(int mu=0; mu<molecule.GetTotalNumberAOs(); mu++){
            for(int nu=0; nu<molecule.GetTotalNumberAOs(); nu++){
               temp += (-1.0*fockMatrix[moI][mu]*fockMatrix[moI][nu] + fockMatrix[moA][mu]*fockMatrix[moA][nu])
                       *(cartesianMatrix[mu][nu][axis] - molecule.GetXyzCOC()[axis]*overlap[mu][nu]);
            }
         }
         value += matrixCIS[from-1][l]*matrixCIS[to-1][l]*temp;
      }
   }
   else if(from == groundState && to != groundState){
      // transition dipole moment from the ground to excited states
      for(int l=0; l<this->matrixCISdimension; l++){
         double temp = 0.0;
         // single excitation from I-th (occupied)MO to A-th (virtual)MO
         int moI = this->GetActiveOccIndex(molecule, l);
         int moA = this->GetActiveVirIndex(molecule, l);
         for(int mu=0; mu<molecule.GetTotalNumberAOs(); mu++){
            for(int nu=0; nu<molecule.GetTotalNumberAOs(); nu++){
               temp += fockMatrix[moA][mu]
                       *(cartesianMatrix[mu][nu][axis] - molecule.GetXyzCOC()[axis]*overlap[mu][nu])
                       *fockMatrix[moI][nu];
            }
         }
         value += this->matrixCIS[to-1][l]*sqrt(2.0)*temp;
      }
   }
   else if(from != groundState && to == groundState){
      // transition dipole moment from the excited to ground states
      for(int l=0; l<this->matrixCISdimension; l++){
         double temp = 0.0;
         // single excitation from I-th (occupied)MO to A-th (virtual)MO
         int moI = this->GetActiveOccIndex(molecule, l);
         int moA = this->GetActiveVirIndex(molecule, l);
         for(int mu=0; mu<molecule.GetTotalNumberAOs(); mu++){
            for(int nu=0; nu<molecule.GetTotalNumberAOs(); nu++){
               temp += fockMatrix[moI][mu]
                       *(cartesianMatrix[mu][nu][axis] - molecule.GetXyzCOC()[axis]*overlap[mu][nu])
                       *fockMatrix[moA][nu];
            }
         }
         value += matrixCIS[from-1][l]*sqrt(2.0)*temp;
      }
   }
   else if(from == to){
      if(from != groundState){
         // dipole moment of the excited state. It is needed that the dipole of ground state has been already calculated!!
         for(int l=0; l<this->matrixCISdimension; l++){
            double temp = 0.0;
            // single excitation from I-th (occupied)MO to A-th (virtual)MO
            int moI = this->GetActiveOccIndex(molecule, l);
            int moA = this->GetActiveVirIndex(molecule, l);
            for(int mu=0; mu<molecule.GetTotalNumberAOs(); mu++){
               for(int nu=0; nu<molecule.GetTotalNumberAOs(); nu++){
                  temp += (-1.0*fockMatrix[moI][mu]*fockMatrix[moI][nu] + fockMatrix[moA][mu]*fockMatrix[moA][nu])
                          *(cartesianMatrix[mu][nu][axis] - molecule.GetXyzCOC()[axis]*overlap[mu][nu]);
               }
            }
            value += matrixCIS[from-1][l]*matrixCIS[to-1][l]*temp;
         }
         value += groundStateDipole[axis];
      }
      else{
         // dipole moment of the ground state
         value = Cndo2::GetElectronicTransitionDipoleMoment(to,
                                                            from,
                                                            axis,
                                                            fockMatrix,
                                                            matrixCIS,
                                                            cartesianMatrix,
                                                            molecule,
                                                            orbitalElectronPopulation,
                                                            overlap,
                                                            NULL);
         }
   }

   return value;
}

void ZindoS::CalcFreeExcitonEnergies(double** freeExcitonEnergiesCIS, 
                                     const Molecule& molecule, 
                                     double const* energiesMO, 
                                     double const* const* matrixCIS,
                                     int matrixCISdimension) const{
   if(Parameters::GetInstance()->RequiresExcitonEnergiesCIS()){
      // malloc or initialize free exciton energies
      if(*freeExcitonEnergiesCIS == NULL){
         MallocerFreer::GetInstance()->Malloc<double>(freeExcitonEnergiesCIS,
                                                      matrixCISdimension);
      }
      else{
         MallocerFreer::GetInstance()->Initialize<double>(*freeExcitonEnergiesCIS, 
                                                          matrixCISdimension);
      }
      // clac free exciton energies
      for(int k=0; k<matrixCISdimension; k++){
         double value = 0.0;
         for(int l=0; l<matrixCISdimension; l++){
            // single excitation from I-th (occupied)MO to A-th (virtual)MO
            int moI = this->GetActiveOccIndex(molecule, l);
            int moA = this->GetActiveVirIndex(molecule, l);
            value += pow(matrixCIS[k][l],2.0)*(energiesMO[moA] - energiesMO[moI]);
         }
         (*freeExcitonEnergiesCIS)[k] = value;
      }
   }
}

void ZindoS::OutputCISResults() const{
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();

   // output cis eigen energies
   this->OutputLog(this->messageExcitedStatesEnergiesTitle);
   double eV2AU = Parameters::GetInstance()->GetEV2AU();
   for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
      this->OutputLog((boost::format("%s\t%d\t%e\t%e\t") 
         % this->messageExcitedStatesEnergies
         % (k+1) 
         % this->excitedEnergies[k]
         % (this->excitedEnergies[k]/eV2AU)).str());

      // sort eigen vector coefficeits of CIS and output
      vector<CISEigenVectorCoefficient> cisEigenVectorCoefficients;
      this->SortCISEigenVectorCoefficients(&cisEigenVectorCoefficients, this->matrixCIS[k]);
      for(int l=0; l<Parameters::GetInstance()->GetNumberPrintCoefficientsCIS(); l++){
         this->OutputLog((boost::format("%e (%d -> %d)\t") % cisEigenVectorCoefficients[l].coefficient
                                                           % cisEigenVectorCoefficients[l].occIndex
                                                           % cisEigenVectorCoefficients[l].virIndex).str());
      }
      this->OutputLog("\n");
   }
   this->OutputLog("\n");

   // output dipole moment
   this->OutputCISDipole();

   // output transition dipole moment
   this->OutputCISTransitionDipole();

   // output exciton energies
   if(Parameters::GetInstance()->RequiresExcitonEnergiesCIS()){
      this->OutputLog(this->messageExcitonEnergiesCIS);
      this->OutputLog(this->messageExcitonEnergiesCISTitle);
      for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
         this->OutputLog((boost::format("%s\t%d\t%e\t%e\t%e\t%e\n") 
                           % this->messageExcitonEnergiesShortCIS
                           % (k+1) 
                           %  this->freeExcitonEnergiesCIS[k]
                           % (this->freeExcitonEnergiesCIS[k]/eV2AU)
                           % (this->excitedEnergies[k]-this->freeExcitonEnergiesCIS[k])
                           %((this->excitedEnergies[k]-this->freeExcitonEnergiesCIS[k])/eV2AU)).str());
      }
   }
   this->OutputLog("\n");

   // output Hole density
   if(Parameters::GetInstance()->RequiresHolePlot()){
      MolDS_base_loggers::DensityLogger* holeDensityLogger = new MolDS_base_loggers::HoleDensityLogger(
                                                                                     *this->molecule, 
                                                                                     this->fockMatrix, 
                                                                                     this->matrixCIS, 
                                                                                     this->theory);
      holeDensityLogger->DrawDensity(*(Parameters::GetInstance()->GetElecIndecesHolePlot()));
      delete holeDensityLogger;
   }

   // output particle density
   if(Parameters::GetInstance()->RequiresParticlePlot()){
      MolDS_base_loggers::DensityLogger* particleDensityLogger = new MolDS_base_loggers::ParticleDensityLogger(
                                                                                         *this->molecule, 
                                                                                         this->fockMatrix, 
                                                                                         this->matrixCIS, 
                                                                                         this->theory);
      particleDensityLogger->DrawDensity(*(Parameters::GetInstance()->GetElecIndecesParticlePlot()));
      delete particleDensityLogger;
   }
}

// this module output ground and excited state dipole moments.
void ZindoS::OutputCISDipole() const{
   double debye2AU = Parameters::GetInstance()->GetDebye2AU();

   // total dipole (electronic + core dipole)
   this->OutputLog(this->messageTotalDipoleMomentsTitle);
   for(int k=0; k<=Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
      double magnitude = 0.0; 
      double temp = 0.0;
      temp += pow(this->electronicTransitionDipoleMoments[k][k][XAxis]+this->coreDipoleMoment[XAxis],2.0);
      temp += pow(this->electronicTransitionDipoleMoments[k][k][YAxis]+this->coreDipoleMoment[YAxis],2.0);
      temp += pow(this->electronicTransitionDipoleMoments[k][k][ZAxis]+this->coreDipoleMoment[ZAxis],2.0);
      magnitude = sqrt(temp);
      this->OutputLog((boost::format("\t%s\t%d\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n") 
         % this->messageTotalDipoleMoment
         % k
         % (this->electronicTransitionDipoleMoments[k][k][XAxis]+this->coreDipoleMoment[XAxis])
         % (this->electronicTransitionDipoleMoments[k][k][YAxis]+this->coreDipoleMoment[YAxis])
         % (this->electronicTransitionDipoleMoments[k][k][ZAxis]+this->coreDipoleMoment[ZAxis])
         % magnitude
         % ((this->electronicTransitionDipoleMoments[k][k][XAxis]+this->coreDipoleMoment[XAxis])/debye2AU)
         % ((this->electronicTransitionDipoleMoments[k][k][YAxis]+this->coreDipoleMoment[YAxis])/debye2AU)
         % ((this->electronicTransitionDipoleMoments[k][k][ZAxis]+this->coreDipoleMoment[ZAxis])/debye2AU)
         % (magnitude/debye2AU)).str());
   }
   this->OutputLog("\n");

   // electronic dipole
   this->OutputLog(this->messageElectronicDipoleMomentsTitle);
   for(int k=0; k<=Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
      double magnitude = 0.0; 
      double temp = 0.0;
      temp += pow(this->electronicTransitionDipoleMoments[k][k][XAxis],2.0);
      temp += pow(this->electronicTransitionDipoleMoments[k][k][YAxis],2.0);
      temp += pow(this->electronicTransitionDipoleMoments[k][k][ZAxis],2.0);
      magnitude = sqrt(temp);
      this->OutputLog((boost::format("\t%s\t%d\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n") 
         % this->messageElectronicDipoleMoment
         % k
         % (this->electronicTransitionDipoleMoments[k][k][XAxis])
         % (this->electronicTransitionDipoleMoments[k][k][YAxis])
         % (this->electronicTransitionDipoleMoments[k][k][ZAxis])
         % magnitude
         % (this->electronicTransitionDipoleMoments[k][k][XAxis]/debye2AU)
         % (this->electronicTransitionDipoleMoments[k][k][YAxis]/debye2AU)
         % (this->electronicTransitionDipoleMoments[k][k][ZAxis]/debye2AU)
         % (magnitude/debye2AU)).str());
   }
   this->OutputLog("\n");

}

void ZindoS::OutputCISTransitionDipole() const{
   double debye2AU = Parameters::GetInstance()->GetDebye2AU();
   int groundState = 0;
   // electronic dipole
   this->OutputLog(this->messageTransitionDipoleMomentsTitle);
   for(int from=0; from<=Parameters::GetInstance()->GetNumberExcitedStatesCIS(); from++){
      if(groundState < from && !Parameters::GetInstance()->RequiresAllTransitionDipoleMomentsCIS()){
         break;
      }
      for(int to=0; to<=Parameters::GetInstance()->GetNumberExcitedStatesCIS(); to++){
         if(from != to){
            double magnitude = 0.0; 
            double temp = 0.0;
            temp += pow(this->electronicTransitionDipoleMoments[to][from][XAxis],2.0);
            temp += pow(this->electronicTransitionDipoleMoments[to][from][YAxis],2.0);
            temp += pow(this->electronicTransitionDipoleMoments[to][from][ZAxis],2.0);
            magnitude = sqrt(temp);
            this->OutputLog((boost::format("\t%s\t%d -> %d\t\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n") 
               % this->messageTransitionDipoleMoment
               % from
               % to
               % (this->electronicTransitionDipoleMoments[to][from][XAxis])
               % (this->electronicTransitionDipoleMoments[to][from][YAxis])
               % (this->electronicTransitionDipoleMoments[to][from][ZAxis])
               % magnitude
               % (this->electronicTransitionDipoleMoments[to][from][XAxis]/debye2AU)
               % (this->electronicTransitionDipoleMoments[to][from][YAxis]/debye2AU)
               % (this->electronicTransitionDipoleMoments[to][from][ZAxis]/debye2AU)
               % (magnitude/debye2AU)).str());
         }
      }
   }
   this->OutputLog("\n");

}
void ZindoS::SortCISEigenVectorCoefficients(vector<CISEigenVectorCoefficient>* cisEigenVectorCoefficients,
                                            double* cisEigenVector) const{
   for(int l=0; l<this->matrixCISdimension; l++){
      // single excitation from I-th (occupied)MO to A-th (virtual)MO
      int moI = this->GetActiveOccIndex(*this->molecule, l);
      int moA = this->GetActiveVirIndex(*this->molecule, l);
      CISEigenVectorCoefficient cisEigenVectorCoefficient = {cisEigenVector[l], moI, moA, k};
      cisEigenVectorCoefficients->push_back(cisEigenVectorCoefficient);
   }
   sort(cisEigenVectorCoefficients->begin(), 
        cisEigenVectorCoefficients->end(), 
        MoreCISEigenVectorCoefficient());
}

void ZindoS::SortSingleExcitationSlaterDeterminants(vector<MoEnergyGap>* moEnergyGaps) const{
   for(int k=0; k<this->matrixCISdimension; k++){
      // single excitation from I-th (occupied)MO to A-th (virtual)MO
      int moI = this->GetActiveOccIndex(*this->molecule, k);
      int moA = this->GetActiveVirIndex(*this->molecule, k);
      MoEnergyGap moEnergyGap = {this->energiesMO[moA]-this->energiesMO[moI], moI, moA, k};
      moEnergyGaps->push_back(moEnergyGap);
   }
   sort(moEnergyGaps->begin(), moEnergyGaps->end(), LessMoEnergyGap());
}

// This method is used for Davidson
void ZindoS::CalcRitzVector(double* ritzVector, 
                            double const* const* expansionVectors, 
                            double const* const* interactionMatrix, 
                            int interactionMatrixDimension, 
                            int ritzVectorIndex) const{
   for(int j=0; j<this->matrixCISdimension; j++){
      ritzVector[j] = 0.0;
      for(int k=0; k<interactionMatrixDimension; k++){
         ritzVector[j] += expansionVectors[j][k]*interactionMatrix[ritzVectorIndex][k];
      }
   }
}

// This method is used for Davidson
void ZindoS::CalcResidualVectorAndNorm(double* residualVector, 
                                       double* norm, 
                                       double const* ritzVector, 
                                       double const* interactionEigenEnergies, 
                                       int residualVectorIndex) const{
   double sqNorm = 0.0;
   for(int j=0; j<this->matrixCISdimension; j++){
      residualVector[j] = interactionEigenEnergies[residualVectorIndex] * ritzVector[j];
      for(int k=0; k<this->matrixCISdimension; k++){
         double value = j<=k ? this->matrixCIS[j][k] : this->matrixCIS[k][j];
         residualVector[j] -= value*ritzVector[k];
      }
      sqNorm += pow(residualVector[j],2.0);
   }
   *norm = sqrt(sqNorm);
}

// This method is used for Davidson
void ZindoS::UpdateExpansionVectors(double** expansionVectors, 
                                    int* notConvergedStates, 
                                    double const* interactionEigenEnergies, 
                                    double const* residualVector,
                                    int interactionMatrixDimension, 
                                    int residualVectorIndex) const{
   double newExpansionVector[this->matrixCISdimension];
   // calculate new expansion vector from residual vector
   for(int j=0; j<this->matrixCISdimension; j++){
      double temp = interactionEigenEnergies[residualVectorIndex]-this->matrixCIS[j][j];
      if(temp == 0.0){
         // prevent dividing by 0.
         temp = pow(10,-100);
      }
      newExpansionVector[j]=pow(temp, -1.0)*residualVector[j];
   }

   // orthonormalize old expansion vectors and new expansion vector
   for(int k=0; k<interactionMatrixDimension+*notConvergedStates; k++){
      double overlap=0.0;
      for(int j=0; j<this->matrixCISdimension; j++){
         overlap += expansionVectors[j][k] * newExpansionVector[j];
      }
      for(int j=0; j<this->matrixCISdimension; j++){
         newExpansionVector[j] -= overlap*expansionVectors[j][k];
      }
   }

   // add new expansion vector to old expansion vectors
   double sqNormNewExpVect = 0.0;
   for(int j=0; j<this->matrixCISdimension; j++){
      sqNormNewExpVect += pow(newExpansionVector[j],2.0);
   }
   double normNewExpVect = sqrt(sqNormNewExpVect);
   for(int j=0; j<this->matrixCISdimension; j++){
      expansionVectors[j][interactionMatrixDimension+*notConvergedStates] 
            = pow(normNewExpVect,-1.0)*newExpansionVector[j];
   }
   *notConvergedStates += 1;
}

// This method is used for Davidson
void ZindoS::CalcInteractionMatrix(double** interactionMatrix, 
                                   double const* const* expansionVectors, 
                                   int interactionMatrixDimension) const{
   stringstream ompErrors;
   #pragma omp parallel for schedule(auto)
   for(int k=0; k<interactionMatrixDimension*interactionMatrixDimension; k++){
      try{
         int i = k/interactionMatrixDimension;
         int j = k%interactionMatrixDimension;
         if(i<=j){
            for(int a=0; a<this->matrixCISdimension; a++){
               for(int b=0; b<this->matrixCISdimension; b++){
                  double value = a<=b ? this->matrixCIS[a][b] : this->matrixCIS[b][a];
                  interactionMatrix[i][j] += expansionVectors[a][i] 
                                            *value
                                            *expansionVectors[b][j];
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

void ZindoS::DoCISDavidson(){
   this->OutputLog(this->messageStartDavidsonCIS);
   int numberOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberVir = Parameters::GetInstance()->GetActiveVirCIS();
   int numberExcitedStates = Parameters::GetInstance()->GetNumberExcitedStatesCIS();
   int maxIter = Parameters::GetInstance()->GetMaxIterationsCIS();
   int maxDim  = Parameters::GetInstance()->GetMaxDimensionsCIS();
   double normTol = Parameters::GetInstance()->GetNormToleranceCIS();
   bool convergeExcitedStates[numberExcitedStates];
   int interactionMatrixDimension;
   bool reachMaxDim;
   bool allConverged;
   int notConvergedStates;
   bool goToDirectCIS;
   double** expansionVectors = NULL;
   double*  ritzVector = NULL;
   double*  residualVector = NULL;

   try{
      MallocerFreer::GetInstance()->Malloc<double>(&expansionVectors, 
                                                   this->matrixCISdimension, 
                                                   maxDim);
      MallocerFreer::GetInstance()->Malloc<double>(&ritzVector, this->matrixCISdimension);
      MallocerFreer::GetInstance()->Malloc<double>(&residualVector, this->matrixCISdimension);
      // sort single excitation slater determinants
      vector<MoEnergyGap> moEnergyGaps;
      this->SortSingleExcitationSlaterDeterminants(&moEnergyGaps);

      // set initial expansion vectors and initial conveged vectors
      for(int k=0; k<numberExcitedStates; k++){
         expansionVectors[moEnergyGaps[k].slaterIndex][k] = 1.0;
         convergeExcitedStates[k] = false;
      }

      interactionMatrixDimension = 0;
      reachMaxDim = false;
      goToDirectCIS = false;
      // Davidson roop
      for(int k=0; k<maxIter; k++){
         this->OutputLog((boost::format("%s%d\n") % this->messageNumIterCIS.c_str() % k ).str());
         // calculate dimension of the interaction matrix 
         // (= number of the expansion vectors).
         for(int i=0; i<numberExcitedStates; i++){
            if(!convergeExcitedStates[i]){
               interactionMatrixDimension += 1;
            }
         }

         // malloc interaction matrix and etc.
         double** interactionMatrix = NULL;
         double* interactionEigenEnergies = NULL;

         try{
            MallocerFreer::GetInstance()->Malloc<double>(&interactionMatrix,
                                                         interactionMatrixDimension, 
                                                         interactionMatrixDimension);
            MallocerFreer::GetInstance()->Malloc<double>(&interactionEigenEnergies,
                                                         interactionMatrixDimension);
            // calculate interaction matrix
            this->CalcInteractionMatrix(interactionMatrix, 
                                        expansionVectors, 
                                        interactionMatrixDimension);

            // diagonalize interaction matrix
            bool calcEigenVectors = true;
            MolDS_wrappers::Lapack::GetInstance()->Dsyevd(interactionMatrix,
                                                          interactionEigenEnergies, 
                                                          interactionMatrixDimension, 
                                                          calcEigenVectors);

            // check convergence of all excited states
            notConvergedStates=0;
            allConverged = true;
            for(int i=0; i<numberExcitedStates; i++){
       
               // calculate i-th ritz vector
               this->CalcRitzVector(ritzVector, 
                                    expansionVectors, 
                                    interactionMatrix, 
                                    interactionMatrixDimension, 
                                    i);

               // calculate i-th residual vector and the norm of the residual vector
               double norm = 0.0;
               this->CalcResidualVectorAndNorm(residualVector, 
                                               &norm, 
                                               ritzVector, 
                                               interactionEigenEnergies, i);

               // output norm of residual vector
               this->OutputLog((boost::format("\t  %d%s%e\n") % (i+1) % this->messageResidualNorm.c_str() % norm ).str());
               if(i == numberExcitedStates-1){
                  this->OutputLog("\n");
               }

               // check tolerance for the norm of the residual vector.
               if(norm < normTol){
                  convergeExcitedStates[i] = true;
               }
               else{
                  convergeExcitedStates[i] = false;
                  allConverged = false;
                  if(interactionMatrixDimension+notConvergedStates == maxDim && maxDim !=this->matrixCISdimension){
                     reachMaxDim = true;
                     break;
                  }
                  else if(interactionMatrixDimension+notConvergedStates == this->matrixCISdimension){
                     goToDirectCIS = true;
                     break;
                  }
            
                  // update expansion vectors
                  this->UpdateExpansionVectors(expansionVectors, 
                                               &notConvergedStates, 
                                               interactionEigenEnergies, 
                                               residualVector,
                                               interactionMatrixDimension, 
                                               i);

               }
            } 

            if(allConverged){
               // copy to cis eigen vector and value
               for(int i=0; i<numberExcitedStates; i++){
                  this->excitedEnergies[i] = interactionEigenEnergies[i];
                  this->CalcRitzVector(ritzVector, 
                                       expansionVectors, 
                                       interactionMatrix, 
                                       interactionMatrixDimension, 
                                       i);
                  for(int j=0; j<this->matrixCISdimension; j++){
                     this->matrixCIS[i][j] = ritzVector[j];
                  }
               }
            }
         }
         catch(MolDSException ex){
            this->FreeDavidsonRoopCISTemporaryMtrices(&interactionMatrix, 
                                                      interactionMatrixDimension, 
                                                      &interactionEigenEnergies);
            throw ex;
         }
         this->FreeDavidsonRoopCISTemporaryMtrices(&interactionMatrix, 
                                                   interactionMatrixDimension, 
                                                   &interactionEigenEnergies);

         // stop the Davidson roop
         if(allConverged){
            this->OutputLog(this->messageDavidsonConverge);
            break;
         }
         else if(!allConverged && goToDirectCIS){
            this->OutputLog(this->messageDavidsonReachCISMatrix);
            this->OutputLog(this->messageDavidsonGoToDirect);
            break;
         }
         else if(!allConverged && reachMaxDim){
            stringstream ss;
            ss << endl;
            ss << this->errorMessageDavidsonNotConverged;
            ss << this->errorMessageDavidsonMaxDim << maxDim << endl;
            throw MolDSException(ss.str());
         }
         else if(!allConverged && k==maxIter-1){
            stringstream ss;
            ss << this->errorMessageDavidsonNotConverged;
            ss << this->errorMessageDavidsonMaxIter << maxIter << endl;
            throw MolDSException(ss.str());
         }

      }// end Davidson roop
   }
   catch(MolDSException ex){
      this->FreeDavidsonCISTemporaryMtrices(&expansionVectors, 
                                            &residualVector, 
                                            &ritzVector);
      throw ex;
   }
   this->FreeDavidsonCISTemporaryMtrices(&expansionVectors, 
                                         &residualVector, 
                                         &ritzVector);

   this->OutputLog(this->messageDoneDavidsonCIS);
   // change algorithm from Davidso to direct
   if(goToDirectCIS){
      this->DoCISDirect();
   }
}

void ZindoS::FreeDavidsonCISTemporaryMtrices(double*** expansionVectors, 
                                             double** residualVector, 
                                             double** ritzVector) const{
   int maxDim  = Parameters::GetInstance()->GetMaxDimensionsCIS();
   MallocerFreer::GetInstance()->Free<double>(expansionVectors, 
                                              this->matrixCISdimension,
                                              maxDim);
   MallocerFreer::GetInstance()->Free<double>(residualVector, this->matrixCISdimension);
   MallocerFreer::GetInstance()->Free<double>(ritzVector, this->matrixCISdimension);
}

void ZindoS::FreeDavidsonRoopCISTemporaryMtrices(double*** interactionMatrix, 
                                                 int interactionMatrixDimension, 
                                                 double** interactionEigenEnergies) const{
   MallocerFreer::GetInstance()->Free<double>(interactionMatrix, 
                                              interactionMatrixDimension,
                                              interactionMatrixDimension);
   MallocerFreer::GetInstance()->Free<double>(interactionEigenEnergies, interactionMatrixDimension);
}

void ZindoS::DoCISDirect(){
   this->OutputLog(this->messageStartDirectCIS);
   bool calcEigenVectors = true;
   MolDS_wrappers::Lapack::GetInstance()->Dsyevd(this->matrixCIS,
                                                 this->excitedEnergies, 
                                                 this->matrixCISdimension, 
                                                 calcEigenVectors);
   this->OutputLog(this->messageDoneDirectCIS);
}

void ZindoS::CalcCISMatrix(double** matrixCIS) const{
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
            int firstAOIndexA;
            int firstAOIndexB;
            int numberAOsA;
            int numberAOsB;
            double gamma;
            double exchange;
            double coulomb;
            OrbitalType orbitalMu;
            OrbitalType orbitalNu;
            // Off diagonal term (right upper)
            if(k<l){
               for(int A=0; A<molecule->GetNumberAtoms(); A++){
                  const Atom& atomA = *molecule->GetAtom(A);
                  firstAOIndexA = atomA.GetFirstAOIndex();
                  numberAOsA = atomA.GetValenceSize();

                  for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
                     orbitalMu = atomA.GetValence(mu-firstAOIndexA);

                     // CNDO term
                     for(int B=A; B<molecule->GetNumberAtoms(); B++){
                        const Atom& atomB = *molecule->GetAtom(B);
                        firstAOIndexB = atomB.GetFirstAOIndex();
                        numberAOsB = atomB.GetValenceSize();

                        for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
                           orbitalNu = atomB.GetValence(nu-firstAOIndexB);

                           if(A<B){
                              gamma = this->GetNishimotoMatagaTwoEleInt(atomA, 
                                                                        orbitalMu, 
                                                                        atomB, 
                                                                        orbitalNu);
                              value += 2.0*gamma*fockMatrix[moA][mu]
                                                *fockMatrix[moI][mu]
                                                *fockMatrix[moJ][nu]
                                                *fockMatrix[moB][nu];
                              value -=     gamma*fockMatrix[moA][mu]
                                                *fockMatrix[moB][mu]
                                                *fockMatrix[moI][nu]
                                                *fockMatrix[moJ][nu];
                              value += 2.0*gamma*fockMatrix[moA][nu]
                                                *fockMatrix[moI][nu]
                                                *fockMatrix[moJ][mu]
                                                *fockMatrix[moB][mu];
                              value -=     gamma*fockMatrix[moA][nu]
                                                *fockMatrix[moB][nu]
                                                *fockMatrix[moI][mu]
                                                *fockMatrix[moJ][mu];
                           }
                           else{
                              gamma = atomA.GetZindoF0ss();
                              value += 2.0*gamma*fockMatrix[moA][mu]
                                                *fockMatrix[moI][mu]
                                                *fockMatrix[moJ][nu]
                                                *fockMatrix[moB][nu];
                              value -=     gamma*fockMatrix[moA][mu]
                                                *fockMatrix[moB][mu]
                                                *fockMatrix[moI][nu]
                                                *fockMatrix[moJ][nu];
                           }  
                        }
                     }

                     // Aditional term for INDO or ZIND/S, see Eq. (10) in [RZ_1973]
                     for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
                        orbitalNu = atomA.GetValence(nu-firstAOIndexA);

                        if(mu!=nu){
                           exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
                           value += 2.0*exchange*fockMatrix[moA][mu]
                                                *fockMatrix[moI][nu]
                                                *fockMatrix[moJ][nu]
                                                *fockMatrix[moB][mu];
                           value += 2.0*exchange*fockMatrix[moA][mu]
                                                *fockMatrix[moI][nu]
                                                *fockMatrix[moJ][mu]
                                                *fockMatrix[moB][nu];
                           value -=     exchange*fockMatrix[moA][mu]
                                                *fockMatrix[moB][nu]
                                                *fockMatrix[moI][nu]
                                                *fockMatrix[moJ][mu];
                           value -=     exchange*fockMatrix[moA][mu]
                                                *fockMatrix[moB][nu]
                                                *fockMatrix[moI][mu]
                                                *fockMatrix[moJ][nu];
                        }

                        coulomb = this->GetCoulombInt(orbitalMu, orbitalNu, atomA);

                        if( (orbitalMu == s || orbitalMu == px || orbitalMu == py || orbitalMu == pz) &&
                            (orbitalNu == s || orbitalNu == px || orbitalNu == py || orbitalNu == pz) ){
                              gamma = atomA.GetZindoF0ss();
                        }
                        else{
                           stringstream ss;
                           ss << this->errorMessageCalcCISMatrix;
                           ss << this->errorMessageAtomType << AtomTypeStr(atomA.GetAtomType()) << "\n";
                           ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalMu) << "\n";
                           ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalNu) << "\n";
                           #pragma omp critical
                           ompErrors << ss.str() << endl ;
                        }   

                        value += 2.0*(coulomb-gamma)*fockMatrix[moA][mu]
                                                    *fockMatrix[moI][mu]
                                                    *fockMatrix[moJ][nu]
                                                    *fockMatrix[moB][nu];
                        value -=     (coulomb-gamma)*fockMatrix[moA][mu]
                                                    *fockMatrix[moB][mu]
                                                    *fockMatrix[moI][nu]
                                                    *fockMatrix[moJ][nu];
                     }
                  }
               }
            }
            // Diagonal term
            else if(k==l){
               value = this->energiesMO[moA] - this->energiesMO[moI];
               for(int A=0; A<molecule->GetNumberAtoms(); A++){
                  const Atom& atomA = *molecule->GetAtom(A);
                  firstAOIndexA = atomA.GetFirstAOIndex();
                  numberAOsA = atomA.GetValenceSize();

                  for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
                     orbitalMu = atomA.GetValence(mu-firstAOIndexA);

                     // CNDO term
                     for(int B=A; B<molecule->GetNumberAtoms(); B++){
                        const Atom& atomB = *molecule->GetAtom(B);
                        firstAOIndexB = atomB.GetFirstAOIndex();
                        numberAOsB = atomB.GetValenceSize();

                        for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
                           orbitalNu = atomB.GetValence(nu-firstAOIndexB);

                           if(A<B){
                              gamma = this->GetNishimotoMatagaTwoEleInt(atomA, 
                                                                        orbitalMu, 
                                                                        atomB, 
                                                                        orbitalNu);
                              value += 2.0*gamma*fockMatrix[moI][mu]
                                                *fockMatrix[moA][mu]
                                                *fockMatrix[moA][nu]
                                                *fockMatrix[moI][nu];
                              value -=     gamma*fockMatrix[moI][mu]
                                                *fockMatrix[moI][mu]
                                                *fockMatrix[moA][nu]
                                                *fockMatrix[moA][nu];
                              value += 2.0*gamma*fockMatrix[moI][nu]
                                                *fockMatrix[moA][nu]
                                                *fockMatrix[moA][mu]
                                                *fockMatrix[moI][mu];
                              value -=     gamma*fockMatrix[moI][nu]
                                                *fockMatrix[moI][nu]
                                                *fockMatrix[moA][mu]
                                                *fockMatrix[moA][mu];
                           }
                           else{
                              gamma = atomA.GetZindoF0ss();
                              value += 2.0*gamma*fockMatrix[moI][mu]
                                                *fockMatrix[moA][mu]
                                                *fockMatrix[moA][nu]
                                                *fockMatrix[moI][nu];
                              value -=     gamma*fockMatrix[moI][mu]
                                                *fockMatrix[moI][mu]
                                                *fockMatrix[moA][nu]
                                                *fockMatrix[moA][nu];
                           }     
                        }
                     }

                     // Aditional term for INDO or ZIND/S, see Eq. (10) in [RZ_1973]
                     for(int nu=firstAOIndexA; nu<firstAOIndexA+numberAOsA; nu++){
                        orbitalNu = atomA.GetValence(nu-firstAOIndexA);

                        if(mu!=nu){
                           exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
                           value += 2.0*exchange*fockMatrix[moI][mu]
                                                *fockMatrix[moA][nu]
                                                *fockMatrix[moA][nu]
                                                *fockMatrix[moI][mu];
                           value += 2.0*exchange*fockMatrix[moI][mu]
                                                *fockMatrix[moA][nu]
                                                *fockMatrix[moA][mu]
                                                *fockMatrix[moI][nu];
                           value -=     exchange*fockMatrix[moI][mu]
                                                *fockMatrix[moI][nu]
                                                *fockMatrix[moA][nu]
                                                *fockMatrix[moA][mu];
                           value -=     exchange*fockMatrix[moI][mu]
                                                *fockMatrix[moI][nu]
                                                *fockMatrix[moA][mu]
                                                *fockMatrix[moA][nu];
                        }

                        coulomb = this->GetCoulombInt(orbitalMu, orbitalNu, atomA);

                        if( (orbitalMu == s || orbitalMu == px || orbitalMu == py || orbitalMu == pz) &&
                            (orbitalNu == s || orbitalNu == px || orbitalNu == py || orbitalNu == pz) ){
                              gamma = atomA.GetZindoF0ss();
                        }
                        else{
                           stringstream ss;
                           ss << this->errorMessageCalcCISMatrix;
                           ss << this->errorMessageAtomType << AtomTypeStr(atomA.GetAtomType()) << "\n";
                           ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalMu) << "\n";
                           ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalNu) << "\n";
                           #pragma omp critical
                           ompErrors << ss.str() << endl ;
                        }   

                        value += 2.0*(coulomb-gamma)*fockMatrix[moI][mu]
                                                    *fockMatrix[moA][mu]
                                                    *fockMatrix[moA][nu]
                                                    *fockMatrix[moI][nu];
                        value -=     (coulomb-gamma)*fockMatrix[moI][mu]
                                                    *fockMatrix[moI][mu]
                                                    *fockMatrix[moA][nu]
                                                    *fockMatrix[moA][nu];
                     }
                  }
               }
            }
            // End of the fast algorith.
            /* 
            // Slow algorith, but this is easy to read. Fast altorithm is also written above.
            value = 2.0*this->GetMolecularIntegralElement(moA, moI, moJ, moB, 
                                                          this->molecule, 
                                                          this->fockMatrix, 
                                                          NULL)
                       -this->GetMolecularIntegralElement(moA, moB, moI, moJ, 
                                                          this->molecule, 
                                                          this->fockMatrix, 
                                                          NULL);
            if(k==l){
               value += this->energiesMO[moA] - this->energiesMO[moI];
            }
            // End of the slow algorith.
            */
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
                                                 % this->messageDoneCalcCISMatrix.c_str() ).str());
}

void ZindoS::CheckMatrixForce(const vector<int>& elecStates){
   // malloc or initialize Force matrix
   if(this->matrixForce == NULL){
      MallocerFreer::GetInstance()->Malloc<double>(&this->matrixForce, 
                                                   elecStates.size(),
                                                   this->molecule->GetNumberAtoms(), 
                                                   CartesianType_end);
      this->matrixForceElecStatesNum = elecStates.size();
   }
   else{
      MallocerFreer::GetInstance()->
      Initialize<double>(this->matrixForce,
                         elecStates.size(),
                         this->molecule->GetNumberAtoms(),
                         CartesianType_end);
   }
}

// Note taht activeOccIndex and activeVirIndex are not MO's number.
// activeOccIndex=0 means HOMO and activeVirIndex=0 means LUMO.
int ZindoS::GetSlaterDeterminantIndex(int activeOccIndex, 
                                      int activeVirIndex) const{
   return Parameters::GetInstance()->GetActiveVirCIS()
         *activeOccIndex
         +activeVirIndex; 
}

// This returns an index of occupied MO. Generally, This index=0 means the lowest energy MO;
int ZindoS::GetActiveOccIndex(const MolDS_base::Molecule& molecule, int matrixCISIndex) const{
   return molecule.GetTotalNumberValenceElectrons()/2 
         -(matrixCISIndex/Parameters::GetInstance()->GetActiveVirCIS()) -1;
}

// This returns an index of virtual MO. Generally, This index=0 means the lowest energy MO;
int ZindoS::GetActiveVirIndex(const MolDS_base::Molecule& molecule, int matrixCISIndex) const{
   return molecule.GetTotalNumberValenceElectrons()/2
         +(matrixCISIndex%Parameters::GetInstance()->GetActiveVirCIS());
}

// elecStates is indeces of the electroinc eigen states.
// The index = 0 means electronic ground state. 
void ZindoS::CalcForce(const vector<int>& elecStates){
   int elecState = elecStates[0];
   int groundState = 0;
   if(elecState != groundState){
      stringstream ss;
      ss << this->errorMessageCalcForceNotGroundState;
      ss << this->errorMessageElecState << elecState << "\n";
      throw MolDSException(ss.str());
   }

   this->CheckMatrixForce(elecStates);
   stringstream ompErrors;
   #pragma omp parallel 
   {
      double*** overlapDer=NULL;
      try{
         MallocerFreer::GetInstance()->Malloc<double>(&overlapDer,
                                                      OrbitalType_end, 
                                                      OrbitalType_end, 
                                                      CartesianType_end);

         #pragma omp for schedule(auto)
         for(int a=0; a<this->molecule->GetNumberAtoms(); a++){
            const Atom& atomA = *molecule->GetAtom(a);
            int firstAOIndexA = atomA.GetFirstAOIndex();
            int numberAOsA = atomA.GetValenceSize();
            double coreRepulsion[CartesianType_end] = {0.0,0.0,0.0};
            double electronicForce1[CartesianType_end] = {0.0,0.0,0.0};
            double electronicForce2[CartesianType_end] = {0.0,0.0,0.0};
            double electronicForce3[CartesianType_end] = {0.0,0.0,0.0};
            for(int b=0; b<this->molecule->GetNumberAtoms(); b++){
               if(a != b){
                  const Atom& atomB = *molecule->GetAtom(b);
                  int firstAOIndexB = atomB.GetFirstAOIndex();
                  int numberAOsB = atomB.GetValenceSize();

                  // calc. first derivative of overlap.
                  this->CalcDiatomicOverlapFirstDerivatives(overlapDer, atomA, atomB);

                  for(int i=0; i<CartesianType_end; i++){
                     coreRepulsion[i] += this->GetDiatomCoreRepulsionFirstDerivative(
                                               a, b, (CartesianType)i);
                     if(Parameters::GetInstance()->RequiresVdWSCF()){
                        coreRepulsion[i] += this->GetDiatomVdWCorrectionFirstDerivative(
                                                  a, b, (CartesianType)i);
                     }

                     electronicForce1[i] += ( atomA.GetCoreCharge()
                                             *atomicElectronPopulation[b]
                                             +atomB.GetCoreCharge()
                                             *atomicElectronPopulation[a])
                                             *this->GetNishimotoMatagaTwoEleIntFirstDerivative(
                                                    atomA, s, atomB, s, (CartesianType)i);
                  }
                  for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
                     OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                     for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
                        OrbitalType orbitalNu = atomB.GetValence(nu-firstAOIndexB);
                        double bondParameter = 0.5*(atomA.GetBondingParameter(
                                                           this->theory, orbitalMu) 
                                                   +atomB.GetBondingParameter(
                                                           this->theory, orbitalNu)); 
                        for(int i=0; i<CartesianType_end; i++){
                           electronicForce2[i] += 2.0*this->orbitalElectronPopulation[mu][nu]
                                                 *bondParameter
                                                 *overlapDer[mu-firstAOIndexA][nu-firstAOIndexB][i];
                           electronicForce3[i] += (this->orbitalElectronPopulation[mu][mu]
                                                  *this->orbitalElectronPopulation[nu][nu]
                                                  -0.5*pow(this->orbitalElectronPopulation[mu][nu],2.0))
                                                  *this->GetNishimotoMatagaTwoEleIntFirstDerivative(
                                                         atomA, orbitalMu, atomB, orbitalNu,
                                                         (CartesianType)i);
                        }
                     }
                  }
               }
            }
            for(int i=0; i<CartesianType_end; i++){
               this->matrixForce[elecState][a][i] = -1.0*(coreRepulsion[i]
                                                         -electronicForce1[i] 
                                                         +electronicForce2[i]
                                                         +electronicForce3[i]);
            }
         }
      }
      catch(MolDSException ex){
         #pragma omp critical
         ompErrors << ex.what() << endl ;
      }
      MallocerFreer::GetInstance()->Free<double>(&overlapDer, 
                                                 OrbitalType_end,
                                                 OrbitalType_end,
                                                 CartesianType_end);
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException(ompErrors.str());
   }
  
   /*
   // Calculate force. First derivative of overlap integral is
   // calculated with GTO expansion technique.
   stringstream ompErrors;
   #pragma omp parallel for schedule(auto)
   for(int a=0; a<this->molecule->GetNumberAtoms(); a++){
      try{
         const Atom& atomA = *molecule->GetAtom(a);
         int firstAOIndexA = atomA.GetFirstAOIndex();
         int numberAOsA = atomA.GetValenceSize();
         for(int i=0; i<CartesianType_end; i++){

            double coreRepulsion = 0.0;
            double electronicForce1 = 0.0;
            double electronicForce2 = 0.0;
            double electronicForce3 = 0.0;
            for(int b=0; b<this->molecule->GetNumberAtoms(); b++){
               if(a != b){
                  const Atom& atomB = *molecule->GetAtom(b);
                  int firstAOIndexB = atomB.GetFirstAOIndex();
                  int numberAOsB = atomB.GetValenceSize();

                  // Calculation of core repusion force
                  coreRepulsion += this->GetDiatomCoreRepulsionFirstDerivative(
                                         a, b, (CartesianType)i);
                  if(Parameters::GetInstance()->RequiresVdWSCF()){
                     coreRepulsion += this->GetDiatomVdWCorrectionFirstDerivative(
                                            a, b, (CartesianType)i);
                  }

                  // Calculate force arise from electronic part.
                  electronicForce1 += ( atomA.GetCoreCharge()*atomicElectronPopulation[b]
                                       +atomB.GetCoreCharge()*atomicElectronPopulation[a])
                                       *this->GetNishimotoMatagaTwoEleIntFirstDerivative
                                             (atomA, s, atomB, s, (CartesianType)i);

                  for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
                     OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                     for(int nu=firstAOIndexB; nu<firstAOIndexB+numberAOsB; nu++){
                        OrbitalType orbitalNu = atomB.GetValence(nu-firstAOIndexB);

                        double bondParameter = 0.5*(atomA.GetBondingParameter(this->theory, 
                                                                              orbitalMu) 
                                                   +atomB.GetBondingParameter(this->theory, 
                                                                              orbitalNu)); 

                        electronicForce2 += 2.0*this->orbitalElectronPopulation[mu][nu]
                                            *bondParameter
                                            *this->GetOverlapElementFirstDerivativeByGTOExpansion
                                                   (atomA, mu-firstAOIndexA, 
                                                    atomB, nu-firstAOIndexB,
                                                    STO6G, (CartesianType)i);

                        electronicForce3 += (this->orbitalElectronPopulation[mu][mu]
                                            *this->orbitalElectronPopulation[nu][nu]
                                            -0.5*pow(this->orbitalElectronPopulation[mu][nu],2.0))
                                            *this->GetNishimotoMatagaTwoEleIntFirstDerivative
                                                   (atomA, orbitalMu, atomB, orbitalNu,
                                                   (CartesianType)i);
                     }
                  }
               }
            }

            this->matrixForce[0][a][i] = -1.0*(coreRepulsion 
                                              -electronicForce1 
                                              +electronicForce2
                                              +electronicForce3);
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
   */
}

}



