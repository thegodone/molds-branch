//************************************************************************//
// Copyright (C) 2011-2013 Mikiya Fujii                                   //
// Copyright (C) 2012-2013 Michihiro Okuyama
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
#include"../base/Uncopyable.h"
#include"../mpi/MpiProcess.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../wrappers/Blas.h"
#include"../wrappers/Lapack.h"
#include"../base/Enums.h"
#include"../base/MathUtilities.h"
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
 *  Main Reference for Zindo is [RZ_1973]
 */

ZindoS::ZindoS() : MolDS_cndo::Cndo2(){
   //protected variables and methods
   this->theory = ZINDOS;
   this->SetMessages();
   this->SetEnableAtomTypes();

   this->zMatrixForceElecStatesNum = 0;
   this->etaMatrixForceElecStatesNum = 0;
   this->zMatrixForce = NULL;
   this->etaMatrixForce = NULL;

   //private variables
   this->nishimotoMatagaMatrix = NULL;
   this->matrixForceElecStatesNum = 0;
   this->nishimotoMatagaParamA = 1.2;
   this->nishimotoMatagaParamB = 2.4;
   this->overlapAOsCorrectionSigma = 1.267;
   this->overlapAOsCorrectionPi = 0.585;
   //this->OutputLog("ZindoS created\n");
}

ZindoS::~ZindoS(){
   if(this->theory==ZINDOS){
      MallocerFreer::GetInstance()->Free<double>(&this->nishimotoMatagaMatrix, 
                                                 this->molecule->GetNumberAtoms(), 
                                                 OrbitalType_end, 
                                                 this->molecule->GetNumberAtoms(), 
                                                 OrbitalType_end);
   }
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
   MallocerFreer::GetInstance()->Free<double>(&this->zMatrixForce, 
                                              this->zMatrixForceElecStatesNum,
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs());
   MallocerFreer::GetInstance()->Free<double>(&this->etaMatrixForce, 
                                              this->etaMatrixForceElecStatesNum,
                                              this->molecule->GetTotalNumberAOs(),
                                              this->molecule->GetTotalNumberAOs());
   if(Parameters::GetInstance()->RequiresMullikenCIS()){
      vector<int>* elecStates = Parameters::GetInstance()->GetElectronicStateIndecesMullikenCIS();
      MallocerFreer::GetInstance()->Free<double>(&this->orbitalElectronPopulationCIS, 
                                                 elecStates->size(),
                                                 this->molecule->GetTotalNumberAOs(),
                                                 this->molecule->GetTotalNumberAOs());
      MallocerFreer::GetInstance()->Free<double>(&this->atomicElectronPopulationCIS, 
                                                 elecStates->size(),
                                                 this->molecule->GetNumberAtoms());
      if(Parameters::GetInstance()->RequiresUnpairedPopCIS()){
         MallocerFreer::GetInstance()->Free<double>(&this->atomicUnpairedPopulationCIS, 
                                                    elecStates->size(),
                                                    this->molecule->GetNumberAtoms());
      }
   }
   //this->OutputLog("ZindoS deleted\n");
}

void ZindoS::SetMolecule(Molecule* molecule){
   Cndo2::SetMolecule(molecule);
   if(this->theory==ZINDOS){
      MallocerFreer::GetInstance()->Malloc<double>(&this->nishimotoMatagaMatrix, 
                                                   this->molecule->GetNumberAtoms(), 
                                                   OrbitalType_end, 
                                                   this->molecule->GetNumberAtoms(), 
                                                   OrbitalType_end);
   }
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
   this->errorMessageGetDiatomCoreRepulsion2ndDerivativeNotImplemented
      = "Error in indo::ZindoS::GetDiatomCoreRepulsion2ndDerivative: Second derivative is not implemented for ZINDO/S.\n";
   this->errorMessageCalcCISMatrix
      = "Error in zindo::ZindoS::CalcCISMatrix: Non available orbital is contained.\n";
   this->errorMessageDavidsonNotConverged =  "Error in zindo::ZindoS::DoCISDavidson: Davidson did not met convergence criterion. \n";
   this->errorMessageDavidsonMaxIter = "Davidson loop reaches max_iter=";
   this->errorMessageDavidsonMaxDim = "Dimension of the expansion vectors reaches max_dim=";
   this->errorMessageElecState = "Electronic State = ";
   this->errorMessageGetElectronicEnergyEnergyNotCalculated
      = "Error in zindo::ZindoS::GetElectronicEnergy: Set electronic state is not calculated by CIS.\n";
   this->errorMessageGetElectronicEnergyNULLCISEnergy 
      = "Error in zindo::ZindoS::GetElectronicEnergy: excitedEnergies is NULL\n";
   this->errorMessageCalcElectronicTransitionDipoleMomentBadState
      = "Error in zindo::ZindoS::CalcElectronicTransitionDipoleMoment: Bad eigen state is set to calculate the transition dipole moment. Note taht state=0 means the ground state and other state = i means the i-th excited state in below.\n";
   this->errorMessageCalcFrequenciesNormalModesBadTheory
      = "Error in zindo::ZindoS::CalcFrequenciesNormalModesBadTheory: ZINDO/S is not supported for frequency (normal mode) analysis.\n";
   this->errorMessageCalcZMatrixForceEtaNull 
      = "Error in zindo::ZindoS::CalcZMatrixForce: Nndo::etaMatrixForce is NULL. Call Mndo::CalcEtaMatrixForce before calling Mndo::CalcZMatrixForce.\n";
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
         if(B != indexAtomA){
            const Atom& atomB = *molecule.GetAtom(B);
            OrbitalType orbitalSigma;
            int sigma;
            int atomBNumberValence = atomB.GetValenceSize();
            double rAB = molecule.GetDistanceAtoms(atomA, atomB);
            for(int i=0; i<atomBNumberValence; i++){
               sigma = i + atomB.GetFirstAOIndex();
               orbitalSigma = atomB.GetValence(i);
               temp += orbitalElectronPopulationDiagPart[sigma]
                      *this->nishimotoMatagaMatrix[indexAtomA][orbitalMu][B][orbitalSigma];
            }
            temp -= atomB.GetCoreCharge() 
                   *this->nishimotoMatagaMatrix[indexAtomA][s][B][s];
         }
      }
      value += temp;
   }

   return value;
}

double ZindoS::GetFockOffDiagElement(const Atom& atomA, 
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
   OrbitalType orbitalMu = atomA.GetValence(mu-atomA.GetFirstAOIndex());
   OrbitalType orbitalNu = atomB.GetValence(nu-atomB.GetFirstAOIndex());
   double bondParameter = 0.5*(atomA.GetBondingParameter(this->theory, orbitalMu) 
                              +atomB.GetBondingParameter(this->theory, orbitalNu)); 

   if(isGuess){
      value = bondParameter*overlapAOs[mu][nu];
   }
   else{
      double coulomb = 0.0;
      double exchange = 0.0;
      if(indexAtomA == indexAtomB){
         coulomb  = this->GetCoulombInt(orbitalMu, orbitalNu, atomA); 
         exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA); 
         value = (1.5*exchange - 0.5*coulomb)*orbitalElectronPopulation[mu][nu];
      }
      else{
         value = bondParameter*overlapAOs[mu][nu];
         value -= 0.5*orbitalElectronPopulation[mu][nu]
                  *this->nishimotoMatagaMatrix[indexAtomA][orbitalMu][indexAtomB][orbitalNu];
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

void ZindoS::CalcTwoElecTwoCore(double****** twoElecTwoCore, 
                              const Molecule& molecule) const{
   this->CalcNishimotoMatagaMatrix(this->nishimotoMatagaMatrix, molecule);
}

// ref. [MN_1957] and (5a) in [AEZ_1986]
double ZindoS::GetNishimotoMatagaTwoEleInt(const Atom& atomA, OrbitalType orbitalA, 
                                           const Atom& atomB, OrbitalType orbitalB) const{
   double r = this->molecule->GetDistanceAtoms(atomA, atomB);
   return this->GetNishimotoMatagaTwoEleInt(atomA, orbitalA, atomB, orbitalB,r);
}

// ref. [MN_1957] and (5a) in [AEZ_1986]
double ZindoS::GetNishimotoMatagaTwoEleInt(const Atom& atomA, OrbitalType orbitalA, 
                                           const Atom& atomB, OrbitalType orbitalB,
                                           const double rAB) const{
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

   double gamma=gammaAA+gammaBB;
   return this->nishimotoMatagaParamA/( rAB+this->nishimotoMatagaParamB/gamma );

}

// First derivative of Nishimoto-Mataga related to the coordinate of atom A.
// For Nishimoto-Mataga, See ZindoS::GetNishimotoMatagaTwoEleInt 
// or ref. [MN_1957] and (5a) in [AEZ_1986]
double ZindoS::GetNishimotoMatagaTwoEleInt1stDerivative(const Atom& atomA, 
                                                        OrbitalType orbitalA, 
                                                        const Atom& atomB, 
                                                        OrbitalType orbitalB,
                                                        CartesianType axisA) const{
   double r = this->molecule->GetDistanceAtoms(atomA, atomB);
   return this->GetNishimotoMatagaTwoEleInt1stDerivative(atomA, orbitalA, atomB, orbitalB, r, axisA);
}

// First derivative of Nishimoto-Mataga related to the coordinate of atom A.
// For Nishimoto-Mataga, See ZindoS::GetNishimotoMatagaTwoEleInt 
// or ref. [MN_1957] and (5a) in [AEZ_1986]
double ZindoS::GetNishimotoMatagaTwoEleInt1stDerivative(const Atom& atomA, 
                                                        OrbitalType orbitalA, 
                                                        const Atom& atomB, 
                                                        OrbitalType orbitalB,
                                                        const double rAB,
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

   double dCartesian = atomA.GetXyz()[axisA] - atomB.GetXyz()[axisA];
   double value = -1.0*dCartesian/rAB;
   value *= this->nishimotoMatagaParamA;
   value *= pow( rAB+this->nishimotoMatagaParamB/(gammaAA+gammaBB) ,-2.0);
   return value;
}

void ZindoS::CalcNishimotoMatagaMatrix(double**** nishimotoMatagaMatrix, const Molecule& molecule) const{
   int totalNumberAtoms = molecule.GetNumberAtoms();
   stringstream ompErrors;
#pragma omp parallel for schedule(auto)
   for(int A=0; A<totalNumberAtoms; A++){
      try{
         const Atom& atomA = *molecule.GetAtom(A);
         int firstAOIndexA = atomA.GetFirstAOIndex();
         int lastAOIndexA  = atomA.GetLastAOIndex();
         for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
            OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
            for(int B=A; B<totalNumberAtoms; B++){
               const Atom& atomB = *molecule.GetAtom(B);
               int firstAOIndexB = atomB.GetFirstAOIndex();
               int lastAOIndexB  = atomB.GetLastAOIndex();
               double rAB = molecule.GetDistanceAtoms(atomA, atomB);
               for(int nu=firstAOIndexB; nu<=lastAOIndexB; nu++){
                  OrbitalType orbitalNu = atomB.GetValence(nu-firstAOIndexB);
                  nishimotoMatagaMatrix[A][orbitalMu][B][orbitalNu] = this->GetNishimotoMatagaTwoEleInt(atomA, 
                                                                                                        orbitalMu, 
                                                                                                        atomB, 
                                                                                                        orbitalNu,
                                                                                                        rAB);
                  if(A!=B){
                     nishimotoMatagaMatrix[B][orbitalNu][A][orbitalMu] = nishimotoMatagaMatrix[A][orbitalMu][B][orbitalNu];
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
}

void ZindoS::CalcDiatomicOverlapAOsInDiatomicFrame(double** diatomicOverlapAOs, 
                                                   const Atom& atomA, 
                                                   const Atom& atomB) const{
   MolDS_cndo::Cndo2::CalcDiatomicOverlapAOsInDiatomicFrame(diatomicOverlapAOs, atomA, atomB);

   // see (4f) in [AEZ_1986]
   diatomicOverlapAOs[pz][pz] *= this->overlapAOsCorrectionSigma;
   diatomicOverlapAOs[py][py] *= this->overlapAOsCorrectionPi;
   diatomicOverlapAOs[px][px] *= this->overlapAOsCorrectionPi;
   
   for(int i=0;i<OrbitalType_end;i++){
      for(int j=0;j<OrbitalType_end;j++){
         //this->OutputLog(boost::format("diatomicOverlapAOs[%d][%d]=%lf\n") % i % j % diatomicOverlapAOs[i][j]);
      }
   }
   
}

// First derivative of (B.40) in J. A. Pople book with bond correction.
void ZindoS::CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(double** diatomicOverlapAOsDeri, 
                                                                const Atom& atomA, 
                                                                const Atom& atomB) const{

   MolDS_cndo::Cndo2::CalcDiatomicOverlapAOs1stDerivativeInDiatomicFrame(diatomicOverlapAOsDeri,atomA, atomB);

   // see (4f) in [AEZ_1986] like as overlapAOs integlral
   diatomicOverlapAOsDeri[pz][pz] *= this->overlapAOsCorrectionSigma;
   diatomicOverlapAOsDeri[py][py] *= this->overlapAOsCorrectionPi;
   diatomicOverlapAOsDeri[px][px] *= this->overlapAOsCorrectionPi;
}

// Second derivative of (B.40) in J. A. Pople book with bond correction.
void ZindoS::CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(double** diatomicOverlapAOs2ndDeri, 
                                                                const Atom& atomA, 
                                                                const Atom& atomB) const{

   MolDS_cndo::Cndo2::CalcDiatomicOverlapAOs2ndDerivativeInDiatomicFrame(
                      diatomicOverlapAOs2ndDeri,atomA, atomB);

   // see (4f) in [AEZ_1986] like as overlapAOs integlral
   diatomicOverlapAOs2ndDeri[pz][pz] *= this->overlapAOsCorrectionSigma;
   diatomicOverlapAOs2ndDeri[py][py] *= this->overlapAOsCorrectionPi;
   diatomicOverlapAOs2ndDeri[px][px] *= this->overlapAOsCorrectionPi;
}

// calculate OverlapSingletSDs matrix between different electronic-structure, S^{SSD}_{ij}.
// i and j are singlet SDs belonging to left and right hand side electronic-structures, respectively.
// The index i=0 means the Hartree-Fock state.
// This overlapSingletSDs are calculated from overlapMOs.
// Note that rhs-electronic-structure is this electronic-structure  
// and lhs-electronic-structure is another electronic-structure.
void ZindoS::CalcOverlapSingletSDsWithAnotherElectronicStructure(double** overlapSingletSDs, 
                                                                 double const* const* overlapMOs) const{
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   MallocerFreer::GetInstance()->Initialize<double>(overlapSingletSDs, 
                                                    this->matrixCISdimension, 
                                                    this->matrixCISdimension);
   double** tmpMatrix1=NULL;
   double** tmpMatrix2=NULL;
   double** tmpMatrix3=NULL;
   double sqrtGroundStateOverlap;
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&tmpMatrix1, numberOcc, numberOcc);
      MallocerFreer::GetInstance()->Malloc<double>(&tmpMatrix2, numberOcc, numberOcc);
      MallocerFreer::GetInstance()->Malloc<double>(&tmpMatrix3, numberOcc, numberOcc);
      // between ground state
      for(int i=0; i<numberOcc; i++){
         for(int j=0; j<numberOcc; j++){
            tmpMatrix1[i][j] = overlapMOs[i][j];
         }
      }
      sqrtGroundStateOverlap = GetDeterminant(tmpMatrix1, numberOcc);
      overlapSingletSDs[0][0] = pow(sqrtGroundStateOverlap,2.0);

      for(int k=0; k<this->matrixCISdimension; k++){
         // single excitation from I-th (occupied)MO to A-th (virtual)MO
         int moI = this->GetActiveOccIndex(*this->molecule, k);
         int moA = this->GetActiveVirIndex(*this->molecule, k);
         for(int l=0; l<this->matrixCISdimension; l++){
            // single excitation from I-th (occupied)MO to A-th (virtual)MO
            int moJ = this->GetActiveOccIndex(*this->molecule, l);
            int moB = this->GetActiveVirIndex(*this->molecule, l);
            for(int i=0; i<numberOcc; i++){
               int destMO = i==moI ? moA : i;
               for(int j=0; j<numberOcc; j++){
                  int sourceMO = j==moJ ? moB : j;
                  tmpMatrix1[i][j] = overlapMOs[destMO][j];
                  tmpMatrix2[i][j] = overlapMOs[i][sourceMO];
                  tmpMatrix3[i][j] = overlapMOs[destMO][sourceMO];
               }
            }
            double det1 = GetDeterminant(tmpMatrix1, numberOcc);
            double det2 = GetDeterminant(tmpMatrix2, numberOcc);
            double det3 = GetDeterminant(tmpMatrix3, numberOcc);
            // from ground state to singet SDs
            overlapSingletSDs[k+1][0]   = sqrtGroundStateOverlap*det1;
            // from singet SDs to ground state
            overlapSingletSDs[0]  [l+1] = sqrtGroundStateOverlap*det2;
            // from singet SDs to singlet SDs
            overlapSingletSDs[k+1][l+1] = sqrtGroundStateOverlap*det3 + det1*det2;
         }
      }
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&tmpMatrix1, numberOcc, numberOcc);
      MallocerFreer::GetInstance()->Free<double>(&tmpMatrix2, numberOcc, numberOcc);
      MallocerFreer::GetInstance()->Free<double>(&tmpMatrix3, numberOcc, numberOcc);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&tmpMatrix1, numberOcc, numberOcc);
   MallocerFreer::GetInstance()->Free<double>(&tmpMatrix2, numberOcc, numberOcc);
   MallocerFreer::GetInstance()->Free<double>(&tmpMatrix3, numberOcc, numberOcc);
}

// calculate overlapESs (ES means eigenstate) matrix between different electronic-structure, S^{ES}_{ij}.
// i and j are singlet SDs belonging to left and right hand side electronic-structures, respectively.
// The index i=0 means the ground state.
// This overlapESs is calculated from the overlapsingletSDs.
// Note that rhs-electronic-structure is this electronic-structure  
// and lhs-electronic-structure is another electronic-structure.
void ZindoS::CalcOverlapESsWithAnotherElectronicStructure(double** overlapESs, 
                                                          double const* const* overlapSingletSDs,
                                                          const MolDS_base::ElectronicStructure& lhsElectronicStructure) const{
   const ElectronicStructure* rhsElectronicStructure = this;
   double const* const* rhsMatrixCIS = this->matrixCIS;
   double const* const* lhsMatrixCIS = lhsElectronicStructure.GetMatrixCIS();
   int dimOverlapSingletSDs = this->matrixCISdimension + 1;
   int dimOverlapESs = Parameters::GetInstance()->GetNumberElectronicStatesNASCO();
   int groundstate = 0;
   MallocerFreer::GetInstance()->Initialize<double>(overlapESs, dimOverlapESs, dimOverlapESs);
   // extended CIS matrix includes groundstate althoug matrixCIS does not include groundstate.
   double** lhsExtendedMatrixCIS=NULL;
   double** rhsExtendedMatrixCIS=NULL;
   double** tmpMatrix=NULL;
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&lhsExtendedMatrixCIS, dimOverlapSingletSDs, dimOverlapSingletSDs);
      MallocerFreer::GetInstance()->Malloc<double>(&rhsExtendedMatrixCIS, dimOverlapSingletSDs, dimOverlapSingletSDs);
      MallocerFreer::GetInstance()->Malloc<double>(&tmpMatrix,            dimOverlapSingletSDs, dimOverlapESs);
      lhsExtendedMatrixCIS[groundstate][groundstate] = 1.0;
      rhsExtendedMatrixCIS[groundstate][groundstate] = 1.0;
      for(int i=1; i<dimOverlapESs; i++){
         for(int j=1; j<dimOverlapSingletSDs; j++){
            rhsExtendedMatrixCIS[i][j] = rhsMatrixCIS[i-1][j-1];
            lhsExtendedMatrixCIS[i][j] = lhsMatrixCIS[i-1][j-1];
         }
      }
      // calc. overlap between eigenstates
      bool isColumnMajorOverlapSingletSDs = false;
      bool isColumnMajorRhsMatrixCIS = true;
      double alpha=1.0;
      double beta=0.0;
      MolDS_wrappers::Blas::GetInstance()->Dgemm(isColumnMajorOverlapSingletSDs,
                                                 isColumnMajorRhsMatrixCIS,
                                                 dimOverlapSingletSDs, dimOverlapESs, dimOverlapSingletSDs,
                                                 alpha,
                                                 overlapSingletSDs,
                                                 rhsExtendedMatrixCIS,
                                                 beta,
                                                 tmpMatrix);
      MolDS_wrappers::Blas::GetInstance()->Dgemm(dimOverlapESs, dimOverlapESs, dimOverlapSingletSDs,
                                                 lhsExtendedMatrixCIS,
                                                 tmpMatrix,
                                                 overlapESs);
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&lhsExtendedMatrixCIS, dimOverlapSingletSDs, dimOverlapSingletSDs);
      MallocerFreer::GetInstance()->Free<double>(&rhsExtendedMatrixCIS, dimOverlapSingletSDs, dimOverlapSingletSDs);
      MallocerFreer::GetInstance()->Free<double>(&tmpMatrix,            dimOverlapSingletSDs, dimOverlapESs);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&lhsExtendedMatrixCIS, dimOverlapSingletSDs, dimOverlapSingletSDs);
   MallocerFreer::GetInstance()->Free<double>(&rhsExtendedMatrixCIS, dimOverlapSingletSDs, dimOverlapSingletSDs);
   MallocerFreer::GetInstance()->Free<double>(&tmpMatrix,            dimOverlapSingletSDs, dimOverlapESs);
}

// The order of moI, moJ, moK, moL is consistent with Eq. (9) in [RZ_1973]
double ZindoS::GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                           const Molecule& molecule, 
                                           double const* const* fockMatrix, 
                                           double const* const* gammaAB) const{
   double value = 0.0;
   double gamma;
   double exchange;
   double coulomb;
   for(int A=0; A<molecule.GetNumberAtoms(); A++){
      const Atom& atomA = *molecule.GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();

      for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
         OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);

         // CNDO term
         for(int B=A; B<molecule.GetNumberAtoms(); B++){
            const Atom& atomB = *molecule.GetAtom(B);
            int firstAOIndexB = atomB.GetFirstAOIndex();
            int lastAOIndexB  = atomB.GetLastAOIndex();

            for(int nu=firstAOIndexB; nu<=lastAOIndexB; nu++){
               OrbitalType orbitalNu = atomB.GetValence(nu-firstAOIndexB);

               if(A<B){
                  gamma = this->nishimotoMatagaMatrix[A][orbitalMu][B][orbitalNu];
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
         for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
            OrbitalType orbitalNu = atomA.GetValence(nu-firstAOIndexA);

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
   this->OutputLog(boost::format("%s%lf%s\n%s") % this->messageOmpElapsedTimeCIS.c_str()
                                                % (ompEndTime - ompStartTime)
                                                % this->messageUnitSec.c_str()
                                                % this->messageDoneCIS.c_str() );
}

void ZindoS::CalcCISProperties(){
//calculate dipole moments and transitiondipolemoment
{
   int totalNumberAOs = this->molecule->GetTotalNumberAOs();
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   int mpiRank = MolDS_mpi::MpiProcess::GetInstance()->GetRank();
   int mpiSize = MolDS_mpi::MpiProcess::GetInstance()->GetSize();
   double*** dipoleMOs = NULL;
   double**  overlapMOs = NULL;
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&dipoleMOs,       CartesianType_end, totalNumberAOs, totalNumberAOs);
      MallocerFreer::GetInstance()->Malloc<double>(&overlapMOs, totalNumberAOs, totalNumberAOs);
      double alpha=1.0;
      double beta =0.0;
      //double ompStartTime = omp_get_wtime();
      MolDS_wrappers::Blas::GetInstance()->Dgemmm(false, false, true, totalNumberAOs, totalNumberAOs, totalNumberAOs, totalNumberAOs,
                                                  alpha, 
                                                  this->fockMatrix,
                                                  this->cartesianMatrix[XAxis],
                                                  this->fockMatrix,
                                                  beta,
                                                  dipoleMOs[XAxis]);
      MolDS_wrappers::Blas::GetInstance()->Dgemmm(false, false, true, totalNumberAOs, totalNumberAOs, totalNumberAOs, totalNumberAOs,
                                                  alpha, 
                                                  this->fockMatrix,
                                                  this->cartesianMatrix[YAxis],
                                                  this->fockMatrix,
                                                  beta,
                                                  dipoleMOs[YAxis]);
      MolDS_wrappers::Blas::GetInstance()->Dgemmm(false, false, true, totalNumberAOs, totalNumberAOs, totalNumberAOs, totalNumberAOs,
                                                  alpha, 
                                                  this->fockMatrix,
                                                  this->cartesianMatrix[ZAxis],
                                                  this->fockMatrix,
                                                  beta,
                                                  dipoleMOs[ZAxis]);
 
      double const* centerOfDipole = this->molecule->GetXyzCOC();
      // set orign of dipole
      MolDS_wrappers::Blas::GetInstance()->Dgemmm(false, false, true, totalNumberAOs, totalNumberAOs, totalNumberAOs, totalNumberAOs,
                                                  alpha, 
                                                  this->fockMatrix,
                                                  this->overlapAOs,
                                                  this->fockMatrix,
                                                  beta,
                                                  overlapMOs);
      MolDS_wrappers::Blas::GetInstance()->Daxpy(totalNumberAOs*totalNumberAOs,
                                                 -centerOfDipole[XAxis], 
                                                 &overlapMOs[0][0],
                                                 &dipoleMOs[XAxis][0][0]);
      MolDS_wrappers::Blas::GetInstance()->Daxpy(totalNumberAOs*totalNumberAOs,
                                                 -centerOfDipole[YAxis], 
                                                 &overlapMOs[0][0],
                                                 &dipoleMOs[YAxis][0][0]);
      MolDS_wrappers::Blas::GetInstance()->Daxpy(totalNumberAOs*totalNumberAOs,
                                                 -centerOfDipole[ZAxis], 
                                                 &overlapMOs[0][0],
                                                 &dipoleMOs[ZAxis][0][0]);
 
 
      // dipole moments of excited states
      //this->CalcElectronicDipoleMomentsExcitedStates(this->electronicTransitionDipoleMoments,
      //                                               this->fockMatrix,
      //                                               this->matrixCIS,
      //                                               this->cartesianMatrix,
      //                                               *this->molecule, 
      //                                               this->orbitalElectronPopulation,
      //                                               this->overlapAOs);
      int groundState = 0;
      for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
         int excitedState = k+1; // (k+1)-th excited state
         if(excitedState%mpiSize != mpiRank){continue;}
         this->electronicTransitionDipoleMoments[excitedState][excitedState][XAxis] = this->electronicTransitionDipoleMoments[groundState][groundState][XAxis];
         this->electronicTransitionDipoleMoments[excitedState][excitedState][YAxis] = this->electronicTransitionDipoleMoments[groundState][groundState][YAxis];
         this->electronicTransitionDipoleMoments[excitedState][excitedState][ZAxis] = this->electronicTransitionDipoleMoments[groundState][groundState][ZAxis];
         double tmpX=0.0, tmpY=0.0, tmpZ=0.0;
#pragma omp parallel for reduction(+:tmpX,tmpY,tmpZ) schedule(auto)
         for(int l=0; l<this->matrixCISdimension; l++){
            // single excitation from I-th (occupied)MO to A-th (virtual)MO
            int moI = this->GetActiveOccIndex(*this->molecule, l);
            int moA = this->GetActiveVirIndex(*this->molecule, l);
            double temp    = matrixCIS[k][l]*matrixCIS[k][l];
            tmpX += temp*(-dipoleMOs[XAxis][moI][moI]+dipoleMOs[XAxis][moA][moA]);
            tmpY += temp*(-dipoleMOs[YAxis][moI][moI]+dipoleMOs[YAxis][moA][moA]);
            tmpZ += temp*(-dipoleMOs[ZAxis][moI][moI]+dipoleMOs[ZAxis][moA][moA]);
         }
         this->electronicTransitionDipoleMoments[excitedState][excitedState][XAxis] += tmpX;
         this->electronicTransitionDipoleMoments[excitedState][excitedState][YAxis] += tmpY;
         this->electronicTransitionDipoleMoments[excitedState][excitedState][ZAxis] += tmpZ;
         
      }
 
      // transition dipole moment
      //this->CalcElectronicTransitionDipoleMoments(this->electronicTransitionDipoleMoments,
      //                                            this->fockMatrix,
      //                                            this->matrixCIS,
      //                                            this->cartesianMatrix,
      //                                            *this->molecule, 
      //                                            this->orbitalElectronPopulation,
      //                                            this->overlapAOs);
      for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
         int excitedState = k+1; // (k+1)-th excited state
         if(excitedState%mpiSize != mpiRank){continue;}
         this->electronicTransitionDipoleMoments[excitedState][groundState][XAxis] = 0.0;
         this->electronicTransitionDipoleMoments[excitedState][groundState][YAxis] = 0.0;
         this->electronicTransitionDipoleMoments[excitedState][groundState][ZAxis] = 0.0;
         double tmpX=0.0, tmpY=0.0, tmpZ=0.0;
#pragma omp parallel for reduction(+:tmpX,tmpY,tmpZ) schedule(auto)
         for(int l=0; l<this->matrixCISdimension; l++){
            // single excitation from I-th (occupied)MO to A-th (virtual)MO
            int moI = this->GetActiveOccIndex(*this->molecule, l);
            int moA = this->GetActiveVirIndex(*this->molecule, l);
            //double temp    = matrixCIS[k][l]*matrixCIS[k][l];
            double tmp    = this->matrixCIS[k][l]*sqrt(2.0);
            tmpX += tmp*dipoleMOs[XAxis][moA][moI];
            tmpY += tmp*dipoleMOs[YAxis][moA][moI];
            tmpZ += tmp*dipoleMOs[ZAxis][moA][moI];
         }
         this->electronicTransitionDipoleMoments[excitedState][groundState][XAxis] += tmpX;
         this->electronicTransitionDipoleMoments[excitedState][groundState][YAxis] += tmpY;
         this->electronicTransitionDipoleMoments[excitedState][groundState][ZAxis] += tmpZ;
      }
      if(Parameters::GetInstance()->RequiresAllTransitionDipoleMomentsCIS()){
         for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
            int departureExcitedState = k+1; // (k+1)-th excited state
            for(int l=k+1; l<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); l++){
               int destinationExcitedState = l+1; // (l+1)-th excited state
               if(destinationExcitedState%mpiSize != mpiRank){continue;}
               this->electronicTransitionDipoleMoments[destinationExcitedState][departureExcitedState][XAxis] = 0.0;
               this->electronicTransitionDipoleMoments[destinationExcitedState][departureExcitedState][YAxis] = 0.0;
               this->electronicTransitionDipoleMoments[destinationExcitedState][departureExcitedState][ZAxis] = 0.0;
               double tmpX=0.0, tmpY=0.0, tmpZ=0.0;
#pragma omp parallel for reduction(+:tmpX,tmpY,tmpZ) schedule(auto)
               for(int l=0; l<this->matrixCISdimension; l++){
                  // single excitation from I-th (occupied)MO to A-th (virtual)MO
                  int moI = this->GetActiveOccIndex(*this->molecule, l);
                  int moA = this->GetActiveVirIndex(*this->molecule, l);
                  double tmp    = matrixCIS[departureExcitedState-1][l]*matrixCIS[destinationExcitedState-1][l];
                  tmpX += tmp*(-dipoleMOs[XAxis][moI][moI]+dipoleMOs[XAxis][moA][moA]);
                  tmpY += tmp*(-dipoleMOs[YAxis][moI][moI]+dipoleMOs[YAxis][moA][moA]);
                  tmpZ += tmp*(-dipoleMOs[ZAxis][moI][moI]+dipoleMOs[ZAxis][moA][moA]);
               }
               this->electronicTransitionDipoleMoments[destinationExcitedState][departureExcitedState][XAxis] += tmpX;
               this->electronicTransitionDipoleMoments[destinationExcitedState][departureExcitedState][YAxis] += tmpY;
               this->electronicTransitionDipoleMoments[destinationExcitedState][departureExcitedState][ZAxis] += tmpZ;
            }
         }
      }
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&dipoleMOs, CartesianType_end, totalNumberAOs, totalNumberAOs);
      MallocerFreer::GetInstance()->Free<double>(&overlapMOs, totalNumberAOs, totalNumberAOs);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&dipoleMOs, CartesianType_end, totalNumberAOs, totalNumberAOs);
   MallocerFreer::GetInstance()->Free<double>(&overlapMOs, totalNumberAOs, totalNumberAOs);

   // communication to collect all matrix data on rank 0
   int numTransported = (Parameters::GetInstance()->GetNumberExcitedStatesCIS()+1)*CartesianType_end;
   if(mpiRank == 0){
      // receive the matrix data from other ranks
      for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS()+1; k++){
         if(k%mpiSize == 0){continue;}
         int source = k%mpiSize;
         int tag = k;
         MolDS_mpi::MpiProcess::GetInstance()->Recv(source, tag, &this->electronicTransitionDipoleMoments[k][0][0], numTransported);
      }
   }
   else{
      // send the matrix data to rank-0
      for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS()+1; k++){
         if(k%mpiSize != mpiRank){continue;}
         int dest = 0;
         int tag = k;
         MolDS_mpi::MpiProcess::GetInstance()->Send(dest, tag, &this->electronicTransitionDipoleMoments[k][0][0], numTransported);
      }
   }

   // right upper part of the matrix is copied from left lower part.
   if(mpiRank == 0 && Parameters::GetInstance()->RequiresAllTransitionDipoleMomentsCIS()){
      for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS()+1; k++){
         for(int l=k+1; l<Parameters::GetInstance()->GetNumberExcitedStatesCIS()+1; l++){
            for(int axis=0; axis<CartesianType_end; axis++){
               this->electronicTransitionDipoleMoments[k][l][axis] 
                  = this->electronicTransitionDipoleMoments[l][k][axis];
            }
         }
      }
   }

   // broadcast all matrix data to all ranks
   numTransported *= (Parameters::GetInstance()->GetNumberExcitedStatesCIS()+1);
   int root=0;
   MolDS_mpi::MpiProcess::GetInstance()->Broadcast(&this->electronicTransitionDipoleMoments[0][0][0], numTransported, root);


}// end of "calculate dipole moments and transitiondipolemoment"



   // free exciton energies
   this->CalcFreeExcitonEnergies(&this->freeExcitonEnergiesCIS, 
                                 *this->molecule, 
                                 this->energiesMO, 
                                 this->matrixCIS,
                                 this->matrixCISdimension);
   
   // orbital electron population
   this->CalcOrbitalElectronPopulationCIS(&this->orbitalElectronPopulationCIS, 
                                          this->orbitalElectronPopulation,
                                          *this->molecule, 
                                          this->fockMatrix,
                                          this->matrixCIS);

   // atomic electron population
   this->CalcAtomicElectronPopulationCIS(&this->atomicElectronPopulationCIS,
                                         this->orbitalElectronPopulationCIS, 
                                         *this->molecule);
   // atomic unpaired electron population
   this->CalcAtomicUnpairedPopulationCIS(&this->atomicUnpairedPopulationCIS,
                                         this->orbitalElectronPopulationCIS, 
                                         *this->molecule);
}

void ZindoS::CalcElectronicDipoleMomentsExcitedStates(double*** electronicTransitionDipoleMoments,
                                                      double const* const* fockMatrix,
                                                      double const* const* matrixCIS,
                                                      double const* const* const* cartesianMatrix,
                                                      const MolDS_base::Molecule& molecule, 
                                                      double const* const* orbitalElectronPopulation,
                                                      double const* const* overlapAOs) const{
   int groundState = 0;
   // dipole moment of excited states
   for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
      int excitedState = k+1; // (k+1)-th excited state
      this->CalcElectronicTransitionDipoleMoment(electronicTransitionDipoleMoments[excitedState][excitedState],
                                                 excitedState,
                                                 excitedState,
                                                 fockMatrix,
                                                 matrixCIS,
                                                 cartesianMatrix,
                                                 molecule,
                                                 orbitalElectronPopulation,
                                                 overlapAOs,
                                                 electronicTransitionDipoleMoments[groundState][groundState]);
   }
}

void ZindoS::CalcElectronicTransitionDipoleMoments(double*** electronicTransitionDipoleMoments,
                                                   double const* const* fockMatrix,
                                                   double const* const* matrixCIS,
                                                   double const* const* const* cartesianMatrix,
                                                   const MolDS_base::Molecule& molecule, 
                                                   double const* const* orbitalElectronPopulation,
                                                   double const* const* overlapAOs) const{
   int groundState = 0;
   // transition dipole moments from ground state to excited states
   for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
      int excitedState = k+1; // (k+1)-th excited state
      this->CalcElectronicTransitionDipoleMoment(electronicTransitionDipoleMoments[excitedState][groundState],
                                                 excitedState,
                                                 groundState,
                                                 fockMatrix,
                                                 matrixCIS,
                                                 cartesianMatrix,
                                                 molecule,
                                                 orbitalElectronPopulation,
                                                 overlapAOs,
                                                 electronicTransitionDipoleMoments[groundState][groundState]);
   }

   if(Parameters::GetInstance()->RequiresAllTransitionDipoleMomentsCIS()){
      // transition dipole moments between excited states
      for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
         int departureExcitedState = k+1; // (k+1)-th excited state
         for(int l=k+1; l<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); l++){
            int destinationExcitedState = l+1; // (l+1)-th excited state
            this->CalcElectronicTransitionDipoleMoment(electronicTransitionDipoleMoments[destinationExcitedState][departureExcitedState],
                                                       destinationExcitedState,
                                                       departureExcitedState,
                                                       fockMatrix,
                                                       matrixCIS,
                                                       cartesianMatrix,
                                                       molecule,
                                                       orbitalElectronPopulation,
                                                       overlapAOs,
                                                       electronicTransitionDipoleMoments[groundState][groundState]);
         }
      }
      for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS()+1; k++){
         for(int l=k+1; l<Parameters::GetInstance()->GetNumberExcitedStatesCIS()+1; l++){
            for(int axis=0; axis<CartesianType_end; axis++){
               electronicTransitionDipoleMoments[k][l][axis] 
                  = electronicTransitionDipoleMoments[l][k][axis];
            }
         }
      }
   }
}

void ZindoS::CalcElectronicTransitionDipoleMoment(double* transitionDipoleMoment,
                                                  int to, int from,
                                                  double const* const* fockMatrix,
                                                  double const* const* matrixCIS,
                                                  double const* const* const* cartesianMatrix,
                                                  const MolDS_base::Molecule& molecule, 
                                                  double const* const* orbitalElectronPopulation,
                                                  double const* const* overlapAOs,
                                                  double const* groundStateDipole) const{
   double valueX = 0.0;
   double valueY = 0.0;
   double valueZ = 0.0;
   double const* xyzCOC = molecule.GetXyzCOC();
   int groundState = 0;
   int totalNumberAOs = molecule.GetTotalNumberAOs();
   stringstream ompErrors;
   if(Parameters::GetInstance()->GetNumberExcitedStatesCIS() < from ||
      Parameters::GetInstance()->GetNumberExcitedStatesCIS() < to ){
      stringstream ss;
      ss << this->errorMessageCalcElectronicTransitionDipoleMomentBadState;
      ss << this->errorMessageFromState << from << endl;
      ss << this->errorMessageToState << to << endl;
      throw MolDSException(ss.str());
   }

   if(from != to){
      if(from != groundState && to != groundState){
         // transition dipole moment between different excited states
#pragma omp parallel for reduction(+:valueX,valueY,valueZ) schedule(auto)
         for(int l=0; l<this->matrixCISdimension; l++){
            try{
               double temp  = 0.0;
               double tempX = 0.0;
               double tempY = 0.0;
               double tempZ = 0.0;
               // single excitation from I-th (occupied)MO to A-th (virtual)MO
               int moI = this->GetActiveOccIndex(molecule, l);
               int moA = this->GetActiveVirIndex(molecule, l);
               for(int mu=0; mu<totalNumberAOs; mu++){
                  for(int nu=0; nu<totalNumberAOs; nu++){
                     temp   = (-1.0*fockMatrix[moI][mu]*fockMatrix[moI][nu] + fockMatrix[moA][mu]*fockMatrix[moA][nu]);
                     tempX += temp*(cartesianMatrix[XAxis][mu][nu] - xyzCOC[XAxis]*overlapAOs[mu][nu]);
                     tempY += temp*(cartesianMatrix[YAxis][mu][nu] - xyzCOC[YAxis]*overlapAOs[mu][nu]);
                     tempZ += temp*(cartesianMatrix[ZAxis][mu][nu] - xyzCOC[ZAxis]*overlapAOs[mu][nu]);
                  }
               }
               temp    = matrixCIS[from-1][l]*matrixCIS[to-1][l];
               valueX += temp*tempX;
               valueY += temp*tempY;
               valueZ += temp*tempZ;
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
         transitionDipoleMoment[XAxis] = valueX;
         transitionDipoleMoment[YAxis] = valueY;
         transitionDipoleMoment[ZAxis] = valueZ;
      }
      else if(from == groundState && to != groundState){
         // transition dipole moment from the ground to excited states
#pragma omp parallel for reduction(+:valueX,valueY,valueZ) schedule(auto)
         for(int l=0; l<this->matrixCISdimension; l++){
            try{
               double temp  = 0.0;
               double tempX = 0.0;
               double tempY = 0.0;
               double tempZ = 0.0;
               // single excitation from I-th (occupied)MO to A-th (virtual)MO
               int moI = this->GetActiveOccIndex(molecule, l);
               int moA = this->GetActiveVirIndex(molecule, l);
               for(int mu=0; mu<totalNumberAOs; mu++){
                  for(int nu=0; nu<totalNumberAOs; nu++){
                     temp   = fockMatrix[moA][mu]*fockMatrix[moI][nu];
                     tempX += temp*(cartesianMatrix[XAxis][mu][nu] - xyzCOC[XAxis]*overlapAOs[mu][nu]);
                     tempY += temp*(cartesianMatrix[YAxis][mu][nu] - xyzCOC[YAxis]*overlapAOs[mu][nu]);
                     tempZ += temp*(cartesianMatrix[ZAxis][mu][nu] - xyzCOC[ZAxis]*overlapAOs[mu][nu]);
                  }
               }
               temp    = this->matrixCIS[to-1][l]*sqrt(2.0);
               valueX += temp*tempX;
               valueY += temp*tempY;
               valueZ += temp*tempZ;
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
         transitionDipoleMoment[XAxis] = valueX;
         transitionDipoleMoment[YAxis] = valueY;
         transitionDipoleMoment[ZAxis] = valueZ;
      }
      else if(from != groundState && to == groundState){
         // transition dipole moment from the excited to ground states
#pragma omp parallel for reduction(+:valueX,valueY,valueZ) schedule(auto)
         for(int l=0; l<this->matrixCISdimension; l++){
            try{
               double temp  = 0.0;
               double tempX = 0.0;
               double tempY = 0.0;
               double tempZ = 0.0;
               // single excitation from I-th (occupied)MO to A-th (virtual)MO
               int moI = this->GetActiveOccIndex(molecule, l);
               int moA = this->GetActiveVirIndex(molecule, l);
               for(int mu=0; mu<totalNumberAOs; mu++){
                  for(int nu=0; nu<totalNumberAOs; nu++){
                     temp   = fockMatrix[moI][mu]*fockMatrix[moA][nu];
                     tempX += temp*(cartesianMatrix[XAxis][mu][nu] - xyzCOC[XAxis]*overlapAOs[mu][nu]);
                     tempY += temp*(cartesianMatrix[YAxis][mu][nu] - xyzCOC[YAxis]*overlapAOs[mu][nu]);
                     tempZ += temp*(cartesianMatrix[ZAxis][mu][nu] - xyzCOC[ZAxis]*overlapAOs[mu][nu]);
                  }
               }
               temp    = matrixCIS[from-1][l]*sqrt(2.0);
               valueX += temp*tempX;
               valueY += temp*tempY;
               valueZ += temp*tempZ;
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
         transitionDipoleMoment[XAxis] = valueX;
         transitionDipoleMoment[YAxis] = valueY;
         transitionDipoleMoment[ZAxis] = valueZ;
      }
   }
   else{
      if(from != groundState){
         // dipole moment of the excited state. It is needed that the dipole of ground state has been already calculated!!
#pragma omp parallel for reduction(+:valueX,valueY,valueZ) schedule(auto)
         for(int l=0; l<this->matrixCISdimension; l++){
            try{
               double temp  = 0.0;
               double tempX = 0.0;
               double tempY = 0.0;
               double tempZ = 0.0;
               // single excitation from I-th (occupied)MO to A-th (virtual)MO
               int moI = this->GetActiveOccIndex(molecule, l);
               int moA = this->GetActiveVirIndex(molecule, l);
               for(int mu=0; mu<totalNumberAOs; mu++){
                  for(int nu=0; nu<totalNumberAOs; nu++){
                     temp   = (-1.0*fockMatrix[moI][mu]*fockMatrix[moI][nu] + fockMatrix[moA][mu]*fockMatrix[moA][nu]);
                     tempX += temp*(cartesianMatrix[XAxis][mu][nu] - xyzCOC[XAxis]*overlapAOs[mu][nu]);
                     tempY += temp*(cartesianMatrix[YAxis][mu][nu] - xyzCOC[YAxis]*overlapAOs[mu][nu]);
                     tempZ += temp*(cartesianMatrix[ZAxis][mu][nu] - xyzCOC[ZAxis]*overlapAOs[mu][nu]);
                  }
               }
               temp    = matrixCIS[from-1][l]*matrixCIS[to-1][l];
               valueX += temp*tempX;
               valueY += temp*tempY;
               valueZ += temp*tempZ;
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
         transitionDipoleMoment[XAxis] = valueX + groundStateDipole[XAxis];
         transitionDipoleMoment[YAxis] = valueY + groundStateDipole[YAxis];
         transitionDipoleMoment[ZAxis] = valueZ + groundStateDipole[ZAxis];
      }
      else{
         // dipole moment of the ground state
         Cndo2::CalcElectronicTransitionDipoleMoment(transitionDipoleMoment,
                                                     to,
                                                     from,
                                                     fockMatrix,
                                                     matrixCIS,
                                                     cartesianMatrix,
                                                     molecule,
                                                     orbitalElectronPopulation,
                                                     overlapAOs,
                                                     NULL);
      }
   }
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

void ZindoS::CalcOrbitalElectronPopulationCIS(double**** orbitalElectronPopulationCIS, 
                                              double const* const* orbitalElectronPopulation, 
                                              const MolDS_base::Molecule& molecule, 
                                              double const* const* fockMatrix,
                                              double const* const* matrixCIS) const{
   if(!Parameters::GetInstance()->RequiresMullikenCIS()){
      return;
   }
   vector<int>* elecStates = Parameters::GetInstance()->GetElectronicStateIndecesMullikenCIS();
   // malloc or initialize free exciton energies
   if(*orbitalElectronPopulationCIS == NULL){
      MallocerFreer::GetInstance()->Malloc<double>(orbitalElectronPopulationCIS, 
                                                   elecStates->size(),
                                                   molecule.GetTotalNumberAOs(),
                                                   molecule.GetTotalNumberAOs());
   }
   else{
      MallocerFreer::GetInstance()->Initialize<double>(*orbitalElectronPopulationCIS, 
                                                       elecStates->size(),
                                                       molecule.GetTotalNumberAOs(),
                                                       molecule.GetTotalNumberAOs());
   }
   // clac orbital electron population
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   for(int k=0; k<elecStates->size(); k++){
      int excitedStateIndex = (*elecStates)[k]-1;
      stringstream ompErrors;
#pragma omp parallel for schedule(auto)
      for(int mu=0; mu<molecule.GetTotalNumberAOs(); mu++){
         try{
            for(int nu=0; nu<molecule.GetTotalNumberAOs(); nu++){
               double value = orbitalElectronPopulation[mu][nu];
               for(int moI=0; moI<numberActiveOcc; moI++){
                  for(int moA=numberOcc; moA<numberOcc+numberActiveVir; moA++){
                     int slaterDeterminantIndex = this->GetSlaterDeterminantIndex(moI,moA);
                     value += pow(matrixCIS[excitedStateIndex][slaterDeterminantIndex],2.0)
                             *(-fockMatrix[moI][mu]*fockMatrix[moI][nu] 
                               +fockMatrix[moA][mu]*fockMatrix[moA][nu]);
                     double tmpVal1=0.0;
                     for(int moB=numberOcc; moB<numberOcc+numberActiveVir; moB++){
                        if(moB==moA) continue;
                        int tmpSDIndex = this->GetSlaterDeterminantIndex(moI,moB);
                        tmpVal1 += matrixCIS[excitedStateIndex][tmpSDIndex]*fockMatrix[moB][nu];
                     }
                     double tmpVal2=0.0;
                     for(int moJ=0; moJ<numberActiveOcc; moJ++){
                        if(moJ==moI) continue;
                        int tmpSDIndex = this->GetSlaterDeterminantIndex(moJ,moA);
                        tmpVal2 += matrixCIS[excitedStateIndex][tmpSDIndex]*fockMatrix[moJ][mu];
                     }
                     value += matrixCIS[excitedStateIndex][slaterDeterminantIndex]
                             *(fockMatrix[moA][mu]*tmpVal1 + fockMatrix[moI][nu]*tmpVal2);
                  }
               }
               (*orbitalElectronPopulationCIS)[k][mu][nu] = value;
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
}

void ZindoS::CalcAtomicElectronPopulationCIS(double*** atomicElectronPopulationCIS,
                                             double const* const* const* orbitalElectronPopulationCIS, 
                                             const Molecule& molecule) const{
   if(!Parameters::GetInstance()->RequiresMullikenCIS()){
      return;
   }
   int totalNumberAtoms = molecule.GetNumberAtoms();
   vector<int>* elecStates = Parameters::GetInstance()->GetElectronicStateIndecesMullikenCIS();
   // malloc or initialize free exciton energies
   if(*atomicElectronPopulationCIS == NULL){
      MallocerFreer::GetInstance()->Malloc<double>(atomicElectronPopulationCIS, 
                                                   elecStates->size(),
                                                   totalNumberAtoms);
   }
   else{
      MallocerFreer::GetInstance()->Initialize<double>(*atomicElectronPopulationCIS,
                                                       elecStates->size(),
                                                       totalNumberAtoms);
   }
   // clac atomic electron population
   for(int k=0; k<elecStates->size(); k++){
      stringstream ompErrors;
#pragma omp parallel for schedule(auto)
      for(int a=0; a<totalNumberAtoms; a++){
         try{
            int firstAOIndex = molecule.GetAtom(a)->GetFirstAOIndex();
            int numberAOs = molecule.GetAtom(a)->GetValenceSize();
            for(int i=firstAOIndex; i<firstAOIndex+numberAOs; i++){
               (*atomicElectronPopulationCIS)[k][a] += orbitalElectronPopulationCIS[k][i][i];
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
}

void ZindoS::CalcAtomicUnpairedPopulationCIS(double*** atomicUnpairedPopulationCIS,
                                             double const* const* const* orbitalElectronPopulationCIS, 
                                             const Molecule& molecule) const{
   if(!Parameters::GetInstance()->RequiresMullikenCIS()){
      return;
   }
   if(!Parameters::GetInstance()->RequiresUnpairedPopCIS()){
      return;
   }
   int totalNumberAtoms = molecule.GetNumberAtoms();
   vector<int>* elecStates = Parameters::GetInstance()->GetElectronicStateIndecesMullikenCIS();
   // malloc or initialize free exciton energies
   if(*atomicUnpairedPopulationCIS == NULL){
      MallocerFreer::GetInstance()->Malloc<double>(atomicUnpairedPopulationCIS, 
                                                   elecStates->size(),
                                                   totalNumberAtoms);
   }
   else{
      MallocerFreer::GetInstance()->Initialize<double>(*atomicUnpairedPopulationCIS,
                                                       elecStates->size(),
                                                       totalNumberAtoms);
   }
   // calc atomic electron population
   for(int k=0; k<elecStates->size(); k++){
      stringstream ompErrors;
#pragma omp parallel for schedule(auto)
      for(int a=0; a<totalNumberAtoms; a++){
         try{
            int firstAOIndex = molecule.GetAtom(a)->GetFirstAOIndex();
            int numberAOs = molecule.GetAtom(a)->GetValenceSize();
            (*atomicUnpairedPopulationCIS)[k][a] = 0.0; 
            for(int i=firstAOIndex; i<firstAOIndex+numberAOs; i++){
               double orbitalSquarePopulation = 0.0; 
               int    totalNumberAOs = molecule.GetTotalNumberAOs(); 
               for(int j=0; j<totalNumberAOs; j++) {
                  orbitalSquarePopulation += orbitalElectronPopulationCIS[k][i][j] * orbitalElectronPopulationCIS[k][j][i]; 
               }
               (*atomicUnpairedPopulationCIS)[k][a] += 2.0 * orbitalElectronPopulationCIS[k][i][i] - orbitalSquarePopulation;
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
}

void ZindoS::OutputCISResults() const{
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();

   // output cis eigen energies
   this->OutputLog(this->messageExcitedStatesEnergiesTitle);
   double eV2AU = Parameters::GetInstance()->GetEV2AU();
   for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
      this->OutputLog(boost::format("%s\t%d\t%e\t%e\t") 
         % this->messageExcitedStatesEnergies
         % (k+1) 
         % this->excitedEnergies[k]
         % (this->excitedEnergies[k]/eV2AU));

      // sort eigen vector coefficeits of CIS and output
      vector<CISEigenVectorCoefficient> cisEigenVectorCoefficients;
      this->SortCISEigenVectorCoefficients(&cisEigenVectorCoefficients, this->matrixCIS[k]);
      for(int l=0; l<Parameters::GetInstance()->GetNumberPrintCoefficientsCIS(); l++){
         this->OutputLog(boost::format("%e (%d -> %d)\t") % cisEigenVectorCoefficients[l].coefficient
                                                          % cisEigenVectorCoefficients[l].occIndex
                                                          % cisEigenVectorCoefficients[l].virIndex);
      }
      this->OutputLog("\n");
   }
   this->OutputLog("\n");

   // output dipole moment
   this->OutputCISDipole();

   // output transition dipole moment
   this->OutputCISTransitionDipole();

   // output mulliken population
   this->OutputCISMulliken();

   // output unpaired electron population
   this->OutputCISUnpairedPop(); 

   // output exciton energies
   if(Parameters::GetInstance()->RequiresExcitonEnergiesCIS()){
      this->OutputLog(this->messageExcitonEnergiesCIS);
      this->OutputLog(this->messageExcitonEnergiesCISTitle);
      for(int k=0; k<Parameters::GetInstance()->GetNumberExcitedStatesCIS(); k++){
         this->OutputLog(boost::format("%s\t%d\t%e\t%e\t%e\t%e\n") 
                           % this->messageExcitonEnergiesShortCIS
                           % (k+1) 
                           %  this->freeExcitonEnergiesCIS[k]
                           % (this->freeExcitonEnergiesCIS[k]/eV2AU)
                           % (this->excitedEnergies[k]-this->freeExcitonEnergiesCIS[k])
                           %((this->excitedEnergies[k]-this->freeExcitonEnergiesCIS[k])/eV2AU));
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
      this->OutputLog(boost::format("\t%s\t%d\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n") 
         % this->messageTotalDipoleMoment
         % k
         % (this->electronicTransitionDipoleMoments[k][k][XAxis]+this->coreDipoleMoment[XAxis])
         % (this->electronicTransitionDipoleMoments[k][k][YAxis]+this->coreDipoleMoment[YAxis])
         % (this->electronicTransitionDipoleMoments[k][k][ZAxis]+this->coreDipoleMoment[ZAxis])
         % magnitude
         % ((this->electronicTransitionDipoleMoments[k][k][XAxis]+this->coreDipoleMoment[XAxis])/debye2AU)
         % ((this->electronicTransitionDipoleMoments[k][k][YAxis]+this->coreDipoleMoment[YAxis])/debye2AU)
         % ((this->electronicTransitionDipoleMoments[k][k][ZAxis]+this->coreDipoleMoment[ZAxis])/debye2AU)
         % (magnitude/debye2AU));
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
      this->OutputLog(boost::format("\t%s\t%d\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n") 
         % this->messageElectronicDipoleMoment
         % k
         % (this->electronicTransitionDipoleMoments[k][k][XAxis])
         % (this->electronicTransitionDipoleMoments[k][k][YAxis])
         % (this->electronicTransitionDipoleMoments[k][k][ZAxis])
         % magnitude
         % (this->electronicTransitionDipoleMoments[k][k][XAxis]/debye2AU)
         % (this->electronicTransitionDipoleMoments[k][k][YAxis]/debye2AU)
         % (this->electronicTransitionDipoleMoments[k][k][ZAxis]/debye2AU)
         % (magnitude/debye2AU));
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
         if(from < to){
            double magnitude = 0.0; 
            double temp = 0.0;
            temp += pow(this->electronicTransitionDipoleMoments[to][from][XAxis],2.0);
            temp += pow(this->electronicTransitionDipoleMoments[to][from][YAxis],2.0);
            temp += pow(this->electronicTransitionDipoleMoments[to][from][ZAxis],2.0);
            magnitude = sqrt(temp);
            this->OutputLog(boost::format("\t%s\t%d -> %d\t\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n") 
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
               % (magnitude/debye2AU));
         }
      }
   }
   this->OutputLog("\n");
}

void ZindoS::OutputCISMulliken() const{
   if(!Parameters::GetInstance()->RequiresMullikenCIS()){
      return;
   }
   int totalNumberAtoms = this->molecule->GetNumberAtoms();
   this->OutputLog(this->messageMullikenAtomsTitle);
   vector<int>* elecStates = Parameters::GetInstance()->GetElectronicStateIndecesMullikenCIS();
   for(int k=0; k<elecStates->size(); k++){
      for(int a=0; a<totalNumberAtoms; a++){
         const Atom& atom = *this->molecule->GetAtom(a);
         this->OutputLog(boost::format("%s\t%d\t%d\t%s\t%e\t%e\n") % this->messageMullikenAtoms
                                                                   % (*elecStates)[k]
                                                                   % a
                                                                   % AtomTypeStr(atom.GetAtomType())
                                                                   % atom.GetCoreCharge()
                                                                   % (atom.GetCoreCharge()-this->atomicElectronPopulationCIS[k][a]));
      }
      this->OutputLog("\n");
   }
}

void ZindoS::OutputCISUnpairedPop() const{
   if(!Parameters::GetInstance()->RequiresMullikenCIS()){
      return;
   }
   if(!Parameters::GetInstance()->RequiresUnpairedPopCIS()){
      return;
   }
   int totalNumberAtoms = this->molecule->GetNumberAtoms();
   this->OutputLog(this->messageUnpairedAtomsTitle);
   vector<int>* elecStates = Parameters::GetInstance()->GetElectronicStateIndecesMullikenCIS();
   for(int k=0; k<elecStates->size(); k++){
      for(int a=0; a<totalNumberAtoms; a++){
         const Atom& atom = *this->molecule->GetAtom(a);
         this->OutputLog(boost::format("%s\t%d\t%d\t%s\t%e\n") % this->messageUnpairedAtoms
                                                               % (*elecStates)[k]
                                                               % a
                                                               % AtomTypeStr(atom.GetAtomType())
                                                               % this->atomicUnpairedPopulationCIS[k][a]);
      }
      this->OutputLog("\n");
   }
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
      double overlapAOs=0.0;
      for(int j=0; j<this->matrixCISdimension; j++){
         overlapAOs += expansionVectors[j][k] * newExpansionVector[j];
      }
      for(int j=0; j<this->matrixCISdimension; j++){
         newExpansionVector[j] -= overlapAOs*expansionVectors[j][k];
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
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
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
      // Davidson loop
      for(int k=0; k<maxIter; k++){
         this->OutputLog(boost::format("%s%d\n") % this->messageNumIterCIS.c_str() % k );
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
               this->OutputLog(boost::format("\t  %d%s%e\n") % (i+1) % this->messageResidualNorm.c_str() % norm );
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

         // stop the Davidson loop
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

      }// end Davidson loop
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

            // Fast algorithm, but this is not easy to read. Slow algorithm is also written below.
            if(k<l){
               // Off diagonal term (right upper)
               matrixCIS[k][l] = this->GetCISOffDiagElement(this->nishimotoMatagaMatrix, *this->molecule, this->fockMatrix, moI, moA, moJ, moB);
            }
            else if(k==l){
               // Diagonal term
               matrixCIS[k][l] = this->GetCISDiagElement(energiesMO, this->nishimotoMatagaMatrix, *this->molecule, this->fockMatrix, moI, moA);
            } 
            // End of the fast algorith.

            /*// Slow algorith, but this is easy to read. Fast altorithm is also written above.
            double value=0.0;
            value = 2.0*this->GetMolecularIntegralElement(moA, moI, moJ, moB, 
                                                          *this->molecule, 
                                                          this->fockMatrix, 
                                                          NULL)
                       -this->GetMolecularIntegralElement(moA, moB, moI, moJ, 
                                                          *this->molecule, 
                                                          this->fockMatrix, 
                                                          NULL);
            if(k==l){
               value += this->energiesMO[moA] - this->energiesMO[moI];
            }
            matrixCIS[k][l] = value;
            // End of the slow algorith. */
         }
         catch(MolDSException ex){
#pragma omp critical
            ex.Serialize(ompErrors);
         }
      } // end of l-loop
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
                                                % this->messageDoneCalcCISMatrix.c_str() );
}

double ZindoS::GetCISDiagElement(double const* energiesMO,
                                 double const* const* const* const* nishimotoMatagaMatrix,
                                 const MolDS_base::Molecule& molecule,
                                 double const* const* fockMatrix, 
                                 int moI,
                                 int moA) const{
   double value    = energiesMO[moA] - energiesMO[moI];
   double gamma    = 0.0;
   double exchange = 0.0;
   double coulomb  = 0.0;
   int totalNumberAtoms = molecule.GetNumberAtoms();
   for(int A=0; A<totalNumberAtoms; A++){
      const Atom& atomA = *molecule.GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();

      for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
         OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
         double tmp1 = fockMatrix[moA][mu]*fockMatrix[moI][mu]; 
         double tmp2 = fockMatrix[moI][mu]*fockMatrix[moI][mu]; 
         double tmp3 = fockMatrix[moA][mu]*fockMatrix[moA][mu];

         // CNDO term
         for(int B=A; B<totalNumberAtoms; B++){
            const Atom& atomB = *molecule.GetAtom(B);
            int firstAOIndexB = atomB.GetFirstAOIndex();
            int lastAOIndexB  = atomB.GetLastAOIndex();

            for(int nu=firstAOIndexB; nu<=lastAOIndexB; nu++){
               OrbitalType orbitalNu = atomB.GetValence(nu-firstAOIndexB);

               if(A<B){
                  gamma = nishimotoMatagaMatrix[A][orbitalMu][B][orbitalNu];
                  value += 4.0*gamma*tmp1
                                    *fockMatrix[moA][nu]
                                    *fockMatrix[moI][nu];
                  value -=     gamma*tmp2
                                    *fockMatrix[moA][nu]
                                    *fockMatrix[moA][nu];
                  value -=     gamma*tmp3
                                    *fockMatrix[moI][nu]
                                    *fockMatrix[moI][nu];
               }
               else{
                  gamma = atomA.GetZindoF0ss();
                  value += 2.0*gamma*tmp1
                                    *fockMatrix[moA][nu]
                                    *fockMatrix[moI][nu];
                  value -=     gamma*tmp2
                                    *fockMatrix[moA][nu]
                                    *fockMatrix[moA][nu];
               }     
            }
         }

         //Aditional term for INDO or ZIND/S, see Eq. (10) in [RZ_1973]
         for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
            OrbitalType orbitalNu = atomA.GetValence(nu-firstAOIndexA);

            if(mu!=nu){
               exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
               value += 2.0*exchange*tmp2
                                    *fockMatrix[moA][nu]
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
               throw MolDSException(ss.str());
            }   
            value += 2.0*(coulomb-gamma)*tmp1
                                        *fockMatrix[moA][nu]
                                        *fockMatrix[moI][nu];
            value -=     (coulomb-gamma)*tmp2
                                        *fockMatrix[moA][nu]
                                        *fockMatrix[moA][nu];
         }
      }
   }
   return value;
}

double ZindoS::GetCISOffDiagElement(double const* const* const* const* nishimotoMatagaMatrix,
                                    const MolDS_base::Molecule& molecule,
                                    double const* const* fockMatrix, 
                                    int moI,
                                    int moA,
                                    int moJ,
                                    int moB) const{
   double value    = 0.0;
   double gamma    = 0.0;
   double exchange = 0.0;
   double coulomb  = 0.0;
   int totalNumberAtoms = molecule.GetNumberAtoms();
   for(int A=0; A<totalNumberAtoms; A++){
      const Atom& atomA = *molecule.GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();

      for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){ 
         OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA); 
         double tmp1 = fockMatrix[moA][mu]*fockMatrix[moI][mu]; 
         double tmp2 = fockMatrix[moA][mu]*fockMatrix[moB][mu]; 
         double tmp3 = fockMatrix[moJ][mu]*fockMatrix[moB][mu];
         double tmp4 = fockMatrix[moI][mu]*fockMatrix[moJ][mu];

         // CNDO term
         for(int B=A; B<totalNumberAtoms; B++){
            const Atom& atomB = *molecule.GetAtom(B);
            int firstAOIndexB = atomB.GetFirstAOIndex(); 
            int lastAOIndexB  = atomB.GetLastAOIndex(); 
            for(int nu=firstAOIndexB; nu<=lastAOIndexB; nu++){
               OrbitalType orbitalNu = atomB.GetValence(nu-firstAOIndexB);
               
               if(A<B){
                  gamma = nishimotoMatagaMatrix[A][orbitalMu][B][orbitalNu];
                  value += 2.0*gamma*tmp1
                                    *fockMatrix[moJ][nu]
                                    *fockMatrix[moB][nu];
                  value -=     gamma*tmp2
                                    *fockMatrix[moI][nu]
                                    *fockMatrix[moJ][nu];
                  value += 2.0*gamma*tmp3
                                    *fockMatrix[moA][nu]
                                    *fockMatrix[moI][nu];
                  value -=     gamma*tmp4
                                    *fockMatrix[moA][nu]
                                    *fockMatrix[moB][nu];
               }
               else{
                  gamma = atomA.GetZindoF0ss();
                  value += 2.0*gamma*tmp1
                                    *fockMatrix[moJ][nu]
                                    *fockMatrix[moB][nu];
                  value -=     gamma*tmp2
                                    *fockMatrix[moI][nu]
                                    *fockMatrix[moJ][nu];
               }  
            }
         }

         // Aditional term for INDO or ZIND/S, see Eq. (10) in [RZ_1973]
         double tmp5 = fockMatrix[moA][mu]*fockMatrix[moB][mu];
         double tmp6 = fockMatrix[moA][mu]*fockMatrix[moJ][mu];
         double tmp7 = fockMatrix[moA][mu]*fockMatrix[moI][mu];
         for(int nu=firstAOIndexA; nu<=lastAOIndexA; nu++){
            OrbitalType orbitalNu = atomA.GetValence(nu-firstAOIndexA);

            if(mu!=nu){
               exchange = this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
               value += 2.0*exchange
                       *(tmp5*fockMatrix[moJ][nu] + tmp6*fockMatrix[moB][nu])
                       *fockMatrix[moI][nu];
               value -= exchange
                       *(tmp6*fockMatrix[moI][nu] + tmp7*fockMatrix[moJ][nu])
                       *fockMatrix[moB][nu];
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
               throw MolDSException(ss.str());
            }   
            value += 2.0*(coulomb-gamma)*tmp7
                                        *fockMatrix[moJ][nu]
                                        *fockMatrix[moB][nu];
            value -=     (coulomb-gamma)*tmp5
                                        *fockMatrix[moI][nu]
                                        *fockMatrix[moJ][nu];
         }
      }
   }
   return value;
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

bool ZindoS::RequiresExcitedStatesForce(const vector<int>& elecStates) const{
   bool requires = true;
   if(elecStates.size()==1 && elecStates[0]==0){
      requires = false;
   }
   return requires;
}

void ZindoS::CheckZMatrixForce(const vector<int>& elecStates){
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

void ZindoS::CheckEtaMatrixForce(const vector<int>& elecStates){
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

void ZindoS::CalcEtaMatrixForce(const vector<int>& elecStates){
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
                  ex.Serialize(ompErrors);
               }
            }
            // Exception throwing for omp-region
            if(!ompErrors.str().empty()){
               throw MolDSException::Deserialize(ompErrors);
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

// see [PT_1996, PT_1997]
void ZindoS::CalcZMatrixForce(const vector<int>& elecStates){
#ifdef MOLDS_DBG
   if(this->etaMatrixForce == NULL){
      throw MolDSException(this->errorMessageCalcZMatrixForceEtaNull);
   }
#endif
   this->CheckZMatrixForce(elecStates); 

   // creat MO-index-pair for Q variables. 
   vector<MoIndexPair> nonRedundantQIndeces;
   vector<MoIndexPair> redundantQIndeces;
   this->CalcActiveSetVariablesQ(&nonRedundantQIndeces, 
                                 &redundantQIndeces,
                                 Parameters::GetInstance()->GetActiveOccCIS(),
                                 Parameters::GetInstance()->GetActiveVirCIS());

   // malloc temporary arrays
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
         if(elecStates[n] <= groundState){continue;}
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
               ex.Serialize(ompErrors);
            }
         }
         // Exception throwing for omp-region
         if(!ompErrors.str().empty()){
            throw MolDSException::Deserialize(ompErrors);
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

// each element (mu, nu) of z matrix.
// see (57) in [PT_1996]
double ZindoS::GetZMatrixForceElement(double const* y,
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

void ZindoS::MallocTempMatrixForZMatrix(double** delta,
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

void ZindoS::FreeTempMatrixForZMatrix(double** delta,
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

// see (9) in [PT_1997]
void ZindoS::CalcDeltaVector(double* delta, int exciteState) const{
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
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
}

// see variable Q-vector in [PT_1996, PT_1997]
void ZindoS::CalcActiveSetVariablesQ(vector<MoIndexPair>* nonRedundantQIndeces, 
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

// see (20) - (23) in [PT_1997]
void ZindoS::CalcQVector(double* q, 
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
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
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
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
   /* 
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
      this->OutputLog(boost::format("q[%d] = %e\n") % i % q[i]);
   }
   for(int i=0; i<redundantQIndeces.size(); i++){
      int r = nonRedundantQIndeces.size() + i;
      this->OutputLog(boost::format("q[%d] = %e\n") % r % q[r]);
   }
   */
}

// see (18) in [PT_1997]
double ZindoS::GetSmallQElement(int moI, 
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
            double rAB = this->molecule->GetDistanceAtoms(atomA, atomB);
            for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
               const OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
               for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
                  const OrbitalType orbitalLambda = atomB.GetValence(lambda-firstAOIndexB);
                  double twoElecInt = 0.0;
                  twoElecInt = this->nishimotoMatagaMatrix[A][orbitalMu][B][orbitalLambda];
                  double temp = 0.0;
                  if(isMoPOcc){
                     int p = numberOcc - (moP+1);
                     temp = 4.0*xiOcc[p][mu]    *eta[lambda][lambda]
                           -1.0*xiOcc[p][lambda]*eta[mu][lambda]
                           -1.0*xiOcc[p][lambda]*eta[mu][lambda];
                     value += twoElecInt*this->fockMatrix[moI][mu]*temp;
                     
                     temp = 4.0*xiOcc[p][lambda]*eta[mu][mu]
                           -1.0*xiOcc[p][mu]    *eta[lambda][mu]
                           -1.0*xiOcc[p][mu]    *eta[lambda][mu];
                     value += twoElecInt*this->fockMatrix[moI][lambda]*temp;
                     
                  }
                  else{
                     int p = moP - numberOcc;
                     temp = 4.0*xiVir[p][mu]    *eta[lambda][lambda]
                           -1.0*xiVir[p][lambda]*eta[lambda][mu]
                           -1.0*xiVir[p][lambda]*eta[lambda][mu];
                     value += twoElecInt*this->fockMatrix[moI][mu]*temp;
                     
                     temp = 4.0*xiVir[p][lambda]*eta[mu][mu]
                           -1.0*xiVir[p][mu]    *eta[mu][lambda]
                           -1.0*xiVir[p][mu]    *eta[mu][lambda];
                     value += twoElecInt*this->fockMatrix[moI][lambda]*temp;
                     
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

void ZindoS::CalcXiMatrices(double** xiOcc, 
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
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
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
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
}

// right hand side of (54) in [PT_1996]      
void ZindoS::CalcAuxiliaryVector(double* y, 
                               double const* q, 
                               double const* const* kRDagerGammaRInv, 
                               const vector<MoIndexPair>& nonRedundantQIndeces, 
                               const vector<MoIndexPair>& redundantQIndeces) const{
   MallocerFreer::GetInstance()->Initialize<double>(
                                 y,
                                 nonRedundantQIndeces.size());
   MolDS_wrappers::Blas::GetInstance()->Dgemv(nonRedundantQIndeces.size(),
                                              redundantQIndeces.size(),
                                              kRDagerGammaRInv,
                                              &(q[nonRedundantQIndeces.size()]),
                                              y);
   stringstream ompErrors;
#pragma omp parallel
#pragma omp single nowait
   for(int i=0; i<nonRedundantQIndeces.size(); i++){
#pragma omp task
      {
         try{
            int moI = nonRedundantQIndeces[i].moI;
            int moJ = nonRedundantQIndeces[i].moJ;
            y[i] += q[i]/this->GetNNRElement(moI, moJ, moI, moJ);
         }
         catch(MolDSException ex){
#pragma omp critical
            ex.Serialize(ompErrors);
         }
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
}

// see (40) and (45) in [PT_1996].
// This method calculates "\Gamma_{NR} - K_{NR}" to solve (54) in [PT_1966]
// Note taht K_{NR} is not calculated.
void ZindoS::CalcGammaNRMinusKNRMatrix(double** gammaNRMinusKNR, const vector<MoIndexPair>& nonRedundantQIndeces) const{
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
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
}

// see (41), (42), and (46) in [PT_1996].
// This method calculates "K_{R}^{\dagger} * Gamma_{R}" matrix, see (41), (42), and (46) to solve (54) in [PT_1996]
// Note taht K_{R}^{\dager} is not calculated.
void ZindoS::CalcKRDagerGammaRInvMatrix(double** kRDagerGammaRInv, 
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
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
}

// see (40) in [PT_1996]
double ZindoS::GetGammaNRElement(int moI, int moJ, int moK, int moL) const{
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
double ZindoS::GetGammaRElement(int moI, int moJ, int moK, int moL) const{
   double value=0.0;
   if(moI==moK && moJ==moL){
      value = moI==moJ ? 1.0 : this->energiesMO[moJ]-this->energiesMO[moI];
   }
   return value;
}

// see (43) in [PT_1996]
double ZindoS::GetNNRElement(int moI, int moJ, int moK, int moL) const{
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
double ZindoS::GetNRElement(int moI, int moJ, int moK, int moL) const{
   double value=0.0;
   if(moI==moK && moJ==moL){
      value = 1.0;
   }
   return value;
}

// see (44) in [PT_1996]
double ZindoS::GetKNRElement(int moI, int moJ, int moK, int moL) const{
   double value=0.0;
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int nI = moI<numberOcc ? 2 : 0;
   int nJ = moJ<numberOcc ? 2 : 0;
   int nK = moK<numberOcc ? 2 : 0;
   int nL = moL<numberOcc ? 2 : 0;
   
   if(nI!=nJ && nK!=nL){
      value = this->GetAuxiliaryKNRKRElement(moI, moJ, moK, moL);
   }
   //See (24) in [DL_1990] about "0.5" multiplied to "GetKNRElement".
   return 0.5*value;
}

// Dager of (45) in [PT_1996]. Note taht the (45) is real number.
double ZindoS::GetKRDagerElement(int moI, int moJ, int moK, int moL) const{
   return this->GetKRElement(moK, moL, moI, moJ);
}

// see (45) in [PT_1996]
double ZindoS::GetKRElement(int moI, int moJ, int moK, int moL) const{
   double value=0.0;
   int numberOcc = this->molecule->GetTotalNumberValenceElectrons()/2;
   int nI = moI<numberOcc ? 2 : 0;
   int nJ = moJ<numberOcc ? 2 : 0;
   int nK = moK<numberOcc ? 2 : 0;
   int nL = moL<numberOcc ? 2 : 0;

   if(nI==nJ && nK!=nL){
      value = this->GetAuxiliaryKNRKRElement(moI, moJ, moK, moL);
   }
   //See (24) in [DL_1990] about "0.5" multiplied to "GetKRElement".
   return 0.5*value;
}

// see common term in eqs. (45) and (46) in [PT_1996],
// that is, 4.0(ij|kl) - (ik|jl) - (il|jk).
double ZindoS::GetAuxiliaryKNRKRElement(int moI, int moJ, int moK, int moL) const{
   double value = 0.0;

   // Fast algorith, but this is not easy to read. 
   // Slow algorithm is alos written below.
   for(int A=0; A<this->molecule->GetNumberAtoms(); A++){
      const Atom& atomA = *this->molecule->GetAtom(A);
      int firstAOIndexA = atomA.GetFirstAOIndex();
      int lastAOIndexA  = atomA.GetLastAOIndex();

      for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
         OrbitalType orbitalMu = atomA.GetValence(mu - firstAOIndexA);
         double tmpM01 = this->fockMatrix[moI][mu]
                        *this->fockMatrix[moJ][mu];
         double tmpM02 = this->fockMatrix[moI][mu]
                        *this->fockMatrix[moL][mu];
         double tmpM03 = this->fockMatrix[moJ][mu]
                        *this->fockMatrix[moK][mu];
         double tmpM04 = this->fockMatrix[moI][mu]
                        *this->fockMatrix[moK][mu];
         double tmpM05 = this->fockMatrix[moJ][mu]
                        *this->fockMatrix[moL][mu];
         double tmpM06 = this->fockMatrix[moK][mu]
                        *this->fockMatrix[moL][mu];

         // A=B && (mu==lambda && nu==sigma for (mu nu|lamba sigma))
         for(int nu=mu+1; nu<=lastAOIndexA; nu++){
            OrbitalType orbitalNu = atomA.GetValence(nu - firstAOIndexA);
            double tmpValue = 0.0;
            tmpValue -= 4.0*tmpM01 
                       *this->fockMatrix[moK][nu]
                       *this->fockMatrix[moL][nu];
            tmpValue += 6.0*tmpM02
                       *this->fockMatrix[moJ][nu]
                       *this->fockMatrix[moK][nu];
            tmpValue += 6.0*tmpM03
                       *this->fockMatrix[moL][nu]
                       *this->fockMatrix[moI][nu];
            tmpValue += 6.0*tmpM04
                       *this->fockMatrix[moJ][nu]
                       *this->fockMatrix[moL][nu];
            tmpValue += 6.0*tmpM05
                       *this->fockMatrix[moI][nu]
                       *this->fockMatrix[moK][nu];
            tmpValue -= 4.0*tmpM06
                       *this->fockMatrix[moI][nu]
                       *this->fockMatrix[moJ][nu];

            double gamma = 0.5*this->GetExchangeInt(orbitalMu, orbitalNu, atomA);
            value += tmpValue*gamma;
         }

         //  (A==B || A!=B) && (mu==nu && lambda==sigma for (mu nu|lamba sigma))
         for(int B=A; B<this->molecule->GetNumberAtoms(); B++){
            const Atom& atomB = *this->molecule->GetAtom(B);
            int firstAOIndexB = atomB.GetFirstAOIndex();
            int lastAOIndexB  = atomB.GetLastAOIndex();

            for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
               OrbitalType orbitalLambda = atomB.GetValence(lambda - firstAOIndexB);
               double tmpValue = 0.0;
               tmpValue += 4.0*tmpM01
                          *this->fockMatrix[moK][lambda]
                          *this->fockMatrix[moL][lambda];
               tmpValue -= tmpM02
                          *this->fockMatrix[moJ][lambda]
                          *this->fockMatrix[moK][lambda];
               tmpValue -= tmpM03
                          *this->fockMatrix[moI][lambda]
                          *this->fockMatrix[moL][lambda];
               tmpValue -= tmpM04
                          *this->fockMatrix[moJ][lambda]
                          *this->fockMatrix[moL][lambda];
               tmpValue -= tmpM05
                          *this->fockMatrix[moI][lambda]
                          *this->fockMatrix[moK][lambda];
               tmpValue += 4.0*tmpM06
                          *this->fockMatrix[moI][lambda]
                          *this->fockMatrix[moJ][lambda];
               double gamma = 0.0;
               if(A!=B){
                  gamma = this->nishimotoMatagaMatrix[A][orbitalMu][B][orbitalLambda];
               }
               else{
                  gamma = 0.5*this->GetCoulombInt(orbitalMu, orbitalLambda, atomA);
               }
               value += tmpValue*gamma;
            }
         }
      }
   }
   // End of the fast algorith.

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

void ZindoS::CalcDiatomicTwoElecTwoCore1stDerivatives(double*** matrix,  
                                                      int indexAtomA, 
                                                      int indexAtomB) const{
   const Atom& atomA = *molecule->GetAtom(indexAtomA);
   const int firstAOIndexA = atomA.GetFirstAOIndex();
   const int lastAOIndexA  = atomA.GetLastAOIndex();
   const Atom& atomB = *molecule->GetAtom(indexAtomB);
   const int firstAOIndexB = atomB.GetFirstAOIndex();
   const int lastAOIndexB  = atomB.GetLastAOIndex();
   const double rAB = this->molecule->GetDistanceAtoms(atomA, atomB);
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      const OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
      for(int nu=firstAOIndexB; nu<=lastAOIndexB; nu++){
         const OrbitalType orbitalNu = atomB.GetValence(nu-firstAOIndexB);
         for(int i=0; i<CartesianType_end; i++){
            matrix[mu-firstAOIndexA][nu-firstAOIndexB][i]
               = this->GetNishimotoMatagaTwoEleInt1stDerivative(
                     atomA, orbitalMu, atomB, orbitalNu, rAB, static_cast<CartesianType>(i));
         }
      }
   }
}                                                      


// elecStates is indeces of the electroinc eigen states.
// The index = 0 means electronic ground state. 
void ZindoS::CalcForce(const vector<int>& elecStates){
   this->CheckMatrixForce(elecStates);
   if(this->RequiresExcitedStatesForce(elecStates)){
      this->CalcEtaMatrixForce(elecStates);
      this->CalcZMatrixForce(elecStates);
   }
   stringstream ompErrors;
#pragma omp parallel 
   {
      double*** diatomicTwoElecTwoCore1stDerivs = NULL;
      double*** diatomicOverlapAOs1stDerivs = NULL;
      try{
         MallocerFreer::GetInstance()->Malloc<double>(&diatomicTwoElecTwoCore1stDerivs,
                                                      OrbitalType_end, 
                                                      OrbitalType_end, 
                                                      CartesianType_end);
         MallocerFreer::GetInstance()->Malloc<double>(&diatomicOverlapAOs1stDerivs,
                                                      OrbitalType_end, 
                                                      OrbitalType_end, 
                                                      CartesianType_end);

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
               double rAB = this->molecule->GetDistanceAtoms(atomA, atomB);

               // calc. first derivative of overlapAOs.
               this->CalcDiatomicOverlapAOs1stDerivatives(diatomicOverlapAOs1stDerivs, atomA, atomB);

               // calc. first derivative of two elec two core interaction by Nishimoto-Mataga
               this->CalcDiatomicTwoElecTwoCore1stDerivatives(diatomicTwoElecTwoCore1stDerivs, a, b);

               double coreRepulsion       [CartesianType_end] = {0.0,0.0,0.0};
               double forceElecCoreAttPart[CartesianType_end] = {0.0,0.0,0.0};
               for(int i=0; i<CartesianType_end; i++){
                  // core repulsion part (ground state)
                  coreRepulsion[i] = this->GetDiatomCoreRepulsion1stDerivative(a, b, static_cast<CartesianType>(i));
                  if(Parameters::GetInstance()->RequiresVdWSCF()){
                     coreRepulsion[i] += this->GetDiatomVdWCorrection1stDerivative(a, b, static_cast<CartesianType>(i));
                  }
                  // electron core attraction part (ground state)                     
                  forceElecCoreAttPart[i] = ( atomA.GetCoreCharge()*atomicElectronPopulation[b]
                                             +atomB.GetCoreCharge()*atomicElectronPopulation[a])
                                           *diatomicTwoElecTwoCore1stDerivs[s][s][i];
               }
               double forceOverlapAOsPart [CartesianType_end] = {0.0,0.0,0.0};
               double forceTwoElecPart    [CartesianType_end] = {0.0,0.0,0.0};
               for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
                  OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                  for(int nu=firstAOIndexB; nu<=lastAOIndexB; nu++){
                     OrbitalType orbitalNu = atomB.GetValence(nu-firstAOIndexB);
                     double bondParameter = 0.5*(atomA.GetBondingParameter(this->theory, orbitalMu) 
                                                +atomB.GetBondingParameter(this->theory, orbitalNu)); 
                     for(int i=0; i<CartesianType_end; i++){
                        // overlapAOs part (ground state)
                        forceOverlapAOsPart[i] += 2.0*this->orbitalElectronPopulation[mu][nu]
                                                      *bondParameter
                                                      *diatomicOverlapAOs1stDerivs[mu-firstAOIndexA][nu-firstAOIndexB][i];
                        // two electron part (ground state)
                        forceTwoElecPart[i] += (this->orbitalElectronPopulation[mu][mu]
                                               *this->orbitalElectronPopulation[nu][nu]
                                               -0.5*pow(this->orbitalElectronPopulation[mu][nu],2.0))
                                              *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA][nu-firstAOIndexB][i];
                     }
                  }
               }
               // sum up contributions from each part (ground state)
#pragma omp critical
               {
                  for(int n=0; n<elecStates.size(); n++){
                     for(int i=0; i<CartesianType_end; i++){
                        this->matrixForce[n][a][i] += -coreRepulsion[i]
                                                      +forceElecCoreAttPart[i]
                                                      -forceOverlapAOsPart[i]
                                                      -forceTwoElecPart[i];
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
                  // electron core attraction part (excited state)
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
               }
            } // end of for(int b)
         }    // end of for(int a)
      }       // end of try
      catch(MolDSException ex){
#pragma omp critical
         ex.Serialize(ompErrors);
      }
      MallocerFreer::GetInstance()->Free<double>(&diatomicTwoElecTwoCore1stDerivs, 
                                                 OrbitalType_end,
                                                 OrbitalType_end,
                                                 CartesianType_end);
      MallocerFreer::GetInstance()->Free<double>(&diatomicOverlapAOs1stDerivs, 
                                                 OrbitalType_end,
                                                 OrbitalType_end,
                                                 CartesianType_end);
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }
  
   /*
   // Calculate force (on the ground state only). 
   // First derivative of overlapAOs integral is
   // calculated with GTO expansion technique.
   stringstream ompErrors;
#pragma omp parallel for schedule(auto)
   for(int a=0; a<this->molecule->GetNumberAtoms(); a++){
      try{
         const Atom& atomA = *molecule->GetAtom(a);
         int firstAOIndexA = atomA.GetFirstAOIndex();
         int lastAOIndexA  = atomA.GetLastAOIndex();
         for(int i=0; i<CartesianType_end; i++){

            double coreRepulsion = 0.0;
            double forceElecCoreAttPart = 0.0;
            double forceOverlapAOsPart = 0.0;
            double forceTwoElecPart = 0.0;
            for(int b=0; b<this->molecule->GetNumberAtoms(); b++){
               if(a != b){
                  const Atom& atomB = *molecule->GetAtom(b);
                  int firstAOIndexB = atomB.GetFirstAOIndex();
                  int lastAOIndexB  = atomB.GetLastAOIndex();

                  // Calculation of core repusion force
                  coreRepulsion += this->GetDiatomCoreRepulsion1stDerivative(
                                         a, b, (CartesianType)i);
                  if(Parameters::GetInstance()->RequiresVdWSCF()){
                     coreRepulsion += this->GetDiatomVdWCorrection1stDerivative(
                                            a, b, (CartesianType)i);
                  }

                  // Calculate force arise from electronic part.
                  forceElecCoreAttPart += ( atomA.GetCoreCharge()*atomicElectronPopulation[b]
                                       +atomB.GetCoreCharge()*atomicElectronPopulation[a])
                                       *this->GetNishimotoMatagaTwoEleInt1stDerivative
                                             (atomA, s, atomB, s, (CartesianType)i);

                  for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
                     OrbitalType orbitalMu = atomA.GetValence(mu-firstAOIndexA);
                     for(int nu=firstAOIndexB; nu<=lastAOIndexB; nu++){
                        OrbitalType orbitalNu = atomB.GetValence(nu-firstAOIndexB);

                        double bondParameter = 0.5*(atomA.GetBondingParameter(this->theory, 
                                                                              orbitalMu) 
                                                   +atomB.GetBondingParameter(this->theory, 
                                                                              orbitalNu)); 

                        forceOverlapAOsPart += 2.0*this->orbitalElectronPopulation[mu][nu]
                                            *bondParameter
                                            *this->GetOverlapAOsElement1stDerivativeByGTOExpansion
                                                   (atomA, mu-firstAOIndexA, 
                                                    atomB, nu-firstAOIndexB,
                                                    STO6G, (CartesianType)i);

                        forceTwoElecPart += (this->orbitalElectronPopulation[mu][mu]
                                            *this->orbitalElectronPopulation[nu][nu]
                                            -0.5*pow(this->orbitalElectronPopulation[mu][nu],2.0))
                                            *this->GetNishimotoMatagaTwoEleInt1stDerivative
                                                   (atomA, orbitalMu, atomB, orbitalNu,
                                                   (CartesianType)i);
                     }
                  }
               }
            }

            this->matrixForce[0][a][i] = -1.0*(coreRepulsion 
                                              -forceElecCoreAttPart 
                                              +forceOverlapAOsPart
                                              +forceTwoElecPart);
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
   */
}

void ZindoS::CalcForceExcitedStaticPart(double* force, 
                                      int elecStateIndex,
                                      int indexAtomA, 
                                      int indexAtomB,
                                      double const* const* const* diatomicTwoElecTwoCore1stDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int lastAOIndexA  = atomA.GetLastAOIndex();
   int lastAOIndexB  = atomB.GetLastAOIndex();
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
         for(int i=0; i<CartesianType_end; i++){
            double temp= 2.0*this->etaMatrixForce[elecStateIndex][mu][mu]
                            *this->etaMatrixForce[elecStateIndex][lambda][lambda]
                        -1.0*this->etaMatrixForce[elecStateIndex][mu][lambda]
                            *this->etaMatrixForce[elecStateIndex][mu][lambda];
            force[i] += temp
                       *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA]
                                                       [lambda-firstAOIndexB]
                                                       [i];
         }
      }
   }
}

void ZindoS::CalcForceExcitedElecCoreAttractionPart(double* force, 
                                                  int elecStateIndex,
                                                  int indexAtomA, 
                                                  int indexAtomB,
                                                  double const* const* const* diatomicTwoElecTwoCore1stDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int lastAOIndexA  = atomA.GetLastAOIndex();
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      for(int i=0; i<CartesianType_end; i++){
         force[i] += this->zMatrixForce[elecStateIndex][mu][mu]
                    *atomB.GetCoreCharge()
                    *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA][s][i];
      }
   }
}

void ZindoS::CalcForceExcitedOverlapAOsPart(double* force, 
                                            int elecStateIndex,
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
         bondParameter *= 0.5;
         for(int i=0; i<CartesianType_end; i++){
            force[i] += -1.0
                       *this->zMatrixForce[elecStateIndex][mu][nu]
                       *bondParameter
                       *diatomicOverlapAOs1stDerivs[mu-firstAOIndexA][nu-firstAOIndexB][i];
         }
      }
   }
}

void ZindoS::CalcForceExcitedTwoElecPart(double* force, 
                                         int elecStateIndex,
                                         int indexAtomA, 
                                         int indexAtomB,
                                         double const* const* const* diatomicTwoElecTwoCore1stDerivs) const{
   const Atom& atomA = *this->molecule->GetAtom(indexAtomA);
   const Atom& atomB = *this->molecule->GetAtom(indexAtomB);
   int firstAOIndexA = atomA.GetFirstAOIndex();
   int firstAOIndexB = atomB.GetFirstAOIndex();
   int lastAOIndexA  = atomA.GetLastAOIndex();
   int lastAOIndexB  = atomB.GetLastAOIndex();
   for(int mu=firstAOIndexA; mu<=lastAOIndexA; mu++){
      for(int lambda=firstAOIndexB; lambda<=lastAOIndexB; lambda++){
         for(int i=0; i<CartesianType_end; i++){
            force[i] -= this->zMatrixForce[elecStateIndex][mu][mu]
                       *this->orbitalElectronPopulation[lambda][lambda]
                       *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA]
                                                       [lambda-firstAOIndexB]
                                                       [i];
            force[i] += 0.50
                       *this->zMatrixForce[elecStateIndex][mu][lambda]
                       *this->orbitalElectronPopulation[mu][lambda]
                       *diatomicTwoElecTwoCore1stDerivs[mu-firstAOIndexA]
                                                       [lambda-firstAOIndexB]
                                                       [i];
         }
      }
   }
}


}

