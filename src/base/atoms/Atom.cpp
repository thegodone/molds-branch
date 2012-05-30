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
#include<new>
#include<string>
#include<vector>
#include<stdexcept>
#include"../PrintController.h"
#include"../MolDSException.h"
#include"../Uncopyable.h"
#include"../Enums.h"
#include"../MathUtilities.h"
#include"../MallocerFreer.h"
#include"../EularAngle.h"
#include"Atom.h"
using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{
Atom::Atom(){
   this->SetMessages();
   this->xyz = NULL;
   this->pxyz = NULL;
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&this->xyz, CartesianType_end);
      MallocerFreer::GetInstance()->Malloc<double>(&this->pxyz, CartesianType_end);
   }
   catch(exception ex){
      MallocerFreer::GetInstance()->Free<double>(&this->xyz, CartesianType_end);
      MallocerFreer::GetInstance()->Free<double>(&this->pxyz, CartesianType_end);
      throw MolDSException(ex.what());
   }
}

Atom::~Atom(){
   MallocerFreer::GetInstance()->Free<double>(&this->xyz, CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(&this->pxyz, CartesianType_end);
   //this->OutputLog("atom deleted\n");
}

void Atom::SetMessages(){
   this->errorMessageOrbitalExponent = "Error in base_atoms::Atom::GetOrbitalExponent: Invalid shelltype or orbitalType.\n";
   this->errorMessageCndo2CoreIntegral = "Error in base_atoms::Atom::GetCoreIntegral: Invalid orbitalType for Cndo2.\n";
   this->errorMessageIndoCoreIntegral = "Error in base_atoms::Atom::GetCoreIntegral: Invalid orbitalType for INDO.\n";
   this->errorMessageMndoCoreIntegral = "Error in base_atoms::Atom::GetMndoCoreINtegral: Invalid orbitalType for MNDO.\n";
   this->errorMessageAm1CoreIntegral = "Error in base_atoms::Atom::GetAm1CoreINtegral: Invalid orbitalType for AM1.\n";
   this->errorMessageAm1DCoreIntegral = "Error in base_atoms::Atom::GetAm1DCoreINtegral: Invalid orbitalType for AM1.\n";
   this->errorMessagePm3CoreIntegral = "Error in base_atoms::Atom::GetPm3CoreINtegral: Invalid orbitalType for PM3.\n";
   this->errorMessagePm3DCoreIntegral = "Error in base_atoms::Atom::GetPm3DCoreINtegral: Invalid orbitalType for PM3-D.\n";
   this->errorMessagePm3PddgCoreIntegral = "Error in base_atoms::Atom::GetPm3PddgCoreINtegral: Invalid orbitalType for PM3/PDDG.\n";
   this->errorMessageIonPot = "Error in base_atoms::Atom::GetZindoIonPot: Invalid orbitalType.\n";
   this->errorMessageAtomType = "\tatom type = ";
   this->errorMessageOrbitalType = "\torbital type = ";
   this->errorMessageShellType = "\tshell type = ";
   this->errorMessageTheoryType = "\tTheory = ";
   this->errorMessageNumberValences = "\tnumber of valences = ";
   this->errorMessageValenceIndex = "\tvalenceIndex = ";
   this->errorMessageGetAtomicBasisValueBadValenceIndex 
      = "Error in molds_atoms::Atom::GetAtomicBasisValue: Bad valenceIndex is set.\n";
   this->errorMessageGetRealAngularPartAOBadValence 
      = "Error in molds_atoms::Atom::GetRealAngularPartAO: Bad valence orbital is set.\n";
   this->errorMessageEffectivPrincipalQuantumNumber = 
      "Error in base::Atom::GetEffectivePrincipalQuantumNumber: invalid shelltype.\n";
   this->errorMessageZindoCoreIntegral = "Error in base_atoms::Atom::GetZindoCoreINtegral: Invalid orbitalType.\n";
   this->errorMessageGetOrbitalExponentBadTheory = "Erro in base_atoms::Atom::GetOrbitalExponent: Bad theory is set.\n";
   this->errorMessageGetBondingParameterBadTheoryBadOrbital = "Error in base_atoms::Atom::GetBondingParameter: Bad Theory of bad orbital is set.\n";
   this->errorMessageGetNddoDerivedParameterDBadDIndex 
      = "Error in base_atoms::Atom::GetNddoDerivedParameterD: Bad index for parameter D(dIndex). Only 0, 1, and 2 are permitted.\n";
   this->errorMessageGetNddoDerivedParameterDBadTheory
      = "Error in base_atoms::Atom::GetNddoDerivedParameterD: Bad theory is set.\n";
   this->errorMessageGetNddoAlphaBadTheory
      = "Error in base_atoms::Atom::GetNddoAlpha: Bad theory is set.\n";
   this->errorMessageDIndex  = "dIndex = ";
   this->errorMessageGetNddoDerivedParameterRhoBadRhoIndex 
      = "Error in base_atoms::Atom::GetNddoDerivedParameterRho: Bad index for parameter rho(rhoIndex). Only 0, 1, and 2 are permitted.\n";
   this->errorMessageRhoIndex = "rhoIndex = ";
   this->errorMessageGetNddoDerivedParameterRhoBadTheory 
      = "Error in base_atoms::Atom::GetNddoDerivedParameterRho: Bad thory is set.\n";
   this->errorMessageGetNddoParameterKBadKIndex 
      = "Error in base_atoms::Atom::GetNddoParameterK: Bad index for parameter K(kIndex). Only 0, 1, 2, and 3 are permitted.\n";
   this->errorMessageGetNddoParameterKBadTheory
      = "Error in base_atoms::Atom::GetNddoParameterK: Bad theory is set.\n";
   this->errorMessageKIndex  = "kIndex = ";
   this->errorMessageGetNddoParameterLBadLIndex 
      = "Error in base_atoms::Atom::GetNddoParameterL: Bad index for parameter L(lIndex). Only 0, 1, 2, and 3 are permitted.\n";
   this->errorMessageGetNddoParameterLBadTheory
      = "Error in base_atoms::Atom::GetNddoParameterL: Bad theory is set.\n";
   this->errorMessageLIndex  = "lIndex = ";
   this->errorMessageGetNddoParameterMBadMIndex 
      = "Error in base_atoms::Atom::GetNddoParameterM: Bad index for parameter M(mIndex). Only 0, 1, 2, and 3 are permitted.\n";
   this->errorMessageGetNddoParameterMBadTheory
      = "Error in base_atoms::Atom::GetNddoParameterM: Bad theory is set.\n";
   this->errorMessageMIndex  = "mIndex = ";
   this->errorMessageGetNddoGssBadTheory 
      = "Error in base_atoms::Atom::GetNddoGss: Bad theory is set.\n";
   this->errorMessageGetNddoGppBadTheory 
      = "Error in base_atoms::Atom::GetNddoGpp Bad theory is set.\n";
   this->errorMessageGetNddoGspBadTheory 
      = "Error in base_atoms::Atom::GetNddoGsp: Bad theory is set.\n";
   this->errorMessageGetNddoGpp2BadTheory 
      = "Error in base_atoms::Atom::GetNddoGpp2: Bad theory is set.\n";
   this->errorMessageGetNddoHspBadTheory 
      = "Error in base_atoms::Atom::GetNddoHsp: Bad theory is set.\n";
   this->errorMessageGetNddoHppBadTheory 
      = "Error in base_atoms::Atom::GetNddoHp: Bad theory is set.\n";
   this->errorMessageGetPm3PddgParameterPaBadPaIndex 
      = "Error in base_atoms::Atom::GetPm3PddgParameterPa: Bad index for parameter Pa(paIndex). Only 0, and 1 are permitted.\n";
   this->errorMessagePaIndex  = "paIndex = ";
   this->errorMessageGetPm3PddgParameterDaBadDaIndex 
      = "Error in base_atoms::Atom::GetPm3PddgParameterDa: Bad index for parameter Da(daIndex). Only 0, and 1 are permitted.\n";
   this->errorMessageDaIndex  = "daIndex = ";
   this->errorMessageGetXyzCoordinatesNull = "Error in base_atoms::Atom::GetXyz: xyz is NULL\n";
   this->errorMessageGetPxyzMomentaNull  = "Error in base_atoms::Atom::GetPxyz: pxyz is NULL\n";
}

AtomType Atom::GetAtomType() const{
   return this->atomType;
}

double Atom::GetAtomicMass() const{
   return this->atomicMass;
}

double Atom::GetCoreMass() const{
   return this->atomicMass - static_cast<double>(this->numberValenceElectrons);
}

double* Atom::GetXyz() const{
   if(this->xyz==NULL){
      stringstream ss;
      ss << this->errorMessageGetXyzCoordinatesNull;
      throw MolDSException(ss.str());
   }
   return this->xyz;
}

double* Atom::GetPxyz() const{
   if(this->pxyz==NULL){
      stringstream ss;
      ss << this->errorMessageGetPxyzMomentaNull;
      throw MolDSException(ss.str());
   }
   return this->pxyz;
}

void Atom::SetXyz(double x, double y, double z) const{
   xyz[0]= x;
   xyz[1]= y;
   xyz[2]= z;
}

void Atom::SetPxyz(double px, double py, double pz) const{
   pxyz[0]= px;
   pxyz[1]= py;
   pxyz[2]= pz;
}

int Atom::GetValenceSize() const{
   return this->valence.size();
}

OrbitalType Atom::GetValence(int index) const{
   return this->valence[index];
}

double Atom::GetVdWCoefficient() const{
   return this->vdWCoefficient;
}

double Atom::GetVdWRadii() const{
   return this->vdWRadii;
}

double Atom::GetAtomicBasisValue(double x, 
                                 double y, 
                                 double z, 
                                 int valenceIndex, 
                                 TheoryType theory) const{
   if(this->valence.size()<=valenceIndex){
      stringstream ss;
      ss << this->errorMessageGetAtomicBasisValueBadValenceIndex;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageNumberValences << this->valence.size() << endl;
      ss << this->errorMessageValenceIndex << valenceIndex << endl;
      throw MolDSException(ss.str());
   }
   double dx = x - this->xyz[XAxis];
   double dy = y - this->xyz[YAxis];
   double dz = z - this->xyz[ZAxis];
   double dr = sqrt( pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0) );
   EularAngle eularAngle(dx, dy, dz);
   double angularPart = this->GetRealAngularPartAO(eularAngle.GetBeta(),
                                                   eularAngle.GetAlpha(),
                                                   this->valence[valenceIndex]);
   double orbitalExponent = this->GetOrbitalExponent(this->valenceShellType,
                                                     this->valence[valenceIndex],
                                                     theory);
   double radialPart = this->GetRadialPartAO(dr, orbitalExponent, this->valenceShellType);
   return angularPart*radialPart;
}

// See (1.74) & (1.72) in J. A. Pople book.
double Atom::GetRadialPartAO(double dr, 
                             double orbitalExponent, 
                             MolDS_base::ShellType shell) const{
   int principalQuantumNumber = static_cast<int>(shell) + 1;
   double temp1 = pow(2.0*orbitalExponent,static_cast<double>(principalQuantumNumber)+0.5);
   double temp2 = pow(Factorial(2*principalQuantumNumber),-0.5);
   return temp1*temp2*pow(dr,principalQuantumNumber-1)*exp(-1.0*orbitalExponent*dr);
}

// See Table 1 in [BFB_1997] or Table 1.2 in J. A. Pople book.
// See Table 1 in [BFB_1997] or p25 in J. A. Pople book for defenitions of theta and phi.
double Atom::GetRealAngularPartAO(double theta, 
                                  double phi, 
                                  OrbitalType orbital) const{
   double value=0.0;
   switch(orbital){
      case s:
         value = pow(4.0*M_PI,-0.5);
         break;
      case py:
         value = pow(3.0/(4.0*M_PI),0.5)*sin(theta)*sin(phi);
         break;
      case pz:
         value = pow(3.0/(4.0*M_PI),0.5)*cos(theta);
         break;
      case px:
         value = pow(3.0/(4.0*M_PI),0.5)*sin(theta)*cos(phi);
         break;
      case dxy:
         value = pow(15.0/(16.0*M_PI),0.5)*pow(sin(theta),2.0)*sin(2.0*phi);
         break;
      case dyz:
         value = pow(15.0/(16.0*M_PI),0.5)*sin(2.0*theta)*sin(phi);
         break;
      case dzz:
         value = pow(5.0/(16.0*M_PI),0.5)*(3.0*pow(cos(theta),2.0) - 1.0);
         break;
      case dzx:
         value = pow(15.0/(16.0*M_PI),0.5)*sin(2.0*theta)*cos(phi);
         break;
      case dxxyy:
         value = pow(15.0/(16.0*M_PI),0.5)*pow(sin(theta),2.0)*cos(2.0*phi);
         break;
      default:
         stringstream ss;
         ss << this->errorMessageGetRealAngularPartAOBadValence;
         ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
         ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
         throw MolDSException(ss.str());
   }
   return value;
}


double Atom::GetBondingParameter(TheoryType theory, OrbitalType orbital) const{

   double value = 0.0;
   if(theory == CNDO2 || theory == INDO){
      value = this->bondingParameter;
   }     
   else if(theory == ZINDOS && ( orbital == s ||
                                 orbital == px ||
                                 orbital == py ||
                                 orbital == pz ) ){
      value = this->zindoBondingParameterS;
   }
   else if(theory == ZINDOS && ( orbital == dxy ||
                                 orbital == dyz ||
                                 orbital == dzz ||
                                 orbital == dzx ||
                                 orbital == dxxyy ) ){
      value = this->zindoBondingParameterD;
   }
   else if(theory == MNDO && orbital == s){
      value = this->mndoBondingParameterS;
   }
   else if(theory == MNDO && ( orbital == px ||
                               orbital == py ||
                               orbital == pz ) ){
      value = this->mndoBondingParameterP;
   }
   else if(theory == AM1 && orbital == s){
      value = this->am1BondingParameterS;
   }
   else if(theory == AM1 && ( orbital == px ||
                              orbital == py ||
                              orbital == pz ) ){
      value = this->am1BondingParameterP;
   }
   else if(theory == AM1D && orbital == s){
      value = this->am1DBondingParameterS;
   }
   else if(theory == AM1D && ( orbital == px ||
                               orbital == py ||
                               orbital == pz ) ){
      value = this->am1DBondingParameterP;
   }
   else if(theory == PM3 && orbital == s){
      value = this->pm3BondingParameterS;
   }
   else if(theory == PM3 && ( orbital == px ||
                              orbital == py ||
                              orbital == pz ) ){
      value = this->pm3BondingParameterP;
   }
   else if(theory == PM3D && orbital == s){
      value = this->pm3DBondingParameterS;
   }
   else if(theory == PM3D && ( orbital == px ||
                               orbital == py ||
                               orbital == pz ) ){
      value = this->pm3DBondingParameterP;
   }
   else if(theory == PM3PDDG && orbital == s){
      value = this->pm3PddgBondingParameterS;
   }
   else if(theory == PM3PDDG && ( orbital == px ||
                                  orbital == py ||
                                  orbital == pz ) ){
      value = this->pm3PddgBondingParameterP;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetBondingParameterBadTheoryBadOrbital;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << "\n";
      throw MolDSException(ss.str());
   }

   return value;

}

double Atom::GetBondingParameter() const{
   return this->GetBondingParameter(CNDO2, s);
}

double Atom::GetCoreCharge() const{
   return this->coreCharge;
}

int Atom::GetFirstAOIndex() const{
   return this->firstAOIndex;
}

void Atom::SetFirstAOIndex(int firstAOIndex){
   this->firstAOIndex = firstAOIndex;
}

ShellType Atom::GetValenceShellType() const{
   return this->valenceShellType;
}

int Atom::GetEffectivePrincipalQuantumNumber(ShellType shellType) const{
   if(shellType == k){
      return 1.0;
   }
   else if(shellType == l){
      return 2.0;
   }
   else if(shellType == m){
      return 3.0;
   }
   else{
      stringstream ss;
      ss << this->errorMessageEffectivPrincipalQuantumNumber;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
      ss << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
      throw MolDSException(ss.str());
   }
}   

int Atom::GetNumberValenceElectrons() const{
   return this->numberValenceElectrons;
}

// (1.73) in J. A. Pople book
double Atom::GetOrbitalExponent(ShellType shellType, 
                                OrbitalType orbitalType, 
                                TheoryType theory) const{
   if(theory == CNDO2 || theory == INDO || theory == ZINDOS){
      if(shellType == k && orbitalType == s){ 
         return this->effectiveNuclearChargeK
               /this->GetEffectivePrincipalQuantumNumber(shellType);
      }   
      else if(shellType == l && (orbitalType == s  || 
                                 orbitalType == px || 
                                 orbitalType == py || 
                                 orbitalType == pz)){
         return this->effectiveNuclearChargeL
               /this->GetEffectivePrincipalQuantumNumber(shellType);
      }   
      else if(shellType == m && (orbitalType == s  || 
                                 orbitalType == px || 
                                 orbitalType == py || 
                                 orbitalType == pz )){
         return this->effectiveNuclearChargeMsp
               /this->GetEffectivePrincipalQuantumNumber(shellType);
      }   
      else if(shellType == m && (orbitalType == dxy  || 
                                 orbitalType == dyz ||
                                 orbitalType == dzz ||
                                 orbitalType == dzx ||
                                 orbitalType == dxxyy)){
         return this->effectiveNuclearChargeMd
               /this->GetEffectivePrincipalQuantumNumber(shellType);
      }   
      else{
         stringstream ss;
         ss << this->errorMessageOrbitalExponent;
         ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
         ss << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
         ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << "\n";
         throw MolDSException(ss.str());
      }   
   }
   else if(theory == MNDO){
      if(orbitalType == s){ 
         return this->mndoOrbitalExponentS;
      }
      else if(orbitalType == px ||
              orbitalType == py ||
              orbitalType == pz){
         return this->mndoOrbitalExponentP;
      }
      else{
         stringstream ss;
         ss << this->errorMessageOrbitalExponent;
         ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
         ss << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
         ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else if(theory == AM1 || theory == AM1D){
      if(orbitalType == s){ 
         return this->am1OrbitalExponentS;
      }
      else if(orbitalType == px ||
              orbitalType == py ||
              orbitalType == pz){
         return this->am1OrbitalExponentP;
      }
      else{
         stringstream ss;
         ss << this->errorMessageOrbitalExponent;
         ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
         ss << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
         ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else if(theory == PM3 || theory == PM3D){
      if(orbitalType == s){ 
         return this->pm3OrbitalExponentS;
      }
      else if(orbitalType == px ||
              orbitalType == py ||
              orbitalType == pz){
         return this->pm3OrbitalExponentP;
      }
      else{
         stringstream ss;
         ss << this->errorMessageOrbitalExponent;
         ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
         ss << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
         ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else if(theory == PM3PDDG){
      if(orbitalType == s){ 
         return this->pm3PddgOrbitalExponentS;
      }
      else if(orbitalType == px ||
              orbitalType == py ||
              orbitalType == pz){
         return this->pm3PddgOrbitalExponentP;
      }
      else{
         stringstream ss;
         ss << this->errorMessageOrbitalExponent;
         ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
         ss << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
         ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetOrbitalExponentBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
}


// Part of Eq. (13) in [BZ_1979]
double Atom::GetZindoJss() const{
   return this->zindoF0ss;
}

// Part of Eq. (13) in [BZ_1979]
double Atom::GetZindoJsp() const{
   // F0ss = F0sp
   return this->zindoF0ss - this->zindoG1sp/6.0;
}

// Part of Eq. (13) in [BZ_1979]
double Atom::GetZindoJsd() const{
   return this->zindoF0sd - this->zindoG2sd/10.0;
}

// Part of Eq. (13) in [BZ_1979]
double Atom::GetZindoJpp() const{
   // F0pp = F0ss
   return this->zindoF0ss - 2.0*this->zindoF2pp/25.0;
}

// Part of Eq. (13) in [BZ_1979]
double Atom::GetZindoJpd() const{
   // F0pd = F0sd
   return this->zindoF0sd - this->zindoG1pd/15.0 - 3.0*this->zindoG3pd/70.0;
}

// Part of Eq. (13) in [BZ_1979]
double Atom::GetZindoJdd() const{
   return this->zindoF0dd - 2.0*(this->zindoF2dd + this->zindoF4dd)/63.0;
}

// (3.72) in J. A. Pople book.
double Atom::GetCndo2CoreIntegral(OrbitalType orbital, double gamma, bool isGuess) const{
   double value = 0.0;
   if(orbital == s){
      value = -1.0*this->imuAmuS;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = -1.0*this->imuAmuP;
   }
   else if(orbital == dxy || 
           orbital == dyz || 
           orbital == dzz || 
           orbital == dzx || 
           orbital == dxxyy ){
      value = -1.0*this->imuAmuD;
   }   
   else{
      stringstream ss;
      ss << this->errorMessageCndo2CoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }
   if(!isGuess){
      value -= (this->coreCharge - 0.5)*gamma;
   }
   return value;
}

// (3.93) - (3.99) in J. A. Pople book.
double Atom::GetIndoCoreIntegral(OrbitalType orbital, double gamma, bool isGuess) const{
   double value = 0.0;
   if(orbital == s){
      value = -1.0*this->imuAmuS;
      if(!isGuess){
         value -= this->indoF0CoefficientS*gamma 
                 +this->indoG1CoefficientS*this->indoG1
                 +this->indoF2CoefficientS*this->indoF2;
      }
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = -1.0*this->imuAmuP;
      if(!isGuess){
         value -= this->indoF0CoefficientP*gamma 
                 +this->indoG1CoefficientP*this->indoG1
                 +this->indoF2CoefficientP*this->indoF2;
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageIndoCoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }
   return value;
}


// Eq. (13) in [BZ_1979]
double Atom::GetZindoCoreIntegral(OrbitalType orbital) const{
   double value=0.0;

   if(orbital == s){
      value = -1.0*this->zindoIonPotS 
              - this->GetZindoJss()*static_cast<double>(this->zindoL-1) 
              - this->GetZindoJsp()*static_cast<double>(this->zindoM)
              - this->GetZindoJsd()*static_cast<double>(this->zindoN);
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = -1.0*this->zindoIonPotP
              - this->GetZindoJpp()*static_cast<double>(this->zindoM-1) 
              - this->GetZindoJsp()*static_cast<double>(this->zindoL)
              - this->GetZindoJpd()*static_cast<double>(this->zindoN);
   }
   else if(orbital == dxy || orbital == dyz || orbital == dzz || orbital == dzx || orbital == dxxyy ){
      value = -1.0*this->zindoIonPotD
              - this->GetZindoJdd()*static_cast<double>(this->zindoN-1) 
              - this->GetZindoJsd()*static_cast<double>(this->zindoL)
              - this->GetZindoJpd()*static_cast<double>(this->zindoM);
   }
   else{
      stringstream ss;
      ss << this->errorMessageZindoCoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}

double Atom::GetMndoCoreIntegral(OrbitalType orbital) const{
   double value=0.0;

   if(orbital == s){
      value = this->mndoCoreintegralS;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = this->mndoCoreintegralP;
   }
   else{
      stringstream ss;
      ss << this->errorMessageMndoCoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}

double Atom::GetAm1CoreIntegral(OrbitalType orbital) const{
   double value=0.0;

   if(orbital == s){
      value = this->am1CoreintegralS;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = this->am1CoreintegralP;
   }
   else{
      stringstream ss;
      ss << this->errorMessageAm1CoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}

double Atom::GetAm1DCoreIntegral(OrbitalType orbital) const{
   double value=0.0;

   if(orbital == s){
      value = this->am1DCoreintegralS;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = this->am1DCoreintegralP;
   }
   else{
      stringstream ss;
      ss << this->errorMessageAm1DCoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}

double Atom::GetPm3CoreIntegral(OrbitalType orbital) const{
   double value=0.0;

   if(orbital == s){
      value = this->pm3CoreintegralS;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = this->pm3CoreintegralP;
   }
   else{
      stringstream ss;
      ss << this->errorMessagePm3CoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}

double Atom::GetPm3DCoreIntegral(OrbitalType orbital) const{
   double value=0.0;

   if(orbital == s){
      value = this->pm3DCoreintegralS;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = this->pm3DCoreintegralP;
   }
   else{
      stringstream ss;
      ss << this->errorMessagePm3DCoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}

double Atom::GetPm3PddgCoreIntegral(OrbitalType orbital) const{
   double value=0.0;

   if(orbital == s){
      value = this->pm3PddgCoreintegralS;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = this->pm3PddgCoreintegralP;
   }
   else{
      stringstream ss;
      ss << this->errorMessagePm3PddgCoreIntegral;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}

double Atom::GetIndoF2() const{
   return this->indoF2;
}

double Atom::GetIndoG1() const{
   return this->indoG1;
}

double Atom::GetNddoAlpha(TheoryType theory) const{
   double value = 0.0;
   if(theory == MNDO){
      value = this->mndoAlpha;
   }
   else if(theory == AM1){
      value = this->am1Alpha;
   }
   else if(theory == AM1D){
      value = this->am1DAlpha;
   }
   else if(theory == PM3){
      value = this->pm3Alpha;
   }
   else if(theory == PM3D){
      value = this->pm3DAlpha;
   }
   else if(theory == PM3PDDG){
      value = this->pm3PddgAlpha;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoAlphaBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
   return value;
}

double Atom::GetNddoDerivedParameterD(TheoryType theory, int dIndex) const{
   if(dIndex == 0 || dIndex == 1 || dIndex == 2){
      if(theory == MNDO){
         return this->mndoDerivedParameterD[dIndex];
      }
      else if(theory == AM1 || theory == AM1D){
         return this->am1DerivedParameterD[dIndex];
      }
      else if(theory == PM3 || theory == PM3D){
         return this->pm3DerivedParameterD[dIndex];
      }
      else if(theory == PM3PDDG){
         return this->pm3PddgDerivedParameterD[dIndex];
      }
      else{
         stringstream ss;
         ss << this->errorMessageGetNddoDerivedParameterDBadTheory;
         ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoDerivedParameterDBadDIndex;
      ss << this->errorMessageDIndex << dIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoDerivedParameterRho(TheoryType theory, int rhoIndex) const{
   if(rhoIndex == 0 || rhoIndex == 1 || rhoIndex == 2){
      if(theory == MNDO){
         return this->mndoDerivedParameterRho[rhoIndex];
      }
      else if(theory == AM1 || theory == AM1D){
         return this->am1DerivedParameterRho[rhoIndex];
      }
      else if(theory == PM3 || theory == PM3D){
         return this->pm3DerivedParameterRho[rhoIndex];
      }
      else if(theory == PM3PDDG){
         return this->pm3PddgDerivedParameterRho[rhoIndex];
      }
      else{
         stringstream ss;
         ss << this->errorMessageGetNddoDerivedParameterDBadTheory;
         ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoDerivedParameterRhoBadRhoIndex;
      ss << this->errorMessageRhoIndex << rhoIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoParameterK(TheoryType theory, int kIndex) const{
   if(kIndex == 0 || kIndex == 1 || kIndex == 2 || kIndex == 3){
      if(theory == AM1 || theory == AM1D){
         return this->am1ParameterK[kIndex];
      }
      else if(theory == PM3 || theory == PM3D){
         return this->pm3ParameterK[kIndex];
      }
      else if(theory == PM3PDDG){
         return this->pm3PddgParameterK[kIndex];
      }
      else{
         stringstream ss;
         ss << this->errorMessageGetNddoParameterKBadTheory;
         ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoParameterKBadKIndex;
      ss << this->errorMessageKIndex << kIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoParameterL(TheoryType theory, int lIndex) const{
   if(lIndex == 0 || lIndex == 1 || lIndex == 2 || lIndex == 3){
      if(theory == AM1 || theory == AM1D){
         return this->am1ParameterL[lIndex];
      }
      else if(theory == PM3 || theory == PM3D){
         return this->pm3ParameterL[lIndex];
      }
      else if(theory == PM3PDDG){
         return this->pm3PddgParameterL[lIndex];
      }
      else{
         stringstream ss;
         ss << this->errorMessageGetNddoParameterLBadTheory;
         ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoParameterLBadLIndex;
      ss << this->errorMessageLIndex << lIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoParameterM(TheoryType theory, int mIndex) const{
   if(mIndex == 0 || mIndex == 1 || mIndex == 2 || mIndex == 3){
      if(theory == AM1 || theory == AM1D){
         return this->am1ParameterM[mIndex];
      }
      else if(theory == PM3 || theory == PM3D){
         return this->pm3ParameterM[mIndex];
      }
      else if(theory == PM3PDDG){
         return this->pm3PddgParameterM[mIndex];
      }
      else{
         stringstream ss;
         ss << this->errorMessageGetNddoParameterMBadTheory;
         ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
         throw MolDSException(ss.str());
      }
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoParameterMBadMIndex;
      ss << this->errorMessageMIndex << mIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetPm3PddgParameterPa(int paIndex) const{
   if(paIndex == 0 || paIndex == 1 ){
      return this->pm3PddgParameterPa[paIndex];
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetPm3PddgParameterPaBadPaIndex;
      ss << this->errorMessagePaIndex << paIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetPm3PddgParameterDa(int daIndex) const{
   if(daIndex == 0 || daIndex == 1 ){
      return this->pm3PddgParameterDa[daIndex];
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetPm3PddgParameterDaBadDaIndex;
      ss << this->errorMessageDaIndex << daIndex << endl;
      throw MolDSException(ss.str());
   }
}

double Atom::GetMndoElecEnergyAtom() const{
   return this->mndoElecEnergyAtom;
}

double Atom::GetMndoHeatsFormAtom() const{
   return this->mndoHeatsFormAtom;
}

double Atom::GetNddoGss(TheoryType theory) const{
   if(theory == MNDO){
      return this->mndoGss;
   }
   else if(theory == AM1 || theory == AM1D){
      return this->am1Gss;
   }
   else if(theory == PM3 || theory == PM3D){
      return this->pm3Gss;
   }
   else if(theory == PM3PDDG){
      return this->pm3Gss;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoGssBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoGpp(TheoryType theory) const{
   if(theory == MNDO){
      return this->mndoGpp;
   }
   else if(theory == AM1 || theory == AM1D){
      return this->am1Gpp;
   }
   else if(theory == PM3 || theory == PM3D){
      return this->pm3Gpp;
   }
   else if(theory == PM3PDDG){
      return this->pm3Gpp;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoGppBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoGsp(TheoryType theory) const{
   if(theory == MNDO){
      return this->mndoGsp;
   }
   else if(theory == AM1 || theory == AM1D){
      return this->am1Gsp;
   }
   else if(theory == PM3 || theory == PM3D){
      return this->pm3Gsp;
   }
   else if(theory == PM3PDDG){
      return this->pm3Gsp;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoGspBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoGpp2(TheoryType theory) const{
   if(theory == MNDO){
      return this->mndoGpp2;
   }
   else if(theory == AM1 || theory == AM1D){
      return this->am1Gpp2;
   }
   else if(theory == PM3 || theory == PM3D){
      return this->pm3Gpp2;
   }
   else if(theory == PM3PDDG){
      return this->pm3Gpp2;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoGpp2BadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
}

double Atom::GetNddoHsp(TheoryType theory) const{
   if(theory == MNDO){
      return this->mndoHsp;
   }
   else if(theory == AM1 || theory == AM1D){
      return this->am1Hsp;
   }
   else if(theory == PM3 || theory == PM3D){
      return this->pm3Hsp;
   }
   else if(theory == PM3PDDG){
      return this->pm3Hsp;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoHspBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
}

// see p17 in [MOPAC_1990]
double Atom::GetNddoHpp(TheoryType theory) const{
   if(theory == MNDO){
      return 0.5*(this->mndoGpp - this->mndoGpp2);
   }
   else if(theory == AM1 || theory == AM1D){
      return 0.5*(this->am1Gpp - this->am1Gpp2);
   }
   else if(theory == PM3 || theory == PM3D){
      return 0.5*(this->pm3Gpp - this->pm3Gpp2);
   }
   else if(theory == PM3PDDG){
      return 0.5*(this->pm3Gpp - this->pm3Gpp2);
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetNddoHppBadTheory;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
      throw MolDSException(ss.str());
   }
}

// Table 1 in ref. [RZ_1976], Table 1 in [AEZ_1986], or Table 1 in [GD_1972]
double Atom::GetZindoF0ss() const{
   return this->zindoF0ss;
}

// Table 1 in [AEZ_1986]
double Atom::GetZindoF0sd() const{
   return this->zindoF0sd;
}

// Table 1 in [AEZ_1986]
double Atom::GetZindoF0dd() const{
   return this->zindoF0dd;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoG1sp() const{
   return this->zindoG1sp;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoF2pp() const{
   return this->zindoF2pp;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoG2sd() const{
   return this->zindoG2sd;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoG1pd() const{
   return this->zindoG1pd;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoF2pd() const{
   return this->zindoF2pd;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoG3pd() const{
   return this->zindoG3pd;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoF2dd() const{
   return this->zindoF2dd;
}

// Table 3 in ref. [BZ_1979]
double Atom::GetZindoF4dd() const{
   return this->zindoF4dd;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF0ssLower() const{
   return this->zindoF0ss;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF0sdLower() const{
   return this->zindoF0sd;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF0ddLower() const{
   return this->zindoF0dd;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoG1spLower() const{
   return this->zindoG1sp/3.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF2ppLower() const{
   return this->zindoF2pp/25.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoG2sdLower() const{
   return this->zindoG2sd/5.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoG1pdLower() const{
   return this->zindoG1pd/15.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF2pdLower() const{
   return this->zindoF2pd/35.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoG3pdLower() const{
   return this->zindoG3pd/245.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF2ddLower() const{
   return this->zindoF2dd/49.0;
}

// Apendix in ref. [BZ_1979]
double Atom::GetZindoF4ddLower() const{
   return this->zindoF4dd/441.0;
}

double Atom::GetCoreIntegral(OrbitalType orbital, 
                             double gamma, 
                             bool isGuess, 
                             TheoryType theory) const{
   double value = 0.0;
   if(theory == CNDO2){
      value = this->GetCndo2CoreIntegral(orbital, gamma, isGuess);
   }
   else if(theory == INDO){
      value = this->GetIndoCoreIntegral(orbital, gamma, isGuess);
   }
   else if(theory == ZINDOS){
      value = this->GetZindoCoreIntegral(orbital);
   }
   else if(theory == MNDO){
      value = this->GetMndoCoreIntegral(orbital);
   }
   else if(theory == AM1){
      value = this->GetAm1CoreIntegral(orbital);
   }
   else if(theory == AM1D){
      value = this->GetAm1DCoreIntegral(orbital);
   }
   else if(theory == PM3){
      value = this->GetPm3CoreIntegral(orbital);
   }
   else if(theory == PM3D){
      value = this->GetPm3DCoreIntegral(orbital);
   }
   else if(theory == PM3PDDG){
      value = this->GetPm3PddgCoreIntegral(orbital);
   }
   return value;
}

double Atom::GetCoreIntegral(OrbitalType orbital, bool isGuess, TheoryType theory) const{
   return this->GetCoreIntegral(orbital, 0.0, isGuess, theory);
}

double Atom::GetZindoIonPot(OrbitalType orbital) const{
   double value=0.0;

   if(orbital == s){
      value = -1.0*this->zindoIonPotS;
   }
   else if(orbital == px || orbital == py || orbital == pz){
      value = -1.0*this->zindoIonPotP;
   }
   else if(orbital == dxy || orbital == dyz || orbital == dzz || orbital == dzx || orbital == dxxyy ){
      value = -1.0*this->zindoIonPotD;
   }
   else{
      stringstream ss;
      ss << this->errorMessageIonPot;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << endl;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbital) << endl;
      throw MolDSException(ss.str());
   }

   return value;
}
}




















