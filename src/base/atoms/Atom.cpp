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
#include<new>
#include<string>
#include<vector>
#include<stdexcept>
#include<boost/format.hpp>
#include"../Enums.h"
#include"../Uncopyable.h"
#include"../PrintController.h"
#include"../MolDSException.h"
#include"../MallocerFreer.h"
#include"../../mpi/MpiInt.h"
#include"../../mpi/MpiProcess.h"
#include"../MathUtilities.h"
#include"../MallocerFreer.h"
#include"../EularAngle.h"
#include"../RealSphericalHarmonicsIndex.h"
#include"Atom.h"
using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{
Atom::Atom(int index){
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
   this->index = index;
}

Atom::~Atom(){
   MallocerFreer::GetInstance()->Free<double>(&this->xyz, CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(&this->pxyz, CartesianType_end);
   int valenceSize = this->valence.size();
   for(int i=0; i<valenceSize; i++){
      delete this->realSphericalHarmonicsIndeces[i];
   }
   this->realSphericalHarmonicsIndeces.clear();
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
   this->errorMessageGetNddoDerivedParameterDBadMultipoleType 
      = "Error in base_atoms::Atom::GetNddoDerivedParameterD: Bad multipole tyep for NDDO derived parameter D is set.\n";
   this->errorMessageGetNddoDerivedParameterDBadTheory
      = "Error in base_atoms::Atom::GetNddoDerivedParameterD: Bad theory is set.\n";
   this->errorMessageGetNddoAlphaBadTheory
      = "Error in base_atoms::Atom::GetNddoAlpha: Bad theory is set.\n";
   this->errorMessageMultipoleType  = "MultipoleType = ";
   this->errorMessageGetNddoDerivedParameterRhoBadMultipoleType 
      = "Error in base_atoms::Atom::GetNddoDerivedParameterRho: Bad multipole type for parameter rho is set.\n";
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
   this->errorMessageSetXyzCoordinatesNull = "Error in base_atoms::Atom::SetXyz: xyz is NULL\n";
   this->errorMessageGetPxyzMomentaNull  = "Error in base_atoms::Atom::GetPxyz: pxyz is NULL\n";
   this->errorMessageSetPxyzMomentaNull  = "Error in base_atoms::Atom::SetPxyz: pxyz is NULL\n";
}

double* Atom::GetXyz() const{
#ifdef MOLDS_DBG
   if(this->xyz==NULL) throw MolDSException(this->errorMessageGetXyzCoordinatesNull);
#endif
   return this->xyz;
}

double* Atom::GetPxyz() const{
#ifdef MOLDS_DBG
   if(this->pxyz==NULL) throw MolDSException(this->errorMessageGetPxyzMomentaNull);
#endif
   return this->pxyz;
}

void Atom::SetXyz(double x, double y, double z) const{
#ifdef MOLDS_DBG
   if(this->xyz==NULL){
      printf("xyz\n\n");
      throw MolDSException("aaa");
      //throw MolDSException(this->errorMessageSetXyzCoordinatesNull);
   }
#endif
   xyz[0]= x; xyz[1]= y; xyz[2]= z;
}

void Atom::SetPxyz(double px, double py, double pz) const{
#ifdef MOLDS_DBG
   if(this->pxyz==NULL) throw MolDSException(this->errorMessageSetPxyzMomentaNull);
#endif
   pxyz[0]= px; pxyz[1]= py; pxyz[2]= pz;
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

// Table 1.4 in J.A. Pople book
int Atom::GetEffectivePrincipalQuantumNumber(ShellType shellType) const{
   if(shellType == kShell){
      return 1.0;
   }
   else if(shellType == lShell){
      return 2.0;
   }
   else if(shellType == mShell){
      return 3.0;
   }
   else if(shellType == nShell){
      return 3.7;
   }
   else{
      stringstream ss;
      ss << this->errorMessageEffectivPrincipalQuantumNumber;
      ss << this->errorMessageAtomType << AtomTypeStr(this->atomType) << "\n";
      ss << this->errorMessageShellType << ShellTypeStr(shellType) << "\n";
      throw MolDSException(ss.str());
   }
}   

// (1.73) in J. A. Pople book
double Atom::GetOrbitalExponent(ShellType shellType, 
                                OrbitalType orbitalType, 
                                TheoryType theory) const{
   if(theory == CNDO2 || theory == INDO || theory == ZINDOS){
      if(shellType == kShell && orbitalType == s){ 
         return this->effectiveNuclearChargeK
               /this->GetEffectivePrincipalQuantumNumber(shellType);
      }   
      else if(shellType == lShell && (orbitalType == s  || 
                                      orbitalType == px || 
                                      orbitalType == py || 
                                      orbitalType == pz)){
         return this->effectiveNuclearChargeL
               /this->GetEffectivePrincipalQuantumNumber(shellType);
      }   
      else if(shellType == mShell && (orbitalType == s  || 
                                      orbitalType == px || 
                                      orbitalType == py || 
                                      orbitalType == pz )){
         return this->effectiveNuclearChargeMsp
               /this->GetEffectivePrincipalQuantumNumber(shellType);
      }   
      else if(shellType == mShell && (orbitalType == dxy  || 
                                      orbitalType == dyz ||
                                      orbitalType == dzz ||
                                      orbitalType == dzx ||
                                      orbitalType == dxxyy)){
         return this->effectiveNuclearChargeMd
               /this->GetEffectivePrincipalQuantumNumber(shellType);
      }   
      else if(shellType == nShell && (orbitalType == s  || 
                                      orbitalType == px || 
                                      orbitalType == py || 
                                      orbitalType == pz )){
         return this->effectiveNuclearChargeNsp
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

double Atom::GetNddoDerivedParameterD(TheoryType theory, MultipoleType multipole) const{
   int dIndex=0;
   switch(multipole){
      case sQ:
         dIndex = 0;
         break;
      case mux:
      case muy:
      case muz:
         dIndex = 1;
         break;
      case Qxx:
      case Qyy:
      case Qzz:
      case Qxz:
      case Qyz:
      case Qxy:
         dIndex = 2;
         break;
      default:
         stringstream ss;
         ss << this->errorMessageGetNddoDerivedParameterDBadMultipoleType;
         ss << this->errorMessageMultipoleType << MultipoleTypeStr(multipole) << endl;
         throw MolDSException(ss.str());
   }

   switch(theory){
      case MNDO:
         return this->mndoDerivedParameterD[dIndex];
         break;
      case AM1:
      case AM1D:
         return this->am1DerivedParameterD[dIndex];
         break;
      case PM3:
      case PM3D:
         return this->pm3DerivedParameterD[dIndex];
         break;
      case PM3PDDG:
         return this->pm3PddgDerivedParameterD[dIndex];
         break;
      default:
         stringstream ss;
         ss << this->errorMessageGetNddoDerivedParameterDBadTheory;
         ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
         throw MolDSException(ss.str());
   }
}

double Atom::GetNddoDerivedParameterRho(TheoryType theory, MultipoleType multipole) const{
   int rhoIndex=0;
   switch(multipole){
      case sQ:
         rhoIndex = 0;
         break;
      case mux:
      case muy:
      case muz:
         rhoIndex = 1;
         break;
      case Qxx:
      case Qyy:
      case Qzz:
      case Qxz:
      case Qyz:
      case Qxy:
         rhoIndex = 2;
         break;
      default:
         stringstream ss;
         ss << this->errorMessageGetNddoDerivedParameterRhoBadMultipoleType;
         ss << this->errorMessageMultipoleType << MultipoleTypeStr(multipole) << endl;
         throw MolDSException(ss.str());
   }

   switch(theory){
      case MNDO:
         return this->mndoDerivedParameterRho[rhoIndex];
         break;
      case AM1:
      case AM1D:
         return this->am1DerivedParameterRho[rhoIndex];
         break;
      case PM3:
      case PM3D:
         return this->pm3DerivedParameterRho[rhoIndex];
         break;
      case PM3PDDG:
         return this->pm3PddgDerivedParameterRho[rhoIndex];
         break;
      default:
         stringstream ss;
         ss << this->errorMessageGetNddoDerivedParameterRhoBadTheory;
         ss << this->errorMessageTheoryType << TheoryTypeStr(theory) << "\n";
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




















