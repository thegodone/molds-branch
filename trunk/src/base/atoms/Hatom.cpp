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
#include<vector>
#include<boost/format.hpp>
#include<boost/shared_ptr.hpp>
#include"../Enums.h"
#include"../Uncopyable.h"
#include"../PrintController.h"
#include"../MolDSException.h"
#include"../MallocerFreer.h"
#include"../../mpi/MpiInt.h"
#include"../../mpi/MpiProcess.h"
#include"../EularAngle.h"
#include"../Parameters.h"
#include"../RealSphericalHarmonicsIndex.h"
#include"Atom.h"
#include"Hatom.h"
using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{

Hatom::Hatom(int index) : Atom(index){
   this->SetAtomicParameters();
}

Hatom::~Hatom(){}

void Hatom::SetAtomicParameters(){
   this->atomType = H;
   this->atomicMass = 1.00794*Parameters::GetInstance()->GetGMolin2AU();
   this->coreCharge = 1.0;
   this->numberValenceElectrons = 1;
   this->valenceShellType = kShell;
   this->valence.push_back(s);
   for(int i=0; i<this->valence.size();i++){
      this->realSphericalHarmonicsIndeces.push_back(new RealSphericalHarmonicsIndex(this->valence[i]));
   }
   this->vdWCoefficient = 0.16*Parameters::GetInstance()->GetJ2AU()
                              *pow(Parameters::GetInstance()->GetNm2AU(),6.0)
                              /Parameters::GetInstance()->GetAvogadro();
   this->vdWRadii = 1.110*Parameters::GetInstance()->GetAngstrom2AU();
   this->bondingParameter = -9.0*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuS = 7.176*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 0.0;
   this->imuAmuD = 0.0;
   this->effectiveNuclearChargeK = 1.2; // see P78 in J. A. Pople book
   this->effectiveNuclearChargeL = 0.0;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   this->indoG1 = 0.0;
   this->indoF2 = 0.0;
   this->indoF0CoefficientS = 0.5;
   this->indoF0CoefficientP = 0.0;
   this->indoG1CoefficientS = 0.0;
   this->indoG1CoefficientP = 0.0;
   this->indoF2CoefficientS = 0.0;
   this->indoF2CoefficientP = 0.0;
   this->zindoBondingParameterS = -12.0*Parameters::GetInstance()->GetEV2AU();
   this->zindoBondingParameterD = 0.0;
   this->zindoF0ss = 12.85 * Parameters::GetInstance()->GetEV2AU();                  
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                  
   this->zindoG1sp = 0.0;                 
   this->zindoF2pp = 0.0;                 
   this->zindoG2sd = 0.0;                 
   this->zindoG1pd = 0.0;                 
   this->zindoF2pd = 0.0;                 
   this->zindoG3pd = 0.0;                 
   this->zindoF2dd = 0.0;                 
   this->zindoF4dd = 0.0;                 
   this->zindoL = 1;
   this->zindoM = 0;
   this->zindoN = 0;
   this->zindoIonPotS = 13.06 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotP = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotD = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralS = -11.906276 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = 0.0;         
   this->mndoOrbitalExponentS = 1.331967;      
   this->mndoOrbitalExponentP = 0.0;      
   this->mndoBondingParameterS = -6.989064 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = 0.0;     
   this->mndoAlpha = 2.544134 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -11.906276 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 52.102 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoGss = 12.848 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp = 0.0 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp = 0.0 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 = 0.0 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp = 0.0 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] = 0.0;    
   this->mndoDerivedParameterD[1] = 0.0;    
   this->mndoDerivedParameterD[2] = 0.0;    
   this->mndoDerivedParameterRho[0] = 0.5/0.4721515374;
   //this->mndoDerivedParameterRho[0] = 0.560345 * Parameters::GetInstance()->GetAngstrom2AU();  
   this->mndoDerivedParameterRho[1] = 0.0;  
   this->mndoDerivedParameterRho[2] = 0.0;  
   this->am1CoreintegralS = -11.396427 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = 0.0;         
   this->am1OrbitalExponentS = 1.188078;      
   this->am1OrbitalExponentP = 0.0;      
   this->am1BondingParameterS = -6.173787 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = 0.0;     
   this->am1Alpha = 2.882324 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 0.0;    
   this->am1DerivedParameterD[2] = 0.0;    
   this->am1DerivedParameterRho[0] = 0.5/0.4721515374;
   this->am1DerivedParameterRho[1] = 0.0;  
   this->am1DerivedParameterRho[2] = 0.0;  
   this->am1ParameterK[0] = 0.122796 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.005090 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] =-0.018336 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.000000 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 2.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 1.20 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 1.80 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 2.10 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1DCoreintegralS = -11.223791 * Parameters::GetInstance()->GetEV2AU();         
   this->am1DCoreintegralP = 0.0 * Parameters::GetInstance()->GetEV2AU();         
   this->am1DBondingParameterS = -6.376265 * Parameters::GetInstance()->GetEV2AU();     
   this->am1DBondingParameterP = 0.0 * Parameters::GetInstance()->GetEV2AU();     
   this->am1DAlpha = 3.577756 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3CoreintegralS = -13.073321 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3CoreintegralP = 0.0;
   this->pm3OrbitalExponentS = 0.967807;      
   this->pm3OrbitalExponentP = 0.0;      
   this->pm3BondingParameterS = -5.626512 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3BondingParameterP = 0.0;
   this->pm3Alpha = 3.356386 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3DerivedParameterD[0] = 0.0;    
   this->pm3DerivedParameterD[1] = 0.0;    
   this->pm3DerivedParameterD[2] = 0.0;    
   this->pm3DerivedParameterRho[0] = 0.5/0.5436727936;
   this->pm3DerivedParameterRho[1] = 0.0;  
   this->pm3DerivedParameterRho[2] = 0.0;  
   this->pm3ParameterK[0] = 1.128750 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[1] =-1.060329 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[2] = 0.0;
   this->pm3ParameterK[3] = 0.0;
   this->pm3ParameterL[0] = 5.096282 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[1] = 6.003788 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[2] = 0.00;
   this->pm3ParameterL[3] = 0.00;
   this->pm3ParameterM[0] = 1.537465 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[1] = 1.570189 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[2] = 0.00;
   this->pm3ParameterM[3] = 0.00;
   this->pm3Gss = 14.794208 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp = 0.0;
   this->pm3Gsp = 0.0;
   this->pm3Gpp2 = 0.0;   
   this->pm3Hsp = 0.0;    
   this->pm3PddgCoreintegralS = -12.893272 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3PddgCoreintegralP = 0.0;
   this->pm3PddgOrbitalExponentS = 0.972786;      
   this->pm3PddgOrbitalExponentP = 0.0;      
   this->pm3PddgBondingParameterS = -6.152654 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3PddgBondingParameterP = 0.0;
   this->pm3PddgAlpha = 3.381686 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3PddgDerivedParameterD[0] = 0.0;    
   this->pm3PddgDerivedParameterD[1] = 0.0;    
   this->pm3PddgDerivedParameterD[2] = 0.0;    
   this->pm3PddgDerivedParameterRho[0] = 0.919616;
   this->pm3PddgDerivedParameterRho[1] = 0.0;  
   this->pm3PddgDerivedParameterRho[2] = 0.0;  
   this->pm3PddgParameterK[0] = 1.122244 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterK[1] =-1.069737 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterK[2] = 0.0;
   this->pm3PddgParameterK[3] = 0.0;
   this->pm3PddgParameterL[0] = 4.707790 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3PddgParameterL[1] = 5.857995 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3PddgParameterL[2] = 0.00;
   this->pm3PddgParameterL[3] = 0.00;
   this->pm3PddgParameterM[0] = 1.547099 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterM[1] = 1.567893 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterM[2] = 0.00;
   this->pm3PddgParameterM[3] = 0.00;
   this->pm3PddgParameterPa[0] = 0.057193 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterPa[1] =-0.034823 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterDa[0] = 0.663395 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterDa[1] = 1.081901 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3DCoreintegralS = -13.054076 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3DCoreintegralP = 0.0 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3DBondingParameterS = -5.628901 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3DBondingParameterP = 0.0 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3DAlpha = 3.417532 / Parameters::GetInstance()->GetAngstrom2AU();        
}
}
