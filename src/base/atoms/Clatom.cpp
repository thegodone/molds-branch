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
#include"Clatom.h"
using namespace std;
using namespace MolDS_base;
namespace MolDS_base_atoms{
Clatom::Clatom(int index) : Atom(index){
   this->SetAtomicParameters();
}

Clatom::~Clatom(){}

void Clatom::SetAtomicParameters(){
   this->atomType = Cl;
   this->atomicMass = 35.453*Parameters::GetInstance()->GetGMolin2AU();
   this->coreCharge = 7.0;
   this->numberValenceElectrons = 7;
   this->valenceShellType = mShell;
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   if(Parameters::GetInstance()->GetCurrentTheory() == CNDO2){
      this->valence.push_back(dxy);
      this->valence.push_back(dyz);
      this->valence.push_back(dzz);
      this->valence.push_back(dzx);
      this->valence.push_back(dxxyy);
   }
   for(int i=0; i<this->valence.size();i++){
      this->realSphericalHarmonicsIndeces.push_back(new RealSphericalHarmonicsIndex(this->valence[i]));
   }
   this->vdWCoefficient = 8.00*Parameters::GetInstance()->GetJ2AU()
                              *pow(Parameters::GetInstance()->GetNm2AU(),6.0)
                              /Parameters::GetInstance()->GetAvogadro();
   this->vdWRadii = 1.820*Parameters::GetInstance()->GetAngstrom2AU();
   this->bondingParameter = -22.330*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuS = 21.591*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 8.708*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.977*Parameters::GetInstance()->GetEV2AU();
   this->effectiveNuclearChargeK = 16.70;
   this->effectiveNuclearChargeL = 12.85;
   if(Parameters::GetInstance()->GetCurrentTheory() == ZINDOS){
      this->effectiveNuclearChargeMsp = 2.130*3.0; // from orca 3.0.1
      this->effectiveNuclearChargeMd  = 0.0;       // not used
   }
   else{
      this->effectiveNuclearChargeMsp = 6.10;
      this->effectiveNuclearChargeMd  = 6.10;
   }
   //this->indoG1 = 0.0;
   //this->indoF2 = 0.0;
   //this->indoF0CoefficientS = 0.0;
   //this->indoF0CoefficientP = 0.0;
   //this->indoG1CoefficientS = 0.0;
   //this->indoG1CoefficientP = 0.0;
   //this->indoF2CoefficientS = 0.0;
   //this->indoF2CoefficientP = 0.0;

   // ORCA parameter 3.0.1 set
   // see "ORCA 2.8"( http://www.thch.uni-bonn.de/tc/orca/ ).
   this->zindoBondingParameterS = -11.0*Parameters::GetInstance()->GetEV2AU();
   this->zindoBondingParameterD =   0.0*Parameters::GetInstance()->GetEV2AU();
   this->zindoF0ss = 11.25 * Parameters::GetInstance()->GetEV2AU(); 
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                 
   this->zindoG1sp = 8.8027 * Parameters::GetInstance()->GetEV2AU();                 
   this->zindoF2pp = 6.4470 * Parameters::GetInstance()->GetEV2AU();
   this->zindoG2sd = 0.0;
   this->zindoG1pd = 0.0;        
   this->zindoF2pd = 0.0;
   this->zindoG3pd = 0.0;
   this->zindoF2dd = 0.0;
   this->zindoF4dd = 0.0;
   // end (ORCA 2.8 parameter set)

   this->zindoL = 2;
   this->zindoM = 5;
   this->zindoN = 0;
   this->zindoIonPotS = 25.23 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotP = 15.03 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotD =  6.00 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralS = -100.227166 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = -77.378667 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoOrbitalExponentS = 3.784645;      
   this->mndoOrbitalExponentP = 2.036263;      
   this->mndoBondingParameterS = -14.262320 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = -14.26320 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha = 2.542201 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -353.137667 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 28.99 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoGss =  15.03 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp =  11.30 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp =  13.16 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 =  9.97 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =   2.42 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0.0;
   this->mndoDerivedParameterD[1] =   0.4986870220;
   this->mndoDerivedParameterD[2] =   0.8217602800;
   this->mndoDerivedParameterRho[0] = 0.5/0.5523379209;
   this->mndoDerivedParameterRho[1] = 0.5/0.8061021276;
   this->mndoDerivedParameterRho[2] = 0.5/0.6053315152;
   this->am1CoreintegralS = -111.613948 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP =  -76.640107 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 3.631376;      
   this->am1OrbitalExponentP = 2.076799;      
   this->am1BondingParameterS = -24.594670 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -14.637216 * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 2.919368 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss  = this->mndoGss;
   this->am1Gpp  = this->mndoGpp;
   this->am1Gsp  = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp  = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 0.5406286370;    
   this->am1DerivedParameterD[2] = 0.8057207525;    
   this->am1DerivedParameterRho[0] = 0.5/0.5523379209;
   this->am1DerivedParameterRho[1] = 0.5/0.7693007940;  
   this->am1DerivedParameterRho[2] = 0.5/0.6133247965;  
   this->am1ParameterK[0] = 0.094243 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.027168 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 4.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 4.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 1.30 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 2.10 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1DCoreintegralS = this->am1CoreintegralS;
   this->am1DCoreintegralP = this->am1CoreintegralP;
   this->am1DBondingParameterS = this->am1BondingParameterS;
   this->am1DBondingParameterP = this->am1BondingParameterP;
   this->am1DAlpha = this->am1DAlpha;
   this->pm3CoreintegralS = -100.626747 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3CoreintegralP =  -53.614396 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3OrbitalExponentS = 2.246210;      
   this->pm3OrbitalExponentP = 2.151010;      
   this->pm3BondingParameterS = -27.528560 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3BondingParameterP = -11.593922 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3Alpha = 2.517296 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3DerivedParameterD[0] = 0.0;    
   this->pm3DerivedParameterD[1] = 0.9175855709;    
   this->pm3DerivedParameterD[2] = 0.7779229539;    
   this->pm3DerivedParameterRho[0] = 0.5/0.5884843036;
   this->pm3DerivedParameterRho[1] = 0.5/0.6814322305;  
   this->pm3DerivedParameterRho[2] = 0.5/0.2074800643;  
   this->pm3ParameterK[0] = -0.171591 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[1] = -0.013458 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[2] = 0.0;
   this->pm3ParameterK[3] = 0.0;
   this->pm3ParameterL[0] = 6.000802 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[1] = 1.966618 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[2] = 0.00;
   this->pm3ParameterL[3] = 0.00;
   this->pm3ParameterM[0] = 1.087502 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[1] = 2.292891 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[2] = 0.00;
   this->pm3ParameterM[3] = 0.00;
   this->pm3Gss  = 16.013601 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp  =  7.522215 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gsp  =  8.048115 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp2 =  7.504154 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Hsp  =  3.481153 * Parameters::GetInstance()->GetEV2AU();   
   //this->pm3PddgCoreintegralS = -43.906366 * Parameters::GetInstance()->GetEV2AU();         
   //this->pm3PddgCoreintegralP = -43.461348 * Parameters::GetInstance()->GetEV2AU();         
   //this->pm3PddgOrbitalExponentS = 1.012002;      
   //this->pm3PddgOrbitalExponentP = 1.876999;      
   //this->pm3PddgBondingParameterS = -2.953912 * Parameters::GetInstance()->GetEV2AU();     
   //this->pm3PddgBondingParameterP = -8.507779 * Parameters::GetInstance()->GetEV2AU();     
   //this->pm3PddgAlpha = 2.539751 / Parameters::GetInstance()->GetAngstrom2AU();        
   //this->pm3PddgDerivedParameterD[0] = 0.0;    
   //this->pm3PddgDerivedParameterD[1] = 1.006989;    
   //this->pm3PddgDerivedParameterD[2] = 0.891487;    
   //this->pm3PddgDerivedParameterRho[0] = 1.517625;
   //this->pm3PddgDerivedParameterRho[1] = 0.711672;  
   //this->pm3PddgDerivedParameterRho[2] = 0.754336;  
   //this->pm3PddgParameterK[0] =-0.330692 * Parameters::GetInstance()->GetEV2AU();
   //this->pm3PddgParameterK[1] = 0.024171 * Parameters::GetInstance()->GetEV2AU();
   //this->pm3PddgParameterK[2] = 0.0;
   //this->pm3PddgParameterK[3] = 0.0;
   //this->pm3PddgParameterL[0] = 6.000000 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   //this->pm3PddgParameterL[1] = 6.000000 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   //this->pm3PddgParameterL[2] = 0.00;
   //this->pm3PddgParameterL[3] = 0.00;
   //this->pm3PddgParameterM[0] = 0.823837 * Parameters::GetInstance()->GetAngstrom2AU();
   //this->pm3PddgParameterM[1] = 2.017756 * Parameters::GetInstance()->GetAngstrom2AU();
   //this->pm3PddgParameterM[2] = 0.00;
   //this->pm3PddgParameterM[3] = 0.00;
   //this->pm3PddgParameterPa[0] = 0.120434 * Parameters::GetInstance()->GetEV2AU();
   //this->pm3PddgParameterPa[1] =-0.002663 * Parameters::GetInstance()->GetEV2AU();
   //this->pm3PddgParameterDa[0] = 0.672870 * Parameters::GetInstance()->GetAngstrom2AU();
   //this->pm3PddgParameterDa[1] = 2.032340 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3DCoreintegralS = this->pm3CoreintegralS;
   this->pm3DCoreintegralP = this->pm3CoreintegralP;
   this->pm3DBondingParameterS = this->pm3BondingParameterS;
   this->pm3DBondingParameterP = this->pm3BondingParameterP;
   this->pm3DAlpha = this->pm3Alpha;
}
}
