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
#include"Natom.h"
using namespace std;
using namespace MolDS_base;
namespace MolDS_base_atoms{
Natom::Natom(int index) : Atom(index){
   this->SetAtomicParameters();
}

Natom::~Natom(){}

void Natom::SetAtomicParameters(){
   this->atomType = N;
   this->atomicMass = 14.00674*Parameters::GetInstance()->GetGMolin2AU();
   this->coreCharge = 5.0;
   this->numberValenceElectrons = 5;
   this->valenceShellType = lShell;
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   for(int i=0; i<this->valence.size();i++){
      this->realSphericalHarmonicsIndeces.push_back(new RealSphericalHarmonicsIndex(this->valence[i]));
   }
   this->vdWCoefficient = 1.11*Parameters::GetInstance()->GetJ2AU()
                              *pow(Parameters::GetInstance()->GetNm2AU(),6.0)
                              /Parameters::GetInstance()->GetAvogadro();
   this->vdWRadii = 1.550*Parameters::GetInstance()->GetAngstrom2AU();
   this->bondingParameter = -25.0*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuS = 19.316*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 7.275*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.0;
   this->effectiveNuclearChargeK = 6.7;
   this->effectiveNuclearChargeL = 3.90;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   this->indoG1 = 0.346029;
   this->indoF2 = 0.219055;
   this->indoF0CoefficientS = (this->coreCharge - 0.5);
   this->indoF0CoefficientP = (this->coreCharge - 0.5);
   this->indoG1CoefficientS = -1.0*(this->coreCharge - 1.5)/6.0;
   this->indoG1CoefficientP = -1.0/3.0;
   this->indoF2CoefficientS = 0.0;
   this->indoF2CoefficientP = -2.0*(this->coreCharge - 2.5)/25.0;
   this->zindoBondingParameterS = -26.0*Parameters::GetInstance()->GetEV2AU();
   this->zindoBondingParameterD = 0.0;
   this->zindoF0ss = 12.01 * Parameters::GetInstance()->GetEV2AU();                  
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                  
   this->zindoG1sp = 72255*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoF2pp = 52100*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoG2sd = 0.0;                 
   this->zindoG1pd = 0.0;                 
   this->zindoF2pd = 0.0;                 
   this->zindoG3pd = 0.0;                 
   this->zindoF2dd = 0.0;                 
   this->zindoF4dd = 0.0;                 
   this->zindoL = 2;
   this->zindoM = 3;
   this->zindoN = 0;
   this->zindoIonPotS = 25.69 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotP = 14.05 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotD = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralS = -71.932122 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = -57.172319 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoOrbitalExponentS = 2.255614;      
   this->mndoOrbitalExponentP = 2.255614;      
   this->mndoBondingParameterS = -20.495758 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = -20.495758 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha = 2.861342 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -202.581201 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 113.00 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoGss =  13.59 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp =  12.98 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp =  12.66 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 = 11.59 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =   3.14 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0.0;
   this->mndoDerivedParameterD[1] =   0.6399036683;
   this->mndoDerivedParameterD[2] =   0.5429762678;
   this->mndoDerivedParameterRho[0] = 0.5/0.4994193177;
   this->mndoDerivedParameterRho[1] = 0.5/0.7843433156;
   this->mndoDerivedParameterRho[2] = 0.5/0.8126295047;
   //this->mndoDerivedParameterD[0] =   0.0;
   //this->mndoDerivedParameterD[1] =   0.338616 * Parameters::GetInstance()->GetAngstrom2AU();    
   //this->mndoDerivedParameterD[2] =   0.287325 * Parameters::GetInstance()->GetAngstrom2AU();    
   //this->mndoDerivedParameterRho[0] = 0.529751 * Parameters::GetInstance()->GetAngstrom2AU();  
   //this->mndoDerivedParameterRho[1] = 0.337322 * Parameters::GetInstance()->GetAngstrom2AU();  
   //this->mndoDerivedParameterRho[2] = 0.324853 * Parameters::GetInstance()->GetAngstrom2AU();  
   this->am1CoreintegralS = -71.860000 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -57.167581 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 2.315410;      
   this->am1OrbitalExponentP = 2.157940;      
   this->am1BondingParameterS = -20.299110 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -18.238666 * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 2.947286 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 0.6433247425;    
   this->am1DerivedParameterD[2] = 0.5675527917;    
   this->am1DerivedParameterRho[0] = 0.5/0.4994193177;
   this->am1DerivedParameterRho[1] = 0.5/0.7820630445;  
   this->am1DerivedParameterRho[2] = 0.5/0.7883351388;  
   this->am1ParameterK[0] = 0.025251 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.028953 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] =-0.005806 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 5.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 2.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 1.50 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 2.10 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 2.40 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1DCoreintegralS = -71.997845 * Parameters::GetInstance()->GetEV2AU();         
   this->am1DCoreintegralP = -57.401718 * Parameters::GetInstance()->GetEV2AU();         
   this->am1DBondingParameterS = -20.092408 * Parameters::GetInstance()->GetEV2AU();     
   this->am1DBondingParameterP = -18.470679 * Parameters::GetInstance()->GetEV2AU();     
   this->am1DAlpha = 2.968737 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3CoreintegralS = -49.335672 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3CoreintegralP = -47.509736 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3OrbitalExponentS = 2.028094;      
   this->pm3OrbitalExponentP = 2.313728;      
   this->pm3BondingParameterS = -14.062521 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3BondingParameterP = -20.043848 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3Alpha = 2.830545 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3DerivedParameterD[0] = 0.0;    
   this->pm3DerivedParameterD[1] = 0.6577005762;    
   this->pm3DerivedParameterD[2] = 0.5293383109;    
   this->pm3DerivedParameterRho[0] = 0.5/0.4374893746;
   this->pm3DerivedParameterRho[1] = 0.5/0.5030877737;  
   this->pm3DerivedParameterRho[2] = 0.5/0.7364801616;  
   this->pm3ParameterK[0] = 1.501674 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[1] =-1.505772 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[2] = 0.0;
   this->pm3ParameterK[3] = 0.0;
   this->pm3ParameterL[0] = 5.901148 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[1] = 6.004658 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[2] = 0.00;
   this->pm3ParameterL[3] = 0.00;
   this->pm3ParameterM[0] = 1.710740 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[1] = 1.716149 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[2] = 0.00;
   this->pm3ParameterM[3] = 0.00;
   this->pm3Gss = 11.904787 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp = 11.754672 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gsp = 7.348565 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp2 = 10.807277 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Hsp = 1.136713 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3PddgCoreintegralS = -49.454546 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3PddgCoreintegralP = -47.757406 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgOrbitalExponentS = 2.035807;      
   this->pm3PddgOrbitalExponentP = 2.324327;      
   this->pm3PddgBondingParameterS = -14.117230 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3PddgBondingParameterP = -19.938509 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgAlpha = 2.849124 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3PddgDerivedParameterD[0] = 0.0;    
   this->pm3PddgDerivedParameterD[1] = 0.654855;    
   this->pm3PddgDerivedParameterD[2] = 0.526924;    
   this->pm3PddgDerivedParameterRho[0] = 1.142818;
   this->pm3PddgDerivedParameterRho[1] = 0.991235;  
   this->pm3PddgDerivedParameterRho[2] = 0.676704;  
   this->pm3PddgParameterK[0] = 1.513320 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterK[1] =-1.511892 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterK[2] = 0.0;
   this->pm3PddgParameterK[3] = 0.0;
   this->pm3PddgParameterL[0] = 5.904394 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3PddgParameterL[1] = 6.030014 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3PddgParameterL[2] = 0.00;
   this->pm3PddgParameterL[3] = 0.00;
   this->pm3PddgParameterM[0] = 1.728376 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterM[1] = 1.734108 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterM[2] = 0.00;
   this->pm3PddgParameterM[3] = 0.00;
   this->pm3PddgParameterPa[0] =-0.003160 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterPa[1] = 0.012501 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterDa[0] = 1.004172 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterDa[1] = 1.516336 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3DCoreintegralS = -49.348460 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3DCoreintegralP = -47.543768 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3DBondingParameterS = -14.068411 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3DBondingParameterP = -20.039292 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3DAlpha = 3.060404 / Parameters::GetInstance()->GetAngstrom2AU();        
}
}
