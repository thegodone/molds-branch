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
#include"Fatom.h"
using namespace std;
using namespace MolDS_base;
namespace MolDS_base_atoms{
Fatom::Fatom(int index) : Atom(index){
   this->SetAtomicParameters();
}

Fatom::~Fatom(){}

void Fatom::SetAtomicParameters(){
   this->atomType = F;
   this->atomicMass = 18.9984032*Parameters::GetInstance()->GetGMolin2AU();
   this->coreCharge = 7.0;
   this->numberValenceElectrons = 7;
   this->valenceShellType = lShell;
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   for(int i=0; i<this->valence.size();i++){
      this->realSphericalHarmonicsIndeces.push_back(new RealSphericalHarmonicsIndex(this->valence[i]));
   }
   this->vdWCoefficient = 0.57*Parameters::GetInstance()->GetJ2AU()
                              *pow(Parameters::GetInstance()->GetNm2AU(),6.0)
                              /Parameters::GetInstance()->GetAvogadro();
   this->vdWRadii = 1.430*Parameters::GetInstance()->GetAngstrom2AU();
   this->bondingParameter = -39.0*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuS = 32.272*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 11.080*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.0;
   this->effectiveNuclearChargeK = 8.7;
   this->effectiveNuclearChargeL = 5.20;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   this->indoG1 = 0.532305;
   this->indoF2 = 0.31580;
   this->indoF0CoefficientS = (this->coreCharge - 0.5);
   this->indoF0CoefficientP = (this->coreCharge - 0.5);
   this->indoG1CoefficientS = -1.0*(this->coreCharge - 1.5)/6.0;
   this->indoG1CoefficientP = -1.0/3.0;
   this->indoF2CoefficientS = 0.0;
   this->indoF2CoefficientP = -2.0*(this->coreCharge - 2.5)/25.0;
   this->zindoBondingParameterS = -44.0*Parameters::GetInstance()->GetEV2AU(); //from orca3.0.1
   this->zindoBondingParameterD = 0.0;
   this->zindoF0ss = 14.00 * Parameters::GetInstance()->GetEV2AU(); // from orca3.0.1
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                  
   this->zindoG1sp = 116828*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoF2pp =  69310*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoG2sd = 0.0;                 
   this->zindoG1pd = 0.0;                 
   this->zindoF2pd = 0.0;                 
   this->zindoG3pd = 0.0;                 
   this->zindoF2dd = 0.0;                 
   this->zindoF4dd = 0.0;                 
   this->zindoL = 2;
   this->zindoM = 5;
   this->zindoN = 0;
   this->zindoIonPotS = 39.39 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotP = 20.86 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotD = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralS = -131.071548 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = -105.782137 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoOrbitalExponentS = 2.848487;      
   this->mndoOrbitalExponentP = 2.848487;      
   this->mndoBondingParameterS = -48.290460 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = -36.508540 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha =  3.419661 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -476.683781 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 18.86 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoGss =  16.92 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp =  16.71 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp =  17.25 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 = 14.91 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =   4.83 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0.0;
   this->mndoDerivedParameterD[1] =   0.5067166088;
   this->mndoDerivedParameterD[2] =   0.4299633003;    
   this->mndoDerivedParameterRho[0] = 0.5/0.6217935876;  
   this->mndoDerivedParameterRho[1] = 0.5/1.0850000098;  
   this->mndoDerivedParameterRho[2] = 0.5/1.0343451703;  
   this->am1CoreintegralS = -136.105579 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -104.889885 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 3.770082;      
   this->am1OrbitalExponentP = 2.494670;      
   this->am1BondingParameterS = -69.590277 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -27.922360 * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 5.517800/ Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 0.4145203025;    
   this->am1DerivedParameterD[2] = 0.4909446425;    
   this->am1DerivedParameterRho[0] = 0.5/0.6217935876;
   this->am1DerivedParameterRho[1] = 0.5/1.2088469198;  
   this->am1DerivedParameterRho[2] = 0.5/0.9449175360;  
   this->am1ParameterK[0] = 0.242079 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.003607 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 4.80 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 4.60 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 0.930 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 1.660 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1DCoreintegralS = this->am1CoreintegralS;
   this->am1DCoreintegralP = this->am1CoreintegralP;
   this->am1DBondingParameterS = this->am1BondingParameterS;
   this->am1DBondingParameterP = this->am1BondingParameterP;
   this->am1DAlpha = this->am1DAlpha;
   this->pm3CoreintegralS = -110.435303 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3CoreintegralP = -105.685047 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3OrbitalExponentS = 4.708555;      
   this->pm3OrbitalExponentP = 2.491178;      
   this->pm3BondingParameterS = -48.405939 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3BondingParameterP = -27.744660 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3Alpha = 3.358921 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3DerivedParameterD[0] = 0.0;    
   this->pm3DerivedParameterD[1] = 0.3125302275;    
   this->pm3DerivedParameterD[2] = 0.4916328225;    
   this->pm3DerivedParameterRho[0] = 0.5/0.3857423305;
   this->pm3DerivedParameterRho[1] = 0.5/0.6768359077;  
   this->pm3DerivedParameterRho[2] = 0.5/0.6119953427;  
   this->pm3ParameterK[0] = -0.012166 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[1] = -0.002852 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[2] = 0.0;
   this->pm3ParameterK[3] = 0.0;
   this->pm3ParameterL[0] = 6.023574 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[1] = 6.003717 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[2] = 0.00;
   this->pm3ParameterL[3] = 0.00;
   this->pm3ParameterM[0] = 1.856859 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[1] = 2.636158 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[2] = 0.00;
   this->pm3ParameterM[3] = 0.00;
   this->pm3Gss = 10.496667 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp = 14.817256 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gsp = 16.073689 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp2= 14.418393 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Hsp =  0.727763 * Parameters::GetInstance()->GetEV2AU();   
   //this->pm3PddgCoreintegralS = -48.241241 * Parameters::GetInstance()->GetEV2AU();         
   //this->pm3PddgCoreintegralP = -36.461256 * Parameters::GetInstance()->GetEV2AU();
   //this->pm3PddgOrbitalExponentS = 1.567864;      
   //this->pm3PddgOrbitalExponentP = 1.846659;      
   //this->pm3PddgBondingParameterS = -11.952818 * Parameters::GetInstance()->GetEV2AU();     
   //this->pm3PddgBondingParameterP =  -9.922411 * Parameters::GetInstance()->GetEV2AU();
   //this->pm3PddgAlpha = 2.725772 / Parameters::GetInstance()->GetAngstrom2AU();        
   //this->pm3PddgDerivedParameterD[0] = 0.0;    
   //this->pm3PddgDerivedParameterD[1] = 0.831413;    
   //this->pm3PddgDerivedParameterD[2] = 0.663222;    
   //this->pm3PddgDerivedParameterRho[0] = 1.214657;
   //this->pm3PddgDerivedParameterRho[1] = 0.848467;  
   //this->pm3PddgDerivedParameterRho[2] = 0.652785;  
   //this->pm3PddgParameterK[0] = 0.048906 * Parameters::GetInstance()->GetEV2AU();
   //this->pm3PddgParameterK[1] = 0.047697 * Parameters::GetInstance()->GetEV2AU();
   //this->pm3PddgParameterK[2] = 0.0;
   //this->pm3PddgParameterK[3] = 0.0;
   //this->pm3PddgParameterL[0] = 5.765340 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   //this->pm3PddgParameterL[1] = 5.973721 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   //this->pm3PddgParameterL[2] = 0.00;
   //this->pm3PddgParameterL[3] = 0.00;
   //this->pm3PddgParameterM[0] = 1.682232 * Parameters::GetInstance()->GetAngstrom2AU();
   //this->pm3PddgParameterM[1] = 0.894406 * Parameters::GetInstance()->GetAngstrom2AU();
   //this->pm3PddgParameterM[2] = 0.00;
   //this->pm3PddgParameterM[3] = 0.00;
   //this->pm3PddgParameterPa[0] =-0.000743 * Parameters::GetInstance()->GetEV2AU();
   //this->pm3PddgParameterPa[1] = 0.000985 * Parameters::GetInstance()->GetEV2AU();
   //this->pm3PddgParameterDa[0] = 0.836915 * Parameters::GetInstance()->GetAngstrom2AU();
   //this->pm3PddgParameterDa[1] = 1.585236 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3DCoreintegralS = this->pm3CoreintegralS;
   this->pm3DCoreintegralP = this->pm3CoreintegralP;
   this->pm3DBondingParameterS = this->pm3BondingParameterS;
   this->pm3DBondingParameterP = this->pm3BondingParameterP;
   this->pm3DAlpha = this->pm3Alpha;
}
}
