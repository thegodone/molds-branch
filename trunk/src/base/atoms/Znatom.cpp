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
#include"Znatom.h"
using namespace std;
using namespace MolDS_base;
namespace MolDS_base_atoms{
Znatom::Znatom(int index) : Atom(index){
   this->SetAtomicParameters();
}

Znatom::~Znatom(){}

void Znatom::SetAtomicParameters(){
   this->atomType = Zn;
   this->atomicMass = 65.38*Parameters::GetInstance()->GetGMolin2AU();
   this->coreCharge = 2.0;
   this->numberValenceElectrons = 2;
   if(Parameters::GetInstance()->GetCurrentTheory() == ZINDOS){
      this->coreCharge = 12.0;            // spd set for ZINDO/S
      this->numberValenceElectrons = 12;  // spd set for ZINDO/S
   }
   this->valenceShellType = nShell;
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   if(Parameters::GetInstance()->GetCurrentTheory() == ZINDOS){
      this->valence.push_back(dxy);
      this->valence.push_back(dyz);
      this->valence.push_back(dzz);
      this->valence.push_back(dzx);
      this->valence.push_back(dxxyy);
   }
   for(int i=0; i<this->valence.size();i++){
      this->realSphericalHarmonicsIndeces.push_back(new RealSphericalHarmonicsIndex(this->valence[i]));
   }
   /*
   this->vdWCoefficient = 1.65*Parameters::GetInstance()->GetJ2AU()
                              *pow(Parameters::GetInstance()->GetNm2AU(),6.0)
                              /Parameters::GetInstance()->GetAvogadro();
   this->vdWRadii = 1.610*Parameters::GetInstance()->GetAngstrom2AU();
   this->bondingParameter = -21.0*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuS = 14.051*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 5.572*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.0;
   */
   this->effectiveNuclearChargeK   = 29.70;
   this->effectiveNuclearChargeL   = 25.85;
   this->effectiveNuclearChargeMsp = 18.75;
   this->effectiveNuclearChargeMd  =  8.85;
   this->effectiveNuclearChargeNsp =  4.35;
   /*
   this->indoG1 = 0.267708;
   this->indoF2 = 0.17372;
   this->indoF0CoefficientS = (this->coreCharge - 0.5);
   this->indoF0CoefficientP = (this->coreCharge - 0.5);
   this->indoG1CoefficientS = -1.0*(this->coreCharge - 1.5)/6.0;
   this->indoG1CoefficientP = -1.0/3.0;
   this->indoF2CoefficientS = 0.0;
   this->indoF2CoefficientP = -2.0*(this->coreCharge - 2.5)/25.0;
   */
   this->zindoBondingParameterS = -10.0*Parameters::GetInstance()->GetEV2AU(); // orca 3.1
   this->zindoBondingParameterD = -34.0*Parameters::GetInstance()->GetEV2AU(); // orca 3.1
   this->zindoF0ss = 7.98  * Parameters::GetInstance()->GetEV2AU();   // ZHKW_1980
   this->zindoF0sd = 9.39  * Parameters::GetInstance()->GetEV2AU();   // ZHKW_1980
   this->zindoF0dd = 14.55 * Parameters::GetInstance()->GetEV2AU();   // ZHKW_1980
   this->zindoG1sp = 20400 * Parameters::GetInstance()->GetKayser2AU(); // BZ_1979
   this->zindoF2pp =  9500 * Parameters::GetInstance()->GetKayser2AU(); // BZ_1979
   this->zindoG2sd = 0.6199 * Parameters::GetInstance()->GetEV2AU();   // orca 3.1
   this->zindoG1pd = 0.8679 * Parameters::GetInstance()->GetEV2AU();   // orca 3.1
   this->zindoF2pd = 1.4878 * Parameters::GetInstance()->GetEV2AU();   // orca 3.1
   this->zindoG3pd = 0.9919 * Parameters::GetInstance()->GetEV2AU();   // orca 3.1
   this->zindoF2dd =11.4063 * Parameters::GetInstance()->GetEV2AU();   // orca 3.1
   this->zindoF4dd = 7.6249 * Parameters::GetInstance()->GetEV2AU();   // orca 3.1
   this->zindoL =  2;
   this->zindoM =  0;
   this->zindoN = 10;
   this->zindoIonPotS = 3.54  * Parameters::GetInstance()->GetEV2AU(); //  BZ_1979
   this->zindoIonPotP = 4.77  * Parameters::GetInstance()->GetEV2AU(); //  BZ_1979
   this->zindoIonPotD = 17.57 * Parameters::GetInstance()->GetEV2AU(); //  BZ_1979
   this->mndoCoreintegralS = -20.8397160 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = -19.6252240 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoOrbitalExponentS = 2.0473590;      
   this->mndoOrbitalExponentP = 1.4609460;      
   this->mndoBondingParameterS = -1.0000 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = -2.0000  * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha = 1.5064570 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -29.879432 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 31.17 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoGss = 11.800000 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp = 13.300000 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp = 11.182018 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2= 12.930520 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =  0.484606 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0.0;
   this->mndoDerivedParameterD[1] =   1.3037825507;
   this->mndoDerivedParameterD[2] =   1.4520183111;    
   this->mndoDerivedParameterRho[0] = 0.5/0.4336385540;  
   this->mndoDerivedParameterRho[1] = 0.5/0.2375857385;
   this->mndoDerivedParameterRho[2] = 0.5/0.2738809734;  
   this->am1CoreintegralS = -21.040008 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -17.655574 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 1.954299;      
   this->am1OrbitalExponentP = 1.372365;      
   this->am1BondingParameterS = -1.997429 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -4.758119 * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 1.484563 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2= this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 1.3581112717;    
   this->am1DerivedParameterD[2] = 1.5457406328;    
   this->am1DerivedParameterRho[0] = 0.5/0.4336385540;
   this->am1DerivedParameterRho[1] = 0.5/0.2317370207;  
   this->am1DerivedParameterRho[2] = 0.5/0.2621118818;  
   this->am1ParameterK[0] = 0.00;
   this->am1ParameterK[1] = 0.00;
   this->am1ParameterK[2] = 0.00;
   this->am1ParameterK[3] = 0.00;
   this->am1ParameterL[0] = 0.00;
   this->am1ParameterL[1] = 0.00;
   this->am1ParameterL[2] = 0.00;
   this->am1ParameterL[3] = 0.00;
   this->am1ParameterM[0] = 0.00;
   this->am1ParameterM[1] = 0.00;
   this->am1ParameterM[2] = 0.00;
   this->am1ParameterM[3] = 0.00;
   /*
   this->am1DCoreintegralS = -52.183798 * Parameters::GetInstance()->GetEV2AU();         
   this->am1DCoreintegralP = -39.368413 * Parameters::GetInstance()->GetEV2AU();         
   this->am1DBondingParameterS = -15.682341 * Parameters::GetInstance()->GetEV2AU();     
   this->am1DBondingParameterP = -7.804762 * Parameters::GetInstance()->GetEV2AU();     
   this->am1DAlpha = 2.625506 / Parameters::GetInstance()->GetAngstrom2AU();        
   */
   this->pm3CoreintegralS = -18.532198 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3CoreintegralP = -11.047409 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3OrbitalExponentS = 1.819989;      
   this->pm3OrbitalExponentP = 1.506922;      
   this->pm3BondingParameterS = -0.715578 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3BondingParameterP = -6.351864 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3Alpha = 1.350126 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3DerivedParameterD[0] = 0.0;    
   this->pm3DerivedParameterD[1] = 1.5005757621;    
   this->pm3DerivedParameterD[2] = 1.4077174157;    
   this->pm3DerivedParameterRho[0] = 0.5/0.3556275661;
   this->pm3DerivedParameterRho[1] = 0.5/0.2375632937;  
   this->pm3DerivedParameterRho[2] = 0.5/0.2661024312;  
   this->pm3ParameterK[0] = -0.111234 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[1] = -0.132370 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[2] = 0.0;
   this->pm3ParameterK[3] = 0.0;
   this->pm3ParameterL[0] = 6.001478 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[1] = 1.995839 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[2] = 0.00;
   this->pm3ParameterL[3] = 0.00;
   this->pm3ParameterM[0] = 1.516032 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[1] = 2.519642 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[2] = 0.00;
   this->pm3ParameterM[3] = 0.00;
   this->pm3Gss = 9.677196 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp = 4.980174 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gsp = 7.736204 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp2= 4.669656 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Hsp = 0.600413 * Parameters::GetInstance()->GetEV2AU();   
   /*
   this->pm3PddgCoreintegralS = -48.241241 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3PddgCoreintegralP = -36.461256 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgOrbitalExponentS = 1.567864;      
   this->pm3PddgOrbitalExponentP = 1.846659;      
   this->pm3PddgBondingParameterS = -11.952818 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3PddgBondingParameterP =  -9.922411 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgAlpha = 2.725772 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3PddgDerivedParameterD[0] = 0.0;    
   this->pm3PddgDerivedParameterD[1] = 0.831413;    
   this->pm3PddgDerivedParameterD[2] = 0.663222;    
   this->pm3PddgDerivedParameterRho[0] = 1.214657;
   this->pm3PddgDerivedParameterRho[1] = 0.848467;  
   this->pm3PddgDerivedParameterRho[2] = 0.652785;  
   this->pm3PddgParameterK[0] = 0.048906 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterK[1] = 0.047697 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterK[2] = 0.0;
   this->pm3PddgParameterK[3] = 0.0;
   this->pm3PddgParameterL[0] = 5.765340 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3PddgParameterL[1] = 5.973721 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3PddgParameterL[2] = 0.00;
   this->pm3PddgParameterL[3] = 0.00;
   this->pm3PddgParameterM[0] = 1.682232 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterM[1] = 0.894406 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterM[2] = 0.00;
   this->pm3PddgParameterM[3] = 0.00;
   this->pm3PddgParameterPa[0] =-0.000743 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterPa[1] = 0.000985 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterDa[0] = 0.836915 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterDa[1] = 1.585236 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3DCoreintegralS = -47.275431 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3DCoreintegralP = -36.268916 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3DBondingParameterS = -11.941466 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3DBondingParameterP = -9.819760 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3DAlpha = 2.721152 / Parameters::GetInstance()->GetAngstrom2AU();        
   */
}
}
