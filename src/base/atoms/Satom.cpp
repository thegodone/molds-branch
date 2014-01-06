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
#include"Satom.h"
using namespace std;
using namespace MolDS_base;
namespace MolDS_base_atoms{
Satom::Satom(int index) : Atom(index){
   this->SetAtomicParameters();
}

Satom::~Satom(){}

void Satom::SetAtomicParameters(){
   this->atomType = S;
   this->atomicMass = 32.066*Parameters::GetInstance()->GetGMolin2AU();
   this->coreCharge = 6.0;
   this->numberValenceElectrons = 6;
   this->valenceShellType = mShell;
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   if(Parameters::GetInstance()->GetCurrentTheory() == CNDO2 || 
      Parameters::GetInstance()->GetCurrentTheory() == INDO){
      this->valence.push_back(dxy);
      this->valence.push_back(dyz);
      this->valence.push_back(dzz);
      this->valence.push_back(dzx);
      this->valence.push_back(dxxyy);
   }
   for(int i=0; i<this->valence.size();i++){
      this->realSphericalHarmonicsIndeces.push_back(new RealSphericalHarmonicsIndex(this->valence[i]));
   }
   this->vdWCoefficient = 10.3*Parameters::GetInstance()->GetJ2AU()
                              *pow(Parameters::GetInstance()->GetNm2AU(),6.0)
                              /Parameters::GetInstance()->GetAvogadro();
   this->vdWRadii = 1.870*Parameters::GetInstance()->GetAngstrom2AU();
   this->bondingParameter = -18.150*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuS = 17.650*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 6.989*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.713*Parameters::GetInstance()->GetEV2AU();
   this->effectiveNuclearChargeK = 15.70;
   this->effectiveNuclearChargeL = 11.85;
   if(Parameters::GetInstance()->GetCurrentTheory() == ZINDOS){
      this->effectiveNuclearChargeMsp = 1.925*3.0;
      this->effectiveNuclearChargeMd = 1.731*3.0;
   }
   else{
      this->effectiveNuclearChargeMsp = 5.45;
      this->effectiveNuclearChargeMd = 5.45;
   }
   this->indoG1 = 0.267708;
   this->indoF2 = 0.17372;
   this->indoF0CoefficientS = 0.0;
   this->indoF0CoefficientP = 0.0;
   this->indoG1CoefficientS = 0.0;
   this->indoG1CoefficientP = 0.0;
   this->indoF2CoefficientS = 0.0;
   this->indoF2CoefficientP = 0.0;

   // ORCA parameter 2.8 set
   // see "ORCA 2.8"( http://www.thch.uni-bonn.de/tc/orca/ ).
   this->zindoBondingParameterS = -15.0*Parameters::GetInstance()->GetEV2AU();
   this->zindoBondingParameterD =   0.0*Parameters::GetInstance()->GetEV2AU();
   this->zindoF0ss = 10.09 * Parameters::GetInstance()->GetEV2AU(); 
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                 
   this->zindoG1sp = 3.0756 * Parameters::GetInstance()->GetEV2AU();                 
   this->zindoF2pp = 4.5377 * Parameters::GetInstance()->GetEV2AU();
   this->zindoG2sd = 0.0;
   this->zindoG1pd = 0.0;        
   this->zindoF2pd = 0.0;
   this->zindoG3pd = 0.0;
   this->zindoF2dd = 0.0;
   this->zindoF4dd = 0.0;
   // end (ORCA 2.8 parameter set)

   /*
   // Parameter set in [HKLWNZ_1982]
   // Take care that F0s are not included in this paper.
   // So, these parameters may be used as ones in [GD_1972]
   this->zindoBondingParameterS = -14.0*Parameters::GetInstance()->GetEV2AU();
   this->zindoBondingParameterD =   4.0*Parameters::GetInstance()->GetEV2AU();
   this->zindoF0ss = 8.96 * Parameters::GetInstance()->GetEV2AU(); // from [GD_1972]
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                 
   this->zindoG1sp = 3.10 * Parameters::GetInstance()->GetEV2AU();                 
   this->zindoF2pp = 4.57 * Parameters::GetInstance()->GetEV2AU();
   this->zindoG2sd = 3.25 * Parameters::GetInstance()->GetEV2AU();
   this->zindoG1pd = 4.31 * Parameters::GetInstance()->GetEV2AU();        
   this->zindoF2pd = 3.45 * Parameters::GetInstance()->GetEV2AU();
   this->zindoG3pd = 2.57 * Parameters::GetInstance()->GetEV2AU();
   this->zindoF2dd = 3.55 * Parameters::GetInstance()->GetEV2AU();
   this->zindoF4dd = 2.31 * Parameters::GetInstance()->GetEV2AU();
   // end(Parameter set in [HKLWNZ_1982])
   */

   /*
   // Parameter set in [BZ_1979]
   // Take care that F0s and bondingParameters are not included in this paper.
   // So, parameters for those in [HKLWNZ_1982] and [GD_1972] may be used.
   // Furthermore, this parameter set are not suitable for spectroscopy.
   this->zindoBondingParameterS = -14.0*Parameters::GetInstance()->GetEV2AU(); // from [HKLWNZ_1982]
   this->zindoBondingParameterD =   4.0*Parameters::GetInstance()->GetEV2AU(); // from [HKLWNZ_1982]
   this->zindoF0ss = 8.96 * Parameters::GetInstance()->GetEV2AU(); // from [GD_1972]
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                 
   this->zindoG1sp = 24807*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoF2pp = 36600*Parameters::GetInstance()->GetKayser2AU();                 
   this->zindoG2sd = 25972*Parameters::GetInstance()->GetKayser2AU();     
   this->zindoG1pd = 34486*Parameters::GetInstance()->GetKayser2AU();        
   this->zindoF2pd = 29173*Parameters::GetInstance()->GetKayser2AU();           
   this->zindoG3pd = 20587*Parameters::GetInstance()->GetKayser2AU();           
   this->zindoF2dd = 28411*Parameters::GetInstance()->GetKayser2AU();           
   this->zindoF4dd = 18529*Parameters::GetInstance()->GetKayser2AU();           
   // end(Parameter set in [BZ_1979]) 
   */

   this->zindoL = 2;
   this->zindoM = 4;
   this->zindoN = 0;
   this->zindoIonPotS = 21.11 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotP = 12.39 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotD = 4.11 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralS = -72.242281 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoCoreintegralP = -56.973207 * Parameters::GetInstance()->GetEV2AU();         
   this->mndoOrbitalExponentS = 2.312962;      
   this->mndoOrbitalExponentP = 2.009146;      
   this->mndoBondingParameterS = -10.761670 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoBondingParameterP = -10.108433 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha = 2.478026 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -226.01239 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom = 66.40 * Parameters::GetInstance()->GetKcalMolin2AU();
   this->mndoGss =  12.88 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp =   9.90 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp =  11.26 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 =  8.83 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =   2.26 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0.0;
   this->mndoDerivedParameterD[1] =   0.9189935137;
   this->mndoDerivedParameterD[2] =   0.8328513971;
   this->mndoDerivedParameterRho[0] = 0.5/0.4733275064;
   this->mndoDerivedParameterRho[1] = 0.5/0.5544352823;
   this->mndoDerivedParameterRho[2] = 0.5/0.5585137839;
   this->am1CoreintegralS = -56.694056 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -48.717049 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 2.366515;      
   this->am1OrbitalExponentP = 1.667263;      
   this->am1BondingParameterS = -3.920566 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -7.905278 * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 2.461648 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss =  11.786329 * Parameters::GetInstance()->GetEV2AU();   
   this->am1Gpp =  10.039308 * Parameters::GetInstance()->GetEV2AU();   
   this->am1Gsp =   8.663127 * Parameters::GetInstance()->GetEV2AU();   
   this->am1Gpp2 =  7.781688 * Parameters::GetInstance()->GetEV2AU();  
   this->am1Hsp =   2.532137 * Parameters::GetInstance()->GetEV2AU();   
   this->am1DerivedParameterD[0] = 0.0;    
   this->am1DerivedParameterD[1] = 0.9004264562;    
   this->am1DerivedParameterD[2] = 1.0036329320;    
   this->am1DerivedParameterRho[0] = 0.5/0.4331361580;
   this->am1DerivedParameterRho[1] = 0.5/0.5906953135;  
   this->am1DerivedParameterRho[2] = 0.5/0.6454793983;  
   this->am1ParameterK[0] =-0.509195 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] =-0.011863 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] = 0.012334 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 4.593691 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 5.865731 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 13.557336 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 0.770665 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 1.503313 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 2.009173 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1DCoreintegralS = -57.235044 * Parameters::GetInstance()->GetEV2AU();         
   this->am1DCoreintegralP = -48.307513 * Parameters::GetInstance()->GetEV2AU();         
   this->am1DBondingParameterS = -3.311308 * Parameters::GetInstance()->GetEV2AU();     
   this->am1DBondingParameterP = -7.256468 * Parameters::GetInstance()->GetEV2AU();     
   this->am1DAlpha = 2.309315 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3CoreintegralS = -49.895371 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3CoreintegralP = -44.392583 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3OrbitalExponentS = 1.891185;      
   this->pm3OrbitalExponentP = 1.658972;      
   this->pm3BondingParameterS = -8.827465 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3BondingParameterP = -8.091415 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3Alpha = 2.269706 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3DerivedParameterD[0] = 0.0;    
   this->pm3DerivedParameterD[1] = 1.1214312500;    
   this->pm3DerivedParameterD[2] = 1.0086487614;    
   this->pm3DerivedParameterRho[0] = 0.5/0.3294428165;
   this->pm3DerivedParameterRho[1] = 0.5/0.6678906502;  
   this->pm3DerivedParameterRho[2] = 0.5/0.6137333700;  
   this->pm3ParameterK[0] = -0.399191 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[1] = -0.054899 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[2] = 0.0;
   this->pm3ParameterK[3] = 0.0;
   this->pm3ParameterL[0] = 6.000669 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[1] = 6.001845 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[2] = 0.00;
   this->pm3ParameterL[3] = 0.00;
   this->pm3ParameterM[0] = 0.962123 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[1] = 1.579944 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[2] = 0.00;
   this->pm3ParameterM[3] = 0.00;
   this->pm3Gss = 8.964667 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp = 9.968164 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gsp = 6.785936 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Gpp2 = 7.970247 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3Hsp = 4.041836 * Parameters::GetInstance()->GetEV2AU();   
   this->pm3PddgCoreintegralS = -43.906366 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3PddgCoreintegralP = -43.461348 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3PddgOrbitalExponentS = 1.012002;      
   this->pm3PddgOrbitalExponentP = 1.876999;      
   this->pm3PddgBondingParameterS = -2.953912 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3PddgBondingParameterP = -8.507779 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3PddgAlpha = 2.539751 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3PddgDerivedParameterD[0] = 0.0;    
   this->pm3PddgDerivedParameterD[1] = 1.006989;    
   this->pm3PddgDerivedParameterD[2] = 0.891487;    
   this->pm3PddgDerivedParameterRho[0] = 1.517625;
   this->pm3PddgDerivedParameterRho[1] = 0.711672;  
   this->pm3PddgDerivedParameterRho[2] = 0.754336;  
   this->pm3PddgParameterK[0] =-0.330692 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterK[1] = 0.024171 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterK[2] = 0.0;
   this->pm3PddgParameterK[3] = 0.0;
   this->pm3PddgParameterL[0] = 6.000000 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3PddgParameterL[1] = 6.000000 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3PddgParameterL[2] = 0.00;
   this->pm3PddgParameterL[3] = 0.00;
   this->pm3PddgParameterM[0] = 0.823837 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterM[1] = 2.017756 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterM[2] = 0.00;
   this->pm3PddgParameterM[3] = 0.00;
   this->pm3PddgParameterPa[0] = 0.120434 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterPa[1] =-0.002663 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterDa[0] = 0.672870 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterDa[1] = 2.032340 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3DCoreintegralS = -50.249536 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3DCoreintegralP = -43.968965 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3DBondingParameterS = -8.397415 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3DBondingParameterP = -7.594232 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3DAlpha = 2.234331 / Parameters::GetInstance()->GetAngstrom2AU();        
}
}
