//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   // 
// 2016 Guillaume GODIN                                                   // 
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
#include"Patom.h"
using namespace std;
using namespace MolDS_base;
namespace MolDS_base_atoms{
Patom::Patom(int index) : Atom(index){
   this->SetAtomicParameters();
}

Patom::~Patom(){}

void Patom::SetAtomicParameters(){
   this->atomType = P;
   this->atomicMass = 30.9737620*Parameters::GetInstance()->GetGMolin2AU();
   this->coreCharge = 5.0;
   this->numberValenceElectrons = 5;
   this->valenceShellType = lShell;
   this->valenceShellType = lShell;
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   // like for Zn
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

   // values from Grimme 2006! computed using those parameters alpha = 20 & s6 = 1.11!
   this->vdWCoefficient = 7.84*Parameters::GetInstance()->GetJ2AU()
                              *pow(Parameters::GetInstance()->GetNm2AU(),6.0)
                              /Parameters::GetInstance()->GetAvogadro();
   this->vdWRadii = 1.705*Parameters::GetInstance()->GetAngstrom2AU();
   // found in pople paper IV table III
   this->bondingParameter = -15.070*Parameters::GetInstance()->GetEV2AU();
   ///  found in pople paper IV table II
   this->imuAmuS = 14.033*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 5.464*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.5;
   // from http://scientificsentence.net/Equations/Quantum/index.php?key=yes&Integer=slater
   this->effectiveNuclearChargeK = 14.70;
   this->effectiveNuclearChargeL = 10.85;
   this->effectiveNuclearChargeMsp = 4.8;
   this->effectiveNuclearChargeMd = 0.0;

// ORCA parameter 3.0.1 set
   // see "ORCA 2.8"( http://www.thch.uni-bonn.de/tc/orca/ ).
   this->zindoBondingParameterS = -15.0*Parameters::GetInstance()->GetEV2AU();
   this->zindoBondingParameterD =   0.0*Parameters::GetInstance()->GetEV2AU();
   this->zindoF0ss = 8.86 * Parameters::GetInstance()->GetEV2AU(); 
   this->zindoF0sd = 0.0;                   
   this->zindoF0dd = 0.0;                 
   this->zindoG1sp = 1.0556 * Parameters::GetInstance()->GetEV2AU();                 
   this->zindoF2pp = 2.9477 * Parameters::GetInstance()->GetEV2AU();
   this->zindoG2sd = 0.0;
   this->zindoG1pd = 0.0;        
   this->zindoF2pd = 0.0;
   this->zindoG3pd = 0.0;
   this->zindoF2dd = 0.0;
   this->zindoF4dd = 0.0;
   // end (ORCA 2.8 parameter set)

   this->zindoL = 2;
   this->zindoM = 3;
   this->zindoN = 0;
   this->zindoIonPotS = 1 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotP = 1 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotD = 1 * Parameters::GetInstance()->GetEV2AU();

   this->mndoHeatsFormAtom =  75.57 * Parameters::GetInstance()->GetKcalMolin2AU(); //Handbook of Bond Dissociation Energies in Organic Compounds
   // from MOPAC 7 ok
   this->mndoGss = 11.5600050 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp = 7.8775890 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp = 5.2374490 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 = 7.3076480 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =  0.7792380 * Parameters::GetInstance()->GetEV2AU();   
   // END OK MOPAC
   // what about this ? in ORCA it's look like
   this->effectiveNuclearChargeK = 8.7;
   this->effectiveNuclearChargeL = 5.20;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   // END what about this ? in ORCA it's look like
   // OK came from MOPAC 7 & ORIGINAL PUBLICATION
   this->am1CoreintegralS = -42.029863 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -34.030709 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 1.9812800;
   this->am1OrbitalExponentP = 1.8751500;
   this->am1BondingParameterS = -6.3537640 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -6.5907090 * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 2.455322/ Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;
   this->am1DerivedParameterD[1] = 1.0452021761;
   this->am1DerivedParameterD[2] = 0.8923659724;
   this->am1DerivedParameterRho[0] = 0.5/0.4248189705;
   this->am1DerivedParameterRho[1] = 0.5/0.3275242214;
   this->am1DerivedParameterRho[2] = 0;
   this->am1ParameterK[0] = -0.031827 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.0184700 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] = 0.0332900 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 6.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 7.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 9.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 1.4743230 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 1.7793540 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 3.0065760 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1DCoreintegralS = this->am1CoreintegralS;
   this->am1DCoreintegralP = this->am1CoreintegralP;
   this->am1DBondingParameterS = this->am1BondingParameterS;
   this->am1DBondingParameterP = this->am1BondingParameterP;
   this->am1DAlpha = this->am1DAlpha;
   this->rm1CoreintegralS = -41.81533184 * Parameters::GetInstance()->GetEV2AU();    
   this->rm1CoreintegralP = -34.38342529 * Parameters::GetInstance()->GetEV2AU();    
   this->rm1OrbitalExponentS = 2.12240118;
   this->rm1OrbitalExponentP = 1.74327954;
   this->rm1BondingParameterS = -6.13514969 * Parameters::GetInstance()->GetEV2AU();
   this->rm1BondingParameterP = -5.94442127 * Parameters::GetInstance()->GetEV2AU();
   this->rm1Alpha = 1.90993294 / Parameters::GetInstance()->GetAngstrom2AU();    
   this->rm1Gss = 11.08059265 * Parameters::GetInstance()->GetEV2AU();
   this->rm1Gpp = 5.68339201 * Parameters::GetInstance()->GetEV2AU();
   this->rm1Gsp = 7.60417563 * Parameters::GetInstance()->GetEV2AU();
   this->rm1Gpp2 =7.40265182 * Parameters::GetInstance()->GetEV2AU();
   this->rm1Hsp = 1.16181792 * Parameters::GetInstance()->GetEV2AU();
   this->rm1ParameterK[0] = -0.41063467 * Parameters::GetInstance()->GetEV2AU();
   this->rm1ParameterK[1] = -0.16299288 * Parameters::GetInstance()->GetEV2AU();
   this->rm1ParameterK[2] = -0.04887125 * Parameters::GetInstance()->GetEV2AU();
   this->rm1ParameterK[3] = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->rm1ParameterL[0] = 6.08752832 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->rm1ParameterL[1] = 7.09472602 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->rm1ParameterL[2] = 8.99979308 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->rm1ParameterL[3] = 0.0 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->rm1ParameterM[0] = 1.31650261 * Parameters::GetInstance()->GetAngstrom2AU();
   this->rm1ParameterM[1] = 1.90721319 * Parameters::GetInstance()->GetAngstrom2AU();
   this->rm1ParameterM[2] = 2.6585778  * Parameters::GetInstance()->GetAngstrom2AU();
   this->rm1ParameterM[3] = 0.0 * Parameters::GetInstance()->GetAngstrom2AU();
   this->rm1DerivedParameterD[0] = 0.0;
   this->rm1DerivedParameterD[1] = 1.0106954774;
   this->rm1DerivedParameterD[2] =  0.9598690369;
   this->rm1DerivedParameterRho[0] = 0.5 / 0.4072010317;
   this->rm1DerivedParameterRho[1] = 0.5 / 0.3932126744;
   this->rm1DerivedParameterRho[2] = 0.0;
   this->pm3CoreintegralS = -40.413096 * Parameters::GetInstance()->GetEV2AU();
   this->pm3CoreintegralP = -29.593052 * Parameters::GetInstance()->GetEV2AU();
   this->pm3OrbitalExponentS = 2.017563;
   this->pm3OrbitalExponentP = 1.504732;
   this->pm3BondingParameterS = -12.615879 * Parameters::GetInstance()->GetEV2AU();
   this->pm3BondingParameterP = -4.160040 * Parameters::GetInstance()->GetEV2AU();
   this->pm3Alpha = 1.940534 / Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3Gss = 7.8016150 * Parameters::GetInstance()->GetEV2AU();
   this->pm3Gpp = 6.6184780 * Parameters::GetInstance()->GetEV2AU();
   this->pm3Gsp = 5.1869490 * Parameters::GetInstance()->GetEV2AU();
   this->pm3Gpp2 =6.0620020 * Parameters::GetInstance()->GetEV2AU();
   this->pm3Hsp = 1.5428090 * Parameters::GetInstance()->GetEV2AU();
   this->pm3DerivedParameterD[0] = 0.0;
   this->pm3DerivedParameterD[1] = 1.0644946619;
   this->pm3DerivedParameterD[2] = 1.1120385910;    
   this->pm3DerivedParameterRho[0] = 0.5/0.2867017837;
   this->pm3DerivedParameterRho[1] = 0.5/0.4309335108;  
   this->pm3DerivedParameterRho[2] = 0.5/0.3732450359;  
   this->pm3ParameterK[0] = -0.611421 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[1] = -0.093935 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[2] = 0.0;
   this->pm3ParameterK[3] = 0.0;
   this->pm3ParameterL[0] = 1.997272 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[1] = 1.99836 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[2] = 0.00;
   this->pm3ParameterL[3] = 0.00;
   this->pm3ParameterM[0] = 0.794624 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[1] = 1.910677 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[2] = 0.00;
   this->pm3ParameterM[3] = 0.00;
   this->pm3PddgCoreintegralS = -37.882113 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgCoreintegralP = -30.312979 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgOrbitalExponentS = 2.395882;
   this->pm3PddgOrbitalExponentP = 1.742213;
   this->pm3PddgBondingParameterS = -12.676297 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3PddgBondingParameterP = -7.093318 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3PddgAlpha = 2.005294 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3PddgDerivedParameterD[0] = 0.0;    
   this->pm3PddgDerivedParameterD[1] = 0.8939785327;    
   this->pm3PddgDerivedParameterD[2] = 0.9604566451;    
   this->pm3PddgDerivedParameterRho[0] = 0.5/0.2867017837;
   this->pm3PddgDerivedParameterRho[1] = 0.5/0.4757949783;  
   this->pm3PddgDerivedParameterRho[2] = 0.5/0.4135905723;  
   this->pm3PddgParameterK[0] =-0.398055 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterK[1] =-0.079653 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterK[2] = 0.0;
   this->pm3PddgParameterK[3] = 0.0;
   this->pm3PddgParameterL[0] = 1.997272 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3PddgParameterL[1] = 1.99836 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3PddgParameterL[2] = 0.00;
   this->pm3PddgParameterL[3] = 0.00;
   this->pm3PddgParameterM[0] = 0.950073 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterM[1] = 2.336959 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterM[2] = 0.00;
   this->pm3PddgParameterM[3] = 0.00;
   this->pm3PddgParameterPa[0] = 0.462741 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterPa[1] =-0.020444 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterDa[0] = 0.714296 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterDa[1] = 2.041209 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3DCoreintegralS = this->pm3CoreintegralS;
   this->pm3DCoreintegralP = this->pm3CoreintegralP;
   this->pm3DBondingParameterS = this->pm3BondingParameterS;
   this->pm3DBondingParameterP = this->pm3BondingParameterP;
   this->pm3DAlpha = this->pm3Alpha;

}
}
