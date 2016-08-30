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
#include"Iatom.h"
using namespace std;
using namespace MolDS_base;
namespace MolDS_base_atoms{
Iatom::Iatom(int index) : Atom(index){
   this->SetAtomicParameters();
}

Iatom::~Iatom(){}

void Iatom::SetAtomicParameters(){
   this->atomType = I;
   this->atomicMass = 126.9045*Parameters::GetInstance()->GetGMolin2AU();
   this->coreCharge = 7.0;
   this->numberValenceElectrons = 7;
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

   // ????
   this->vdWCoefficient = 0.57*Parameters::GetInstance()->GetJ2AU()
                              *pow(Parameters::GetInstance()->GetNm2AU(),6.0)
                              /Parameters::GetInstance()->GetAvogadro();
   this->vdWRadii = 1.430*Parameters::GetInstance()->GetAngstrom2AU();
   this->bondingParameter = -39.0*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuS = 32.272*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 11.080*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.0;
   // END ????
   // what about this ? in ORCA it's look like
   this->effectiveNuclearChargeK = 8.7;
   this->effectiveNuclearChargeL = 5.20;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   // END what about this ? in ORCA it's look like
   // from MOPAC 7 + Extension of the PDDG/PM3 and PDDG/MNDO Semiempirical Molecular Orbital Methods to the Halogens
   // J Comput Chem 25: 138 –150, 2004
   this->mndoCoreintegralS = -100.003054 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralP = -74.611469 * Parameters::GetInstance()->GetEV2AU();
   this->mndoOrbitalExponentS = 2.272961;
   this->mndoOrbitalExponentP = 2.169498;
   this->mndoBondingParameterS = -7.414451 * Parameters::GetInstance()->GetEV2AU();
   this->mndoBondingParameterP = -6.196781 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha = 2.445705 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -340.598360 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom =  25.516 * Parameters::GetInstance()->GetKcalMolin2AU(); //Handbook of Bond Dissociation Energies in Organic Compounds
   this->mndoGss =  15.0404486 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp =  11.1477837 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp =  13.0565580 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 =  9.9140907 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =   2.4563820 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0;
   this->mndoDerivedParameterD[1] =   0;
   this->mndoDerivedParameterD[2] =   0;
   this->mndoDerivedParameterRho[0] = 0.5/1;
   this->mndoDerivedParameterRho[1] = 0.5/1;
   this->mndoDerivedParameterRho[2] = 0.5/1;
   // END OK MOPAC & J Comput Chem 25: 138 –150, 2004
   // OK came from MOPAC 7 & ORIGINAL PUBLICATION

   this->am1CoreintegralS = -103.589663 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -74.429997 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 2.102858;
   this->am1OrbitalExponentP = 2.161153;
   this->am1BondingParameterS = -8.443327 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -6.323405 * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 2.299424/ Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   // need to check publication there
   this->am1DerivedParameterD[0] = 0.0;
   this->am1DerivedParameterD[1] = 0.8458104;
   this->am1DerivedParameterD[2] = 1.0407133;
   this->am1DerivedParameterRho[0] = 0.5/0.5526071;
   this->am1DerivedParameterRho[1] = 0.5/0.6024598;
   this->am1DerivedParameterRho[2] = 0.5/0.5307555;
   // end need  to check publicaiton there
   this->am1ParameterK[0] = 0.004361 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.015706 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 2.30 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 3.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 1.80 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 2.24 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();

   this->rm1CoreintegralS = -74.89997837 * Parameters::GetInstance()->GetEV2AU();    
   this->rm1CoreintegralP = -51.41023805 * Parameters::GetInstance()->GetEV2AU();    
   this->rm1OrbitalExponentS =  2.53003753;
   this->rm1OrbitalExponentP =  2.31738678;
   this->rm1BondingParameterS = -4.19316149 * Parameters::GetInstance()->GetEV2AU();
   this->rm1BondingParameterP =  -4.40038412 * Parameters::GetInstance()->GetEV2AU();
   this->rm1Alpha = 2.14157092  / Parameters::GetInstance()->GetAngstrom2AU();    
   this->rm1Gss = 19.99974131 * Parameters::GetInstance()->GetEV2AU();
   this->rm1Gpp = 7.68957672  * Parameters::GetInstance()->GetEV2AU();
   this->rm1Gsp = 7.30488343  * Parameters::GetInstance()->GetEV2AU();
   this->rm1Gpp2 = 6.85424614  * Parameters::GetInstance()->GetEV2AU();
   this->rm1Hsp = 1.4160294   * Parameters::GetInstance()->GetEV2AU();
   this->rm1DerivedParameterD[0] =  0;    
   this->rm1DerivedParameterD[1] =  1.2963424986;    
   this->rm1DerivedParameterD[2] =  1.1085963456;    
   this->rm1DerivedParameterRho[0] = 0.5/0.7349710934;     
   this->rm1DerivedParameterRho[1] = 0.5/0.3717314674;     
   this->rm1DerivedParameterRho[2] = 0.5/0.4248407280;     
   this->rm1ParameterK[0] = -0.08147724 * Parameters::GetInstance()->GetEV2AU();
   this->rm1ParameterK[1] = 0.05914991  * Parameters::GetInstance()->GetEV2AU();
   this->rm1ParameterK[2] = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->rm1ParameterK[3] = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->rm1ParameterL[0] = 1.56065072 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->rm1ParameterL[1] = 5.76111270  / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->rm1ParameterL[2] = 0.0 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->rm1ParameterL[3] = 0.0 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->rm1ParameterM[0] = 2.00002063 * Parameters::GetInstance()->GetAngstrom2AU();
   this->rm1ParameterM[1] = 2.204888 * Parameters::GetInstance()->GetAngstrom2AU();
   this->rm1ParameterM[2] = 0.0 * Parameters::GetInstance()->GetAngstrom2AU();
   this->rm1ParameterM[3] = 0.0 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3CoreintegralS = -96.454040 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3CoreintegralP =  -61.091580 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3OrbitalExponentS = 7.001010;      
   this->pm3OrbitalExponentP = 2.454350;      
   this->pm3BondingParameterS = -14.494230 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3BondingParameterP = -5.894700 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3Alpha = 1.990190 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3DerivedParameterD[0] = 0.0;    
   this->pm3DerivedParameterD[1] = 0;    
   this->pm3DerivedParameterD[2] = 0;    
   this->pm3DerivedParameterRho[0] = 0.5/1;
   this->pm3DerivedParameterRho[1] = 0.5/1;  
   this->pm3DerivedParameterRho[2] = 0.5/1;  
   this->pm3ParameterK[0] = -0.131480 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[1] = -0.036900 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[2] = 0.0;
   this->pm3ParameterK[3] = 0.0;
   this->pm3ParameterL[0] = 5.206420 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[1] = 6.010120 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[2] = 0.00;
   this->pm3ParameterL[3] = 0.00;
   this->pm3ParameterM[0] = 1.748820 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[1] = 2.710370 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[2] = 0.00;
   this->pm3ParameterM[3] = 0.00;
   this->pm3PddgCoreintegralS = -97.664174 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3PddgCoreintegralP = -61.167137 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3PddgOrbitalExponentS = 5.062801;      
   this->pm3PddgOrbitalExponentP = 2.417757;      
   this->pm3PddgBondingParameterS = -16.592621 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3PddgBondingParameterP = -6.599816 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3PddgAlpha = 1.978170 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3PddgDerivedParameterD[0] = 0.0;    
   this->pm3PddgDerivedParameterD[1] = 0;    
   this->pm3PddgDerivedParameterD[2] = 0;    
   this->pm3PddgDerivedParameterRho[0] = 0.5/1;
   this->pm3PddgDerivedParameterRho[1] = 0.5/1;  
   this->pm3PddgDerivedParameterRho[2] = 0.5/1;  
   this->pm3PddgParameterK[0] =-0.136003 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterK[1] =-0.037287 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterK[2] = 0.0;
   this->pm3PddgParameterK[3] = 0.0;
   this->pm3PddgParameterL[0] = 3.852912 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3PddgParameterL[1] = 5.229264 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3PddgParameterL[2] = 0.00;
   this->pm3PddgParameterL[3] = 0.00;
   this->pm3PddgParameterM[0] = 1.697455 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterM[1] = 2.768669 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterM[2] = 0.00;
   this->pm3PddgParameterM[3] = 0.00;
   this->pm3PddgParameterPa[0] = 0.012901 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterPa[1] =-0.012825 * Parameters::GetInstance()->GetEV2AU();
   this->pm3PddgParameterDa[0] = 1.994299 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3PddgParameterDa[1] = 2.263417 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3DCoreintegralS = this->pm3CoreintegralS;
   this->pm3DCoreintegralP = this->pm3CoreintegralP;
   this->pm3DBondingParameterS = this->pm3BondingParameterS;
   this->pm3DBondingParameterP = this->pm3BondingParameterP;
   this->pm3DAlpha = this->pm3Alpha;

}
}
