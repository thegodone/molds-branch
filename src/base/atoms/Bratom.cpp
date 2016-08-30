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
#include"Bratom.h"
using namespace std;
using namespace MolDS_base;
namespace MolDS_base_atoms{
Bratom::Bratom(int index) : Atom(index){
   this->SetAtomicParameters();
}

Bratom::~Bratom(){}

void Bratom::SetAtomicParameters(){
   this->atomType = Br;
   this->atomicMass = 79.90400*Parameters::GetInstance()->GetGMolin2AU();
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
   this->mndoCoreintegralS = -99.986441 * Parameters::GetInstance()->GetEV2AU();
   this->mndoCoreintegralP = -75.671308 * Parameters::GetInstance()->GetEV2AU();
   this->mndoOrbitalExponentS = 3.854302;
   this->mndoOrbitalExponentP = 2.199209;
   this->mndoBondingParameterS = -8.917107 * Parameters::GetInstance()->GetEV2AU();
   this->mndoBondingParameterP = -9.943740 * Parameters::GetInstance()->GetEV2AU();     
   this->mndoAlpha = 2.445705 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->mndoElecEnergyAtom = -346.681250 * Parameters::GetInstance()->GetEV2AU();        
   this->mndoHeatsFormAtom =  26.735 * Parameters::GetInstance()->GetKcalMolin2AU(); //Handbook of Bond Dissociation Energies in Organic Compounds
   this->mndoGss = 15.03643948 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp = 11.27632539 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGsp = 13.03468242 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoGpp2 = 9.85442552 * Parameters::GetInstance()->GetEV2AU();  
   this->mndoHsp =  2.45586832 * Parameters::GetInstance()->GetEV2AU();   
   this->mndoDerivedParameterD[0] =   0;
   this->mndoDerivedParameterD[1] =   0;
   this->mndoDerivedParameterD[2] =   0;    
   this->mndoDerivedParameterRho[0] = 0.5/1;  
   this->mndoDerivedParameterRho[1] = 0.5/1;
   this->mndoDerivedParameterRho[2] = 0.5/1;
   // END OK MOPAC & J Comput Chem 25: 138 –150, 2004
   // OK came from MOPAC 7 & ORIGINAL PUBLICATION
   this->am1CoreintegralS = -104.656063 * Parameters::GetInstance()->GetEV2AU();         
   this->am1CoreintegralP = -74.930052 * Parameters::GetInstance()->GetEV2AU();         
   this->am1OrbitalExponentS = 3.064133;
   this->am1OrbitalExponentP = 2.038333;
   this->am1BondingParameterS = -19.399880 * Parameters::GetInstance()->GetEV2AU();     
   this->am1BondingParameterP = -8.957195 * Parameters::GetInstance()->GetEV2AU();     
   this->am1Alpha = 2.576546/ Parameters::GetInstance()->GetAngstrom2AU();        
   this->am1Gss = this->mndoGss;
   this->am1Gpp = this->mndoGpp;
   this->am1Gsp = this->mndoGsp;
   this->am1Gpp2 = this->mndoGpp2;
   this->am1Hsp = this->mndoHsp;
   this->am1DerivedParameterD[0] = 0.0;
   this->am1DerivedParameterD[1] = 0.8458103844;
   this->am1DerivedParameterD[2] = 1.0407133396;
   this->am1DerivedParameterRho[0] = 0.5/0.5525745655;
   this->am1DerivedParameterRho[1] = 0.5/0.6024436099;
   this->am1DerivedParameterRho[2] = 0.5/0.5307443828;
   this->am1ParameterK[0] = 0.066685 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[1] = 0.025568 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[2] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterK[3] = 0.00 * Parameters::GetInstance()->GetEV2AU();
   this->am1ParameterL[0] = 4.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[1] = 4.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[2] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterL[3] = 0.00 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->am1ParameterM[0] = 1.50 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[1] = 2.30 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[2] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->am1ParameterM[3] = 0.00 * Parameters::GetInstance()->GetAngstrom2AU();
   this->rm1CoreintegralS = -113.4839818 * Parameters::GetInstance()->GetEV2AU();    
   this->rm1CoreintegralP = -76.18720023 * Parameters::GetInstance()->GetEV2AU();    
   this->rm1OrbitalExponentS = 5.73157215;
   this->rm1OrbitalExponentP = 2.03147582;
   this->rm1BondingParameterS = -1.34139841 * Parameters::GetInstance()->GetEV2AU();
   this->rm1BondingParameterP =   -8.20225991 * Parameters::GetInstance()->GetEV2AU();
   this->rm1Alpha = 2.86710532 / Parameters::GetInstance()->GetAngstrom2AU();    
   this->rm1Gss = 17.11563074 * Parameters::GetInstance()->GetEV2AU();
   this->rm1Gpp = 15.62419251 * Parameters::GetInstance()->GetEV2AU();
   this->rm1Gsp = 10.73546293 * Parameters::GetInstance()->GetEV2AU();
   this->rm1Gpp2 = 8.86056199  * Parameters::GetInstance()->GetEV2AU();
   this->rm1Hsp =  2.23512762  * Parameters::GetInstance()->GetEV2AU();
   this->rm1DerivedParameterD[0] = 0;    
   this->rm1DerivedParameterD[1] = 0.2099005325;    
   this->rm1DerivedParameterD[2] = 1.0442262333;    
   this->rm1DerivedParameterRho[0] = 0.5 / 0.6289828275;     
   this->rm1DerivedParameterRho[1] = 0.5 / 1.3165461497;     
   this->rm1DerivedParameterRho[2] = 0.5 / 1.0457055880;     
   this->rm1ParameterK[1] = -0.92731247 * Parameters::GetInstance()->GetEV2AU();
   this->rm1ParameterK[2] = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->rm1ParameterK[3] = 0.0 * Parameters::GetInstance()->GetEV2AU();
   this->rm1ParameterL[0] = 4.28484191 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->rm1ParameterL[1] = 4.5400591  / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->rm1ParameterL[2] = 0.0 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->rm1ParameterL[3] = 0.0 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->rm1ParameterM[0] = 2.00019696 * Parameters::GetInstance()->GetAngstrom2AU();
   this->rm1ParameterM[1] = 2.01617695 * Parameters::GetInstance()->GetAngstrom2AU();
   this->rm1ParameterM[2] = 0.0 * Parameters::GetInstance()->GetAngstrom2AU();
   this->rm1ParameterM[3] = 0.0 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3CoreintegralS = -116.619310 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3CoreintegralP =  -74.227130 * Parameters::GetInstance()->GetEV2AU();         
   this->pm3OrbitalExponentS = 5.348460;      
   this->pm3OrbitalExponentP = 2.127590;      
   this->pm3BondingParameterS = -31.171340 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3BondingParameterP = -6.814010 * Parameters::GetInstance()->GetEV2AU();     
   this->pm3Alpha = 2.511840 / Parameters::GetInstance()->GetAngstrom2AU();        
   this->pm3DerivedParameterD[0] = 0.0;    
   this->pm3DerivedParameterD[1] = 0;    
   this->pm3DerivedParameterD[2] = 0;    
   this->pm3DerivedParameterRho[0] = 0.5/1;
   this->pm3DerivedParameterRho[1] = 0.5/1;  
   this->pm3DerivedParameterRho[2] = 0.5/1;  
   this->pm3ParameterK[0] = 0.960460 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[1] = -0.954920 * Parameters::GetInstance()->GetEV2AU();
   this->pm3ParameterK[2] = 0.0;
   this->pm3ParameterK[3] = 0.0;
   this->pm3ParameterL[0] = 5.976510 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[1] = 5.944700 / pow(Parameters::GetInstance()->GetAngstrom2AU(),2.0);
   this->pm3ParameterL[2] = 0.00;
   this->pm3ParameterL[3] = 0.00;
   this->pm3ParameterM[0] = 5.944700 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[1] = 2.328140 * Parameters::GetInstance()->GetAngstrom2AU();
   this->pm3ParameterM[2] = 0.00;
   this->pm3ParameterM[3] = 0.00;

}
}
