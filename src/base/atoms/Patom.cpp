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
   this->am1DerivedParameterD[1] = 0;
   this->am1DerivedParameterD[2] = 0;
   this->am1DerivedParameterRho[0] = 0.5/1;
   this->am1DerivedParameterRho[1] = 0.5/1;
   this->am1DerivedParameterRho[2] = 0.1;
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
   this->rm1DerivedParameterD[0] = 0;    
   this->rm1DerivedParameterD[1] = 0;    
   this->rm1DerivedParameterD[2] = 0;    
   this->rm1DerivedParameterRho[0] = 0.5 / 1;     
   this->rm1DerivedParameterRho[1] = 0.5 / 1;     
   this->rm1DerivedParameterRho[2] = 0.5 / 1;     

}
}
