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
#include"Liatom.h"
using namespace std;
using namespace MolDS_base;

namespace MolDS_base_atoms{
Liatom::Liatom(int index) : Atom(index){
   this->SetAtomicParameters();
}

Liatom::~Liatom(){}

void Liatom::SetAtomicParameters(){
   this->atomType = Li;
   this->atomicMass = 6.941*Parameters::GetInstance()->GetGMolin2AU();
   this->coreCharge = 1.0;
   this->numberValenceElectrons = 1;
   this->valenceShellType = lShell;
   this->valence.push_back(s);
   this->valence.push_back(py);
   this->valence.push_back(pz);
   this->valence.push_back(px);
   for(int i=0; i<this->valence.size();i++){
      this->realSphericalHarmonicsIndeces.push_back(new RealSphericalHarmonicsIndex(this->valence[i]));
   }
   this->bondingParameter = -9.0*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuS = 3.106*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuP = 1.258*Parameters::GetInstance()->GetEV2AU();
   this->imuAmuD = 0.0;
   this->effectiveNuclearChargeK = 2.7;
   this->effectiveNuclearChargeL = 1.3;
   this->effectiveNuclearChargeMsp = 0.0;
   this->effectiveNuclearChargeMd = 0.0;
   this->indoG1 = 0.092012;
   this->indoF2 = 0.049865;
   this->indoF0CoefficientS = 0.5;
   this->indoF0CoefficientP = 0.5;
   this->indoG1CoefficientS = 0.0;
   this->indoG1CoefficientP = -1.0/12.0;
   this->indoF2CoefficientS = 0.0;
   this->indoF2CoefficientP = 0.0;
   this->zindoBondingParameterS = 0.0;
   this->zindoBondingParameterD = 0.0;
   this->zindoF0ss = 0.0;
   this->zindoF0sd = 0.0;        
   this->zindoF0dd = 0.0;              
   this->zindoG1sp = 20194*Parameters::GetInstance()->GetKayser2AU();           
   this->zindoF2pp = 10944*Parameters::GetInstance()->GetKayser2AU();        
   this->zindoG2sd = 0.0;                
   this->zindoG1pd = 0.0;                 
   this->zindoF2pd = 0.0;                
   this->zindoG3pd = 0.0;             
   this->zindoF2dd = 0.0;             
   this->zindoF4dd = 0.0;       
   this->zindoL = 1;
   this->zindoM = 0;
   this->zindoN = 0;
   this->zindoIonPotS = 5.39 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotP = 3.54 * Parameters::GetInstance()->GetEV2AU();
   this->zindoIonPotD = 0.0 * Parameters::GetInstance()->GetEV2AU();
}
}
