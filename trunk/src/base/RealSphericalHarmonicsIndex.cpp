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
#include<string>
#include<stdexcept>
#include"Enums.h"
#include"MolDSException.h"
#include"RealSphericalHarmonicsIndex.h"
using namespace std;
namespace MolDS_base{

RealSphericalHarmonicsIndex::RealSphericalHarmonicsIndex(){}

RealSphericalHarmonicsIndex::RealSphericalHarmonicsIndex(OrbitalType orbitalType){
   string errorMessageInvalidOrbital = "Error in base::RealSphericalHarmonicIndex::RealSphericalHarmonicIndex: invalid orbitalType. Indicated orbitalType is not prepared\n";

   if(orbitalType == s){
      this->l = 0;
      this->m = 0;
   }
   else if(orbitalType == py){
      this->l = 1;
      this->m = -1;
   }
   else if(orbitalType == pz){
      this->l = 1;
      this->m = 0;
   }
   else if(orbitalType == px){
      this->l = 1;
      this->m = 1;
   }
   else if(orbitalType == dxy){
      this->l = 2;
      this->m = -2;
   }
   else if(orbitalType == dyz){
      this->l = 2;
      this->m = -1;
   }
   else if(orbitalType == dzz){
      this->l = 2;
      this->m = 0;
   }
   else if(orbitalType == dzx){
      this->l = 2;
      this->m = 1;
   }
   else if(orbitalType == dxxyy){
      this->l = 2;
      this->m = 2;
   }
   else{
      stringstream ss;
      ss << errorMessageInvalidOrbital;
      throw MolDSException(ss.str());
   }

}

RealSphericalHarmonicsIndex::RealSphericalHarmonicsIndex(int l, int m){
   this->l = l;
   this->m = m;
}

}
