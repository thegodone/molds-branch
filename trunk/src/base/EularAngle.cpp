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
#include<sstream>
#include<string>
#include<math.h>
#include<stdexcept>
#include"MolDSException.h"
#include"EularAngle.h"
using namespace std;
namespace MolDS_base{

EularAngle::EularAngle(){
   this->SetMessage();
   // e the [BFB_1997] for defenitions of alpha, beta, gamma;
   this->alpha = 0.0; // (= "phi" in P25 in J. A. Pople book)
   this->beta  = 0.0; // (= "theta" in P25 in J. A. Pople book)
   this->gamma = 0.0;
}

EularAngle::EularAngle(double x, double y, double z){
   this->SetMessage();
   double r = 0.0;

   // calc. beta
   if(x==0.0 && y==0.0 && z==0.0){
      stringstream ss;
      ss << this->errorMessageInvalidXYZ;
      throw MolDSException(ss.str());
   }

   r = sqrt( pow(x, 2.0) + pow(y, 2.0) + pow(z, 2.0) );
   this->beta = acos(z/r);

   // calc. alpha
   if(x==0.0 && y==0.0){
      this->alpha = 0.0;
   }
   else{ 
      r = sqrt( pow(x, 2.0) + pow(y, 2.0) );
      this->alpha = atan2(y/r, x/r);
   }

   // set gamma
   this->gamma = 0.0;
   
}

EularAngle::EularAngle(double* angles){
   this->SetMessage();
   this->alpha = angles[0];
   this->beta  = angles[1];
   this->gamma = angles[2];
}

void EularAngle::SetMessage(){
   this->errorMessageInvalidXYZ="Error in base::EularAngle: Invalid coordinates. x=y=z=0.\n";
}



}
