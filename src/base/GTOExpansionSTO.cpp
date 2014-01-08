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
#include<string>
#include<stdexcept>
#include"Enums.h"
#include"Uncopyable.h"
#include"MolDSException.h"
#include"GTOExpansionSTO.h"
using namespace std;
namespace MolDS_base{

GTOExpansionSTO* GTOExpansionSTO::gTOExpansionSTO = NULL;

GTOExpansionSTO::GTOExpansionSTO(){
   this->SetCoefficientsExponents();
   this->errorMessageGetCoefficientNonValidOrbital 
      = "Error in base::GTOExpansionSTO::GetCoefficient: Non available orbital is contained.\n";
   this->errorMessageGetExponentNonValidOrbital 
      = "Error in base::GTOExpansionSTO::GetExponent: Non available orbital is contained.\n";
   this->errorMessageOrbitalType = "\torbital type = ";
   this->errorMessageSTOnGType = "\tSTOnG type = ";
}

GTOExpansionSTO::~GTOExpansionSTO(){
}

GTOExpansionSTO* GTOExpansionSTO::GetInstance(){
   if(gTOExpansionSTO == NULL){
      gTOExpansionSTO = new GTOExpansionSTO();
   }
   return gTOExpansionSTO;
}

void GTOExpansionSTO::DeleteInstance(){
   if(gTOExpansionSTO != NULL){
      delete gTOExpansionSTO; 
   }
   gTOExpansionSTO = NULL;
}

double GTOExpansionSTO::GetExponent(STOnGType stonG, ShellType shellType, OrbitalType orbitalType, int index) const{

   AzimuthalType azimuthalType;
   if(orbitalType == s){
      azimuthalType = sAzimuthal;
   }
   else if(orbitalType == px || orbitalType == py ||orbitalType == pz ){
      azimuthalType = pAzimuthal;
   }
   else if(orbitalType == dxy || orbitalType == dyz ||orbitalType == dzz || orbitalType == dzx ||orbitalType == dxxyy ){
      azimuthalType = dAzimuthal;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetExponentNonValidOrbital;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << endl;
      ss << this->errorMessageSTOnGType << STOnGTypeStr(stonG) << endl;
      throw MolDSException(ss.str());
   }

   return this->exponents[stonG][shellType][azimuthalType][index];
}

double GTOExpansionSTO::GetCoefficient(STOnGType stonG, ShellType shellType, OrbitalType orbitalType, int index) const{

   AzimuthalType azimuthalType;
   if(orbitalType == s){
      azimuthalType = sAzimuthal;
   }
   else if(orbitalType == px || orbitalType == py ||orbitalType == pz ){
      azimuthalType = pAzimuthal;
   }
   else if(orbitalType == dxy || orbitalType == dyz ||orbitalType == dzz || orbitalType == dzx ||orbitalType == dxxyy ){
      azimuthalType = dAzimuthal;
   }
   else{
      stringstream ss;
      ss << this->errorMessageGetCoefficientNonValidOrbital;
      ss << this->errorMessageOrbitalType << OrbitalTypeStr(orbitalType) << endl;
      ss << this->errorMessageSTOnGType << STOnGTypeStr(stonG) << endl;
      throw MolDSException(ss.str());
   }

   return this->coefficients[stonG][shellType][azimuthalType][index];
}

//  see Table I and II in [S_1970]
void GTOExpansionSTO::SetCoefficientsExponents(){

   //STO-1G, 
   {   
      // 1s
      exponents[STO1G][kShell][sAzimuthal][0] = 2.709498091e-1;   coefficients[STO1G][kShell][sAzimuthal][0] = 1.0000;
      // 2s
      exponents[STO1G][lShell][sAzimuthal][0] = 1.012151084e-1;   coefficients[STO1G][lShell][sAzimuthal][0] = 1.0000;
      // 2p
      exponents[STO1G][lShell][pAzimuthal][0] = 1.759666885e-1;   coefficients[STO1G][lShell][pAzimuthal][0] = 1.0000;
      // 3s
      exponents[STO1G][mShell][sAzimuthal][0] = 5.296881757e-2;   coefficients[STO1G][mShell][sAzimuthal][0] = 1.0000;
      // 3p
      exponents[STO1G][mShell][pAzimuthal][0] = 9.113614253e-2;   coefficients[STO1G][mShell][pAzimuthal][0] = 1.0000;
      // 3d
      exponents[STO1G][mShell][dAzimuthal][0] = 1.302270363e-1;   coefficients[STO1G][mShell][dAzimuthal][0] = 1.0000;
      // 4s
      exponents[STO1G][nShell][sAzimuthal][0] = 3.264600274e-2;   coefficients[STO1G][nShell][sAzimuthal][0] = 1.0000;
      // 4p
      exponents[STO1G][nShell][pAzimuthal][0] = 5.578350235e-2;   coefficients[STO1G][nShell][pAzimuthal][0] = 1.0000;
      // 4d
      exponents[STO1G][nShell][dAzimuthal][0] = 7.941656339e-2;   coefficients[STO1G][nShell][dAzimuthal][0] = 1.0000;
   }

   //STO-2G, 
   {
      // 1s
      exponents[STO2G][kShell][sAzimuthal][0] = 8.518186635e-1;   coefficients[STO2G][kShell][sAzimuthal][0] = 4.301284983e-1; 
      exponents[STO2G][kShell][sAzimuthal][1] = 1.516232927e-1;   coefficients[STO2G][kShell][sAzimuthal][1] = 6.789135305e-1; 
      // 2s
      exponents[STO2G][lShell][sAzimuthal][0] = 1.292278611e-1;   coefficients[STO2G][lShell][sAzimuthal][0] = 7.470867124e-1; 
      exponents[STO2G][lShell][sAzimuthal][1] = 4.908584205e-2;   coefficients[STO2G][lShell][sAzimuthal][1] = 2.855980556e-1; 
      // 2p
      exponents[STO2G][lShell][pAzimuthal][0] = 4.323908358e-1;   coefficients[STO2G][lShell][pAzimuthal][0] = 4.522627513e-1; 
      exponents[STO2G][lShell][pAzimuthal][1] = 1.069439065e-1;   coefficients[STO2G][lShell][pAzimuthal][1] = 6.713122642e-1; 
      // 3s
      exponents[STO2G][mShell][sAzimuthal][0] = 6.694095822e-1;   coefficients[STO2G][mShell][sAzimuthal][0] =-1.529645716e-1; 
      exponents[STO2G][mShell][sAzimuthal][1] = 5.837135094e-2;   coefficients[STO2G][mShell][sAzimuthal][1] = 1.051370110; 
      // 3p
      exponents[STO2G][mShell][pAzimuthal][0] = 1.458620964e-1;   coefficients[STO2G][mShell][pAzimuthal][0] = 5.349653144e-1; 
      exponents[STO2G][mShell][pAzimuthal][1] = 5.664210742e-2;   coefficients[STO2G][mShell][pAzimuthal][1] = 5.299607212e-1; 
      // 3d
      exponents[STO2G][mShell][dAzimuthal][0] = 2.777427345e-1;   coefficients[STO2G][mShell][dAzimuthal][0] = 4.666137923e-1; 
      exponents[STO2G][mShell][dAzimuthal][1] = 8.336507714e-2;   coefficients[STO2G][mShell][dAzimuthal][1] = 6.644706516e-1; 
      // 4s
      exponents[STO2G][nShell][sAzimuthal][0] = 2.441785453e-1;   coefficients[STO2G][nShell][sAzimuthal][0] =-3.046656896e-1; 
      exponents[STO2G][nShell][sAzimuthal][1] = 4.051097664e-2;   coefficients[STO2G][nShell][sAzimuthal][1] = 1.146877294e00; 
      // 4p                                                               
      exponents[STO2G][nShell][pAzimuthal][0] = 6.190052680e-2;   coefficients[STO2G][nShell][pAzimuthal][0] = 8.743116767e-1; 
      exponents[STO2G][nShell][pAzimuthal][1] = 2.648418407e-2;   coefficients[STO2G][nShell][pAzimuthal][1] = 1.513640107e-1; 
      // 4d                                                               
      exponents[STO2G][nShell][dAzimuthal][0] = 1.330958892e-1;   coefficients[STO2G][nShell][dAzimuthal][0] = 4.932764167e-1; 
      exponents[STO2G][nShell][dAzimuthal][1] = 5.272119659e-2;   coefficients[STO2G][nShell][dAzimuthal][1] = 5.918727866e-1; 
   }

   //STO-3G, 
   {
      // 1s
      exponents[STO3G][kShell][sAzimuthal][0] = 2.227660584e00;   coefficients[STO3G][kShell][sAzimuthal][0] = 1.543289673e-1; 
      exponents[STO3G][kShell][sAzimuthal][1] = 4.057711562e-1;   coefficients[STO3G][kShell][sAzimuthal][1] = 5.353281423e-1; 
      exponents[STO3G][kShell][sAzimuthal][2] = 1.098175104e-1;   coefficients[STO3G][kShell][sAzimuthal][2] = 4.446345422e-1; 
      // 2s
      exponents[STO3G][lShell][sAzimuthal][0] = 2.581578398e00;   coefficients[STO3G][lShell][sAzimuthal][0] =-5.994474934e-2; 
      exponents[STO3G][lShell][sAzimuthal][1] = 1.567622104e-1;   coefficients[STO3G][lShell][sAzimuthal][1] = 5.960385398e-1; 
      exponents[STO3G][lShell][sAzimuthal][2] = 6.018332272e-2;   coefficients[STO3G][lShell][sAzimuthal][2] = 4.581786291e-1; 
      // 2p
      exponents[STO3G][lShell][pAzimuthal][0] = 9.192379002e-1;   coefficients[STO3G][lShell][pAzimuthal][0] = 1.623948553e-1; 
      exponents[STO3G][lShell][pAzimuthal][1] = 2.359194503e-1;   coefficients[STO3G][lShell][pAzimuthal][1] = 5.661708862e-1; 
      exponents[STO3G][lShell][pAzimuthal][2] = 8.009805746e-2;   coefficients[STO3G][lShell][pAzimuthal][2] = 4.223071752e-1; 
      // 3s
      exponents[STO3G][mShell][sAzimuthal][0] = 5.641487709e-1;   coefficients[STO3G][mShell][sAzimuthal][0] =-1.782577972e-1; 
      exponents[STO3G][mShell][sAzimuthal][1] = 6.924421391e-2;   coefficients[STO3G][mShell][sAzimuthal][1] = 8.612761663e-1; 
      exponents[STO3G][mShell][sAzimuthal][2] = 3.269529097e-2;   coefficients[STO3G][mShell][sAzimuthal][2] = 2.261841969e-1; 
      // 3p
      exponents[STO3G][mShell][pAzimuthal][0] = 2.692880368e00;   coefficients[STO3G][mShell][pAzimuthal][0] =-1.061945788e-2; 
      exponents[STO3G][mShell][pAzimuthal][1] = 1.489359592e-1;   coefficients[STO3G][mShell][pAzimuthal][1] = 5.218564264e-1; 
      exponents[STO3G][mShell][pAzimuthal][2] = 5.739585040e-2;   coefficients[STO3G][mShell][pAzimuthal][2] = 5.450015143e-1; 
      // 3d
      exponents[STO3G][mShell][dAzimuthal][0] = 5.229112225e-1;   coefficients[STO3G][mShell][dAzimuthal][0] = 1.686596060e-1; 
      exponents[STO3G][mShell][dAzimuthal][1] = 1.639595876e-1;   coefficients[STO3G][mShell][dAzimuthal][1] = 5.847984817e-1; 
      exponents[STO3G][mShell][dAzimuthal][2] = 6.386630021e-2;   coefficients[STO3G][mShell][dAzimuthal][2] = 4.056779523e-1; 
      // 4s
      exponents[STO3G][nShell][sAzimuthal][0] = 2.267938753e-1;   coefficients[STO3G][nShell][sAzimuthal][0] =-3.349048323e-1; 
      exponents[STO3G][nShell][sAzimuthal][1] = 4.448178019e-2;   coefficients[STO3G][nShell][sAzimuthal][1] = 1.056744667e00; 
      exponents[STO3G][nShell][sAzimuthal][2] = 2.195294664e-2;   coefficients[STO3G][nShell][sAzimuthal][2] = 1.256661680e-1; 
      // 4p                                                               
      exponents[STO3G][nShell][pAzimuthal][0] = 4.859692220e-1;   coefficients[STO3G][nShell][pAzimuthal][0] =-6.147823411e-2; 
      exponents[STO3G][nShell][pAzimuthal][1] = 7.430216918e-2;   coefficients[STO3G][nShell][pAzimuthal][1] = 6.604172234e-1; 
      exponents[STO3G][nShell][pAzimuthal][2] = 3.653340923e-2;   coefficients[STO3G][nShell][pAzimuthal][2] = 3.932639495e-1; 
      // 4d                                                               
      exponents[STO3G][nShell][dAzimuthal][0] = 1.777717219e-1;   coefficients[STO3G][nShell][dAzimuthal][0] = 2.308552718e-1; 
      exponents[STO3G][nShell][dAzimuthal][1] = 8.040647350e-2;   coefficients[STO3G][nShell][dAzimuthal][1] = 6.042409177e-1; 
      exponents[STO3G][nShell][dAzimuthal][2] = 3.949855551e-2;   coefficients[STO3G][nShell][dAzimuthal][2] = 2.595768926e-1; 
   }

   //STO-4G, 
   {
      // 1s
      exponents[STO4G][kShell][sAzimuthal][0] = 5.216844534e00;   coefficients[STO4G][kShell][sAzimuthal][0] = 5.675242080e-2; 
      exponents[STO4G][kShell][sAzimuthal][1] = 9.546182760e-1;   coefficients[STO4G][kShell][sAzimuthal][1] = 2.601413550e-1; 
      exponents[STO4G][kShell][sAzimuthal][2] = 2.652034102e-1;   coefficients[STO4G][kShell][sAzimuthal][2] = 5.328461143e-1; 
      exponents[STO4G][kShell][sAzimuthal][3] = 8.801862774e-2;   coefficients[STO4G][kShell][sAzimuthal][3] = 2.916254405e-1; 
      // 2s
      exponents[STO4G][lShell][sAzimuthal][0] = 1.161525551e01;   coefficients[STO4G][lShell][sAzimuthal][0] =-1.198411747e-2; 
      exponents[STO4G][lShell][sAzimuthal][1] = 2.000243111e00;   coefficients[STO4G][lShell][sAzimuthal][1] =-5.472052539e-2; 
      exponents[STO4G][lShell][sAzimuthal][2] = 1.607280687e-1;   coefficients[STO4G][lShell][sAzimuthal][2] = 5.805587176e-1; 
      exponents[STO4G][lShell][sAzimuthal][3] = 6.125744532e-2;   coefficients[STO4G][lShell][sAzimuthal][3] = 4.770079976e-1; 
      // 2p
      exponents[STO4G][lShell][pAzimuthal][0] = 1.798260992e00;   coefficients[STO4G][lShell][pAzimuthal][0] = 5.713170255e-2; 
      exponents[STO4G][lShell][pAzimuthal][1] = 4.662622228e-1;   coefficients[STO4G][lShell][pAzimuthal][1] = 2.857455515e-1; 
      exponents[STO4G][lShell][pAzimuthal][2] = 1.643718620e-1;   coefficients[STO4G][lShell][pAzimuthal][2] = 5.517873105e-1; 
      exponents[STO4G][lShell][pAzimuthal][3] = 6.543927065e-2;   coefficients[STO4G][lShell][pAzimuthal][3] = 2.632314924e-1; 
      // 3s
      exponents[STO4G][mShell][sAzimuthal][0] = 1.513265591e00;   coefficients[STO4G][mShell][sAzimuthal][0] =-3.295496352e-2; 
      exponents[STO4G][mShell][sAzimuthal][1] = 4.262497508e-1;   coefficients[STO4G][mShell][sAzimuthal][1] =-1.724516959e-1; 
      exponents[STO4G][mShell][sAzimuthal][2] = 7.643320863e-2;   coefficients[STO4G][mShell][sAzimuthal][2] = 7.518511194e-1; 
      exponents[STO4G][mShell][sAzimuthal][3] = 3.760545063e-2;   coefficients[STO4G][mShell][sAzimuthal][3] = 3.589627317e-1; 
      // 3p
      exponents[STO4G][mShell][pAzimuthal][0] = 1.853180239e00;   coefficients[STO4G][mShell][pAzimuthal][0] =-1.434249391e-2; 
      exponents[STO4G][mShell][pAzimuthal][1] = 1.915075719e-1;   coefficients[STO4G][mShell][pAzimuthal][1] = 2.755177589e-1; 
      exponents[STO4G][mShell][pAzimuthal][2] = 8.655487938e-2;   coefficients[STO4G][mShell][pAzimuthal][2] = 5.846750879e-1; 
      exponents[STO4G][mShell][pAzimuthal][3] = 4.184253862e-2;   coefficients[STO4G][mShell][pAzimuthal][3] = 2.144986514e-1; 
      // 3d
      exponents[STO4G][mShell][dAzimuthal][0] = 9.185846715e-1;   coefficients[STO4G][mShell][dAzimuthal][0] = 5.799057705e-2; 
      exponents[STO4G][mShell][dAzimuthal][1] = 2.920461109e-1;   coefficients[STO4G][mShell][dAzimuthal][1] = 3.045581349e-1; 
      exponents[STO4G][mShell][dAzimuthal][2] = 1.187568890e-1;   coefficients[STO4G][mShell][dAzimuthal][2] = 5.601358038e-1; 
      exponents[STO4G][mShell][dAzimuthal][3] = 5.286755896e-2;   coefficients[STO4G][mShell][dAzimuthal][3] = 2.432423313e-1; 
      // 4s
      exponents[STO4G][nShell][sAzimuthal][0] = 3.242212833e-1;   coefficients[STO4G][nShell][sAzimuthal][0] =-1.120682822e-1; 
      exponents[STO4G][nShell][sAzimuthal][1] = 1.663217177e-1;   coefficients[STO4G][nShell][sAzimuthal][1] =-2.845426863e-1; 
      exponents[STO4G][nShell][sAzimuthal][2] = 5.081097451e-2;   coefficients[STO4G][nShell][sAzimuthal][2] = 8.909873788e-1; 
      exponents[STO4G][nShell][sAzimuthal][3] = 2.829066600e-2;   coefficients[STO4G][nShell][sAzimuthal][3] = 3.517811205e-1; 
      // 4p                                                               
      exponents[STO4G][nShell][pAzimuthal][0] = 1.492607880e00;   coefficients[STO4G][nShell][pAzimuthal][0] =-6.035216774e-3; 
      exponents[STO4G][nShell][pAzimuthal][1] = 4.327619272e-1;   coefficients[STO4G][nShell][pAzimuthal][1] =-6.013310874e-2; 
      exponents[STO4G][nShell][pAzimuthal][2] = 7.553156064e-2;   coefficients[STO4G][nShell][pAzimuthal][2] = 6.451518200e-1; 
      exponents[STO4G][nShell][pAzimuthal][3] = 3.706272183e-2;   coefficients[STO4G][nShell][pAzimuthal][3] = 4.117923820e-1; 
      // 4d                                                               
      exponents[STO4G][nShell][dAzimuthal][0] = 1.995825422e00;   coefficients[STO4G][nShell][dAzimuthal][0] =-2.816702620e-3; 
      exponents[STO4G][nShell][dAzimuthal][1] = 1.823461280e-1;   coefficients[STO4G][nShell][dAzimuthal][1] = 2.177095871e-1; 
      exponents[STO4G][nShell][dAzimuthal][2] = 8.197240896e-2;   coefficients[STO4G][nShell][dAzimuthal][2] = 6.058047348e-1; 
      exponents[STO4G][nShell][dAzimuthal][3] = 4.000634951e-2;   coefficients[STO4G][nShell][dAzimuthal][3] = 2.717811257e-1; 
   }

   //STO-5G, 
   {
      // 1s
      exponents[STO5G][kShell][sAzimuthal][0] = 1.130563696e01;   coefficients[STO5G][kShell][sAzimuthal][0] = 2.214055312e-2; 
      exponents[STO5G][kShell][sAzimuthal][1] = 2.071728178e00;   coefficients[STO5G][kShell][sAzimuthal][1] = 1.135411520e-1; 
      exponents[STO5G][kShell][sAzimuthal][2] = 5.786484833e-1;   coefficients[STO5G][kShell][sAzimuthal][2] = 3.318161484e-1; 
      exponents[STO5G][kShell][sAzimuthal][3] = 1.975724573e-1;   coefficients[STO5G][kShell][sAzimuthal][3] = 4.825700713e-1; 
      exponents[STO5G][kShell][sAzimuthal][4] = 7.445271746e-2;   coefficients[STO5G][kShell][sAzimuthal][4] = 1.935721966e-1; 
      // 2s
      exponents[STO5G][lShell][sAzimuthal][0] = 8.984956862e00;   coefficients[STO5G][lShell][sAzimuthal][0] =-1.596349096e-2; 
      exponents[STO5G][lShell][sAzimuthal][1] = 1.673710636e00;   coefficients[STO5G][lShell][sAzimuthal][1] =-5.685884883e-2; 
      exponents[STO5G][lShell][sAzimuthal][2] = 1.944726668e-1;   coefficients[STO5G][lShell][sAzimuthal][2] = 3.698265599e-1; 
      exponents[STO5G][lShell][sAzimuthal][3] = 8.806345634e-2;   coefficients[STO5G][lShell][sAzimuthal][3] = 5.480512593e-1; 
      exponents[STO5G][lShell][sAzimuthal][4] = 4.249068522e-2;   coefficients[STO5G][lShell][sAzimuthal][4] = 1.472634893e-1; 
      // 2p
      exponents[STO5G][lShell][pAzimuthal][0] = 3.320386533e00;   coefficients[STO5G][lShell][pAzimuthal][0] = 2.079051117e-2; 
      exponents[STO5G][lShell][pAzimuthal][1] = 8.643257633e-1;   coefficients[STO5G][lShell][pAzimuthal][1] = 1.235472099e-1; 
      exponents[STO5G][lShell][pAzimuthal][2] = 3.079819284e-1;   coefficients[STO5G][lShell][pAzimuthal][2] = 3.667738986e-1; 
      exponents[STO5G][lShell][pAzimuthal][3] = 1.273309895e-1;   coefficients[STO5G][lShell][pAzimuthal][3] = 4.834930290e-1; 
      exponents[STO5G][lShell][pAzimuthal][4] = 5.606243164e-2;   coefficients[STO5G][lShell][pAzimuthal][4] = 1.653444074e-1; 
      // 3s
      exponents[STO5G][mShell][sAzimuthal][0] = 4.275877914e00;   coefficients[STO5G][mShell][sAzimuthal][0] =-3.920358850e-3; 
      exponents[STO5G][mShell][sAzimuthal][1] = 1.132409433e00;   coefficients[STO5G][mShell][sAzimuthal][1] =-4.168430506e-2; 
      exponents[STO5G][mShell][sAzimuthal][2] = 4.016256968e-1;   coefficients[STO5G][mShell][sAzimuthal][2] =-1.637440990e-1; 
      exponents[STO5G][mShell][sAzimuthal][3] = 7.732370620e-2;   coefficients[STO5G][mShell][sAzimuthal][3] = 7.419373723e-1; 
      exponents[STO5G][mShell][sAzimuthal][4] = 3.800708627e-2;   coefficients[STO5G][mShell][sAzimuthal][4] = 3.724364929e-1; 
      // 3p
      exponents[STO5G][mShell][pAzimuthal][0] = 6.466803859e00;   coefficients[STO5G][mShell][pAzimuthal][0] =-2.329023747e-3; 
      exponents[STO5G][mShell][pAzimuthal][1] = 1.555914802e00;   coefficients[STO5G][mShell][pAzimuthal][1] =-1.357395221e-2; 
      exponents[STO5G][mShell][pAzimuthal][2] = 1.955925255e-1;   coefficients[STO5G][mShell][pAzimuthal][2] = 2.632185383e-1; 
      exponents[STO5G][mShell][pAzimuthal][3] = 8.809647701e-2;   coefficients[STO5G][mShell][pAzimuthal][3] = 5.880427024e-1; 
      exponents[STO5G][mShell][pAzimuthal][4] = 4.234835707e-2;   coefficients[STO5G][mShell][pAzimuthal][4] = 2.242794445e-1; 
      // 3d
      exponents[STO5G][mShell][dAzimuthal][0] = 1.539033958e00;   coefficients[STO5G][mShell][dAzimuthal][0] = 2.020869128e-2; 
      exponents[STO5G][mShell][dAzimuthal][1] = 4.922090297e-1;   coefficients[STO5G][mShell][dAzimuthal][1] = 1.321157923e-1; 
      exponents[STO5G][mShell][dAzimuthal][2] = 2.029756928e-1;   coefficients[STO5G][mShell][dAzimuthal][2] = 3.911240346e-1; 
      exponents[STO5G][mShell][dAzimuthal][3] = 9.424112917e-2;   coefficients[STO5G][mShell][dAzimuthal][3] = 4.779609701e-1; 
      exponents[STO5G][mShell][dAzimuthal][4] = 4.569058269e-2;   coefficients[STO5G][mShell][dAzimuthal][4] = 1.463662294e-1; 
      // 4s
      exponents[STO5G][nShell][sAzimuthal][0] = 2.980263783e00;   coefficients[STO5G][nShell][sAzimuthal][0] = 1.513948997e-3; 
      exponents[STO5G][nShell][sAzimuthal][1] = 3.792228833e-1;   coefficients[STO5G][nShell][sAzimuthal][1] =-7.316801518e-2; 
      exponents[STO5G][nShell][sAzimuthal][2] = 1.789717224e-1;   coefficients[STO5G][nShell][sAzimuthal][2] =-3.143703799e-1; 
      exponents[STO5G][nShell][sAzimuthal][3] = 5.002110360e-2;   coefficients[STO5G][nShell][sAzimuthal][3] = 9.032615169e-1; 
      exponents[STO5G][nShell][sAzimuthal][4] = 2.789361681e-2;   coefficients[STO5G][nShell][sAzimuthal][4] = 3.294210848e-1; 
      // 4p                                                               
      exponents[STO5G][nShell][pAzimuthal][0] = 1.091977298e00;   coefficients[STO5G][nShell][pAzimuthal][0] =-1.143929558e-2; 
      exponents[STO5G][nShell][pAzimuthal][1] = 3.719985051e-1;   coefficients[STO5G][nShell][pAzimuthal][1] =-6.322651538e-2; 
      exponents[STO5G][nShell][pAzimuthal][2] = 8.590019352e-2;   coefficients[STO5G][nShell][pAzimuthal][2] = 4.398907721e-1; 
      exponents[STO5G][nShell][pAzimuthal][3] = 4.786503860e-2;   coefficients[STO5G][nShell][pAzimuthal][3] = 5.245859166e-1; 
      exponents[STO5G][nShell][pAzimuthal][4] = 2.730479990e-2;   coefficients[STO5G][nShell][pAzimuthal][4] = 1.017072253e-1; 
      // 4d                                                               
      exponents[STO5G][nShell][dAzimuthal][0] = 1.522122079e00;   coefficients[STO5G][nShell][dAzimuthal][0] =-3.673711876e-3; 
      exponents[STO5G][nShell][dAzimuthal][1] = 2.173041823e-1;   coefficients[STO5G][nShell][dAzimuthal][1] = 1.167122499e-1; 
      exponents[STO5G][nShell][dAzimuthal][2] = 1.084876577e-1;   coefficients[STO5G][nShell][dAzimuthal][2] = 4.216476416e-1; 
      exponents[STO5G][nShell][dAzimuthal][3] = 5.836797641e-2;   coefficients[STO5G][nShell][dAzimuthal][3] = 4.547673415e-1; 
      exponents[STO5G][nShell][dAzimuthal][4] = 3.206682246e-2;   coefficients[STO5G][nShell][dAzimuthal][4] = 1.037803318e-1; 
   }

   //STO-6G, 
   {
      // 1s
      exponents[STO6G][kShell][sAzimuthal][0] = 2.310303149e01;   coefficients[STO6G][kShell][sAzimuthal][0] = 9.163596280e-3; 
      exponents[STO6G][kShell][sAzimuthal][1] = 4.235915534e00;   coefficients[STO6G][kShell][sAzimuthal][1] = 4.936149294e-2; 
      exponents[STO6G][kShell][sAzimuthal][2] = 1.185056519e00;   coefficients[STO6G][kShell][sAzimuthal][2] = 1.685383049e-1; 
      exponents[STO6G][kShell][sAzimuthal][3] = 4.070988982e-1;   coefficients[STO6G][kShell][sAzimuthal][3] = 3.705627997e-1; 
      exponents[STO6G][kShell][sAzimuthal][4] = 1.580884151e-1;   coefficients[STO6G][kShell][sAzimuthal][4] = 4.164915298e-1; 
      exponents[STO6G][kShell][sAzimuthal][5] = 6.510953954e-2;   coefficients[STO6G][kShell][sAzimuthal][5] = 1.303340841e-1; 
      // 2s
      exponents[STO6G][lShell][sAzimuthal][0] = 2.768496241e01;   coefficients[STO6G][lShell][sAzimuthal][0] =-4.151277819e-3; 
      exponents[STO6G][lShell][sAzimuthal][1] = 5.077140627e00;   coefficients[STO6G][lShell][sAzimuthal][1] =-2.067024148e-2; 
      exponents[STO6G][lShell][sAzimuthal][2] = 1.426786050e00;   coefficients[STO6G][lShell][sAzimuthal][2] =-5.150303337e-2; 
      exponents[STO6G][lShell][sAzimuthal][3] = 2.040335729e-1;   coefficients[STO6G][lShell][sAzimuthal][3] = 3.346271174e-1; 
      exponents[STO6G][lShell][sAzimuthal][4] = 9.260298399e-2;   coefficients[STO6G][lShell][sAzimuthal][4] = 5.621061301e-1; 
      exponents[STO6G][lShell][sAzimuthal][5] = 4.416183978e-2;   coefficients[STO6G][lShell][sAzimuthal][5] = 1.712994697e-1; 
      // 2p
      exponents[STO6G][lShell][pAzimuthal][0] = 5.868285913e00;   coefficients[STO6G][lShell][pAzimuthal][0] = 7.924233646e-3; 
      exponents[STO6G][lShell][pAzimuthal][1] = 1.530329631e00;   coefficients[STO6G][lShell][pAzimuthal][1] = 5.144104825e-2; 
      exponents[STO6G][lShell][pAzimuthal][2] = 5.475665231e-1;   coefficients[STO6G][lShell][pAzimuthal][2] = 1.898400060e-1; 
      exponents[STO6G][lShell][pAzimuthal][3] = 2.288932733e-1;   coefficients[STO6G][lShell][pAzimuthal][3] = 4.049863191e-1; 
      exponents[STO6G][lShell][pAzimuthal][4] = 1.046655969e-1;   coefficients[STO6G][lShell][pAzimuthal][4] = 4.012362861e-1; 
      exponents[STO6G][lShell][pAzimuthal][5] = 4.948220127e-2;   coefficients[STO6G][lShell][pAzimuthal][5] = 1.051855189e-1; 
      // 3s
      exponents[STO6G][mShell][sAzimuthal][0] = 3.273031938e00;   coefficients[STO6G][mShell][sAzimuthal][0] =-6.775596947e-3; 
      exponents[STO6G][mShell][sAzimuthal][1] = 9.200611311e-1;   coefficients[STO6G][mShell][sAzimuthal][1] =-5.639325779e-2; 
      exponents[STO6G][mShell][sAzimuthal][2] = 3.593349765e-1;   coefficients[STO6G][mShell][sAzimuthal][2] =-1.587656086e-1; 
      exponents[STO6G][mShell][sAzimuthal][3] = 8.636686991e-2;   coefficients[STO6G][mShell][sAzimuthal][3] = 5.534527651e-1; 
      exponents[STO6G][mShell][sAzimuthal][4] = 4.797373812e-2;   coefficients[STO6G][mShell][sAzimuthal][4] = 5.015351020e-1; 
      exponents[STO6G][mShell][sAzimuthal][5] = 2.724741144e-2;   coefficients[STO6G][mShell][sAzimuthal][5] = 7.223633674e-2; 
      // 3p
      exponents[STO6G][mShell][pAzimuthal][0] = 5.077973607e00;   coefficients[STO6G][mShell][pAzimuthal][0] =-3.329929840e-3; 
      exponents[STO6G][mShell][pAzimuthal][1] = 1.340786940e00;   coefficients[STO6G][mShell][pAzimuthal][1] =-1.419488340e-2; 
      exponents[STO6G][mShell][pAzimuthal][2] = 2.248434849e-1;   coefficients[STO6G][mShell][pAzimuthal][2] = 1.639395770e-1; 
      exponents[STO6G][mShell][pAzimuthal][3] = 1.131741848e-1;   coefficients[STO6G][mShell][pAzimuthal][3] = 4.485358256e-1; 
      exponents[STO6G][mShell][pAzimuthal][4] = 6.076408893e-2;   coefficients[STO6G][mShell][pAzimuthal][4] = 3.908813050e-1; 
      exponents[STO6G][mShell][pAzimuthal][5] = 3.315424265e-2;   coefficients[STO6G][mShell][pAzimuthal][5] = 7.411456232e-2; 
      // 3d
      exponents[STO6G][mShell][dAzimuthal][0] = 2.488296923e00;   coefficients[STO6G][mShell][dAzimuthal][0] = 7.283828112e-3; 
      exponents[STO6G][mShell][dAzimuthal][1] = 7.981487853e-1;   coefficients[STO6G][mShell][dAzimuthal][1] = 5.386799363e-2; 
      exponents[STO6G][mShell][dAzimuthal][2] = 3.311327490e-1;   coefficients[STO6G][mShell][dAzimuthal][2] = 2.072139149e-1; 
      exponents[STO6G][mShell][dAzimuthal][3] = 1.559114463e-1;   coefficients[STO6G][mShell][dAzimuthal][3] = 4.266269092e-1; 
      exponents[STO6G][mShell][dAzimuthal][4] = 7.877734732e-2;   coefficients[STO6G][mShell][dAzimuthal][4] = 3.843100204e-1; 
      exponents[STO6G][mShell][dAzimuthal][5] = 4.058484363e-2;   coefficients[STO6G][mShell][dAzimuthal][5] = 8.902827546e-2; 
      // 4s
      exponents[STO6G][nShell][sAzimuthal][0] = 3.232838646e00;   coefficients[STO6G][nShell][sAzimuthal][0] = 1.374817488e-3; 
      exponents[STO6G][nShell][sAzimuthal][1] = 3.605788802e-1;   coefficients[STO6G][nShell][sAzimuthal][1] =-8.666390043e-2; 
      exponents[STO6G][nShell][sAzimuthal][2] = 1.717905487e-1;   coefficients[STO6G][nShell][sAzimuthal][2] =-3.130627309e-1; 
      exponents[STO6G][nShell][sAzimuthal][3] = 5.277666487e-2;   coefficients[STO6G][nShell][sAzimuthal][3] = 7.812787397e-1; 
      exponents[STO6G][nShell][sAzimuthal][4] = 3.163400284e-2;   coefficients[STO6G][nShell][sAzimuthal][4] = 4.389247988e-1; 
      exponents[STO6G][nShell][sAzimuthal][5] = 1.874093091e-2;   coefficients[STO6G][nShell][sAzimuthal][5] = 2.487178756e-2; 
      // 4p                                                               
      exponents[STO6G][nShell][pAzimuthal][0] = 2.389722618e00;   coefficients[STO6G][nShell][pAzimuthal][0] =-1.665913575e-3; 
      exponents[STO6G][nShell][pAzimuthal][1] = 7.960947826e-1;   coefficients[STO6G][nShell][pAzimuthal][1] =-1.657464971e-2; 
      exponents[STO6G][nShell][pAzimuthal][2] = 3.415541380e-1;   coefficients[STO6G][nShell][pAzimuthal][2] =-5.958513378e-2; 
      exponents[STO6G][nShell][pAzimuthal][3] = 8.847434525e-2;   coefficients[STO6G][nShell][pAzimuthal][3] = 4.053115554e-1; 
      exponents[STO6G][nShell][pAzimuthal][4] = 4.958248334e-2;   coefficients[STO6G][nShell][pAzimuthal][4] = 5.433958189e-1; 
      exponents[STO6G][nShell][pAzimuthal][5] = 2.816929784e-2;   coefficients[STO6G][nShell][pAzimuthal][5] = 1.204970491e-1; 
      // 4d                                                               
      exponents[STO6G][nShell][dAzimuthal][0] = 4.634239420e00;   coefficients[STO6G][nShell][dAzimuthal][0] =-4.749842876e-4; 
      exponents[STO6G][nShell][dAzimuthal][1] = 1.341648295e00;   coefficients[STO6G][nShell][dAzimuthal][1] =-3.566777891e-3; 
      exponents[STO6G][nShell][dAzimuthal][2] = 2.209593028e-1;   coefficients[STO6G][nShell][dAzimuthal][2] = 1.108670481e-1; 
      exponents[STO6G][nShell][dAzimuthal][3] = 1.101467943e-1;   coefficients[STO6G][nShell][dAzimuthal][3] = 4.159646930e-1; 
      exponents[STO6G][nShell][dAzimuthal][4] = 5.904190370e-2;   coefficients[STO6G][nShell][dAzimuthal][4] = 4.621672517e-1; 
      exponents[STO6G][nShell][dAzimuthal][5] = 3.232628887e-2;   coefficients[STO6G][nShell][dAzimuthal][5] = 1.081250196e-1; 
   }

}

}





