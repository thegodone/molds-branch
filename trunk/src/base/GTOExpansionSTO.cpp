//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
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

   //STO-1G, k-shell
   {   
      // 1s
      exponents[STO1G][k][sAzimuthal][0] = 2.709498091e-1;   coefficients[STO1G][k][sAzimuthal][0] = 1.0000;
   }

   //STO-1G, l-shell
   {
      // 2s
      exponents[STO1G][l][sAzimuthal][0] = 1.012151084e-1;   coefficients[STO1G][l][sAzimuthal][0] = 1.0000;
      // 2p
      exponents[STO1G][l][pAzimuthal][0] = 1.759666885e-1;   coefficients[STO1G][l][pAzimuthal][0] = 1.0000;
   }

   //STO-1G, m-shell
   {
      // 3s
      exponents[STO1G][m][sAzimuthal][0] = 5.296881757e-2;   coefficients[STO1G][m][sAzimuthal][0] = 1.0000;
      // 3p
      exponents[STO1G][m][pAzimuthal][0] = 9.113614253e-2;   coefficients[STO1G][m][pAzimuthal][0] = 1.0000;
      // 3d
      exponents[STO1G][m][dAzimuthal][0] = 1.302270363e-1;   coefficients[STO1G][m][dAzimuthal][0] = 1.0000;
   }

   //STO-2G, k-shell
   {
      // 1s
      exponents[STO2G][k][sAzimuthal][0] = 8.518186635e-1;   coefficients[STO2G][k][sAzimuthal][0] = 4.301284983e-1; 
      exponents[STO2G][k][sAzimuthal][1] = 1.516232927e-1;   coefficients[STO2G][k][sAzimuthal][1] = 6.789135305e-1; 
   }   

   //STO-2G, l-shell
   {
      // 2s
      exponents[STO2G][l][sAzimuthal][0] = 1.292278611e-1;   coefficients[STO2G][l][sAzimuthal][0] = 7.470867124e-1; 
      exponents[STO2G][l][sAzimuthal][1] = 4.908584205e-2;   coefficients[STO2G][l][sAzimuthal][1] = 2.855980556e-1; 
      // 2p
      exponents[STO2G][l][pAzimuthal][0] = 4.323908358e-1;   coefficients[STO2G][l][pAzimuthal][0] = 4.522627513e-1; 
      exponents[STO2G][l][pAzimuthal][1] = 1.069439065e-1;   coefficients[STO2G][l][pAzimuthal][1] = 6.713122642e-1; 
   }   

   //STO-2G, m-shell
   {
      // 3s
      exponents[STO2G][m][sAzimuthal][0] = 6.694095822e-1;   coefficients[STO2G][m][sAzimuthal][0] =-1.529645716e-1; 
      exponents[STO2G][m][sAzimuthal][1] = 5.837135094e-2;   coefficients[STO2G][m][sAzimuthal][1] = 1.051370110; 
      // 3p
      exponents[STO2G][m][pAzimuthal][0] = 1.458620964e-1;   coefficients[STO2G][m][pAzimuthal][0] = 5.349653144e-1; 
      exponents[STO2G][m][pAzimuthal][1] = 5.664210742e-2;   coefficients[STO2G][m][pAzimuthal][1] = 5.299607212e-1; 
      // 3d
      exponents[STO2G][m][dAzimuthal][0] = 2.777427345e-1;   coefficients[STO2G][m][dAzimuthal][0] = 4.666137923e-1; 
      exponents[STO2G][m][dAzimuthal][1] = 8.336507714e-2;   coefficients[STO2G][m][dAzimuthal][1] = 6.644706516e-1; 
   }

   //STO-3G, k-shell
   {
      // 1s
      exponents[STO3G][k][sAzimuthal][0] = 2.227660584e00;   coefficients[STO3G][k][sAzimuthal][0] = 1.543289673e-1; 
      exponents[STO3G][k][sAzimuthal][1] = 4.057711562e-1;   coefficients[STO3G][k][sAzimuthal][1] = 5.353281423e-1; 
      exponents[STO3G][k][sAzimuthal][2] = 1.098175104e-1;   coefficients[STO3G][k][sAzimuthal][2] = 4.446345422e-1; 
   }   

   //STO-3G, l-shell
   {
      // 2s
      exponents[STO3G][l][sAzimuthal][0] = 2.581578398e00;   coefficients[STO3G][l][sAzimuthal][0] =-5.994474934e-2; 
      exponents[STO3G][l][sAzimuthal][1] = 1.567622104e-1;   coefficients[STO3G][l][sAzimuthal][1] = 5.960385398e-1; 
      exponents[STO3G][l][sAzimuthal][2] = 6.018332272e-2;   coefficients[STO3G][l][sAzimuthal][2] = 4.581786291e-1; 
      // 2p
      exponents[STO3G][l][pAzimuthal][0] = 9.192379002e-1;   coefficients[STO3G][l][pAzimuthal][0] = 1.623948553e-1; 
      exponents[STO3G][l][pAzimuthal][1] = 2.359194503e-1;   coefficients[STO3G][l][pAzimuthal][1] = 5.661708862e-1; 
      exponents[STO3G][l][pAzimuthal][2] = 8.009805746e-2;   coefficients[STO3G][l][pAzimuthal][2] = 4.223071752e-1; 
   }   

   //STO-3G, m-shell
   {
      // 3s
      exponents[STO3G][m][sAzimuthal][0] = 5.641487709e-1;   coefficients[STO3G][m][sAzimuthal][0] =-1.782577972e-1; 
      exponents[STO3G][m][sAzimuthal][1] = 6.924421391e-2;   coefficients[STO3G][m][sAzimuthal][1] = 8.612761663e-1; 
      exponents[STO3G][m][sAzimuthal][2] = 3.269529097e-2;   coefficients[STO3G][m][sAzimuthal][2] = 2.261841969e-1; 
      // 3p
      exponents[STO3G][m][pAzimuthal][0] = 2.692880368e00;   coefficients[STO3G][m][pAzimuthal][0] =-1.061945788e-2; 
      exponents[STO3G][m][pAzimuthal][1] = 1.489359592e-1;   coefficients[STO3G][m][pAzimuthal][1] = 5.218564264e-1; 
      exponents[STO3G][m][pAzimuthal][2] = 5.739585040e-2;   coefficients[STO3G][m][pAzimuthal][2] = 5.450015143e-1; 
      // 3d
      exponents[STO3G][m][dAzimuthal][0] = 5.229112225e-1;   coefficients[STO3G][m][dAzimuthal][0] = 1.686596060e-1; 
      exponents[STO3G][m][dAzimuthal][1] = 1.639595876e-1;   coefficients[STO3G][m][dAzimuthal][1] = 5.847984817e-1; 
      exponents[STO3G][m][dAzimuthal][2] = 6.386630021e-2;   coefficients[STO3G][m][dAzimuthal][2] = 4.056779523e-1; 
   }

   //STO-4G, k-shell
   {
      // 1s
      exponents[STO4G][k][sAzimuthal][0] = 5.216844534e00;   coefficients[STO4G][k][sAzimuthal][0] = 5.675242080e-2; 
      exponents[STO4G][k][sAzimuthal][1] = 9.546182760e-1;   coefficients[STO4G][k][sAzimuthal][1] = 2.601413550e-1; 
      exponents[STO4G][k][sAzimuthal][2] = 2.652034102e-1;   coefficients[STO4G][k][sAzimuthal][2] = 5.328461143e-1; 
      exponents[STO4G][k][sAzimuthal][3] = 8.801862774e-2;   coefficients[STO4G][k][sAzimuthal][3] = 2.916254405e-1; 
   }   

   //STO-4G, l-shell
   {
      // 2s
      exponents[STO4G][l][sAzimuthal][0] = 1.161525551e01;   coefficients[STO4G][l][sAzimuthal][0] =-1.198411747e-2; 
      exponents[STO4G][l][sAzimuthal][1] = 2.000243111e00;   coefficients[STO4G][l][sAzimuthal][1] =-5.472052539e-2; 
      exponents[STO4G][l][sAzimuthal][2] = 1.607280687e-1;   coefficients[STO4G][l][sAzimuthal][2] = 5.805587176e-1; 
      exponents[STO4G][l][sAzimuthal][3] = 6.125744532e-2;   coefficients[STO4G][l][sAzimuthal][3] = 4.770079976e-1; 
      // 2p
      exponents[STO4G][l][pAzimuthal][0] = 1.798260992e00;   coefficients[STO4G][l][pAzimuthal][0] = 5.713170255e-2; 
      exponents[STO4G][l][pAzimuthal][1] = 4.662622228e-1;   coefficients[STO4G][l][pAzimuthal][1] = 2.857455515e-1; 
      exponents[STO4G][l][pAzimuthal][2] = 1.643718620e-1;   coefficients[STO4G][l][pAzimuthal][2] = 5.517873105e-1; 
      exponents[STO4G][l][pAzimuthal][3] = 6.543927065e-2;   coefficients[STO4G][l][pAzimuthal][3] = 2.632314924e-1; 
   }   

   //STO-4G, m-shell
   {
      // 3s
      exponents[STO4G][m][sAzimuthal][0] = 1.513265591e00;   coefficients[STO4G][m][sAzimuthal][0] =-3.295496352e-2; 
      exponents[STO4G][m][sAzimuthal][1] = 4.262497508e-1;   coefficients[STO4G][m][sAzimuthal][1] =-1.724516959e-1; 
      exponents[STO4G][m][sAzimuthal][2] = 7.643320863e-2;   coefficients[STO4G][m][sAzimuthal][2] = 7.518511194e-1; 
      exponents[STO4G][m][sAzimuthal][3] = 3.760545063e-2;   coefficients[STO4G][m][sAzimuthal][3] = 3.589627317e-1; 
      // 3p
      exponents[STO4G][m][pAzimuthal][0] = 1.853180239e00;   coefficients[STO4G][m][pAzimuthal][0] =-1.434249391e-2; 
      exponents[STO4G][m][pAzimuthal][1] = 1.915075719e-1;   coefficients[STO4G][m][pAzimuthal][1] = 2.755177589e-1; 
      exponents[STO4G][m][pAzimuthal][2] = 8.655487938e-2;   coefficients[STO4G][m][pAzimuthal][2] = 5.846750879e-1; 
      exponents[STO4G][m][pAzimuthal][3] = 4.184253862e-2;   coefficients[STO4G][m][pAzimuthal][3] = 2.144986514e-1; 
      // 3d
      exponents[STO4G][m][dAzimuthal][0] = 9.185846715e-1;   coefficients[STO4G][m][dAzimuthal][0] = 5.799057705e-2; 
      exponents[STO4G][m][dAzimuthal][1] = 2.920461109e-1;   coefficients[STO4G][m][dAzimuthal][1] = 3.045581349e-1; 
      exponents[STO4G][m][dAzimuthal][2] = 1.187568890e-1;   coefficients[STO4G][m][dAzimuthal][2] = 5.601358038e-1; 
      exponents[STO4G][m][dAzimuthal][3] = 5.286755896e-2;   coefficients[STO4G][m][dAzimuthal][3] = 2.432423313e-1; 
   }

   //STO-5G, k-shell
   {
      // 1s
      exponents[STO5G][k][sAzimuthal][0] = 1.130563696e01;   coefficients[STO5G][k][sAzimuthal][0] = 2.214055312e-2; 
      exponents[STO5G][k][sAzimuthal][1] = 2.071728178e00;   coefficients[STO5G][k][sAzimuthal][1] = 1.135411520e-1; 
      exponents[STO5G][k][sAzimuthal][2] = 5.786484833e-1;   coefficients[STO5G][k][sAzimuthal][2] = 3.318161484e-1; 
      exponents[STO5G][k][sAzimuthal][3] = 1.975724573e-1;   coefficients[STO5G][k][sAzimuthal][3] = 4.825700713e-1; 
      exponents[STO5G][k][sAzimuthal][4] = 7.445271746e-2;   coefficients[STO5G][k][sAzimuthal][4] = 1.935721966e-1; 
   }   

   //STO-5G, l-shell
   {
      // 2s
      exponents[STO5G][l][sAzimuthal][0] = 8.984956862e00;   coefficients[STO5G][l][sAzimuthal][0] =-1.596349096e-2; 
      exponents[STO5G][l][sAzimuthal][1] = 1.673710636e00;   coefficients[STO5G][l][sAzimuthal][1] =-5.685884883e-2; 
      exponents[STO5G][l][sAzimuthal][2] = 1.944726668e-1;   coefficients[STO5G][l][sAzimuthal][2] = 3.698265599e-1; 
      exponents[STO5G][l][sAzimuthal][3] = 8.806345634e-2;   coefficients[STO5G][l][sAzimuthal][3] = 5.480512593e-1; 
      exponents[STO5G][l][sAzimuthal][4] = 4.249068522e-2;   coefficients[STO5G][l][sAzimuthal][4] = 1.472634893e-1; 
      // 2p
      exponents[STO5G][l][pAzimuthal][0] = 3.320386533e00;   coefficients[STO5G][l][pAzimuthal][0] = 2.079051117e-2; 
      exponents[STO5G][l][pAzimuthal][1] = 8.643257633e-1;   coefficients[STO5G][l][pAzimuthal][1] = 1.235472099e-1; 
      exponents[STO5G][l][pAzimuthal][2] = 3.079819284e-1;   coefficients[STO5G][l][pAzimuthal][2] = 3.667738986e-1; 
      exponents[STO5G][l][pAzimuthal][3] = 1.273309895e-1;   coefficients[STO5G][l][pAzimuthal][3] = 4.834930290e-1; 
      exponents[STO5G][l][pAzimuthal][4] = 5.606243164e-2;   coefficients[STO5G][l][pAzimuthal][4] = 1.653444074e-1; 
   }   

   //STO-5G, m-shell
   {
      // 3s
      exponents[STO5G][m][sAzimuthal][0] = 4.275877914e00;   coefficients[STO5G][m][sAzimuthal][0] =-3.920358850e-3; 
      exponents[STO5G][m][sAzimuthal][1] = 1.132409433e00;   coefficients[STO5G][m][sAzimuthal][1] =-4.168430506e-2; 
      exponents[STO5G][m][sAzimuthal][2] = 4.016256968e-1;   coefficients[STO5G][m][sAzimuthal][2] =-1.637440990e-1; 
      exponents[STO5G][m][sAzimuthal][3] = 7.732370620e-2;   coefficients[STO5G][m][sAzimuthal][3] = 7.419373723e-1; 
      exponents[STO5G][m][sAzimuthal][4] = 3.800708627e-2;   coefficients[STO5G][m][sAzimuthal][4] = 3.724364929e-1; 
      // 3p
      exponents[STO5G][m][pAzimuthal][0] = 6.466803859e00;   coefficients[STO5G][m][pAzimuthal][0] =-2.329023747e-3; 
      exponents[STO5G][m][pAzimuthal][1] = 1.555914802e00;   coefficients[STO5G][m][pAzimuthal][1] =-1.357395221e-2; 
      exponents[STO5G][m][pAzimuthal][2] = 1.955925255e-1;   coefficients[STO5G][m][pAzimuthal][2] = 2.632185383e-1; 
      exponents[STO5G][m][pAzimuthal][3] = 8.809647701e-2;   coefficients[STO5G][m][pAzimuthal][3] = 5.880427024e-1; 
      exponents[STO5G][m][pAzimuthal][4] = 4.234835707e-2;   coefficients[STO5G][m][pAzimuthal][4] = 2.242794445e-1; 
      // 3d
      exponents[STO5G][m][dAzimuthal][0] = 1.539033958e00;   coefficients[STO5G][m][dAzimuthal][0] = 2.020869128e-2; 
      exponents[STO5G][m][dAzimuthal][1] = 4.922090297e-1;   coefficients[STO5G][m][dAzimuthal][1] = 1.321157923e-1; 
      exponents[STO5G][m][dAzimuthal][2] = 2.029756928e-1;   coefficients[STO5G][m][dAzimuthal][2] = 3.911240346e-1; 
      exponents[STO5G][m][dAzimuthal][3] = 9.424112917e-2;   coefficients[STO5G][m][dAzimuthal][3] = 4.779609701e-1; 
      exponents[STO5G][m][dAzimuthal][4] = 4.569058269e-2;   coefficients[STO5G][m][dAzimuthal][4] = 1.463662294e-1; 
   }

   //STO-6G, k-shell
   {
      // 1s
      exponents[STO6G][k][sAzimuthal][0] = 2.310303149e01;   coefficients[STO6G][k][sAzimuthal][0] = 9.163596280e-3; 
      exponents[STO6G][k][sAzimuthal][1] = 4.235915534e00;   coefficients[STO6G][k][sAzimuthal][1] = 4.936149294e-2; 
      exponents[STO6G][k][sAzimuthal][2] = 1.185056519e00;   coefficients[STO6G][k][sAzimuthal][2] = 1.685383049e-1; 
      exponents[STO6G][k][sAzimuthal][3] = 4.070988982e-1;   coefficients[STO6G][k][sAzimuthal][3] = 3.705627997e-1; 
      exponents[STO6G][k][sAzimuthal][4] = 1.580884151e-1;   coefficients[STO6G][k][sAzimuthal][4] = 4.164915298e-1; 
      exponents[STO6G][k][sAzimuthal][5] = 6.510953954e-2;   coefficients[STO6G][k][sAzimuthal][5] = 1.303340841e-1; 
   }   

   //STO-6G, l-shell
   {
      // 2s
      exponents[STO6G][l][sAzimuthal][0] = 2.768496241e01;   coefficients[STO6G][l][sAzimuthal][0] =-4.151277819e-3; 
      exponents[STO6G][l][sAzimuthal][1] = 5.077140627e00;   coefficients[STO6G][l][sAzimuthal][1] =-2.067024148e-2; 
      exponents[STO6G][l][sAzimuthal][2] = 1.426786050e00;   coefficients[STO6G][l][sAzimuthal][2] =-5.150303337e-2; 
      exponents[STO6G][l][sAzimuthal][3] = 2.040335729e-1;   coefficients[STO6G][l][sAzimuthal][3] = 3.346271174e-1; 
      exponents[STO6G][l][sAzimuthal][4] = 9.260298399e-2;   coefficients[STO6G][l][sAzimuthal][4] = 5.621061301e-1; 
      exponents[STO6G][l][sAzimuthal][5] = 4.416183978e-2;   coefficients[STO6G][l][sAzimuthal][5] = 1.712994697e-1; 
      // 2p
      exponents[STO6G][l][pAzimuthal][0] = 5.868285913e00;   coefficients[STO6G][l][pAzimuthal][0] = 7.924233646e-3; 
      exponents[STO6G][l][pAzimuthal][1] = 1.530329631e00;   coefficients[STO6G][l][pAzimuthal][1] = 5.144104825e-2; 
      exponents[STO6G][l][pAzimuthal][2] = 5.475665231e-1;   coefficients[STO6G][l][pAzimuthal][2] = 1.898400060e-1; 
      exponents[STO6G][l][pAzimuthal][3] = 2.288932733e-1;   coefficients[STO6G][l][pAzimuthal][3] = 4.049863191e-1; 
      exponents[STO6G][l][pAzimuthal][4] = 1.046655969e-1;   coefficients[STO6G][l][pAzimuthal][4] = 4.012362861e-1; 
      exponents[STO6G][l][pAzimuthal][5] = 4.948220127e-2;   coefficients[STO6G][l][pAzimuthal][5] = 1.051855189e-1; 
   }   

   //STO-6G, m-shell
   {
      // 3s
      exponents[STO6G][m][sAzimuthal][0] = 3.273031938e00;   coefficients[STO6G][m][sAzimuthal][0] =-6.775596947e-3; 
      exponents[STO6G][m][sAzimuthal][1] = 9.200611311e-1;   coefficients[STO6G][m][sAzimuthal][1] =-5.639325779e-2; 
      exponents[STO6G][m][sAzimuthal][2] = 3.593349765e-1;   coefficients[STO6G][m][sAzimuthal][2] =-1.587656086e-1; 
      exponents[STO6G][m][sAzimuthal][3] = 8.636686991e-2;   coefficients[STO6G][m][sAzimuthal][3] = 5.534527651e-1; 
      exponents[STO6G][m][sAzimuthal][4] = 4.797373812e-2;   coefficients[STO6G][m][sAzimuthal][4] = 5.015351020e-1; 
      exponents[STO6G][m][sAzimuthal][5] = 2.724741144e-2;   coefficients[STO6G][m][sAzimuthal][5] = 7.223633674e-2; 
      // 3p
      exponents[STO6G][m][pAzimuthal][0] = 5.077973607e00;   coefficients[STO6G][m][pAzimuthal][0] =-3.329929840e-3; 
      exponents[STO6G][m][pAzimuthal][1] = 1.340786940e00;   coefficients[STO6G][m][pAzimuthal][1] =-1.419488340e-2; 
      exponents[STO6G][m][pAzimuthal][2] = 2.248434849e-1;   coefficients[STO6G][m][pAzimuthal][2] = 1.639395770e-1; 
      exponents[STO6G][m][pAzimuthal][3] = 1.131741848e-1;   coefficients[STO6G][m][pAzimuthal][3] = 4.485358256e-1; 
      exponents[STO6G][m][pAzimuthal][4] = 6.076408893e-2;   coefficients[STO6G][m][pAzimuthal][4] = 3.908813050e-1; 
      exponents[STO6G][m][pAzimuthal][5] = 3.315424265e-2;   coefficients[STO6G][m][pAzimuthal][5] = 7.411456232e-2; 
      // 3d
      exponents[STO6G][m][dAzimuthal][0] = 2.488296923e00;   coefficients[STO6G][m][dAzimuthal][0] = 7.283828112e-3; 
      exponents[STO6G][m][dAzimuthal][1] = 7.981487853e-1;   coefficients[STO6G][m][dAzimuthal][1] = 5.386799363e-2; 
      exponents[STO6G][m][dAzimuthal][2] = 3.311327490e-1;   coefficients[STO6G][m][dAzimuthal][2] = 2.072139149e-1; 
      exponents[STO6G][m][dAzimuthal][3] = 1.559114463e-1;   coefficients[STO6G][m][dAzimuthal][3] = 4.266269092e-1; 
      exponents[STO6G][m][dAzimuthal][4] = 7.877734732e-2;   coefficients[STO6G][m][dAzimuthal][4] = 3.843100204e-1; 
      exponents[STO6G][m][dAzimuthal][5] = 4.058484363e-2;   coefficients[STO6G][m][dAzimuthal][5] = 8.902827546e-2; 
   }

}

}





