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
#include"MolDSException.h"
#include"Uncopyable.h"
#include"Enums.h"
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
      exponents[STO1G][k][sAzimuthal][0] = 2.709498091*pow(10.0,-1.0);   coefficients[STO1G][k][sAzimuthal][0] = 1.0000;
   }

   //STO-1G, l-shell
   {
      // 2s
      exponents[STO1G][l][sAzimuthal][0] = 1.012151084*pow(10.0,-1.0);   coefficients[STO1G][l][sAzimuthal][0] = 1.0000;
      // 2p
      exponents[STO1G][l][pAzimuthal][0] = 1.759666885*pow(10.0,-1.0);   coefficients[STO1G][l][pAzimuthal][0] = 1.0000;
   }

   //STO-1G, m-shell
   {
      // 3s
      exponents[STO1G][m][sAzimuthal][0] = 5.296881757*pow(10.0,-2.0);   coefficients[STO1G][m][sAzimuthal][0] = 1.0000;
      // 3p
      exponents[STO1G][m][pAzimuthal][0] = 9.113614253*pow(10.0,-2.0);   coefficients[STO1G][m][pAzimuthal][0] = 1.0000;
      // 3d
      exponents[STO1G][m][dAzimuthal][0] = 1.302270363*pow(10.0,-1.0);   coefficients[STO1G][m][dAzimuthal][0] = 1.0000;
   }

   //STO-2G, k-shell
   {
      // 1s
      exponents[STO2G][k][sAzimuthal][0] = 8.518186635*pow(10.0,-1.0);   coefficients[STO2G][k][sAzimuthal][0] = 4.301284983*pow(10,-1.0); 
      exponents[STO2G][k][sAzimuthal][1] = 1.516232927*pow(10.0,-1.0);   coefficients[STO2G][k][sAzimuthal][1] = 6.789135305*pow(10,-1.0); 
   }   

   //STO-2G, l-shell
   {
      // 2s
      exponents[STO2G][l][sAzimuthal][0] = 1.292278611*pow(10.0,-1.0);   coefficients[STO2G][l][sAzimuthal][0] = 7.470867124*pow(10,-1.0); 
      exponents[STO2G][l][sAzimuthal][1] = 4.908584205*pow(10.0,-2.0);   coefficients[STO2G][l][sAzimuthal][1] = 2.855980556*pow(10,-1.0); 
      // 2p
      exponents[STO2G][l][pAzimuthal][0] = 4.323908358*pow(10.0,-1.0);   coefficients[STO2G][l][pAzimuthal][0] = 4.522627513*pow(10,-1.0); 
      exponents[STO2G][l][pAzimuthal][1] = 1.069439065*pow(10.0,-1.0);   coefficients[STO2G][l][pAzimuthal][1] = 6.713122642*pow(10,-1.0); 
   }   

   //STO-2G, m-shell
   {
      // 3s
      exponents[STO2G][m][sAzimuthal][0] = 6.694095822*pow(10.0,-1.0);   coefficients[STO2G][m][sAzimuthal][0] =-1.529645716*pow(10,-1.0); 
      exponents[STO2G][m][sAzimuthal][1] = 5.837135094*pow(10.0,-2.0);   coefficients[STO2G][m][sAzimuthal][1] = 1.051370110; 
      // 3p
      exponents[STO2G][m][pAzimuthal][0] = 1.458620964*pow(10.0,-1.0);   coefficients[STO2G][m][pAzimuthal][0] = 5.349653144*pow(10,-1.0); 
      exponents[STO2G][m][pAzimuthal][1] = 5.664210742*pow(10.0,-2.0);   coefficients[STO2G][m][pAzimuthal][1] = 5.299607212*pow(10,-1.0); 
      // 3d
      exponents[STO2G][m][dAzimuthal][0] = 2.777427345*pow(10.0,-1.0);   coefficients[STO2G][m][dAzimuthal][0] = 4.666137923*pow(10,-1.0); 
      exponents[STO2G][m][dAzimuthal][1] = 8.336507714*pow(10.0,-2.0);   coefficients[STO2G][m][dAzimuthal][1] = 6.644706516*pow(10,-1.0); 
   }

   //STO-3G, k-shell
   {
      // 1s
      exponents[STO3G][k][sAzimuthal][0] = 2.227660584;                  coefficients[STO3G][k][sAzimuthal][0] = 1.543289673*pow(10,-1.0); 
      exponents[STO3G][k][sAzimuthal][1] = 4.057711562*pow(10.0,-1.0);   coefficients[STO3G][k][sAzimuthal][1] = 5.353281423*pow(10,-1.0); 
      exponents[STO3G][k][sAzimuthal][2] = 1.098175104*pow(10.0,-1.0);   coefficients[STO3G][k][sAzimuthal][2] = 4.446345422*pow(10,-1.0); 
   }   

   //STO-3G, l-shell
   {
      // 2s
      exponents[STO3G][l][sAzimuthal][0] = 2.581578398;                  coefficients[STO3G][l][sAzimuthal][0] =-5.994474934*pow(10,-2.0); 
      exponents[STO3G][l][sAzimuthal][1] = 1.567622104*pow(10.0,-1.0);   coefficients[STO3G][l][sAzimuthal][1] = 5.960385398*pow(10,-1.0); 
      exponents[STO3G][l][sAzimuthal][2] = 6.018332272*pow(10.0,-2.0);   coefficients[STO3G][l][sAzimuthal][2] = 4.581786291*pow(10,-1.0); 
      // 2p
      exponents[STO3G][l][pAzimuthal][0] = 9.192379002*pow(10.0,-1.0);   coefficients[STO3G][l][pAzimuthal][0] = 1.623948553*pow(10,-1.0); 
      exponents[STO3G][l][pAzimuthal][1] = 2.359194503*pow(10.0,-1.0);   coefficients[STO3G][l][pAzimuthal][1] = 5.661708862*pow(10,-1.0); 
      exponents[STO3G][l][pAzimuthal][2] = 8.009805746*pow(10.0,-2.0);   coefficients[STO3G][l][pAzimuthal][2] = 4.223071752*pow(10,-1.0); 
   }   

   //STO-3G, m-shell
   {
      // 3s
      exponents[STO3G][m][sAzimuthal][0] = 5.641487709*pow(10.0,-1.0);   coefficients[STO3G][m][sAzimuthal][0] =-1.782577972*pow(10,-1.0); 
      exponents[STO3G][m][sAzimuthal][1] = 6.924421391*pow(10.0,-2.0);   coefficients[STO3G][m][sAzimuthal][1] = 8.612761663*pow(10,-1.0); 
      exponents[STO3G][m][sAzimuthal][2] = 3.269529097*pow(10.0,-2.0);   coefficients[STO3G][m][sAzimuthal][2] = 2.261841969*pow(10,-1.0); 
      // 3p
      exponents[STO3G][m][pAzimuthal][0] = 2.692880368;                  coefficients[STO3G][m][pAzimuthal][0] =-1.061945788*pow(10,-2.0); 
      exponents[STO3G][m][pAzimuthal][1] = 1.489359592*pow(10.0,-1.0);   coefficients[STO3G][m][pAzimuthal][1] = 5.218564264*pow(10,-1.0); 
      exponents[STO3G][m][pAzimuthal][2] = 5.739585040*pow(10.0,-2.0);   coefficients[STO3G][m][pAzimuthal][2] = 5.450015143*pow(10,-1.0); 
      // 3d
      exponents[STO3G][m][dAzimuthal][0] = 5.229112225*pow(10.0,-1.0);   coefficients[STO3G][m][dAzimuthal][0] = 1.686596060*pow(10,-1.0); 
      exponents[STO3G][m][dAzimuthal][1] = 1.639595876*pow(10.0,-1.0);   coefficients[STO3G][m][dAzimuthal][1] = 5.847984817*pow(10,-1.0); 
      exponents[STO3G][m][dAzimuthal][2] = 6.386630021*pow(10.0,-2.0);   coefficients[STO3G][m][dAzimuthal][2] = 4.056779523*pow(10,-1.0); 
   }

   //STO-4G, k-shell
   {
      // 1s
      exponents[STO4G][k][sAzimuthal][0] = 5.216844534;                  coefficients[STO4G][k][sAzimuthal][0] = 5.675242080*pow(10,-2.0); 
      exponents[STO4G][k][sAzimuthal][1] = 9.546182760*pow(10.0,-1.0);   coefficients[STO4G][k][sAzimuthal][1] = 2.601413550*pow(10,-1.0); 
      exponents[STO4G][k][sAzimuthal][2] = 2.652034102*pow(10.0,-1.0);   coefficients[STO4G][k][sAzimuthal][2] = 5.328461143*pow(10,-1.0); 
      exponents[STO4G][k][sAzimuthal][3] = 8.801862774*pow(10.0,-2.0);   coefficients[STO4G][k][sAzimuthal][3] = 2.916254405*pow(10,-1.0); 
   }   

   //STO-4G, l-shell
   {
      // 2s
      exponents[STO4G][l][sAzimuthal][0] = 1.161525551*pow(10.0, 1.0);   coefficients[STO4G][l][sAzimuthal][0] =-1.198411747*pow(10,-2.0); 
      exponents[STO4G][l][sAzimuthal][1] = 2.000243111;                  coefficients[STO4G][l][sAzimuthal][1] =-5.472052539*pow(10,-2.0); 
      exponents[STO4G][l][sAzimuthal][2] = 1.607280687*pow(10.0,-1.0);   coefficients[STO4G][l][sAzimuthal][2] = 5.805587176*pow(10,-1.0); 
      exponents[STO4G][l][sAzimuthal][3] = 6.125744532*pow(10.0,-2.0);   coefficients[STO4G][l][sAzimuthal][3] = 4.770079976*pow(10,-1.0); 
      // 2p
      exponents[STO4G][l][pAzimuthal][0] = 1.798260992;                  coefficients[STO4G][l][pAzimuthal][0] = 5.713170255*pow(10,-2.0); 
      exponents[STO4G][l][pAzimuthal][1] = 4.662622228*pow(10.0,-1.0);   coefficients[STO4G][l][pAzimuthal][1] = 2.857455515*pow(10,-1.0); 
      exponents[STO4G][l][pAzimuthal][2] = 1.643718620*pow(10.0,-1.0);   coefficients[STO4G][l][pAzimuthal][2] = 5.517873105*pow(10,-1.0); 
      exponents[STO4G][l][pAzimuthal][3] = 6.543927065*pow(10.0,-2.0);   coefficients[STO4G][l][pAzimuthal][3] = 2.632314924*pow(10,-1.0); 
   }   

   //STO-4G, m-shell
   {
      // 3s
      exponents[STO4G][m][sAzimuthal][0] = 1.513265591;                  coefficients[STO4G][m][sAzimuthal][0] =-3.295496352*pow(10,-2.0); 
      exponents[STO4G][m][sAzimuthal][1] = 4.262497508*pow(10.0,-1.0);   coefficients[STO4G][m][sAzimuthal][1] =-1.724516959*pow(10,-1.0); 
      exponents[STO4G][m][sAzimuthal][2] = 7.643320863*pow(10.0,-2.0);   coefficients[STO4G][m][sAzimuthal][2] = 7.518511194*pow(10,-1.0); 
      exponents[STO4G][m][sAzimuthal][3] = 3.760545063*pow(10.0,-2.0);   coefficients[STO4G][m][sAzimuthal][3] = 3.589627317*pow(10,-1.0); 
      // 3p
      exponents[STO4G][m][pAzimuthal][0] = 1.853180239;                  coefficients[STO4G][m][pAzimuthal][0] =-1.434249391*pow(10,-2.0); 
      exponents[STO4G][m][pAzimuthal][1] = 1.915075719*pow(10.0,-1.0);   coefficients[STO4G][m][pAzimuthal][1] = 2.755177589*pow(10,-1.0); 
      exponents[STO4G][m][pAzimuthal][2] = 8.655487938*pow(10.0,-2.0);   coefficients[STO4G][m][pAzimuthal][2] = 5.846750879*pow(10,-1.0); 
      exponents[STO4G][m][pAzimuthal][3] = 4.184253862*pow(10.0,-2.0);   coefficients[STO4G][m][pAzimuthal][3] = 2.144986514*pow(10,-1.0); 
      // 3d
      exponents[STO4G][m][dAzimuthal][0] = 9.185846715*pow(10.0,-1.0);   coefficients[STO4G][m][dAzimuthal][0] = 5.799057705*pow(10,-2.0); 
      exponents[STO4G][m][dAzimuthal][1] = 2.920461109*pow(10.0,-1.0);   coefficients[STO4G][m][dAzimuthal][1] = 3.045581349*pow(10,-1.0); 
      exponents[STO4G][m][dAzimuthal][2] = 1.187568890*pow(10.0,-1.0);   coefficients[STO4G][m][dAzimuthal][2] = 5.601358038*pow(10,-1.0); 
      exponents[STO4G][m][dAzimuthal][3] = 5.286755896*pow(10.0,-2.0);   coefficients[STO4G][m][dAzimuthal][3] = 2.432423313*pow(10,-1.0); 
   }

   //STO-5G, k-shell
   {
      // 1s
      exponents[STO5G][k][sAzimuthal][0] = 1.130563696*pow(10.0, 1.0);   coefficients[STO5G][k][sAzimuthal][0] = 2.214055312*pow(10,-2.0); 
      exponents[STO5G][k][sAzimuthal][1] = 2.071728178;                  coefficients[STO5G][k][sAzimuthal][1] = 1.135411520*pow(10,-1.0); 
      exponents[STO5G][k][sAzimuthal][2] = 5.786484833*pow(10.0,-1.0);   coefficients[STO5G][k][sAzimuthal][2] = 3.318161484*pow(10,-1.0); 
      exponents[STO5G][k][sAzimuthal][3] = 1.975724573*pow(10.0,-1.0);   coefficients[STO5G][k][sAzimuthal][3] = 4.825700713*pow(10,-1.0); 
      exponents[STO5G][k][sAzimuthal][4] = 7.445271746*pow(10.0,-2.0);   coefficients[STO5G][k][sAzimuthal][4] = 1.935721966*pow(10,-1.0); 
   }   

   //STO-5G, l-shell
   {
      // 2s
      exponents[STO5G][l][sAzimuthal][0] = 8.984956862;                  coefficients[STO5G][l][sAzimuthal][0] =-1.596349096*pow(10,-2.0); 
      exponents[STO5G][l][sAzimuthal][1] = 1.673710636;                  coefficients[STO5G][l][sAzimuthal][1] =-5.685884883*pow(10,-2.0); 
      exponents[STO5G][l][sAzimuthal][2] = 1.944726668*pow(10.0,-1.0);   coefficients[STO5G][l][sAzimuthal][2] = 3.698265599*pow(10,-1.0); 
      exponents[STO5G][l][sAzimuthal][3] = 8.806345634*pow(10.0,-2.0);   coefficients[STO5G][l][sAzimuthal][3] = 5.480512593*pow(10,-1.0); 
      exponents[STO5G][l][sAzimuthal][4] = 4.249068522*pow(10.0,-2.0);   coefficients[STO5G][l][sAzimuthal][4] = 1.472634893*pow(10,-1.0); 
      // 2p
      exponents[STO5G][l][pAzimuthal][0] = 3.320386533;                  coefficients[STO5G][l][pAzimuthal][0] = 2.079051117*pow(10,-2.0); 
      exponents[STO5G][l][pAzimuthal][1] = 8.643257633*pow(10.0,-1.0);   coefficients[STO5G][l][pAzimuthal][1] = 1.235472099*pow(10,-1.0); 
      exponents[STO5G][l][pAzimuthal][2] = 3.079819284*pow(10.0,-1.0);   coefficients[STO5G][l][pAzimuthal][2] = 3.667738986*pow(10,-1.0); 
      exponents[STO5G][l][pAzimuthal][3] = 1.273309895*pow(10.0,-1.0);   coefficients[STO5G][l][pAzimuthal][3] = 4.834930290*pow(10,-1.0); 
      exponents[STO5G][l][pAzimuthal][4] = 5.606243164*pow(10.0,-2.0);   coefficients[STO5G][l][pAzimuthal][4] = 1.653444074*pow(10,-1.0); 
   }   

   //STO-5G, m-shell
   {
      // 3s
      exponents[STO5G][m][sAzimuthal][0] = 4.275877914;                  coefficients[STO5G][m][sAzimuthal][0] =-3.920358850*pow(10,-3.0); 
      exponents[STO5G][m][sAzimuthal][1] = 1.132409433;                  coefficients[STO5G][m][sAzimuthal][1] =-4.168430506*pow(10,-2.0); 
      exponents[STO5G][m][sAzimuthal][2] = 4.016256968*pow(10.0,-1.0);   coefficients[STO5G][m][sAzimuthal][2] =-1.637440990*pow(10,-1.0); 
      exponents[STO5G][m][sAzimuthal][3] = 7.732370620*pow(10.0,-2.0);   coefficients[STO5G][m][sAzimuthal][3] = 7.419373723*pow(10,-1.0); 
      exponents[STO5G][m][sAzimuthal][4] = 3.800708627*pow(10.0,-2.0);   coefficients[STO5G][m][sAzimuthal][4] = 3.724364929*pow(10,-1.0); 
      // 3p
      exponents[STO5G][m][pAzimuthal][0] = 6.466803859;                  coefficients[STO5G][m][pAzimuthal][0] =-2.329023747*pow(10,-3.0); 
      exponents[STO5G][m][pAzimuthal][1] = 1.555914802;                  coefficients[STO5G][m][pAzimuthal][1] =-1.357395221*pow(10,-2.0); 
      exponents[STO5G][m][pAzimuthal][2] = 1.955925255*pow(10.0,-1.0);   coefficients[STO5G][m][pAzimuthal][2] = 2.632185383*pow(10,-1.0); 
      exponents[STO5G][m][pAzimuthal][3] = 8.809647701*pow(10.0,-2.0);   coefficients[STO5G][m][pAzimuthal][3] = 5.880427024*pow(10,-1.0); 
      exponents[STO5G][m][pAzimuthal][4] = 4.234835707*pow(10.0,-2.0);   coefficients[STO5G][m][pAzimuthal][4] = 2.242794445*pow(10,-1.0); 
      // 3d
      exponents[STO5G][m][dAzimuthal][0] = 1.539033958;                  coefficients[STO5G][m][dAzimuthal][0] = 2.020869128*pow(10,-2.0); 
      exponents[STO5G][m][dAzimuthal][1] = 4.922090297*pow(10.0,-1.0);   coefficients[STO5G][m][dAzimuthal][1] = 1.321157923*pow(10,-1.0); 
      exponents[STO5G][m][dAzimuthal][2] = 2.029756928*pow(10.0,-1.0);   coefficients[STO5G][m][dAzimuthal][2] = 3.911240346*pow(10,-1.0); 
      exponents[STO5G][m][dAzimuthal][3] = 9.424112917*pow(10.0,-2.0);   coefficients[STO5G][m][dAzimuthal][3] = 4.779609701*pow(10,-1.0); 
      exponents[STO5G][m][dAzimuthal][4] = 4.569058269*pow(10.0,-2.0);   coefficients[STO5G][m][dAzimuthal][4] = 1.463662294*pow(10,-1.0); 
   }

   //STO-6G, k-shell
   {
      // 1s
      exponents[STO6G][k][sAzimuthal][0] = 2.310303149*pow(10.0, 1.0);   coefficients[STO6G][k][sAzimuthal][0] = 9.163596280*pow(10,-3.0); 
      exponents[STO6G][k][sAzimuthal][1] = 4.235915534;                  coefficients[STO6G][k][sAzimuthal][1] = 4.936149294*pow(10,-2.0); 
      exponents[STO6G][k][sAzimuthal][2] = 1.185056519;                  coefficients[STO6G][k][sAzimuthal][2] = 1.685383049*pow(10,-1.0); 
      exponents[STO6G][k][sAzimuthal][3] = 4.070988982*pow(10.0,-1.0);   coefficients[STO6G][k][sAzimuthal][3] = 3.705627997*pow(10,-1.0); 
      exponents[STO6G][k][sAzimuthal][4] = 1.580884151*pow(10.0,-1.0);   coefficients[STO6G][k][sAzimuthal][4] = 4.164915298*pow(10,-1.0); 
      exponents[STO6G][k][sAzimuthal][5] = 6.510953954*pow(10.0,-2.0);   coefficients[STO6G][k][sAzimuthal][5] = 1.303340841*pow(10,-1.0); 
   }   

   //STO-6G, l-shell
   {
      // 2s
      exponents[STO6G][l][sAzimuthal][0] = 2.768496241*pow(10.0, 1.0);   coefficients[STO6G][l][sAzimuthal][0] =-4.151277819*pow(10,-3.0); 
      exponents[STO6G][l][sAzimuthal][1] = 5.077140627;                  coefficients[STO6G][l][sAzimuthal][1] =-2.067024148*pow(10,-2.0); 
      exponents[STO6G][l][sAzimuthal][2] = 1.426786050;                  coefficients[STO6G][l][sAzimuthal][2] =-5.150303337*pow(10,-2.0); 
      exponents[STO6G][l][sAzimuthal][3] = 2.040335729*pow(10.0,-1.0);   coefficients[STO6G][l][sAzimuthal][3] = 3.346271174*pow(10,-1.0); 
      exponents[STO6G][l][sAzimuthal][4] = 9.260298399*pow(10.0,-2.0);   coefficients[STO6G][l][sAzimuthal][4] = 5.621061301*pow(10,-1.0); 
      exponents[STO6G][l][sAzimuthal][5] = 4.416183978*pow(10.0,-2.0);   coefficients[STO6G][l][sAzimuthal][5] = 1.712994697*pow(10,-1.0); 
      // 2p
      exponents[STO6G][l][pAzimuthal][0] = 5.868285913;                  coefficients[STO6G][l][pAzimuthal][0] = 7.924233646*pow(10,-3.0); 
      exponents[STO6G][l][pAzimuthal][1] = 1.530329631;                  coefficients[STO6G][l][pAzimuthal][1] = 5.144104825*pow(10,-2.0); 
      exponents[STO6G][l][pAzimuthal][2] = 5.475665231*pow(10.0,-1.0);   coefficients[STO6G][l][pAzimuthal][2] = 1.898400060*pow(10,-1.0); 
      exponents[STO6G][l][pAzimuthal][3] = 2.288932733*pow(10.0,-1.0);   coefficients[STO6G][l][pAzimuthal][3] = 4.049863191*pow(10,-1.0); 
      exponents[STO6G][l][pAzimuthal][4] = 1.046655969*pow(10.0,-1.0);   coefficients[STO6G][l][pAzimuthal][4] = 4.012362861*pow(10,-1.0); 
      exponents[STO6G][l][pAzimuthal][5] = 4.948220127*pow(10.0,-2.0);   coefficients[STO6G][l][pAzimuthal][5] = 1.051855189*pow(10,-1.0); 
   }   

   //STO-6G, m-shell
   {
      // 3s
      exponents[STO6G][m][sAzimuthal][0] = 3.273031938;                  coefficients[STO6G][m][sAzimuthal][0] =-6.775596947*pow(10,-3.0); 
      exponents[STO6G][m][sAzimuthal][1] = 9.200611311*pow(10.0,-1.0);   coefficients[STO6G][m][sAzimuthal][1] =-5.639325779*pow(10,-2.0); 
      exponents[STO6G][m][sAzimuthal][2] = 3.593349765*pow(10.0,-1.0);   coefficients[STO6G][m][sAzimuthal][2] =-1.587656086*pow(10,-1.0); 
      exponents[STO6G][m][sAzimuthal][3] = 8.636686991*pow(10.0,-2.0);   coefficients[STO6G][m][sAzimuthal][3] = 5.534527651*pow(10,-1.0); 
      exponents[STO6G][m][sAzimuthal][4] = 4.797373812*pow(10.0,-2.0);   coefficients[STO6G][m][sAzimuthal][4] = 5.015351020*pow(10,-1.0); 
      exponents[STO6G][m][sAzimuthal][5] = 2.724741144*pow(10.0,-2.0);   coefficients[STO6G][m][sAzimuthal][5] = 7.223633674*pow(10,-2.0); 
      // 3p
      exponents[STO6G][m][pAzimuthal][0] = 5.077973607;                  coefficients[STO6G][m][pAzimuthal][0] =-3.329929840*pow(10,-3.0); 
      exponents[STO6G][m][pAzimuthal][1] = 1.340786940;                  coefficients[STO6G][m][pAzimuthal][1] =-1.419488340*pow(10,-2.0); 
      exponents[STO6G][m][pAzimuthal][2] = 2.248434849*pow(10.0,-1.0);   coefficients[STO6G][m][pAzimuthal][2] = 1.639395770*pow(10,-1.0); 
      exponents[STO6G][m][pAzimuthal][3] = 1.131741848*pow(10.0,-1.0);   coefficients[STO6G][m][pAzimuthal][3] = 4.485358256*pow(10,-1.0); 
      exponents[STO6G][m][pAzimuthal][4] = 6.076408893*pow(10.0,-2.0);   coefficients[STO6G][m][pAzimuthal][4] = 3.908813050*pow(10,-1.0); 
      exponents[STO6G][m][pAzimuthal][5] = 3.315424265*pow(10.0,-2.0);   coefficients[STO6G][m][pAzimuthal][5] = 7.411456232*pow(10,-2.0); 
      // 3d
      exponents[STO6G][m][dAzimuthal][0] = 2.488296923;                  coefficients[STO6G][m][dAzimuthal][0] = 7.283828112*pow(10,-3.0); 
      exponents[STO6G][m][dAzimuthal][1] = 7.981487853*pow(10.0,-1.0);   coefficients[STO6G][m][dAzimuthal][1] = 5.386799363*pow(10,-2.0); 
      exponents[STO6G][m][dAzimuthal][2] = 3.311327490*pow(10.0,-1.0);   coefficients[STO6G][m][dAzimuthal][2] = 2.072139149*pow(10,-1.0); 
      exponents[STO6G][m][dAzimuthal][3] = 1.559114463*pow(10.0,-1.0);   coefficients[STO6G][m][dAzimuthal][3] = 4.266269092*pow(10,-1.0); 
      exponents[STO6G][m][dAzimuthal][4] = 7.877734732*pow(10.0,-2.0);   coefficients[STO6G][m][dAzimuthal][4] = 3.843100204*pow(10,-1.0); 
      exponents[STO6G][m][dAzimuthal][5] = 4.058484363*pow(10.0,-2.0);   coefficients[STO6G][m][dAzimuthal][5] = 8.902827546*pow(10,-2.0); 
   }

}

}





