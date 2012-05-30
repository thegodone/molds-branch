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
#include<time.h>
#include<list>
#include<vector>
#include<stdexcept>
using namespace std;

/*******************************************************************
 * This program calculates reduced parameters for NDDO-series (MNDO, AM1, PM3, and etc).
 * See p20 & 21 in [MOPAC_1990] for implemeneted procedures and notations.
 * Note that iterative equations in p20 & 21 in [MOPAC_1990] are wrong.
 * The correct iterative equations are shown in Mopac HP:
 *    Therory > Semiempirical theory > NDDO two-electron two-center integrals
 *    in http://openmopac.net/manual/index.html
 *
 * Methods calculating D1 and D2 should be selected according to the periodic table.
 * That is, see commented outed code to calculate D1 and D2 in this file.
 *
 * Note that eV2AU should be equal to MolDS_base::Parameters::eV2AU.
 *
 * refferences
 * [MOPAC_1990] J. J. P. Stewart, J. Computer-Aided Molecular Design 4, 1 (1990)
*********************************************************************/

long double GetAForAD(long double AD, long double D1){
   long double a=0.0;
   a = 0.5*AD
      -0.5/sqrt(4.0*pow(D1,2.0)+pow(AD,-2.0));
   return a;
}

long double GetAForAQ(long double AQ, long double D2){
   long double a=0.0;
   a = 0.25*AQ
      -0.50/sqrt(4.0*pow(D2,2.0)+pow(AQ,-2.0)) 
      +0.25/sqrt(8.0*pow(D2,2.0)+pow(AQ,-2.0)) ;
   return a;
}

int main(){
   // notation is [MOPAC1970]
   // all valuable should be in atomic units.
   long double D1=0.0;
   long double D2=0.0;
   long double AM=0.0;
   long double AD=0.0;
   long double AQ=0.0;
   long double AD_old=0.0;
   long double AD_old2=0.0;
   long double AQ_old=0.0;
   long double AQ_old2=0.0;

   long double orbitalExponentS=1.891185;
   long double orbitalExponentP=1.658972;

   double eV2AU = 0.03674903;
   long double Gss = 8.964667 * eV2AU;
   long double Gpp = 9.968164 * eV2AU;
   long double Gsp = 6.785936 * eV2AU;
   long double Gpp2= 7.970247 * eV2AU;
   long double Hsp = 4.041836 * eV2AU;
   long double Hpp = 0.5*(Gpp - Gpp2);

   // output prepared parameters
   printf("=====  NDDO parameters =====\n");
   printf("orbital exponent S in [a.u.] = %.10lf\n",(double)orbitalExponentS);
   printf("orbital exponent P in [a.u.] = %.10lf\n",(double)orbitalExponentP);
   printf("Gss in [a.u.] = %.10lf\n",(double)Gss);
   printf("Gpp in [a.u.] = %.10lf\n",(double)Gpp);
   printf("Gsp in [a.u.] = %.10lf\n",(double)Gsp);
   printf("Gpp2 in [a.u.] = %.10lf\n",(double)Gpp2);
   printf("Hsp in [a.u.] = %.10lf\n",(double)Hsp);
   printf("Hpp = 0.5*(Gpp - Gpp2)\n\n\n");

   // calculateion and output derived parameters
   printf("=====  NDDO derived parameters =====\n");

   /*
   // Calc. D1 for n=2 (C, N, O, and etc.)
   D1 = 5.0*pow(3.0,-0.5)
       *pow(4.0*orbitalExponentS*orbitalExponentP,2.5)
       /pow(orbitalExponentS+orbitalExponentP,6.0);
   printf("D1 in [a.u.] = %.10lf\n",(double)D1);
   */

   // Calc. D1 for n=3 (S and etc.)
   D1 = 7.0*pow(3.0,-0.5)
       *pow(4.0*orbitalExponentS*orbitalExponentP,3.5)
       /pow(orbitalExponentS+orbitalExponentP,8.0);
   printf("D1 in [a.u.] = %.10lf\n",(double)D1);

   /*
   // Calc. D2 for n=2  (C, N, O, and etc.)
   D2 = pow(1.5,0.5)/orbitalExponentP;
   printf("D2 in [a.u.] = %.10lf\n",(double)D2);
   */

   // Calc. D2 for n=3  (S and etc.)
   D2 = pow(2.8,0.5)/orbitalExponentP;
   printf("D2 in [a.u.] = %.10lf\n",(double)D2);

   // Calc. AM
   AM = Gss;
   printf("AM in [a.u.] = %.10lf\n\n",(double)AM);

   // Calc. AD
   AD_old2 = 1.0;
   AD_old  = pow(Hsp/pow(D1,2.0),1.0/3.0);
   for(int n=0;n<10;n++){
      long double a_old2=GetAForAD(AD_old2,D1);
      long double a_old =GetAForAD(AD_old, D1);
      AD = AD_old2 + (AD_old - AD_old2)*(Hsp-a_old2)/(a_old - a_old2);
      AD_old2 = AD_old;
      AD_old = AD;
      printf("iter=%d\tAD in [a.u.] = %.10lf\n",n,(double)AD);
   }
   cout <<endl;

   // Calc. AQ
   AQ_old2 = 1.0;
   AQ_old = pow(Hpp/(3.0*pow(D2,2.0)),1.0/5.0);
   for(int n=0;n<15;n++){
      long double a_old2=GetAForAQ(AQ_old2,D2);
      long double a_old =GetAForAQ(AQ_old, D2);
      AQ = AQ_old2 + (AQ_old - AQ_old2)*(Hpp-a_old2)/(a_old - a_old2);
      AQ_old2 = AQ_old;
      AQ_old = AQ;
      printf("iter=%d\tAQ in [a.u.] = %.10lf\n",n,(double)AQ);
   }
   return 0;
}














