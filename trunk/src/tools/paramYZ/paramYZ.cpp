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
#include<iostream>
#include<math.h>
#include<string>
#include<time.h>
#include<stdexcept>
#include"../../base/Enums.h"
//#undef INCLUDED_ENUMS
//#define RENUMSTR_BODY 1
//#include"../../base/Enums.h"
#include"../../base/MolDSException.h"
#include"../../base/MathUtilities.h"

using namespace std;
using namespace MolDS_base;


/*******************************************************************
This program claculates Z[na-1][nb-1][k] and Y[na-1][nb-1][la][lb][m][i][j].

Z[na][nb][k] corresponds Z_{k\lammda} for given n_{a} and n_{b} 
in (B.28) in Pople book.

Y[na][nb][la][lb][m][i][j] corresponds Y_{ij\lammda} 
for given n_{a}, n_{b}, l_{a}, l_{b}, and m in (B.20) in Pople book.
*********************************************************************/
int main(){
   try{
   int I = 2*ShellType_end+1;
   int J = 2*ShellType_end+1;
   double C[ShellType_end][ShellType_end][ShellType_end];
   double Y[ShellType_end+1][ShellType_end+1][ShellType_end][ShellType_end][ShellType_end][I][J];
   double Z[2*ShellType_end][2*ShellType_end][4*ShellType_end-1];


   for(int i=0;i<ShellType_end;i++){
      for(int j=0;j<ShellType_end;j++){
         for(int k=0;k<ShellType_end;k++){
            C[i][j][k] = 0.0;
         }
      }
   }

   for(int i=0;i<ShellType_end+1;i++){
      for(int j=0;j<ShellType_end+1;j++){
         for(int k=0;k<ShellType_end;k++){
            for(int l=0;l<ShellType_end;l++){
               for(int m=0;m<ShellType_end;m++){
                  for(int n=0;n<I;n++){
                     for(int o=0;o<J;o++){
                        Y[i][j][k][l][m][n][o] = 0.0;
                     }
                  }
               }
            }
         }
      }
   }

   for(int i=0;i<2*ShellType_end;i++){
      for(int j=0;j<2*ShellType_end;j++){
         for(int k=0;k<4*ShellType_end-1;k++){
            Z[i][j][k] = 0.0;
         }
      }
   }


   // set C
   if(4 < ShellType_end){
      cout << "Error in tools/paramYZ/paramYZ.cpp: parameter C is prepared upto N-Shell.\n";
      exit(EXIT_FAILURE);
   }
   C[0][0][0] =  8.0;
   C[1][0][1] =  8.0;
   C[1][1][0] =  4.0;
   C[2][0][0] = -4.0;
   C[2][0][2] = 12.0;
   C[2][1][1] = 12.0;
   C[2][2][0] =  4.0;
   C[3][1][0] = -6.0;
   C[3][1][2] = 30.0;
   C[3][2][1] = 20.0;
   C[3][3][0] =  5.0;


   // calculate Z
   printf("Calculate Z[%d][%d][%d]\n", 2*ShellType_end, 2*ShellType_end, 4*ShellType_end-1);
   double valueZ = 0.0;
   double tempZ = 0.0;
   for(int na=0; na<2*ShellType_end; na++){
      for(int nb=0; nb<2*ShellType_end; nb++){
         for(int k=0; k<4*ShellType_end-1; k++){
            valueZ = 0.0; 
            tempZ = 0.0;

            for(int i=0; i<=na; i++){
               for(int j=0; j<=nb; j++){
                  if(i+j == k){
                     tempZ = pow(-1.0, nb-j);
                     tempZ *= Conbination(na, i);
                     tempZ *= Conbination(nb, j);
                     valueZ += tempZ;
                  }
               }
            }

            Z[na][nb][k] = valueZ;
            printf("\t%lf, \n",valueZ);
         }
      }
   }

   /*
   printf("Z diff check!!!\n");
   for(int na=0; na<2*ShellType_end; na++){
      for(int nb=0; nb<2*ShellType_end; nb++){
         for(int k=0; k<4*ShellType_end-1; k++){
            if(Z[na][nb][k] - pow(-1.0, na+nb-k)*Z[nb][na][k] > 0){
               printf("diff: %lf %lf\n",Z[na][nb][k],pow(-1.0, na+nb-k)*Z[nb][na][k]);
            }
         }
      }
   }
   */
   cout << "\n\n\n";


   // calculate Y
   printf("Calculate Y[%d][%d][%d][%d][%d][%d][%d]\n", ShellType_end+1, ShellType_end+1, ShellType_end, ShellType_end, ShellType_end, I, J);
   double valueY = 0.0;
   double tempY = 0.0;
   for(int na=0; na<=ShellType_end; na++){
      for(int nb=0; nb<=ShellType_end; nb++){
         for(int la=0; la<ShellType_end; la++){
            for(int lb=0; lb<ShellType_end; lb++){
               for(int m=0; m<ShellType_end; m++){
                  for(int i=0; i<I; i++){
                     for(int j=0; j<J; j++){
                        valueY = 0.0;
                        tempY = 0.0;

                        if(0<na && 0<nb){
                        for(int u=0; u<=la-m; u++){
                           for(int v=0; v<=lb-m; v++){
                              if(0<fabs(C[la][m][u]) && 0<fabs(C[lb][m][v])){
                                 for(int a=0; a<=m; a++){
                                    for(int b=0; b<=m; b++){
                                       for(int c=0; c<=u; c++){
                                          for(int d=0; d<=v; d++){
                                             for(int e=0; e<=na-m-u; e++){
                                                for(int f=0; f<=nb-m-v; f++){

                                                   if(i == 2*a+c+d+e+f && j== 2*b+c+d+na+nb-2*m-u-v-e-f){
                                                      tempY = C[la][m][u]*C[lb][m][v];
                                                      tempY *= pow(-1.0, m-a+b+d+nb-m-v-f);
                                                      tempY *= pow(-1.0, v);
                                                      tempY *= Conbination(m, a);
                                                      tempY *= Conbination(m, b);
                                                      tempY *= Conbination(u, c);
                                                      tempY *= Conbination(v, d);
                                                      tempY *= Conbination(na-m-u, e);
                                                      tempY *= Conbination(nb-m-v, f);
                                                      valueY += tempY;
                                                   }

                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                        }

                        Y[na][nb][la][lb][m][i][j] = valueY;
                        printf("\t%lf, \n",valueY);
                        //printf("Y[%d][%d][%d][%d][%d][%d][%d]=%lf\n",na,nb,la,lb,m,i,j,Y[na][nb][la][lb][m][i][j]);
                     }
                  }
               }
            }
         }
      }
   }

  
   /*
   printf("Y diff check!!!\n");
   for(int na=1; na<=ShellType_end; na++){
      for(int nb=1; nb<=ShellType_end; nb++){
         for(int la=0; la<ShellType_end; la++){
            for(int lb=0; lb<ShellType_end; lb++){
               for(int m=0; m<ShellType_end; m++){
                  for(int i=0;i<I;i++){
                     for(int j=0;j<J;j++){

                           if(fabs( Y[na][nb][la][lb][m][i][j] - Y[nb][na][lb][la][m][i][j]*pow(-1.0,j)    )>0){
                              printf("diffY %lf %lf\n",
                              Y[na][nb][la][lb][m][i][j], 
                              Y[nb][na][lb][la][m][i][j]*pow(-1.0,j));
                           }
                     }
                  }

               }
            }
         }
      }
   }
   */
   }
   catch(MolDSException ex){
      cout << ex.what();
   }
   return 0;
}














