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
#include<stdexcept>
#include<boost/format.hpp>
#include<boost/math/special_functions/factorials.hpp>
#include"Enums.h"
#include"Uncopyable.h"
#include"PrintController.h"
#include"MolDSException.h"
#include"MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
#include"../wrappers/Blas.h"
#include"../wrappers/Lapack.h"
#include"MathUtilities.h"
using namespace std;

namespace MolDS_base{

// n!
int Factorial(int n){
   if(n<0){
      stringstream ss;
      ss << "Error in base::MathUtility::Factorial: n<0 \n";
      throw MolDSException(ss.str());
   }
   return static_cast<int>(boost::math::factorial<double>(n));
}

// nCk
int Conbination(int n, int k){
   if(n < 0){ 
      stringstream ss;
      ss << "Error in base::MathUtility::Conbination: n<0 \n";
      throw MolDSException(ss.str());
   }
   else if(k < 0){ 
      stringstream ss;
      ss << "Error in base::MathUtility::Conbination: k<0 \n";
      throw MolDSException(ss.str());
   }
   else if(n < k){ 
      stringstream ss;
      ss << "Error in base::MathUtility::Conbination: n<k \n";
      throw MolDSException(ss.str());
   }
   else{
      return Factorial(n)/(Factorial(k)*Factorial(n-k));
   }
}

// rotating matrix
void CalcRotatingMatrix(double matrix[][3], double theta, CartesianType cartesianType){
   if(cartesianType == XAxis){
      matrix[0][0] = 1.0;
      matrix[0][1] = 0.0;
      matrix[0][2] = 0.0;

      matrix[1][0] = 0.0;
      matrix[1][1] = cos(theta);
      matrix[1][2] = sin(theta);

      matrix[2][0] = 0.0;
      matrix[2][1] = -sin(theta);
      matrix[2][2] = cos(theta);
   }
   else if(cartesianType == YAxis){
      matrix[0][0] = cos(theta);
      matrix[0][1] = 0.0;
      matrix[0][2] = -sin(theta);

      matrix[1][0] = 0.0;
      matrix[1][1] = 1.0;
      matrix[1][2] = 0.0;

      matrix[2][0] = sin(theta);
      matrix[2][1] = 0.0;
      matrix[2][2] = cos(theta);
   }
   else if(cartesianType == ZAxis){
      matrix[0][0] = cos(theta);
      matrix[0][1] = sin(theta);
      matrix[0][2] = 0.0;

      matrix[1][0] = -sin(theta);
      matrix[1][1] = cos(theta);
      matrix[1][2] = 0.0;

      matrix[2][0] = 0.0;
      matrix[2][1] = 0.0;
      matrix[2][2] = 1.0;
   }
   else{
      stringstream ss;
      ss << "Error in base::MathUtility::CalcRotatingMatrix: invalid cartesianType \n";
      throw MolDSException(ss.str());
   }
}

// calculate determinant of the matrix
double GetDeterminant(double** matrix, int dim){
   double determinant=1.0;
#ifdef __FCC_VERSION
   int*      ipiv=NULL;
   MallocerFreer::GetInstance()->Malloc<int>(&ipiv, dim);
#else
   intptr_t* ipiv=NULL;
   MallocerFreer::GetInstance()->Malloc<intptr_t>(&ipiv, dim);
#endif
   MolDS_wrappers::Lapack::GetInstance()->Dgetrf(matrix, ipiv, dim, dim);
   for(int i=0; i<dim; i++){
      determinant*=matrix[i][i];
      if(ipiv[i] != i-1){
         determinant *= -1.0;
      }
   }
#ifdef __FCC_VERSION
   MallocerFreer::GetInstance()->Free<int>(&ipiv, dim);
#else
   MallocerFreer::GetInstance()->Free<intptr_t>(&ipiv, dim);
#endif
   return determinant;
}

void NormalizeVector(double* vector, int dim){
   double naiseki = MolDS_wrappers::Blas::GetInstance()->Ddot(dim, vector, vector);
   double norm = sqrt(naiseki);
   for(int i=0; i<dim; i++){
      vector[i] /= norm;
   }
}

void OrthoNormalizeSchmidt(double** vectors, int dim, int numVector){
   NormalizeVector(&vectors[0][0], dim);
   for(int i=1; i<numVector; i++){
      for(int j=0; j<i; j++){
         double naiseki = MolDS_wrappers::Blas::GetInstance()->Ddot(dim, &vectors[i][0], &vectors[j][0]);
         MolDS_wrappers::Blas::GetInstance()->Daxpy(dim, -naiseki, &vectors[j][0], &vectors[i][0]);
      }
      NormalizeVector(&vectors[i][0], dim);
   }
}

}
 
