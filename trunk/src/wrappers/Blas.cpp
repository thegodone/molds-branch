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
#include<boost/format.hpp>
#include"mkl.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/Uncopyable.h"
#include"Blas.h"
using namespace std;
using namespace MolDS_base;

namespace MolDS_wrappers{
Blas* Blas::blas = NULL;

Blas::Blas(){
}

Blas::~Blas(){
}

Blas* Blas::GetInstance(){
   if(blas == NULL){
      blas = new Blas();
      //this->OutputLog("Blas created.\n\n");
   }
   return blas;
}

void Blas::DeleteInstance(){
   if(blas != NULL){
      delete blas;
      //this->OutputLog("Blas deleted\n\n");
   }
   blas = NULL;
}

// vectorY = vectorX
//    vectorX: n-vector
//    vectorY: n-vector
void Blas::Dcopy(int n,
                 double const* vectorX,
                 double *      vectorY)const{
   int incrementX=1;
   int incrementY=1;
   this->Dcopy(n, vectorX, incrementX, vectorY, incrementY);
}

// vectorY = vectorX
//    vectorX: n-vector
//    vectorY: n-vector
void Blas::Dcopy(int n,
                 double const* vectorX, int incrementX,
                 double*       vectorY, int incrementY) const{
   dcopy(&n, vectorX, &incrementX, vectorY, &incrementY);
}

// vectorY = alpha*vectorX + vectorY
//    vectorX: n-vector
//    vectorY: n-vector
void Blas::Daxpy(int n, double alpha,
           double const* vectorX,
           double*       vectorY) const{
   int incrementX=1;
   int incrementY=1;
   this->Daxpy(n, alpha, vectorX, incrementX, vectorY, incrementY);
}

// vectorY = alpha*vectorX + vectorY
//    vectorX: n-vector
//    vectorY: n-vector
void Blas::Daxpy(int n, double alpha,
           double const* vectorX, int incrementX,
           double*       vectorY, int incrementY) const{
   daxpy(&n, &alpha, vectorX, &incrementX, vectorY, &incrementY);
}

// vectorY = matrixA*vectorX
//    matrixA: m*n-matrix (matrixA[m][n] in row-major (C/C++ style))
//    vectorX: n-vector
//    vectorY: m-vector
void Blas::Dgemv(int m, int n,
                 double const* const* matrixA,
                 double const* vectorX,
                 double*       vectorY) const{
   bool isColumnMajorMatrixA = false; // because, in general, C/C++ style is row-major.
   int incrementX=1;
   int incrementY=1;
   double alpha  =1.0; 
   double beta   =0.0; 
   this->Dgemv(isColumnMajorMatrixA, m, n, alpha, matrixA, vectorX, incrementX, beta, vectorY, incrementY);
}

// vectorY = alpha*matrixA*vectorX + beta*vectorY
//    matrixA: m*n-matrix
//    vectorX: n-vector
//    vectorY: m-vector
void Blas::Dgemv(bool isColumnMajorMatrixA,
                 int m, int n,
                 double alpha,
                 double const* const* matrixA,
                 double const* vectorX ,
                 int incrementX,
                 double beta,
                 double* vectorY,
                 int incrementY) const{
   double const* a = &matrixA[0][0];
   char transA;
   if(isColumnMajorMatrixA){
      transA = 'N';
   }
   else{
      transA = 'T';
      swap(m,n);
   }
   int lda = m;
   dgemv(&transA, &m, &n, &alpha, a, &lda, vectorX, &incrementX, &beta, vectorY, &incrementY);
}

// matrixC = matrixA*matrixB
//    matrixA: m*k-matrix (matrixA[m][k] in row-major (C/C++ style))
//    matrixB: k*n-matrix (matrixB[k][n] in row-major (C/C++ style))
//    matrixC: m*n-matrix (matrixC[m][n] in row-major (C/C++ style))
void Blas::Dgemm(int m, int n, int k, 
                 double const* const* matrixA, 
                 double const* const* matrixB, 
                 double**             matrixC) const{
   bool isColumnMajorMatrixA = false; // because, in general, C/C++ style is row-major.
   bool isColumnMajorMatrixB = false; // because, in general, C/C++ style is row-major.
   double alpha=1.0;
   double beta =0.0;
   this->Dgemm(isColumnMajorMatrixA, isColumnMajorMatrixB, m, n, k, alpha, matrixA, matrixB, beta, matrixC);
}

// matrixC = alpha*matrixA*matrixB + beta*matrixC
//    matrixA: m*k-matrix 
//    matrixB: k*n-matrix
//    matrixC: m*n-matrix (matrixC[m][n] in row-major (C/C++ style))
void Blas::Dgemm(bool isColumnMajorMatrixA, 
                 bool isColumnMajorMatrixB, 
                 int m, int n, int k,  
                 double alpha,
                 double const* const* matrixA,
                 double const* const* matrixB,
                 double beta,
                 double** matrixC) const{
   double const* a = &matrixA[0][0];
   double const* b = &matrixB[0][0];
   double*       c = &matrixC[0][0];

   char transA;
   int lda;
   if(isColumnMajorMatrixA){
      transA = 'N'; //ka=k
      lda = m;
   }
   else{
      transA = 'T'; //ka=m
      lda = k;
   }

   char transB;
   int ldb;
   if(isColumnMajorMatrixB){
      transB = 'N'; //kb=n
      ldb = k;
   }
   else{
      transB = 'T'; //kb=k
      ldb = n;
   }

   double* tmpC;
   tmpC = (double*)mkl_malloc( sizeof(double)*m*n, 16 );
   for(int i=0; i<m; i++){
      for(int j=0; j<n; j++){
         tmpC[i+j*m] = matrixC[i][j];
      }
   }
   int ldc = m;
   //call blas
   dgemm(&transA, &transB, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, tmpC, &ldc);
   for(int i=0; i<m; i++){
      for(int j=0; j<n; j++){
         matrixC[i][j] = tmpC[i+j*m];
      }
   }
   mkl_free(tmpC);
}

}
