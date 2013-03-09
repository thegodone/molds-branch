//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
// Copyright (C) 2012-2012 Katsuhiko Nishimra                             // 
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
#ifdef __INTEL_COMPILER
#include"mkl.h"
#else
#include"cblas.h"
#endif
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/Uncopyable.h"
#include"../base/MallocerFreer.h"
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
void Blas::Dcopy(molds_blas_int n,
                 double const* vectorX,
                 double *      vectorY)const{
   molds_blas_int incrementX=1;
   molds_blas_int incrementY=1;
   this->Dcopy(n, vectorX, incrementX, vectorY, incrementY);
}

// vectorY = vectorX
//    vectorX: n-vector
//    vectorY: n-vector
void Blas::Dcopy(molds_blas_int n,
                 double const* vectorX, molds_blas_int incrementX,
                 double*       vectorY, molds_blas_int incrementY) const{
   double* x = const_cast<double*>(&vectorX[0]);
   cblas_dcopy(n, x, incrementX, vectorY, incrementY);
}

// vectorY = alpha*vectorX + vectorY
//    vectorX: n-vector
//    vectorY: n-vector
void Blas::Daxpy(molds_blas_int n, double alpha,
           double const* vectorX,
           double*       vectorY) const{
   molds_blas_int incrementX=1;
   molds_blas_int incrementY=1;
   this->Daxpy(n, alpha, vectorX, incrementX, vectorY, incrementY);
}

// vectorY = alpha*vectorX + vectorY
//    vectorX: n-vector
//    vectorY: n-vector
void Blas::Daxpy(molds_blas_int n, double alpha,
           double const* vectorX, molds_blas_int incrementX,
           double*       vectorY, molds_blas_int incrementY) const{
   double* x = const_cast<double*>(&vectorX[0]);
   cblas_daxpy(n, alpha, x, incrementX, vectorY, incrementY);
}

// returns vectorX^T*vectorY
//    vectorX: n-vector
//    vectorY: n-vector
double Blas::Ddot(molds_blas_int n,
            double const* vectorX,
            double const* vectorY) const{
   molds_blas_int incrementX=1;
   molds_blas_int incrementY=1;
   return this->Ddot(n, vectorX, incrementX, vectorY, incrementY);
}

// returns vectorX^T*vectorY
//    vectorX: n-vector
//    vectorY: n-vector
double Blas::Ddot(molds_blas_int n,
            double const* vectorX, molds_blas_int incrementX,
            double const* vectorY, molds_blas_int incrementY)const{
   double* x=const_cast<double*>(vectorX),
         * y=const_cast<double*>(vectorY);
   return cblas_ddot(n, x, incrementX, y, incrementY);
}

// vectorY = matrixA*vectorX
//    matrixA: m*n-matrix (matrixA[m][n] in row-major (C/C++ style))
//    vectorX: n-vector
//    vectorY: m-vector
void Blas::Dgemv(molds_blas_int m, molds_blas_int n,
                 double const* const* matrixA,
                 double const* vectorX,
                 double*       vectorY) const{
   bool isColumnMajorMatrixA = false; // because, in general, C/C++ style is row-major.
   molds_blas_int incrementX=1;
   molds_blas_int incrementY=1;
   double alpha  =1.0; 
   double beta   =0.0; 
   this->Dgemv(isColumnMajorMatrixA, m, n, alpha, matrixA, vectorX, incrementX, beta, vectorY, incrementY);
}

// vectorY = alpha*matrixA*vectorX + beta*vectorY
//    matrixA: m*n-matrix
//    vectorX: n-vector
//    vectorY: m-vector
void Blas::Dgemv(bool isColumnMajorMatrixA,
                 molds_blas_int m, molds_blas_int n,
                 double alpha,
                 double const* const* matrixA,
                 double const* vectorX ,
                 molds_blas_int incrementX,
                 double beta,
                 double* vectorY,
                 molds_blas_int incrementY) const{
   double* a = const_cast<double*>(&matrixA[0][0]);
   double* x = const_cast<double*>(&vectorX[0]);
   CBLAS_TRANSPOSE transA;
   if(isColumnMajorMatrixA){
      transA = CblasNoTrans;
   }
   else{
      transA = CblasTrans;
      swap(m,n);
   }
   int lda = m;
   cblas_dgemv(CblasColMajor, transA, m, n, alpha, a, lda, x, incrementX, beta, vectorY, incrementY);
}

// vectorY = matrixA*vectorX
//    matrixA: n*n-matrix,symmetric (Use the upper triangular part)
//    vectorX: n-vector
//    vectorY: n-vector
void Blas::Dsymv(molds_blas_int n,
           double const* const* matrixA,
           double const* vectorX,
           double*       vectorY) const{
   bool isColumnMajorMatrixA = false; // because, in general, C/C++ style is row-major.
   molds_blas_int incrementX=1, incrementY=1;
   double alpha=1.0, beta=0.0;
   this->Dsymv(n, alpha, matrixA, vectorX, incrementX, beta, vectorY, incrementY);
}

// vectorY = alpha*matrixA*vectorX + beta*vectorY
//    matrixA: n*n-matrix,symmetric (Use the upper triangular part)
//    vectorX: n-vector
//    vectorY: n-vector
void Blas::Dsymv(molds_blas_int n, double alpha,
           double const* const* matrixA,
           double const* vectorX, molds_blas_int incrementX,
           double beta,
           double*       vectorY, molds_blas_int incrementY) const{
   double* a = const_cast<double*>(&matrixA[0][0]);
   double* x = const_cast<double*>(&vectorX[0]);
   CBLAS_UPLO uploA=CblasUpper;
   int lda = n;
   cblas_dsymv(CblasRowMajor, uploA, n, alpha, a, lda, x, incrementX, beta, vectorY, incrementY);
}

// matrixA = alpha*vectorX*vectorX^T + matrixA
//    matrixA: n*n-matrix,symmetric (Use the upper triangular part, and copy it to the lower part.)
//    vectorX: n-matrix
void Blas::Dsyr(molds_blas_int n, double alpha,
          double const* vectorX,
          double ** matrixA)const{
   molds_blas_int incrementX=1;
   this->Dsyr(n, alpha, vectorX, incrementX, matrixA);
}

void Blas::Dsyr(molds_blas_int n, double alpha,
          double const* vectorX, molds_blas_int incrementX,
          double ** matrixA)const{
   double* a = &matrixA[0][0];
   double* x = const_cast<double*>(&vectorX[0]);
   CBLAS_UPLO uploA=CblasUpper;
   int lda = n;
   cblas_dsyr(CblasRowMajor, uploA, n, alpha, x, incrementX, a, lda);
#pragma omp parallel for schedule(auto)
   for(molds_blas_int i=0;i<n;i++){
      for(molds_blas_int j=i+1;j<n;j++){
         matrixA[j][i] = matrixA[i][j];
      }
   }
}

// matrixC = matrixA*matrixB
//    matrixA: m*k-matrix (matrixA[m][k] in row-major (C/C++ style))
//    matrixB: k*n-matrix (matrixB[k][n] in row-major (C/C++ style))
//    matrixC: m*n-matrix (matrixC[m][n] in row-major (C/C++ style))
void Blas::Dgemm(molds_blas_int m, molds_blas_int n, molds_blas_int k, 
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
                 molds_blas_int m, molds_blas_int n, molds_blas_int k,  
                 double alpha,
                 double const* const* matrixA,
                 double const* const* matrixB,
                 double beta,
                 double** matrixC) const{
   double* a = const_cast<double*>(&matrixA[0][0]);
   double* b = const_cast<double*>(&matrixB[0][0]);
   double*       c = &matrixC[0][0];

   int lda;
   CBLAS_TRANSPOSE transA;
   if(isColumnMajorMatrixA){
      transA = CblasNoTrans;
      lda = m;
   }
   else{
      transA = CblasTrans;
      lda = k;
   }

   int ldb;
   CBLAS_TRANSPOSE transB;
   if(isColumnMajorMatrixB){
      transB = CblasNoTrans;
      ldb = k;
   }
   else{
      transB = CblasTrans;
      ldb = n;
   }

   double* tmpC;
#ifdef __INTEL_COMPILER
   tmpC = (double*)mkl_malloc( sizeof(double)*m*n, 16 );
#else
   tmpC = (double*)malloc( sizeof(double)*m*n);
#endif
   for(molds_blas_int i=0; i<m; i++){
      for(molds_blas_int j=0; j<n; j++){
         tmpC[i+j*m] = matrixC[i][j];
      }
   }
   molds_blas_int ldc = m;
   //call blas
   cblas_dgemm(CblasColMajor, transA, transB, m, n, k, alpha, a, lda, b, ldb, beta, tmpC, ldc);
   for(molds_blas_int i=0; i<m; i++){
      for(molds_blas_int j=0; j<n; j++){
         matrixC[i][j] = tmpC[i+j*m];
      }
   }
#ifdef __INTEL_COMPILER
   mkl_free(tmpC);
#else
   free(tmpC);
#endif
}

// matrixD = matrixA*matrixB*matrixC
//    matrixA: m*k-matrix (matrixA[m][k] in row-major (C/C++ style))
//    matrixB: k*l-matrix (matrixB[k][n] in row-major (C/C++ style))
//    matrixC: k*n-matrix (matrixC[k][n] in row-major (C/C++ style))
//    matrixD: m*n-matrix (matrixC[m][n] in row-major (C/C++ style))
void Blas::Dgemmm(molds_blas_int m, molds_blas_int n, molds_blas_int k, molds_blas_int l,
                  double const* const* matrixA, 
                  double const* const* matrixB, 
                  double const* const* matrixC, 
                  double**             matrixD) const{
   bool isColumnMajorMatrixA = false; // because, in general, C/C++ style is row-major.
   bool isColumnMajorMatrixB = false; // because, in general, C/C++ style is row-major.
   bool isColumnMajorMatrixC = false; // because, in general, C/C++ style is row-major.
   double alpha=1.0;
   double beta =0.0;
   this->Dgemmm(isColumnMajorMatrixA, isColumnMajorMatrixB, isColumnMajorMatrixC, m, n, k, l, alpha, matrixA, matrixB, matrixC, beta, matrixD);
}

// matrixD = alpha*matrixA*matrixB*matrixC + beta*matrixD
//    matrixA: m*k-matrix 
//    matrixB: k*l-matrix
//    matrixC: l*n-matrix
//    matrixD: m*n-matrix (matrixC[m][n] in row-major (C/C++ style))
void Blas::Dgemmm(bool isColumnMajorMatrixA,
                  bool isColumnMajorMatrixB, 
                  bool isColumnMajorMatrixC, 
                  molds_blas_int m, molds_blas_int n, molds_blas_int k, molds_blas_int l,
                  double alpha,
                  double const* const* matrixA,
                  double const* const* matrixB,
                  double const* const* matrixC,
                  double beta,
                  double** matrixD) const{
   
   double** matrixBC = NULL;
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&matrixBC, k, n); 
      this->Dgemmm(isColumnMajorMatrixA, isColumnMajorMatrixB, isColumnMajorMatrixC,
                   m, n, k, l, 
                   alpha,
                   matrixA,
                   matrixB,
                   matrixC,
                   beta,
                   matrixD,
                   matrixBC);
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&matrixBC, k, n); 
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&matrixBC, k, n); 
}

// matrixD = alpha*matrixA*matrixB*matrixC + beta*matrixD
//    matrixA: m*k-matrix 
//    matrixB: k*l-matrix
//    matrixC: l*n-matrix
//    matrixD: m*n-matrix (matrixC[m][n] in row-major (C/C++ style))
//       tmpMatrixBC is temporary calculated matrix in row-major, (C/C++ style) 
//       tmpMatrixBC = matrixB*matrixC
void Blas::Dgemmm(bool isColumnMajorMatrixA,
                  bool isColumnMajorMatrixB, 
                  bool isColumnMajorMatrixC, 
                  molds_blas_int m, molds_blas_int n, molds_blas_int k, molds_blas_int l,
                  double alpha,
                  double const* const* matrixA,
                  double const* const* matrixB,
                  double const* const* matrixC,
                  double beta,
                  double** matrixD,
                  double** tmpMatrixBC) const{
   
   double alphaBC = 1.0;
   double betaBC  = 0.0;
   bool isColumnMajorMatrixBC = false;
   this->Dgemm(isColumnMajorMatrixB, isColumnMajorMatrixC,  k, n, l, alphaBC, matrixB, matrixC,     betaBC, tmpMatrixBC);
   this->Dgemm(isColumnMajorMatrixA, isColumnMajorMatrixBC, m, n, k, alpha,   matrixA, tmpMatrixBC, beta,   matrixD );
}

// matrixC = matrixA*matrixA^T
//    matrixA: n*k-matrix
//    matrixC: n*n-matrix,symmetric (Use the upper triangular part, and copy it to the lower part.)
void Blas::Dsyrk(molds_blas_int n, molds_blas_int k,
                 double const* const* matrixA,
                 double ** matrixC)const{
   bool isMatrixAColumnMajor = false;
   bool isMatrixATransposed = false;
   bool isLowerTriangularPartMatrixCUsed = false;
   double alpha = 1.0 , beta = 0.0;
   this->Dsyrk(n, k, isMatrixAColumnMajor, isMatrixATransposed, isLowerTriangularPartMatrixCUsed, alpha, matrixA, beta, matrixC);
}

// matrixC = alpha*matrixA*matrixA^T + beta*matrixC (isMatrixATransposed==false)
//  or
// matrixC = alpha*matrixA^T*matrixA + beta*matrixC (isMatrixATransposed==true)
//    matrixA: n*k-matrix (isMatrixATransposed==false) or k*n-matrix (isMatrixATransposed==true)
//    matrixC: n*n-matrix,symmetric (Use the upper triangular part, and copy it to the lower part.)
void Blas::Dsyrk(molds_blas_int n, molds_blas_int k,
                 bool isMatrixAColumnMajor,
                 bool isMatrixATransposed,
                 bool isLowerTriangularPartMatrixCUsed,
                 double alpha, double const* const* matrixA,
                 double beta,  double ** matrixC)const{
   double* c = &matrixC[0][0];
   double* a = const_cast<double*>(&matrixA[0][0]);
   CBLAS_ORDER orderA = isMatrixAColumnMajor ? CblasColMajor : CblasRowMajor;
   CBLAS_UPLO uploC= isLowerTriangularPartMatrixCUsed ? CblasLower : CblasUpper;
   CBLAS_TRANSPOSE transA= isMatrixATransposed ? CblasTrans : CblasNoTrans;
   molds_blas_int lda = &matrixA[1][0] - &matrixA[0][0];
   molds_blas_int ldc = &matrixC[1][0] - &matrixC[0][0];
   cblas_dsyrk(orderA, uploC, transA, n, k, alpha, a, lda, beta, c, ldc);
#pragma omp parallel for schedule(auto)
   for(molds_blas_int i=0;i<n;i++){
      for(molds_blas_int j=i+1;j<n;j++){
         if(isLowerTriangularPartMatrixCUsed){
            matrixC[i][j] = matrixC[j][i];
         }
         else{
            matrixC[j][i] = matrixC[i][j];
         }
      }
   }
}

}
