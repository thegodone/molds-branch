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
#include<boost/format.hpp>
#include"../base/Enums.h"
#include"../base/Uncopyable.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
#include"Lapack.h"

#ifdef __INTEL_COMPILER
   #include"mkl.h"
   #include"mkl_lapacke.h"
#elif defined __FCC_VERSION
   #include"lapacke.h"
#else
   #if ( __WORDSIZE == 32 )
   #else
      #define HAVE_LAPACK_CONFIG_H
      #define LAPACK_ILP64
   #endif
   #include"lapacke.h"
#endif

#ifdef __INTEL_COMPILER
#define MOLDS_LAPACK_malloc(a,b) mkl_malloc(a,b)
#define MOLDS_LAPACK_free(a) mkl_free(a)
#else
#define MOLDS_LAPACK_malloc(a,b) malloc(a)
#define MOLDS_LAPACK_free(a) free(a)
#endif


using namespace std;
using namespace MolDS_base;

namespace MolDS_wrappers{
Lapack* Lapack::lapack = NULL;

Lapack::Lapack(){
   this->errorMessageDsyevdInfo = "Error in wrappers::Lapack::Dsyevd: info != 0: info = ";
   this->errorMessageDsyevdSize = "Error in wrappers::Lapack::Dsyevd: size of matirx < 1\n";
   this->errorMessageDsysvInfo = "Error in wrappers::Lapack::Dsysv: info != 0: info = ";
   this->errorMessageDsysvSize = "Error in wrappers::Lapack::Dsysv: size of matirx < 1\n";
   this->errorMessageDgetrsInfo = "Error in wrappers::Lapack::Dgetrs: info != 0: info = ";
   this->errorMessageDgetrsSize = "Error in wrappers::Lapack::Dgetrs: size of matirx < 1\n";
   this->errorMessageDgetrfInfo = "Error in wrappers::Lapack::Dgetrf: info != 0: info = ";
}

Lapack::~Lapack(){
}

Lapack* Lapack::GetInstance(){
   if(lapack == NULL){
      lapack = new Lapack();
      //this->OutputLog("Lapack created.\n\n");
   }
   return lapack;
}

void Lapack::DeleteInstance(){
   if(lapack != NULL){
      delete lapack;
      //this->OutputLog("Lapack deleted\n\n");
   }
   lapack = NULL;
}


/***
 *
 * Eigenvalue and eigenvector of a real symmetirc matrix are calculated:
 *   - i-th eigenvalue will be stored in eigenValues[i].
 *   - i-th eigenvector will be stored as (matirx[i][0], matirx[i][1], matirx[i][2], ....).
 *
 * ***/
molds_lapack_int Lapack::Dsyevd(double** matrix, double* eigenValues, molds_lapack_int size, bool calcEigenVectors){
   molds_lapack_int info = 0;
   char uplo = 'U';
   molds_lapack_int lda = size;
   // set job type
   char job;
   if(calcEigenVectors){
      job = 'V';
   }
   else{
      job = 'N';
   }

   // call Lapack
   info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, job, uplo, size, &matrix[0][0], lda, eigenValues);

   // make i-th row i-the eigenvector
   double** tmpMatrix=NULL;
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&tmpMatrix, size, size);
      for(molds_lapack_int i = 0; i < size; i++){
         for(molds_lapack_int j = 0; j < size; j++){
            tmpMatrix[j][i] = matrix[i][j];
         }
      }
      for(molds_lapack_int i = 0; i < size; i++){
         for(molds_lapack_int j = 0; j < size; j++){
            matrix[i][j] = tmpMatrix[i][j];
         }
      }
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&tmpMatrix, size, size);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&tmpMatrix, size, size);

   // adjust phase of eigenvectors
   for(molds_lapack_int i=0;i<size;i++){
      double tmp = 0.0;
      for(molds_lapack_int j=0;j<size;j++){
         tmp += matrix[i][j];
      }
      if(tmp<0){
         for(molds_lapack_int j=0;j<size;j++){
            matrix[i][j]*=-1.0;
         }
      }
   }   

   if(info != 0){
      stringstream ss;
      ss << errorMessageDsyevdInfo;
      ss << info << endl;
      MolDSException ex(ss.str());
      ex.SetKeyValue<int>(LapackInfo, info);
      throw ex;
   }
   return info;
}

/***
 *
 * "matrix*X=b" is solved, then we get X by this method.
 * The X will be stored in b.
 * The matrix will be overwriten by this method.
 *
 */
molds_lapack_int Lapack::Dsysv(double** matrix, double* b, molds_lapack_int size){
   molds_lapack_int info = 0;
   char uplo             = 'U';
   molds_lapack_int nrhs = 1;
   molds_lapack_int lda  = size;
   molds_lapack_int ldb  = nrhs;
   molds_lapack_int* ipiv;

   if(size < 1 ){
      stringstream ss;
      ss << errorMessageDsysvSize;
      throw MolDSException(ss.str());
   }

   // malloc
   ipiv = (molds_lapack_int*)MOLDS_LAPACK_malloc( sizeof(molds_lapack_int)*2*size, 16 );

   // call Lapack
   info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, uplo, size, nrhs, &matrix[0][0], lda, ipiv, b, ldb);

   // free
   MOLDS_LAPACK_free(ipiv);
  
   if(info != 0){
      stringstream ss;
      ss << errorMessageDsysvInfo;
      ss << info << endl;
      MolDSException ex(ss.str());
      ex.SetKeyValue<int>(LapackInfo, info);
      throw ex;
   }
   return info;
}

/***
 *
 * "matrix*X[i]=b[i] (i=0, 1, ... , nrhs-1) is solved, then we get X[i] by this method.
 * The X[i] will be stored in b[i], namely
 * the b[i][j] will be j-th element of i-th solution, b[i].
 * Besides, the matrix will be overwriten by this method.
 *
 */
molds_lapack_int Lapack::Dgetrs(double** matrix, double** b, molds_lapack_int size, molds_lapack_int nrhs) const{
   molds_lapack_int info = 0;
   char trans = 'N';
   molds_lapack_int lda = size;
   molds_lapack_int ldb = nrhs;
   double* tmpB;
   molds_lapack_int* ipiv;

   if(size < 1 ){
      stringstream ss;
      ss << errorMessageDgetrsSize;
      throw MolDSException(ss.str());
   }

   try{
      // malloc
      ipiv = (molds_lapack_int*)MOLDS_LAPACK_malloc( sizeof(molds_lapack_int)*2*size, 16 );
      tmpB = (double*)MOLDS_LAPACK_malloc( sizeof(double)*size*nrhs, 16 );
      // matrix b should be transposed
      for(molds_lapack_int i = 0; i < nrhs; i++){
         for(molds_lapack_int j = 0; j < size; j++){
            tmpB[j*nrhs+i] = b[i][j];
         }
      }
      this->Dgetrf(&matrix[0][0], ipiv, size, size);
      // call Lapack
      info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, trans, size, nrhs, &matrix[0][0], lda, ipiv, tmpB, ldb);
      for(molds_lapack_int i = 0; i < nrhs; i++){
         for(molds_lapack_int j = 0; j < size; j++){
            b[i][j] = tmpB[j*nrhs+i];
         }
      }
   }
   catch(MolDSException ex){
      // free
      MOLDS_LAPACK_free(tmpB);
      MOLDS_LAPACK_free(ipiv);
      throw ex;
   }
   // free
   MOLDS_LAPACK_free(ipiv);
   MOLDS_LAPACK_free(tmpB);
  
   if(info != 0){
      stringstream ss;
      ss << errorMessageDgetrsInfo;
      ss << info << endl;
      throw MolDSException(ss.str());
   }
   return info;
}

// Argument "matrix" is sizeM*sizeN matrix.
// Argument "matrix" will be LU-decomposed.
molds_lapack_int Lapack::Dgetrf(double** matrix, molds_lapack_int sizeM, molds_lapack_int sizeN) const{
   molds_lapack_int* ipiv = (molds_lapack_int*)MOLDS_LAPACK_malloc( sizeof(molds_lapack_int)*2*sizeM,        16 );
   this->Dgetrf(&matrix[0][0], ipiv, sizeM, sizeN);
   MOLDS_LAPACK_free(ipiv);
   molds_lapack_int info = 0;
   return info;
}

// Argument "matrix" is sizeM*sizeN matrix in Row-major (C/C++ style)
// Argument "matrix" will be LU-decomposed.
molds_lapack_int Lapack::Dgetrf(double** matrix, molds_lapack_int* ipiv, molds_lapack_int sizeM, molds_lapack_int sizeN) const{
   this->Dgetrf(&matrix[0][0], ipiv, sizeM, sizeN);
   molds_lapack_int info = 0;
   return info;
}

// Argument "matrix" is sizeM*sizeN matrix.
// The each element of "matrix" should be stored in 1-dimensional vecotre with Row major (C/C++ style).
molds_lapack_int Lapack::Dgetrf(double* matrix, molds_lapack_int* ipiv, molds_lapack_int sizeM, molds_lapack_int sizeN) const{
   molds_lapack_int info = 0;
   molds_lapack_int lda = sizeM;
   // call Lapack
   info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, sizeM, sizeN, matrix, lda, ipiv);
   if(info != 0){
      stringstream ss;
      ss << errorMessageDgetrfInfo;
      ss << info << endl;
      MolDSException ex(ss.str());
      ex.SetKeyValue<int>(LapackInfo, info);
      throw ex;
   }
   return info;
}
}
