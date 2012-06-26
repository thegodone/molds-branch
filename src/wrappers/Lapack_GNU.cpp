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
#include"lapacke.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/Uncopyable.h"
#include"Lapack.h"
using namespace std;
using namespace MolDS_base;

namespace MolDS_wrappers{
Lapack* Lapack::lapack = NULL;

Lapack::Lapack(){
   this->calculatedDsysvBlockSize = false;
   this->dsysvBlockSize = 64;
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
 * return notice
 *    i-th eigen value is eigenValues[i].
 *    i-th eigen vector is (matirx[i][0], matirx[i][1], matirx[i][2], ....).
 *
 * ***/
int Lapack::Dsyevd(double** matrix, double* eigenValues, int size, bool calcEigenVectors){
   int info = 0;
   int k = 0;
   int lwork;
   int liwork;
   char job;
   char uplo = 'U';
   int lda = size;
   double* convertedMatrix;
   double* tempEigenValues;
   double* work;
   int* iwork;

   // set job type
   if(calcEigenVectors){
      job = 'V';
   }
   else{
      job = 'N';
   }

   // calc. lwork and liwork
   if(size < 1 ){
      stringstream ss;
      ss << errorMessageDsyevdSize;
      throw MolDSException(ss.str());
   }
   else if(size == 1){
      lwork = 1;
      liwork = 1;
   }
   else if(1 < size && job == 'N'){
      lwork = 2*size + 1;
      liwork = 2;
   }
   else{
      // calc. k
      double temp = log((double)size)/log(2.0);
      if( (double)((int)temp) < temp ){
         k = (int)temp + 1;
      }
      else{
         k = (int)temp;
      }
      lwork = 3*size*size + (5+2*k)*size + 1;
      liwork = 5*size + 3;
   }

   // malloc
   work = (double*)LAPACKE_malloc( sizeof(double)*lwork );
   iwork = (int*)LAPACKE_malloc( sizeof(int)*liwork );
   convertedMatrix = (double*)LAPACKE_malloc( sizeof(double)*size*size );
   tempEigenValues = (double*)LAPACKE_malloc( sizeof(double)*size );

   for(int i = 0; i < size; i++){
      for(int j = i; j < size; j++){
         convertedMatrix[i+j*size] = matrix[i][j];
      }
   }

   // call Lapack
   info = LAPACKE_dsyevd_work(LAPACK_COL_MAJOR, job, uplo, size, convertedMatrix, lda, tempEigenValues, work, lwork, iwork, liwork);

   for(int i = 0; i < size; i++){
      for(int j = 0; j < size; j++){
         matrix[i][j] = convertedMatrix[j+i*size];  //i-th row is i-th eigen vector
         //matrix[j][i] = convertedMatrix[j+i*size];  //i-th column is i-th eigen vector
      }
   }

   for(int i=0;i<size;i++){
      double temp = 0.0;
      for(int j=0;j<size;j++){
         temp += matrix[i][j];
      }
      if(temp<0){
         for(int j=0;j<size;j++){
            matrix[i][j]*=-1.0;
         }
      }
   }   

   for(int i = 0; i < size; i++){
      eigenValues[i] = tempEigenValues[i];
   }
   //this->OutputLog((boost::format("size=%d lwork=%d liwork=%d k=%d info=%d\n") % size % lwork % liwork % k % info ).str());

   // free
   LAPACKE_free(work);
   LAPACKE_free(iwork);
   LAPACKE_free(convertedMatrix);
   LAPACKE_free(tempEigenValues);
  
   if(info != 0){
      stringstream ss;
      ss << errorMessageDsyevdInfo;
      ss << info << endl;
      throw MolDSException(ss.str());
   }
   return info;
}

/***
 *
 * Slove matrix*X=b, then we get X by this method.
 * The X is stored in b.
 *
 */
int Lapack::Dsysv(double const* const* matrix, double* b, int size){
   int info = 0;
   int lwork;
   char uplo = 'U';
   int lda = size;
   int ldb = size;
   int nrhs = 1;
   double* convertedMatrix;
   double* work;
   double* tempB;
   int* ipiv;

   if(size < 1 ){
      stringstream ss;
      ss << errorMessageDsysvSize;
      throw MolDSException(ss.str());
   }

   // malloc
   ipiv = (int*)LAPACKE_malloc( sizeof(int)*2*size );
   convertedMatrix = (double*)LAPACKE_malloc( sizeof(double)*size*size );
   tempB = (double*)LAPACKE_malloc( sizeof(double)*size );

   for(int i = 0; i < size; i++){
      for(int j = i; j < size; j++){
         convertedMatrix[i+j*size] = matrix[i][j];
      }
   }
   for(int i = 0; i < size; i++){
      tempB[i] = b[i];
   }

   // calc. lwork
   #pragma omp critical
   {
      if(!this->calculatedDsysvBlockSize){
         lwork = -1;
         double tempWork[3]={0.0, 0.0, 0.0};
         info = LAPACKE_dsysv_work(LAPACK_COL_MAJOR, uplo, size, nrhs, convertedMatrix, lda, ipiv, tempB, ldb, tempWork, lwork);
         this->calculatedDsysvBlockSize = true;
         this->dsysvBlockSize = tempWork[0]/size;
      }
   }
   info = 0;
   lwork = this->dsysvBlockSize*size;
   work = (double*)LAPACKE_malloc( sizeof(double)*lwork );

   // call Lapack
   info = LAPACKE_dsysv_work(LAPACK_COL_MAJOR, uplo, size, nrhs, convertedMatrix, lda, ipiv, tempB, ldb, work, lwork);
   for(int i = 0; i < size; i++){
      b[i] = tempB[i];
   }

   // free
   LAPACKE_free(convertedMatrix);
   LAPACKE_free(ipiv);
   LAPACKE_free(work);
   LAPACKE_free(tempB);
  
   if(info != 0){
      stringstream ss;
      ss << errorMessageDsysvInfo;
      ss << info << endl;
      throw MolDSException(ss.str());
   }
   return info;
}

/***
 *
 * Slove matrix*X[i]=b[i] (i=0, 1, ... , nrhs-1), then we get X[i] by this method.
 * The X[i] is stored in b[i].
 * b[i][j] is j-th element of i-th solution, b[i].
 *
 */
int Lapack::Dgetrs(double const* const* matrix, double** b, int size, int nrhs) const{
   int info = 0;
   char trans = 'N';
   int lda = size;
   int ldb = size;
   double* convertedMatrix;
   double* convertedB;
   int* ipiv;

   if(size < 1 ){
      stringstream ss;
      ss << errorMessageDgetrsSize;
      throw MolDSException(ss.str());
   }

   try{
      // malloc
      ipiv = (int*)LAPACKE_malloc( sizeof(int)*2*size);
      convertedMatrix = (double*)LAPACKE_malloc( sizeof(double)*size*size);
      convertedB = (double*)LAPACKE_malloc( sizeof(double)*nrhs*size);
      for(int i = 0; i < size; i++){
         for(int j = 0; j < size; j++){
            convertedMatrix[i+j*size] = matrix[i][j];
         }
      }
      for(int i = 0; i < nrhs; i++){
         for(int j = 0; j < size; j++){
            convertedB[j+i*size] = b[i][j];
         }
      }
      this->Dgetrf(convertedMatrix, ipiv, size, size);
      info = LAPACKE_dgetrs_work(LAPACK_COL_MAJOR, trans, size, nrhs, convertedMatrix, lda, ipiv, convertedB, ldb);
      for(int i = 0; i < nrhs; i++){
         for(int j = 0; j < size; j++){
            b[i][j] = convertedB[j+i*size];
         }
      }
   }
   catch(MolDSException ex){
      // free
      LAPACKE_free(convertedMatrix);
      LAPACKE_free(convertedB);
      LAPACKE_free(ipiv);
      throw ex;
   }
   // free
   LAPACKE_free(convertedMatrix);
   LAPACKE_free(convertedB);
   LAPACKE_free(ipiv);
  
   if(info != 0){
      stringstream ss;
      ss << errorMessageDgetrsInfo;
      ss << info << endl;
      throw MolDSException(ss.str());
   }
   return info;
}

// Argument "matrix" means sizeM * sizeN matrix.
// The each element of "matrix" should be stored in 1-dimensional vecotre with column major (Fortran type).
int Lapack::Dgetrf(double* matrix, int* ipiv, int sizeM, int sizeN) const{
   int info = 0;
   int lda = sizeM;
   info = LAPACKE_dgetrf_work(LAPACK_COL_MAJOR, sizeM, sizeN, matrix, lda, ipiv);
   if(info != 0){
      stringstream ss;
      ss << errorMessageDgetrfInfo;
      ss << info << endl;
      throw MolDSException(ss.str());
   }
   return info;
}
}
