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
#include"../config.h"
#include"../base/Enums.h"
#include"../base/Uncopyable.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
#include"ScaLapack.h"

#ifdef __INTEL_COMPILER
#define MOLDS_SCALAPACK_malloc(a,b) mkl_malloc(a,b)
#define MOLDS_SCALAPACK_free(a) mkl_free(a)
#else
#define MOLDS_SCALAPACK_malloc(a,b) malloc(a)
#define MOLDS_SCALAPACK_free(a) free(a)
#endif

using namespace std;
using namespace MolDS_base;

namespace MolDS_wrappers{

#ifdef __FCC_VERSION
extern "C" {
   extern void                descinit_(molds_scalapack_int* desc, molds_scalapack_int* m, molds_scalapack_int* n, molds_scalapack_int* mb, molds_scalapack_int* nb, molds_scalapack_int* irsrc, molds_scalapack_int* icsrc, molds_scalapack_int* ictxt, molds_scalapack_int* lld, molds_scalapack_int* info);
   extern void                blac_gridexit_(molds_scalapack_int* iContext);
   extern molds_scalapack_int numroc_(molds_scalapack_int* n, molds_scalapack_int* nb, molds_scalapack_int* iproc, molds_scalapack_int* isrcproc, molds_scalapack_int* nprocs);
   extern void                sl_init_(molds_scalapack_int* iContext, molds_scalapack_int* npRow, molds_scalapack_int* npCol); 
   extern void                blacs_pinfo_(molds_scalapack_int *mypnum, molds_scalapack_int *nprocs);
   extern void                blacs_gridinfo_(molds_scalapack_int *ConTxt, molds_scalapack_int *nprow, molds_scalapack_int *npcol, molds_scalapack_int *myrow, molds_scalapack_int *mycol);
   extern void                blacs_gridexit_(molds_scalapack_int *ConTxt);
   extern void                pdsyevd_(char *jobz, char *uplo, molds_scalapack_int *n, double *a, molds_scalapack_int *ia, molds_scalapack_int *ja, molds_scalapack_int *desca, double *w, double *z, molds_scalapack_int *iz, molds_scalapack_int *jz, molds_scalapack_int *descz, double *work, molds_scalapack_int *lwork, molds_scalapack_int *iwork, molds_scalapack_int *liwork, molds_scalapack_int *info);
}

#ifdef F77_WITH_NO_UNDERSCORE
#define   descinit_        descinit
#define   pdsyevd_         pdsyevd
#define   blacs_gridinfo_  blacs_gridinfo
#define   blacs_gridexit_  blacs_gridexit
#define   blacs_pinfo_     blacs_pinfo
#define   sl_init_         sl_init
#define   numroc_          numroc
#endif

#endif

ScaLapack* ScaLapack::scaLapack = NULL;
ScaLapack::ScaLapack(){
   this->errorMessagePdsyevdInfo         = "Error in wrappers::ScaLapack::Pdsyevd: info != 0: info = ";
   this->errorMessagePdsyevdSize         = "Error in wrappers::ScaLapack::Pdsyevd: size of matirx < 1\n";
   this->errorMessagePdsyevdNotSupported = "Error in wrappers::ScaLapack::Pdsyevd: ScaLapack is not supported on the current system. ScaLapack is supported only on FX10.\n";
}

ScaLapack::~ScaLapack(){
}

ScaLapack* ScaLapack::GetInstance(){
   if(scaLapack == NULL){
      scaLapack = new ScaLapack();
      //this->OutputLog("ScaLapack created.\n\n");
   }
   return scaLapack;
}

void ScaLapack::DeleteInstance(){
   if(scaLapack != NULL){
      delete scaLapack;
      //this->OutputLog("ScaLapack deleted\n\n");
   }
   scaLapack = NULL;
}

/***
 *  *
 *  * Eigenvalue and eigenvector of a real symmetirc matrix are calculated:
 *  *   - i-th eigenvalue will be stored in eigenValues[i].
 *  *   - i-th eigenvector will be stored as (matirx[i][0], matirx[i][1], matirx[i][2], ....).
 *  *
 *  * ***/
molds_scalapack_int ScaLapack::Pdsyevd(double** matrix, double* eigenValues, molds_scalapack_int size, bool calcEigenVectors){
   molds_scalapack_int info = 0;
#ifdef __FCC_VERSION
   if(size < 1 ){
      stringstream ss;
      ss << errorMessagePdsyevdSize;
      MolDSException ex(ss.str());
      throw ex;
   }

   char job;
   char uplo = 'U';
   double* localMatrix;
   double* localEigenVector;
   double* tempEigenValues;
   // set job type
   if(calcEigenVectors){
      job = 'V';
   }
   else{
      job = 'N';
   }

   //tmporal values
   molds_scalapack_int intOne  = 1;
   molds_scalapack_int intZero = 0;
   molds_scalapack_int myRow   = intZero;
   molds_scalapack_int myCol   = intZero;
   molds_scalapack_int iContext= intZero;
   molds_scalapack_int rsrc    = intZero;
   molds_scalapack_int csrc    = intZero;
   molds_scalapack_int mpiRank = intZero;
   molds_scalapack_int mpiSize = intZero;

   // initialize blacs and scalapack 
   blacs_pinfo_(&mpiRank, &mpiSize);
   molds_scalapack_int squareMpiSize    = lround(sqrt(static_cast<double>(mpiSize)));
   if(mpiSize < squareMpiSize*squareMpiSize){squareMpiSize-=1;}
   molds_scalapack_int npRow            = squareMpiSize;
   molds_scalapack_int npCol            = squareMpiSize;
   molds_scalapack_int blockSizeDefault = 128; 
   molds_scalapack_int blockSize        = blockSizeDefault;
   while(size/2 < blockSize*npRow){
      blockSize /= 2;
      if(blockSize < 1){
         blockSize = 1;
         break;
      }
   }
   sl_init_(&iContext, &npRow, &npCol); 
   blacs_gridinfo_(&iContext, &npRow, &npCol, &myRow, &myCol);
   
   // calculate size of local matrix on each node
   molds_scalapack_int mp   = numroc_(&size, &blockSize, &myRow, &rsrc, &npRow);
   molds_scalapack_int nq   = numroc_(&size, &blockSize, &myCol, &csrc, &npCol);
   localMatrix              = (double*)MOLDS_SCALAPACK_malloc( sizeof(double)*mp*nq, 16 );
   localEigenVector         = (double*)MOLDS_SCALAPACK_malloc( sizeof(double)*mp*nq, 16 );
   tempEigenValues          = (double*)MOLDS_SCALAPACK_malloc( sizeof(double)*size, 16 );
   molds_scalapack_int lldA = max(intOne,mp);
   molds_scalapack_int lldZ = max(intOne,mp);
   molds_scalapack_int descA[9];
   molds_scalapack_int descZ[9];
   descinit_(descA, &size, &size, &blockSize, &blockSize, &intZero, &intZero, &iContext, &lldA,       &info);
   descinit_(descZ, &size, &size, &blockSize, &blockSize, &intZero, &intZero, &iContext, &lldZ,       &info);
   // distribute data to local matrix on each node accoding to the block syclic decomposition
#pragma omp for schedule(dynamic, MOLDS_OMP_DYNAMIC_CHUNK_SIZE)
   for(molds_scalapack_int i=0; i<mp; i++){
      for(molds_scalapack_int j=0; j<nq; j++){
         molds_scalapack_int globalI = (myRow + (i/blockSize)*npRow)*blockSize + i%blockSize;
         molds_scalapack_int globalJ = (myCol + (j/blockSize)*npCol)*blockSize + j%blockSize;
         localMatrix[j*mp+i] = matrix[globalI][globalJ];
      }
   }

   // calculate working array space for pdsyevd
   molds_scalapack_int lwork       = -1;
   molds_scalapack_int trilwmin    = 3*size + std::max(blockSize*(mp+1), 3*blockSize);
   molds_scalapack_int lworkMin    = std::max(1+6*size+2*mp*nq, trilwmin) + 2*size + 1;
   molds_scalapack_int liwork      = 1;
   molds_scalapack_int liworkMin   = 7*size + 8*npCol + 2;
   double              tmpWork[3]  = {0.0, 0.0, 0.0};
   molds_scalapack_int tmpIwork[3] = {0, 0, 0};
   molds_scalapack_int ia          = intOne;
   molds_scalapack_int ja          = intOne;
   molds_scalapack_int iz          = intOne;
   molds_scalapack_int jz          = intOne;

   pdsyevd_(&job, &uplo, &size, localMatrix, &ia, &ja, descA, tempEigenValues, localEigenVector, &iz, &jz ,descZ, tmpWork, &lwork, tmpIwork, &liwork, &info);

   // cll scalapack (pdsyevd)
   info = 0;
   double*               work = NULL;
   molds_scalapack_int* iwork = NULL;
   lwork                      = std::max(lworkMin,  static_cast<molds_scalapack_int>(tmpWork[0])+1);
   liwork                     = std::max(liworkMin, tmpIwork[0]);
   work                       = (double*)MOLDS_SCALAPACK_malloc( sizeof(double)*lwork, 16 );
   iwork                      = (molds_scalapack_int*)MOLDS_SCALAPACK_malloc( sizeof(molds_scalapack_int)*liwork, 16 );
   pdsyevd_(&job, &uplo, &size, localMatrix, &ia, &ja, descA, eigenValues, localEigenVector, &iz, &jz ,descZ, work, &lwork, iwork, &liwork, &info);

   // gathter block cyclic decomposed eigenvector to the global data
   MallocerFreer::GetInstance()->Initialize<double>(matrix, size, size);
#pragma omp for schedule(dynamic, MOLDS_OMP_DYNAMIC_CHUNK_SIZE)
   for(molds_scalapack_int j=0; j<nq; j++){
      molds_scalapack_int globalJ = (myCol + (j/blockSize)*npCol)*blockSize + j%blockSize;
      for(molds_scalapack_int i=0; i<mp; i++){
         molds_scalapack_int globalI = (myRow + (i/blockSize)*npRow)*blockSize + i%blockSize;
         matrix[globalJ][globalI] = localEigenVector[j*mp+i]; // globalj-th row    is globalj-th eigen vector
         //matrix[globalI][globalJ] = localEigenVector[j*mp+i]; // globalj-th column is globalj-th eigen vector
      }
   }
   MolDS_mpi::MpiProcess::GetInstance()->AllReduce(&matrix[0][0], size*size, std::plus<double>());

#pragma omp for schedule(dynamic, MOLDS_OMP_DYNAMIC_CHUNK_SIZE)
   for(molds_scalapack_int i=0;i<size;i++){
      double temp = 0.0;
      for(molds_scalapack_int j=0;j<size;j++){
         temp += matrix[i][j];
      }
      if(temp<0){
         for(molds_scalapack_int j=0;j<size;j++){
            matrix[i][j]*=-1.0;
         }
      }
   }   

   // free
   MOLDS_SCALAPACK_free(work);
   MOLDS_SCALAPACK_free(iwork);
   MOLDS_SCALAPACK_free(localMatrix);
   MOLDS_SCALAPACK_free(localEigenVector);
   MOLDS_SCALAPACK_free(tempEigenValues);

   blacs_gridexit_(&iContext);

#else
   stringstream ss;
   ss << errorMessagePdsyevdNotSupported;
   MolDSException ex(ss.str());
   throw ex;
#endif
   if(info != 0){
      stringstream ss;
      ss << errorMessagePdsyevdInfo;
      ss << info << endl;
      MolDSException ex(ss.str());
      ex.SetKeyValue<int>(LapackInfo, info);
      throw ex;
   }
   return info;
}


}
