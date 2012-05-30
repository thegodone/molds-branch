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
#ifndef INCLUDED_MALLOCERFREER
#define INCLUDED_MALLOCERFREER
namespace MolDS_base{

// MallocerFreer is singleton
class MallocerFreer: public PrintController, private Uncopyable{
public:
   static MallocerFreer* GetInstance();
   static void DeleteInstance();

   //1d
   template<typename T> void Malloc(T** matrix, int size1) const{
      if(*matrix!=NULL){
         return;
      }
      double wannaMalloc = static_cast<double>(size1*sizeof(T));
      this->CheckLimitHeap(wannaMalloc);

      *matrix = new T[size1];
      MallocerFreer::AddCurrentMalloced(wannaMalloc);
      this->Initialize<T>(*matrix, size1);
   }

   template<typename T> void Initialize(T* matrix, int size1) const{
      for(int i=0;i<size1;i++){
         matrix[i] = (T)0;
      }
   }

   template<typename T> void Free(T** matrix, int size1) const{
      if(*matrix==NULL){
         return;
      }
      delete [] *matrix;
      MallocerFreer::SubtCurrentMalloced(static_cast<double>(size1*sizeof(T)));
      *matrix = NULL;
   }

   //2d
   template<typename T> void Malloc(T*** matrix, int size1, int size2) const{
      if(*matrix!=NULL){
         return;
      }
      double wannaMalloc = static_cast<double>(size1*size2*sizeof(T));
      this->CheckLimitHeap(wannaMalloc);

      *matrix = new T*[size1];
      if(*matrix==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int i=0;i<size1;i++) {
         (*matrix)[i] = new T[size2];
         if ((*matrix)[i]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
      }
      MallocerFreer::AddCurrentMalloced(wannaMalloc);
      this->Initialize<T>(*matrix, size1, size2);
   }

   template<typename T> void Initialize(T** matrix, int size1, int size2) const{
      for(int i=0;i<size1;i++){
         for(int j=0;j<size2;j++){
            matrix[i][j] = (T)0.0;
         }
      }
   }

   template<typename T> void Free(T*** matrix, int size1, int size2) const{
      if(*matrix==NULL){
         return;
      }
      for(int i=0;i<size1;i++){
         delete [] (*matrix)[i];
      }
      delete [] *matrix;
      MallocerFreer::SubtCurrentMalloced(static_cast<double>(size1*size2*sizeof(T)));
      *matrix = NULL;
   }

   // 3d
   template<typename T> void Malloc(T**** matrix, int size1, int size2, int size3) const{
      if(*matrix!=NULL){
         return;
      }
      double wannaMalloc = static_cast<double>(size1*size2*size3*sizeof(T));
      this->CheckLimitHeap(wannaMalloc);

      *matrix = new T**[size1];
      if(*matrix==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int i=0;i<size1;i++) {
         (*matrix)[i] = new T*[size2];
         if((*matrix)[i]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         for(int j=0;j<size2;j++){
            (*matrix)[i][j] = new T[size3];
            if((*matrix)[i][j]==NULL){
               throw MolDSException(this->errorMessageMallocFailure);
            }
         }
      }
      MallocerFreer::AddCurrentMalloced(wannaMalloc);
      this->Initialize<T>(*matrix, size1, size2, size3);
   }

   template<typename T> void Initialize(T*** matrix, int size1, int size2, int size3) const{
      for(int i=0;i<size1;i++) {
         for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
               matrix[i][j][k] = (T)0.0;
            }
         }
      }  
   }

   template<typename T> void Free(T**** matrix, int size1, int size2, int size3) const{
      if(*matrix==NULL){
         return;
      }
      for (int i=0;i<size1;i++) {
         for (int j=0;j<size2;j++) {
            delete [] (*matrix)[i][j];
         }
         delete [] (*matrix)[i];
      }
      delete [] *matrix;
      MallocerFreer::SubtCurrentMalloced(static_cast<double>(size1*size2*size3*sizeof(T)));
      *matrix = NULL;
   }

   //4d
   template<typename T> void Malloc(T***** matrix, int size1, int size2, int size3, int size4) const{
      if(*matrix!=NULL){
         return;
      }
      double wannaMalloc = static_cast<double>(size1*size2*size3*size4*sizeof(T));
      this->CheckLimitHeap(wannaMalloc);

      *matrix = new T***[size1];
      if(*matrix==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int i=0;i<size1;i++) {
         (*matrix)[i] = new T**[size2];
         if((*matrix)[i]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         for(int j=0;j<size2;j++){
            (*matrix)[i][j] = new T*[size3];
            if((*matrix)[i][j]==NULL){
               throw MolDSException(this->errorMessageMallocFailure);
            }
            for(int k=0;k<size3;k++){
               (*matrix)[i][j][k] = new T[size4];
               if((*matrix)[i][j][k]==NULL){
                  throw MolDSException(this->errorMessageMallocFailure);
               }
            }
         }
      }
      MallocerFreer::AddCurrentMalloced(wannaMalloc);
      this->Initialize<T>(*matrix, size1, size2, size3, size4);
   }

   template<typename T> void Initialize(T**** matrix, int size1, int size2, int size3, int size4) const{
      for(int i=0;i<size1;i++) {
         for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
               for(int l=0;l<size4;l++){
                  matrix[i][j][k][l] = (T)0.0;
               }
            }
         }
      }
   }

   template<typename T> void Free(T***** matrix, int size1, int size2, int size3, int size4) const{
      if(*matrix==NULL){
         return;
      }
      for (int i=0;i<size1;i++) {
         for (int j=0;j<size2;j++) {
            for (int k=0;k<size3;k++) {
               delete [] (*matrix)[i][j][k];
            }
            delete [] (*matrix)[i][j];
         }
         delete [] (*matrix)[i];
      }
      delete [] *matrix;
      MallocerFreer::SubtCurrentMalloced(static_cast<double>(size1*size2*size3*size4*sizeof(T)));
      *matrix = NULL;
   }

   //5d
   template<typename T> void Malloc(T****** matrix, int size1, int size2, int size3, int size4, int size5) const{
      if(*matrix!=NULL){
         return;
      }
      double wannaMalloc = static_cast<double>(size1*size2*size3*size4*size5*sizeof(T));
      this->CheckLimitHeap(wannaMalloc);

      *matrix = new T****[size1];
      if(*matrix==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int i=0;i<size1;i++) {
         (*matrix)[i] = new T***[size2];
         if((*matrix)[i]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         for(int j=0;j<size2;j++){
            (*matrix)[i][j] = new T**[size3];
            if((*matrix)[i][j]==NULL){
               throw MolDSException(this->errorMessageMallocFailure);
            }
            for(int k=0;k<size3;k++){
               (*matrix)[i][j][k] = new T*[size4];
               if((*matrix)[i][j][k]==NULL){
                  throw MolDSException(this->errorMessageMallocFailure);
               }
               for(int l=0;l<size4;l++){
                  (*matrix)[i][j][k][l] = new T[size5];
                  if((*matrix)[i][j][k][l]==NULL){
                     throw MolDSException(this->errorMessageMallocFailure);
                  }
               }
            }
         }
      }
      MallocerFreer::AddCurrentMalloced(wannaMalloc);
      this->Initialize<T>(*matrix, size1, size2, size3, size4, size5);
   }

   template<typename T> void Initialize(T***** matrix, int size1, int size2, int size3, int size4, int size5) const{
      for(int i=0;i<size1;i++) {
         for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
               for(int l=0;l<size4;l++){
                  for(int m=0;m<size5;m++){
                     matrix[i][j][k][l][m] = 0.0;
                  }
               }
            }
         }
      }
   }

   template<typename T> void Free(T****** matrix, int size1, int size2, int size3, int size4, int size5) const{
      if(*matrix==NULL){
         return;
      }
      for (int i=0;i<size1;i++) {
         for (int j=0;j<size2;j++) {
            for (int k=0;k<size3;k++) {
               for (int l=0;l<size4;l++) {
                  delete [] (*matrix)[i][j][k][l];
               }
               delete [] (*matrix)[i][j][k];
            }
            delete [] (*matrix)[i][j];
         }
         delete [] (*matrix)[i];
      }
      delete [] *matrix;
      MallocerFreer::SubtCurrentMalloced(static_cast<double>(size1*size2*size3*size4*size5*sizeof(T)));
      *matrix = NULL;
   }

   //6d
   template<typename T> void Malloc(T******* matrix, int size1, int size2, int size3, int size4, int size5, int size6) const{
      if(*matrix!=NULL){
         return;
      }
      double wannaMalloc = static_cast<double>(size1*size2*size3*size4*size5*size6*sizeof(T));
      this->CheckLimitHeap(wannaMalloc);

      *matrix = new T*****[size1];
      if(*matrix==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int i=0;i<size1;i++) {
         (*matrix)[i] = new T****[size2];
         if((*matrix)[i]==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         for(int j=0;j<size2;j++){
            (*matrix)[i][j] = new T***[size3];
            if((*matrix)[i][j]==NULL){
               throw MolDSException(this->errorMessageMallocFailure);
            }
            for(int k=0;k<size3;k++){
               (*matrix)[i][j][k] = new T**[size4];
               if((*matrix)[i][j][k]==NULL){
                  throw MolDSException(this->errorMessageMallocFailure);
               }
               for(int l=0;l<size4;l++){
                  (*matrix)[i][j][k][l] = new T*[size5];
                  if((*matrix)[i][j][k][l]==NULL){
                     throw MolDSException(this->errorMessageMallocFailure);
                  }
                  for(int m=0;m<size5;m++){
                     (*matrix)[i][j][k][l][m] = new T[size6];
                     if((*matrix)[i][j][k][l][m]==NULL){
                        throw MolDSException(this->errorMessageMallocFailure);
                     }
                  }
               }
            }
         }
      }
      MallocerFreer::AddCurrentMalloced(wannaMalloc);
      this->Initialize<T>(*matrix, size1, size2, size3, size4, size5, size6);
   }

   template<typename T> void Initialize(T****** matrix, int size1, int size2, int size3, int size4, int size5, int size6) const{
      for(int i=0;i<size1;i++) {
         for(int j=0;j<size2;j++){
            for(int k=0;k<size3;k++){
               for(int l=0;l<size4;l++){
                  for(int m=0;m<size5;m++){
                     for(int n=0;n<size6;n++){
                        matrix[i][j][k][l][m][n] = 0.0;
                     }
                  }
               }
            }
         }
      }
   }

   template<typename T> void Free(T******* matrix, int size1, int size2, int size3, int size4, int size5, int size6) const{
      if(*matrix==NULL){
         return;
      }
      for (int i=0;i<size1;i++) {
         for (int j=0;j<size2;j++) {
            for (int k=0;k<size3;k++) {
               for (int l=0;l<size4;l++) {
                  for (int m=0;m<size5;m++) {
                     delete [] (*matrix)[i][j][k][l][m];
                  }
                  delete [] (*matrix)[i][j][k][l];
               }
               delete [] (*matrix)[i][j][k];
            }
            delete [] (*matrix)[i][j];
         }
         delete [] (*matrix)[i];
      }
      delete [] *matrix;
      MallocerFreer::SubtCurrentMalloced(static_cast<double>(size1*size2*size3*size4*size5*size6*sizeof(T)));
      *matrix = NULL;
   }
private:
   MallocerFreer();
   ~MallocerFreer();
   static MallocerFreer* mallocerFreer;
   static double currentMalloced;
   static double maxMalloced;
   static void AddCurrentMalloced(double amount);
   static void SubtCurrentMalloced(double amount);
   std::string errorMessageMallocFailure;
   std::string errorMessageReachHeapLimit;
   std::string messageMemoryUsage;
   std::string messageMemoryCurrentHeap;
   std::string messageMemoryMaxHeap;
   std::string messageMByte;
   void OutputMemoryUsage() const;
   void CheckLimitHeap(double wannaMalloc) const;
};
}
#endif
