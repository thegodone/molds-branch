//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
// Copyright (C) 2012-2012 Katushiko Nishimra                             // 
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

      // Continuous allocation is necessary for matrix to vector conversion.
      *matrix = new T*[size1];
      if(*matrix==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
			T *buf = new T[size1*size2];
      if(buf==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      for(int i=0;i<size1;i++) {
         (*matrix)[i] = &buf[i*size2];
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
      delete [] (*matrix)[0];
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

      T *p1d=NULL, **p2d=NULL, ***p3d=NULL;
      try{
         p1d = new T[size1*size2*size3];
         if(p1d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p2d = new T*[size1*size2];
         if(p2d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p3d = new T**[size1];
         if(p3d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }

         for(int i=0;i<size1;i++){
            p3d[i] = &p2d[i*size2];
            for(int j=0;j<size2;j++){
               p3d[i][j] = &p1d[i*size2*size3+j*size3];
            }
         }
         *matrix = p3d;
         MallocerFreer::AddCurrentMalloced(wannaMalloc);
         this->Initialize<T>(*matrix, size1, size2, size3);
      }
      catch(MolDSException ex){
         if(p1d!=NULL){
            delete[] p1d;
         }
         if(p2d!=NULL){
            delete[] p2d;
         }
         if(p3d!=NULL){
            delete[] p3d;
         }
         throw ex;
      }
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
      T *p1d=NULL, **p2d=NULL, ***p3d=NULL;
      p3d = *matrix;
      p2d = p3d[0];
      p1d = p2d[0];
      delete [] p3d;
      delete [] p2d;
      delete [] p1d;
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

      T *p1d=NULL, **p2d=NULL, ***p3d=NULL, ****p4d=NULL;
      try{
         p1d = new T[size1*size2*size3*size4];
         if(p1d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p2d = new T*[size1*size2*size3];
         if(p2d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p3d = new T**[size1*size2];
         if(p3d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p4d = new T***[size1];
         if(p4d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }

         for(int i=0;i<size1;i++){
            p4d[i] = &p3d[i*size2];
            for(int j=0;j<size2;j++){
               p4d[i][j] = &p2d[i*size2*size3+j*size3];
               for(int k=0;k<size3;k++){
                  p4d[i][j][k] = &p1d[i*size2*size3*size4+j*size3*size4+k*size4];
               }
            }
         }
         *matrix = p4d;
         MallocerFreer::AddCurrentMalloced(wannaMalloc);
         this->Initialize<T>(*matrix, size1, size2, size3, size4);
      }
      catch(MolDSException ex){
         if(p1d!=NULL){
            delete[] p1d;
         }
         if(p2d!=NULL){
            delete[] p2d;
         }
         if(p3d!=NULL){
            delete[] p3d;
         }
         if(p4d!=NULL){
            delete[] p4d;
         }
         throw ex;
      }
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
      T *p1d=NULL, **p2d=NULL, ***p3d=NULL,****p4d=NULL;
      p4d = *matrix;
      p3d = p4d[0];
      p2d = p3d[0];
      p1d = p2d[0];
      delete [] p4d;
      delete [] p3d;
      delete [] p2d;
      delete [] p1d;
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

      T *p1d=NULL, **p2d=NULL, ***p3d=NULL, ****p4d=NULL, *****p5d=NULL;
      try{
         p1d = new T[size1*size2*size3*size4*size5];
         if(p1d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p2d = new T*[size1*size2*size3*size4];
         if(p2d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p3d = new T**[size1*size2*size3];
         if(p3d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p4d = new T***[size1*size2];
         if(p4d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p5d = new T****[size1];
         if(p5d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }

         for(int i=0;i<size1;i++){
            p5d[i] = &p4d[i*size2];
            for(int j=0;j<size2;j++){
               p5d[i][j] = &p3d[i*size2*size3+j*size3];
               for(int k=0;k<size3;k++){
                  p5d[i][j][k] = &p2d[i*size2*size3*size4+j*size3*size4+k*size4];
                  for(int l=0;l<size4;l++){
                     p5d[i][j][k][l] = &p1d[i*size2*size3*size4*size5+
                                            j*size3*size4*size5+
                                            k*size4*size5+
                                            l*size5];
                  }
               }
            }
         }
         *matrix = p5d;
         MallocerFreer::AddCurrentMalloced(wannaMalloc);
         this->Initialize<T>(*matrix, size1, size2, size3, size4, size5);
      }
      catch(MolDSException ex){
         if(p1d!=NULL){
            delete[] p1d;
         }
         if(p2d!=NULL){
            delete[] p2d;
         }
         if(p3d!=NULL){
            delete[] p3d;
         }
         if(p4d!=NULL){
            delete[] p4d;
         }
         if(p5d!=NULL){
            delete[] p5d;
         }
         throw ex;
      }
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

      T *p1d=NULL, **p2d=NULL, ***p3d=NULL,****p4d=NULL, *****p5d=NULL;
      p5d = *matrix;
      p4d = p5d[0];
      p3d = p4d[0];
      p2d = p3d[0];
      p1d = p2d[0];
      delete [] p5d;
      delete [] p4d;
      delete [] p3d;
      delete [] p2d;
      delete [] p1d;
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

      T *p1d=NULL, **p2d=NULL, ***p3d=NULL, ****p4d=NULL, *****p5d=NULL, ******p6d=NULL;
      try{
         p1d = new T[size1*size2*size3*size4*size5*size6];
         if(p1d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p2d = new T*[size1*size2*size3*size4*size5];
         if(p2d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p3d = new T**[size1*size2*size3*size4];
         if(p3d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p4d = new T***[size1*size2*size3];
         if(p4d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p5d = new T****[size1*size2];
         if(p5d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p6d = new T*****[size1];
         if(p6d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }

         for(int i=0;i<size1;i++){
            p6d[i] = &p5d[i*size2];
            for(int j=0;j<size2;j++){
               p6d[i][j] = &p4d[i*size2*size3+j*size3];
               for(int k=0;k<size3;k++){
                  p6d[i][j][k] = &p3d[i*size2*size3*size4+j*size3*size4+k*size4];
                  for(int l=0;l<size4;l++){
                     p6d[i][j][k][l] = &p2d[i*size2*size3*size4*size5+
                                            j*size3*size4*size5+
                                            k*size4*size5+
                                            l*size5];
                     for(int m=0;m<size5;m++){
                        p6d[i][j][k][l][m] = &p1d[i*size2*size3*size4*size5*size6+
                                                  j*size3*size4*size5*size6+
                                                  k*size4*size5*size6+
                                                  l*size5*size6+
                                                  m*size6];
                     }
                  }
               }
            }
         }
         *matrix = p6d;
         MallocerFreer::AddCurrentMalloced(wannaMalloc);
         this->Initialize<T>(*matrix, size1, size2, size3, size4, size5, size6);
      }
      catch(MolDSException ex){
         if(p1d!=NULL){
            delete[] p1d;
         }
         if(p2d!=NULL){
            delete[] p2d;
         }
         if(p3d!=NULL){
            delete[] p3d;
         }
         if(p4d!=NULL){
            delete[] p4d;
         }
         if(p5d!=NULL){
            delete[] p5d;
         }
         if(p6d!=NULL){
            delete[] p6d;
         }
         throw ex;
      }
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
      T *p1d=NULL, **p2d=NULL, ***p3d=NULL,****p4d=NULL, *****p5d=NULL, ******p6d=NULL;
      p6d = *matrix;
      p5d = p6d[0];
      p4d = p5d[0];
      p3d = p4d[0];
      p2d = p3d[0];
      p1d = p2d[0];
      delete [] p6d;
      delete [] p5d;
      delete [] p4d;
      delete [] p3d;
      delete [] p2d;
      delete [] p1d;
      MallocerFreer::SubtCurrentMalloced(static_cast<double>(size1*size2*size3*size4*size5*sizeof(T)));
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
