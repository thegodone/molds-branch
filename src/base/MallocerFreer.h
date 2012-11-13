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
   template<typename T> void Malloc(T** matrix, size_t size1) const{
      if(*matrix!=NULL){
         return;
      }
      double requiredMalloc = this->GetMemoryAmount<T>(size1);
      this->CheckLimitHeap(requiredMalloc);

      *matrix = new T[size1];
      if(*matrix==NULL){
         throw MolDSException(this->errorMessageMallocFailure);
      }
      MallocerFreer::AddCurrentMalloced(requiredMalloc);
      this->Initialize<T>(*matrix, size1);
   }

   template<typename T> void Initialize(T* matrix, size_t size1) const{
      for(size_t i=0;i<size1;i++){
         matrix[i] = (T)0;
      }
   }

   template<typename T> void Free(T** matrix, size_t size1) const{
      if(*matrix==NULL){
         return;
      }
      delete [] *matrix;
      double freedMalloc = this->GetMemoryAmount<T>(size1);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   template<typename T> double GetMemoryAmount(size_t size1) const{
      double amount = static_cast<double>(sizeof(T))
                     *static_cast<double>(size1);
      return amount;
   }

   // 2d
   template<typename T> void Malloc(T*** matrix, size_t size1, size_t size2) const{
      if(*matrix!=NULL){
         return;
      }
      double requiredMalloc = this->GetMemoryAmount<T>(size1, size2);
      this->CheckLimitHeap(requiredMalloc);

      T *p1d=NULL, **p2d=NULL;
      try{
         p1d = new T[size1*size2];
         if(p1d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p2d = new T*[size1];
         if(p2d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }

         for(size_t i=0;i<size1;i++){
            p2d[i] = &p1d[i*size2];
         }
         *matrix = p2d;
         MallocerFreer::AddCurrentMalloced(requiredMalloc);
         this->Initialize<T>(*matrix, size1, size2);
      }
      catch(MolDSException ex){
         if(p1d!=NULL){
            delete[] p1d;
         }
         if(p2d!=NULL){
            delete[] p2d;
         }
         throw ex;
      }
   }

   template<typename T> void Initialize(T** matrix, size_t size1, size_t size2) const{
      for(size_t i=0;i<size1;i++){
         for(size_t j=0;j<size2;j++){
            matrix[i][j] = (T)0.0;
         }
      }
   }

   template<typename T> void Free(T*** matrix, size_t size1, size_t size2) const{
      if(*matrix==NULL){
         return;
      }
      delete [] (*matrix)[0];
      delete [] *matrix;
      double freedMalloc = this->GetMemoryAmount<T>(size1, size2);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   template<typename T> double GetMemoryAmount(size_t size1, size_t size2) const{
      double amount = static_cast<double>(sizeof(T))
                     *static_cast<double>(size1)
                     *static_cast<double>(size2);
      amount += this->GetMemoryAmount<T*>(size1);
      return amount;
   }

   // 3d
   template<typename T> void Malloc(T**** matrix, size_t size1, size_t size2, size_t size3) const{
      if(*matrix!=NULL){
         return;
      }
      double requiredMalloc = this->GetMemoryAmount<T>(size1, size2, size3);
      this->CheckLimitHeap(requiredMalloc);

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

         for(size_t i=0;i<size1;i++){
            p3d[i] = &p2d[i*size2];
            for(size_t j=0;j<size2;j++){
               p3d[i][j] = &p1d[i*size2*size3+j*size3];
            }
         }
         *matrix = p3d;
         MallocerFreer::AddCurrentMalloced(requiredMalloc);
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

   template<typename T> void Initialize(T*** matrix, size_t size1, size_t size2, size_t size3) const{
      for(size_t i=0;i<size1;i++) {
         for(size_t j=0;j<size2;j++){
            for(size_t k=0;k<size3;k++){
               matrix[i][j][k] = (T)0.0;
            }
         }
      }  
   }

   template<typename T> void Free(T**** matrix, size_t size1, size_t size2, size_t size3) const{
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
      double freedMalloc = this->GetMemoryAmount<T>(size1, size2, size3);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   template<typename T> double GetMemoryAmount(size_t size1, size_t size2, size_t size3) const{
      double amount = static_cast<double>(sizeof(T))
                     *static_cast<double>(size1)
                     *static_cast<double>(size2)
                     *static_cast<double>(size3);
      amount += this->GetMemoryAmount<T*>(size1, size2);
      return amount;
   }

   //4d
   template<typename T> void Malloc(T***** matrix, size_t size1, size_t size2, size_t size3, size_t size4) const{
      if(*matrix!=NULL){
         return;
      }
      double requiredMalloc = this->GetMemoryAmount<T>(size1, size2, size3, size4);
      this->CheckLimitHeap(requiredMalloc);

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

         for(size_t i=0;i<size1;i++){
            p4d[i] = &p3d[i*size2];
            for(size_t j=0;j<size2;j++){
               p4d[i][j] = &p2d[i*size2*size3+j*size3];
               for(size_t k=0;k<size3;k++){
                  p4d[i][j][k] = &p1d[i*size2*size3*size4+j*size3*size4+k*size4];
               }
            }
         }
         *matrix = p4d;
         MallocerFreer::AddCurrentMalloced(requiredMalloc);
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

   template<typename T> void Initialize(T**** matrix, size_t size1, size_t size2, size_t size3, size_t size4) const{
      for(size_t i=0;i<size1;i++) {
         for(size_t j=0;j<size2;j++){
            for(size_t k=0;k<size3;k++){
               for(size_t l=0;l<size4;l++){
                  matrix[i][j][k][l] = (T)0.0;
               }
            }
         }
      }
   }

   template<typename T> void Free(T***** matrix, size_t size1, size_t size2, size_t size3, size_t size4) const{
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
      double freedMalloc = this->GetMemoryAmount<T>(size1, size2, size3, size4);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   template<typename T> double GetMemoryAmount(size_t size1, size_t size2, size_t size3, size_t size4) const{
      double amount = static_cast<double>(sizeof(T))
                     *static_cast<double>(size1)
                     *static_cast<double>(size2)
                     *static_cast<double>(size3)
                     *static_cast<double>(size4);
      amount += this->GetMemoryAmount<T*>(size1, size2, size3);
      return amount;
   }

   //5d
   template<typename T> void Malloc(T****** matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5) const{
      if(*matrix!=NULL){
         return;
      }
      double requiredMalloc = this->GetMemoryAmount<T>(size1, size2, size3, size4, size5);
      this->CheckLimitHeap(requiredMalloc);

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

         for(size_t i=0;i<size1;i++){
            p5d[i] = &p4d[i*size2];
            for(size_t j=0;j<size2;j++){
               p5d[i][j] = &p3d[i*size2*size3+j*size3];
               for(size_t k=0;k<size3;k++){
                  p5d[i][j][k] = &p2d[i*size2*size3*size4+j*size3*size4+k*size4];
                  for(size_t l=0;l<size4;l++){
                     p5d[i][j][k][l] = &p1d[i*size2*size3*size4*size5+
                                            j*size3*size4*size5+
                                            k*size4*size5+
                                            l*size5];
                  }
               }
            }
         }
         *matrix = p5d;
         MallocerFreer::AddCurrentMalloced(requiredMalloc);
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

   template<typename T> void Initialize(T***** matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5) const{
      for(size_t i=0;i<size1;i++) {
         for(size_t j=0;j<size2;j++){
            for(size_t k=0;k<size3;k++){
               for(size_t l=0;l<size4;l++){
                  for(size_t m=0;m<size5;m++){
                     matrix[i][j][k][l][m] = 0.0;
                  }
               }
            }
         }
      }
   }

   template<typename T> void Free(T****** matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5) const{
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
      double freedMalloc = this->GetMemoryAmount<T>(size1, size2, size3, size4, size5);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   template<typename T> double GetMemoryAmount(size_t size1, size_t size2, size_t size3, size_t size4, size_t size5) const{
      double amount = static_cast<double>(sizeof(T))
                     *static_cast<double>(size1)
                     *static_cast<double>(size2)
                     *static_cast<double>(size3)
                     *static_cast<double>(size4)
                     *static_cast<double>(size5);
      amount += this->GetMemoryAmount<T*>(size1, size2, size3, size4);
      return amount;
   }

   //6d
   template<typename T> void Malloc(T******* matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6) const{
      if(*matrix!=NULL){
         return;
      }
      double requiredMalloc = this->GetMemoryAmount<T>(size1, size2, size3, size4, size5, size6);
      this->CheckLimitHeap(requiredMalloc);

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

         for(size_t i=0;i<size1;i++){
            p6d[i] = &p5d[i*size2];
            for(size_t j=0;j<size2;j++){
               p6d[i][j] = &p4d[i*size2*size3+j*size3];
               for(size_t k=0;k<size3;k++){
                  p6d[i][j][k] = &p3d[i*size2*size3*size4+j*size3*size4+k*size4];
                  for(size_t l=0;l<size4;l++){
                     p6d[i][j][k][l] = &p2d[i*size2*size3*size4*size5+
                                            j*size3*size4*size5+
                                            k*size4*size5+
                                            l*size5];
                     for(size_t m=0;m<size5;m++){
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
         MallocerFreer::AddCurrentMalloced(requiredMalloc);
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

   template<typename T> void Initialize(T****** matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6) const{
      for(size_t i=0;i<size1;i++) {
         for(size_t j=0;j<size2;j++){
            for(size_t k=0;k<size3;k++){
               for(size_t l=0;l<size4;l++){
                  for(size_t m=0;m<size5;m++){
                     for(size_t n=0;n<size6;n++){
                        matrix[i][j][k][l][m][n] = 0.0;
                     }
                  }
               }
            }
         }
      }
   }

   template<typename T> void Free(T******* matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6) const{
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
      double freedMalloc = this->GetMemoryAmount<T>(size1, size2, size3, size4, size5, size6);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   template<typename T> double GetMemoryAmount(size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6) const{
      double amount = static_cast<double>(sizeof(T))
                     *static_cast<double>(size1)
                     *static_cast<double>(size2)
                     *static_cast<double>(size3)
                     *static_cast<double>(size4)
                     *static_cast<double>(size5)
                     *static_cast<double>(size6);
      amount += this->GetMemoryAmount<T*>(size1, size2, size3, size4, size5);
      return amount;
   }

   //7d
   template<typename T> void Malloc(T******** matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6, size_t size7) const{
      if(*matrix!=NULL){
         return;
      }
      double requiredMalloc = this->GetMemoryAmount<T>(size1, size2, size3, size4, size5, size6, size7);
      this->CheckLimitHeap(requiredMalloc);

      T *p1d=NULL, **p2d=NULL, ***p3d=NULL, ****p4d=NULL, *****p5d=NULL, ******p6d=NULL, *******p7d=NULL;
      try{
         p1d = new T[size1*size2*size3*size4*size5*size6*size7];
         if(p1d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p2d = new T*[size1*size2*size3*size4*size5*size6];
         if(p2d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p3d = new T**[size1*size2*size3*size4*size5];
         if(p3d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p4d = new T***[size1*size2*size3*size4];
         if(p4d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p5d = new T****[size1*size2*size3];
         if(p5d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p6d = new T*****[size1*size2];
         if(p6d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }
         p7d = new T******[size1];
         if(p7d==NULL){
            throw MolDSException(this->errorMessageMallocFailure);
         }

         for(size_t i=0;i<size1;i++){
            p7d[i] = &p6d[i*size2];
            for(size_t j=0;j<size2;j++){
               p7d[i][j] = &p5d[i*size2*size3+j*size3];
               for(size_t k=0;k<size3;k++){
                  p7d[i][j][k] = &p4d[i*size2*size3*size4+j*size3*size4+k*size4];
                  for(size_t l=0;l<size4;l++){
                     p7d[i][j][k][l] = &p3d[i*size2*size3*size4*size5+
                                            j*size3*size4*size5+
                                            k*size4*size5+
                                            l*size5];
                     for(size_t m=0;m<size5;m++){
                        p7d[i][j][k][l][m] = &p2d[i*size2*size3*size4*size5*size6+
                                                  j*size3*size4*size5*size6+
                                                  k*size4*size5*size6+
                                                  l*size5*size6+
                                                  m*size6];
                        for(size_t n=0;n<size6;n++){
                           p7d[i][j][k][l][m][n] = &p1d[i*size2*size3*size4*size5*size6*size7+
                                                        j*size3*size4*size5*size6*size7+
                                                        k*size4*size5*size6*size7+
                                                        l*size5*size6*size7+
                                                        m*size6*size7+
                                                        n*size7];
                     }
                     }
                  }
               }
            }
         }
         *matrix = p7d;
         MallocerFreer::AddCurrentMalloced(requiredMalloc);
         this->Initialize<T>(*matrix, size1, size2, size3, size4, size5, size6, size7);
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
         if(p7d!=NULL){
            delete[] p7d;
         }
         throw ex;
      }
   }

   template<typename T> void Initialize(T******* matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6, size_t size7) const{
      for(size_t i=0;i<size1;i++) {
         for(size_t j=0;j<size2;j++){
            for(size_t k=0;k<size3;k++){
               for(size_t l=0;l<size4;l++){
                  for(size_t m=0;m<size5;m++){
                     for(size_t n=0;n<size6;n++){
                        for(size_t o=0;o<size7;o++){
                           matrix[i][j][k][l][m][n][o] = 0.0;
                        }
                     }
                  }
               }
            }
         }
      }
   }

   template<typename T> void Free(T******** matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6, size_t size7) const{
      if(*matrix==NULL){
         return;
      }
      T *p1d=NULL, **p2d=NULL, ***p3d=NULL,****p4d=NULL, *****p5d=NULL, ******p6d=NULL, *******p7d=NULL;
      p7d = *matrix;
      p6d = p7d[0];
      p5d = p6d[0];
      p4d = p5d[0];
      p3d = p4d[0];
      p2d = p3d[0];
      p1d = p2d[0];
      delete [] p7d;
      delete [] p6d;
      delete [] p5d;
      delete [] p4d;
      delete [] p3d;
      delete [] p2d;
      delete [] p1d;
      double freedMalloc = this->GetMemoryAmount<T>(size1, size2, size3, size4, size5, size6, size7);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   template<typename T> double GetMemoryAmount(size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6, size_t size7) const{
      double amount = static_cast<double>(sizeof(T))
                     *static_cast<double>(size1)
                     *static_cast<double>(size2)
                     *static_cast<double>(size3)
                     *static_cast<double>(size4)
                     *static_cast<double>(size5)
                     *static_cast<double>(size6)
                     *static_cast<double>(size7);
      amount += this->GetMemoryAmount<T*>(size1, size2, size3, size4, size5, size6);
      return amount;
   }
private:
   MallocerFreer();
   ~MallocerFreer();
   static MallocerFreer* mallocerFreer;
   static double currentMalloced;
   static double maxMalloced;
   static const double byte2MByte;
   static void AddCurrentMalloced(double amount);
   std::string errorMessageMallocFailure;
   std::string errorMessageReachHeapLimit;
   std::string messageMemoryUsage;
   std::string messageMemoryLeakedHeap;
   std::string messageMemoryMaxHeap;
   std::string messageMemoryCurrentHeap;
   std::string messageMemoryRequiredHeap;
   std::string messageMemoryLimitHeap;
   std::string messageMByte;
   void OutputMemoryUsage() const;
   void CheckLimitHeap(double requiredMalloc) const;
};
}
#endif
