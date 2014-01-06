//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   //
// Copyright (C) 2012-2014 Katushiko Nishimra                             // 
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
      if(*matrix!=NULL) return;
      if(size1<=0)      return;
      double requiredMalloc = this->GetMemoryAmount<T>(size1);
      this->CheckLimitHeap(requiredMalloc);
      *matrix = new T[size1]();
      if(*matrix==NULL) throw MolDSException(this->errorMessageMallocFailure);
      MallocerFreer::AddCurrentMalloced(requiredMalloc);
   }

   template<typename T> void Initialize(T* matrix, size_t size1) const{
      for(size_t i=0;i<size1;i++){
         matrix[i] = (T)0;
      }
   }

   template<typename T> void Free(T** matrix, size_t size1) const{
      if(*matrix==NULL) return;
      delete [] *matrix;
      double freedMalloc = this->GetMemoryAmount<T>(size1);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   template<typename T> double GetMemoryAmount(size_t size1) const{
      return static_cast<double>(sizeof(T))*static_cast<double>(size1);
   }

   // 2d
   template<typename T> void Malloc(T*** matrix, size_t size1, size_t size2) const{
      if(*matrix!=NULL) return;
      if(size1*size2<=0)      return;
      double requiredMalloc = this->GetMemoryAmount<T*>(size1);
      this->CheckLimitHeap(requiredMalloc);

      T *p1d=NULL, **p2d=NULL;
      try{
         this->Malloc<T>(&p1d, size1*size2);
         p2d = new T*[size1];
         if(p2d==NULL) throw MolDSException(this->errorMessageMallocFailure);

         for(size_t i=0;i<size1;i++){p2d[i] = &p1d[i*size2];}

         MallocerFreer::AddCurrentMalloced(requiredMalloc);
         *matrix = p2d;
      }
      catch(MolDSException ex){
         this->Free<T>(&p1d, size1*size2);
         if(p2d!=NULL) delete[] p2d;
         throw ex;
      }
   }

   template<typename T> void Initialize(T** matrix, size_t size1, size_t size2) const{
      for(size_t i=0;i<size1;i++){
         this->Initialize<T>(matrix[i], size2);
      }
   }

   template<typename T> void Free(T*** matrix, size_t size1, size_t size2) const{
      if(*matrix==NULL) return;
      T *p1d=NULL, **p2d=NULL;
      p2d = *matrix;
      p1d = p2d[0];
      delete [] p2d;
      this->Free<T>(&p1d, size1*size2);
      double freedMalloc = this->GetMemoryAmount<T*>(size1);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   // 3d
   template<typename T> void Malloc(T**** matrix, size_t size1, size_t size2, size_t size3) const{
      if(*matrix!=NULL) return;
      if(size1*size2*size3<=0) return;
      double requiredMalloc = this->GetMemoryAmount<T**>(size1);
      this->CheckLimitHeap(requiredMalloc);

      T **p2d=NULL, ***p3d=NULL;
      try{
         this->Malloc<T>(&p2d, size1*size2, size3);
         p3d = new T**[size1];
         if(p3d==NULL) throw MolDSException(this->errorMessageMallocFailure);

         for(size_t i=0;i<size1;i++){p3d[i]    = &p2d[i*size2];
         for(size_t j=0;j<size2;j++){p3d[i][j] = &p2d[i*size2][j*size3];
         }}

         MallocerFreer::AddCurrentMalloced(requiredMalloc);
         *matrix = p3d;
      }
      catch(MolDSException ex){
         this->Free<T>(&p2d, size1*size2, size3);
         if(p3d!=NULL) delete[] p3d;
         throw ex;
      }
   }

   template<typename T> void Initialize(T*** matrix, size_t size1, size_t size2, size_t size3) const{
      for(size_t i=0;i<size1;i++) {
         this->Initialize<T>(matrix[i], size2, size3);
      }  
   }

   template<typename T> void Free(T**** matrix, size_t size1, size_t size2, size_t size3) const{
      if(*matrix==NULL) return;
      T **p2d=NULL, ***p3d=NULL;
      p3d = *matrix;
      p2d = p3d[0];
      delete [] p3d;
      this->Free<T>(&p2d, size1*size2, size3);
      double freedMalloc = this->GetMemoryAmount<T**>(size1);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   //4d
   template<typename T> void Malloc(T***** matrix, size_t size1, size_t size2, size_t size3, size_t size4) const{
      if(*matrix!=NULL) return;
      if(size1*size2*size3*size4<=0) return;
      double requiredMalloc = this->GetMemoryAmount<T***>(size1);
      this->CheckLimitHeap(requiredMalloc);

      T ***p3d=NULL, ****p4d=NULL;
      try{
         this->Malloc<T>(&p3d, size1*size2, size3, size4);
         p4d = new T***[size1];
         if(p4d==NULL) throw MolDSException(this->errorMessageMallocFailure);

         for(size_t i=0;i<size1;i++){p4d[i]       = &p3d[i*size2];
         for(size_t j=0;j<size2;j++){p4d[i][j]    = &p3d[i*size2][j*size3];
         for(size_t k=0;k<size3;k++){p4d[i][j][k] = &p3d[i*size2][j*size3][k*size4];
         }}}

         MallocerFreer::AddCurrentMalloced(requiredMalloc);
         *matrix = p4d;
      }
      catch(MolDSException ex){
         this->Free<T>(&p3d, size1*size2, size3, size4);
         if(p4d!=NULL) delete[] p4d;
         throw ex;
      }
   }

   template<typename T> void Initialize(T**** matrix, size_t size1, size_t size2, size_t size3, size_t size4) const{
      for(size_t i=0;i<size1;i++) {
         this->Initialize<T>(matrix[i], size2, size3, size4);
      }
   }

   template<typename T> void Free(T***** matrix, size_t size1, size_t size2, size_t size3, size_t size4) const{
      if(*matrix==NULL) return;
      T *p1d=NULL, **p2d=NULL, ***p3d=NULL,****p4d=NULL;
      p4d = *matrix;
      p3d = p4d[0];
      delete [] p4d;
      this->Free<T>(&p3d, size1*size2, size3, size4);
      double freedMalloc = this->GetMemoryAmount<T***>(size1);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   //5d
   template<typename T> void Malloc(T****** matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5) const{
      if(*matrix!=NULL) return;
      if(size1*size2*size3*size4*size5<=0) return;
      double requiredMalloc = this->GetMemoryAmount<T****>(size1);
      this->CheckLimitHeap(requiredMalloc);

      T ****p4d=NULL, *****p5d=NULL;
      try{
         this->Malloc<T>(&p4d, size1*size2, size3, size4, size5);
         p5d = new T****[size1];
         if(p5d==NULL) throw MolDSException(this->errorMessageMallocFailure);

         for(size_t i=0;i<size1;i++){p5d[i]          = &p4d[i*size2];
         for(size_t j=0;j<size2;j++){p5d[i][j]       = &p4d[i*size2][j*size3];
         for(size_t k=0;k<size3;k++){p5d[i][j][k]    = &p4d[i*size2][j*size3][k*size4];
         for(size_t l=0;l<size4;l++){p5d[i][j][k][l] = &p4d[i*size2][j*size3][k*size4][l*size5];
         }}}}

         MallocerFreer::AddCurrentMalloced(requiredMalloc);
         *matrix = p5d;
      }
      catch(MolDSException ex){
         this->Free<T>(&p4d, size1*size2, size3, size4, size5);
         if(p5d!=NULL) delete[] p5d;
         throw ex;
      }
   }

   template<typename T> void Initialize(T***** matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5) const{
      for(size_t i=0;i<size1;i++) {
         this->Initialize<T>(matrix[i], size2, size3, size4, size5);
      }
   }

   template<typename T> void Free(T****** matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5) const{
      if(*matrix==NULL) return;
      T ****p4d=NULL, *****p5d=NULL;
      p5d = *matrix;
      p4d = p5d[0];
      delete [] p5d;
      this->Free<T>(&p4d, size1*size2, size3, size4, size5);
      double freedMalloc = this->GetMemoryAmount<T****>(size1);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   //6d
   template<typename T> void Malloc(T******* matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6) const{
      if(*matrix!=NULL) return;
      if(size1*size2*size3*size4*size5*size6<=0) return;
      double requiredMalloc = this->GetMemoryAmount<T*****>(size1);
      this->CheckLimitHeap(requiredMalloc);

      T *****p5d=NULL, ******p6d=NULL;
      try{
         this->Malloc<T>(&p5d, size1*size2, size3, size4, size5, size6);
         p6d = new T*****[size1];
         if(p6d==NULL) throw MolDSException(this->errorMessageMallocFailure);

         for(size_t i=0;i<size1;i++){p6d[i]             = &p5d[i*size2];
         for(size_t j=0;j<size2;j++){p6d[i][j]          = &p5d[i*size2][j*size3];
         for(size_t k=0;k<size3;k++){p6d[i][j][k]       = &p5d[i*size2][j*size3][k*size4];
         for(size_t l=0;l<size4;l++){p6d[i][j][k][l]    = &p5d[i*size2][j*size3][k*size4][l*size5];
         for(size_t m=0;m<size5;m++){p6d[i][j][k][l][m] = &p5d[i*size2][j*size3][k*size4][l*size5][m*size6];
         }}}}}

         MallocerFreer::AddCurrentMalloced(requiredMalloc);
         *matrix = p6d;
      }
      catch(MolDSException ex){
         this->Free<T>(&p5d, size1*size2, size3, size4, size5, size6);
         if(p6d!=NULL) delete[] p6d;
         throw ex;
      }
   }

   template<typename T> void Initialize(T****** matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6) const{
      for(size_t i=0;i<size1;i++) {
         this->Initialize<T>(matrix[i], size2, size3, size4, size5, size6);
      }
   }

   template<typename T> void Free(T******* matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6) const{
      if(*matrix==NULL) return;
      T *****p5d=NULL, ******p6d=NULL;
      p6d = *matrix;
      p5d = p6d[0];
      delete [] p6d;
      this->Free<T>(&p5d, size1*size2, size3, size4, size5, size6);
      double freedMalloc = this->GetMemoryAmount<T*****>(size1);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
   }

   //7d
   template<typename T> void Malloc(T******** matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6, size_t size7) const{
      if(*matrix!=NULL) return;
      if(size1*size2*size3*size4*size5*size6*size7<=0) return;
      double requiredMalloc = this->GetMemoryAmount<T******>(size1);
      this->CheckLimitHeap(requiredMalloc);

      T ******p6d=NULL, *******p7d=NULL;
      try{
         this->Malloc<T>(&p6d, size1*size2, size3, size4, size5, size6, size7);
         p7d = new T******[size1];
         if(p7d==NULL) throw MolDSException(this->errorMessageMallocFailure);

         for(size_t i=0;i<size1;i++){p7d[i]                = &p6d[i*size2];
         for(size_t j=0;j<size2;j++){p7d[i][j]             = &p6d[i*size2][j*size3];
         for(size_t k=0;k<size3;k++){p7d[i][j][k]          = &p6d[i*size2][j*size3][k*size4];
         for(size_t l=0;l<size4;l++){p7d[i][j][k][l]       = &p6d[i*size2][j*size3][k*size4][l*size5];
         for(size_t m=0;m<size5;m++){p7d[i][j][k][l][m]    = &p6d[i*size2][j*size3][k*size4][l*size5][m*size6];
         for(size_t n=0;n<size6;n++){p7d[i][j][k][l][m][n] = &p6d[i*size2][j*size3][k*size4][l*size5][m*size6][n*size7];
         }}}}}}

         MallocerFreer::AddCurrentMalloced(requiredMalloc);
         *matrix = p7d;
      }
      catch(MolDSException ex){
         this->Free<T>(&p6d, size1*size2, size3, size4, size5, size6, size7);
         if(p7d!=NULL) delete[] p7d;
         throw ex;
      }
   }

   template<typename T> void Initialize(T******* matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6, size_t size7) const{
      for(size_t i=0;i<size1;i++) {
         this->Initialize<T>(matrix[i], size2, size3, size4, size5, size6, size7);
      }
   }

   template<typename T> void Free(T******** matrix, size_t size1, size_t size2, size_t size3, size_t size4, size_t size5, size_t size6, size_t size7) const{
      if(*matrix==NULL) return;
      T ******p6d=NULL, *******p7d=NULL;
      p7d = *matrix;
      p6d = p7d[0];
      delete [] p7d;
      this->Free<T>(&p6d, size1*size2, size3, size4, size5, size6, size7);
      double freedMalloc = this->GetMemoryAmount<T******>(size1);
      MallocerFreer::AddCurrentMalloced(-1.0*freedMalloc);
      *matrix = NULL;
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
