//************************************************************************//
// Copyright (C) 2011-2013 Mikiya Fujii                                   // 
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
#ifndef INCLUDED_ASYNCCOMMUNICATOR
#define INCLUDED_ASYNCCOMMUNICATOR
#include<boost/thread.hpp>
#include<boost/thread/condition.hpp>
#include<boost/bind.hpp>
#define NON_USED 0
namespace MolDS_mpi{

class AsyncCommunicator{
public:
   AsyncCommunicator();
   ~AsyncCommunicator();
   template<typename T> void Run(int passingTimes){
      int mpiRank = MolDS_mpi::MpiProcess::GetInstance()->GetRank();
      while(0<passingTimes){
         sleep(0.1);
         boost::mutex::scoped_lock lk(this->stateGuard);
         try{
            DataInfo dInfo = this->dataQueue.FrontPop();
            if(dInfo.mpiFuncType == MolDS_base::Send){
               MolDS_mpi::MpiProcess::GetInstance()->Send(dInfo.dest,
                                                          dInfo.tag,
                                                          reinterpret_cast<T*>(dInfo.vectorPtr), 
                                                          dInfo.num);
            }
            else if(dInfo.mpiFuncType == MolDS_base::Recv){
               MolDS_mpi::MpiProcess::GetInstance()->Recv(dInfo.source,
                                                          dInfo.tag,
                                                          reinterpret_cast<T*>(dInfo.vectorPtr), 
                                                          dInfo.num);
            }
            else if(dInfo.mpiFuncType == MolDS_base::Broadcast){
               MolDS_mpi::MpiProcess::GetInstance()->Broadcast(reinterpret_cast<T*>(dInfo.vectorPtr), 
                                                               dInfo.num, 
                                                               dInfo.source);
            }
            else{
               std::stringstream ss;
               ss << "non valid mpi function type\n";
               MolDS_base::MolDSException ex(ss.str());
               throw ex;
            }
            this->stateChange.notify_all();
            passingTimes--;
         }
         catch(MolDS_base::MolDSException ex){
            if(ex.HasKey(MolDS_base::EmptyQueue)){
               this->stateChange.wait(lk);
               continue;
            }
            else{
               throw ex;
            }
         }
      }
   }

   template<typename T> void SetSentVector(T* vector, 
                                           intptr_t num, 
                                           int dest,
                                           int tag){
      int source = NON_USED;
      MolDS_base::MpiFunctionType mpiFuncType = MolDS_base::Send;
      this->SetVector(vector, num, source, dest, tag, mpiFuncType);
   }

   template<typename T> void SetRecvedVector(T* vector, 
                                             intptr_t num, 
                                             int source, 
                                             int tag){
      int dest   = NON_USED;
      MolDS_base::MpiFunctionType mpiFuncType = MolDS_base::Recv;
      this->SetVector(vector, num, source, dest, tag, mpiFuncType);
   }

   template<typename T> void SetBroadcastedVector(T* vector, intptr_t num, int root){
      int source = root;
      int dest   = NON_USED;
      int tag    = NON_USED;
      MolDS_base::MpiFunctionType mpiFuncType = MolDS_base::Broadcast;
      this->SetVector(vector, num, source, dest, tag, mpiFuncType);
   }

private:
   struct DataInfo{intptr_t vectorPtr; 
                   intptr_t num; 
                   int source; 
                   int dest; 
                   int tag;
                   MolDS_base::MpiFunctionType mpiFuncType;};
   boost::mutex     stateGuard;
   boost::condition stateChange;
   MolDS_base_containers::ThreadSafeQueue<DataInfo> dataQueue;
   template<typename T> void SetVector(T* vector, 
                                       intptr_t num, 
                                       int source, 
                                       int dest, 
                                       int tag,
                                       MolDS_base::MpiFunctionType mpiFuncType){
      boost::mutex::scoped_lock lk(this->stateGuard);
      DataInfo dInfo = {reinterpret_cast<intptr_t>(vector), num, source, dest, tag, mpiFuncType};
      this->dataQueue.Push(dInfo);
      this->stateChange.notify_all();
   }
};

}
#endif

