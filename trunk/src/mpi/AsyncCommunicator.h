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
   template<typename T> void Run(){
      int mpiRank = MolDS_mpi::MpiProcess::GetInstance()->GetRank();
      while(true){
         boost::mutex::scoped_lock lk(this->stateGuard);
         try{
            MessageInfo mInfo = this->messageQueue.FrontPop();
            if(mInfo.mpiFuncType == MolDS_base::Send){
               MolDS_mpi::MpiProcess::GetInstance()->Send(mInfo.dest,
                                                          mInfo.tag,
                                                          reinterpret_cast<T*>(mInfo.vectorPtr), 
                                                          mInfo.num);
            }
            else if(mInfo.mpiFuncType == MolDS_base::Recv){
               MolDS_mpi::MpiProcess::GetInstance()->Recv(mInfo.source,
                                                          mInfo.tag,
                                                          reinterpret_cast<T*>(mInfo.vectorPtr), 
                                                          mInfo.num);
            }
            else if(mInfo.mpiFuncType == MolDS_base::Broadcast){
               MolDS_mpi::MpiProcess::GetInstance()->Broadcast(reinterpret_cast<T*>(mInfo.vectorPtr), 
                                                               mInfo.num, 
                                                               mInfo.source);
            }
            else{
               std::stringstream ss;
               ss << "non valid mpi function type\n";
               MolDS_base::MolDSException ex(ss.str());
               throw ex;
            }
            this->stateChange.notify_all();
         }
         catch(MolDS_base::MolDSException ex){
            if(ex.HasKey(MolDS_base::EmptyQueue && this->hasAllMessagesSet)){
               break;
            }
            else if(ex.HasKey(MolDS_base::EmptyQueue && !this->hasAllMessagesSet)){
               this->stateChange.wait(lk);
               continue;
            }
            else{
               throw ex;
            }
         }
      }
   }

   template<typename T> void SetSentMessage(T* vector, 
                                            molds_mpi_int num, 
                                            int dest,
                                            int tag){
      int source = NON_USED;
      MolDS_base::MpiFunctionType mpiFuncType = MolDS_base::Send;
      this->SetMessage(vector, num, source, dest, tag, mpiFuncType);
   }

   template<typename T> void SetRecvedMessage(T* vector, 
                                              molds_mpi_int num, 
                                              int source, 
                                              int tag){
      int dest   = NON_USED;
      MolDS_base::MpiFunctionType mpiFuncType = MolDS_base::Recv;
      this->SetMessage(vector, num, source, dest, tag, mpiFuncType);
   }

   template<typename T> void SetBroadcastedMessage(T* vector, molds_mpi_int num, int root){
      int source = root;
      int dest   = NON_USED;
      int tag    = NON_USED;
      MolDS_base::MpiFunctionType mpiFuncType = MolDS_base::Broadcast;
      this->SetMessage(vector, num, source, dest, tag, mpiFuncType);
   }

   void Finalize();

private:
   struct MessageInfo{intptr_t vectorPtr; 
                      molds_mpi_int num; 
                      int source; 
                      int dest; 
                      int tag;
                      MolDS_base::MpiFunctionType mpiFuncType;};
   boost::mutex     stateGuard;
   boost::condition stateChange;
   bool hasAllMessagesSet;
   MolDS_base_containers::ThreadSafeQueue<MessageInfo> messageQueue;
   template<typename T> void SetMessage(T* vector, 
                                        molds_mpi_int num, 
                                        int source, 
                                        int dest, 
                                        int tag,
                                        MolDS_base::MpiFunctionType mpiFuncType){
      boost::mutex::scoped_lock lk(this->stateGuard);
      MessageInfo mInfo = {reinterpret_cast<intptr_t>(vector), num, source, dest, tag, mpiFuncType};
      this->messageQueue.Push(mInfo);
      this->stateChange.notify_all();
   }
};

}
#endif

