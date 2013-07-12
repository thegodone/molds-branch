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
#ifndef INCLUDED_MPIPROCESS
#define INCLUDED_MPIPROCESS
#include<limits.h>
#include<boost/mpi.hpp>
namespace MolDS_mpi{
// MpiProcess is singleton
class MpiProcess: private MolDS_base::Uncopyable{
public:
   static void        CreateInstance(int argc, char *argv[]);
   static void        DeleteInstance();
   static MpiProcess* GetInstance();
   int GetRank() const{return this->communicator->rank();}
   int GetSize() const{return this->communicator->size();}
   template<typename T> void Send(int dest, int tag, const T* values, intptr_t num) const{
      std::vector<Chunk> chunks;
      this->SplitMessage2Chunks(chunks, tag, values, num);
      for(intptr_t i=0; i<chunks.size(); i++){
         this->communicator->send(dest, chunks[i].tag, &values[chunks[i].first], chunks[i].num);
      }
   }
   template<typename T> void Recv(int source, int tag, T* values, intptr_t num) const{
      std::vector<Chunk> chunks;
      this->SplitMessage2Chunks(chunks, tag, values, num);
      for(intptr_t i=0; i<chunks.size(); i++){
         this->communicator->recv(source, chunks[i].tag, &values[chunks[i].first], chunks[i].num);
      }
   }
   template<typename T> void Broadcast(T* values, intptr_t num, int root) const{
      std::vector<Chunk> chunks;
      intptr_t tag=0;
      this->SplitMessage2Chunks(chunks, tag, values, num);
      for(intptr_t i=0; i<chunks.size(); i++){
         broadcast(*this->communicator, &values[chunks[i].first], chunks[i].num, root);
      }
   }
private:
   static MpiProcess* mpiProcess;
   MpiProcess();
   MpiProcess(int argc, char *argv[]);
   ~MpiProcess();
   boost::mpi::environment*  environment;
   boost::mpi::communicator* communicator;
   struct Chunk{int tag; intptr_t first; int num;};
   double messageLimit;
   template<typename T> void SplitMessage2Chunks(std::vector<Chunk>& chunks, const int origianlTag, T* values, intptr_t num) const{
      if(this->messageLimit < static_cast<double>(sizeof(T))*static_cast<double>(num) ){
         int  elementsLimit = static_cast<intptr_t>(messageLimit/sizeof(T));
         intptr_t numChunks = num/elementsLimit;
         int      remaining = num%elementsLimit;
         if(0 < remaining){
            numChunks++;
         }
         int tagBase = origianlTag*numChunks;
         for(intptr_t i=0; i<numChunks; i++){
            int tag = tagBase+i;
            Chunk chunk = {tag, i*elementsLimit, elementsLimit};
            chunks.push_back(chunk);
         }
         if(0 < remaining){
            chunks[numChunks-1].num = remaining;
         }
      }
      else{
         Chunk chunk = {origianlTag, 0, num};
         chunks.push_back(chunk);
      }
   }
};

}
#endif

