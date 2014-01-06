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
#ifndef INCLUDED_MPIPROCESS
#define INCLUDED_MPIPROCESS
#include<limits.h>
#include<omp.h>
#include<boost/mpi.hpp>
namespace MolDS_mpi{
// MpiProcess is singleton
class MpiProcess: private MolDS_base::Uncopyable{
public:
   static void        CreateInstance(int argc, char *argv[]);
   static void        DeleteInstance();
   static MpiProcess* GetInstance();
   int GetHeadRank() const{return 0;}
   int GetRank() const{return this->communicator->rank();}
   int GetSize() const{return this->communicator->size();}
   //template<typename T> void Send(int dest, int tag, const T* values, molds_mpi_int num) const{
   template<typename T> void Send(int dest, int tag, const T* values, molds_mpi_int num) {
      double startTime=0.0;
      double endTime=0.0;
      std::vector<Chunk> chunks;
      this->SplitMessage2Chunks(chunks, tag, values, num);
      for(molds_mpi_int i=0; i<chunks.size(); i++){
         startTime = omp_get_wtime();
         this->communicator->send(dest, chunks[i].tag, &values[chunks[i].first], chunks[i].num);
         endTime = omp_get_wtime();
         this->mpiConsumingTime     += endTime - startTime;
         this->mpiConsumingTimeSend += endTime - startTime;
      }
   }
   //template<typename T> void Recv(int source, int tag, T* values, molds_mpi_int num) const{
   template<typename T> void Recv(int source, int tag, T* values, molds_mpi_int num) {
      double startTime=0.0;
      double endTime=0.0;
      std::vector<Chunk> chunks;
      this->SplitMessage2Chunks(chunks, tag, values, num);
      for(molds_mpi_int i=0; i<chunks.size(); i++){
         startTime = omp_get_wtime();
         this->communicator->recv(source, chunks[i].tag, &values[chunks[i].first], chunks[i].num);
         endTime = omp_get_wtime();
         this->mpiConsumingTime     += endTime - startTime;
         this->mpiConsumingTimeRecv += endTime - startTime;
      }
   }
   //template<typename T> void Broadcast(T* values, molds_mpi_int num, int root) const{
   template<typename T> void Broadcast(T* values, molds_mpi_int num, int root){
      double startTime=0.0;
      double endTime=0.0;
      std::vector<Chunk> chunks;
      molds_mpi_int tag=0;
      this->SplitMessage2Chunks(chunks, tag, values, num);
      for(molds_mpi_int i=0; i<chunks.size(); i++){
         startTime = omp_get_wtime();
         broadcast(*this->communicator, &values[chunks[i].first], chunks[i].num, root);
         endTime = omp_get_wtime();
         this->mpiConsumingTime         += endTime - startTime;
         this->mpiConsumingTimeBrodCast += endTime - startTime;
      }
   }
   template<typename T, typename Op> void Reduce(const T* inValues, molds_mpi_int num, T* outValues, Op op, int root) const{
      std::vector<Chunk> chunks;
      molds_mpi_int tag=0;
      this->SplitMessage2Chunks(chunks, tag, inValues, num);
      for(molds_mpi_int i=0; i<chunks.size(); i++){
         reduce(*this->communicator, &inValues[chunks[i].first], chunks[i].num, &outValues[chunks[i].first], op, root);
      }
   }
   //template<typename T, typename Op> void AllReduce(const T* inValues, molds_mpi_int num, T* outValues, Op op) const{
   template<typename T, typename Op> void AllReduce(const T* inValues, molds_mpi_int num, T* outValues, Op op){
      double startTime=0.0;
      double endTime=0.0;
      std::vector<Chunk> chunks;
      molds_mpi_int tag=0;
      this->SplitMessage2Chunks(chunks, tag, inValues, num);
      for(molds_mpi_int i=0; i<chunks.size(); i++){
         startTime = omp_get_wtime();
         all_reduce(*this->communicator, &inValues[chunks[i].first], chunks[i].num, &outValues[chunks[i].first], op);
         endTime = omp_get_wtime();
         this->mpiConsumingTime          += endTime - startTime;
         this->mpiConsumingTimeAllReduce += endTime - startTime;
      }
   }
   //template<typename T, typename Op> void AllReduce(T* values, molds_mpi_int num, Op op) const{
   template<typename T, typename Op> void AllReduce(T* values, molds_mpi_int num, Op op){
      double* tmpValues=NULL;
      try{
         MolDS_base::MallocerFreer::GetInstance()->Malloc<double>(&tmpValues, num);
         this->AllReduce(values, num, tmpValues, op);
         for(molds_mpi_int i=0; i<num; i++){
            values[i] = tmpValues[i];
         }
      }
      catch(MolDS_base::MolDSException ex){
         MolDS_base::MallocerFreer::GetInstance()->Free<double>(&tmpValues, num);
         throw ex;
      }
      MolDS_base::MallocerFreer::GetInstance()->Free<double>(&tmpValues, num);
   }
   void Barrier();
private:
   static MpiProcess* mpiProcess;
   MpiProcess();
   MpiProcess(int argc, char *argv[]);
   ~MpiProcess();
   static std::string errorMessageCreateInstanceDuplicate;
   static std::string errorMessageGetInstanceNULL;
   std::string errorMessageSplitMessageElemLimNegative;
   std::string errorMessageSplitMessageNumChnkNegative;
   std::string errorMessageSplitMessageTagBaseNegative;
   std::string errorMessageSplitMessageRemainingNegative;
   boost::mpi::environment*  environment;
   boost::mpi::communicator* communicator;
   double messageLimit;
   struct Chunk{int tag; molds_mpi_int first; int num;};
   void SetMessages();
   template<typename T> void SplitMessage2Chunks(std::vector<Chunk>& chunks, const int origianlTag, T* values, molds_mpi_int num) const{
      if(this->messageLimit < static_cast<double>(sizeof(T))*static_cast<double>(num) ){
         int elementsLimit = static_cast<int>(messageLimit/sizeof(T));
         int numChunks     = num/elementsLimit;
         int remaining     = num%elementsLimit;
         if(0 < remaining){
            numChunks++;
         }
         int tagBase = origianlTag*numChunks;
         if(elementsLimit < 0){
            std::stringstream ss;
            ss << this->errorMessageSplitMessageElemLimNegative << elementsLimit << std::endl;
            MolDS_base::MolDSException ex(ss.str());
            throw ex;
         }
         if(numChunks < 0){
            std::stringstream ss;
            ss << this->errorMessageSplitMessageNumChnkNegative << numChunks << std::endl;
            MolDS_base::MolDSException ex(ss.str());
            throw ex;
         }
         if(remaining < 0){
            std::stringstream ss;
            ss << this->errorMessageSplitMessageRemainingNegative << remaining << std::endl;
            MolDS_base::MolDSException ex(ss.str());
            throw ex;
         }
         if(tagBase < 0){
            std::stringstream ss;
            ss << this->errorMessageSplitMessageTagBaseNegative << tagBase << std::endl;
            MolDS_base::MolDSException ex(ss.str());
            throw ex;
         }
         for(int i=0; i<numChunks; i++){
            int tag             = tagBase+i;
            molds_mpi_int first = i*static_cast<molds_mpi_int>(elementsLimit);
            Chunk chunk         = {tag, first, elementsLimit};
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
   double mpiConsumingTime;
   double mpiConsumingTimeSend;
   double mpiConsumingTimeRecv;
   double mpiConsumingTimeBrodCast;
   double mpiConsumingTimeAllReduce;
};

}
#endif

