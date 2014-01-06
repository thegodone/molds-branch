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
#include<omp.h>
#include<boost/format.hpp>
#include"../base/Uncopyable.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/MallocerFreer.h"
#include"MpiInt.h"
#include"MpiProcess.h"
using namespace std;
namespace MolDS_mpi{

MpiProcess* MpiProcess::mpiProcess = NULL;
string      MpiProcess::errorMessageCreateInstanceDuplicate
               = "Error in mpi::MpiProcess::CreateInstance: mpiProcess has been already created, namely duplication error.\n";
string      MpiProcess::errorMessageGetInstanceNULL
               = "Error in mpi::MpiProcess::GetInstance: mpiProcess is NULL.\n";

MpiProcess::MpiProcess(){
}

MpiProcess::MpiProcess(int argc, char *argv[]){
   this->environment  = new boost::mpi::environment(argc, argv);
   this->communicator = new boost::mpi::communicator();
   this->messageLimit = INT_MAX;
   this->mpiConsumingTime=0.0;
   this->mpiConsumingTimeSend=0.0;
   this->mpiConsumingTimeRecv=0.0;
   this->mpiConsumingTimeBrodCast=0.0;
   this->mpiConsumingTimeAllReduce=0.0;
   this->SetMessages();
}

MpiProcess::~MpiProcess(){
   /*
   int rank = this->GetRank();
   printf("\nrnk:%d mpiconsumingtime          = %e [s]\n",rank, this->mpiConsumingTime);
   printf("\nrnk:%d mpiconsumingtimeSend      = %e [s]\n",rank, this->mpiConsumingTimeSend);
   printf("\nrnk:%d mpiconsumingtimeRecv      = %e [s]\n",rank, this->mpiConsumingTimeRecv);
   printf("\nrnk:%d mpiconsumingtimeBroadcast = %e [s]\n",rank, this->mpiConsumingTimeBrodCast);
   printf("\nrnk:%d mpiconsumingtimeAllReduce = %e [s]\n",rank, this->mpiConsumingTimeAllReduce);
   */
   delete this->environment;
   delete this->communicator;
}

void MpiProcess::CreateInstance(int argc, char *argv[]){
   if(mpiProcess != NULL){
      std::stringstream ss;
      ss << errorMessageCreateInstanceDuplicate;
      MolDS_base::MolDSException ex(ss.str());
      throw ex;
   }
   mpiProcess = new MpiProcess(argc, argv);
}

void MpiProcess::DeleteInstance(){
   if(mpiProcess != NULL){
      delete mpiProcess; 
   }
   mpiProcess = NULL;
}

MpiProcess* MpiProcess::GetInstance(){
   if(mpiProcess == NULL){
      std::stringstream ss;
      ss << errorMessageGetInstanceNULL;
      MolDS_base::MolDSException ex(ss.str());
      throw ex;
   }
   return mpiProcess;
}

void MpiProcess::Barrier(){this->communicator->barrier();}

void MpiProcess::SetMessages(){
   this->errorMessageSplitMessageElemLimNegative
      = "Error in mpi::MpiProcess::SplitMessage2Chunks: elementsLimit is negative. \nelementsLimit=";
   this->errorMessageSplitMessageNumChnkNegative
      = "Error in mpi::MpiProcess::SplitMessage2Chunks: numChunks is negative. \nnumChunks=";
   this->errorMessageSplitMessageTagBaseNegative
      = "Error in mpi::MpiProcess::SplitMessage2Chunks: tagBase is negative. \ntagBase=";
   this->errorMessageSplitMessageRemainingNegative
      = "Error in mpi::MpiProcess::SplitMessage2Chunks: remaining is negative. \nremaining=";
}

}





