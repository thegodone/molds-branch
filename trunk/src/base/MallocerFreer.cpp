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
#include<vector>
#include<stdexcept>
#include<boost/format.hpp>
#include"Enums.h"
#include"Uncopyable.h"
#include"PrintController.h"
#include"MolDSException.h"
#include"MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
#include"EularAngle.h"
#include"Parameters.h"
using namespace std;

namespace MolDS_base{

MallocerFreer* MallocerFreer::mallocerFreer = NULL;
double MallocerFreer::currentMalloced = 0.0;
double MallocerFreer::maxMalloced = 0.0;
const double MallocerFreer::byte2MByte = 1e-06;

MallocerFreer::MallocerFreer(){
   this->errorMessageMallocFailure = "Error in base::MallocerFreer::Malloc: Malloc failure...\n";
   this->errorMessageReachHeapLimit = "Error in base::MallocerFreer::Malloc: Reaches limit of heap. Change the \"limit_heap\" option in the memory-directive, machine you using, or your study!!!\n";
   this->messageMemoryUsage =        "\tSummary for memory usage:\n";
   this->messageMemoryMaxHeap =      "\t\tMax Heap: ";
   this->messageMemoryLeakedHeap =   "\t\tCurrent Heap(Leaked): ";
   this->messageMemoryCurrentHeap =  "\t\tCurrent Heap:  ";
   this->messageMemoryRequiredHeap = "\t\tRequired Heap: ";
   this->messageMemoryLimitHeap =    "\t\tHeap Limit:    ";
   this->messageMByte = "[MB].\n";
}

MallocerFreer::~MallocerFreer(){
   this->OutputMemoryUsage();
}

void MallocerFreer::CheckLimitHeap(double requiredMalloc) const{
   double limit = Parameters::GetInstance()->GetLimitHeapMemory();
   if(limit < (MallocerFreer::currentMalloced + requiredMalloc)*MallocerFreer::byte2MByte){
      stringstream ss;
      ss << this->errorMessageReachHeapLimit;
      ss << this->messageMemoryLimitHeap    << limit                                                    << this->messageMByte;
      ss << this->messageMemoryCurrentHeap  << MallocerFreer::currentMalloced*MallocerFreer::byte2MByte << this->messageMByte;
      ss << this->messageMemoryRequiredHeap << requiredMalloc                *MallocerFreer::byte2MByte << this->messageMByte;
      throw MolDSException(ss.str());
   }
}

void MallocerFreer::OutputMemoryUsage() const{
   this->OutputLog(this->messageMemoryUsage);
   this->OutputLog(boost::format("%s%lf%s") % this->messageMemoryMaxHeap.c_str() 
                                            % (MallocerFreer::maxMalloced*MallocerFreer::byte2MByte)
                                            % this->messageMByte.c_str());
   this->OutputLog(boost::format("%s%lf%s") % this->messageMemoryLeakedHeap.c_str() 
                                            % (MallocerFreer::currentMalloced*MallocerFreer::byte2MByte)
                                            % this->messageMByte.c_str());
}

MallocerFreer* MallocerFreer::GetInstance(){
   if(mallocerFreer == NULL){
      mallocerFreer = new MallocerFreer();
   }
   return mallocerFreer;
}

void MallocerFreer::DeleteInstance(){
   if(mallocerFreer != NULL){
      delete mallocerFreer;
   }
   mallocerFreer = NULL;
}

void MallocerFreer::AddCurrentMalloced(double amount){
   #pragma omp atomic 
   MallocerFreer::currentMalloced += amount;
   if(0 < amount){
      #pragma omp flush(MallocerFreer::maxMalloced, MallocerFreer::currentMalloced)
      if(MallocerFreer::maxMalloced < MallocerFreer::currentMalloced){
         #pragma omp critical
         {   
            if(MallocerFreer::maxMalloced < MallocerFreer::currentMalloced){
               MallocerFreer::maxMalloced = MallocerFreer::currentMalloced;
            }   
         }   
      }   
   }   
}

}
