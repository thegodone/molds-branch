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
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<math.h>
#include<string>
#include<vector>
#include<stdexcept>
#include<boost/format.hpp>
#include"PrintController.h"
#include"MolDSException.h"
#include"Uncopyable.h"
#include"Enums.h"
#include"MallocerFreer.h"
#include"EularAngle.h"
#include"Parameters.h"
using namespace std;

namespace MolDS_base{

MallocerFreer* MallocerFreer::mallocerFreer = NULL;
double MallocerFreer::currentMalloced = 0.0;
double MallocerFreer::maxMalloced = 0.0;

MallocerFreer::MallocerFreer(){
   this->errorMessageMallocFailure = "Error in base::MallocFreer::Malloc: Malloc failure...\n";
   this->errorMessageReachHeapLimit = "Error in base::MallocFreer::Malloc: Reaches limit of heap. Change the \"limit_heap\" option in the memory-directive, machine you using, or your study!!!\n";
   this->messageMemoryUsage = "\tSummary for memory usage:\n";
   this->messageMemoryMaxHeap = "\t\tMax Heap: ";
   this->messageMemoryCurrentHeap = "\t\tCurrent Heap(Leaked): ";
   this->messageMByte = "[MB].\n";
}

MallocerFreer::~MallocerFreer(){
   this->OutputMemoryUsage();
}

void MallocerFreer::CheckLimitHeap(double wannaMalloc) const{
   double limit = Parameters::GetInstance()->GetLimitHeapMemory();
   if(limit < (MallocerFreer::currentMalloced + wannaMalloc)/pow(10.0,6.0)){
      throw MolDSException(this->errorMessageReachHeapLimit);
   }
}

void MallocerFreer::OutputMemoryUsage() const{
   this->OutputLog(this->messageMemoryUsage);
   this->OutputLog((boost::format("%s%lf%s") % this->messageMemoryMaxHeap.c_str() 
                                             % (MallocerFreer::maxMalloced/pow(10.0,6.0))
                                             % this->messageMByte.c_str()).str());
   this->OutputLog((boost::format("%s%lf%s") % this->messageMemoryCurrentHeap.c_str() 
                                             % (MallocerFreer::currentMalloced/pow(10.0,6.0))
                                             % this->messageMByte.c_str()).str());
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
   #pragma omp critical
   {
      MallocerFreer::currentMalloced += amount;
      if(MallocerFreer::maxMalloced < MallocerFreer::currentMalloced){
         MallocerFreer::maxMalloced = MallocerFreer::currentMalloced;
      }
   }
}

void MallocerFreer::SubtCurrentMalloced(double amount){
   #pragma omp critical
   MallocerFreer::currentMalloced -= amount;
}

}
