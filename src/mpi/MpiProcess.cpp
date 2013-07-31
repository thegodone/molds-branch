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
#include<sstream>
#include<math.h>
#include<string>
#include<stdexcept>
#include<boost/format.hpp>
#include"../base/Uncopyable.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/MallocerFreer.h"
#include"MpiProcess.h"
using namespace std;
namespace MolDS_mpi{

MpiProcess* MpiProcess::mpiProcess = NULL;

MpiProcess::MpiProcess(){
}

MpiProcess::MpiProcess(int argc, char *argv[]){
   this->environment  = new boost::mpi::environment(argc, argv);
   this->communicator = new boost::mpi::communicator();
   this->messageLimit = INT_MAX;
}

MpiProcess::~MpiProcess(){
   delete this->environment;
   delete this->communicator;
}

void MpiProcess::CreateInstance(int argc, char *argv[]){
   if(mpiProcess != NULL){
      // ToDo: error
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
      //mpiProcess = new MpiProcess();
      // ToDo: error
   }
   return mpiProcess;
}

}





