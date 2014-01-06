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
#include<boost/format.hpp>
#include"../base/Enums.h"
#include"../base/Uncopyable.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/containers/ThreadSafeQueue.h"
#include"../base/MallocerFreer.h"
#include"MpiInt.h"
#include"MpiProcess.h"
#include"AsyncCommunicator.h"
using namespace std;
namespace MolDS_mpi{
AsyncCommunicator::AsyncCommunicator(){
   this->hasAllMessagesSet=false;
}
AsyncCommunicator::~AsyncCommunicator(){}
void AsyncCommunicator::Finalize(){
   boost::mutex::scoped_lock lk(this->stateGuard);
   this->hasAllMessagesSet = true;
   this->stateChange.notify_all();
}
}




