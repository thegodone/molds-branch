//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   // 
// Copyright (C) 2011-2014 Katsuhiko Nishimra                             // 
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
#include<sstream>
#include<string>
#include<iostream>
#include<cctype>
#include<boost/format.hpp>
#include"Uncopyable.h"
#include"PrintController.h"
#include"MolDSException.h"
#include"MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
using namespace std;
namespace MolDS_base{

PrintController::PrintController(){
   this->canOutputLogs = true;
   //this->OutputLog("printController is created.\n");
}

PrintController::~PrintController(){
   //this->OutputLog("printController is destructed.\n";
}

void PrintController::OutputLog(string log) const{
   if(this->canOutputLogs){
      int mpiHeadRank = MolDS_mpi::MpiProcess::GetInstance()->GetHeadRank();
#ifndef MOLDS_DBG
      if(MolDS_mpi::MpiProcess::GetInstance()->GetRank()!=mpiHeadRank){return;}
#endif
      bool endl = false;
      string::reverse_iterator iter;
      for(iter = log.rbegin(); iter != log.rend(); iter++){
         if(*iter == '\n'){
            string::iterator fwditer = iter.base();
            log.erase(--fwditer);
            endl = true;
            break;
         }
         else if(*iter != '\0'){
            break;
         }
      }
      cout << log;
      if(endl){cout << std::endl;}
   }
}
}
