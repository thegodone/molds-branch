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
#include<execinfo.h>
#include<string>
#include<stdexcept>
#include<boost/format.hpp>
#include"MolDSException.h"
using namespace std;
namespace MolDS_base{
MolDSException::MolDSException(string cause) : domain_error(cause), backtraceSize(0){
   this->GetBacktrace(80);
}

MolDSException::MolDSException(const boost::format& cause) : domain_error(cause.str()), backtraceSize(0){
   this->GetBacktrace(80);
}

void MolDSException::GetBacktrace(int bufsize){
   backtracePtr.reset(new void*[bufsize]);
   this->backtraceSize = backtrace(this->backtracePtr.get(), bufsize);
   if(this->backtraceSize==bufsize){
      GetBacktrace(bufsize*2);
   }
}

void MolDSException::PrintBacktrace(){
   if(this->backtraceSize <= 0){
      return;
   }
   cout<<"backtrace:" << endl;
   backtrace_symbols_fd(this->backtracePtr.get(), this->backtraceSize, 1);
}
}

