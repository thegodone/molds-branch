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
#include<sstream>
#include<string>
#include<stdexcept>
#include<iostream>
#include<boost/format.hpp>
#include<boost/serialization/map.hpp>
#include<boost/archive/text_iarchive.hpp>
#include<boost/archive/text_oarchive.hpp>
#include"MolDSException.h"
#include"Enums.h"
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

template<>
int MolDSException::GetKeyValue(int key){
   return intKeyValueMap[key];
}

/*
template<>
other MolDSException::GetKeyValue(int key){
   return otherKeyValueMap[key];
}
*/

template<>
void MolDSException::SetKeyValue(int key, int value){
   intKeyValueMap[key]=value;
}

/*
template<>
void MolDSException::SetKeyValue(int key, other value){
   otherKeyValueMap.insert(otherKeyValueMap::pair(key,value));
}
*/

bool MolDSException::HasKey(int key){
   if(!(intKeyValueMap.find(key)==intKeyValueMap.end())) return true;
   //if(otherKeyValueMap.find(key)!=otherKeyValueMap::end) return true;
   return false;
}

const char* MolDSException::what() const throw(){
   static string str;
   stringstream ss;
   ss << domain_error::what() << "key value pairs:";
   for(intKeyValueMap_t::const_iterator i = intKeyValueMap.begin(); i != intKeyValueMap.end(); i++){
      ss << endl << '\t' << ExceptionKeyStr(i->first) << ":" << i->second;
   }
   str = ss.str();
   return str.c_str();
}

template<class Archive>
void MolDSException::serialize(Archive& ar, const unsigned int ver){
   ar & intKeyValueMap;
   // ar & otherKeyValueMap;

   ar & backtraceSize;
   std::cerr << "backtraceSize:" << backtraceSize << std::endl;
   if(!Archive::is_saving::value){
      backtracePtr.reset(new void*[backtraceSize]);
   }
   for(int i; i<backtraceSize; i++){
      if(Archive::is_saving::value){
         intptr_t p = reinterpret_cast<intptr_t>(backtracePtr[i]);
         ar & p;
         std::cerr << "in: " << p << std::endl;
      }
      else{
         intptr_t p;
         ar & p;
         std::cerr << "out:" << p << std::endl;
         backtracePtr[i]=reinterpret_cast<void*>(p);
      }
   }
}

void MolDSException::Serialize(std::ostream& os){
   boost::archive::text_oarchive oa(os);
   std::string what = domain_error::what();
   std::cerr << "what:" << what << std::endl;
   oa << what;
   oa << (*this);
}

MolDSException MolDSException::Deserialize(std::istream& is){
   boost::archive::text_iarchive ia(is);
   std::string what;
   ia >> what;
   MolDSException e(what);
   ia >> e;
   return e;
}
}
