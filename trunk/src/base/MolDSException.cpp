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
#include<execinfo.h>
#include<sstream>
#include<string>
#include<stdexcept>
#include<iostream>
#include<cctype>
#include<boost/format.hpp>
#include<boost/serialization/map.hpp>
#include<boost/archive/text_iarchive.hpp>
#include<boost/archive/text_oarchive.hpp>
#include<boost/scoped_ptr.hpp>
#include"Enums.h"
#include"MolDSException.h"
using namespace std;

namespace boost { namespace serialization {
template<class Archive>
void save_construct_data(Archive& ar, const MolDS_base::MolDSException* t, const unsigned int){
   string what(t->What());
   ar << what;
}

template<class Archive>
void load_construct_data(Archive& ar, MolDS_base::MolDSException* t, const unsigned int){
   string what;
   ar >> what;
   ::new(t)MolDS_base::MolDSException(what.c_str());
}
}}

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
      this->GetBacktrace(bufsize*2);
   }
}

template<>
int MolDSException::GetKeyValue(int key){
   return this->intKeyValueMap[key];
}

/*
template<>
other MolDSException::GetKeyValue(int key){
   return this->otherKeyValueMap[key];
}
*/

template<>
void MolDSException::SetKeyValue(int key, int value){
   this->intKeyValueMap[key]=value;
}

/*
template<>
void MolDSException::SetKeyValue(int key, other value){
   otherKeyValueMap.insert(otherKeyValueMap::pair(key,value));
}
*/

bool MolDSException::HasKey(int key){
   if(!(this->intKeyValueMap.find(key)==this->intKeyValueMap.end())) return true;
   //if(this->otherKeyValueMap.find(key)!=this->otherKeyValueMap::end) return true;
   return false;
}

const char* MolDSException::what() const throw(){
   static string str;
   stringstream ss;
   ss << domain_error::what();

   ss << "\nkey value pairs:";
   for(intKeyValueMap_t::const_iterator i = this->intKeyValueMap.begin();
       i != this->intKeyValueMap.end(); i++){
      ss << endl << '\t' << ExceptionKeyStr(i->first) << ":" << i->second;
   }
   char** backtraceSymbols = backtrace_symbols(this->backtracePtr.get(),
                                               this->backtraceSize);
   if(backtraceSymbols == NULL){
      ss << "\nCannot Get Backtraces!";
   }
   else{
      ss << "\nbacktrace:";
      for(int i = 0; i < this->backtraceSize; i++){
         ss << "\n\t" << backtraceSymbols[i];
      }
   }
   free(backtraceSymbols);

   if(this->nextException.get() != NULL){
      ss << "\n--- Next Exception ---\n";
      ss << this->nextException->what();
   }
   str = ss.str();
   return str.c_str();
}

template<class Archive>
void MolDSException::serialize(Archive& ar, const unsigned int ver){
   ar & this->intKeyValueMap;
   // ar & this->otherKeyValueMap;

   ar & this->backtraceSize;
   if(!Archive::is_saving::value){
      this->backtracePtr.reset(new void*[this->backtraceSize]);
   }
   for(int i = 0; i < this->backtraceSize; i++){
      if(Archive::is_saving::value){
         intptr_t p = reinterpret_cast<intptr_t>(this->backtracePtr[i]);
         ar & p;
      }
      else{
         intptr_t p;
         ar & p;
         this->backtracePtr[i]=reinterpret_cast<void*>(p);
      }
   }

   bool hasnext = this->nextException.get() != NULL;
   ar & hasnext;
   MolDSException* pe = NULL;
   if(Archive::is_saving::value){
      pe = this->nextException.get();
   }
   if(hasnext){
      ar & pe;
   }
   if(!Archive::is_saving::value){
      this->nextException.reset(pe);
   }
}

void MolDSException::Serialize(std::ostream& os){
   boost::archive::text_oarchive oa(os);
   oa << this;
}

MolDSException MolDSException::Deserialize(std::istream& is){
   boost::archive::text_iarchive ia(is);
   MolDSException *p = NULL;
   boost::scoped_ptr<MolDSException> sp(p);
   ia >> p;
   sp.reset(p);

   while(isspace(is.peek())){
      is.get();
   }
   while(!is.eof()){
      try{
         boost::archive::text_iarchive ia(is);
         MolDSException* pnext = NULL;
         ia >> pnext;
         p->LastException()->nextException.reset(pnext);
         p = pnext;
         while(isspace(is.peek())){
            is.get();
         }
      }
      catch(...){
         p->LastException()->nextException.reset();
         break;
      }
   }

   return *sp.get();
}
}
