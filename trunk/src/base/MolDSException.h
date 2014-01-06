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
#ifndef INCLUDED_MOLDSEXCEPTION
#define INCLUDED_MOLDSEXCEPTION
#include<map>
#include<boost/shared_ptr.hpp>
#include<boost/shared_array.hpp>
#include<boost/serialization/access.hpp>
namespace MolDS_base{
class MolDSException : public std::domain_error {
public:
   explicit MolDSException(std::string cause);
#ifdef BOOST_FORMAT_HPP
   MolDSException(const boost::format& cause);
#endif
   ~MolDSException() throw(){};
   template <class T>
   T GetKeyValue(int key);
   template <class T>
   void SetKeyValue(int key, T value);
   bool HasKey(int key);
   const MolDSException* NextException() const{return this->nextException.get();}
   virtual const char* what() const throw();
   const char* What() const throw(){return domain_error::what();}
   void Serialize(std::ostream& os);
   static MolDSException Deserialize(std::istream& is);
private:
   void GetBacktrace(int bufsize);
   size_t backtraceSize;
   boost::shared_array<void*> backtracePtr;

   typedef std::map<int, int> intKeyValueMap_t;
   intKeyValueMap_t intKeyValueMap;
   //typedef std::map<int, other> otherKeyValueMap_t;
   //otherKeyValueMap_t otherKeyValueMap;

   boost::shared_ptr<MolDSException> nextException;
   MolDSException* LastException(){
      if(this->nextException.get()==NULL){
         return this;
      }
      else{
         return this->nextException->LastException();
      }
   }

   friend class boost::serialization::access;
   template<class Archive>
   void serialize(Archive& ar, const unsigned int ver);
};
}
#endif





