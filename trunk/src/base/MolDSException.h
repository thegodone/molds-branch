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
#ifndef INCLUDED_MOLDSEXCEPTION
#define INCLUDED_MOLDSEXCEPTION
#include<map>
#include<boost/shared_array.hpp>
namespace MolDS_base{
class MolDSException : public std::domain_error {
public:
   explicit MolDSException(std::string cause);
#ifdef BOOST_FORMAT_HPP
   MolDSException(const boost::format& cause);
#endif
   ~MolDSException() throw(){};
   void PrintBacktrace();
   template <class T>
   T GetKeyValue(int key);
   template <class T>
   void SetKeyValue(int key, T value);
   bool HasKey(int key);
   virtual const char* what() const throw();
private:
   void GetBacktrace(int bufsize);
   size_t backtraceSize;
   boost::shared_array<void*> backtracePtr;
   typedef std::map<int, int> intKeyValueMap_t;
   intKeyValueMap_t intKeyValueMap;
   //typedef std::map<int, other> otherKeyValueMap_t;
   //otherKeyValueMap_t otherKeyValueMap;
};
}
#endif





