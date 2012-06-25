//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
// Copyright (C) 2012-2012 Katsuhiko Nishimra                             // 
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
#ifndef INCLUDED_PRINTCONTROLLER
#define INCLUDED_PRINTCONTROLLER
namespace MolDS_base{

class PrintController{
public:
   PrintController(){
      this->canOutputLogs = true;
      //this->OutputLog("printController is created.\n");
   }
   virtual ~PrintController(){
      //this->OutputLog("printController is destructed.\n";
   }
   bool CanOutputLogs() const{
      return this->canOutputLogs;
   }
   void SetCanOutputLogs(bool canOutputLogs){
      this->canOutputLogs = canOutputLogs;
   }
protected:
   void OutputLog(std::string log) const{
      if(this->canOutputLogs){
         bool endl = false;
         std::string::reverse_iterator iter;
         for(iter = log.rbegin(); iter != log.rend(); iter++){
            if(*iter == '\n'){
               std::string::iterator fwditer = iter.base();
               log.erase(--fwditer);
               endl = true;
               break;
            }
            else if(*iter != '\0'){
               break;
            }
         }
         std::cout << log;
         if(endl){std::cout << std::endl;}
      }
   }
private:
   bool canOutputLogs;
};
}
#endif