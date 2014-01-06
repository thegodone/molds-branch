//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   // 
// Copyright (C) 2012-2014 Katsuhiko Nishimra                             // 
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
   PrintController();
   virtual ~PrintController();
   bool CanOutputLogs() const               {return this->canOutputLogs;}
   void SetCanOutputLogs(bool canOutputLogs){this->canOutputLogs = canOutputLogs;}
protected:
   void OutputLog(const boost::format& log) const{this->OutputLog(log.str());}
   void OutputLog(std::string log) const;
private:
   bool canOutputLogs;
};
}
#endif
