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
#ifndef INCLUDED_PM3
#define INCLUDED_PM3
namespace MolDS_pm3{

/***
 *  Main References for PM3 are [S_1989, S_1989-2, S_1991, S_2004, S_2007]
 */
class Pm3 : public MolDS_am1::Am1{
public:
   Pm3();
   virtual ~Pm3();
protected:
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
private:
};

}
#endif



