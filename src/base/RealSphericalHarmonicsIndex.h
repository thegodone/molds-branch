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
#ifndef INCLUDED_REAL_SHPERICAL_HARMONICS_INDEX
#define INCLUDED_REAL_SHPERICAL_HARMONICS_INDEX
namespace MolDS_base{

// l and m correspond to l and m in Eq.(15) in J. Comp. Chem. 20, 383(1999)
class RealSphericalHarmonicsIndex {
public:
   RealSphericalHarmonicsIndex(int l, int m);
   explicit RealSphericalHarmonicsIndex(MolDS_base::OrbitalType  orbitalType);
   int GetL() const{return this->l;}
   int GetM() const{return this->m;}
private:
   int l;
   int m;
   RealSphericalHarmonicsIndex();
};

}
#endif
