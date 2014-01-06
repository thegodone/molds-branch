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
#ifndef INCLUDED_ATOMFACTORY
#define INCLUDED_ATOMFACTORY
namespace MolDS_base_factories{

// AtomFactory is singleton
class AtomFactory{
public:
   static MolDS_base_atoms::Atom* Create(MolDS_base::AtomType atomType,
                                         int index,
                                         double x,
                                         double y,
                                         double z,
                                         double px,
                                         double py,
                                         double pz);
   static MolDS_base_atoms::Atom* Create(MolDS_base::AtomType atomType,
                                         int index,
                                         double x,
                                         double y,
                                         double z,
                                         double px,
                                         double py,
                                         double pz,
                                         double charge);
   static MolDS_base_atoms::Atom* Create(MolDS_base::AtomType atomType,
                                         int index,
                                         double x,
                                         double y,
                                         double z);
   static MolDS_base_atoms::Atom* Create(MolDS_base::AtomType atomType,
                                         int index,
                                         double x,
                                         double y,
                                         double z,
                                         double charge);
private:
   AtomFactory(); 
   ~AtomFactory(); 
   static std::string errorMessageNotEnableAtom;
   static std::string errorMessageNotEnvironmentalPointCharge;
   static std::string errorMessageAtomType;
};

}
#endif





