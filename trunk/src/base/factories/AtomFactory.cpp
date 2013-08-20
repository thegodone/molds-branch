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
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<vector>
#include<stdexcept>
#include<boost/format.hpp>
#include"../Enums.h"
#include"../Uncopyable.h"
#include"../PrintController.h"
#include"../MolDSException.h"
#include"../MallocerFreer.h"
#include"../../mpi/MpiProcess.h"
#include"../RealSphericalHarmonicsIndex.h"
#include"../atoms/Atom.h"
#include"../atoms/Hatom.h"
#include"../atoms/Liatom.h"
#include"../atoms/Catom.h"
#include"../atoms/Natom.h"
#include"../atoms/Oatom.h"
#include"../atoms/Satom.h"
#include"AtomFactory.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_base_factories{

string AtomFactory::errorMessageNotEnableAtom = "Error in base::AtomFactory::Create: Not Enable AtomType is set.";
string AtomFactory::errorMessageAtomType = "\tatom type = ";

Atom* AtomFactory::Create(AtomType atomType, int index, double x, double y, double z, double px, double py, double pz){
   Atom* atom=NULL;
   if(atomType == H){
      atom = new Hatom(index);
   }
   else if(atomType == Li){
      atom = new Liatom(index);
   }
   else if(atomType == C){
      atom = new Catom(index);
   }
   else if(atomType == N){
      atom = new Natom(index);
   }
   else if(atomType == O){
      atom = new Oatom(index);
   }
   else if(atomType == S){
      atom = new Satom(index);
   }
   else{
      stringstream ss;
      ss << AtomFactory::errorMessageNotEnableAtom << endl;
      ss << AtomFactory::errorMessageAtomType << AtomTypeStr(atomType) << endl;
      throw MolDSException(ss.str());
   }
   atom->SetXyz(x, y, z);
   atom->SetPxyz(px, py, pz);
   return atom;
}

Atom* AtomFactory::Create(AtomType atomType, int index, double x, double y, double z){
   double px=0.0;
   double py=0.0;
   double pz=0.0;
   return AtomFactory::Create(atomType, index, x, y, z, px, py, pz);
}
}





