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
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<vector>
#include<stdexcept>
#include<boost/shared_ptr.hpp>
#include<boost/format.hpp>
#include"../Enums.h"
#include"../Uncopyable.h"
#include"../PrintController.h"
#include"../MolDSException.h"
#include"../MallocerFreer.h"
#include"../../mpi/MpiInt.h"
#include"../../mpi/MpiProcess.h"
#include"../EularAngle.h"
#include"../Parameters.h"
#include"../RealSphericalHarmonicsIndex.h"
#include"../atoms/Atom.h"
#include"../Molecule.h"
#include"../ElectronicStructure.h"
#include"../constraints/Constraint.h"
#include"../constraints/SpaceFixedAtoms.h"
#include"../constraints/NonConstraint.h"
#include"ConstraintFactory.h"
using namespace std;
using namespace MolDS_base;
namespace MolDS_base_factories{

MolDS_base_constraints::Constraint* ConstraintFactory::Create(const Molecule& molecule,
                                                             boost::shared_ptr<ElectronicStructure> electronicStructure){
   MolDS_base_constraints::Constraint* c=NULL;
   if(Parameters::GetInstance()->RequiresSpaceFixedAtomsOptimization()){
      c = new MolDS_base_constraints::SpaceFixedAtoms(&molecule, electronicStructure);
   }
   else{
      c = new MolDS_base_constraints::NonConstraint(&molecule, electronicStructure);
   }
   c->SetConstraintCondition();
   return c;
}

}





