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
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<math.h>
#include<sstream>
#include<vector>
#include<stdexcept>
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
#include"AtomFactory.h"
#include"../Molecule.h"
#include"../ElectronicStructure.h"
#include"../../cndo/Cndo2.h"
#include"../../indo/Indo.h"
#include"../../zindo/ZindoS.h"
#include"../../mndo/Mndo.h"
#include"../../am1/Am1.h"
#include"../../am1/Am1D.h"
#include"../../pm3/Pm3.h"
#include"../../pm3/Pm3D.h"
#include"../../pm3/Pm3Pddg.h"
#include"ElectronicStructureFactory.h"
using namespace std;
using namespace MolDS_base;
namespace MolDS_base_factories{
string ElectronicStructureFactory::errorMessageNotEnableTheory 
      = "Error in base::ElectronicStructureFactory::Create: Not Enable TheoryType is set.";
string ElectronicStructureFactory::errorMessageTheoryType = "\ttheory type = ";

ElectronicStructure* ElectronicStructureFactory::Create(TheoryType theoryType){
   ElectronicStructure* electronicStructure=NULL;
   if(theoryType == CNDO2 ){
      electronicStructure = new MolDS_cndo::Cndo2();
   }
   else if(theoryType == INDO ){
      electronicStructure = new MolDS_indo::Indo();
   }
   else if(theoryType == ZINDOS ){
      electronicStructure = new MolDS_zindo::ZindoS();
   }
   else if(theoryType == MNDO ){
      electronicStructure = new MolDS_mndo::Mndo();
   }
   else if(theoryType == AM1 ){
      electronicStructure = new MolDS_am1::Am1();
   }
   else if(theoryType == AM1D ){
      electronicStructure = new MolDS_am1::Am1D();
   }
   else if(theoryType == PM3 ){
      electronicStructure = new MolDS_pm3::Pm3();
   }
   else if(theoryType == PM3D ){
      electronicStructure = new MolDS_pm3::Pm3D();
   }
   else if(theoryType == PM3PDDG ){
      electronicStructure = new MolDS_pm3::Pm3Pddg();
   }
   else{
      stringstream ss;
      ss << ElectronicStructureFactory::errorMessageNotEnableTheory << endl;
      ss << ElectronicStructureFactory::errorMessageTheoryType << TheoryTypeStr(theoryType) << endl;
      throw MolDSException(ss.str());
   }
   return electronicStructure;
}

ElectronicStructure* ElectronicStructureFactory::Create(){
   return ElectronicStructureFactory::Create(Parameters::GetInstance()->GetCurrentTheory());
}
}





