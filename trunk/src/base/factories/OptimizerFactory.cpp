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
#include"../../optimization/Optimizer.h"
#include"../../optimization/ConjugateGradient.h"
#include"../../optimization/BFGS.h"
#include"../../optimization/GEDIIS.h"
#include"../../optimization/SteepestDescent.h"
#include"OptimizerFactory.h"
using namespace std;
using namespace MolDS_base;
namespace MolDS_base_factories{
string OptimizerFactory::errorMessageNotEnableOptimizationMethod
      = "Error in base_factories::OptimizerFactory::Create: Not Enable OptimizationMethodType is set.";
string OptimizerFactory::errorMessageOptimizationMethodType = "\tMethod = ";

MolDS_optimization::Optimizer* OptimizerFactory::Create(OptimizationMethodType methodType){
   MolDS_optimization::Optimizer* optimizer=NULL;
   if(methodType == ConjugateGradientMethod ){
      optimizer = new MolDS_optimization::ConjugateGradient();
   }
   else if(methodType == BFGSMethod ){
      optimizer = new MolDS_optimization::BFGS();
   }
   else if(methodType == GEDIISMethod ){
      optimizer = new MolDS_optimization::GEDIIS();
   }
   else if(methodType == SteepestDescentMethod ){
      optimizer = new MolDS_optimization::SteepestDescent();
   }
   else{
      stringstream ss;
      ss << OptimizerFactory::errorMessageNotEnableOptimizationMethod << endl;
      ss << OptimizerFactory::errorMessageOptimizationMethodType << OptimizationMethodTypeStr(methodType) << endl;
      throw MolDSException(ss.str());
   }
   return optimizer;
}

MolDS_optimization::Optimizer* OptimizerFactory::Create(){
   return OptimizerFactory::Create(Parameters::GetInstance()->GetMethodOptimization());
}
}





