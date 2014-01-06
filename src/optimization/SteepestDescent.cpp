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
#include<math.h>
#include<string>
#include<vector>
#include<stdexcept>
#include<boost/shared_ptr.hpp>
#include<boost/format.hpp>
#include"../base/Enums.h"
#include"../base/Uncopyable.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
#include"../base/EularAngle.h"
#include"../base/Parameters.h"
#include"../base/RealSphericalHarmonicsIndex.h"
#include"../base/atoms/Atom.h"
#include"../base/Molecule.h"
#include"../base/ElectronicStructure.h"
#include"Optimizer.h"
#include"SteepestDescent.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_optimization{
SteepestDescent::SteepestDescent(){
   this->SetMessages();
   //this->OutputLog("SteepestDescent created\n");
}

SteepestDescent::~SteepestDescent(){
   //this->OutputLog("SteepestDescent deleted\n");
}

void SteepestDescent::SetMessages(){
   Optimizer::SetMessages();
   this->errorMessageNotEnebleTheoryType  
      = "Error in optimization::SteepestDescent::CheckEnableTheoryType: Non available theory is set.\n";
   this->errorMessageGeometyrOptimizationNotConverged 
      = "Error in optimization::SteepestDescent::Optimize: Optimization did not met convergence criterion.\n";
   this->messageStartSteepestDescentStep = "\n==========  START: Steepest Descent step ";
}

void SteepestDescent::SearchMinimum(boost::shared_ptr<ElectronicStructure> electronicStructure,
                                    Molecule& molecule,
                                    double* lineSearchedEnergy,
                                    bool* obtainesOptimizedStructure) const{
   int    elecState            = Parameters::GetInstance()->GetElectronicStateIndexOptimization();
   double dt                   = Parameters::GetInstance()->GetTimeWidthOptimization();
   int    totalSteps           = Parameters::GetInstance()->GetTotalStepsOptimization();
   double maxGradientThreshold = Parameters::GetInstance()->GetMaxGradientOptimization();
   double rmsGradientThreshold = Parameters::GetInstance()->GetRmsGradientOptimization();
   double lineSearchCurrentEnergy   = 0.0;
   double lineSearchInitialEnergy   = 0.0;
   double const* const* matrixForce = NULL;

   // initial calculation
   bool requireGuess = true;
   this->UpdateElectronicStructure(electronicStructure, molecule, requireGuess, this->CanOutputLogs());

   requireGuess = false;
   matrixForce = electronicStructure->GetForce(elecState);
   lineSearchCurrentEnergy = electronicStructure->GetElectronicEnergy(elecState);
   for(int s=0; s<totalSteps; s++){
      this->OutputLog(boost::format("%s%d\n\n") % this->messageStartSteepestDescentStep.c_str() % (s+1));
      lineSearchInitialEnergy = lineSearchCurrentEnergy;

      // do line search
      this->LineSearch(electronicStructure, molecule, lineSearchCurrentEnergy, matrixForce, elecState, dt);

      // update force
      matrixForce = electronicStructure->GetForce(elecState);

      // check convergence
      if(this->SatisfiesConvergenceCriterion(matrixForce, 
                                             molecule,
                                             lineSearchInitialEnergy, 
                                             lineSearchCurrentEnergy,
                                             maxGradientThreshold, 
                                             rmsGradientThreshold)){
         *obtainesOptimizedStructure = true;
         break;
      }

   }
   *lineSearchedEnergy = lineSearchCurrentEnergy;
}
}



