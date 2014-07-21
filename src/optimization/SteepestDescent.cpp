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
#include"../base/constraints/Constraint.h"
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

void SteepestDescent::CalcNextStepGeometry(Molecule &molecule,
                                           OptimizerState& state,
                                           boost::shared_ptr<ElectronicStructure> electronicStructure,
                                           const int elecState,
                                           const double dt) const{
   state.SetInitialEnergy(state.GetCurrentEnergy());

   this->LineSearch(electronicStructure, molecule, state.GetCurrentEnergyRef(), state.GetMatrixForce(), elecState, dt);
}

void SteepestDescent::UpdateState(OptimizerState& state) const{
   this->UpdateSearchDirection(state, state.GetElectronicStructure(), state.GetMolecule(), state.GetConstraint(), state.GetElecState());
}

void SteepestDescent::UpdateSearchDirection(OptimizerState& state,
                                            boost::shared_ptr<ElectronicStructure> electronicStructure,
                                            const MolDS_base::Molecule& molecule,
                                            boost::shared_ptr<MolDS_base_constraints::Constraint> constraint,
                                            int elecState) const{
      // update force
      state.SetMatrixForce(constraint->GetForce(elecState));
}
}
