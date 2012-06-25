//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   //
// Copyright (C) 2012-2012 Katsuhiko Nishimra                             //
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
#include<string.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<vector>
#include<stdexcept>
#include<boost/shared_ptr.hpp>
#include<boost/format.hpp>
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/Uncopyable.h"
#include"../base/Enums.h"
#include"../base/EularAngle.h"
#include"../base/Parameters.h"
#include"../base/atoms/Atom.h"
#include"../base/Molecule.h"
#include"../base/ElectronicStructure.h"
#include"../base/MallocerFreer.h"
#include"Optimizer.h"
#include"BFGS.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_optimization{
BFGS::BFGS(){
   this->SetMessages();
   //this->OutputLog("BFGS created\n");
}

BFGS::~BFGS(){
   //this->OutputLog("BFGS deleted\n");
}

void BFGS::SetMessages(){
   Optimizer::SetMessages();
   this->errorMessageNotEnebleTheoryType
      = "Error in optimization::BFGS::CheckEnableTheoryType: Non available theory is set.\n";
   this->errorMessageGeometyrOptimizationNotConverged
      = "Error in optimization::BFGS::Optimize: Optimization did not met convergence criterion.\n";
   this->messageStartBFGSStep = "\n==========  START: BFGS step ";
}

void BFGS::SearchMinimum(boost::shared_ptr<ElectronicStructure> electronicStructure,
                         Molecule& molecule,
                         double* lineSearchedEnergy,
                         bool* obtainesOptimizedStructure) const {
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexOptimization();
   double dt = Parameters::GetInstance()->GetTimeWidthOptimization();
   int totalSteps = Parameters::GetInstance()->GetTotalStepsOptimization();
   double maxGradientThreshold = Parameters::GetInstance()->GetMaxGradientOptimization();
   double rmsGradientThreshold = Parameters::GetInstance()->GetRmsGradientOptimization();
   double lineSearchCurrentEnergy = 0.0;
   double lineSearchInitialEnergy = 0.0;
   double** matrixForce = NULL;
   double* vectorForce = NULL;
   const int dimension = molecule.GetNumberAtoms()*CartesianType_end;
   double matrixHessianInverse[dimension][dimension];

   // initialize Hessian with unit matrix
   for(int i=0; i<dimension; i++){
      for(int j=0; j<dimension; j++){
         matrixHessianInverse[i][j] = i==j ? 1.0 : 0.0;
      }
   }

   // initial calculation
   bool requireGuess = true;
   this->UpdateElectronicStructure(electronicStructure, molecule, requireGuess, this->CanOutputLogs());
   lineSearchCurrentEnergy = electronicStructure->GetElectronicEnergy(elecState);

   requireGuess = false;
   matrixForce = electronicStructure->GetForce(elecState);
   vectorForce = &matrixForce[0][0];

   for(int s=0; s<totalSteps; s++){
      // Store old Force data
      double vectorOldForce[dimension];
      for(int i =0;i < dimension; i++){
         vectorOldForce[i] = vectorForce[i];
      }

      //Store old coordinates for calculating
      double vectorDisplacement[dimension];
      for(int i=0;i<molecule.GetNumberAtoms();i++){
         const Atom* atom = molecule.GetAtom(i);
         for(int j=0;j<CartesianType_end;j++){
            vectorDisplacement[i*CartesianType_end+j] = - atom->GetXyz()[j];
         }
      }

      // Calculate new search direction
      double vectorDirection[dimension];
      double *matrixDirection[molecule.GetNumberAtoms()];
      for(int i=0;i<molecule.GetNumberAtoms();i++){
         matrixDirection[i]=&vectorDirection[i*CartesianType_end];
      }
      for(int i=0;i<dimension;i++){
         vectorDirection[i] = 0;
         for(int j=0;j<dimension;j++){
            vectorDirection[i] += matrixHessianInverse[i][j]*vectorForce[j];
         }
      }

      // Store initial energy
      lineSearchInitialEnergy = lineSearchCurrentEnergy;

      // do line search
      this->LineSearch(electronicStructure, molecule, lineSearchCurrentEnergy, matrixDirection, elecState, dt);
      matrixForce = electronicStructure->GetForce(elecState);
      vectorForce = &matrixForce[0][0];

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

      //Calculate displacement
      for(int i=0;i<molecule.GetNumberAtoms();i++){
         const Atom* atom = molecule.GetAtom(i);
         for(int j=0;j<CartesianType_end;j++){
            vectorDisplacement[i*CartesianType_end+j] += atom->GetXyz()[j];
         }
      }
      // Update Hessian
      double matrixOldHessianInverse[dimension][dimension];
      for(int i=0;i<dimension;i++){
         for(int j=0;j<dimension;j++){
            matrixOldHessianInverse[i][j]=matrixHessianInverse[i][j];
         }
      }
      double sy = 0;
      double yby = 0;
      for(int i=0; i<dimension; i++){
         sy += vectorDisplacement[i] * (vectorForce[i] - vectorOldForce[i]);
         double byi = 0;
         for(int j=0; j<dimension; j++){
            byi += matrixHessianInverse[i][j] * (vectorForce[j] - vectorOldForce[j]);
         }
         yby += (vectorForce[i] - vectorOldForce[i]) * byi;
      }
      double syybysysy = (sy+yby)/sy/sy;
      for(int i=0;i<dimension;i++){
         for(int j=0;j<dimension;j++){
            matrixHessianInverse[i][j] = matrixOldHessianInverse[i][j] + syybysysy *vectorDisplacement[i] * vectorDisplacement[j];
            double byssyb = 0;
            for(int k=0;k<dimension;k++){
               byssyb += matrixOldHessianInverse[i][k] * (vectorForce[k] - vectorOldForce[k]) * vectorDisplacement[j];
               byssyb += vectorDisplacement[i] * (vectorForce[k] - vectorOldForce[k]) * matrixOldHessianInverse[k][j];
            }
            matrixHessianInverse[i][j] -= byssyb / sy;
         }
      }
   }
   *lineSearchedEnergy = lineSearchCurrentEnergy;
}
}