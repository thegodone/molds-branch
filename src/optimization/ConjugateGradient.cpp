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
#include"../base/MallocerFreer.h"
#include"../base/EularAngle.h"
#include"../base/Parameters.h"
#include"../base/atoms/Atom.h"
#include"../base/Molecule.h"
#include"../base/ElectronicStructure.h"
#include"Optimizer.h"
#include"ConjugateGradient.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_optimization{
ConjugateGradient::ConjugateGradient(){
   this->SetMessages();
   //this->OutputLog("ConjugateGradient created\n");
}

ConjugateGradient::~ConjugateGradient(){
   //this->OutputLog("ConjugateGradient deleted\n");
}

void ConjugateGradient::SetMessages(){
   Optimizer::SetMessages();
   this->errorMessageNotEnebleTheoryType  
      = "Error in optimization::ConjugateGradient::CheckEnableTheoryType: Non available theory is set.\n";
   this->errorMessageGeometyrOptimizationNotConverged 
      = "Error in optimization::ConjugateGradient::Optimize: Optimization did not met convergence criterion.\n";
   this->messageStartConjugateGradientStep = "\n==========  START: Conjugate gradient step ";
}

void ConjugateGradient::SearchMinimum(boost::shared_ptr<ElectronicStructure> electronicStructure,
                                      Molecule& molecule,
                                      double* lineSearchedEnergy,
                                      bool* obtainesOptimizedStructure) const{
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexOptimization();
   double dt = Parameters::GetInstance()->GetTimeWidthOptimization();
   int totalSteps = Parameters::GetInstance()->GetTotalStepsOptimization();
   double maxGradientThreshold = Parameters::GetInstance()->GetMaxGradientOptimization();
   double rmsGradientThreshold = Parameters::GetInstance()->GetRmsGradientOptimization();
   double lineSearchCurrentEnergy = 0.0;
   double lineSearchInitialEnergy = 0.0;
   double** matrixForce = NULL;
   double** oldMatrixForce = NULL;
   double** matrixSearchDirection = NULL;

   // initial calculation
   bool requireGuess = true;
   this->UpdateElectronicStructure(electronicStructure, molecule, requireGuess, this->CanOutputLogs());
   lineSearchCurrentEnergy = electronicStructure->GetElectronicEnergy(elecState);

   requireGuess = false;
   matrixForce = electronicStructure->GetForce(elecState);
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&oldMatrixForce, molecule.GetNumberAtoms(), CartesianType_end);
      MallocerFreer::GetInstance()->Malloc<double>(&matrixSearchDirection, molecule.GetNumberAtoms(), CartesianType_end);
      for(int a=0;a<molecule.GetNumberAtoms();a++){
         for(int i=0; i<CartesianType_end; i++){
            matrixSearchDirection[a][i] = matrixForce[a][i];
         }
      }
      
      // conugate gradient roop
      for(int s=0; s<totalSteps; s++){
         this->OutputLog((boost::format("%s%d\n\n") % this->messageStartConjugateGradientStep.c_str() % (s+1)).str());
         lineSearchInitialEnergy = lineSearchCurrentEnergy;

         // line search roop
         bool tempCanOutputLogs = false;
         int lineSearchSteps = 0;
         double lineSearchOldEnergy = lineSearchCurrentEnergy;
         while(lineSearchCurrentEnergy <= lineSearchOldEnergy){
            this->UpdateMolecularCoordinates(molecule, matrixSearchDirection, dt);
            this->UpdateElectronicStructure(electronicStructure, molecule, requireGuess, tempCanOutputLogs);
            lineSearchOldEnergy = lineSearchCurrentEnergy;
            lineSearchCurrentEnergy = electronicStructure->GetElectronicEnergy(elecState);
            lineSearchSteps++;
         }

         // final state of line search
         this->OutputLog((boost::format("%s%d\n\n") % this->messageLineSearchSteps.c_str() % lineSearchSteps).str());
         this->UpdateMolecularCoordinates(molecule, matrixSearchDirection, -0.5*dt);
         this->UpdateElectronicStructure(electronicStructure, molecule, requireGuess, tempCanOutputLogs);
         this->OutputMoleculeElectronicStructure(electronicStructure, molecule, this->CanOutputLogs());

         // update matrixSearchDirection
         this->UpdateSearchDirection(&matrixForce, oldMatrixForce, matrixSearchDirection, electronicStructure, molecule, elecState);

         lineSearchCurrentEnergy = electronicStructure->GetElectronicEnergy(elecState);

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
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free<double>(&oldMatrixForce, molecule.GetNumberAtoms(), CartesianType_end);
      MallocerFreer::GetInstance()->Free<double>(&matrixSearchDirection, molecule.GetNumberAtoms(), CartesianType_end);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free<double>(&oldMatrixForce, molecule.GetNumberAtoms(), CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(&matrixSearchDirection, molecule.GetNumberAtoms(), CartesianType_end);
   *lineSearchedEnergy = lineSearchCurrentEnergy;
}

void ConjugateGradient::UpdateSearchDirection(double*** matrixForce, 
                                              double** oldMatrixForce, 
                                              double** matrixSearchDirection,
                                              boost::shared_ptr<ElectronicStructure> electronicStructure, 
                                              const MolDS_base::Molecule& molecule,
                                              int elecState) const{
   for(int a=0;a<molecule.GetNumberAtoms();a++){
      for(int i=0; i<CartesianType_end; i++){
         oldMatrixForce[a][i] = (*matrixForce)[a][i];
      }
   }
   *matrixForce = electronicStructure->GetForce(elecState);
   double beta=0.0;
   double temp=0.0;
   for(int a=0;a<molecule.GetNumberAtoms();a++){
      for(int i=0; i<CartesianType_end; i++){
         temp += pow(oldMatrixForce[a][i],2.0);
         beta += ((*matrixForce)[a][i] - oldMatrixForce[a][i])*(*matrixForce)[a][i];
      }
   }
   beta /= temp;
   for(int a=0;a<molecule.GetNumberAtoms();a++){
      for(int i=0; i<CartesianType_end; i++){
         matrixSearchDirection[a][i] *= beta;
         matrixSearchDirection[a][i] += (*matrixForce)[a][i];
      }
   }
}


}



