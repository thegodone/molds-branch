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
#include"../base/MallocerFreer.h"
#include"../base/EularAngle.h"
#include"../base/Parameters.h"
#include"../base/atoms/Atom.h"
#include"../base/Molecule.h"
#include"../base/ElectronicStructure.h"
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
   double** matrixHessian = NULL;
   double* vectorOldForce = NULL;
   double* vectorDirection = NULL;
   double** matrixDirection = NULL;
   double*  P    = NULL; // P_k in eq. 14 on [SJTO_1983]
   double*  K    = NULL; // K_k in eq. 15 on [SJTO_1983]
   double** PP   = NULL; // P_k P_k^T at second term in RHS of Eq. (13) in [SJTO_1983]
   double*  HK   = NULL; // H_k K_k at third term on RHS of Eq. (13) in [SJTO_1983]
   double** HKKH = NULL; // H_k K_k K_k^T H_k at third term on RHS of Eq. (13) in [SJTO_1983]

   try{
      // initialize Hessian with unit matrix
      MallocerFreer::GetInstance()->Malloc(&matrixHessian, dimension, dimension);
      for(int i=0; i<dimension; i++){
         for(int j=0; j<dimension; j++){
            matrixHessian[i][j] = i==j ? 1.0 : 0.0;
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
         MallocerFreer::GetInstance()->Malloc(&vectorOldForce, dimension);
         for(int i =0;i < dimension; i++){
            vectorOldForce[i] = vectorForce[i];
         }

         //Store old coordinates for calculating
         MallocerFreer::GetInstance()->Malloc(&K, dimension);
         for(int i=0;i<molecule.GetNumberAtoms();i++){
            const Atom* atom = molecule.GetAtom(i);
            for(int j=0;j<CartesianType_end;j++){
               K[i*CartesianType_end+j] = - atom->GetXyz()[j];
            }
         }

         // Calculate new search direction
         MallocerFreer::GetInstance()->Malloc(&matrixDirection, molecule.GetNumberAtoms(), CartesianType_end);
         vectorDirection = &matrixDirection[0][0];
         for(int i=0;i<dimension;i++){
            vectorDirection[i] = 0;
            for(int j=0;j<dimension;j++){
               vectorDirection[i] += matrixHessian[i][j]*vectorForce[j];
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

         //Calculate displacement (K_k at Eq. (15) in [SJTO_1983])
         for(int i=0;i<molecule.GetNumberAtoms();i++){
            const Atom* atom = molecule.GetAtom(i);
            for(int j=0;j<CartesianType_end;j++){
               K[i*CartesianType_end+j] += atom->GetXyz()[j];
            }
         }
         // Update Hessian
         MallocerFreer::GetInstance()->Malloc(&P, dimension);
         for(int i=0; i<dimension; i++){
            // initialize P_k according to Eq. (14) in [SJTO_1983]
            P[i] = vectorForce[i] - vectorOldForce[i];
         }
         double PK = 0; // P_k^T K_k at second term at RHS of Eq. (13) in [SJTO_1983]
         for(int i=0; i<dimension;i++){
            PK += P[i] * K[i];
         }
         //P_k P_k^T at second term in RHS of Eq. (13) in [SJTO_1983]
         MallocerFreer::GetInstance()->Malloc(&PP, dimension, dimension);
         for(int i=0; i<dimension;i++){
            for(int j=0;j<dimension;j++){
               PP[i][j] = P[i] * P[j];
            }
         }
         //H_k K_k at third term on RHS of Eq. (13) in [SJTO_1983]
         MallocerFreer::GetInstance()->Malloc(&HK, dimension);
         for(int i=0; i<dimension; i++){
            HK[i] = 0;
            for(int j=0; j<dimension; j++){
               HK[i] += matrixHessian[i][j] * K[j];
            }
         }
         double KHK = 0; //K_k^T H_k K_k at third term on RHS of Eq. (13) in [SJTO_1983]
         for(int i=0;i<dimension;i++){
            KHK += K[i]*HK[i];
         }
         //H_k K_k K_k^T H_k at third term on RHS of Eq. (13) in [SJTO_1983]
         MallocerFreer::GetInstance()->Malloc(&HKKH, dimension,dimension);
         for(int i=0;i<dimension;i++){
            for(int j=0;j<dimension;j++){
               HKKH[i][j] = 0;
               for(int k=0;k<dimension;k++){
                  HKKH[i][j] += HK[i]*K[k]*matrixHessian[k][j];
               }
            }
         }
         // Calculate H_k+1 according to Eq. (13) in [SJTO_1983]
         for(int i=0;i<dimension;i++){
            for(int j=0;j<dimension;j++){
               matrixHessian[i][j]+= PP[i][j] / PK
                                   - HKKH[i][j] / KHK;
            }
         }

      }
      *lineSearchedEnergy = lineSearchCurrentEnergy;
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&matrixHessian, dimension, dimension);
      MallocerFreer::GetInstance()->Free(&vectorOldForce, dimension);
      MallocerFreer::GetInstance()->Free(&K, dimension);
      MallocerFreer::GetInstance()->Free(&matrixDirection, molecule.GetNumberAtoms(), CartesianType_end);
      MallocerFreer::GetInstance()->Free(&P, dimension);
      MallocerFreer::GetInstance()->Free(&PP, dimension, dimension);
      MallocerFreer::GetInstance()->Free(&HK, dimension);
      MallocerFreer::GetInstance()->Free(&HKKH, dimension,dimension);
   }
   MallocerFreer::GetInstance()->Free(&matrixHessian, dimension, dimension);
   MallocerFreer::GetInstance()->Free(&vectorOldForce, dimension);
   MallocerFreer::GetInstance()->Free(&K, dimension);
   MallocerFreer::GetInstance()->Free(&matrixDirection, molecule.GetNumberAtoms(), CartesianType_end);
   MallocerFreer::GetInstance()->Free(&P, dimension);
   MallocerFreer::GetInstance()->Free(&PP, dimension, dimension);
   MallocerFreer::GetInstance()->Free(&HK, dimension);
   MallocerFreer::GetInstance()->Free(&HKKH, dimension,dimension);
}
}
