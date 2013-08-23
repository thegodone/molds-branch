//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   //
// Copyright (C) 2012-2013 Katsuhiko Nishimra                             //
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
#include"../base/Enums.h"
#include"../base/Uncopyable.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/MallocerFreer.h"
#include"../mpi/MpiProcess.h"
#include"../wrappers/Blas.h"
#include"../wrappers/Lapack.h"
#include"../base/EularAngle.h"
#include"../base/Parameters.h"
#include"../base/atoms/Atom.h"
#include"../base/Molecule.h"
#include"../base/ElectronicStructure.h"
#include"Optimizer.h"
#include"BFGS.h"
#include"GEDIIS.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;

namespace MolDS_optimization{
GEDIIS::GEDIIS(){
   this->SetMessages();
   //this->OutputLog("BFGS created\n");
}

GEDIIS::~GEDIIS(){
   //this->OutputLog("BFGS deleted\n");
}

void GEDIIS::SetMessages(){
   BFGS::SetMessages();

   this->errorMessageGeometyrOptimizationNotConverged
      = "Error in optimization::GEDIIS::Optimize: Optimization did not met convergence criterion.\n";
   this->messageStartGEDIISStep
      = "\n==========  START: GEDIIS step ";
}

void GEDIIS::SearchMinimum(boost::shared_ptr<ElectronicStructure> electronicStructure,
                           Molecule& molecule,
                           double* lineSearchedEnergy,
                           bool* obtainesOptimizedStructure) const {
   throw MolDSException("GEDIIS not yet implemented");
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexOptimization();
   double dt = Parameters::GetInstance()->GetTimeWidthOptimization();
   int totalSteps = Parameters::GetInstance()->GetTotalStepsOptimization();
   double maxGradientThreshold = Parameters::GetInstance()->GetMaxGradientOptimization();
   double rmsGradientThreshold = Parameters::GetInstance()->GetRmsGradientOptimization();
   double lineSearchCurrentEnergy   = 0.0;
   double lineSearchInitialEnergy   = 0.0;
   double const* const* matrixForce = NULL;
   double const* vectorForce        = NULL;
   const int dimension = molecule.GetNumberAtoms()*CartesianType_end;
   double** matrixHessian        = NULL;
   double*  vectorOldForce       = NULL;
   double*  vectorStep           = NULL;
   double** matrixStep           = NULL;
   double** matrixOldCoordinates = NULL;
   double*  vectorOldCoordinates = NULL;
   double** matrixDisplacement   = NULL;
   double       trustRadius      = Parameters::GetInstance()->GetInitialTrustRadiusOptimization();
   const double maxNormStep      = Parameters::GetInstance()->GetMaxNormStepOptimization();

   try{
      // initialize Hessian with unit matrix
      MallocerFreer::GetInstance()->Malloc(&matrixHessian, dimension, dimension);
      const double one = 1;
      MolDS_wrappers::Blas::GetInstance()->Dcopy(dimension, &one, 0, &matrixHessian[0][0], dimension+1);

      // initial calculation
      bool requireGuess = true;
      this->UpdateElectronicStructure(electronicStructure, molecule, requireGuess, this->CanOutputLogs());
      lineSearchCurrentEnergy = electronicStructure->GetElectronicEnergy(elecState);

      requireGuess = false;
      matrixForce = electronicStructure->GetForce(elecState);
      vectorForce = &matrixForce[0][0];

      for(int s=0; s<totalSteps; s++){
         this->OutputLog(boost::format("%s%d\n\n") % this->messageStartGEDIISStep % (s+1));

         // Store old Force data
         MallocerFreer::GetInstance()->Malloc(&vectorOldForce, dimension);
         MolDS_wrappers::Blas::GetInstance()->Dcopy(dimension, vectorForce, vectorOldForce);

         this->StoreMolecularGeometry(matrixOldCoordinates, molecule);

         // Level shift Hessian redundant modes
         this->ShiftHessianRedundantMode(matrixHessian, molecule);

         // Limit the trustRadius to maxNormStep
         trustRadius=min(trustRadius,maxNormStep);

         //Calculate RFO step
         MallocerFreer::GetInstance()->Malloc(&matrixStep, molecule.GetNumberAtoms(), CartesianType_end);
         vectorStep = &matrixStep[0][0];
         this->CalcRFOStep(vectorStep, matrixHessian, vectorForce, trustRadius, dimension);

         double approximateChange = this->ApproximateEnergyChange(dimension, matrixHessian, vectorForce, vectorStep);

         // Take a RFO step
         bool doLineSearch = false;
         bool tempCanOutputLogs = false;
         lineSearchInitialEnergy = lineSearchCurrentEnergy;
         if(doLineSearch){
            this->LineSearch(electronicStructure, molecule, lineSearchCurrentEnergy, matrixStep, elecState, dt);
         }
         else{
            this->UpdateMolecularCoordinates(molecule, matrixStep);

            // Broadcast to all processes
            int root = MolDS_mpi::MpiProcess::GetInstance()->GetHeadRank();
            molecule.BroadcastConfigurationToAllProcesses(root);

            this->UpdateElectronicStructure(electronicStructure, molecule, requireGuess, tempCanOutputLogs);
            lineSearchCurrentEnergy = electronicStructure->GetElectronicEnergy(elecState);
         }
         this->OutputMoleculeElectronicStructure(electronicStructure, molecule, this->CanOutputLogs());

         this->UpdateTrustRadius(trustRadius, approximateChange, lineSearchInitialEnergy, lineSearchCurrentEnergy);

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

         if(lineSearchCurrentEnergy > lineSearchInitialEnergy){
            this->OutputLog(this->messageHillClimbing);
            this->RollbackMolecularGeometry(molecule, matrixOldCoordinates);
            lineSearchCurrentEnergy = lineSearchInitialEnergy;
         }

         //Calculate displacement (K_k at Eq. (15) in [SJTO_1983])
         this->CalcDisplacement(matrixDisplacement, matrixOldCoordinates, molecule);

         matrixForce = electronicStructure->GetForce(elecState);
         vectorForce = &matrixForce[0][0];

         // Update Hessian
         this->UpdateHessian(matrixHessian, dimension, vectorForce, vectorOldForce, &matrixDisplacement[0][0]);

      }
      *lineSearchedEnergy = lineSearchCurrentEnergy;
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&matrixHessian, dimension, dimension);
      MallocerFreer::GetInstance()->Free(&vectorOldForce, dimension);
      MallocerFreer::GetInstance()->Free(&matrixStep, molecule.GetNumberAtoms(), CartesianType_end);
      MallocerFreer::GetInstance()->Free(&matrixDisplacement, molecule.GetNumberAtoms(), CartesianType_end);
      MallocerFreer::GetInstance()->Free(&matrixOldCoordinates, molecule.GetNumberAtoms(), CartesianType_end);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free(&matrixHessian, dimension, dimension);
   MallocerFreer::GetInstance()->Free(&vectorOldForce, dimension);
   MallocerFreer::GetInstance()->Free(&matrixStep, molecule.GetNumberAtoms(), CartesianType_end);
   MallocerFreer::GetInstance()->Free(&matrixDisplacement, molecule.GetNumberAtoms(), CartesianType_end);
   MallocerFreer::GetInstance()->Free(&matrixOldCoordinates, molecule.GetNumberAtoms(), CartesianType_end);
}

}
