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
#include"../wrappers/Lapack.h"
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
   double** matrixHessian          = NULL;
   double** matrixAugmentedHessian = NULL;
   double*  vectorOldForce         = NULL;
   double*  vectorStep             = NULL;
   double** matrixStep             = NULL;
   double*  vectorEigenValues      = NULL;
   double** matrixDisplacement     = NULL;
   double*  P    = NULL; // P_k in eq. 14 on [SJTO_1983]
   double*  K    = NULL; // K_k in eq. 15 on [SJTO_1983]
   double** PP   = NULL; // P_k P_k^T at second term in RHS of Eq. (13) in [SJTO_1983]
   double*  HK   = NULL; // H_k K_k at third term on RHS of Eq. (13) in [SJTO_1983]
   double** HKKH = NULL; // H_k K_k K_k^T H_k at third term on RHS of Eq. (13) in [SJTO_1983]
   const double maxNormStep = 0.1;

   try{
      // initialize Hessian with unit matrix
      MallocerFreer::GetInstance()->Malloc(&matrixHessian, dimension, dimension);
      for(int i=0; i<dimension; i++){
         for(int j=i; j<dimension; j++){
            matrixHessian[i][j] = matrixHessian[j][i] = i==j ? 1.0 : 0.0;
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
         this->OutputLog((boost::format("%s%d\n\n") % this->messageStartBFGSStep % (s+1)).str());

         // Store old Force data
         MallocerFreer::GetInstance()->Malloc(&vectorOldForce, dimension);
         for(int i =0;i < dimension; i++){
            vectorOldForce[i] = vectorForce[i];
         }

         //Store old coordinates for calculating
         MallocerFreer::GetInstance()->Malloc(&matrixDisplacement, molecule.GetNumberAtoms(), CartesianType_end);
         K = &matrixDisplacement[0][0];
         for(int i=0;i<molecule.GetNumberAtoms();i++){
            const Atom*   atom = molecule.GetAtom(i);
            const double* xyz  = atom->GetXyz();
            for(int j=0;j<CartesianType_end;j++){
               matrixDisplacement[i][j] = -xyz[j];
            }
         }

         // Calculate Hessian eigenvalues
         double** tmpmatrix = NULL;
         double*  tmpvector = NULL;
         try{
            MallocerFreer::GetInstance()->Malloc(&tmpmatrix, dimension, dimension);
            for(int i=0;i<dimension;i++){
               for(int j=i;j<dimension;j++){
                  tmpmatrix[i][j] = tmpmatrix[j][i] = matrixHessian[i][j];
               }
            }
            bool calcEigenVectors = false;
            MallocerFreer::GetInstance()->Malloc(&tmpvector, dimension);
            MolDS_wrappers::Lapack::GetInstance()->Dsyevd(&tmpmatrix[0],
                                                          &tmpvector[0],
                                                          dimension,
                                                          calcEigenVectors);
            this->OutputLog("Eigenvalues of the hessian:");
            for(int i=0;i<dimension;i++){
               if((i%6) == 0){
                  this->OutputLog("\n");
                  this->OutputLog((boost::format("%e")%tmpvector[i]).str());
               }
               else{
                  this->OutputLog((boost::format(",\t%e")%tmpvector[i]).str());
               }
            }
            this->OutputLog("\n");
         }
         catch(MolDSException ex){
            MallocerFreer::GetInstance()->Free(&tmpmatrix, dimension, dimension);
            MallocerFreer::GetInstance()->Free(&tmpvector, dimension);
            throw ex;
         }
         MallocerFreer::GetInstance()->Free(&tmpmatrix, dimension, dimension);
         MallocerFreer::GetInstance()->Free(&tmpvector, dimension);

         // Prepare the augmented Hessian
         // See Eq. (4) in [EPW_1997]
         MallocerFreer::GetInstance()->Malloc(&matrixAugmentedHessian, dimension+1,dimension+1);
         for(int i=0;i<dimension;i++){
            for(int j=i;j<dimension;j++){
               // H_k in Eq. (4) in [EPW_1997]
               matrixAugmentedHessian[i][j] = matrixAugmentedHessian[j][i] = matrixHessian[i][j];
            }
         }
         // g_k and g_k^t in Eq. (4) in [EPW_1997]
         for(int i=0;i<dimension;i++){
            // note: gradient = -1 * force
            matrixAugmentedHessian[i][dimension] =
               matrixAugmentedHessian[dimension][i] = - vectorForce[i];
         }
         // 0 in Eq. (4) in [EPW_1997]
         matrixAugmentedHessian[dimension][dimension] = 0;

         // Solve eigenvalue problem on the augmented Hessian
         // See Eq. (4) in [EPW_1997]
         MallocerFreer::GetInstance()->Malloc(&vectorEigenValues, dimension+1);
         //TODO: calculate eigenvalues first then calculate only an eigenvector needed
         bool calcEigenVectors = true;
         MolDS_wrappers::Lapack::GetInstance()->Dsyevd(&matrixAugmentedHessian[0],
                                                       &vectorEigenValues[0],
                                                       dimension+1,
                                                       calcEigenVectors);

         // Select a RFO step as the eigenvector whose eivenvalue is the lowest
         MallocerFreer::GetInstance()->Malloc(&matrixStep, molecule.GetNumberAtoms(), CartesianType_end);
         vectorStep = &matrixStep[0][0];
         for(int i=0;i<dimension;i++){
            // Scale last element of eigenvector to 1 because
            // [vectorStep, 1] is the eigenvector of augmented Hessian.
            // See Eq. (4) in [EPW_1997].
            vectorStep[i] = matrixAugmentedHessian[0][i] / matrixAugmentedHessian[0][dimension];
         }
         //
         // Calculate size of the RFO step
         double normStep = 0;
         for(int i=0;i<dimension;i++){
            normStep += vectorStep[i] * vectorStep[i];
         }
         normStep = sqrt(normStep);

         this->OutputLog((boost::format("Lowest eigenvalue of the augmented Hessian     = %f\n") % vectorEigenValues[0]).str());
         this->OutputLog((boost::format("2nd lowest eigenvalue of the augmented Hessian = %f\n") % vectorEigenValues[1]).str());
         this->OutputLog((boost::format("3rd lowest eigenvalue of the augmented Hessian = %f\n") % vectorEigenValues[2]).str());
         this->OutputLog((boost::format("Calculated RFO step size                       = %f\n") % normStep).str());

         double r; // correctness of the step
         do{

            // Limit the step size to maxNormStep
            if(normStep > maxNormStep){
               for(int i=0;i<dimension;i++){
                  vectorStep[i] *= maxNormStep/normStep;
               }
               normStep = maxNormStep;
               this->OutputLog((boost::format("RFO step size is limited to %f\n") % normStep).str());
            }

            // Calculate approximate change of energy
            double approximateChange = 0;
            for(int i=0;i<dimension;i++){
               approximateChange -= vectorForce[i] * vectorStep[i];
               for(int j=0;j<dimension;j++){
                  approximateChange += vectorStep[i] * matrixHessian[i][j] * vectorStep[j] / 2;
               }
            }

            // Take a RFO step
            bool doLineSearch = false;
            bool tempCanOutputLogs = true;
            lineSearchInitialEnergy = lineSearchCurrentEnergy;
            if(doLineSearch){
               this->LineSearch(electronicStructure, molecule, lineSearchCurrentEnergy, matrixStep, elecState, dt);
            }
            else{
               this->UpdateMolecularCoordinates(molecule, matrixStep);
               this->UpdateElectronicStructure(electronicStructure, molecule, requireGuess, tempCanOutputLogs);
               lineSearchCurrentEnergy = electronicStructure->GetElectronicEnergy(elecState);
            }

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

            // Calculate the correctness of the approximation
            r = (lineSearchCurrentEnergy - lineSearchInitialEnergy)/approximateChange;
            molecule.SetCanOutputLogs(true);
            this->OutputLog((boost::format("actual energy change          = %e\n") % (lineSearchCurrentEnergy-lineSearchInitialEnergy)).str());
            this->OutputLog((boost::format("expected energy change        = %e\n") % approximateChange).str());
            this->OutputLog((boost::format("actual/expected energy change = %f\n") % r).str());
            molecule.SetCanOutputLogs(tempCanOutputLogs);

            // Update the trust radius
            if(r < 0)
            {
               // Rollback molecular geometry
               for(int i=0;i<molecule.GetNumberAtoms();i++){
                  const Atom* atom = molecule.GetAtom(i);
                  double*     xyz  = atom->GetXyz();
                  for(int j=0;j<CartesianType_end;j++){
                     xyz[j] = -matrixDisplacement[i][j];
                  }
               }
               lineSearchCurrentEnergy = lineSearchInitialEnergy;
               // and rerun with smaller trust radius
               maxNormStep /= 4;
            }
            else if(r<0.25){
               maxNormStep /= 4;
            }
            else if(r<0.75){
               // keep trust radius
            }
            else if(r<2){
               maxNormStep *= 2;
            }
            else{
               maxNormStep /= 2;
            }
         }while(r < 0);


         // Update Hessian
         //Calculate displacement (K_k at Eq. (15) in [SJTO_1983])
         for(int i=0;i<molecule.GetNumberAtoms();i++){
            const Atom*   atom = molecule.GetAtom(i);
            const double* xyz  = atom->GetXyz();
            for(int j=0;j<CartesianType_end;j++){
               matrixDisplacement[i][j] += xyz[j];
            }
         }
         matrixForce = electronicStructure->GetForce(elecState);
         vectorForce = &matrixForce[0][0];
         MallocerFreer::GetInstance()->Malloc(&P, dimension);
         for(int i=0; i<dimension; i++){
            // initialize P_k according to Eq. (14) in [SJTO_1983]
            // note: gradient = -1 * force
            P[i] = - (vectorForce[i] - vectorOldForce[i]);
         }
         double PK = 0; // P_k^T K_k at second term at RHS of Eq. (13) in [SJTO_1983]
         for(int i=0; i<dimension;i++){
            PK += P[i] * K[i];
         }
         //P_k P_k^T at second term in RHS of Eq. (13) in [SJTO_1983]
         MallocerFreer::GetInstance()->Malloc(&PP, dimension, dimension);
         for(int i=0; i<dimension;i++){
            for(int j=i;j<dimension;j++){
               PP[i][j] = PP[j][i] = P[i] * P[j];
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
            for(int j=i;j<dimension;j++){
               HKKH[i][j] = 0;
               for(int k=0;k<dimension;k++){
                  HKKH[i][j] += HK[i]*K[k]*matrixHessian[k][j];
               }
               HKKH[j][i] = HKKH[i][j];
            }
         }
         // Calculate H_k+1 according to Eq. (13) in [SJTO_1983]
         for(int i=0;i<dimension;i++){
            for(int j=i;j<dimension;j++){
               matrixHessian[i][j]+= PP[i][j] / PK
                                   - HKKH[i][j] / KHK;
               matrixHessian[j][i] = matrixHessian[i][j];
            }
         }
      }
      *lineSearchedEnergy = lineSearchCurrentEnergy;
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&matrixHessian, dimension, dimension);
      MallocerFreer::GetInstance()->Free(&vectorOldForce, dimension);
      MallocerFreer::GetInstance()->Free(&matrixStep, molecule.GetNumberAtoms(), CartesianType_end);
      MallocerFreer::GetInstance()->Free(&matrixDisplacement, molecule.GetNumberAtoms(), CartesianType_end);
      MallocerFreer::GetInstance()->Free(&matrixAugmentedHessian, dimension+1, dimension+1);
      MallocerFreer::GetInstance()->Free(&vectorEigenValues, dimension+1);
      MallocerFreer::GetInstance()->Free(&P, dimension);
      MallocerFreer::GetInstance()->Free(&PP, dimension, dimension);
      MallocerFreer::GetInstance()->Free(&HK, dimension);
      MallocerFreer::GetInstance()->Free(&HKKH, dimension,dimension);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free(&matrixHessian, dimension, dimension);
   MallocerFreer::GetInstance()->Free(&vectorOldForce, dimension);
   MallocerFreer::GetInstance()->Free(&matrixStep, molecule.GetNumberAtoms(), CartesianType_end);
   MallocerFreer::GetInstance()->Free(&matrixDisplacement, molecule.GetNumberAtoms(), CartesianType_end);
   MallocerFreer::GetInstance()->Free(&matrixAugmentedHessian, dimension+1, dimension+1);
   MallocerFreer::GetInstance()->Free(&vectorEigenValues, dimension+1);
   MallocerFreer::GetInstance()->Free(&P, dimension);
   MallocerFreer::GetInstance()->Free(&PP, dimension, dimension);
   MallocerFreer::GetInstance()->Free(&HK, dimension);
   MallocerFreer::GetInstance()->Free(&HKKH, dimension,dimension);
}
}
