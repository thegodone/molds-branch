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
   double** matrixHessian        = NULL;
   double*  vectorOldForce       = NULL;
   double*  vectorStep           = NULL;
   double** matrixStep           = NULL;
   double** matrixOldCoordinates = NULL;
   double*  vectorOldCoordinates = NULL;
   double** matrixDisplacement   = NULL;
   double maxNormStep = 0.3;
   const double largeEigenvalue = 1.0e3;

   try{
      // initialize Hessian
      MallocerFreer::GetInstance()->Malloc(&matrixHessian, dimension, dimension);
      for(int i=0; i<dimension; i++){
         for(int j=i; j<dimension; j++){
            // Make translational mode eigenvalue to be largeEigenvalue
            matrixHessian[i][j] = matrixHessian[j][i] = (i%3==j%3 ? (largeEigenvalue-1)/molecule.GetNumberAtoms() : 0.0)
                                                      + (i==j ? 1:0);
         }
      }
      // Make rotational mode eigenvalue to be large eigenvalue
      MallocerFreer::GetInstance()->Malloc(&matrixOldCoordinates, molecule.GetNumberAtoms(), CartesianType_end);
      vectorOldCoordinates = &matrixOldCoordinates[0][0];
      for(int c=0;c<CartesianType_end;c++){
         MallocerFreer::GetInstance()->Initialize<double>(matrixOldCoordinates,
                                                          molecule.GetNumberAtoms(),
                                                          CartesianType_end);
         // Extract the rotational mode orthogonal to the translational mode
         double dotProductTransRot = 0; // Dot product of rotational and translational mode
         for(int n=0;n<molecule.GetNumberAtoms();n++){
            const Atom*   atom = molecule.GetAtom(n);
            const double* xyz  = atom->GetXyz();
            matrixOldCoordinates[n][c] = xyz[c];
            dotProductTransRot += xyz[c] / sqrt(molecule.GetNumberAtoms()*1.0);
         }
         double rotNorm2 = 0; // Norm of the extracted rotational mode
         for(int n=0;n<molecule.GetNumberAtoms();n++){
            // Orthogonalization
            matrixOldCoordinates[n][c] -= dotProductTransRot / sqrt(molecule.GetNumberAtoms());
            rotNorm2 += matrixOldCoordinates[n][c] * matrixOldCoordinates[n][c];
         }
         for(int i=0; i<dimension; i++){
            for(int j=i; j<dimension; j++){
               matrixHessian[i][j] += vectorOldCoordinates[i]
                                    * vectorOldCoordinates[j]
                                    / rotNorm2 * (largeEigenvalue-1);
               matrixHessian[j][i] = matrixHessian[i][j];
            }
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

         //Store old coordinates
         for(int i=0;i<molecule.GetNumberAtoms();i++){
            const Atom*   atom = molecule.GetAtom(i);
            const double* xyz  = atom->GetXyz();
            for(int j=0;j<CartesianType_end;j++){
               matrixOldCoordinates[i][j] = xyz[j];
            }
         }

         // Calculate Hessian eigenvalues
         // TODO: Decide whether this routine is necessary
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
                  this->OutputLog((boost::format("\n%e")%tmpvector[i]).str());
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

         //Calculate RFO step
         MallocerFreer::GetInstance()->Malloc(&matrixStep, molecule.GetNumberAtoms(), CartesianType_end);
         vectorStep = &matrixStep[0][0];
         this->CalcRFOStep(vectorStep, matrixHessian, vectorForce, maxNormStep, dimension);

         // Calculate approximate change of energy using
         // [2/2] Pade approximant
         // See Eq. (2) in [BB_1998]
         double approximateChangeNumerator   = 0;
         double approximateChangeDenominator = 1;
         for(int i=0;i<dimension;i++){
            approximateChangeNumerator -= vectorForce[i] * vectorStep[i];
            approximateChangeDenominator += vectorStep[i] * vectorStep[i];
            for(int j=0;j<dimension;j++){
               approximateChangeNumerator += vectorStep[i] * matrixHessian[i][j] * vectorStep[j] / 2;
            }
         }
         double approximateChange = approximateChangeNumerator / approximateChangeDenominator;

         // Take a RFO step
         bool doLineSearch = false;
         bool tempCanOutputLogs = false;
         lineSearchInitialEnergy = lineSearchCurrentEnergy;
         if(doLineSearch){
            this->LineSearch(electronicStructure, molecule, lineSearchCurrentEnergy, matrixStep, elecState, dt);
         }
         else{
            this->UpdateMolecularCoordinates(molecule, matrixStep);
            this->UpdateElectronicStructure(electronicStructure, molecule, requireGuess, tempCanOutputLogs);
            lineSearchCurrentEnergy = electronicStructure->GetElectronicEnergy(elecState);
         }
         this->OutputMoleculeElectronicStructure(electronicStructure, molecule, this->CanOutputLogs());

         // Calculate the correctness of the approximation
         double r = (lineSearchCurrentEnergy - lineSearchInitialEnergy)
                  / approximateChange; // correctness of the step
         bool aproxcheckCanOutputLogs = true;
         tempCanOutputLogs = molecule.CanOutputLogs();
         molecule.SetCanOutputLogs(aproxcheckCanOutputLogs);
         this->OutputLog((boost::format("\n"
                                        "actual energy change          = %e\n"
                                        "expected energy change        = %e\n"
                                        "actual/expected energy change = %f\n")
                                        % (lineSearchCurrentEnergy-lineSearchInitialEnergy)
                                        % approximateChange % r).str());
         molecule.SetCanOutputLogs(tempCanOutputLogs);

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

         // Update the trust radius
         if(r < 0)
         {
            // Rollback molecular geometry
            bool tempCanOutputLogs = molecule.CanOutputLogs();
            bool rollbackCanOutputLogs = true;
            molecule.SetCanOutputLogs(rollbackCanOutputLogs);
            this->OutputLog("Detected hill climbing.\n"
                            "Rolling back molecular geometry.\n");
            for(int i=0;i<molecule.GetNumberAtoms();i++){
               const Atom* atom = molecule.GetAtom(i);
               double*     xyz  = atom->GetXyz();
               for(int j=0;j<CartesianType_end;j++){
                  xyz[j] = matrixOldCoordinates[i][j];
               }
            }
            lineSearchCurrentEnergy = lineSearchInitialEnergy;
            molecule.SetCanOutputLogs(tempCanOutputLogs);
            // and rerun with smaller trust radius
            // without updating Hessian
            maxNormStep /= 4;
            continue;
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
         //Calculate displacement (K_k at Eq. (15) in [SJTO_1983])
         MallocerFreer::GetInstance()->Malloc(&matrixDisplacement, molecule.GetNumberAtoms(), CartesianType_end);
         for(int i=0;i<molecule.GetNumberAtoms();i++){
            const Atom*   atom = molecule.GetAtom(i);
            const double* xyz  = atom->GetXyz();
            for(int j=0;j<CartesianType_end;j++){
               matrixDisplacement[i][j] = xyz[j] - matrixOldCoordinates[i][j];
            }
         }
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

void BFGS::CalcRFOStep(double* vectorStep,
                       double const* const* matrixHessian,
                       double const* vectorForce,
                       const double maxNormStep,
                       const int dimension) const{
   double** matrixAugmentedHessian = NULL;
   double*  vectorEigenValues      = NULL;
   double normStep = 0;
   try{
      double alpha = 1;
      do{
         // Prepare the modified augmented Hessian
         // See Eq. (7) in [EPW_1997]
         MallocerFreer::GetInstance()->Malloc(&matrixAugmentedHessian, dimension+1,dimension+1);
         for(int i=0;i<dimension;i++){
            for(int j=i;j<dimension;j++){
               // H_k/alpha in Eq. (7) in [EPW_1997]
               matrixAugmentedHessian[i][j] = matrixAugmentedHessian[j][i] = matrixHessian[i][j] / alpha;
            }
         }
         // g_k and g_k^t in Eq. (7) in [EPW_1997]
         for(int i=0;i<dimension;i++){
            // note: gradient = -1 * force
            matrixAugmentedHessian[i][dimension] =
               matrixAugmentedHessian[dimension][i] = - vectorForce[i];
         }
         // 0 in Eq. (7) in [EPW_1997]
         matrixAugmentedHessian[dimension][dimension] = 0;

         // Solve eigenvalue problem on the augmented Hessian
         // See Eq. (7) in [EPW_1997]
         MallocerFreer::GetInstance()->Malloc(&vectorEigenValues, dimension+1);
         //TODO: calculate eigenvalues first then calculate only an eigenvector needed
         bool calcEigenVectors = true;
         MolDS_wrappers::Lapack::GetInstance()->Dsyevd(&matrixAugmentedHessian[0],
                                                       &vectorEigenValues[0],
                                                       dimension+1,
                                                       calcEigenVectors);

         // Select a RFO step as the eigenvector whose eivenvalue is the lowest
         for(int i=0;i<dimension;i++){
            // Scale last element of eigenvector to 1/alpha because
            // [vectorStep, 1] is the eigenvector of augmented Hessian.
            // See Eq. (4) in [EPW_1997].
            vectorStep[i] = matrixAugmentedHessian[0][i] / matrixAugmentedHessian[0][dimension] / alpha;
         }
         //
         // Calculate size of the RFO step
         normStep = 0;
         for(int i=0;i<dimension;i++){
            normStep += vectorStep[i] * vectorStep[i];
         }
         normStep = sqrt(normStep);

         this->OutputLog((boost::format("Lowest eigenvalue of the augmented Hessian     = %f\n") % vectorEigenValues[0]).str());
         this->OutputLog((boost::format("2nd lowest eigenvalue of the augmented Hessian = %f\n") % vectorEigenValues[1]).str());
         this->OutputLog((boost::format("3rd lowest eigenvalue of the augmented Hessian = %f\n") % vectorEigenValues[2]).str());
         this->OutputLog((boost::format("Calculated RFO step size                       = %f\n") % normStep).str());

         this->OutputLog((boost::format("Trust radius is %f\n") % maxNormStep).str());
         // Limit the step size to maxNormStep
         if(normStep > maxNormStep){
            alpha *= normStep / maxNormStep * 1.1; // 1.1 is speed up factor
            this->OutputLog((boost::format("Scaling factor is increased to %e.\n") % alpha).str());
            this->OutputLog("Recalculating RFO step...\n");
         }
      }while(normStep > maxNormStep);
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&matrixAugmentedHessian, dimension+1, dimension+1);
      MallocerFreer::GetInstance()->Free(&vectorEigenValues, dimension+1);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free(&matrixAugmentedHessian, dimension+1, dimension+1);
   MallocerFreer::GetInstance()->Free(&vectorEigenValues, dimension+1);
}

void BFGS::UpdateHessian(double **matrixHessian,
                         const int dimension,
                         double const* vectorForce,
                         double const* vectorOldForce,
                         double const* vectorDisplacement) const{
   double const* const K = &vectorDisplacement[0]; // K_k in eq. 15 on [SJTO_1983]
   double*  P    = NULL;                           // P_k in eq. 14 on [SJTO_1983]
   double** PP   = NULL;                           // P_k P_k^T at second term in RHS of Eq. (13) in [SJTO_1983]
   double*  HK   = NULL;                           // H_k K_k at third term on RHS of Eq. (13) in [SJTO_1983]
   double** HKKH = NULL;                           // H_k K_k K_k^T H_k at third term on RHS of Eq. (13) in [SJTO_1983]
   try{
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
            // H_k K_k = (K_k^T H_k)^T because H_k^T = H_k
            HKKH[i][j] = HKKH[j][i] = HK[i] * HK[j];
         }
      }

      // Calculate H_k+1 according to Eq. (13) in [SJTO_1983]
      for(int i=0;i<dimension;i++){
         for(int j=i;j<dimension;j++){
            matrixHessian[i][j]+= (PP[i][j] / PK
                                - HKKH[i][j] / KHK);
            matrixHessian[j][i] = matrixHessian[i][j];
         }
      }
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&P, dimension);
      MallocerFreer::GetInstance()->Free(&PP, dimension, dimension);
      MallocerFreer::GetInstance()->Free(&HK, dimension);
      MallocerFreer::GetInstance()->Free(&HKKH, dimension,dimension);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free(&P, dimension);
   MallocerFreer::GetInstance()->Free(&PP, dimension, dimension);
   MallocerFreer::GetInstance()->Free(&HK, dimension);
   MallocerFreer::GetInstance()->Free(&HKKH, dimension,dimension);
}
}
