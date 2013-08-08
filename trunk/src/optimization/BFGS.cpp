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
#include"../base/Uncopyable.h"
#include"../base/PrintController.h"
#include"../base/MolDSException.h"
#include"../base/MallocerFreer.h"
#include"../mpi/MpiProcess.h"
#include"../wrappers/Blas.h"
#include"../wrappers/Lapack.h"
#include"../base/Enums.h"
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
   this->errorMessageNaNInRFOStep
      = "Error in optimization::BFGS::Optimize: RFO step has gone NaN. (lambda * s[%d] = %e, lambda = %e, alpha = %e)\n";

   this->messageStartBFGSStep
      = "\n==========  START: BFGS step ";
   this->messageHillClimbing =
      "Detected hill climbing.\n"
      "Rolling back molecular geometry.\n";
   this->messageRecalculateRFOStep
      = "Recalculating RFO step...\n";
   this->messageRawHessianEigenvalues
      = "Eigenvalues of the raw Hessian:";
   this->messageShiftedHessianEigenvalues
      = "Eigenvalues of the level shifted hessian:";

   this->formatEnergyChangeComparison =
      "\n"
      "actual energy change          = %e\n"
      "expected energy change        = %e\n"
      "actual/expected energy change = %f\n";
   this->formatLowestHessianEigenvalue
      = "Lowest eigenvalue of the augmented Hessian     = %f\n";
   this->format2ndLowestHessianEigenvalue
      = "2nd lowest eigenvalue of the augmented Hessian = %f\n";
   this->format3rdLowestHessianEigenvalue
      = "3rd lowest eigenvalue of the augmented Hessian = %f\n";
   this->formatRFOStepSize
      = "Calculated RFO step size                       = %f\n";
   this->formatTrustRadiusIs
      = "Trust radius is %f\n";
   this->formatIncreaseScalingFactor
      = "Scaling factor is increased to %e.\n";
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
         this->OutputLog(boost::format("%s%d\n\n") % this->messageStartBFGSStep % (s+1));

         // Store old Force data
         MallocerFreer::GetInstance()->Malloc(&vectorOldForce, dimension);
         for(int i =0;i < dimension; i++){
            vectorOldForce[i] = vectorForce[i];
         }

         //Store old coordinates
         MallocerFreer::GetInstance()->Malloc(&matrixOldCoordinates, molecule.GetNumberAtoms(), CartesianType_end);
         for(int i=0;i<molecule.GetNumberAtoms();i++){
            const Atom*   atom = molecule.GetAtom(i);
            const double* xyz  = atom->GetXyz();
            for(int j=0;j<CartesianType_end;j++){
               matrixOldCoordinates[i][j] = xyz[j];
            }
         }
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

void BFGS::CalcRFOStep(double* vectorStep,
                       double const* const* matrixHessian,
                       double const* vectorForce,
                       const double trustRadius,
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
            // See Eq. (7) in [EPW_1997].
            vectorStep[i] = matrixAugmentedHessian[0][i] / matrixAugmentedHessian[0][dimension] / alpha;
            if(isnan(vectorStep[i])){
               throw MolDSException(boost::format(this->errorMessageNaNInRFOStep)
                     % i % matrixAugmentedHessian[0][i] % matrixAugmentedHessian[0][dimension] % alpha);
            }
         }
         //
         // Calculate size of the RFO step
         normStep = 0;
         for(int i=0;i<dimension;i++){
            normStep += vectorStep[i] * vectorStep[i];
         }
         normStep = sqrt(normStep);

         this->OutputLog(boost::format(this->formatLowestHessianEigenvalue)    % vectorEigenValues[0]);
         this->OutputLog(boost::format(this->format2ndLowestHessianEigenvalue) % vectorEigenValues[1]);
         this->OutputLog(boost::format(this->format3rdLowestHessianEigenvalue) % vectorEigenValues[2]);
         this->OutputLog(boost::format(this->formatRFOStepSize)                % normStep);
         this->OutputLog(boost::format(this->formatTrustRadiusIs)              % trustRadius);

         // Limit the step size to trustRadius
         if(normStep > trustRadius){
            alpha *= normStep / trustRadius * 1.1; // 1.1 is speed up factor
            this->OutputLog(boost::format(this->formatIncreaseScalingFactor) % alpha);
            this->OutputLog(this->messageRecalculateRFOStep);
         }
      }while(normStep > trustRadius);
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
   double *P  = NULL;                              // P_k in eq. 14 on [SJTO_1983]
   double *HK = NULL;                              // H_k K_k at third term on RHS of Eq. (13) in [SJTO_1983]
   try{
      MallocerFreer::GetInstance()->Malloc(&P, dimension);
      MallocerFreer::GetInstance()->Malloc(&HK, dimension);
      double KHK = 0;
      double PK = 0;
      // initialize P_k according to Eq. (14) in [SJTO_1983]
      // note: gradient = -1 * force
      MolDS_wrappers::Blas::GetInstance()->Dcopy(dimension, vectorOldForce, P);
      MolDS_wrappers::Blas::GetInstance()->Daxpy(dimension, -1.0, vectorForce, P);

      // P_k^T K_k at second term at RHS of Eq. (13) in [SJTO_1983]
      PK = MolDS_wrappers::Blas::GetInstance()->Ddot(dimension, P, K);

      //H_k K_k at third term on RHS of Eq. (13) in [SJTO_1983]
      MolDS_wrappers::Blas::GetInstance()->Dsymv(dimension, matrixHessian, K, HK);

      //K_k^T H_k K_k at third term on RHS of Eq. (13) in [SJTO_1983]
      KHK = MolDS_wrappers::Blas::GetInstance()->Ddot(dimension, K, HK);

      // Calculate H_k+1 according to Eq. (13) in [SJTO_1983]
      // Add second term in RHS of Eq. (13) in [SJTO_1983]
      MolDS_wrappers::Blas::GetInstance()->Dsyr(dimension,  1.0/PK,  P,  matrixHessian);
      // Add third term in RHS of Eq. (13) in [SJTO_1983]
      MolDS_wrappers::Blas::GetInstance()->Dsyr(dimension, -1.0/KHK, HK, matrixHessian);
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&P, dimension);
      MallocerFreer::GetInstance()->Free(&HK, dimension);
      throw ex;
   }
   MallocerFreer::GetInstance()->Free(&P, dimension);
   MallocerFreer::GetInstance()->Free(&HK, dimension);
}

// Level shift eigenvalues of redandant modes to largeEigenvalue

void BFGS::ShiftHessianRedundantMode(double** matrixHessian,
                                     const Molecule& molecule) const{
   const double one                      = 1;
   const double largeEigenvalue          = 1.0e3;
   const int    numAtoms                 = molecule.GetNumberAtoms();
   const int    dimension                = numAtoms *CartesianType_end;
   const int    numTranslationalModes    = 3;
   const int    numRotationalModes       = 3;
   int          numRedundantModes        = numTranslationalModes + numRotationalModes;
   double** vectorsHessianModes          = NULL;
   double*  vectorHessianEigenValues     = NULL;
   double** matrixesRedundantModes[]     = {NULL, NULL, NULL, NULL, NULL, NULL};
   double*  vectorsRedundantModes[]      = {NULL, NULL, NULL, NULL, NULL, NULL};
   double** matrixProjection             = NULL;
   double*  vectorProjectedRedundantMode = NULL;
   double** matrixShiftedHessianBuffer   = NULL;
   const double matrixesRotationalModeGenerators[numRotationalModes]
                                                [CartesianType_end]
                                                [CartesianType_end]
                = {{{0,  0, 0}, {0, 0, -1}, { 0, 1, 0}},
                   {{0,  0, 1}, {0, 0,  0}, {-1, 0, 0}},
                   {{0, -1, 0}, {1, 0,  0}, { 0, 0, 0}}};

   try{
      // Prepare translational modes
      for(int c=0; c<numTranslationalModes;c++){
         MallocerFreer::GetInstance()->Malloc(&matrixesRedundantModes[c], numAtoms, CartesianType_end);
         vectorsRedundantModes[c] = &matrixesRedundantModes[c][0][0];
      }
      for(int c=0; c<numTranslationalModes;c++){
#pragma omp parallel for schedule(auto)
         for(int n=0;n<numAtoms;n++){
            for(int d=0;d<CartesianType_end;d++){
               matrixesRedundantModes[c][n][d] = c==d? 1.0 : 0.0;
            }
         }
      }
      // Prepare rotational modes
      for(int c=0; c<numRotationalModes;c++){
         MallocerFreer::GetInstance()->Malloc(&matrixesRedundantModes[c+numTranslationalModes], numAtoms, CartesianType_end);
         vectorsRedundantModes[c+numTranslationalModes] = &matrixesRedundantModes[c+numTranslationalModes][0][0];
      }
      for(int c=0; c<numRotationalModes;c++){
#pragma omp parallel for schedule(auto)
         for(int n=0;n<numAtoms;n++){
            const double* xyz = molecule.GetAtom(n)->GetXyz();
            for(int d=0;d<CartesianType_end;d++){
               matrixesRedundantModes[c+numTranslationalModes][n][d] = 0.0;
               for(int e=0;e<CartesianType_end;e++){
                  matrixesRedundantModes[c+numTranslationalModes][n][d] += matrixesRotationalModeGenerators[c][d][e] * xyz[e];
               }
            }
         }
      }
      //Prepare the identity matrix I
      MallocerFreer::GetInstance()->Malloc(&matrixProjection, dimension, dimension);
      MolDS_wrappers::Blas::GetInstance()->Dcopy(dimension, &one, 0, &matrixProjection[0][0], dimension+1);

      //Prepare the projection matrix P = I - sum_i u_i u_i^T, that projects redundant modes to 0 vector
      MallocerFreer::GetInstance()->Malloc(&vectorProjectedRedundantMode, dimension);
      for(int c=0; c<numRedundantModes; c++){
         double normSquare = 0;
         MolDS_wrappers::Blas::GetInstance()->Dsymv(dimension, matrixProjection, vectorsRedundantModes[c], vectorProjectedRedundantMode);
         normSquare = MolDS_wrappers::Blas::GetInstance()->Ddot(dimension, vectorProjectedRedundantMode, vectorProjectedRedundantMode);
         MolDS_wrappers::Blas::GetInstance()->Dsyr(dimension, -1.0/normSquare, vectorProjectedRedundantMode, matrixProjection);
      }

      //// Diagonalize hessian
      //MallocerFreer::GetInstance()->Malloc(&vectorHessianEigenValues, dimension);
      //MallocerFreer::GetInstance()->Malloc(&vectorsHessianModes, dimension, dimension);
      //MolDS_wrappers::Blas::GetInstance()->Dcopy(dimension*dimension, &matrixHessian[0][0], &vectorsHessianModes[0][0]);
      //bool calcEigenVectors = true;
      //MolDS_wrappers::Lapack::GetInstance()->Dsyevd(&vectorsHessianModes[0],
      //                                              &vectorHessianEigenValues[0],
      //                                              dimension,
      //                                              calcEigenVectors);
      //
      //// Output eigenvalues of the raw Hessian to the log
      //this->OutputLog(this->messageRawHessianEigenvalues);
      //for(int i=0;i<dimension;i++){
      //   if((i%6) == 0){
      //      this->OutputLog(boost::format("\n%e")%vectorHessianEigenValues[i]);
      //   }
      //   else{
      //      this->OutputLog(boost::format(",\t%e")%vectorHessianEigenValues[i]);
      //   }
      //}
      //this->OutputLog("\n");

      // Project Hessian H' = P H P
      MallocerFreer::GetInstance()->Malloc(&matrixShiftedHessianBuffer, dimension, dimension);
      // TODO: Use dsymm instead of dgemm
      MolDS_wrappers::Blas::GetInstance()->Dgemm(dimension, dimension, dimension,
            matrixProjection            , matrixHessian   , matrixShiftedHessianBuffer);
      MolDS_wrappers::Blas::GetInstance()->Dgemm(dimension, dimension, dimension,
            matrixShiftedHessianBuffer, matrixProjection, matrixHessian);

      // Shift eigenvalues for redundant mode by adding L sum_i u_i*u_i^T = L ( I - P )= -L (P - I)
      MolDS_wrappers::Blas::GetInstance()->Daxpy(dimension, -1, &one, 0, &matrixProjection[0][0], dimension + 1); // (P - I)
      MolDS_wrappers::Blas::GetInstance()->Daxpy(dimension*dimension, -largeEigenvalue, &matrixProjection[0][0], &matrixHessian[0][0]);

      //// Diagonalize shifted hessian
      //MolDS_wrappers::Blas::GetInstance()->Dcopy(dimension*dimension, &matrixHessian[0][0], &vectorsHessianModes[0][0]);
      //calcEigenVectors = true;
      //MolDS_wrappers::Lapack::GetInstance()->Dsyevd(&vectorsHessianModes[0],
      //                                              &vectorHessianEigenValues[0],
      //                                              dimension,
      //                                              calcEigenVectors);
      //
      //// Output eigenvalues of the shifted Hessian to the log
      //this->OutputLog(this->messageShiftedHessianEigenvalues);
      //for(int i=0;i<dimension;i++){
      //   if((i%6) == 0){
      //      this->OutputLog(boost::format("\n%e")%vectorHessianEigenValues[i]);
      //   }
      //   else{
      //      this->OutputLog(boost::format(",\t%e")%vectorHessianEigenValues[i]);
      //   }
      //}
      //this->OutputLog("\n");
   }
   catch(MolDSException ex)
   {
      for(int i=0;i<numRedundantModes;i++){
         MallocerFreer::GetInstance()->Free(&matrixesRedundantModes[i], numAtoms, CartesianType_end);
      }
      MallocerFreer::GetInstance()->Free(&matrixProjection, dimension, dimension);
      MallocerFreer::GetInstance()->Free(&vectorProjectedRedundantMode, dimension);
      MallocerFreer::GetInstance()->Free(&matrixShiftedHessianBuffer, dimension, dimension);
      MallocerFreer::GetInstance()->Free(&vectorHessianEigenValues, dimension);
      MallocerFreer::GetInstance()->Free(&vectorsHessianModes, dimension, dimension);
   }
   for(int i=0;i<numRedundantModes;i++){
      MallocerFreer::GetInstance()->Free(&matrixesRedundantModes[i], numAtoms, CartesianType_end);
   }
   MallocerFreer::GetInstance()->Free(&matrixProjection, dimension, dimension);
   MallocerFreer::GetInstance()->Free(&vectorProjectedRedundantMode, dimension);
   MallocerFreer::GetInstance()->Free(&matrixShiftedHessianBuffer, dimension, dimension);
   MallocerFreer::GetInstance()->Free(&vectorHessianEigenValues, dimension);
   MallocerFreer::GetInstance()->Free(&vectorsHessianModes, dimension, dimension);
}

double BFGS::ApproximateEnergyChange(int dimension,
                                     double const* const* matrixHessian,
                                     double const* vectorForce,
                                     double const* vectorStep) const{
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
   return approximateChangeNumerator / approximateChangeDenominator;
}

void BFGS::UpdateTrustRadius(double &trustRadius,
                       double approximateEnergyChange,
                       double initialEnergy,
                       double currentEnergy)const{
   // Calculate the correctness of the approximation
   double r = (currentEnergy - initialEnergy)
            / approximateEnergyChange;

   this->OutputLog(boost::format(this->formatEnergyChangeComparison)
                   % (currentEnergy-initialEnergy) % approximateEnergyChange % r);

   if(r < 0)
   {
      trustRadius /= 4;
   }
   else if(r<0.25){
      trustRadius /= 4;
   }
   else if(r<0.75){
      // keep trust radius
   }
   else if(r<2){
      trustRadius *= 2;
   }
   else{
      trustRadius /= 2;
   }
}

void BFGS::RollbackMolecularGeometry(MolDS_base::Molecule& molecule,
                                     double const* const* matrixOldCoordinates)const{
   // Rollback molecular geometry
   bool tempCanOutputLogs = molecule.CanOutputLogs();
   bool rollbackCanOutputLogs = true;
   molecule.SetCanOutputLogs(rollbackCanOutputLogs);
   this->OutputLog(this->messageHillClimbing);
   for(int i=0;i<molecule.GetNumberAtoms();i++){
      const Atom* atom = molecule.GetAtom(i);
      double*     xyz  = atom->GetXyz();
      for(int j=0;j<CartesianType_end;j++){
         xyz[j] = matrixOldCoordinates[i][j];
      }
   }
   molecule.SetCanOutputLogs(tempCanOutputLogs);
}

void BFGS::CalcDisplacement(double      *      *& matrixDisplacement,
                            double const* const*  matrixOldCoordinates,
                            const MolDS_base::Molecule& molecule)const{
   //Calculate displacement (K_k at Eq. (15) in [SJTO_1983])
   MallocerFreer::GetInstance()->Malloc(&matrixDisplacement, molecule.GetNumberAtoms(), CartesianType_end);
   for(int i=0;i<molecule.GetNumberAtoms();i++){
      const Atom*   atom = molecule.GetAtom(i);
      const double* xyz  = atom->GetXyz();
      for(int j=0;j<CartesianType_end;j++){
         matrixDisplacement[i][j] = xyz[j] - matrixOldCoordinates[i][j];
      }
   }
}
}
