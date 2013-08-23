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
#include"../base/RealSphericalHarmonicsIndex.h"
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

GEDIIS::GEDIISHistory::GEDIISHistory(){
   this->SetMessages();
}

GEDIIS::GEDIISHistory::~GEDIISHistory(){
   for(entryList_t::iterator i = this->entryList.begin(); i != this->entryList.end(); i++){
     delete *i;
   }
}

void GEDIIS::GEDIISHistory::SetMessages(){
   this->errorMessageNegativeGEDIISCoefficient
      = "GEDIIS coefficients contains negative value.";
   this->errorMessageNotSufficientHistory
      = "GEDIIS history is not sufficient.";
}

void GEDIIS::GEDIISHistory::AddEntry(double energy,
                                     const MolDS_base::Molecule& molecule,
                                     double const* const* matrixForce){
   this->entryList.push_back(new Entry(energy, molecule, matrixForce));
}

void GEDIIS::GEDIISHistory::DiscardEntries(){
   this->entryList.clear();
}

GEDIIS::GEDIISHistory::Entry::Entry(double energy,
                                    const MolDS_base::Molecule& molecule,
                                    double const* const* matrixForce):
   energy(energy),numAtoms(molecule.GetNumberAtoms()),matrixCoordinate(NULL),matrixForce(NULL) {
   MallocerFreer::GetInstance()->Malloc(&this->matrixCoordinate, this->numAtoms, CartesianType_end);
   MallocerFreer::GetInstance()->Malloc(&this->matrixForce,      this->numAtoms, CartesianType_end);
#pragma omp parallel for schedule(auto)
   for(int i = 0; i < this->numAtoms; i++){
      const Atom*   atom = molecule.GetAtom(i);
      const double* xyz  = atom->GetXyz();
      for(int j = 0; j < CartesianType_end; j++){
         this->matrixCoordinate[i][j] = xyz[j];
         this->matrixForce[i][j]      = matrixForce[i][j];
      }
   }
}

GEDIIS::GEDIISHistory::Entry::~Entry(){
   MallocerFreer::GetInstance()->Free(&this->matrixCoordinate, this->numAtoms, CartesianType_end);
   MallocerFreer::GetInstance()->Free(&this->matrixForce,      this->numAtoms, CartesianType_end);
}

void GEDIIS::GEDIISHistory::SolveGEDIISEquation(double* gediisEnergy, double** matrixCoordinate, double** matrixForce){
   double**  gediisMatrix = NULL;
   double*   gediisCoeffs = NULL;
   double*   bufForce     = NULL;
   double*   bufCoord     = NULL;
   const int numCoeffs    = this->entryList.size();
   const int size         = numCoeffs + 1;
   const int numAtoms     = this->entryList.front()->GetNumberAtoms();
   const int dimension    = numAtoms * CartesianType_end;
   typedef entryList_t::iterator iter;

   if(numCoeffs <= 1){
      MolDSException ex(this->errorMessageNotSufficientHistory);
      ex.SetKeyValue<int>(GEDIISErrorID, GEDIISNotSufficientHistory);
      throw ex;
   }

   MallocerFreer::GetInstance()->Malloc(&gediisMatrix, size, size);
   MallocerFreer::GetInstance()->Malloc(&gediisCoeffs, size);
   MallocerFreer::GetInstance()->Malloc(&bufForce, dimension);
   MallocerFreer::GetInstance()->Malloc(&bufCoord, dimension);
   try{
      iter it1 = this->entryList.begin();
      for(int i = 0; it1 != this->entryList.end(); it1++,i++){
         const Entry* entry1 = *it1;
         gediisCoeffs[i] = entry1->GetEnergy();
         iter it2 = it1;
         for(int j = i; it2 != this->entryList.end(); it2++, j++){
            const Entry* entry2 = *it2;
            MolDS_wrappers::Blas::GetInstance()->Dcopy(dimension,       &entry1->GetForce()[0][0],      &bufForce[0]);
            MolDS_wrappers::Blas::GetInstance()->Dcopy(dimension,       &entry1->GetCoordinate()[0][0], &bufCoord[0]);
            MolDS_wrappers::Blas::GetInstance()->Daxpy(dimension, -1.0, &entry2->GetForce()[0][0],      &bufForce[0]);
            MolDS_wrappers::Blas::GetInstance()->Daxpy(dimension, -1.0, &entry2->GetCoordinate()[0][0], &bufCoord[0]);
            gediisMatrix[i][j] = gediisMatrix[j][i] =
               - MolDS_wrappers::Blas::GetInstance()->Ddot(dimension, bufCoord, bufForce);
         }
         gediisMatrix[i][size-1] = gediisMatrix[size-1][i] = 1;
      }
      gediisMatrix[size-1][size-1] = 0;
      gediisCoeffs[size-1]         = 1;

      MolDS_wrappers::Lapack::GetInstance()->Dsysv(gediisMatrix, gediisCoeffs, size);

      MallocerFreer::GetInstance()->Initialize(matrixCoordinate, numAtoms, CartesianType_end);
      it1 = this->entryList.begin();
      for(int i = 0; it1 != this->entryList.end(); it1++,i++){
         if(gediisCoeffs[i]<0){
//            delete *it1;
//            this->entryList.erase(it1);
            MolDSException ex(this->errorMessageNegativeGEDIISCoefficient);
            ex.SetKeyValue<int>(GEDIISErrorID, GEDIISNegativeCoefficient);
            throw ex;
         }
         MolDS_wrappers::Blas::GetInstance()->Daxpy(dimension, gediisCoeffs[i], &(*it1)->GetCoordinate()[0][0], &matrixCoordinate[0][0]);
         MolDS_wrappers::Blas::GetInstance()->Daxpy(dimension, gediisCoeffs[i], &(*it1)->GetForce()[0][0],      &matrixForce[0][0]);
      }
      *gediisEnergy = gediisCoeffs[numCoeffs];
   }
   catch(MolDSException ex){
      MallocerFreer::GetInstance()->Free(&gediisMatrix, size, size);
      MallocerFreer::GetInstance()->Free(&gediisCoeffs, size);
      MallocerFreer::GetInstance()->Free(&bufForce, dimension);
      MallocerFreer::GetInstance()->Free(&bufCoord, dimension);
      if(ex.HasKey(LapackInfo)){
         ex.SetKeyValue<int>(GEDIISErrorID, GEDIISLapackInfo);
      }
      throw ex;
   }
   MallocerFreer::GetInstance()->Free(&gediisMatrix, size, size);
   MallocerFreer::GetInstance()->Free(&gediisCoeffs, size);
   MallocerFreer::GetInstance()->Free(&bufForce, dimension);
   MallocerFreer::GetInstance()->Free(&bufCoord, dimension);
}

}
