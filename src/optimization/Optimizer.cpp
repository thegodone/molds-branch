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
#include"../config.h"
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
#include"../base/factories/ElectronicStructureFactory.h"
#include"Optimizer.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
using namespace MolDS_base_factories;

namespace MolDS_optimization{
Optimizer::Optimizer(){
   this->SetEnableTheoryTypes();
   //this->OutputLog("Optimizer created\n");
}

Optimizer::~Optimizer(){
   //this->OutputLog("Optimizer deleted\n");
}

void Optimizer::Optimize(Molecule& molecule){
   this->OutputLog(this->messageStartGeometryOptimization);
   this->ClearMolecularMomenta(molecule);

   // malloc electornic structure
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   this->CheckEnableTheoryType(theory);
   boost::shared_ptr<ElectronicStructure> electronicStructure(ElectronicStructureFactory::Create());
   electronicStructure->SetMolecule(&molecule);
   electronicStructure->SetCanOutputLogs(this->CanOutputLogs());
   molecule.SetCanOutputLogs(this->CanOutputLogs());

   // Search Minimum
   double lineSearchedEnergy = 0.0;
   bool obtainesOptimizedStructure = false;
   this->SearchMinimum(electronicStructure, molecule, &lineSearchedEnergy, &obtainesOptimizedStructure);
  
   // Not converged
   if(!obtainesOptimizedStructure){
      int totalSteps = Parameters::GetInstance()->GetTotalStepsOptimization();
      stringstream ss;
      ss << this->errorMessageGeometyrOptimizationNotConverged;
      ss << this->errorMessageTotalSteps << totalSteps << endl;
      throw MolDSException(ss.str());
   }
   this->OutputLog(this->messageEndGeometryOptimization);
}

void Optimizer::SetMessages(){
   this->errorMessageTheoryType = "\ttheory type = ";
   this->errorMessageTotalSteps = "\tTotal steps = ";
   this->messageGeometyrOptimizationMetConvergence 
      = "\t\tGeometry otimization met convergence criterion(^^b\n\n\n";
   this->messageStartGeometryOptimization = "**********  START: Geometry optimization  **********\n";
   this->messageEndGeometryOptimization =   "**********  DONE: Geometry optimization  **********\n";
   this->messageReducedTimeWidth = "dt is reduced to ";
   this->messageLineSearchSteps = "\tNumber of Line search steps: ";
   this->messageOptimizationLog = "\t====== Optimization Logs ======\n";
   this->messageEnergyDifference = "\tEnergy difference: ";
   this->messageMaxGradient = "\tMax gradient: ";
   this->messageRmsGradient = "\tRms gradient: ";
   this->messageAu = "[a.u.]";
}

void Optimizer::SetEnableTheoryTypes(){
   this->enableTheoryTypes.clear();
   this->enableTheoryTypes.push_back(ZINDOS);
   this->enableTheoryTypes.push_back(MNDO);
   this->enableTheoryTypes.push_back(AM1);
   this->enableTheoryTypes.push_back(AM1D);
   this->enableTheoryTypes.push_back(PM3);
   this->enableTheoryTypes.push_back(PM3D);
   this->enableTheoryTypes.push_back(PM3PDDG);
}

void Optimizer::CheckEnableTheoryType(TheoryType theoryType) const{
   bool isEnable = false;
   for(int i=0; i<this->enableTheoryTypes.size();i++){
      if(theoryType == this->enableTheoryTypes[i]){
         isEnable = true;
         break;
      }
   }
   if(!isEnable){
      stringstream ss;
      ss << this->errorMessageNotEnebleTheoryType;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theoryType) << endl;
      throw MolDSException(ss.str());
   }
}

void Optimizer::ClearMolecularMomenta(Molecule& molecule) const{
#pragma omp parallel for schedule(dynamic, MOLDS_OMP_DYNAMIC_CHUNK_SIZE) 
   for(int a=0; a<molecule.GetAtomVect().size(); a++){
      const Atom* atom = molecule.GetAtomVect()[a];
      atom->SetPxyz(0.0, 0.0, 0.0);
   }
}

void Optimizer::UpdateMolecularCoordinates(Molecule& molecule, double const* const* matrixForce, double dt) const{
#pragma omp parallel for schedule(dynamic, MOLDS_OMP_DYNAMIC_CHUNK_SIZE) 
   for(int a=0; a<molecule.GetAtomVect().size(); a++){
      const Atom* atom = molecule.GetAtomVect()[a];
      double coreMass = atom->GetCoreMass();
      for(int i=0; i<CartesianType_end; i++){
         atom->GetXyz()[i] += dt*matrixForce[a][i]/coreMass;
      }
   }
   molecule.CalcBasicsConfiguration();
}

void Optimizer::UpdateMolecularCoordinates(Molecule& molecule, double const* const* matrixForce) const{
#pragma omp parallel for schedule(dynamic, MOLDS_OMP_DYNAMIC_CHUNK_SIZE)
   for(int a=0; a<molecule.GetAtomVect().size(); a++){
      const Atom* atom = molecule.GetAtomVect()[a];
      for(int i=0; i<CartesianType_end; i++){
         atom->GetXyz()[i] += matrixForce[a][i];
      }
   }
   molecule.CalcBasicsConfiguration();
}

void Optimizer::UpdateElectronicStructure(boost::shared_ptr<ElectronicStructure> electronicStructure, 
                                          Molecule& molecule,
                                          bool requireGuess, 
                                          bool canOutputLogs) const{
   electronicStructure->SetCanOutputLogs(canOutputLogs);
   molecule.SetCanOutputLogs(canOutputLogs);
   electronicStructure->DoSCF(requireGuess);
   if(Parameters::GetInstance()->RequiresCIS()){
      electronicStructure->DoCIS();
   }
}

void Optimizer::OutputMoleculeElectronicStructure(boost::shared_ptr<ElectronicStructure> electronicStructure, 
                                                  Molecule& molecule,
                                                  bool canOutputLogs) const{
   // output molecular configuration
   molecule.SetCanOutputLogs(canOutputLogs);
   molecule.OutputConfiguration();
   molecule.OutputXyzCOM();
   molecule.OutputXyzCOC();

   // output electornic structure
   electronicStructure->SetCanOutputLogs(canOutputLogs);
   electronicStructure->OutputSCFResults();
   if(Parameters::GetInstance()->RequiresCIS()){
      electronicStructure->OutputCISResults();
   }
}

void Optimizer::LineSearch(boost::shared_ptr<ElectronicStructure> electronicStructure,
                           MolDS_base::Molecule& molecule,
                           double& lineSearchCurrentEnergy,
                           double const* const* matrixForce,
                           int elecState,
                           double dt) const{
   bool tempCanOutputLogs = false;
   int lineSearchSteps = 0;
   double lineSearchOldEnergy = lineSearchCurrentEnergy;

   bool requireGuess = false;
   while(lineSearchCurrentEnergy <= lineSearchOldEnergy){
      this->UpdateMolecularCoordinates(molecule, matrixForce, dt);
      this->UpdateElectronicStructure(electronicStructure, molecule, requireGuess, tempCanOutputLogs);
      lineSearchOldEnergy = lineSearchCurrentEnergy;
      lineSearchCurrentEnergy = electronicStructure->GetElectronicEnergy(elecState);
      lineSearchSteps++;
   }

   // final state of line search
   this->OutputLog(boost::format("%s%d\n\n") % this->messageLineSearchSteps.c_str() % lineSearchSteps);
   this->UpdateMolecularCoordinates(molecule, matrixForce, -0.5*dt);

   // Broadcast to all processes
   int root = MolDS_mpi::MpiProcess::GetInstance()->GetHeadRank();
   molecule.BroadcastConfigurationToAllProcesses(root);

   // update and output electronic structure
   this->UpdateElectronicStructure(electronicStructure, molecule, requireGuess, tempCanOutputLogs);
   this->OutputMoleculeElectronicStructure(electronicStructure, molecule, this->CanOutputLogs());

   lineSearchCurrentEnergy = electronicStructure->GetElectronicEnergy(elecState);
}

bool Optimizer::SatisfiesConvergenceCriterion(double const* const* matrixForce, 
                                              const MolDS_base::Molecule& molecule,
                                              double oldEnergy,
                                              double currentEnergy,
                                              double maxGradientThreshold,
                                              double rmsGradientThreshold) const{
   bool satisfies = false;
   double maxGradient = 0.0;
   double sumSqureGradient = 0.0;
   double energyDifference = currentEnergy - oldEnergy;
   for(int a=0; a<molecule.GetAtomVect().size(); a++){
      for(int i=0; i<CartesianType_end; i++){
         if(maxGradient<fabs(matrixForce[a][i])){
            maxGradient = fabs(matrixForce[a][i]);
         }
         sumSqureGradient += pow(matrixForce[a][i],2.0);
      }
   }
   sumSqureGradient /= static_cast<double>(molecule.GetAtomVect().size()*CartesianType_end);
   double rmsGradient = sqrt(sumSqureGradient);

   // output logs
   this->OutputLog("\n");
   this->OutputLog(this->messageOptimizationLog);
   this->OutputLog(boost::format("%s %e %s\n") % this->messageEnergyDifference.c_str() 
                                               % energyDifference 
                                               % this->messageAu.c_str());
   this->OutputLog(boost::format("%s %e %s\n") % this->messageMaxGradient.c_str() 
                                               % maxGradient 
                                               % this->messageAu.c_str());
   this->OutputLog(boost::format("%s %e %s\n") % this->messageRmsGradient.c_str() 
                                               % rmsGradient 
                                               % this->messageAu.c_str());
   this->OutputLog("\n\n");
  
   // judge convergence
   if(maxGradient < maxGradientThreshold && rmsGradient < rmsGradientThreshold && energyDifference < 0){
      this->OutputLog(this->messageGeometyrOptimizationMetConvergence);
      satisfies = true;
   }
   return satisfies;
}

}



