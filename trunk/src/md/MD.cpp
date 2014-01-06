//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   //
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
#include"MD.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
using namespace MolDS_base_factories;

namespace MolDS_md{
MD::MD(){
   this->molecule = NULL;
   this->SetMessages();
   this->SetEnableTheoryTypes();
   //this->OutputLog("MD created \n");
}

MD::~MD(){
   //this->OutputLog("MD deleted\n");
}

void MD::SetMolecule(Molecule* molecule){
   // check enable electonic theory
   this->molecule = molecule;
}

void MD::DoMD(){
   this->OutputLog(this->messageStartMD);

   // malloc electornic structure
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   this->CheckEnableTheoryType(theory);
   boost::shared_ptr<ElectronicStructure> electronicStructure(ElectronicStructureFactory::Create());
   electronicStructure->SetMolecule(this->molecule);
   electronicStructure->SetCanOutputLogs(this->CanOutputLogs());
   this->molecule->SetCanOutputLogs(this->CanOutputLogs());

   int totalSteps       = Parameters::GetInstance()->GetTotalStepsMD();
   int elecState        = Parameters::GetInstance()->GetElectronicStateIndexMD();
   double dt            = Parameters::GetInstance()->GetTimeWidthMD();
   double time          = 0.0;
   bool requireGuess    = false;
   double initialEnergy = 0.0;
   double const* const* matrixForce = NULL;

   // initial calculation
   electronicStructure->DoSCF();
   if(Parameters::GetInstance()->RequiresCIS()){
      electronicStructure->DoCIS();
   }
   matrixForce = electronicStructure->GetForce(elecState);

   // output initial conditions
   this->OutputLog(this->messageinitialConditionMD);
   initialEnergy = this->OutputEnergies(electronicStructure);
   this->OutputLog("\n");
   this->molecule->OutputConfiguration();
   this->molecule->OutputXyzCOM();
   this->molecule->OutputXyzCOC();
   this->molecule->OutputMomenta();

   for(int s=0; s<totalSteps; s++){
      this->OutputLog(boost::format("%s%d\n") % this->messageStartStepMD.c_str() % (s+1) );

      // update momenta & coordinates
      this->UpdateMomenta    (*this->molecule, matrixForce, dt);
      this->UpdateCoordinates(*this->molecule, dt);

      // update electronic structure
      electronicStructure->DoSCF(requireGuess);
      if(Parameters::GetInstance()->RequiresCIS()){
         electronicStructure->DoCIS();
      }

      // update force
      matrixForce = electronicStructure->GetForce(elecState);

      // update momenta
      this->UpdateMomenta(*this->molecule, matrixForce, dt);

      // Broadcast to all processes
      int root = MolDS_mpi::MpiProcess::GetInstance()->GetHeadRank();
      this->molecule->BroadcastPhaseSpacePointToAllProcesses(root);

      // output results
      this->OutputEnergies(electronicStructure, initialEnergy);
      this->molecule->OutputConfiguration();
      this->molecule->OutputXyzCOM();
      this->molecule->OutputXyzCOC();
      this->molecule->OutputMomenta();
      this->OutputLog(boost::format("%s%lf\n") % this->messageTime.c_str() 
                                               % (dt*static_cast<double>(s+1)/Parameters::GetInstance()->GetFs2AU()));
      this->OutputLog(boost::format("%s%d\n") % this->messageEndStepMD.c_str() % (s+1) );
   }

   this->OutputLog(this->messageEndMD);
}

void MD::UpdateMomenta(const Molecule& molecule, double const* const* matrixForce, double dt) const{
#pragma omp parallel for schedule(dynamic, MOLDS_OMP_DYNAMIC_CHUNK_SIZE)
   for(int a=0; a<molecule.GetAtomVect().size(); a++){
      Atom* atom = molecule.GetAtomVect()[a];
      for(int i=0; i<CartesianType_end; i++){
         atom->GetPxyz()[i] += 0.5*dt*(matrixForce[a][i]);
      }
   }
}

void MD::UpdateCoordinates(Molecule& molecule, double dt) const{
#pragma omp parallel for schedule(dynamic, MOLDS_OMP_DYNAMIC_CHUNK_SIZE)
      for(int a=0; a<molecule.GetAtomVect().size(); a++){
         Atom* atom = molecule.GetAtomVect()[a];
         double coreMass = atom->GetCoreMass();
         for(int i=0; i<CartesianType_end; i++){
            atom->GetXyz()[i] += dt*atom->GetPxyz()[i]/coreMass;
         }
      }
      molecule.CalcBasicsConfiguration();
}

void MD::SetMessages(){
   this->errorMessageTheoryType = "\ttheory type = ";
   this->errorMessageNotEnebleTheoryType  
      = "Error in md::MD::CheckEnableTheoryType: Non available theory is set.\n";
   this->messageStartMD = "**********  START: Molecular dynamics  **********\n";
   this->messageEndMD = "**********  DONE: Molecular dynamics  **********\n";
   this->messageinitialConditionMD = "\n\t========= Initial conditions \n";
   this->messageStartStepMD = "\n\t========== START: MD step ";
   this->messageEndStepMD =     "\t========== DONE: MD step ";
   this->messageEnergies = "\tEnergies:\n";
   this->messageEnergiesTitle = "\t\t|\tkind\t\t\t| [a.u.] | [eV] | \n";
   this->messageCoreKineticEnergy =   "Core kinetic:     ";
   this->messageCoreRepulsionEnergy = "Core repulsion:   ";
   this->messageVdWCorrectionEnergy = "VdW correction:   ";
   this->messageElectronicEnergy = "Electronic\n\t\t(inc. core rep.):";
   this->messageElectronicEnergyVdW = "Electronic\n\t\t(inc. core rep. and vdW):";
   this->messageTotalEnergy =         "Total:            ";
   this->messageErrorEnergy =         "Error:            ";
   this->messageTime = "\tTime in [fs]: ";
}

double MD::OutputEnergies(boost::shared_ptr<ElectronicStructure> electronicStructure){
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexMD();
   double eV2AU  = Parameters::GetInstance()->GetEV2AU();
   double coreKineticEnergy = 0.0;
   for(int a=0; a<this->molecule->GetAtomVect().size(); a++){
      Atom* atom = this->molecule->GetAtomVect()[a];
      double coreMass = atom->GetCoreMass();
      for(int i=0; i<CartesianType_end; i++){
         coreKineticEnergy += 0.5*pow(atom->GetPxyz()[i],2.0)/coreMass;
      }
   }  
   // output energies:
   this->OutputLog(this->messageEnergies);
   this->OutputLog(this->messageEnergiesTitle);
   this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageCoreKineticEnergy.c_str()
                                                     % coreKineticEnergy
                                                     % (coreKineticEnergy/eV2AU));
   this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageCoreRepulsionEnergy.c_str()
                                                     % electronicStructure->GetCoreRepulsionEnergy()
                                                     % (electronicStructure->GetCoreRepulsionEnergy()/eV2AU));
   if(Parameters::GetInstance()->RequiresVdWSCF()){
      this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageVdWCorrectionEnergy.c_str()
                                                        % electronicStructure->GetVdWCorrectionEnergy()
                                                        % (electronicStructure->GetVdWCorrectionEnergy()/eV2AU));
      this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageElectronicEnergyVdW.c_str()
                                                        % electronicStructure->GetElectronicEnergy(elecState)
                                                        % (electronicStructure->GetElectronicEnergy(elecState)/eV2AU));
   }
   else{
      this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageElectronicEnergy.c_str()
                                                        % electronicStructure->GetElectronicEnergy(elecState)
                                                        % (electronicStructure->GetElectronicEnergy(elecState)/eV2AU));
   }
   this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageTotalEnergy.c_str()
                                                     % (coreKineticEnergy + electronicStructure->GetElectronicEnergy(elecState))
                                                     % ((coreKineticEnergy + electronicStructure->GetElectronicEnergy(elecState))/eV2AU));
   return (coreKineticEnergy + electronicStructure->GetElectronicEnergy(elecState));
}

void MD::OutputEnergies(boost::shared_ptr<ElectronicStructure> electronicStructure, 
                        double initialEnergy){
   double energy = this->OutputEnergies(electronicStructure);
   double eV2AU = Parameters::GetInstance()->GetEV2AU();
   this->OutputLog(boost::format("\t\t%s\t%e\t%e\n\n") % this->messageErrorEnergy.c_str()
                                                       % (initialEnergy - energy)
                                                       % ((initialEnergy - energy)/eV2AU));
}

void MD::SetEnableTheoryTypes(){
   this->enableTheoryTypes.clear();
   this->enableTheoryTypes.push_back(ZINDOS);
   this->enableTheoryTypes.push_back(MNDO);
   this->enableTheoryTypes.push_back(AM1);
   this->enableTheoryTypes.push_back(AM1D);
   this->enableTheoryTypes.push_back(PM3);
   this->enableTheoryTypes.push_back(PM3D);
   this->enableTheoryTypes.push_back(PM3PDDG);
}

void MD::CheckEnableTheoryType(TheoryType theoryType){

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

}



