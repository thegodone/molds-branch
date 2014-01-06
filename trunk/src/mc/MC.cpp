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
#include<boost/random.hpp>
#include<boost/format.hpp>
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
#include"MC.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
using namespace MolDS_base_factories;

namespace MolDS_mc{
MC::MC(){
   this->molecule = NULL;
   this->SetMessages();
   //this->OutputLog("MC created \n");
}

MC::~MC(){
   //this->OutputLog("MC deleted\n");
}

void MC::SetMolecule(Molecule* molecule){
   this->molecule = molecule;
}

void MC::SetMessages(){
   this->messageStartMC = "**********  START: Monte Carlo  **********\n";
   this->messageEndMC = "**********  DONE: Monte Carlo  **********\n\n";
   this->messageinitialConditionMC = "\n\t========= Initial conditions \n";
   this->messageStartStepMC = "\n\t========== START: MC step ";
   this->messageEndStepMC =     "\t========== DONE: MC step ";
   this->messageEnergies = "\tEnergies:\n";
   this->messageEnergiesTitle = "\t\t|\tkind\t\t\t| [a.u.] | [eV] | \n";
   this->messageCoreRepulsionEnergy = "Core repulsion:   ";
   this->messageVdWCorrectionEnergy = "VdW correction:   ";
   this->messageElectronicEnergy = "Electronic\n\t\t(inc. core rep.):";
   this->messageElectronicEnergyVdW = "Electronic\n\t\t(inc. core rep. and vdW):";
   this->messageTransitionRate = "\tTransition Rate: ";
}

void MC::DoMC(){
   int totalSteps = Parameters::GetInstance()->GetTotalStepsMC();
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexMC();
   double temperature = Parameters::GetInstance()->GetTemperatureMC();
   double stepWidth = Parameters::GetInstance()->GetStepWidthMC();
   unsigned long seed = Parameters::GetInstance()->GetSeedMC();
   this->DoMC(totalSteps, elecState, temperature, stepWidth, seed);
}

void MC::DoMC(int totalSteps, int elecState, double temperature, double stepWidth, unsigned long seed){
   this->OutputLog(this->messageStartMC);
   double transitionRate = 0.0;
   // create real random generator
	boost::mt19937 realGenerator(seed);
	boost::uniform_real<> range(0, 1);
	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > realRand( realGenerator, range );

   // prepare trial molecule and electronic structure pointa
   Molecule trialMolecule(*this->molecule);
   boost::shared_ptr<ElectronicStructure> electronicStructure2(ElectronicStructureFactory::Create());
   ElectronicStructure* trialES = electronicStructure2.get();
   trialES->SetMolecule(&trialMolecule);
   trialES->SetCanOutputLogs(this->CanOutputLogs());
   trialMolecule.SetCanOutputLogs(this->CanOutputLogs());

   // initial calculation
   boost::shared_ptr<ElectronicStructure> electronicStructure1(ElectronicStructureFactory::Create());
   ElectronicStructure* currentES = electronicStructure1.get();
   currentES->SetMolecule(this->molecule);
   currentES->SetCanOutputLogs(this->CanOutputLogs());
   this->molecule->SetCanOutputLogs(this->CanOutputLogs());
   currentES->DoSCF();
   if(Parameters::GetInstance()->RequiresCIS()){
      currentES->DoCIS();
   }
   this->OutputLog(this->messageinitialConditionMC);
   this->OutputMolecule(*currentES, *this->molecule, elecState);

   // Monte Carlo loop 
   for(int s=0; s<totalSteps; s++){
      this->OutputLog(boost::format("%s%d\n\n") % this->messageStartStepMC.c_str() % (s+1) );
      // create trial molecule
      this->CreateTrialConfiguration(&trialMolecule, *this->molecule, &realRand, stepWidth);
     
      // calculate trilal electronic structure
      bool requireGuess = (s==0) ? true : false;
      trialES->DoSCF(requireGuess);
      if(Parameters::GetInstance()->RequiresCIS()){
         trialES->DoCIS();
      }

      // which Electronic Structure is used?
      if(UsesTrial(*currentES, *trialES, elecState, &realRand, temperature)){
         this->molecule->SynchronizeConfigurationTo(trialMolecule);
         swap(currentES, trialES);
         currentES->SetMolecule(this->molecule);
         trialES->SetMolecule(&trialMolecule);
         transitionRate += 1.0;
      }
      else{
         trialMolecule.SynchronizeConfigurationTo(*this->molecule);
      }

      // Broadcast to all processes
      int root = MolDS_mpi::MpiProcess::GetInstance()->GetHeadRank();
      this->molecule->BroadcastConfigurationToAllProcesses(root);
      trialMolecule.BroadcastConfigurationToAllProcesses(root);
      
      // output molecular states
      this->OutputMolecule(*currentES, *this->molecule, elecState);
      this->OutputLog(boost::format("%s%d\n\n") % this->messageEndStepMC.c_str() % (s+1) );
   }
   this->OutputLog(boost::format("%s%lf\n\n") % this->messageTransitionRate.c_str() % (transitionRate/static_cast<double>(totalSteps)) );
   this->OutputLog(this->messageEndMC);
}

void MC::CreateTrialConfiguration(Molecule* trial,
                                  const Molecule& current,
                                  boost::variate_generator<
                                     boost::mt19937&,
                                     boost::uniform_real<>
                                  > (*realRand),
                                  double stepWidth) const{
   // disturb an atom in trial molecule
   int movedAtomIndex = static_cast<int>((*realRand)()*this->molecule->GetAtomVect().size());
   const Atom& reffAtom = *current.GetAtomVect()[movedAtomIndex];
   Atom* trialAtom = trial->GetAtomVect()[movedAtomIndex];
   double dr[CartesianType_end] = {0.0, 0.0, 0.0};
   for(int i=0; i<CartesianType_end; i++){
      dr[i] = stepWidth*(2.0*(*realRand)() -1.0);
      trialAtom->GetXyz()[i] = reffAtom.GetXyz()[i] + dr[i];
   }

   // shift all atoms in trial molecule, namely shift the center of core.
   double trialAtomCoreMass = reffAtom.GetCoreMass();
   double totalCoreMass = current.GetTotalCoreMass();
   double coreCenterShift[CartesianType_end] = {0.0, 0.0, 0.0};
   for(int i=0; i<CartesianType_end; i++){
      coreCenterShift[i] = dr[i]*trialAtomCoreMass/totalCoreMass;
   }
   for(int a=0; a<current.GetAtomVect().size(); a++){
      Atom* trialAtom = trial->GetAtomVect()[a];
      for(int i=0; i<CartesianType_end; i++){
         trialAtom->GetXyz()[i] -= coreCenterShift[i];
      }
   }
   trial->CalcBasicsConfiguration();
}

bool MC::UsesTrial(const ElectronicStructure& currentES, 
                   const ElectronicStructure& trialES,
                   int elecState,
                   boost::variate_generator<
                     boost::mt19937&,
                     boost::uniform_real<>
                   > (*realRand),
                   double temperature) const{
   double currentElecEne = currentES.GetElectronicEnergy(elecState);
   double trialElecEne = trialES.GetElectronicEnergy(elecState);
   double deltaElecEne = trialElecEne - currentElecEne;
   if(deltaElecEne <= 0.0){
      return true;
   }
   else{
      double kB = Parameters::GetInstance()->GetBoltzmann();
      double p = exp(-1.0*deltaElecEne/(kB*temperature));
      double random = (*realRand)();
      if(p>random){
         return true;
      }
      else{
         return false;
      }
   }
}

void MC::OutputMolecule(const ElectronicStructure& electronicStructure,
                        const Molecule& molecule,
                        int elecState) const{
   this->OutputEnergies(electronicStructure, elecState);
   this->OutputLog("\n");
   molecule.OutputConfiguration();
   molecule.OutputXyzCOC();
}

void MC::OutputEnergies(const MolDS_base::ElectronicStructure& electronicStructure,
                        int elecState) const{
   double eV2AU = Parameters::GetInstance()->GetEV2AU();
   this->OutputLog(this->messageEnergies);
   this->OutputLog(this->messageEnergiesTitle);
   this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageCoreRepulsionEnergy.c_str()
                                                     % electronicStructure.GetCoreRepulsionEnergy()
                                                     % (electronicStructure.GetCoreRepulsionEnergy()/eV2AU));
   if(Parameters::GetInstance()->RequiresVdWSCF()){
      this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageVdWCorrectionEnergy.c_str()
                                                        % electronicStructure.GetVdWCorrectionEnergy()
                                                        % (electronicStructure.GetVdWCorrectionEnergy()/eV2AU));
      this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageElectronicEnergyVdW.c_str()
                                                        % electronicStructure.GetElectronicEnergy(elecState)
                                                        % (electronicStructure.GetElectronicEnergy(elecState)/eV2AU));
   }
   else{
      this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageElectronicEnergy.c_str()
                                                        % electronicStructure.GetElectronicEnergy(elecState)
                                                        % (electronicStructure.GetElectronicEnergy(elecState)/eV2AU));
   }
}

}

