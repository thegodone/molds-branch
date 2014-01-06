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
#include"../mc/MC.h"
#include"RPMD.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
using namespace MolDS_base_factories;

namespace MolDS_rpmd{
RPMD::RPMD(){
   this->SetMessages();
   this->SetEnableTheoryTypes();
   //this->OutputLog("RPMD created \n");
}

RPMD::~RPMD(){
   //this->OutputLog("RPMD deleted\n");
}

void RPMD::CreateBeads(vector<boost::shared_ptr<Molecule> >& molecularBeads,
                       vector<boost::shared_ptr<ElectronicStructure> >& electronicStructureBeads,
                       const Molecule& refferenceMolecule,
                       int numBeads){
   for(int b=0; b<numBeads; b++){
      // create molecular beads
      boost::shared_ptr<Molecule> molecule(new Molecule());
      *molecule = refferenceMolecule;
      molecularBeads.push_back(molecule);
      // create electronic structure beads
      boost::shared_ptr<ElectronicStructure> electronicStructure(ElectronicStructureFactory::Create());
      electronicStructure->SetMolecule(molecule.get());
      electronicStructureBeads.push_back(electronicStructure);
   }
}

void RPMD::UpdateElectronicStructure(const std::vector<boost::shared_ptr<ElectronicStructure> >& electronicStructureBeads){
   int numBeads = electronicStructureBeads.size();
   for(int b=0; b<numBeads; b++){
      this->OutputLog(boost::format("%s%d\n") % this->messageBeadsNum.c_str() % b);
      electronicStructureBeads[b]->DoSCF();
      if(Parameters::GetInstance()->RequiresCIS()){
         electronicStructureBeads[b]->DoCIS();
      }
   }
}

// elecState=0 means ground state.
void RPMD::UpdateMomenta(const vector<boost::shared_ptr<Molecule> >& molecularBeads,
                         const std::vector<boost::shared_ptr<ElectronicStructure> >& electronicStructureBeads,
                         int elecState,
                         double dt,
                         double temperature){
   double kB = Parameters::GetInstance()->GetBoltzmann();
   int numBeads = molecularBeads.size();
   int numAtom  = molecularBeads[0]->GetAtomVect().size();
   for(int b=0; b<numBeads; b++){
      int preB  = b==0 ? numBeads-1 : b-1;
      int postB = b==numBeads-1 ? 0 : b+1;
      double const* const* electronicForceMatrix 
         = electronicStructureBeads[b]->GetForce(elecState);;
      for(int a=0; a<numAtom; a++){
         Atom* atom      = molecularBeads[b]->GetAtomVect()[a];
         Atom* preAtom   = molecularBeads[preB]->GetAtomVect()[a];
         Atom* postAtom  = molecularBeads[postB]->GetAtomVect()[a];
         double coreMass = atom->GetCoreMass();
         for(int i=0; i<CartesianType_end; i++){
            double beadsForce = -1.0*coreMass*pow(kB*temperature*static_cast<double>(numBeads),2.0)
                               *(2.0*atom->GetXyz()[i] - preAtom->GetXyz()[i] - postAtom->GetXyz()[i]);
            double force = beadsForce + electronicForceMatrix[a][i];
            atom->GetPxyz()[i] += 0.5*dt*(force);
         }
      }
   }
}

void RPMD::UpdateCoordinates(const vector<boost::shared_ptr<Molecule> >& molecularBeads,
                             double dt){
   int numBeads = molecularBeads.size();
   int numAtom = molecularBeads[0]->GetAtomVect().size();
   for(int b=0; b<numBeads; b++){
      for(int a=0; a<numAtom; a++){
         Atom* atom = molecularBeads[b]->GetAtomVect()[a];
         double coreMass = atom->GetCoreMass();
         for(int i=0; i<CartesianType_end; i++){
            atom->GetXyz()[i] += dt*atom->GetPxyz()[i]/coreMass;
         }
      }
      molecularBeads[b]->CalcBasicsConfiguration();
   }
}

void RPMD::BroadcastPhaseSpacepointsToAllProcesses(std::vector<boost::shared_ptr<MolDS_base::Molecule> >& molecularBeads, int root) const{
   int numBeads = molecularBeads.size();
   for(int b=0; b<numBeads; b++){
      molecularBeads[b]->BroadcastPhaseSpacePointToAllProcesses(root);
   }
}

void RPMD::FluctuateBeads(const vector<boost::shared_ptr<Molecule> >& molecularBeads,
                          int elecState,
                          double temperature,
                          unsigned long seed){
   int numBeads = molecularBeads.size();
   double stepWidth = 0.01;
   for(int b=0; b<numBeads; b++){
      boost::shared_ptr<MolDS_mc::MC> mc(new MolDS_mc::MC());
      Molecule* molecule = molecularBeads[b].get();
      molecule->SetCanOutputLogs(false);
      mc->SetMolecule(molecule);
      mc->SetCanOutputLogs(false);
      mc->DoMC(molecule->GetAtomVect().size(), elecState, temperature, stepWidth, seed+b);
   }
}

void RPMD::DoRPMD(const Molecule& refferenceMolecule){
   this->OutputLog(this->messageStartRPMD);

   // validate theory
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   int elecState = Parameters::GetInstance()->GetElectronicStateIndexRPMD();
   this->CheckEnableTheoryType(theory, elecState);

   double temperature = Parameters::GetInstance()->GetTemperatureRPMD();
   unsigned long seed = Parameters::GetInstance()->GetSeedRPMD();
   int totalSteps     = Parameters::GetInstance()->GetTotalStepsRPMD();
   double dt          = Parameters::GetInstance()->GetTimeWidthRPMD();
   double kB          = Parameters::GetInstance()->GetBoltzmann();
   int numBeads       = Parameters::GetInstance()->GetNumberBeadsRPMD();
   int numAtom        = refferenceMolecule.GetAtomVect().size();

   // create Beads
   vector<boost::shared_ptr<Molecule> > molecularBeads;
   vector<boost::shared_ptr<ElectronicStructure> > electronicStructureBeads;
   this->CreateBeads(molecularBeads, electronicStructureBeads, refferenceMolecule, numBeads);

   // initialize Beads fluctuations (quantum fluctuations)
   this->FluctuateBeads(molecularBeads, elecState, temperature, seed);

   // initialize electronic states of each bead.
   this->OutputLog(this->messageStartInitialRPMD);
   this->UpdateElectronicStructure(electronicStructureBeads);
   this->OutputLog(this->messageinitialConditionRPMD);
   double initialEnergy = this->OutputEnergies(molecularBeads,
                                               electronicStructureBeads,
                                               elecState,
                                               temperature);
   this->OutputLog(this->messageEndInitialRPMD);

   // time step loop
   for(int s=0; s<totalSteps; s++){
      this->OutputLog(boost::format("%s%d\n") % this->messageStartStepRPMD.c_str() % (s+1));
      // update momenta
      this->UpdateMomenta(molecularBeads, electronicStructureBeads, elecState, dt, temperature);

      // update coordinates
      this->UpdateCoordinates(molecularBeads, dt);
      
      // update electronic structure
      this->UpdateElectronicStructure(electronicStructureBeads);

      // update momenta
      this->UpdateMomenta(molecularBeads, electronicStructureBeads, elecState, dt, temperature);

      // Broadcast coordinates and momenta of beads to all processes
      int root = MolDS_mpi::MpiProcess::GetInstance()->GetHeadRank();
      this->BroadcastPhaseSpacepointsToAllProcesses(molecularBeads, root);

      // output energy
      this->OutputEnergies(molecularBeads, 
                           electronicStructureBeads,
                           elecState,
                           temperature,
                           initialEnergy);

      this->OutputLog(boost::format("%s%d\n") % this->messageEndStepRPMD.c_str() % (s+1));
   }
      this->OutputLog(this->messageEndRPMD);
}

void RPMD::SetMessages(){
   this->errorMessageTheoryType = "\ttheory type = ";
   this->errorMessageElecState = "\tElectronic state = ";
   this->errorMessageNotEnebleTheoryType  
      = "Error in rpmd::RPMD::CheckEnableTheoryType: Non available theory is set.\n";
   this->messageStartRPMD 
      = "**********  START: Ring Polymer Molecular dynamics  **********\n";
   this->messageEndRPMD 
      = "**********  DONE: Ring Polymer Molecular dynamics  **********\n";
   this->messageStartInitialRPMD = "\n**********  START: Initial calculation of electronic structure of each bead  ********* \n";
   this->messageEndInitialRPMD =   "\n**********  DONE: Initial calculation of electronic structure of each bead   ********* \n";
   this->messageinitialConditionRPMD = "\n=========  Initial conditions of the beads  ==========\n";
   this->messageStartStepRPMD =    "\n==========  START: RPMD step ";
   this->messageEndStepRPMD =        "\n==========  DONE: RPMD step ";
   this->messageBeadsNum = "----------  Beads number ";
   this->messageEnergies = "\tEnergies:\n";
   this->messageEnergiesTitle = "\t\t|\tkind\t\t\t| [a.u.] | [eV] | \n";
   this->messageBeadsKineticEnergy =   "Beads kinetic     ";
   this->messageBeadsHarmonicEnergy = "Beads harmonic   ";
   this->messageElecStateEnergy = "Electronic\n\t\t(inc. core rep.)";
   this->messageTotalEnergy =         "Total            ";
   this->messageErrorEnergy =         "Error            ";
   this->messageTime = "\tTime in [fs]: ";
}

double RPMD::OutputEnergies(const vector<boost::shared_ptr<Molecule> >& molecularBeads,
                            const std::vector<boost::shared_ptr<ElectronicStructure> >& electronicStructureBeads,
                            int elecState,
                            double temperature){
   int numBeads = molecularBeads.size();
   int numAtom = molecularBeads[0]->GetAtomVect().size();
   double beadsKineticEnergy = 0.0;;
   for(int b=0; b<numBeads; b++){
      double coreKineticEnergy = 0.0;
      for(int a=0; a<numAtom; a++){
         Atom* atom = molecularBeads[b]->GetAtomVect()[a];
         double coreMass = atom->GetCoreMass();
         for(int i=0; i<CartesianType_end; i++){
            coreKineticEnergy += 0.5*pow(atom->GetPxyz()[i],2.0)/coreMass;
         }
      }
      beadsKineticEnergy += coreKineticEnergy;
   }  

   double kB = Parameters::GetInstance()->GetBoltzmann();
   double beadsHarmonicEnergy = 0.0;
   for(int b=0; b<numBeads; b++){
      double harmonicEnergy = 0.0;
      int preB = b==0 ? numBeads-1 : b-1;
      for(int a=0; a<numAtom; a++){
         Atom* atom    = molecularBeads[b]->GetAtomVect()[a];
         Atom* preAtom = molecularBeads[preB]->GetAtomVect()[a];
         double coreMass = atom->GetCoreMass();
         double dx = atom->GetXyz()[XAxis] - preAtom->GetXyz()[XAxis];
         double dy = atom->GetXyz()[YAxis] - preAtom->GetXyz()[YAxis];
         double dz = atom->GetXyz()[ZAxis] - preAtom->GetXyz()[ZAxis];
         harmonicEnergy += coreMass*(pow(dx,2.0)+pow(dy,2.0)+pow(dz,2.0));
      }
      harmonicEnergy *= 0.5*pow(kB*temperature*static_cast<double>(numBeads),2.0);
      beadsHarmonicEnergy += harmonicEnergy;
   }

   double elecStateEnergy = 0.0;
   for(int b=0; b<numBeads; b++){
      elecStateEnergy += electronicStructureBeads[b]->GetElectronicEnergy(elecState);
   }

   double totalEnergy = beadsKineticEnergy + beadsHarmonicEnergy + elecStateEnergy;
   double eV2AU = Parameters::GetInstance()->GetEV2AU();
   // output energies:
   this->OutputLog(this->messageEnergies);
   this->OutputLog(this->messageEnergiesTitle);
   this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageBeadsKineticEnergy.c_str()
                                                     % beadsKineticEnergy
                                                     % (beadsKineticEnergy/eV2AU));
   this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageBeadsHarmonicEnergy.c_str()
                                                     % beadsHarmonicEnergy
                                                     % (beadsHarmonicEnergy/eV2AU));
   this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageElecStateEnergy.c_str()
                                                     % elecStateEnergy
                                                     % (elecStateEnergy/eV2AU));
   this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageTotalEnergy.c_str()
                                                     % totalEnergy
                                                     % (totalEnergy/eV2AU));
   return totalEnergy;
}

void RPMD::OutputEnergies(const vector<boost::shared_ptr<Molecule> >& molecularBeads,
                          const std::vector<boost::shared_ptr<ElectronicStructure> >& electronicStructureBeads,
                          int elecState,
                          double temperature,
                          double initialEnergy){
   double energy = this->OutputEnergies(molecularBeads, electronicStructureBeads, elecState, temperature);
   this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageErrorEnergy.c_str()
                                                     % (initialEnergy - energy)
                                                     % ((initialEnergy - energy)/Parameters::GetInstance()->GetEV2AU()));
}

void RPMD::SetEnableTheoryTypes(){
   // ground state
   this->enableGroundStateTheoryTypes.clear();
   this->enableGroundStateTheoryTypes.push_back(ZINDOS);
   this->enableGroundStateTheoryTypes.push_back(MNDO);
   this->enableGroundStateTheoryTypes.push_back(AM1);
   this->enableGroundStateTheoryTypes.push_back(AM1D);
   this->enableGroundStateTheoryTypes.push_back(PM3);
   this->enableGroundStateTheoryTypes.push_back(PM3D);
   this->enableGroundStateTheoryTypes.push_back(PM3PDDG);

   // excited state
   this->enableExcitedStateTheoryTypes.clear();
   this->enableExcitedStateTheoryTypes.push_back(ZINDOS);
   this->enableExcitedStateTheoryTypes.push_back(MNDO);
   this->enableExcitedStateTheoryTypes.push_back(AM1);
   this->enableExcitedStateTheoryTypes.push_back(AM1D);
   this->enableExcitedStateTheoryTypes.push_back(PM3);
   this->enableExcitedStateTheoryTypes.push_back(PM3D);
   this->enableExcitedStateTheoryTypes.push_back(PM3PDDG);
}

void RPMD::CheckEnableTheoryType(TheoryType theoryType, int elecState){

   bool isEnable = false;
   int groundState = 0;
   if(elecState == groundState){
      for(int i=0; i<this->enableGroundStateTheoryTypes.size();i++){
         if(theoryType == this->enableGroundStateTheoryTypes[i]){
            isEnable = true;
            break;
         }
      }
   }
   else{
      for(int i=0; i<this->enableExcitedStateTheoryTypes.size();i++){
         if(theoryType == this->enableExcitedStateTheoryTypes[i]){
            isEnable = true;
            break;
         }
      }
   }
   if(!isEnable){
      stringstream ss;
      ss << this->errorMessageNotEnebleTheoryType;
      ss << this->errorMessageTheoryType << TheoryTypeStr(theoryType) << endl;
      ss << this->errorMessageElecState << elecState << endl;
      throw MolDSException(ss.str());
   }
}

}



