//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
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
#include"../base/factories/ElectronicStructureFactory.h"
#include"NASCO.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
using namespace MolDS_base_factories;

namespace MolDS_nasco{
NASCO::NASCO(){
   this->SetMessages();
   this->SetEnableTheoryTypes();
   //this->OutputLog("NASCO created \n");
}

NASCO::~NASCO(){
   //this->OutputLog("NASCO deleted\n");
}

void NASCO::DoNASCO(Molecule& molecule){
   this->OutputLog(this->messageStartNASCO);
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   this->CheckEnableTheoryType(theory);

   // create tmp molecule and electronic structure
   Molecule tmpMolecule(molecule);

   // malloc electornic structure 
   boost::shared_ptr<ElectronicStructure> electronicStructure1(ElectronicStructureFactory::Create());
   ElectronicStructure* currentES = electronicStructure1.get();
   currentES->SetMolecule(&molecule);
   currentES->SetCanOutputLogs(this->CanOutputLogs());
   molecule.SetCanOutputLogs(this->CanOutputLogs());

   // create temporary electronic structure
   boost::shared_ptr<ElectronicStructure> electronicStructure2(ElectronicStructureFactory::Create());
   ElectronicStructure* tmpES = electronicStructure2.get();
   tmpES->SetMolecule(&tmpMolecule);
   tmpES->SetCanOutputLogs(this->CanOutputLogs());
   tmpMolecule.SetCanOutputLogs(this->CanOutputLogs());

   // create real random generator
	boost::mt19937 realGenerator(Parameters::GetInstance()->GetSeedNASCO());
	boost::uniform_real<> range(0, 1);
	boost::variate_generator<boost::mt19937&, boost::uniform_real<> > realRand( realGenerator, range );

   const int totalSteps         = Parameters::GetInstance()->GetTotalStepsNASCO();
   const double dt              = Parameters::GetInstance()->GetTimeWidthNASCO();
   const int numElecStatesNASCO = Parameters::GetInstance()->GetNumberElectronicStatesNASCO();
   int elecState                = 0; //Parameters::GetInstance()->GetElectronicStateIndexNASCO();
   double time                  = 0.0;
   bool requireGuess            = false;
   double** matrixForce         = NULL;
   double initialEnergy         = 0.0;

   // initial calculation
   currentES->DoSCF();
   currentES->DoCIS();
   matrixForce = currentES->GetForce(elecState);

   // output initial conditions
   this->OutputLog(this->messageinitialConditionNASCO);
   initialEnergy = this->OutputEnergies(*currentES, molecule);
   this->OutputLog("\n");
   molecule.OutputConfiguration();
   molecule.OutputXyzCOM();
   molecule.OutputXyzCOC();
   molecule.OutputMomenta();

   // malloc ovelap AOs, MOs, ESs between differentMolecules
   double** overlapAOs = NULL;
   double** overlapMOs = NULL;
   double** overlapESs = NULL;

   try{
      this->MallocOverlapsDifferentMolecules(&overlapAOs, &overlapMOs, &overlapESs, molecule);
      for(int s=0; s<totalSteps; s++){
         this->OutputLog(boost::format("%s%d\n") % this->messageStartStepNASCO.c_str() % (s+1) );

         // update momenta
         this->UpdateMomenta(molecule, matrixForce, dt);

         // update coordinates
         for(int a=0; a<molecule.GetNumberAtoms(); a++){
            Atom* atom = molecule.GetAtom(a);
            Atom* tmpAtom  = tmpMolecule.GetAtom(a);
            double coreMass = tmpAtom->GetCoreMass();
            for(int i=0; i<CartesianType_end; i++){
               tmpAtom->GetXyz()[i] =  atom->GetXyz()[i] + dt*atom->GetPxyz()[i]/coreMass;
            }
         }
         tmpMolecule.CalcXyzCOM();
         tmpMolecule.CalcXyzCOC();

         // update electronic structure
         requireGuess = (s==0);
         tmpES->DoSCF(requireGuess);
         tmpES->DoCIS();

         // update force
         matrixForce = tmpES->GetForce(elecState);

         // update momenta
         this->UpdateMomenta(molecule, matrixForce, dt);

         // calculate overlaps
         currentES->CalcOverlapAOsWithAnotherConfiguration(overlapAOs, tmpMolecule);
         currentES->CalcOverlapMOsWithAnotherElectronicStructure(overlapMOs, overlapAOs, *tmpES);
         cout << "overlapMOs" << endl;
         for(int i=0; i<molecule.GetTotalNumberAOs(); i++){
            for(int j=0; j<molecule.GetTotalNumberAOs(); j++){
               printf("%e\t",overlapMOs[i][j]);
            }
            cout << endl;
         }
            
         // Synchronous molecular configuration and electronic states
         this->SynchronousMolecularConfiguration(molecule, tmpMolecule);
         swap(currentES, tmpES);
         currentES->SetMolecule(&molecule);
         tmpES->SetMolecule(&tmpMolecule);

         // output results
         this->OutputEnergies(*currentES, initialEnergy, molecule);
         molecule.OutputConfiguration();
         molecule.OutputXyzCOM();
         molecule.OutputXyzCOC();
         molecule.OutputMomenta();
         this->OutputLog(boost::format("%s%lf\n") % this->messageTime.c_str() 
                                                  % (dt*static_cast<double>(s+1)/Parameters::GetInstance()->GetFs2AU()));
         this->OutputLog(boost::format("%s%d\n") % this->messageEndStepNASCO.c_str() % (s+1) );
      }
   }
   catch(MolDSException ex){
      this->FreeOverlapsDifferentMolecules(&overlapAOs, &overlapMOs, &overlapESs, molecule);
      throw ex;
   }
   this->FreeOverlapsDifferentMolecules(&overlapAOs, &overlapMOs, &overlapESs, molecule);
   this->OutputLog(this->messageEndNASCO);
}

void NASCO::UpdateMomenta(Molecule& molecule, double const* const* matrixForce, double dt) const{
   for(int a=0; a<molecule.GetNumberAtoms(); a++){
      Atom* atom = molecule.GetAtom(a);
      for(int i=0; i<CartesianType_end; i++){
         atom->GetPxyz()[i] += 0.5*dt*(matrixForce[a][i]);
      }
   }
}

void NASCO::SynchronousMolecularConfiguration(Molecule& target, 
                                              const Molecule& refference) const{
   for(int a=0; a<target.GetNumberAtoms(); a++){
      Atom& targetAtom = *target.GetAtom(a);
      const Atom& refferenceAtom = *refference.GetAtom(a);
      for(int i=0; i<CartesianType_end; i++){
         targetAtom.GetXyz()[i] = refferenceAtom.GetXyz()[i];
      }
   }
   target.CalcXyzCOC();
   target.CalcXyzCOM();
}

void NASCO::SetMessages(){
   this->errorMessageTheoryType = "\ttheory type = ";
   this->errorMessageNotEnebleTheoryType  
      = "Error in nasco::NASCO::CheckEnableTheoryType: Non available theory is set.\n";
   this->messageStartNASCO = "**********  START: Molecular dynamics  **********\n";
   this->messageEndNASCO = "**********  DONE: Molecular dynamics  **********\n";
   this->messageinitialConditionNASCO = "\n\t========= Initial conditions \n";
   this->messageStartStepNASCO = "\n\t========== START: NASCO step ";
   this->messageEndStepNASCO =     "\t========== DONE: NASCO step ";
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

double NASCO::OutputEnergies(const ElectronicStructure& electronicStructure, const Molecule& molecule){
   int elecState = 0; //Parameters::GetInstance()->GetElectronicStateIndexNASCO();
   double eV2AU = Parameters::GetInstance()->GetEV2AU();
   double coreKineticEnergy = 0.0;
   for(int a=0; a<molecule.GetNumberAtoms(); a++){
      Atom* atom = molecule.GetAtom(a);
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
   this->OutputLog(boost::format("\t\t%s\t%e\t%e\n") % this->messageTotalEnergy.c_str()
                                                     % (coreKineticEnergy + electronicStructure.GetElectronicEnergy(elecState))
                                                     % ((coreKineticEnergy + electronicStructure.GetElectronicEnergy(elecState))/eV2AU));
   return (coreKineticEnergy + electronicStructure.GetElectronicEnergy(elecState));
}

void NASCO::OutputEnergies(const ElectronicStructure& electronicStructure, double initialEnergy, const Molecule& molecule){
   double energy = this->OutputEnergies(electronicStructure, molecule);
   double eV2AU = Parameters::GetInstance()->GetEV2AU();
   this->OutputLog(boost::format("\t\t%s\t%e\t%e\n\n") % this->messageErrorEnergy.c_str()
                                                       % (initialEnergy - energy)
                                                       % ((initialEnergy - energy)/eV2AU));
}

void NASCO::SetEnableTheoryTypes(){
   this->enableTheoryTypes.clear();
   this->enableTheoryTypes.push_back(MNDO);
   this->enableTheoryTypes.push_back(AM1);
   this->enableTheoryTypes.push_back(AM1D);
   this->enableTheoryTypes.push_back(PM3);
   this->enableTheoryTypes.push_back(PM3D);
   this->enableTheoryTypes.push_back(PM3PDDG);
}

void NASCO::CheckEnableTheoryType(TheoryType theoryType){

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

void NASCO::MallocOverlapsDifferentMolecules(double*** overlapAOs,
                                             double*** overlapMOs, 
                                             double*** overlapESs, 
                                             const Molecule& molecule) const{
   int dimOverlapAOs = molecule.GetTotalNumberAOs();
   int dimOverlapMOs = dimOverlapAOs;
   int dimOverlapESs = Parameters::GetInstance()->GetNumberElectronicStatesNASCO();
   MallocerFreer::GetInstance()->Malloc<double>(overlapAOs, dimOverlapAOs, dimOverlapAOs);
   MallocerFreer::GetInstance()->Malloc<double>(overlapMOs, dimOverlapMOs, dimOverlapMOs);
   MallocerFreer::GetInstance()->Malloc<double>(overlapESs, dimOverlapESs, dimOverlapESs);
}

void NASCO::FreeOverlapsDifferentMolecules(double*** overlapAOs,
                                           double*** overlapMOs, 
                                           double*** overlapESs, 
                                           const MolDS_base::Molecule& molecule) const{
   int dimOverlapAOs = molecule.GetTotalNumberAOs();
   int dimOverlapMOs = dimOverlapAOs;
   int dimOverlapESs = Parameters::GetInstance()->GetNumberElectronicStatesNASCO();
   MallocerFreer::GetInstance()->Free<double>(overlapAOs, dimOverlapAOs, dimOverlapAOs);
   MallocerFreer::GetInstance()->Free<double>(overlapMOs, dimOverlapMOs, dimOverlapMOs);
   MallocerFreer::GetInstance()->Free<double>(overlapESs, dimOverlapESs, dimOverlapESs);
}

}



