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
#include<complex>
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
#include"Ehrenfest.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
using namespace MolDS_base_factories;

namespace MolDS_ehrenfest{
Ehrenfest::Ehrenfest(){
   this->molecule = NULL;
   this->SetMessages();
   this->SetEnableTheoryTypes();
   this->elecStates = NULL;
   this->superpositionCoeff=NULL;
   //this->OutputLog("Ehrenfest created \n");
}

Ehrenfest::~Ehrenfest(){
   this->elecStates->clear();
   if(this->elecStates != NULL){
      delete this->elecStates;
   }
   if(this->superpositionCoeff != NULL){
      int highestElecState = Parameters::GetInstance()->GetHighestElectronicStateIndexEhrenfest();
      MallocerFreer::GetInstance()->Free< std::complex<double> >(&this->superpositionCoeff, highestElecState+1);
   }
   //this->OutputLog("Ehrenfest deleted\n");
}

void Ehrenfest::SetMolecule(Molecule* molecule){
   // check enable electonic theory
   this->molecule = molecule;
   this->elecStates = new vector<int>;
   int highestElecState = Parameters::GetInstance()->GetHighestElectronicStateIndexEhrenfest();
   MallocerFreer::GetInstance()->Malloc< std::complex<double> >(&this->superpositionCoeff, highestElecState+1);
}

void Ehrenfest::DoEhrenfest(){
   this->OutputLog(this->messageStartEhrenfest);

   // validate
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   this->CheckEnableTheoryType(theory);

   // declare variables and malloc arrays
   int                         atomNum               = this->molecule->GetAtomVect().size();
   int                         totalSteps            = Parameters::GetInstance()->GetTotalStepsEhrenfest();
   int                         iniElecState          = Parameters::GetInstance()->GetInitialElectronicStateIndexEhrenfest();
   int                         highestElecState      = Parameters::GetInstance()->GetHighestElectronicStateIndexEhrenfest();
   int                         lowestElecState       = Parameters::GetInstance()->GetLowestElectronicStateIndexEhrenfest();
   double                      dt                    = Parameters::GetInstance()->GetTimeWidthEhrenfest();
   int                         dtRatio               = 1;
   double                      dtElec                = dt/static_cast<double>(dtRatio);
   double                      time                  = 0.0;
   bool                        requireGuess          = false;
   double                      initialEnergy         = 0.0;
   double                      elecNorm              = 0.0;
   double const* const* const* matrixForce           = NULL;
   double**                    matrixMeanForce       = NULL;
   double**                    overlapAOs            = NULL;
   double**                    overlapMOs            = NULL;
   double**                    overlapSDs            = NULL;
   double**                    overlapESs            = NULL;
   std::complex<double>*       tmpSuperpositionCoeff = NULL;
   try{
      this->MallocOverlapsDifferentMolecules(&overlapAOs, &overlapMOs, &overlapSDs, &overlapESs, &matrixMeanForce, &tmpSuperpositionCoeff, *this->molecule);
      for(int i=lowestElecState; i<=highestElecState; i++){
         this->elecStates->push_back(i);
      }

      // prepare trial molecule and electronic structure pointa
      Molecule trialMolecule(*this->molecule);
      boost::shared_ptr<ElectronicStructure> electronicStructure2(ElectronicStructureFactory::Create());
      ElectronicStructure* trialES = electronicStructure2.get();
      trialES->SetMolecule(&trialMolecule);
      trialES->SetCanOutputLogs(this->CanOutputLogs());
      trialMolecule.SetCanOutputLogs(this->CanOutputLogs());

      // malloc electornic structure
      boost::shared_ptr<ElectronicStructure> electronicStructure1(ElectronicStructureFactory::Create());
      ElectronicStructure* currentES = electronicStructure1.get();
      currentES->SetMolecule(this->molecule);
      currentES->SetCanOutputLogs(this->CanOutputLogs());
      this->molecule->SetCanOutputLogs(this->CanOutputLogs());

      // initial calculation
      currentES->DoSCF();
      if(Parameters::GetInstance()->RequiresCIS()){
         currentES->DoCIS();
      }
      matrixForce = electronicStructure1->GetForce(*this->elecStates);
      this->superpositionCoeff[iniElecState] = std::complex<double>(1.0,0.0);
      elecNorm = this->GetElecNorm();
      this->CalcMeanForce(matrixMeanForce, matrixForce, elecNorm);

      // output initial conditions
      this->OutputLog(this->messageinitialConditionEhrenfest);
      initialEnergy = this->OutputEnergies(*currentES, elecNorm);
      this->OutputLog("\n");
      this->molecule->OutputConfiguration();
      this->molecule->OutputXyzCOM();
      this->molecule->OutputXyzCOC();
      this->molecule->OutputMomenta();

      for(int s=0; s<totalSteps; s++){
         this->OutputLog(boost::format("%s%d\n") % this->messageStartStepEhrenfest.c_str() % (s+1) );

         // update momenta & coordinates
         this->UpdateMomenta    (trialMolecule, matrixMeanForce, dt);
         this->UpdateCoordinates(trialMolecule, dt);

         // calculate trial electronic structure
         requireGuess = (s==0) ? true : false;
         trialES->DoSCF(requireGuess);
         if(Parameters::GetInstance()->RequiresCIS()){
            trialES->DoCIS();
         }

         // update force
         matrixForce = trialES->GetForce(*this->elecStates);
         this->CalcMeanForce(matrixMeanForce, matrixForce, elecNorm);

         // update momenta
         this->UpdateMomenta(trialMolecule, matrixMeanForce, dt);

         // update electronic superpositioncoeff
         currentES->CalcOverlapAOsWithAnotherConfiguration             (overlapAOs, trialMolecule);
         currentES->CalcOverlapMOsWithAnotherElectronicStructure       (overlapMOs, overlapAOs, *trialES);
         currentES->CalcOverlapSingletSDsWithAnotherElectronicStructure(overlapSDs, overlapMOs);
         currentES->CalcOverlapESsWithAnotherElectronicStructure       (overlapESs, overlapSDs, *trialES);
         for(int es=0; es<dtRatio; es++){
            for(int lState=lowestElecState; lState<=highestElecState; lState++){
               double phase = 0.0;
               tmpSuperpositionCoeff[lState] = std::complex<double>(0.0,0.0);
               for(int rState=lowestElecState; rState<=highestElecState; rState++){
                  phase = -1.0*currentES->GetElectronicEnergy(rState)*dtElec;
                  std::complex<double> timePropagator = std::polar(1.0, phase);
                  tmpSuperpositionCoeff[lState] += overlapESs[lState][rState]
                                                  *timePropagator
                                                  *this->superpositionCoeff[rState];
               }
            }
            for(int state=lowestElecState; state<=highestElecState; state++){
               this->superpositionCoeff[state] = tmpSuperpositionCoeff[state];
            }
         }
         elecNorm = this->GetElecNorm();

         // swap electronicstructures
         this->molecule->SynchronizePhaseSpacePointTo(trialMolecule);
         swap(currentES, trialES);
         currentES->SetMolecule(this->molecule);
         trialES->SetMolecule(&trialMolecule);

         // Broadcast to all processes
         int root = MolDS_mpi::MpiProcess::GetInstance()->GetHeadRank();
         this->molecule->BroadcastPhaseSpacePointToAllProcesses(root);
         MolDS_mpi::MpiProcess::GetInstance()->Broadcast(this->superpositionCoeff, highestElecState+1, root);

         // output results
         this->OutputEnergies(*currentES, initialEnergy, elecNorm);
         this->molecule->OutputConfiguration();
         this->molecule->OutputXyzCOM();
         this->molecule->OutputXyzCOC();
         this->molecule->OutputMomenta();
         this->OutputLog(boost::format("%s%lf\n") % this->messageTime.c_str() 
                                                  % (dt*static_cast<double>(s+1)/Parameters::GetInstance()->GetFs2AU()));
         this->OutputLog(boost::format("%s%d\n") % this->messageEndStepEhrenfest.c_str() % (s+1) );
      }
   }
   catch(MolDSException ex){
      this->FreeOverlapsDifferentMolecules(&overlapAOs, &overlapMOs, &overlapSDs, &overlapESs, &matrixMeanForce, &tmpSuperpositionCoeff, *this->molecule);
      throw ex;
   }
   this->FreeOverlapsDifferentMolecules(&overlapAOs, &overlapMOs, &overlapSDs, &overlapESs, &matrixMeanForce, &tmpSuperpositionCoeff, *this->molecule);
   this->OutputLog(this->messageEndEhrenfest);
}

double Ehrenfest::GetElecNorm(){
   double elecNorm=0.0;
   int    atomNum = this->molecule->GetAtomVect().size();
   int    highestElecState = Parameters::GetInstance()->GetHighestElectronicStateIndexEhrenfest();
   int    lowestElecState  = Parameters::GetInstance()->GetLowestElectronicStateIndexEhrenfest();
   for(int state=lowestElecState; state<=highestElecState; state++){
      double popu = abs(this->superpositionCoeff[state]);
      elecNorm += popu;
   }
   return elecNorm;
}

void Ehrenfest::CalcMeanForce(double** matrixMeanForce, double const* const* const* matrixForce, double elecNorm){
   int atomNum = this->molecule->GetAtomVect().size();
   int highestElecState = Parameters::GetInstance()->GetHighestElectronicStateIndexEhrenfest();
   int lowestElecState  = Parameters::GetInstance()->GetLowestElectronicStateIndexEhrenfest();
   MallocerFreer::GetInstance()->Initialize<double>(matrixMeanForce, atomNum, CartesianType_end);
   for(int state=lowestElecState; state<=highestElecState; state++){
      double popu = abs(this->superpositionCoeff[state])/elecNorm;
      for(int a=0; a<atomNum; a++){
         for(int axis=0; axis<CartesianType_end; axis++){
            matrixMeanForce[a][axis] += popu*matrixForce[state][a][axis];
         }
      }
   }
}

void Ehrenfest::UpdateMomenta(const Molecule& molecule, double const* const* matrixForce, double dt) const{
#pragma omp parallel for schedule(dynamic, MOLDS_OMP_DYNAMIC_CHUNK_SIZE)
   for(int a=0; a<molecule.GetAtomVect().size(); a++){
      Atom* atom = molecule.GetAtomVect()[a];
      for(int i=0; i<CartesianType_end; i++){
         atom->GetPxyz()[i] += 0.5*dt*(matrixForce[a][i]);
      }
   }
}

void Ehrenfest::UpdateCoordinates(Molecule& molecule, double dt) const{
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

void Ehrenfest::SetMessages(){
   this->errorMessageTheoryType = "\ttheory type = ";
   this->errorMessageNotEnebleTheoryType  
      = "Error in ehrenfest::Ehrenfest::CheckEnableTheoryType: Non available theory is set.\n";
   this->messageStartEhrenfest = "**********  START: Ehrenfest dynamics  **********\n";
   this->messageEndEhrenfest = "**********  DONE: Ehrenfest dynamics  **********\n";
   this->messageinitialConditionEhrenfest = "\n\t========= Initial conditions \n";
   this->messageStartStepEhrenfest = "\n\t========== START: Ehrenfest step ";
   this->messageEndStepEhrenfest =     "\t========== DONE: Ehrenfest step ";
   this->messageEnergies            = "\t\tEnergies:\n\t\t(Note that following electronic energies include core repulsion energies)\n\n";
   this->messageEnergiesVdW         = "\t\tEnergies:\n\t\t(Note that following electronic energies include core repulsion and VdW energies)\n\n";
   this->messageElecPopuTitle       = "\t\t|  i-th eigenstate      | total (norm) |      0      |      1      | ...\n";
   this->messageElecPopulations     = "\t\tElectronic populations:  ";
   this->messageElecEigenEnergyTitle= "\t\t\t|       i-th eigenstate          |      0      |      1      | ...\n";
   this->messageElecEigenEnergyAU   = "\t\t\tElectronic eiginenergies [a.u.]:  ";
   this->messageElecEigenEnergyEV   = "\t\t\tElectronic eiginenergies [eV]:    ";
   this->messageEnergiesTitle       = "\t\t\t|  kind           | [a.u.] | [eV] | \n";
   this->messageCoreKineticEnergy   = "\t\t\tCore kinetic:     ";
   this->messageCoreRepulsionEnergy = "\t\t\tCore repulsion:   ";
   this->messageVdWCorrectionEnergy = "\t\t\tVdW correction:   ";
   this->messageElectronicEnergy    = "\t\t\tMean electronic:  ";
   this->messageTotalEnergy         = "\t\t\tTotal:            ";
   this->messageErrorEnergy         = "\t\t\tError:            ";
   this->messageTime = "\tTime in [fs]: ";
}

double Ehrenfest::OutputEnergies(const MolDS_base::ElectronicStructure& electronicStructure, double elecNorm){
   //int elecState = Parameters::GetInstance()->GetElectronicStateIndexEhrenfest();
   double eV2AU  = Parameters::GetInstance()->GetEV2AU();
   double coreKineticEnergy = 0.0;
   for(int a=0; a<this->molecule->GetAtomVect().size(); a++){
      Atom* atom = this->molecule->GetAtomVect()[a];
      double coreMass = atom->GetCoreMass();
      for(int i=0; i<CartesianType_end; i++){
         coreKineticEnergy += 0.5*pow(atom->GetPxyz()[i],2.0)/coreMass;
      }
   }  

   stringstream ssEigenEneInAU;
   stringstream ssEigenEneInEV;
   stringstream ssPopulations;
   ssEigenEneInAU << this->messageElecEigenEnergyTitle;
   ssEigenEneInAU << this->messageElecEigenEnergyAU;
   ssEigenEneInEV << this->messageElecEigenEnergyEV;
   ssPopulations  << this->messageElecPopuTitle;
   ssPopulations  << this->messageElecPopulations;
      ssPopulations  << (boost::format("%e\t") % elecNorm).str();

   double meanElecEnergy=0.0;
   int    highestElecState = Parameters::GetInstance()->GetHighestElectronicStateIndexEhrenfest();
   int    lowestElecState  = Parameters::GetInstance()->GetLowestElectronicStateIndexEhrenfest();
   for(int state=0; state<=highestElecState; state++){
      double popu = abs(this->superpositionCoeff[state]);
      double ene  = electronicStructure.GetElectronicEnergy(state);
      meanElecEnergy += popu*ene;
      ssEigenEneInAU << (boost::format("%e\t") % ene).str();
      ssEigenEneInEV << (boost::format("%e\t") % (ene/eV2AU)).str();
      ssPopulations  << (boost::format("%e\t") % abs(this->superpositionCoeff[state])).str();
   }
   meanElecEnergy /= elecNorm;
   ssEigenEneInAU << endl;
   ssEigenEneInEV << endl << endl;
   ssPopulations  << endl << endl;

   // output electronic population
   this->OutputLog(ssPopulations.str());

   // output energies:
   if(Parameters::GetInstance()->RequiresVdWSCF()){
      this->OutputLog(this->messageEnergiesVdW);
   }
   else{
      this->OutputLog(this->messageEnergies);
   }
   this->OutputLog(ssEigenEneInAU.str());
   this->OutputLog(ssEigenEneInEV.str());
   this->OutputLog(this->messageEnergiesTitle);
   this->OutputLog(boost::format("%s\t%e\t%e\n") % this->messageCoreKineticEnergy.c_str()
                                                 % coreKineticEnergy
                                                 % (coreKineticEnergy/eV2AU));
   this->OutputLog(boost::format("%s\t%e\t%e\n") % this->messageCoreRepulsionEnergy.c_str()
                                                 % electronicStructure.GetCoreRepulsionEnergy()
                                                 % (electronicStructure.GetCoreRepulsionEnergy()/eV2AU));
   if(Parameters::GetInstance()->RequiresVdWSCF()){
      this->OutputLog(boost::format("%s\t%e\t%e\n") % this->messageVdWCorrectionEnergy.c_str()
                                                    % electronicStructure.GetVdWCorrectionEnergy()
                                                    % (electronicStructure.GetVdWCorrectionEnergy()/eV2AU));
   }
   this->OutputLog(boost::format("%s\t%e\t%e\n") % this->messageElectronicEnergy.c_str()
                                                 % meanElecEnergy 
                                                 % (meanElecEnergy/eV2AU));
   this->OutputLog(boost::format("%s\t%e\t%e\n") % this->messageTotalEnergy.c_str()
                                                 % ( coreKineticEnergy + meanElecEnergy)
                                                 % ((coreKineticEnergy + meanElecEnergy)/eV2AU));
   return (coreKineticEnergy + meanElecEnergy);
}

void Ehrenfest::OutputEnergies(const MolDS_base::ElectronicStructure& electronicStructure, double initialEnergy, double elecNorm){
   double energy = this->OutputEnergies(electronicStructure, elecNorm);
   double eV2AU = Parameters::GetInstance()->GetEV2AU();
   this->OutputLog(boost::format("%s\t%e\t%e\n\n") % this->messageErrorEnergy.c_str()
                                                   % (initialEnergy - energy)
                                                   % ((initialEnergy - energy)/eV2AU));
}

void Ehrenfest::SetEnableTheoryTypes(){
   this->enableTheoryTypes.clear();
   this->enableTheoryTypes.push_back(ZINDOS);
   this->enableTheoryTypes.push_back(MNDO);
   this->enableTheoryTypes.push_back(AM1);
   this->enableTheoryTypes.push_back(AM1D);
   this->enableTheoryTypes.push_back(PM3);
   this->enableTheoryTypes.push_back(PM3D);
   this->enableTheoryTypes.push_back(PM3PDDG);
}

void Ehrenfest::CheckEnableTheoryType(TheoryType theoryType){

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

void Ehrenfest::MallocOverlapsDifferentMolecules(double*** overlapAOs,
                                                 double*** overlapMOs, 
                                                 double*** overlapSDs, 
                                                 double*** overlapESs, 
                                                 double*** matrixMeanForce,
                                                 complex<double>** tmpSuperpositionCoeff,
                                                 const Molecule& molecule) const{
   int dimOverlapAOs    = molecule.GetTotalNumberAOs();
   int dimOverlapMOs    = dimOverlapAOs;
   int dimOverlapSDs    = Parameters::GetInstance()->GetActiveOccCIS()
                          *Parameters::GetInstance()->GetActiveVirCIS()
                          +1;
   int dimOverlapESs    = dimOverlapSDs;
   int atomNum          = this->molecule->GetAtomVect().size();
   int highestElecState = Parameters::GetInstance()->GetHighestElectronicStateIndexEhrenfest();
   MallocerFreer::GetInstance()->Malloc<double>(overlapAOs, dimOverlapAOs, dimOverlapAOs);
   MallocerFreer::GetInstance()->Malloc<double>(overlapMOs, dimOverlapMOs, dimOverlapMOs);
   MallocerFreer::GetInstance()->Malloc<double>(overlapSDs, dimOverlapSDs, dimOverlapSDs);
   MallocerFreer::GetInstance()->Malloc<double>(overlapESs, dimOverlapESs, dimOverlapESs);
   MallocerFreer::GetInstance()->Malloc<double>(matrixMeanForce, atomNum, CartesianType_end);
   MallocerFreer::GetInstance()->Malloc< std::complex<double> >(tmpSuperpositionCoeff, highestElecState+1);
}

void Ehrenfest::FreeOverlapsDifferentMolecules(double*** overlapAOs,
                                               double*** overlapMOs, 
                                               double*** overlapSDs, 
                                               double*** overlapESs, 
                                               double*** matrixMeanForce,
                                               complex<double>** tmpSuperpositionCoeff,
                                               const MolDS_base::Molecule& molecule) const{
   int dimOverlapAOs    = molecule.GetTotalNumberAOs();
   int dimOverlapMOs    = dimOverlapAOs;
   int dimOverlapSDs    = Parameters::GetInstance()->GetActiveOccCIS()
                          *Parameters::GetInstance()->GetActiveVirCIS()
                          +1;
   int dimOverlapESs    = dimOverlapSDs;
   int atomNum          = this->molecule->GetAtomVect().size();
   int highestElecState = Parameters::GetInstance()->GetHighestElectronicStateIndexEhrenfest();
   MallocerFreer::GetInstance()->Free<double>(overlapAOs, dimOverlapAOs, dimOverlapAOs);
   MallocerFreer::GetInstance()->Free<double>(overlapMOs, dimOverlapMOs, dimOverlapMOs);
   MallocerFreer::GetInstance()->Free<double>(overlapSDs, dimOverlapSDs, dimOverlapSDs);
   MallocerFreer::GetInstance()->Free<double>(overlapESs, dimOverlapESs, dimOverlapESs);
   MallocerFreer::GetInstance()->Free<double>(matrixMeanForce, atomNum, CartesianType_end);
   MallocerFreer::GetInstance()->Free< std::complex<double> >(tmpSuperpositionCoeff, highestElecState+1);
}

}



