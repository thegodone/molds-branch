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
#include<time.h>
#include<list>
#include<vector>
#include<stdexcept>
#include<omp.h>
#include<boost/shared_ptr.hpp>
#include<boost/random.hpp>
#include"PrintController.h"
#include"MolDSException.h"
#include"Uncopyable.h"
#include"../wrappers/Lapack.h"
#include"Utilities.h"
#include"Enums.h"
#include"MallocerFreer.h"
#include"EularAngle.h"
#include"Parameters.h"
#include"atoms/Atom.h"
#include"factories/AtomFactory.h"
#include"Molecule.h"
#include"InputParser.h"
#include"GTOExpansionSTO.h"
#include"ElectronicStructure.h"
#include"factories/ElectronicStructureFactory.h"
#include"../md/MD.h"
#include"../mc/MC.h"
#include"../rpmd/RPMD.h"
#include"../optimization/Optimizer.h"
#include"factories/OptimizerFactory.h"
#include"MolDS.h"
using namespace std;
using namespace MolDS_base_factories;
namespace MolDS_base{
void MolDS::Run() const{
   // Welcome Messages
   this->OutputLog(Utilities::GetWelcomeMessage());
   
   //timer set
   time_t startTime;
   time(&startTime);
   clock_t startTick = clock();
   double ompStartTime = omp_get_wtime();

   Molecule* molecule = NULL;
   bool runningNormally = true;
   try{
      // declare 
      MallocerFreer::GetInstance();
      InputParser::GetInstance();
      molecule = new Molecule();
      Parameters::GetInstance();
      MolDS_wrappers::Lapack::GetInstance();
      GTOExpansionSTO::GetInstance();
      // Parse input
      InputParser::GetInstance()->Parse(molecule);
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      runningNormally = false;
   }
   
   // once electronic structure calculation
   if(runningNormally && Parameters::GetInstance()->GetCurrentSimulation() == Once){
      this->CalculateElectronicStructureOnce(molecule, &runningNormally);
   }

   // MD
   else if(runningNormally && Parameters::GetInstance()->GetCurrentSimulation() == MD){
      this->DoMD(molecule, &runningNormally);
   }

   // MC
   else if(runningNormally && Parameters::GetInstance()->GetCurrentSimulation() == MC){
      this->DoMC(molecule, &runningNormally);
   }

   // RPMD
   else if(runningNormally && Parameters::GetInstance()->GetCurrentSimulation() == RPMD){
      this->DoRPMD(molecule, &runningNormally);
   }

   // Optimization
   else if(runningNormally && Parameters::GetInstance()->GetCurrentSimulation() == Optimization){
      this->OptimizeGeometry(molecule, &runningNormally);
   }

   // Diagonalize Inertia Tensor
   else if(runningNormally && Parameters::GetInstance()->GetCurrentSimulation() == PrincipalAxes ){
      this->DiagonalizePrincipalAxes(molecule, &runningNormally);
   }

   // Translate molecule
   else if(runningNormally && Parameters::GetInstance()->GetCurrentSimulation() == Translate){
      this->TranslateMolecule(molecule, &runningNormally);
   }

   // Rotate molecule
   else if(runningNormally && Parameters::GetInstance()->GetCurrentSimulation() == Rotate){
      this->RotateMolecule(molecule, &runningNormally);
   }

   //Free 
   GTOExpansionSTO::DeleteInstance();
   MolDS_wrappers::Lapack::DeleteInstance(); 
   Parameters::DeleteInstance();
   delete molecule;
   InputParser::DeleteInstance();
   MallocerFreer::DeleteInstance();

   // Farewell Messages
   this->OutputLog(Utilities::GetFarewellMessage(startTime, startTick, ompStartTime, runningNormally));
}

void MolDS::CalculateElectronicStructureOnce(Molecule* molecule, bool* runningNormally) const{
   try{
      boost::shared_ptr<ElectronicStructure> electronicStructure(ElectronicStructureFactory::Create());
      electronicStructure->SetMolecule(molecule);
      electronicStructure->DoSCF();
      if(Parameters::GetInstance()->RequiresCIS()){
         electronicStructure->DoCIS();
      }
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      *runningNormally = false;
   }
}

void MolDS::DoMC(Molecule* molecule, bool* runningNormally) const{
   try{
      boost::shared_ptr<MolDS_mc::MC> mc(new MolDS_mc::MC());
      mc->SetMolecule(molecule);
      mc->DoMC();
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      *runningNormally = false;
   }
}

void MolDS::DoMD(Molecule* molecule, bool* runningNormally) const{
   try{
      boost::shared_ptr<MolDS_md::MD> md(new MolDS_md::MD());
      md->SetMolecule(molecule);
      md->DoMD();
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      *runningNormally = false;
   }
}

void MolDS::DoRPMD(Molecule* molecule, bool* runningNormally) const{
   try{
      boost::shared_ptr<MolDS_rpmd::RPMD> rpmd(new MolDS_rpmd::RPMD());
      rpmd->DoRPMD(*molecule);
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      *runningNormally = false;
   }
}

void MolDS::OptimizeGeometry(Molecule* molecule, bool* runningNormally) const{
   try{
      boost::shared_ptr<MolDS_optimization::Optimizer> optimizer(OptimizerFactory::Create());
      optimizer->Optimize(*molecule);
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      *runningNormally = false;
   }
}

void MolDS::DiagonalizePrincipalAxes(Molecule* molecule, bool* runningNormally) const{
   try{
      molecule->CalcPrincipalAxes();
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      *runningNormally = false;
   }
}

void MolDS::TranslateMolecule(Molecule* molecule, bool* runningNormally) const{
   try{
      molecule->Translate();
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      *runningNormally = false;
   }
}

void MolDS::RotateMolecule(Molecule* molecule, bool* runningNormally) const{
   try{
      molecule->Rotate();
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      *runningNormally = false;
   }
}

}

