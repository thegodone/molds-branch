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
#include"../wrappers/Blas.h"
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
#include"../nasco/NASCO.h"
#include"../optimization/Optimizer.h"
#include"factories/OptimizerFactory.h"
#include"MolDS.h"
using namespace std;
using namespace MolDS_base_factories;
namespace MolDS_base{
void MolDS::Run(int argc, char *argv[]){
   bool runsNormally(true);
   Molecule* molecule(NULL);
   
   // timer and initialize
   try{
      this->Initialize();
      molecule = new Molecule();
      InputParser::GetInstance()->Parse(molecule, argc, argv);
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      ex.PrintBacktrace();
      runsNormally = false;
   }

   // once electronic structure calculation
   if(runsNormally && Parameters::GetInstance()->GetCurrentSimulation() == Once){
      this->CalculateElectronicStructureOnce(molecule, &runsNormally);
   }

   // MD
   else if(runsNormally && Parameters::GetInstance()->GetCurrentSimulation() == MD){
      this->DoMD(molecule, &runsNormally);
   }

   // MC
   else if(runsNormally && Parameters::GetInstance()->GetCurrentSimulation() == MC){
      this->DoMC(molecule, &runsNormally);
   }

   // RPMD
   else if(runsNormally && Parameters::GetInstance()->GetCurrentSimulation() == RPMD){
      this->DoRPMD(molecule, &runsNormally);
   }

   // NASCO
   else if(runsNormally && Parameters::GetInstance()->GetCurrentSimulation() == NASCO){
      this->DoNASCO(molecule, &runsNormally);
   }

   // Optimization
   else if(runsNormally && Parameters::GetInstance()->GetCurrentSimulation() == Optimization){
      this->OptimizeGeometry(molecule, &runsNormally);
   }

   // Diagonalize Inertia Tensor
   else if(runsNormally && Parameters::GetInstance()->GetCurrentSimulation() == PrincipalAxes ){
      this->DiagonalizePrincipalAxes(molecule, &runsNormally);
   }

   // Translate molecule
   else if(runsNormally && Parameters::GetInstance()->GetCurrentSimulation() == Translate){
      this->TranslateMolecule(molecule, &runsNormally);
   }

   // Rotate molecule
   else if(runsNormally && Parameters::GetInstance()->GetCurrentSimulation() == Rotate){
      this->RotateMolecule(molecule, &runsNormally);
   }

   delete molecule;
   this->Finalize(runsNormally);
}

void MolDS::Initialize(){
   // Welcome Messages
   this->OutputLog(Utilities::GetWelcomeMessage());
   
   //timer set
   time(&this->startTime);
   this->startTick = clock();
   this->ompStartTime = omp_get_wtime();

   // declare 
   MallocerFreer::GetInstance();
   InputParser::GetInstance();
   Parameters::GetInstance();
   MolDS_wrappers::Blas::GetInstance();
   MolDS_wrappers::Lapack::GetInstance();
   GTOExpansionSTO::GetInstance();
}

void MolDS::Finalize(bool runsNormally) const{
   //Free 
   GTOExpansionSTO::DeleteInstance();
   MolDS_wrappers::Lapack::DeleteInstance(); 
   MolDS_wrappers::Blas::DeleteInstance(); 
   Parameters::DeleteInstance();
   InputParser::DeleteInstance();
   MallocerFreer::DeleteInstance();

   // Farewell Messages
   this->OutputLog(Utilities::GetFarewellMessage(this->startTime, this->startTick, this->ompStartTime, runsNormally));
}

void MolDS::CalculateElectronicStructureOnce(Molecule* molecule, bool* runsNormally) const{
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
      ex.PrintBacktrace();
      *runsNormally = false;
   }
}

void MolDS::DoMC(Molecule* molecule, bool* runsNormally) const{
   try{
      boost::shared_ptr<MolDS_mc::MC> mc(new MolDS_mc::MC());
      mc->SetMolecule(molecule);
      mc->DoMC();
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      ex.PrintBacktrace();
      *runsNormally = false;
   }
}

void MolDS::DoMD(Molecule* molecule, bool* runsNormally) const{
   try{
      boost::shared_ptr<MolDS_md::MD> md(new MolDS_md::MD());
      md->SetMolecule(molecule);
      md->DoMD();
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      ex.PrintBacktrace();
      *runsNormally = false;
   }
}

void MolDS::DoRPMD(Molecule* molecule, bool* runsNormally) const{
   try{
      boost::shared_ptr<MolDS_rpmd::RPMD> rpmd(new MolDS_rpmd::RPMD());
      rpmd->DoRPMD(*molecule);
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      ex.PrintBacktrace();
      *runsNormally = false;
   }
}

void MolDS::DoNASCO(Molecule* molecule, bool* runsNormally) const{
   try{
      boost::shared_ptr<MolDS_nasco::NASCO> nasco(new MolDS_nasco::NASCO());
      nasco->DoNASCO(*molecule);
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      ex.PrintBacktrace();
      *runsNormally = false;
   }
}

void MolDS::OptimizeGeometry(Molecule* molecule, bool* runsNormally) const{
   try{
      boost::shared_ptr<MolDS_optimization::Optimizer> optimizer(OptimizerFactory::Create());
      optimizer->Optimize(*molecule);
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      ex.PrintBacktrace();
      *runsNormally = false;
   }
}

void MolDS::DiagonalizePrincipalAxes(Molecule* molecule, bool* runsNormally) const{
   try{
      molecule->CalcPrincipalAxes();
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      ex.PrintBacktrace();
      *runsNormally = false;
   }
}

void MolDS::TranslateMolecule(Molecule* molecule, bool* runsNormally) const{
   try{
      molecule->Translate();
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      ex.PrintBacktrace();
      *runsNormally = false;
   }
}

void MolDS::RotateMolecule(Molecule* molecule, bool* runsNormally) const{
   try{
      molecule->Rotate();
   }
   catch(MolDSException ex){
      this->OutputLog(boost::format("%s\n") % ex.what());
      ex.PrintBacktrace();
      *runsNormally = false;
   }
}

}

