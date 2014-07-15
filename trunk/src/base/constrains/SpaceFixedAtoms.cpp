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
#include<fstream>
#include<string>
#include<math.h>
#include<vector>
#include<stdexcept>
#include<omp.h>
#include<boost/shared_ptr.hpp>
#include<boost/format.hpp>
#include"../../config.h"
#include"../Enums.h"
#include"../Uncopyable.h"
#include"../PrintController.h"
#include"../MolDSException.h"
#include"../MallocerFreer.h"
#include"../../mpi/MpiInt.h"
#include"../../mpi/MpiProcess.h"
#include"../Utilities.h"
#include"../EularAngle.h"
#include"../Parameters.h"
#include"../RealSphericalHarmonicsIndex.h"
#include"../atoms/Atom.h"
#include"../Molecule.h"
#include"../ElectronicStructure.h"
#include"Constrain.h"
#include"SpaceFixedAtoms.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_base_constrains{

SpaceFixedAtoms::SpaceFixedAtoms(const MolDS_base::Molecule* molecule,
                                 const boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure)
                                 : Constrain(molecule, electronicStructure){
   //this->OutputLog("SpaceFixedAtoms created\n");
}

SpaceFixedAtoms::~SpaceFixedAtoms(){
   if(this->constrainedMatrixForce!=NULL){
      MallocerFreer::GetInstance()->Free<double>(&this->constrainedMatrixForce, 
                                                 this->molecule->GetAtomVect().size(),
                                                 CartesianType_end);
      //this->OutputLog("SpaceFixedAtoms:constrainedMatrixForce freed\n");
   }
   if(this->refMolecule!=NULL){
      delete this->refMolecule;
      //this->OutputLog("refMolecule deleted\n");
   }
   //this->OutputLog("SpaceFixedAtoms deleted\n");
}

void SpaceFixedAtoms::SetConstrainCondition(){
   int atomNum = this->molecule->GetAtomVect().size();
   MallocerFreer::GetInstance()->Malloc<double>(&this->constrainedMatrixForce, 
                                                atomNum,
                                                CartesianType_end);
   this->refMolecule = new MolDS_base::Molecule(*this->molecule);

   if(Parameters::GetInstance()->RequiresSpaceFixedAtomsOptimization()){
      // create unique index list of fixed atoms
      const vector<AtomIndexPair>* indexPairs = Parameters::GetInstance()->GetSpaceFixedAtomIndexPairsOptimization();
      for(vector<AtomIndexPair>::const_iterator itr=indexPairs->begin();itr!=indexPairs->end();++itr){
         for(int j=itr->firstAtomIndex; j<=itr->lastAtomIndex; j++){
            this->fixedAtomIndeces.push_back(j);
         }
      }
      std::sort(this->fixedAtomIndeces.begin(), this->fixedAtomIndeces.end()); 
      vector<int>::iterator unq = std::unique(this->fixedAtomIndeces.begin(), this->fixedAtomIndeces.end());
      this->fixedAtomIndeces.erase(unq, this->fixedAtomIndeces.end());
   }
}

double const* const* SpaceFixedAtoms::GetForce(int elecState){
   double const* const* matrixForce = this->electronicStructure->GetForce(elecState);
   int atomNum = this->molecule->GetAtomVect().size();
#pragma omp parallel for schedule(dynamic, MOLDS_OMP_DYNAMIC_CHUNK_SIZE) 
   for(int a=0; a<atomNum; a++){
      for(int i=0; i<CartesianType_end; i++){
         constrainedMatrixForce[a][i] = matrixForce[a][i];
      }
   }
   if(Parameters::GetInstance()->RequiresSpaceFixedAtomsOptimization()){
      for(vector<int>::iterator itr=this->fixedAtomIndeces.begin();itr!=this->fixedAtomIndeces.end();++itr){
         for(int i=0; i<CartesianType_end; i++){
            constrainedMatrixForce[*itr][i] = 0.0;
         }
      }
   }
   return this->constrainedMatrixForce;
}


}
