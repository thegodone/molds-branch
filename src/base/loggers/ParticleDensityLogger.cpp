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
#include<string.h>
#include<math.h>
#include<vector>
#include<stdexcept>
#include<omp.h>
#include<boost/format.hpp>
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
#include"DensityLogger.h"
#include"ParticleDensityLogger.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_base_loggers{

ParticleDensityLogger::ParticleDensityLogger(const Molecule& molecule, 
                                             double const* const* fockMatrix, 
                                             double const* const* cisMatrix, 
                                             TheoryType theory) 
                     : DensityLogger(molecule, fockMatrix, cisMatrix, theory){
   this->SetMessages();
   //this->OutputLog("Particle density logger is created.\n");
}

ParticleDensityLogger::ParticleDensityLogger(){
   //this->OutputLog("Particle density logger is created.\n");
}

ParticleDensityLogger::~ParticleDensityLogger(){
   //this->OutputLog("Particle density logger is deleted.\n");
}

void ParticleDensityLogger::SetMessages(){
   DensityLogger::SetMessages();
   this->errorMessageCISMatrixNULL
      = "Error in base::logger::ParticleDensityPlot::DrawDensity: CIS Matrix is NULL.\n";
   this->errorMessageFockMatrixNULL
      = "Error in base::logger::ParticleDensityPlot::DrawDensity: Fock Matrix is NULL.\n";
   this->messageCubeHeaderComment1 = "MolDS cube file (in atomic units) for Particle density.\n";
   this->messageStartDensityPlot = "\t== START: Particle density plot ==\n";
   this->messageEndDensityPlot = "\t== DONE: Particle density plot ==\n\n";
   this->messageOmpElapsedTimeDensityPlot = "\t\tElapsed time(omp) for the Particle density plot = ";
}

double ParticleDensityLogger::GetDensityValue(int elecStateIndex, 
                                              double const* const* const* const* activeOccMOs,
                                              double const* const* const* const* activeVirMOs,
                                              double const* const* cisMatrix, 
                                              int ix, 
                                              int iy, 
                                              int iz)  const{
   double density = 0.0;
   int excitedStateIndex = elecStateIndex-1;
   int numberActiveOcc = Parameters::GetInstance()->GetActiveOccCIS();
   int numberActiveVir = Parameters::GetInstance()->GetActiveVirCIS();
   for(int i=0; i<numberActiveOcc; i++){
      for(int a=0; a<numberActiveVir; a++){
         for(int b=0; b<numberActiveVir; b++){
            int slaterDeterminatIndexIA = i*numberActiveVir + a;
            int slaterDeterminatIndexIB = i*numberActiveVir + b;
            double moAValue = activeVirMOs[a][ix][iy][iz];
            double moBValue = activeVirMOs[b][ix][iy][iz];
            density += moAValue*cisMatrix[excitedStateIndex][slaterDeterminatIndexIA]
                      *moBValue*cisMatrix[excitedStateIndex][slaterDeterminatIndexIB];
         }
      }
   }
   return density;
}

string ParticleDensityLogger::GetFileName(int elecStateIndex, int digit) const{
   stringstream fileName;
   fileName << Parameters::GetInstance()->GetFileNamePrefixParticlePlot();
   fileName << Utilities::Num2String(elecStateIndex,digit);
   fileName << this->stringCubeExtension;
   return fileName.str();
}

}
