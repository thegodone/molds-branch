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
#include"MOLogger.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_base_loggers{

MOLogger::MOLogger(const Molecule& molecule, 
                   double const* const* fockMatrix, 
                   TheoryType theory){
   this->molecule = &molecule;
   this->fockMatrix = fockMatrix;
   this->theory = theory;
   this->SetMessages();
}

MOLogger::MOLogger(){
}

void MOLogger::SetMessages(){
   this->stringCubeExtension = ".cube";
   this->errorMessageFockMatrixNULL
      = "Error in base::logger::MOLogger::DrawMO: Fock Matrix is NULL.\n";
   this->messageCubeHeaderComment1 = "MolDS cube file (in atomic units).\n";
   this->messageCubeHeaderComment2 = "outer loop:x, middle loop:y, inner loop:z\n";
   this->messageStartMOPlot = "\t== START: MO Plot ==\n";
   this->messageEndMOPlot = "\t== DONE: MO Plot ==\n\n";
   this->messageSkippedMOIndex = "\t\tBad MO-index is skipped. The skipped MO-index: ";
   this->messageOmpElapsedTimeMOPlot = "\t\tElapsed time(omp) for the MO plot = ";
   this->messageUnitSec = "[s].";
}

void MOLogger::DrawMO(int moIndex){
   vector<int> moIndeces;
   moIndeces.push_back(moIndex);
   this->DrawMO(moIndeces);
}

void MOLogger::DrawMO(vector<int> moIndeces){
   this->MatricesNullCheck();
   this->OutputLog(this->messageStartMOPlot);
   double ompStartTime = omp_get_wtime();

   // set frame basics
   double dx=0.0, dy=0.0, dz=0.0;
   double origin[CartesianType_end] = {0.0, 0.0, 0.0};
   this->CalcGridDisplacement(&dx, &dy, &dz);
   this->CalcOrigin(origin);

   // MO output 
   stringstream ompErrors;
#pragma omp parallel for schedule(dynamic, MOLDS_OMP_DYNAMIC_CHUNK_SIZE) 
   for(int i=0; i<moIndeces.size(); i++){
      try{
         // validate mo number
         if(this->molecule->GetTotalNumberAOs() <= moIndeces[i]){
            this->OutputLog(boost::format("%s%d\n") % this->messageSkippedMOIndex.c_str() 
                                                    % moIndeces[i]) ;
            continue;
         }
      
         // open the cube file
         int digit = 5;
         string fileName = this->GetFileName(moIndeces[i], digit);
         ofstream ofs(fileName.c_str());
      
         // output feader and molecule to the cube file
         this->OutputHeaderToFile(ofs, origin, dx, dy, dz);
         this->OutputMoleculeToFile(ofs, *this->molecule);
      
         // output grid data to the cube file
         int lineBreakCounter=0;
         for(int ix=0; ix<Parameters::GetInstance()->GetGridNumberMOPlot()[XAxis]; ix++){
            double x = origin[XAxis] + dx*static_cast<double>(ix);
            for(int iy=0; iy<Parameters::GetInstance()->GetGridNumberMOPlot()[YAxis]; iy++){
               double y = origin[YAxis] + dy*static_cast<double>(iy);
               for(int iz=0; iz<Parameters::GetInstance()->GetGridNumberMOPlot()[ZAxis]; iz++){
                  double z = origin[ZAxis] + dz*static_cast<double>(iz);
      
                  double moValue = this->GetMoValue(moIndeces[i], *this->molecule, this->fockMatrix, x, y, z);
                  ofs << boost::format("\t%e") % moValue ;
                  lineBreakCounter++;
                  if(lineBreakCounter%6==0){
                     ofs << endl;
                     lineBreakCounter=0;
                  }
      
               }
            }
         }
      }
      catch(MolDSException ex){
#pragma omp critical
         ex.Serialize(ompErrors);
      }
   }
   // Exception throwing for omp-region
   if(!ompErrors.str().empty()){
      throw MolDSException::Deserialize(ompErrors);
   }

   double ompEndTime = omp_get_wtime();
   this->OutputLog(boost::format("%s%lf%s\n%s") % this->messageOmpElapsedTimeMOPlot.c_str()
                                                % (ompEndTime - ompStartTime)
                                                % this->messageUnitSec.c_str()
                                                % this->messageEndMOPlot.c_str());
}

string MOLogger::GetFileName(int moIndex, int digit) const{
   stringstream fileName;
   fileName << Parameters::GetInstance()->GetFileNamePrefixMOPlot();
   fileName << Utilities::Num2String(moIndex,digit);
   fileName << this->stringCubeExtension;
   return fileName.str();
}

void MOLogger::OutputHeaderToFile(ofstream& ofs, double const* origin, double dx, double dy, double dz) const{
   int gridNumber[CartesianType_end] = {Parameters::GetInstance()->GetGridNumberMOPlot()[XAxis], 
                                        Parameters::GetInstance()->GetGridNumberMOPlot()[YAxis],
                                        Parameters::GetInstance()->GetGridNumberMOPlot()[ZAxis]};
   char data[1000] = "";
   // output header to the cube file
   ofs << this->messageCubeHeaderComment1;
   ofs << this->messageCubeHeaderComment2;
   sprintf(data,"\t%ld\t%e\t%e\t%e\n", this->molecule->GetAtomVect().size(),
                                       origin[XAxis], 
                                       origin[YAxis], 
                                       origin[ZAxis]);
   ofs << string(data);
   memset(data,0,sizeof(data));
   sprintf(data,"\t%d\t%e\t%e\t%e\n", gridNumber[XAxis], dx, 0.0, 0.0);
   ofs << string(data);
   memset(data,0,sizeof(data));
   sprintf(data,"\t%d\t%e\t%e\t%e\n", gridNumber[YAxis], 0.0, dy, 0.0);
   ofs << string(data);
   memset(data,0,sizeof(data));
   sprintf(data,"\t%d\t%e\t%e\t%e\n", gridNumber[ZAxis], 0.0, 0.0, dz);
   ofs << string(data);
}

void MOLogger::OutputMoleculeToFile(ofstream& ofs, const Molecule& molecule) const{
   char data[1000] = "";
   // output molecule to the cube file
   for(int a=0; a<molecule.GetAtomVect().size(); a++){
      Atom* atomA = molecule.GetAtomVect()[a];
      memset(data,0,sizeof(data));
      sprintf(data,"\t%d\t%d\t%e\t%e\t%e\n", atomA->GetAtomType()+1, 
                                       atomA->GetNumberValenceElectrons(),
                                       atomA->GetXyz()[XAxis],
                                       atomA->GetXyz()[YAxis],
                                       atomA->GetXyz()[ZAxis]);
      ofs << string(data);
   }
}

void MOLogger::CalcOrigin(double* origin) const{
   for(int i=0; i<CartesianType_end; i++){
      origin[i] = this->molecule->GetXyzCOC()[i];
      origin[i] -= 0.5*Parameters::GetInstance()->GetFrameLengthMOPlot()[i];
   }
}

void MOLogger::MatricesNullCheck() const{
   // NULL check
   if(this->fockMatrix == NULL){
      stringstream ss; 
      ss << this->errorMessageFockMatrixNULL;
      throw MolDSException(ss.str());
   }   
}

void MOLogger::CalcGridDisplacement(double* dx, double* dy, double* dz) const{
   *dx = Parameters::GetInstance()->GetFrameLengthMOPlot()[XAxis]
        /static_cast<double>(Parameters::GetInstance()->GetGridNumberMOPlot()[XAxis]);
   *dy = Parameters::GetInstance()->GetFrameLengthMOPlot()[YAxis]
        /static_cast<double>(Parameters::GetInstance()->GetGridNumberMOPlot()[YAxis]);
   *dz = Parameters::GetInstance()->GetFrameLengthMOPlot()[ZAxis]
        /static_cast<double>(Parameters::GetInstance()->GetGridNumberMOPlot()[ZAxis]);
}

double MOLogger::GetMoValue(int moIndex, 
                            const MolDS_base::Molecule& molecule, 
                            double const* const* fockMatrix, 
                            double x, 
                            double y, 
                            double z) const{
   double moValue = 0.0;
   for(int a=0; a<molecule.GetAtomVect().size(); a++){
      Atom* atomA = molecule.GetAtomVect()[a];
      int firstAOIndexA = atomA->GetFirstAOIndex();
      int numberAOsA = atomA->GetValenceSize();
      for(int mu=firstAOIndexA; mu<firstAOIndexA+numberAOsA; mu++){
         double aoValue = atomA->GetAtomicBasisValue(x,
                                                     y,
                                                     z,
                                                     mu-firstAOIndexA,
                                                     this->theory);
         moValue += fockMatrix[moIndex][mu]*aoValue;
      }
   }
   return moValue;
}

}
