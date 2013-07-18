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
#include<fstream>
#include<string>
#include<string.h>
#include<math.h>
#include<vector>
#include<stdexcept>
#include<omp.h>
#include<boost/format.hpp>
#include"../Uncopyable.h"
#include"../../mpi/MpiProcess.h"
#include"../PrintController.h"
#include"../MolDSException.h"
#include"../Utilities.h"
#include"../Enums.h"
#include"../MallocerFreer.h"
#include"../EularAngle.h"
#include"../Parameters.h"
#include"../atoms/Atom.h"
#include"../Molecule.h"
#include"DensityLogger.h"
using namespace std;
using namespace MolDS_base;
using namespace MolDS_base_atoms;
namespace MolDS_base_loggers{

DensityLogger::DensityLogger(const Molecule& molecule, 
                       double const* const* fockMatrix, 
                       double const* const* cisMatrix, 
                   TheoryType theory){
   this->molecule = &molecule;
   this->fockMatrix = fockMatrix;
   this->cisMatrix = cisMatrix;
   this->theory = theory;
   //this->OutputLog("Density logger is created.\n");
}

DensityLogger::DensityLogger(){
   //this->OutputLog("Density loger is created.\n");
}

DensityLogger::~DensityLogger(){
   //this->OutputLog("Density logger is deleted.\n");
}

void DensityLogger::SetMessages(){
   this->messageCubeHeaderComment2 = "outer loop:x, middle loop:y, inner loop:z\n";
   this->messageSkippedElecStateIndex = "\t\tBad electronic state is skipped. The skipped electronic state: ";
   this->messageUnitSec = "[s].";
   this->stringCubeExtension = ".cube";
}

void DensityLogger::DrawDensity(int elecStateIndex) const{
   vector<int> elecStateIndeces;
   elecStateIndeces.push_back(elecStateIndex);
   this->DrawDensity(elecStateIndeces);
}

void DensityLogger::DrawDensity(vector<int> elecStateIndeces) const{
   this->MatricesNullCheck();
   this->OutputLog(this->messageStartDensityPlot);
   double ompStartTime = omp_get_wtime();
   double**** activeOccMOs=NULL;
   double**** activeVirMOs=NULL;

   // set frame basics
   double dx=0.0, dy=0.0, dz=0.0;
   double origin[CartesianType_end] = {0.0, 0.0, 0.0};
   this->CalcGridDisplacement(&dx, &dy, &dz);
   this->CalcOrigin(origin);

   try{
      // calc. active MOs
      this->MallocTemporaryActiveMOs(&activeOccMOs, &activeVirMOs);
      this->CalcActiveMOs(activeOccMOs, 
                          activeVirMOs,
                          dx, dy, dz,
                          origin,
                          *this->molecule,
                          this->fockMatrix,
                          this->cisMatrix);

      // density output 
      stringstream ompErrors;
#pragma omp parallel for schedule(auto) 
      for(int n=0; n<elecStateIndeces.size(); n++){
         try{
            // validate electronic state
            int groundState = 0;
            if(Parameters::GetInstance()->GetNumberExcitedStatesCIS() < elecStateIndeces[n] || 
               groundState == elecStateIndeces[n]){
               this->OutputLog(boost::format("%s%d\n") % this->messageSkippedElecStateIndex.c_str() 
                                                       % elecStateIndeces[n]) ;
               continue;
            }
         
            // open the cube file
            int digit=5;
            string fileName = this->GetFileName(elecStateIndeces[n], digit);
            ofstream ofs(fileName.c_str());
         
            // output feader and molecule to the cube file
            this->OutputHeaderToFile(ofs, origin, dx, dy, dz);
            this->OutputMoleculeToFile(ofs, *this->molecule);
         
            // output grid data to the cube file
            int lineBreakCounter=0;
            for(int ix=0; ix<this->GetGridNumber()[XAxis]; ix++){
               double x = origin[XAxis] + dx*static_cast<double>(ix);
               for(int iy=0; iy<this->GetGridNumber()[YAxis]; iy++){
                  double y = origin[YAxis] + dy*static_cast<double>(iy);
                  for(int iz=0; iz<this->GetGridNumber()[ZAxis]; iz++){
                     double z = origin[ZAxis] + dz*static_cast<double>(iz);
         
                     double density = this->GetDensityValue(elecStateIndeces[n], 
                                                            activeOccMOs,
                                                            activeVirMOs,
                                                            this->cisMatrix, 
                                                            ix, iy, iz);
                     ofs << boost::format("\t%e") % density ;
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
   }
   catch(MolDSException ex){
      this->FreeTemporaryActiveMOs(&activeOccMOs, &activeVirMOs);
      throw ex;
   }
   this->FreeTemporaryActiveMOs(&activeOccMOs, &activeVirMOs);

   double ompEndTime = omp_get_wtime();
   this->OutputLog(boost::format("%s%lf%s\n%s") % this->messageOmpElapsedTimeDensityPlot.c_str()
                                                % (ompEndTime - ompStartTime)
                                                % this->messageUnitSec.c_str()
                                                % this->messageEndDensityPlot.c_str());
}

void DensityLogger::CalcActiveMOs(double**** activeOccMOs, 
                                  double**** activeVirMOs,
                                  double dx, double dy, double dz,
                                  double const* origin,
                                  const MolDS_base::Molecule& molecule,
                                  double const* const* fockMatrix,
                                  double const* const* cisMatrix) const{
   int numberOcc = molecule.GetTotalNumberValenceElectrons()/2;
   for(int i=0; i<Parameters::GetInstance()->GetActiveOccCIS(); i++){
      int moI = numberOcc - (i+1);
      for(int ix=0; ix<this->GetGridNumber()[XAxis]; ix++){
         double x = origin[XAxis] + dx*static_cast<double>(ix);
         for(int iy=0; iy<this->GetGridNumber()[YAxis]; iy++){
            double y = origin[YAxis] + dy*static_cast<double>(iy);
            for(int iz=0; iz<this->GetGridNumber()[ZAxis]; iz++){
               double z = origin[ZAxis] + dz*static_cast<double>(iz);
               activeOccMOs[i][ix][iy][iz] = this->GetMOValue(moI, molecule, fockMatrix, x, y, z);
            }
         }
      }
   }
   for(int a=0; a<Parameters::GetInstance()->GetActiveVirCIS(); a++){
      int moA = numberOcc + a;
      for(int ix=0; ix<this->GetGridNumber()[XAxis]; ix++){
         double x = origin[XAxis] + dx*static_cast<double>(ix);
         for(int iy=0; iy<this->GetGridNumber()[YAxis]; iy++){
            double y = origin[YAxis] + dy*static_cast<double>(iy);
            for(int iz=0; iz<this->GetGridNumber()[ZAxis]; iz++){
               double z = origin[ZAxis] + dz*static_cast<double>(iz);
               activeVirMOs[a][ix][iy][iz] = this->GetMOValue(moA, molecule, fockMatrix, x, y, z);
            }
         }
      }
   }
}

void DensityLogger::MallocTemporaryActiveMOs(double***** activeOccMOs, double***** activeVirMOs) const{
   MallocerFreer::GetInstance()->Malloc<double>(activeOccMOs, 
                                                Parameters::GetInstance()->GetActiveOccCIS(),
                                                this->GetGridNumber()[XAxis],
                                                this->GetGridNumber()[YAxis],
                                                this->GetGridNumber()[ZAxis]);
   MallocerFreer::GetInstance()->Malloc<double>(activeVirMOs, 
                                                Parameters::GetInstance()->GetActiveVirCIS(),
                                                this->GetGridNumber()[XAxis],
                                                this->GetGridNumber()[YAxis],
                                                this->GetGridNumber()[ZAxis]);
}

void DensityLogger::FreeTemporaryActiveMOs(double***** activeOccMOs, double***** activeVirMOs) const{
   MallocerFreer::GetInstance()->Free<double>(activeOccMOs, 
                                              Parameters::GetInstance()->GetActiveOccCIS(),
                                              this->GetGridNumber()[XAxis],
                                              this->GetGridNumber()[YAxis],
                                              this->GetGridNumber()[ZAxis]);
   MallocerFreer::GetInstance()->Free<double>(activeVirMOs, 
                                              Parameters::GetInstance()->GetActiveVirCIS(),
                                              this->GetGridNumber()[XAxis],
                                              this->GetGridNumber()[YAxis],
                                              this->GetGridNumber()[ZAxis]);
}

void DensityLogger::CalcGridDisplacement(double* dx, double* dy, double* dz) const{
   *dx = this->GetFrameLength()[XAxis]/static_cast<double>(this->GetGridNumber()[XAxis]);
   *dy = this->GetFrameLength()[YAxis]/static_cast<double>(this->GetGridNumber()[YAxis]);
   *dz = this->GetFrameLength()[ZAxis]/static_cast<double>(this->GetGridNumber()[ZAxis]);
}

void DensityLogger::CalcOrigin(double* origin) const{
   for(int i=0; i<CartesianType_end; i++){
      origin[i] = this->molecule->GetXyzCOC()[i];
      origin[i] -= 0.5*this->GetFrameLength()[i];
   }
}


double DensityLogger::GetMOValue(int moIndex, 
                                 const MolDS_base::Molecule& molecule, 
                                 double const* const* forckMatrix,
                                 double x, double y, double z) const{
   double moValue = 0.0;
   for(int a=0; a<this->molecule->GetNumberAtoms(); a++){
      Atom* atomA = this->molecule->GetAtom(a);
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

void DensityLogger::OutputHeaderToFile(ofstream& ofs, double const* origin, double dx, double dy, double dz) const{
   int gridNumber[CartesianType_end] = {this->GetGridNumber()[XAxis], 
                                        this->GetGridNumber()[YAxis],
                                        this->GetGridNumber()[ZAxis]};
   char data[1000] = "";
   // output header to the cube file
   ofs << this->messageCubeHeaderComment1;
   ofs << this->messageCubeHeaderComment2;
   sprintf(data,"\t%d\t%e\t%e\t%e\n", this->molecule->GetNumberAtoms(),
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

void DensityLogger::OutputMoleculeToFile(ofstream& ofs, const Molecule& molecule) const{
   char data[1000] = "";
   // output molecule to the cube file
   for(int a=0; a<molecule.GetNumberAtoms(); a++){
      const Atom& atomA = *molecule.GetAtom(a);
      memset(data,0,sizeof(data));
      sprintf(data,"\t%d\t%d\t%e\t%e\t%e\n", atomA.GetAtomType()+1, 
                                       atomA.GetNumberValenceElectrons(),
                                       atomA.GetXyz()[XAxis],
                                       atomA.GetXyz()[YAxis],
                                       atomA.GetXyz()[ZAxis]);
      ofs << string(data);
   }
}

void DensityLogger::MatricesNullCheck() const{
   // NULL check
   if(this->cisMatrix == NULL){
      stringstream ss;
      ss << this->errorMessageCISMatrixNULL;
      throw MolDSException(ss.str());
   }
   if(this->fockMatrix == NULL){
      stringstream ss;
      ss << this->errorMessageFockMatrixNULL;
      throw MolDSException(ss.str());
   }
}

}
