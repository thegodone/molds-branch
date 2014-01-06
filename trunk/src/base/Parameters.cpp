//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   //
// Copyright (C) 2012-2014 Katsuhiko Nishimra                             // 
// Copyright (C) 2012-2014 Michihiro Okuyama
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
#include<boost/format.hpp>
#include"Enums.h"
#include"Uncopyable.h"
#include"PrintController.h"
#include"MolDSException.h"
#include"MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
#include"MallocerFreer.h"
#include"EularAngle.h"
#include"Parameters.h"
using namespace std;
namespace MolDS_base{

Parameters* Parameters::parameters = NULL;
// Physical constsnts
const double Parameters::eV2AU         = 0.03674903;
const double Parameters::j2AU          = 1.0e18/4.35974394;
const double Parameters::kcalMolin2AU  = 0.00159360175;
const double Parameters::angstrom2AU   = 1.0/0.5291772;
const double Parameters::nm2AU         = 10.0*Parameters::angstrom2AU;
const double Parameters::kayser2AU     = 4.556336e-6;
const double Parameters::fs2AU         = 1.0/(2.418884326505e-2);
const double Parameters::gMolin2AU     = 1.0e5/(6.0221415*9.1095);
const double Parameters::degree2Radian = M_PI / 180.0;
const double Parameters::boltzmann     = 3.166791e-6;
const double Parameters::avogadro      = 6.0221415e23;
const double Parameters::debye2AU      = 0.393430191;

// constant
const double Parameters::vdWScalingFactorSCFPM3DAM1D = 1.40;
const double Parameters::vdWDampingFactorSCFPM3DAM1D = 23.0;


Parameters::Parameters(){
   this->SetDefaultValues();
   this->SetMessages();
   this->indecesMOPlot           = NULL;
   this->elecIndecesHolePlot     = NULL;
   this->elecIndecesParticlePlot = NULL;
   this->electronicStateIndecesMullikenCIS = NULL;
}

Parameters::~Parameters(){
   MallocerFreer::GetInstance()->Free<double>(&this->inertiaTensorOrigin,
                                              CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(&this->rotatingOrigin,
                                              CartesianType_end);
   if(this->indecesMOPlot != NULL){
      delete this->indecesMOPlot;
      this->indecesMOPlot = NULL;
      //this->OutputLog("indecesMOPlot deleted\n");
   }
   if(this->elecIndecesHolePlot != NULL){
      delete this->elecIndecesHolePlot;
      this->elecIndecesHolePlot = NULL;
      //this->OutputLog("elecIndecesHolePlot deleted\n");
   }
   if(this->elecIndecesParticlePlot != NULL){
      delete this->elecIndecesParticlePlot;
      this->elecIndecesParticlePlot = NULL;
      //this->OutputLog("elecIndecesParticlePlot deleted\n");
   }
   if(this->electronicStateIndecesMullikenCIS != NULL){
      delete this->electronicStateIndecesMullikenCIS;
      this->electronicStateIndecesMullikenCIS= NULL;
      //this->OutputLog("electronicStateIndecesMullikenCIS deleted\n");
   }
}

Parameters* Parameters::GetInstance(){
   if(parameters == NULL){
      parameters = new Parameters();
   }
   return parameters;
}

void Parameters::DeleteInstance(){
   if(parameters != NULL){
      delete parameters; 
   }
   parameters = NULL;
}

void Parameters::SetDefaultValues(){
   this->currentSimulation = Once;
   this->currentTheory = CNDO2;
   // SCF
   this->thresholdSCF        = 1.0e-8;
   this->maxIterationsSCF    = 100;
   this->dampingThreshSCF    = 1.0;
   this->dampingWeightSCF    = 0.8;
   this->diisNumErrorVectSCF = 5;
   this->diisStartErrorSCF   = 1.0e-2;
   this->diisEndErrorSCF     = 1.0e-8;
   this->requiresVdWSCF      = false;
   this->vdWScalingFactorSCF = 1.40;
   this->vdWDampingFactorSCF = 23.0;
   // MOPlot
   this->fileNamePrefixMOPlot     = "MO_";
   this->gridNumberMOPlot[XAxis]  = 25;
   this->gridNumberMOPlot[YAxis]  = 25;
   this->gridNumberMOPlot[ZAxis]  = 25;
   this->frameLengthMOPlot[XAxis] = 20.0;
   this->frameLengthMOPlot[YAxis] = 20.0;
   this->frameLengthMOPlot[ZAxis] = 20.0;
   // HolePlot
   this->fileNamePrefixHolePlot     = "hole_";
   this->gridNumberHolePlot[XAxis]  = 25;
   this->gridNumberHolePlot[YAxis]  = 25;
   this->gridNumberHolePlot[ZAxis]  = 25;
   this->frameLengthHolePlot[XAxis] = 20.0;
   this->frameLengthHolePlot[YAxis] = 20.0;
   this->frameLengthHolePlot[ZAxis] = 20.0;
   // ParticlePlot
   this->fileNamePrefixParticlePlot     = "particle_";
   this->gridNumberParticlePlot[XAxis]  = 25;
   this->gridNumberParticlePlot[YAxis]  = 25;
   this->gridNumberParticlePlot[ZAxis]  = 25;
   this->frameLengthParticlePlot[XAxis] = 20.0;
   this->frameLengthParticlePlot[YAxis] = 20.0;
   this->frameLengthParticlePlot[ZAxis] = 20.0;
   // Translation
   this->translatingDifference[0] = 0.0;
   this->translatingDifference[1] = 0.0;
   this->translatingDifference[2] = 0.0;
   // Principal axes
   this->inertiaTensorOrigin = NULL;
   // Rotation
   this->rotatingOrigin  = NULL;
   this->rotatingAxis[0] = 0.0;
   this->rotatingAxis[1] = 0.0;
   this->rotatingAxis[2] = 1.0;
   this->rotatingType    = Axis;
   this->rotatingEularAngles.SetAlpha(0.0);
   this->rotatingEularAngles.SetBeta(0.0);
   this->rotatingEularAngles.SetGamma(0.0);
   // CIS
   this->activeOccCIS                          = 10;
   this->activeVirCIS                          = 10;
   this->numberExcitedStatesCIS                = 5;
   this->requiresCIS                           = false;
   this->isDavidsonCIS                         = true;
   this->maxIterationsCIS                      = 100;
   this->maxDimensionsCIS                      = 100;
   this->normToleranceCIS                      = 1.0e-6;
   this->numberPrintCoefficientsCIS            = 1;
   this->requiresExcitonEnergiesCIS            = false;
   this->requiresAllTransitionDipoleMomentsCIS = false;
   this->requiresUnpairedPopCIS                = false;
   // Memory
   this->limitHeapMemory = 256;
   // MD
   this->electronicStateIndexMD = 0;
   this->totalStepsMD           = 10;
   this->timeWidthMD            = 0.1*this->fs2AU;
   // MC
   this->electronicStateIndexMC = 0;
   this->totalStepsMC           = 10;
   this->stepWidthMC            = 0.05*this->angstrom2AU;
   this->temperatureMC          = 300;
   this->seedMC                 = static_cast<unsigned long>(time(0));
   // RPMD
   this->electronicStateIndexRPMD   = 0;
   this->numberElectronicStatesRPMD = 1;
   this->totalStepsRPMD             = 10;
   this->timeWidthRPMD              = 0.1*this->fs2AU;
   this->temperatureRPMD            = 300;
   this->numberBeadsRPMD            = 10;
   this->seedRPMD                   = static_cast<unsigned long>(time(0));
   // NASCO
   this->numberElectronicStatesNASCO = 3;
   this->initialElectronicStateNASCO = 0;
   this->totalStepsNASCO             = 10;
   this->timeWidthNASCO              = 0.1*this->fs2AU;
   this->seedNASCO                   = static_cast<unsigned long>(time(0));
   // Optimization 
   this->methodOptimization               = ConjugateGradientMethod;
   this->totalStepsOptimization           = 50;
   this->electronicStateIndexOptimization = 0;
   this->maxGradientOptimization          = 0.00045;
   this->rmsGradientOptimization          = 0.00030;
   this->timeWidthOptimization            = 50.0*this->fs2AU;
   this->initialTrustRadiusOptimization   = 0.3;
   this->maxNormStepOptimization          = 0.3;
   // Frequencies
   this->requiresFrequencies             = false;
   this->electronicStateIndexFrequencies = 0;
}

void Parameters::SetMessages(){
   this->errorMessageGetIndecesMOPlotNull
      = "Error in base::Parameters::GetIndecesMOPlot: indecesMOPlot is NULL.\n";
   this->errorMessageGetIndecesHolePlotNull
      = "Error in base::Parameters::GetIndecesHolePlot: elecIndecesHolePlot is NULL.\n";
   this->errorMessageGetIndecesParticlePlotNull
      = "Error in base::Parameters::GetIndecesParticlePlot: elecIndecesParticlePlot is NULL.\n";
   this->errorMessageGetElectronicStateIndecesMullikenCISNull
      = "Error in base::Parameters::GetElectronicStateIndecesMullikenCIS: electronicStateIndecesMullikenCIS is NULL.\n";
}

// methods for translation
void Parameters::SetTranslatingDifference(double x, double y, double z){
   this->translatingDifference[0] = x;
   this->translatingDifference[1] = y;
   this->translatingDifference[2] = z;
}

// methods for principal axes
void Parameters::SetInertiaTensorOrigin(double x, double y, double z){
   MallocerFreer::GetInstance()->Malloc<double>(&this->inertiaTensorOrigin, CartesianType_end);
   this->inertiaTensorOrigin[0] = x;
   this->inertiaTensorOrigin[1] = y;
   this->inertiaTensorOrigin[2] = z;
}

// methods for rotation
void Parameters::SetRotatingOrigin(double x, double y, double z){
   MallocerFreer::GetInstance()->Malloc<double>(&this->rotatingOrigin, CartesianType_end);
   this->rotatingOrigin[0] = x;
   this->rotatingOrigin[1] = y;
   this->rotatingOrigin[2] = z;
}

void Parameters::SetRotatingAxis(double x, double y, double z){
   this->rotatingAxis[0] = x;
   this->rotatingAxis[1] = y;
   this->rotatingAxis[2] = z;
}

void Parameters::SetRotatingEularAngles(double alpha, double beta, double gamma){
   this->rotatingEularAngles.SetAlpha(alpha);
   this->rotatingEularAngles.SetBeta(beta);
   this->rotatingEularAngles.SetGamma(gamma);
}

// methods for MOPlot
vector<int>* Parameters::GetIndecesMOPlot() const{
#ifdef MOLDS_DBG
   if(this->indecesMOPlot==NULL) throw MolDSException(this->errorMessageGetIndecesMOPlotNull);
#endif
   return this->indecesMOPlot;
}

void Parameters::AddIndexMOPlot(int moIndex){
   if(this->indecesMOPlot==NULL){
      this->indecesMOPlot = new vector<int>;
   }
   this->indecesMOPlot->push_back(moIndex);
}

void Parameters::SetGridNumberMOPlot(int Nx, int Ny, int Nz){
   this->gridNumberMOPlot[XAxis] = Nx;
   this->gridNumberMOPlot[YAxis] = Ny;
   this->gridNumberMOPlot[ZAxis] = Nz;
}

void Parameters::SetFrameLengthMOPlot(double lx, double ly, double lz){
   this->frameLengthMOPlot[XAxis] = lx;
   this->frameLengthMOPlot[YAxis] = ly;
   this->frameLengthMOPlot[ZAxis] = lz;
}

// methods for HolePlot
vector<int>* Parameters::GetElecIndecesHolePlot() const{
#ifdef MOLDS_DBG
   if(this->elecIndecesHolePlot==NULL) throw MolDSException(this->errorMessageGetIndecesHolePlotNull);
#endif
   return this->elecIndecesHolePlot;
}

void Parameters::AddElecIndexHolePlot(int elecIndex){
   if(this->elecIndecesHolePlot==NULL){
      this->elecIndecesHolePlot = new vector<int>;
   }
   this->elecIndecesHolePlot->push_back(elecIndex);
}

void Parameters::SetGridNumberHolePlot(int Nx, int Ny, int Nz){
   this->gridNumberHolePlot[XAxis] = Nx;
   this->gridNumberHolePlot[YAxis] = Ny;
   this->gridNumberHolePlot[ZAxis] = Nz;
}

void Parameters::SetFrameLengthHolePlot(double lx, double ly, double lz){
   this->frameLengthHolePlot[XAxis] = lx;
   this->frameLengthHolePlot[YAxis] = ly;
   this->frameLengthHolePlot[ZAxis] = lz;
}

// methods for ParticlePlot
const vector<int>* Parameters::GetElecIndecesParticlePlot() const{
#ifdef MOLDS_DBG
   if(this->elecIndecesParticlePlot==NULL) throw MolDSException(this->errorMessageGetIndecesParticlePlotNull);
#endif
   return this->elecIndecesParticlePlot;
}

void Parameters::AddElecIndexParticlePlot(int elecIndex){
   if(this->elecIndecesParticlePlot==NULL){
      this->elecIndecesParticlePlot = new vector<int>;
   }
   this->elecIndecesParticlePlot->push_back(elecIndex);
}

void Parameters::SetGridNumberParticlePlot(int Nx, int Ny, int Nz){
   this->gridNumberParticlePlot[XAxis] = Nx;
   this->gridNumberParticlePlot[YAxis] = Ny;
   this->gridNumberParticlePlot[ZAxis] = Nz;
}

void Parameters::SetFrameLengthParticlePlot(double lx, double ly, double lz){
   this->frameLengthParticlePlot[XAxis] = lx;
   this->frameLengthParticlePlot[YAxis] = ly;
   this->frameLengthParticlePlot[ZAxis] = lz;
}

// methods for CIS
vector<int>* Parameters::GetElectronicStateIndecesMullikenCIS() const{
#ifdef MOLDS_DBG
   if(this->electronicStateIndecesMullikenCIS==NULL) throw MolDSException(this->errorMessageGetElectronicStateIndecesMullikenCISNull);
#endif
   return this->electronicStateIndecesMullikenCIS;
}

void Parameters::AddElectronicStateIndexMullikenCIS(int electronicStateIndex){
   if(this->electronicStateIndecesMullikenCIS==NULL){
      this->electronicStateIndecesMullikenCIS = new vector<int>;
   }
   this->electronicStateIndecesMullikenCIS->push_back(electronicStateIndex);
}

bool Parameters::RequiresMullikenCIS() const{
   return (this->electronicStateIndecesMullikenCIS!=NULL && 
           0<this->electronicStateIndecesMullikenCIS->size());
}

}





