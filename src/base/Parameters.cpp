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
#include<boost/format.hpp>
#include"PrintController.h"
#include"MolDSException.h"
#include"Uncopyable.h"
#include"Enums.h"
#include"MallocerFreer.h"
#include"EularAngle.h"
#include"Parameters.h"
using namespace std;
namespace MolDS_base{

Parameters* Parameters::parameters = NULL;
// Physical constsnts
const double Parameters::eV2AU         = 0.03674903;
const double Parameters::j2AU          = pow(10.0,18.0)/4.35974394;
const double Parameters::kcalMolin2AU  = 0.00159360175;
const double Parameters::angstrom2AU   = 1.0/0.5291772;
const double Parameters::nm2AU         = 10.0*Parameters::angstrom2AU;
const double Parameters::kayser2AU     = 4.556336*pow(10.0,-6.0);
const double Parameters::fs2AU         = 1.0/(2.418884326505*pow(10.0,-2.0));
const double Parameters::gMolin2AU     = pow(10.0,5.0)/(6.0221415*9.1095);
const double Parameters::degree2Radian = M_PI / 180.0;
const double Parameters::boltzmann     = 3.166791*pow(10.0,-6.0);
const double Parameters::avogadro      = 6.0221415*pow(10.0, 23.0);
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
   this->thresholdSCF        = pow(10.0, -8.0);
   this->maxIterationsSCF    = 100;
   this->dampingThreshSCF    = 1.0;
   this->dampingWeightSCF    = 0.8;
   this->diisNumErrorVectSCF = 5;
   this->diisStartErrorSCF   = pow(10.0, -2.0);
   this->diisEndErrorSCF     = pow(10.0, -8.0);
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
   this->normToleranceCIS                      = pow(10.0, -6.0);
   this->numberPrintCoefficientsCIS            = 1;
   this->requiresExcitonEnergiesCIS            = false;
   this->requiresAllTransitionDipoleMomentsCIS = false;
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
}

SimulationType Parameters::GetCurrentSimulation() const{
   return this->currentSimulation;
}

void Parameters::SetCurrentSimulation(SimulationType simulation){
   this->currentSimulation = simulation;
}

TheoryType Parameters::GetCurrentTheory() const{
   return this->currentTheory;
}

void Parameters::SetCurrentTheory(TheoryType theory){
   this->currentTheory = theory;
}

// methods for SCF
double Parameters::GetThresholdSCF() const{
   return this->thresholdSCF;
}

void Parameters::SetThresholdSCF(double thresholdSCF){
   this->thresholdSCF = thresholdSCF;
}

int Parameters::GetMaxIterationsSCF() const{
   return this->maxIterationsSCF;
}

void Parameters::SetMaxIterationsSCF(int maxIterationsSCF){
   this->maxIterationsSCF = maxIterationsSCF;
}

double Parameters::GetDampingThreshSCF() const{
   return this->dampingThreshSCF;
}

void Parameters::SetDampingThreshSCF(double dampingThreshSCF){
   this->dampingThreshSCF = dampingThreshSCF;
}

double Parameters::GetDampingWeightSCF() const{
   return this->dampingWeightSCF;
}

void Parameters::SetDampingWeightSCF(double dampingWeightSCF){
   this->dampingWeightSCF = dampingWeightSCF;
}

int Parameters::GetDiisNumErrorVectSCF() const{
   return this->diisNumErrorVectSCF;
}

void Parameters::SetDiisNumErrorVectSCF(int diisNumErrorVectSCF){
   this->diisNumErrorVectSCF = diisNumErrorVectSCF;
}

double Parameters::GetDiisStartErrorSCF() const{
   return this->diisStartErrorSCF;
}

void Parameters::SetDiisStartErrorSCF(double diisStartErrorSCF){
   this->diisStartErrorSCF = diisStartErrorSCF;
}

double Parameters::GetDiisEndErrorSCF() const{
   return this->diisEndErrorSCF;
}

void Parameters::SetDiisEndErrorSCF(double diisEndErrorSCF){
   this->diisEndErrorSCF = diisEndErrorSCF;
}

bool Parameters::RequiresVdWSCF() const{
   return this->requiresVdWSCF;
}

void Parameters::SetRequiresVdWSCF(bool requiresVdWSCF){
   this->requiresVdWSCF = requiresVdWSCF;
}

double Parameters::GetVdWScalingFactorSCF() const{
   return this->vdWScalingFactorSCF;
}

void Parameters::SetVdWScalingFactorSCF(){
   this->vdWScalingFactorSCF = this->vdWScalingFactorSCFPM3DAM1D;
}

void Parameters::SetVdWScalingFactorSCF(double vdWScalingFactorSCF){
   this->vdWScalingFactorSCF = vdWScalingFactorSCF;
}

double Parameters::GetVdWDampingFactorSCF() const{
   return this->vdWDampingFactorSCF;
}

void Parameters::SetVdWDampingFactorSCF(){
   this->vdWDampingFactorSCF = this->vdWDampingFactorSCFPM3DAM1D;
}

void Parameters::SetVdWDampingFactorSCF(double vdWScalingFactorSCF){
   this->vdWDampingFactorSCF = vdWScalingFactorSCF;
}

// methods for physical constants
double Parameters::GetEV2AU() const{
   return this->eV2AU;
}

double Parameters::GetJ2AU() const{
   return this->j2AU;
}

double Parameters::GetKcalMolin2AU() const{
   return this->kcalMolin2AU;
}

double Parameters::GetAngstrom2AU() const{
   return this->angstrom2AU;
}

double Parameters::GetNm2AU() const{
   return this->nm2AU;
}

double Parameters::GetKayser2AU() const{
   return this->kayser2AU;
}

double Parameters::GetGMolin2AU() const{
   return this->gMolin2AU;
}

double Parameters::GetDegree2Radian() const{
   return this->degree2Radian;
}

double Parameters::GetFs2AU() const{
   return this->fs2AU;
}

double Parameters::GetBoltzmann() const{
   return this->boltzmann;
}

double Parameters::GetAvogadro() const{
   return this->avogadro;
}

double Parameters::GetDebye2AU() const{
   return this->debye2AU;
}

// methods for translation
void Parameters::SetTranslatingDifference(double x, double y, double z){
   this->translatingDifference[0] = x;
   this->translatingDifference[1] = y;
   this->translatingDifference[2] = z;
}

double* Parameters::GetTranslatingDifference() const{
   return (double*)this->translatingDifference;
}

// methods for principal axes
void Parameters::SetInertiaTensorOrigin(double x, double y, double z){
   MallocerFreer::GetInstance()->Malloc<double>(&this->inertiaTensorOrigin, CartesianType_end);
   this->inertiaTensorOrigin[0] = x;
   this->inertiaTensorOrigin[1] = y;
   this->inertiaTensorOrigin[2] = z;
}

double* Parameters::GetInertiaTensorOrigin() const{
   return (double*)this->inertiaTensorOrigin;
}

// methods for rotation
void Parameters::SetRotatingOrigin(double x, double y, double z){
   MallocerFreer::GetInstance()->Malloc<double>(&this->rotatingOrigin, CartesianType_end);
   this->rotatingOrigin[0] = x;
   this->rotatingOrigin[1] = y;
   this->rotatingOrigin[2] = z;
}

double* Parameters::GetRotatingOrigin() const{
   return (double*)this->rotatingOrigin;
}

void Parameters::SetRotatingType(RotatingType rotatingType){
   this->rotatingType = rotatingType;
}

RotatingType Parameters::GetRotatingType() const{
   return this->rotatingType;
}

void Parameters::SetRotatingAxis(double x, double y, double z){
   this->rotatingAxis[0] = x;
   this->rotatingAxis[1] = y;
   this->rotatingAxis[2] = z;

}

double* Parameters::GetRotatingAxis() const{
   return (double*)this->rotatingAxis;
}

void Parameters::SetRotatingAngle(double rotatingAngle){
   this->rotatingAngle = rotatingAngle;
}

double Parameters::GetRotatingAngle() const{
   return this->rotatingAngle;
}

void Parameters::SetRotatingEularAngles(double alpha, double beta, double gamma){
   this->rotatingEularAngles.SetAlpha(alpha);
   this->rotatingEularAngles.SetBeta(beta);
   this->rotatingEularAngles.SetGamma(gamma);
}

EularAngle Parameters::GetRotatingEularAngles() const{
   return this->rotatingEularAngles;
}

// methods for MOPlot
vector<int>* Parameters::GetIndecesMOPlot() const{
   if(this->indecesMOPlot==NULL){
      stringstream ss;
      ss << this->errorMessageGetIndecesMOPlotNull; 
      throw MolDSException(ss.str());
   }
   return this->indecesMOPlot;
}

void Parameters::AddIndexMOPlot(int moIndex){
   if(this->indecesMOPlot==NULL){
      this->indecesMOPlot = new vector<int>;
   }
   this->indecesMOPlot->push_back(moIndex);
}

bool Parameters::RequiresMOPlot() const{
   return (this->indecesMOPlot!=NULL && 0<this->indecesMOPlot->size());
}

string Parameters::GetFileNamePrefixMOPlot() const{
   return this->fileNamePrefixMOPlot;
}

void Parameters::SetFileNamePrefixMOPlot(string fileNamePrefixMOPlot){
   this->fileNamePrefixMOPlot = fileNamePrefixMOPlot;
}

int* Parameters::GetGridNumberMOPlot() const{
   return (int*)this->gridNumberMOPlot;
}

void Parameters::SetGridNumberMOPlot(int Nx, int Ny, int Nz){
   this->gridNumberMOPlot[XAxis] = Nx;
   this->gridNumberMOPlot[YAxis] = Ny;
   this->gridNumberMOPlot[ZAxis] = Nz;
}

double* Parameters::GetFrameLengthMOPlot() const{
   return (double*)this->frameLengthMOPlot;
}

void Parameters::SetFrameLengthMOPlot(double lx, double ly, double lz){
   this->frameLengthMOPlot[XAxis] = lx;
   this->frameLengthMOPlot[YAxis] = ly;
   this->frameLengthMOPlot[ZAxis] = lz;
}

// methods for HolePlot
vector<int>* Parameters::GetElecIndecesHolePlot() const{
   if(this->elecIndecesHolePlot==NULL){
      stringstream ss;
      ss << this->errorMessageGetIndecesHolePlotNull; 
      throw MolDSException(ss.str());
   }
   return this->elecIndecesHolePlot;
}

void Parameters::AddElecIndexHolePlot(int elecIndex){
   if(this->elecIndecesHolePlot==NULL){
      this->elecIndecesHolePlot = new vector<int>;
   }
   this->elecIndecesHolePlot->push_back(elecIndex);
}

bool Parameters::RequiresHolePlot() const{
   return (this->elecIndecesHolePlot!=NULL && 0<this->elecIndecesHolePlot->size());
}

string Parameters::GetFileNamePrefixHolePlot() const{
   return this->fileNamePrefixHolePlot;
}

void Parameters::SetFileNamePrefixHolePlot(string fileNamePrefixHolePlot){
   this->fileNamePrefixHolePlot = fileNamePrefixHolePlot;
}

int* Parameters::GetGridNumberHolePlot() const{
   return (int*)this->gridNumberHolePlot;
}

void Parameters::SetGridNumberHolePlot(int Nx, int Ny, int Nz){
   this->gridNumberHolePlot[XAxis] = Nx;
   this->gridNumberHolePlot[YAxis] = Ny;
   this->gridNumberHolePlot[ZAxis] = Nz;
}

double* Parameters::GetFrameLengthHolePlot() const{
   return (double*)this->frameLengthHolePlot;
}

void Parameters::SetFrameLengthHolePlot(double lx, double ly, double lz){
   this->frameLengthHolePlot[XAxis] = lx;
   this->frameLengthHolePlot[YAxis] = ly;
   this->frameLengthHolePlot[ZAxis] = lz;
}

// methods for ParticlePlot
vector<int>* Parameters::GetElecIndecesParticlePlot() const{
   if(this->elecIndecesParticlePlot==NULL){
      stringstream ss;
      ss << this->errorMessageGetIndecesParticlePlotNull; 
      throw MolDSException(ss.str());
   }
   return this->elecIndecesParticlePlot;
}

void Parameters::AddElecIndexParticlePlot(int elecIndex){
   if(this->elecIndecesParticlePlot==NULL){
      this->elecIndecesParticlePlot = new vector<int>;
   }
   this->elecIndecesParticlePlot->push_back(elecIndex);
}

bool Parameters::RequiresParticlePlot() const{
   return (this->elecIndecesParticlePlot!=NULL && 0<this->elecIndecesParticlePlot->size());
}

string Parameters::GetFileNamePrefixParticlePlot() const{
   return this->fileNamePrefixParticlePlot;
}

void Parameters::SetFileNamePrefixParticlePlot(string fileNamePrefixParticlePlot){
   this->fileNamePrefixParticlePlot = fileNamePrefixParticlePlot;
}

int* Parameters::GetGridNumberParticlePlot() const{
   return (int*)this->gridNumberParticlePlot;
}

void Parameters::SetGridNumberParticlePlot(int Nx, int Ny, int Nz){
   this->gridNumberParticlePlot[XAxis] = Nx;
   this->gridNumberParticlePlot[YAxis] = Ny;
   this->gridNumberParticlePlot[ZAxis] = Nz;
}

double* Parameters::GetFrameLengthParticlePlot() const{
   return (double*)this->frameLengthParticlePlot;
}

void Parameters::SetFrameLengthParticlePlot(double lx, double ly, double lz){
   this->frameLengthParticlePlot[XAxis] = lx;
   this->frameLengthParticlePlot[YAxis] = ly;
   this->frameLengthParticlePlot[ZAxis] = lz;
}

// methods for CIS
int Parameters::GetActiveOccCIS() const{
   return this->activeOccCIS;
}
   
void Parameters::SetActiveOccCIS(int activeOccCIS){
   this->activeOccCIS = activeOccCIS;
}

int Parameters::GetActiveVirCIS() const{
   return this->activeVirCIS;
}

void Parameters::SetActiveVirCIS(int activeVirCIS){
   this->activeVirCIS = activeVirCIS;
}

int Parameters::GetNumberExcitedStatesCIS() const{
   return this->numberExcitedStatesCIS;
}

void Parameters::SetNumberExcitedStatesCIS(int nStates){
   this->numberExcitedStatesCIS = nStates;
}

bool Parameters::RequiresCIS() const{
   return this->requiresCIS;
}

void Parameters::SetRequiresCIS(bool requiresCIS){
   this->requiresCIS = requiresCIS;
}

bool Parameters::IsDavidsonCIS() const{
   return this->isDavidsonCIS;
}

void Parameters::SetIsDavidsonCIS(bool isDavidsonCIS){
   this->isDavidsonCIS = isDavidsonCIS;
}

int Parameters::GetMaxIterationsCIS() const{
   return this->maxIterationsCIS;
}

void Parameters::SetMaxIterationsCIS(int maxIterationsCIS){
   this->maxIterationsCIS = maxIterationsCIS;
}

int Parameters::GetMaxDimensionsCIS() const{
   return this->maxDimensionsCIS;
}

void Parameters::SetMaxDimensionsCIS(int maxDimensionsCIS){
   this->maxDimensionsCIS = maxDimensionsCIS;
}

double Parameters::GetNormToleranceCIS() const{
   return this->normToleranceCIS;
}

void Parameters::SetNormToleranceCIS(double normToleranceCIS){
   this->normToleranceCIS = normToleranceCIS;
}

int Parameters::GetNumberPrintCoefficientsCIS() const{
   return this->numberPrintCoefficientsCIS;
}
   
void Parameters::SetNumberPrintCoefficientsCIS(int numberPrintCoefficientsCIS){
   this->numberPrintCoefficientsCIS = numberPrintCoefficientsCIS;
}

bool Parameters::RequiresExcitonEnergiesCIS() const{
   return this->requiresExcitonEnergiesCIS;
}

void Parameters::SetRequiresExcitonEnergiesCIS(bool requiresExcitonEnergiesCIS){
   this->requiresExcitonEnergiesCIS = requiresExcitonEnergiesCIS;
}

bool Parameters::RequiresAllTransitionDipoleMomentsCIS() const{
   return this->requiresAllTransitionDipoleMomentsCIS;
}

void Parameters::SetRequiresAllTransitionDipoleMomentsCIS(bool requiresAllTransitionDipoleMomentsCIS){
   this->requiresAllTransitionDipoleMomentsCIS = requiresAllTransitionDipoleMomentsCIS;
}

// methods for Memory
double Parameters::GetLimitHeapMemory() const{
   return this->limitHeapMemory;
}

void Parameters::SetLimitHeapMemory(double limitHeapMemory){
   this->limitHeapMemory = limitHeapMemory;
}

// methods for MD
int Parameters::GetElectronicStateIndexMD() const{
   return this->electronicStateIndexMD;
}

void Parameters::SetElectronicStateIndexMD(int electronicStateIndex){
   this->electronicStateIndexMD = electronicStateIndex;
}

int Parameters::GetTotalStepsMD() const{
   return this->totalStepsMD;
}

void Parameters::SetTotalStepsMD(int totalSteps){
   this->totalStepsMD = totalSteps;
}

double Parameters::GetTimeWidthMD() const{
   return this->timeWidthMD;
}

void Parameters::SetTimeWidthMD(double timeWidth){
   this->timeWidthMD = timeWidth;
}

// methods for MC
int Parameters::GetElectronicStateIndexMC() const{
   return this->electronicStateIndexMC;
}

void Parameters::SetElectronicStateIndexMC(int electronicStateIndex){
   this->electronicStateIndexMC = electronicStateIndex;
}

int Parameters::GetTotalStepsMC() const{
   return this->totalStepsMC;
}

void Parameters::SetTotalStepsMC(int totalSteps){
   this->totalStepsMC = totalSteps;
}

double Parameters::GetTemperatureMC() const{
   return this->temperatureMC;
}

void Parameters::SetTemperatureMC(double temperature){
   this->temperatureMC = temperature;
}

double Parameters::GetStepWidthMC() const{
   return this->stepWidthMC;
}

void Parameters::SetStepWidthMC(double stepWidth){
   this->stepWidthMC = stepWidth;
}

unsigned long Parameters::GetSeedMC() const{
   return this->seedMC;
}

void Parameters::SetSeedMC(unsigned long seed){
   this->seedMC = seed;
}

// methods for RPMD
int Parameters::GetElectronicStateIndexRPMD() const{
   return this->electronicStateIndexRPMD;
}

void Parameters::SetElectronicStateIndexRPMD(int electronicStateIndex){
   this->electronicStateIndexRPMD = electronicStateIndex;
}

int Parameters::GetNumberElectronicStatesRPMD() const{
   return this->numberElectronicStatesRPMD;
}

void Parameters::SetNumberElectronicStatesRPMD(int numberElectronicStates){
   this->numberElectronicStatesRPMD = numberElectronicStates;
}

int Parameters::GetTotalStepsRPMD() const{
   return this->totalStepsRPMD;
}

void Parameters::SetTotalStepsRPMD(int totalSteps){
   this->totalStepsRPMD = totalSteps;
}

double Parameters::GetTemperatureRPMD() const{
   return this->temperatureRPMD;
}

void Parameters::SetTemperatureRPMD(double temperature){
   this->temperatureRPMD = temperature;
}

double Parameters::GetTimeWidthRPMD() const{
   return this->timeWidthRPMD;
}

void Parameters::SetTimeWidthRPMD(double timeWidth){
   this->timeWidthRPMD = timeWidth;
}

int Parameters::GetNumberBeadsRPMD() const{
   return this->numberBeadsRPMD;
}

void Parameters::SetNumberBeadsRPMD(int numberBeads){
   this->numberBeadsRPMD = numberBeads;
}

unsigned long Parameters::GetSeedRPMD() const{
   return this->seedRPMD;
}

void Parameters::SetSeedRPMD(unsigned long seed){
   this->seedRPMD = seed;
}

// methods for NASCO
int Parameters::GetTotalStepsNASCO() const{
   return this->totalStepsNASCO;
}

void Parameters::SetTotalStepsNASCO(int totalSteps){
   this->totalStepsNASCO = totalSteps;
}

int Parameters::GetNumberElectronicStatesNASCO() const{
   return this->numberElectronicStatesNASCO;
}

void Parameters::SetNumberElectronicStatesNASCO(int numberElectronicStates){
   this->numberElectronicStatesNASCO = numberElectronicStates;
}

double Parameters::GetTimeWidthNASCO() const{
   return this->timeWidthNASCO;
}

void Parameters::SetTimeWidthNASCO(double timeWidth){
   this->timeWidthNASCO = timeWidth;
}

unsigned long Parameters::GetSeedNASCO() const{
   return this->seedNASCO;
}

void Parameters::SetSeedNASCO(unsigned long seed){
   this->seedNASCO = seed;
}

// Optimization
OptimizationMethodType Parameters::GetMethodOptimization() const{
   return this->methodOptimization;
}

void Parameters::SetMethodOptimization(OptimizationMethodType method){
   this->methodOptimization = method;
}

int Parameters::GetTotalStepsOptimization() const{
   return this->totalStepsOptimization;
}

void Parameters::SetTotalStepsOptimization(int steps){
   this->totalStepsOptimization = steps;
}

int Parameters::GetElectronicStateIndexOptimization() const{
   return this->electronicStateIndexOptimization;
}

void Parameters::SetElectronicStateIndexOptimization(int electronicStateIndex){
   this->electronicStateIndexOptimization = electronicStateIndex;
}

double Parameters::GetMaxGradientOptimization() const{
   return this->maxGradientOptimization;
}

void Parameters::SetMaxGradientOptimization(double maxGradient){
   this->maxGradientOptimization = maxGradient;
}

double Parameters::GetRmsGradientOptimization() const{
   return this->rmsGradientOptimization;
}

void Parameters::SetRmsGradientOptimization(double rmsGradient){
   this->rmsGradientOptimization = rmsGradient;
}

double Parameters::GetTimeWidthOptimization() const{
   return this->timeWidthOptimization;
}

void Parameters::SetTimeWidthOptimization(double timeWidth){
   this->timeWidthOptimization = timeWidth;
}

// Frequencies
bool Parameters::RequiresFrequencies() const{
   return this->requiresFrequencies;
}

void Parameters::SetRequiresFrequencies(bool requiresFrequencies){
   this->requiresFrequencies = requiresFrequencies;
}

int Parameters::GetElectronicStateIndexFrequencies() const{
   return this->electronicStateIndexFrequencies;
}

void Parameters::SetElectronicStateIndexFrequencies(int electronicStateIndex){
   this->electronicStateIndexFrequencies = electronicStateIndex;
}
}





