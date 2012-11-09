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
#ifndef INCLUDED_PARAMETERS
#define INCLUDED_PARAMETERS

namespace MolDS_base{

// Parameters is singleton
class Parameters: public PrintController, private Uncopyable{
public:
   static Parameters* GetInstance();
   static void DeleteInstance();

   SimulationType GetCurrentSimulation() const;
   void SetCurrentSimulation(SimulationType simulation);
   TheoryType GetCurrentTheory() const;
   void SetCurrentTheory(TheoryType theory);
   // Pysical constants
   double GetEV2AU() const;
   double GetJ2AU() const;
   double GetKcalMolin2AU() const;
   double GetAngstrom2AU() const;
   double GetNm2AU() const;
   double GetKayser2AU() const;
   double GetGMolin2AU() const;
   double GetDegree2Radian() const;
   double GetFs2AU() const;
   double GetBoltzmann() const;
   double GetAvogadro() const;
   double GetDebye2AU() const;
   // SCF
   double GetThresholdSCF() const;
   void   SetThresholdSCF(double thresholdSCF);
   int    GetMaxIterationsSCF() const;
   void   SetMaxIterationsSCF(int maxIterationsSCF);
   double GetDampingThreshSCF() const;
   void   SetDampingThreshSCF(double dampingThreshSCF);
   double GetDampingWeightSCF() const;
   void   SetDampingWeightSCF(double dampingWeightSCF);
   int    GetDiisNumErrorVectSCF() const;
   void   SetDiisNumErrorVectSCF(int diisNumErrorVectSCF);
   double GetDiisStartErrorSCF() const;
   void   SetDiisStartErrorSCF(double diisStartErrorSCF);
   double GetDiisEndErrorSCF() const;
   void   SetDiisEndErrorSCF(double diisEndErrorSCF);
   bool   RequiresVdWSCF() const;
   void   SetRequiresVdWSCF(bool requiresVdWSCF);
   double GetVdWScalingFactorSCF() const;
   void   SetVdWScalingFactorSCF();
   void   SetVdWScalingFactorSCF(double vdWScalingFactorSCF);
   double GetVdWDampingFactorSCF() const;
   void   SetVdWDampingFactorSCF();
   void   SetVdWDampingFactorSCF(double vdWDampingFactorSCF);
   // MOPlot
   std::string       GetFileNamePrefixMOPlot() const;
   void              SetFileNamePrefixMOPlot(std::string fileNamePrefixMOPlot);
   int*              GetGridNumberMOPlot() const;
   void              SetGridNumberMOPlot(int Nx, int Ny, int Nz);
   double*           GetFrameLengthMOPlot() const;
   void              SetFrameLengthMOPlot(double lx, double ly, double lz);
   std::vector<int>* GetIndecesMOPlot() const;
   void              AddIndexMOPlot(int moIndex);
   bool              RequiresMOPlot() const;
   // HoleEPlot
   std::string       GetFileNamePrefixHolePlot() const;
   void              SetFileNamePrefixHolePlot(std::string fileNamePrefixHolePlot);
   int*              GetGridNumberHolePlot() const;
   void              SetGridNumberHolePlot(int Nx, int Ny, int Nz);
   double*           GetFrameLengthHolePlot() const;
   void              SetFrameLengthHolePlot(double lx, double ly, double lz);
   std::vector<int>* GetElecIndecesHolePlot() const;
   void              AddElecIndexHolePlot(int elecIndex);
   bool              RequiresHolePlot() const;
   // ParticlePlot
   std::string       GetFileNamePrefixParticlePlot() const;
   void              SetFileNamePrefixParticlePlot(std::string fileNamePrefixParticlePlot);
   int*              GetGridNumberParticlePlot() const;
   void              SetGridNumberParticlePlot(int Nx, int Ny, int Nz);
   double*           GetFrameLengthParticlePlot() const;
   void              SetFrameLengthParticlePlot(double lx, double ly, double lz);
   std::vector<int>* GetElecIndecesParticlePlot() const;
   void              AddElecIndexParticlePlot(int elecIndex);
   bool              RequiresParticlePlot() const;
   // Translation
   void    SetTranslatingDifference(double x, double y, double z);
   double* GetTranslatingDifference() const;
   // Principal axes
   void    SetInertiaTensorOrigin(double x, double y, double z);
   double* GetInertiaTensorOrigin() const;
   // Rotation
   void         SetRotatingOrigin(double x, double y, double z);
   double*      GetRotatingOrigin() const;
   void         SetRotatingType(RotatingType rotatingType);
   RotatingType GetRotatingType() const;
   void         SetRotatingAxis(double x, double y, double z);
   double*      GetRotatingAxis() const;
   void         SetRotatingAngle(double rotatingAngle);
   double       GetRotatingAngle() const;
   void         SetRotatingEularAngles(double alpha, double beta, double gamma);
   EularAngle   GetRotatingEularAngles() const;
   // CIS
   int    GetActiveOccCIS() const;
   void   SetActiveOccCIS(int activeOccCIS);
   int    GetActiveVirCIS() const;
   void   SetActiveVirCIS(int activeOccCIS);
   int    GetNumberExcitedStatesCIS() const;
   void   SetNumberExcitedStatesCIS(int nStates);
   bool   RequiresCIS() const;
   void   SetRequiresCIS(bool requiresCIS);
   bool   IsDavidsonCIS() const;
   void   SetIsDavidsonCIS(bool isDavidsonCIS);
   int    GetMaxIterationsCIS() const;
   void   SetMaxIterationsCIS(int maxIterationsCIS);
   int    GetMaxDimensionsCIS() const;
   void   SetMaxDimensionsCIS(int maxDimensionsCIS);
   double GetNormToleranceCIS() const;
   void   SetNormToleranceCIS(double normToleranceCIS);
   int    GetNumberPrintCoefficientsCIS() const;
   void   SetNumberPrintCoefficientsCIS(int numberPrintCoefficientsCIS);
   bool   RequiresExcitonEnergiesCIS() const;
   void   SetRequiresExcitonEnergiesCIS(bool requiresExcitonEnergiesCIS);
   bool   RequiresAllTransitionDipoleMomentsCIS() const;
   void   SetRequiresAllTransitionDipoleMomentsCIS(bool requiresAllTransitionDipoleMomentsCIS);
   // Memory
   double GetLimitHeapMemory() const;
   void   SetLimitHeapMemory(double limitHeap);
   // MD
   int    GetElectronicStateIndexMD() const;
   void   SetElectronicStateIndexMD(int electronicStateIndex);
   int    GetTotalStepsMD() const;
   void   SetTotalStepsMD(int totalSteps);
   double GetTimeWidthMD() const;
   void   SetTimeWidthMD(double timeWidth);
   // MC
   int           GetElectronicStateIndexMC() const;
   void          SetElectronicStateIndexMC(int electronicStateIndex);
   int           GetTotalStepsMC() const;
   void          SetTotalStepsMC(int totalSteps);
   double        GetTemperatureMC() const;
   void          SetTemperatureMC(double temperature);
   double        GetStepWidthMC() const;
   void          SetStepWidthMC(double stepWidth);
   unsigned long GetSeedMC() const;
   void          SetSeedMC(unsigned long seed);
   // RPMD
   int           GetElectronicStateIndexRPMD() const;
   void          SetElectronicStateIndexRPMD(int electronicStateIndex);
   int           GetNumberElectronicStatesRPMD() const;
   void          SetNumberElectronicStatesRPMD(int NumberElectronicStates);
   int           GetTotalStepsRPMD() const;
   void          SetTotalStepsRPMD(int totalSteps);
   double        GetTemperatureRPMD() const;
   void          SetTemperatureRPMD(double temperature);
   double        GetTimeWidthRPMD() const;
   void          SetTimeWidthRPMD(double stepWidth);
   int           GetNumberBeadsRPMD() const;
   void          SetNumberBeadsRPMD(int numberBeads);
   unsigned long GetSeedRPMD() const;
   void          SetSeedRPMD(unsigned long seed);
   // NASCO
   int           GetTotalStepsNASCO() const;
   void          SetTotalStepsNASCO(int totalSteps);
   int           GetNumberElectronicStatesNASCO() const;
   void          SetNumberElectronicStatesNASCO(int NumberElectronicStates);
   int           GetInitialElectronicStateNASCO() const;
   void          SetInitialElectronicStateNASCO(int initialElectronicState);
   unsigned long GetSeedNASCO() const;
   void          SetSeedNASCO(unsigned long seed);
   double        GetTimeWidthNASCO() const;
   void          SetTimeWidthNASCO(double stepWidth);
   // Optimization
   OptimizationMethodType GetMethodOptimization() const;
   void                   SetMethodOptimization(OptimizationMethodType method);
   int                    GetTotalStepsOptimization() const;
   void                   SetTotalStepsOptimization(int totalSteps);
   int                    GetElectronicStateIndexOptimization() const;
   void                   SetElectronicStateIndexOptimization(int electronicStateIndex);
   double                 GetMaxGradientOptimization() const;
   void                   SetMaxGradientOptimization(double maxGradient);
   double                 GetRmsGradientOptimization() const;
   void                   SetRmsGradientOptimization(double rmsGradient);
   double                 GetTimeWidthOptimization() const;
   void                   SetTimeWidthOptimization(double timeWidth);
   // Frequencies 
   bool RequiresFrequencies() const;
   void SetRequiresFrequencies(bool requiresFrequencies);
   int  GetElectronicStateIndexFrequencies() const;
   void SetElectronicStateIndexFrequencies(int electronicStateIndex);

private:
   static Parameters* parameters;
   Parameters();
   ~Parameters();
   std::string errorMessageGetIndecesMOPlotNull;
   std::string errorMessageGetIndecesHolePlotNull;
   std::string errorMessageGetIndecesParticlePlotNull;
   SimulationType currentSimulation;
   TheoryType currentTheory;
   // Physical constants
   static const double eV2AU;
   static const double j2AU;
   static const double kcalMolin2AU;
   static const double angstrom2AU;
   static const double nm2AU;
   static const double kayser2AU;
   static const double gMolin2AU;
   static const double degree2Radian;
   static const double fs2AU;
   static const double boltzmann;
   static const double avogadro;
   static const double debye2AU;
   // SCF
   double thresholdSCF;
   int    maxIterationsSCF;
   double dampingThreshSCF;
   double dampingWeightSCF;
   int    diisNumErrorVectSCF;
   double diisStartErrorSCF;
   double diisEndErrorSCF;
   bool   requiresVdWSCF;
   double vdWScalingFactorSCF;
   double vdWDampingFactorSCF;
   static const double vdWScalingFactorSCFPM3DAM1D;
   static const double vdWDampingFactorSCFPM3DAM1D;
   // MOPlot
   std::string       fileNamePrefixMOPlot;
   int               gridNumberMOPlot[CartesianType_end];
   double            frameLengthMOPlot[CartesianType_end];
   std::vector<int>* indecesMOPlot;
   // HolePlot
   std::string       fileNamePrefixHolePlot;
   int               gridNumberHolePlot[CartesianType_end];
   double            frameLengthHolePlot[CartesianType_end];
   std::vector<int>* elecIndecesHolePlot;
   // ParticlePlot
   std::string       fileNamePrefixParticlePlot;
   int               gridNumberParticlePlot[CartesianType_end];
   double            frameLengthParticlePlot[CartesianType_end];
   std::vector<int>* elecIndecesParticlePlot;
   // Translation
   double translatingDifference[3];
   // Principal axes
   double* inertiaTensorOrigin;
   // Rotation
   double*      rotatingOrigin;
   double       rotatingAxis[3];
   double       rotatingAngle;
   RotatingType rotatingType;
   EularAngle   rotatingEularAngles;
   // CIS
   int    activeOccCIS;
   int    activeVirCIS;
   int    numberExcitedStatesCIS;
   int    maxIterationsCIS;
   int    maxDimensionsCIS;
   double normToleranceCIS;
   bool   requiresCIS;
   bool   isDavidsonCIS;
   int    numberPrintCoefficientsCIS;
   bool   requiresExcitonEnergiesCIS;
   bool   requiresAllTransitionDipoleMomentsCIS;
   // Memory
   double limitHeapMemory;
   // MD
   int    electronicStateIndexMD;
   int    totalStepsMD;
   double timeWidthMD;
   // MC
   int           electronicStateIndexMC;
   int           totalStepsMC;
   double        temperatureMC;
   double        stepWidthMC;
   unsigned long seedMC;
   // RPMD
   int           electronicStateIndexRPMD;
   int           numberElectronicStatesRPMD;
   int           totalStepsRPMD;
   double        temperatureRPMD;
   double        timeWidthRPMD;
   int           numberBeadsRPMD;
   unsigned long seedRPMD;
   // NASCO
   int           totalStepsNASCO;
   int           numberElectronicStatesNASCO;
   int           initialElectronicStateNASCO;
   double        timeWidthNASCO;
   unsigned long seedNASCO;
   // Optimization
   OptimizationMethodType methodOptimization;
   int                    totalStepsOptimization;
   int                    electronicStateIndexOptimization;
   double                 maxGradientOptimization;
   double                 rmsGradientOptimization;
   double                 timeWidthOptimization;
   // Frequencies
   bool requiresFrequencies;
   int  electronicStateIndexFrequencies;
   // Other
   void SetDefaultValues();
   void SetMessages();
};

}
#endif





