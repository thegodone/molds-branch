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
#ifndef INCLUDED_PARAMETERS
#define INCLUDED_PARAMETERS

namespace MolDS_base{

// Parameters is singleton
class Parameters: public PrintController, private Uncopyable{
public:
   static Parameters* GetInstance();
   static void DeleteInstance();
   inline SimulationType GetCurrentSimulation() const{return this->currentSimulation;}
   inline void           SetCurrentSimulation(SimulationType simulation) {this->currentSimulation = simulation;}
   inline TheoryType GetCurrentTheory() const{return this->currentTheory;}
   inline void       SetCurrentTheory(TheoryType theory){this->currentTheory = theory;}
   // Pysical constants
   inline double GetEV2AU() const        {return this->eV2AU;}
   inline double GetJ2AU() const         {return this->j2AU;}
   inline double GetKcalMolin2AU() const {return this->kcalMolin2AU;}
   inline double GetAngstrom2AU() const  {return this->angstrom2AU;}
   inline double GetNm2AU() const        {return this->nm2AU;}
   inline double GetKayser2AU() const    {return this->kayser2AU;}
   inline double GetGMolin2AU() const    {return this->gMolin2AU;}
   inline double GetDegree2Radian() const{return this->degree2Radian;}
   inline double GetFs2AU() const        {return this->fs2AU;}
   inline double GetBoltzmann() const    {return this->boltzmann;}
   inline double GetAvogadro() const     {return this->avogadro;}
   inline double GetDebye2AU() const     {return this->debye2AU;}
   // SCF
   inline double GetThresholdSCF() const               {return this->thresholdSCF;}
   inline void   SetThresholdSCF(double threshold)     {this->thresholdSCF = threshold;}
   inline int    GetMaxIterationsSCF() const           {return this->maxIterationsSCF;}
   inline void   SetMaxIterationsSCF(int maxIter)      {this->maxIterationsSCF = maxIter;}
   inline double GetDampingThreshSCF() const           {return this->dampingThreshSCF;}
   inline void   SetDampingThreshSCF(double dThresh)   {this->dampingThreshSCF = dThresh;}
   inline double GetDampingWeightSCF() const           {return this->dampingWeightSCF;}
   inline void   SetDampingWeightSCF(double dWeight)   {this->dampingWeightSCF = dWeight;}
   inline int    GetDiisNumErrorVectSCF() const        {return this->diisNumErrorVectSCF;}
   inline void   SetDiisNumErrorVectSCF(int numEVect)  {this->diisNumErrorVectSCF = numEVect;}
   inline double GetDiisStartErrorSCF() const          {return this->diisStartErrorSCF;}
   inline void   SetDiisStartErrorSCF(double sError)   {this->diisStartErrorSCF = sError;}
   inline double GetDiisEndErrorSCF() const            {return this->diisEndErrorSCF;}
   inline void   SetDiisEndErrorSCF(double eError)     {this->diisEndErrorSCF = eError;}
   inline bool   RequiresVdWSCF() const                {return this->requiresVdWSCF;}
   inline void   SetRequiresVdWSCF(bool requires)      {this->requiresVdWSCF = requires;}
   inline double GetVdWScalingFactorSCF() const        {return this->vdWScalingFactorSCF;}
   inline void   SetVdWScalingFactorSCF()              {this->vdWScalingFactorSCF = this->vdWScalingFactorSCFPM3DAM1D;}
   inline void   SetVdWScalingFactorSCF(double vdWScal){this->vdWScalingFactorSCF = vdWScal;}
   inline double GetVdWDampingFactorSCF() const        {return this->vdWDampingFactorSCF;}
   inline void   SetVdWDampingFactorSCF()              {this->vdWDampingFactorSCF = this->vdWDampingFactorSCFPM3DAM1D;}
   inline void   SetVdWDampingFactorSCF(double vdWDamp){this->vdWDampingFactorSCF = vdWDamp;}
   // MOPlot
   inline bool          RequiresMOPlot() const                     {return (this->indecesMOPlot!=NULL && 0<this->indecesMOPlot->size());}
   inline std::string   GetFileNamePrefixMOPlot() const            {return this->fileNamePrefixMOPlot;}
   inline void          SetFileNamePrefixMOPlot(std::string prefix){this->fileNamePrefixMOPlot = prefix;}
   inline const int*    GetGridNumberMOPlot() const                {return (int*)this->gridNumberMOPlot;}
   void                 SetGridNumberMOPlot(int Nx, int Ny, int Nz);
   inline const double* GetFrameLengthMOPlot() const               {return (double*)this->frameLengthMOPlot;}
   void                 SetFrameLengthMOPlot(double lx, double ly, double lz);
   std::vector<int>*    GetIndecesMOPlot() const;
   void                 AddIndexMOPlot(int moIndex);
   // HoleEPlot
   inline bool          RequiresHolePlot() const                     {return (this->elecIndecesHolePlot!=NULL && 0<this->elecIndecesHolePlot->size());}
   inline std::string   GetFileNamePrefixHolePlot() const            {return this->fileNamePrefixHolePlot;}
   void                 SetFileNamePrefixHolePlot(std::string prefix){this->fileNamePrefixHolePlot = prefix;}
   inline const int*    GetGridNumberHolePlot() const                {return (int*)this->gridNumberHolePlot;}
   void                 SetGridNumberHolePlot(int Nx, int Ny, int Nz);
   inline const double* GetFrameLengthHolePlot() const               {return (double*)this->frameLengthHolePlot;}
   void                 SetFrameLengthHolePlot(double lx, double ly, double lz);
   std::vector<int>*    GetElecIndecesHolePlot() const;
   void                 AddElecIndexHolePlot(int elecIndex);
   // ParticlePlot
   inline bool             RequiresParticlePlot() const                     {return (this->elecIndecesParticlePlot!=NULL && 0<this->elecIndecesParticlePlot->size());}
   const std::vector<int>* GetElecIndecesParticlePlot() const;
   void                    AddElecIndexParticlePlot(int elecIndex);
   inline std::string      GetFileNamePrefixParticlePlot() const            {return this->fileNamePrefixParticlePlot;}
   inline void             SetFileNamePrefixParticlePlot(std::string prefix){this->fileNamePrefixParticlePlot = prefix;}
   inline const int*       GetGridNumberParticlePlot() const                {return (int*)this->gridNumberParticlePlot;}
   void                    SetGridNumberParticlePlot(int Nx, int Ny, int Nz);
   inline const double*    GetFrameLengthParticlePlot() const               {return (double*)this->frameLengthParticlePlot;}
   void                    SetFrameLengthParticlePlot(double lx, double ly, double lz);
   // Translation
   inline const double* GetTranslatingDifference() const{return (double*)this->translatingDifference;}
   void                 SetTranslatingDifference(double x, double y, double z);
   // Principal axes
   inline const double* GetInertiaTensorOrigin() const{return (double*)this->inertiaTensorOrigin;}
   void                 SetInertiaTensorOrigin(double x, double y, double z);
   // Rotation
   inline const double* GetRotatingOrigin() const                 {return (double*)this->rotatingOrigin;}
   void                 SetRotatingOrigin(double x, double y, double z);
   inline void          SetRotatingType(RotatingType rotatingType){this->rotatingType = rotatingType;}
   inline RotatingType  GetRotatingType() const                   {return this->rotatingType;}
   void                 SetRotatingAxis(double x, double y, double z);
   inline const double* GetRotatingAxis() const                   {return (double*)this->rotatingAxis;}
   inline double        GetRotatingAngle() const                  {return this->rotatingAngle;}
   inline void          SetRotatingAngle(double rotatingAngle)    {this->rotatingAngle = rotatingAngle;}
   void                 SetRotatingEularAngles(double alpha, double beta, double gamma);
   inline EularAngle    GetRotatingEularAngles() const            {return this->rotatingEularAngles;}
   // CIS
   inline int        GetActiveOccCIS() const                                {return this->activeOccCIS;}
   inline void       SetActiveOccCIS(int activeOccCIS)                      {this->activeOccCIS = activeOccCIS;}
   inline int        GetActiveVirCIS() const                                {return this->activeVirCIS;}
   inline void       SetActiveVirCIS(int activeVirCIS)                      {this->activeVirCIS = activeVirCIS;}
   inline int        GetNumberExcitedStatesCIS() const                      {return this->numberExcitedStatesCIS;}
   inline void       SetNumberExcitedStatesCIS(int nStates)                 {this->numberExcitedStatesCIS = nStates;}
   inline bool       RequiresCIS() const                                    {return this->requiresCIS;}
   inline void       SetRequiresCIS(bool requiresCIS)                       {this->requiresCIS = requiresCIS;}
   inline bool       IsDavidsonCIS() const                                  {return this->isDavidsonCIS;}
   inline void       SetIsDavidsonCIS(bool isDavidsonCIS)                   {this->isDavidsonCIS = isDavidsonCIS;}
   inline int        GetMaxIterationsCIS() const                            {return this->maxIterationsCIS;}
   inline void       SetMaxIterationsCIS(int maxIter)                       {this->maxIterationsCIS = maxIter;}
   inline int        GetMaxDimensionsCIS() const                            {return this->maxDimensionsCIS;}
   inline void       SetMaxDimensionsCIS(int maxDim)                        {this->maxDimensionsCIS = maxDim;}
   inline double     GetNormToleranceCIS() const                            {return this->normToleranceCIS;}
   inline void       SetNormToleranceCIS(double normTol)                    {this->normToleranceCIS = normTol;}
   inline int        GetNumberPrintCoefficientsCIS() const                  {return this->numberPrintCoefficientsCIS;}
   inline void       SetNumberPrintCoefficientsCIS(int number)              {this->numberPrintCoefficientsCIS = number;}
   inline bool       RequiresExcitonEnergiesCIS() const                     {return this->requiresExcitonEnergiesCIS;}
   inline void       SetRequiresExcitonEnergiesCIS(bool requires)           {this->requiresExcitonEnergiesCIS = requires;}
   inline bool       RequiresAllTransitionDipoleMomentsCIS() const          {return this->requiresAllTransitionDipoleMomentsCIS;}
   inline void       SetRequiresAllTransitionDipoleMomentsCIS(bool requires){this->requiresAllTransitionDipoleMomentsCIS = requires;}
   std::vector<int>* GetElectronicStateIndecesMullikenCIS() const;
   void              AddElectronicStateIndexMullikenCIS(int electronicStateIndex);
   bool              RequiresMullikenCIS() const;
   inline bool       RequiresUnpairedPopCIS() const                         {return this->requiresUnpairedPopCIS;}
   inline void       SetRequiresUnpairedPopCIS(bool requires)               {this->requiresUnpairedPopCIS = requires;}
   // Memory
   double GetLimitHeapMemory() const          {return this->limitHeapMemory;}
   void   SetLimitHeapMemory(double limitHeap){this->limitHeapMemory = limitHeap;}
   // MD
   int    GetElectronicStateIndexMD() const   {return this->electronicStateIndexMD;}
   void   SetElectronicStateIndexMD(int index){this->electronicStateIndexMD = index;}
   int    GetTotalStepsMD() const             {return this->totalStepsMD;}
   void   SetTotalStepsMD(int steps)          {this->totalStepsMD = steps;}
   double GetTimeWidthMD() const              {return this->timeWidthMD;}
   void   SetTimeWidthMD(double dt)           {this->timeWidthMD = dt;}
   // MC
   int           GetElectronicStateIndexMC() const{return this->electronicStateIndexMC;}
   void          SetElectronicStateIndexMC(int i) {this->electronicStateIndexMC = i;}
   int           GetTotalStepsMC() const          {return this->totalStepsMC;}
   void          SetTotalStepsMC(int steps)       {this->totalStepsMC = steps;}
   double        GetTemperatureMC() const         {return this->temperatureMC;}
   void          SetTemperatureMC(double t)       {this->temperatureMC = t;}
   double        GetStepWidthMC() const           {return this->stepWidthMC;}
   void          SetStepWidthMC(double dr)        {this->stepWidthMC = dr;}
   unsigned long GetSeedMC() const                {return this->seedMC;}
   void          SetSeedMC(unsigned long seed)    {this->seedMC = seed;}
   // RPMD
   int           GetElectronicStateIndexRPMD() const  {return this->electronicStateIndexRPMD;}
   void          SetElectronicStateIndexRPMD(int i)   {this->electronicStateIndexRPMD = i;}
   int           GetNumberElectronicStatesRPMD() const{return this->numberElectronicStatesRPMD;}
   void          SetNumberElectronicStatesRPMD(int n) {this->numberElectronicStatesRPMD = n;}
   int           GetTotalStepsRPMD() const            {return this->totalStepsRPMD;}
   void          SetTotalStepsRPMD(int steps)         {this->totalStepsRPMD = steps;}
   double        GetTemperatureRPMD() const           {return this->temperatureRPMD;}
   void          SetTemperatureRPMD(double t)         {this->temperatureRPMD = t;}
   double        GetTimeWidthRPMD() const             {return this->timeWidthRPMD;}
   void          SetTimeWidthRPMD(double dr)          {this->timeWidthRPMD = dr;}
   int           GetNumberBeadsRPMD() const           {return this->numberBeadsRPMD;}
   void          SetNumberBeadsRPMD(int b)            {this->numberBeadsRPMD = b;}
   unsigned long GetSeedRPMD() const                  {return this->seedRPMD;}
   void          SetSeedRPMD(unsigned long seed)      {this->seedRPMD = seed;}
   // NASCO
   int           GetTotalStepsNASCO() const            {return this->totalStepsNASCO;}
   void          SetTotalStepsNASCO(int steps)         {this->totalStepsNASCO = steps;}
   int           GetNumberElectronicStatesNASCO() const{return this->numberElectronicStatesNASCO;}
   void          SetNumberElectronicStatesNASCO(int n) {this->numberElectronicStatesNASCO = n;}
   int           GetInitialElectronicStateNASCO() const{return this->initialElectronicStateNASCO;}
   void          SetInitialElectronicStateNASCO(int i) {this->initialElectronicStateNASCO = i;}
   unsigned long GetSeedNASCO() const                  {return this->seedNASCO;}
   void          SetSeedNASCO(unsigned long seed)      {this->seedNASCO = seed;}
   double        GetTimeWidthNASCO() const             {return this->timeWidthNASCO;}
   void          SetTimeWidthNASCO(double dt)          {this->timeWidthNASCO = dt;}
   // Optimization
   OptimizationMethodType GetMethodOptimization() const                  {return this->methodOptimization;}
   void                   SetMethodOptimization(OptimizationMethodType m){this->methodOptimization = m;}
   int                    GetTotalStepsOptimization() const              {return this->totalStepsOptimization;}
   void                   SetTotalStepsOptimization(int steps)           {this->totalStepsOptimization = steps;}
   int                    GetElectronicStateIndexOptimization() const    {return this->electronicStateIndexOptimization;}
   void                   SetElectronicStateIndexOptimization(int i)     {this->electronicStateIndexOptimization = i;}
   double                 GetMaxGradientOptimization() const             {return this->maxGradientOptimization;}
   void                   SetMaxGradientOptimization(double m)           {this->maxGradientOptimization = m;}
   double                 GetRmsGradientOptimization() const             {return this->rmsGradientOptimization;}
   void                   SetRmsGradientOptimization(double r)           {this->rmsGradientOptimization = r;}
   double                 GetTimeWidthOptimization() const               {return this->timeWidthOptimization;}
   void                   SetTimeWidthOptimization(double dt)            {this->timeWidthOptimization = dt;}
   double                 GetInitialTrustRadiusOptimization() const      {return this->initialTrustRadiusOptimization;}
   void                   SetInitialTrustRadiusOptimization(double r)    {this->initialTrustRadiusOptimization = r;}
   double                 GetMaxNormStepOptimization() const             {return this->maxNormStepOptimization;}
   void                   SetMaxNormStepOptimization(double n)           {this->maxNormStepOptimization = n;}
   // Frequencies 
   bool RequiresFrequencies() const               {return this->requiresFrequencies;}
   void SetRequiresFrequencies(bool b)            {this->requiresFrequencies = b;}
   int  GetElectronicStateIndexFrequencies() const{return this->electronicStateIndexFrequencies;}
   void SetElectronicStateIndexFrequencies(int i) {this->electronicStateIndexFrequencies = i;}

private:
   static Parameters* parameters;
   Parameters();
   ~Parameters();
   std::string errorMessageGetIndecesMOPlotNull;
   std::string errorMessageGetIndecesHolePlotNull;
   std::string errorMessageGetIndecesParticlePlotNull;
   std::string errorMessageGetElectronicStateIndecesMullikenCISNull;
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
   int               activeOccCIS;
   int               activeVirCIS;
   int               numberExcitedStatesCIS;
   int               maxIterationsCIS;
   int               maxDimensionsCIS;
   double            normToleranceCIS;
   bool              requiresCIS;
   bool              isDavidsonCIS;
   int               numberPrintCoefficientsCIS;
   bool              requiresExcitonEnergiesCIS;
   bool              requiresAllTransitionDipoleMomentsCIS;
   std::vector<int>* electronicStateIndecesMullikenCIS;
   bool              requiresUnpairedPopCIS;
   // Memory
   double limitHeapMemory; // in [MB]
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
   double                 initialTrustRadiusOptimization;
   double                 maxNormStepOptimization;
   // Frequencies
   bool requiresFrequencies;
   int  electronicStateIndexFrequencies;
   // Other
   void SetDefaultValues();
   void SetMessages();
};

}
#endif





