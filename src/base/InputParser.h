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
#ifndef INCLUDED_INPUT_PARSER
#define INCLUDED_INPUT_PARSER
namespace MolDS_base{

// InputParser is singleton
class InputParser: public PrintController, private Uncopyable{
public:
   static InputParser* GetInstance();
   static void DeleteInstance();
   void Parse(Molecule* molecule) const;
private:
   static InputParser* inputParser;
   InputParser();
   ~InputParser();
   void SetMessages();
   std::string errorMessageNonValidExcitedStatesMD;
   std::string errorMessageNonValidExcitedStatesMC;
   std::string errorMessageNonValidExcitedStatesRPMD;
   std::string errorMessageNonValidExcitedStatesOptimization;
   std::string errorMessageElecState;
   std::string errorMessageTheory;
   std::string errorMessageNumberExcitedStateCIS;
   std::string messageStartParseInput;
   std::string messageDoneParseInput;
   std::string messageTotalNumberAOs;
   std::string messageTotalNumberAtoms;
   std::string messageTotalNumberValenceElectrons;
   std::string messageInputTerms;
   // SCF
   std::string messageScfConditions;
   std::string messageScfMaxIterations;
   std::string messageScfRmsDensity;
   std::string messageScfDampingThresh;
   std::string messageScfDampingWeight;
   std::string messageScfDiisNumErrorVect;
   std::string messageScfDiisStartError;
   std::string messageScfDiisEndError;
   std::string messageScfVdW;
   std::string messageScfVdWScalingFactor;
   std::string messageScfVdWDampingFactor;
   // CIS
   std::string messageCisConditions;
   std::string messageCisNumberActiveOcc;
   std::string messageCisNumberActiveVir;
   std::string messageCisNumberExcitedStates;
   std::string messageCisDavidson;
   std::string messageCisNormTolerance;
   std::string messageCisMaxIterations;
   std::string messageCisMaxDimensions;
   std::string messageCisExcitonEnergies;
   std::string messageCisAllTransitionDipoleMoments;
   std::string messageCisNumPrintCoefficients;
   // Memory
   std::string messageMemoryConditions;
   std::string messageMemoryLimitHeap;
   std::string messageMemoryMB;
   // MD
   std::string messageMdConditions;
   std::string messageMdTotalSteps;
   std::string messageMdElecState;
   std::string messageMdTimeWidth;
   // MC
   std::string messageMcConditions;
   std::string messageMcTotalSteps;
   std::string messageMcElecState;
   std::string messageMcStepWidth;
   std::string messageMcTemperature;
   std::string messageMcSeed;
   // RPMD
   std::string messageRpmdConditions;
   std::string messageRpmdTotalSteps;
   std::string messageRpmdElecState;
   std::string messageRpmdNumElecStates;
   std::string messageRpmdTimeWidth;
   std::string messageRpmdTemperature;
   std::string messageRpmdNumBeads;
   std::string messageRpmdSeed;
   // Optimization
   std::string messageOptimizationConditions;
   std::string messageOptimizationMethod;
   std::string messageOptimizationTotalSteps;
   std::string messageOptimizationElecState;
   std::string messageOptimizationTimeWidth;
   std::string messageOptimizationRmsGradient;
   std::string messageOptimizationMaxGradient;
   // MOPlot
   std::string messageMOPlotConditions;
   std::string messageMOPlotIndex;
   std::string messageMOPlotGridNumber;
   std::string messageMOPlotFrameLength;
   std::string messageMOPlotFilePrefix;
   // HolePlot
   std::string messageHolePlotConditions;
   std::string messageHolePlotElecIndex;
   std::string messageHolePlotGridNumber;
   std::string messageHolePlotFrameLength;
   std::string messageHolePlotFilePrefix;
   // ParticlePlot
   std::string messageParticlePlotConditions;
   std::string messageParticlePlotElecIndex;
   std::string messageParticlePlotGridNumber;
   std::string messageParticlePlotFrameLength;
   std::string messageParticlePlotFilePrefix;
   // unit 
   std::string messageFs;
   std::string messageK;
   std::string messageAngst;
   // others
   std::string stringYES;
   std::string stringNO;
   std::string stringSpace;
   std::string stringTab;
   std::string stringCommentOut;
   // Theory
   std::string stringTheory;
   std::string stringTheoryEnd;
   std::string stringTheoryCNDO2;
   std::string stringTheoryINDO;
   std::string stringTheoryZINDOS;
   std::string stringTheoryMNDO;
   std::string stringTheoryAM1;
   std::string stringTheoryAM1D;
   std::string stringTheoryPM3;
   std::string stringTheoryPM3D;
   std::string stringTheoryPM3PDDG;
   // geometry
   std::string stringGeometry;
   std::string stringGeometryEnd;
   // SCF
   std::string stringScf;
   std::string stringScfEnd;
   std::string stringScfMaxIter;
   std::string stringScfRmsDensity;
   std::string stringScfDampingThresh;
   std::string stringScfDampingWeight;
   std::string stringScfDiisNumErrorVect;
   std::string stringScfDiisStartError;
   std::string stringScfDiisEndError;
   std::string stringScfVdW;
   std::string stringScfVdWScalingFactor;
   std::string stringScfVdWDampingFactor;
   // MOPlot
   std::string stringMO;
   std::string stringMOPlot;
   std::string stringMOPlotEnd;
   std::string stringMOPlotGridNumber;
   std::string stringMOPlotFrameLength;
   std::string stringMOPlotFilePrefix;
   // HolePlot
   std::string stringHolePlot;
   std::string stringHolePlotEnd;
   std::string stringHolePlotElecIndex;
   std::string stringHolePlotGridNumber;
   std::string stringHolePlotFrameLength;
   std::string stringHolePlotFilePrefix;
   // ParticlePlot
   std::string stringParticlePlot;
   std::string stringParticlePlotEnd;
   std::string stringParticlePlotElecIndex;
   std::string stringParticlePlotGridNumber;
   std::string stringParticlePlotFrameLength;
   std::string stringParticlePlotFilePrefix;
   // Principal axes
   std::string stringInertiaTensor;
   std::string stringInertiaTensorEnd;
   std::string stringInertiaTensorOrigin;
   // Rotation
   std::string stringRotate;
   std::string stringRotateEnd;
   std::string stringRotatingOrigin;
   std::string stringRotatingAxis;
   std::string stringRotatingAngle;
   std::string stringRotatingAngles;
   std::string stringRotatingType;
   std::string stringRotatingTypeAxis;
   std::string stringRotatingTypeEularAngle;
   // Translation
   std::string stringTranslate;
   std::string stringTranslateEnd;
   std::string stringTranslatingDifference;
   // CIS
   std::string stringCIS;
   std::string stringCISEnd;
   std::string stringCISActiveOcc;
   std::string stringCISActiveVir;
   std::string stringCISNStates;
   std::string stringCISDavidson;
   std::string stringCISMaxIter;
   std::string stringCISMaxDimensions;
   std::string stringCISNormTolerance;
   std::string stringCISExcitonEnergies;
   std::string stringCISAllTransitionDipoleMoments;
   std::string stringCISNumPrintCoefficients;
   // Memory
   std::string stringMemory;
   std::string stringMemoryEnd;
   std::string stringMemoryLimitHeap;
   // MD
   std::string stringMD;
   std::string stringMDEnd;
   std::string stringMDTotalSteps;
   std::string stringMDElecState;
   std::string stringMDTimeWidth;
   // MC
   std::string stringMC;
   std::string stringMCEnd;
   std::string stringMCTotalSteps;
   std::string stringMCElecState;
   std::string stringMCStepWidth;
   std::string stringMCTemperature;
   std::string stringMCSeed;
   // RPMD
   std::string stringRPMD;
   std::string stringRPMDEnd;
   std::string stringRPMDTotalSteps;
   std::string stringRPMDElecState;
   std::string stringRPMDNumElecStates;
   std::string stringRPMDTimeWidth;
   std::string stringRPMDTemperature;
   std::string stringRPMDNumBeads;
   std::string stringRPMDSeed;
   // Optimization
   std::string stringOptimization;
   std::string stringOptimizationEnd;
   std::string stringOptimizationMethod;
   std::string stringOptimizationBFGS;
   std::string stringOptimizationConjugateGradient;
   std::string stringOptimizationSteepestDescent;
   std::string stringOptimizationTotalSteps;
   std::string stringOptimizationElecState;
   std::string stringOptimizationMaxGradient;
   std::string stringOptimizationRmsGradient;
   std::string stringOptimizationTimeWidth;
   void CalcMolecularBasics(Molecule* molecule) const;
   void ValidateVdWConditions() const;
   void ValidateCisConditions(const Molecule& molecule) const;
   void ValidateMdConditions(const Molecule& molecule) const;
   void ValidateMcConditions(const Molecule& molecule) const;
   void ValidateRpmdConditions(const Molecule& molecule) const;
   void ValidateOptimizationConditions(const Molecule& molecule) const;
   void OutputMolecularBasics(Molecule* molecule) const;
   void OutputScfConditions() const;
   void OutputMemoryConditions() const;
   void OutputCisConditions() const;
   void OutputMdConditions() const;
   void OutputMcConditions() const;
   void OutputRpmdConditions() const;
   void OutputOptimizationConditions() const;
   void OutputMOPlotConditions() const;
   void OutputHolePlotConditions() const;
   void OutputParticlePlotConditions() const;
   void OutputInputTerms(std::vector<std::string> inputTerms) const;
   bool IsCommentOut(std::string str) const;
   std::vector<std::string> GetInputTerms() const;
   int ParseMolecularGeometry(Molecule* molecule, std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseTheory(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsSCF(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsPrincipalAxes(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsTranslation(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsRotation(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsMOPlot(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsHolePlot(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsParticlePlot(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsCIS(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsMC(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsMD(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsRPMD(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsOptimization(std::vector<std::string>* inputTerms, int parseIndex) const;
   int ParseConditionsMemory(std::vector<std::string>* inputTerms, int parseIndex) const;
};

}
#endif





