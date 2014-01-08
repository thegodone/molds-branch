//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   // 
// Copyright (C) 2012-2014 Katsuhiko Nishimra                             // 
// Copyright (C) 2012-2014 Michihiro Okuyama                              //
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
#include<fstream>
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
#include"Utilities.h"
#include"EularAngle.h"
#include"Parameters.h"
#include"RealSphericalHarmonicsIndex.h"
#include"atoms/Atom.h"
#include"factories/AtomFactory.h"
#include"Molecule.h"
#include"InputParser.h"
using namespace std;
using namespace MolDS_base_atoms;
using namespace MolDS_base_factories;
namespace MolDS_base{

InputParser* InputParser::inputParser = NULL;

InputParser::InputParser(){
   this->SetMessages();
}

InputParser::~InputParser(){
}

InputParser* InputParser::GetInstance(){
   if(inputParser == NULL){
      inputParser = new InputParser();
   }
   return inputParser;
}

void InputParser::DeleteInstance(){
   if(inputParser != NULL){
      delete inputParser; 
   }
   inputParser = NULL;
}

void InputParser::SetMessages(){
   this->errorMessageInputFileEmpty
      = "Error in base::InputParser::GetInputTerms: Input file is empty.\n"; 
   this->errorMessageNotFoundInputFile
      = "Error in base::InputParser::StoreInputTermsFromFile: Not found.\n"; 
   this->errorMessageNonValidTheoriesEpc
      = "Error in base::InputParser::ValidateEpcConditions: Theory you set is not supported for EPC. The supported theories are MNDO-sefies(MNDO, AM1, AM1D, PM3, PM3D, PDDG/PM3) only.\n";
   this->errorMessageNonValidTheoriesMD
      = "Error in base::InputParser::ValidateMdConditions: Theory you set is not supported for MD.\n";
   this->errorMessageNonValidExcitedStatesMD
      = "Error in base::InputParser::ValidateMdConditions: Excited state on which MD runs or CIS condition are wrong.\n";
   this->errorMessageNonValidExcitedStatesMC
      = "Error in base::InputParser::ValidateMcConditions: Excited state on which MC runs or CIS condition are wrong.\n";
   this->errorMessageNonValidTheoriesRPMD
      = "Error in base::InputParser::ValidateRpmdConditions: heory you set is not supported for RMPD.\n";
   this->errorMessageNonValidExcitedStatesRPMD
      = "Error in base::InputParser::ValidateRpmdConditions: Excited state on which RPMD runs or CIS condition are wrong.\n";
   this->errorMessageNonValidTheoriesNASCO
      = "Error in base::InputParser::ValidateNascoConditions: Theory you set is not supported for NASCO.\n";
   this->errorMessageNonValidNumberExcitedStatesNASCO
      = "Error in base::InputParser::ValidateNascoConditions: The Number of electronic states of NASCO should be not over the number of CIS excited states plus 1.\n";
   this->errorMessageNonValidInitialElectronicStateNASCO
      = "Error in base::InputParser::ValidateNascoConditions: The initial electronic states for NASCO should be set to one of the electronic eigenstates used in NASCO.\n";
   this->errorMessageNonValidTheoriesOptimization
      = "Error in base::InputParser::ValidateOptimizationConditions: heory you set is not supported for optimization.\n";
   this->errorMessageNonValidExcitedStatesOptimization
      = "Error in base::InputParser::ValidateOptimizationConditions: Excited state on which optimization is carried out or CIS condition are wrong.\n";
   this->errorMessageNonValidElectronicStateFrequencies
      = "Error in base::InputParser::ValidateFrequenciesConditions: Excited states are not supported for the frequencies (normal modes) analysis.\n";
   this->errorMessageNonValidTheoryFrequencies
      = "Error in base::InputParser::ValidateFrequenciesConditions: CNDO2, INDO, and ZINDO/S are supported for the frequencies (normal modes) analysis.\n";
   this->errorMessageElecState = "Electronic eigenstate: ";
   this->errorMessageInputFile = "Inputfile: "; 
   this->errorMessageTheory = "Theory: ";
   this->errorMessageNumberExcitedStateCIS = "Number of CIS excited states: ";
   this->errorMessageNumberElectronicStatesNASCO = "Number of electronic states for NASCO: ";
   this->errorMessageInitialElectronicStateNASCO = "Initial electronic state for NASCO: ";
   this->messageStartParseInput = "**********  START: Parse input  **********\n";
   this->messageDoneParseInput =  "**********  DONE: Parse input  ***********\n\n\n";
   this->messageTotalNumberAOs = "\tTotal number of valence AOs: ";
   this->messageTotalNumberAtoms = "\tTotal number of atoms: ";
   this->messageTotalNumberValenceElectrons = "\tTotal number of valence electrons: ";
   this->messageInputTerms = "Input terms:\n";

   // SCF
   this->messageScfConditions       = "\tSCF conditions:\n";
   this->messageScfMaxIterations    = "\t\tMax iterations: ";
   this->messageScfRmsDensity       = "\t\tRMS density: ";
   this->messageScfDampingThresh    = "\t\tDamping threshold: ";
   this->messageScfDampingWeight    = "\t\tDamping weight: ";
   this->messageScfDiisNumErrorVect = "\t\tDIIS number of error vectors: ";
   this->messageScfDiisStartError   = "\t\tDIIS starting error: ";
   this->messageScfDiisEndError     = "\t\tDIIS ending error: ";
   this->messageScfVdW              = "\t\tvan der Waals (vdW) correction: ";
   this->messageScfVdWScalingFactor = "\t\tvdW corr. scaling factor (s6): ";
   this->messageScfVdWDampingFactor = "\t\tvdW corr. damping factor (d): ";

   // CIS
   this->messageCisConditions                 = "\tCIS conditions:\n";
   this->messageCisNumberActiveOcc            = "\t\tNumber of active Occ.: ";
   this->messageCisNumberActiveVir            = "\t\tNumber of active Vir.: ";
   this->messageCisNumberExcitedStates        = "\t\tNumber of excited states: ";
   this->messageCisDavidson                   = "\t\tCIS-Davidson: ";
   this->messageCisNormTolerance              = "\t\tNorm tolerance for the residual of the Davidson: ";
   this->messageCisMaxIterations              = "\t\tMax iterations for the Davidson: ";
   this->messageCisMaxDimensions              = "\t\tMax dimensions for the Davidson: ";
   this->messageCisExcitonEnergies            = "\t\tExciton energies: ";
   this->messageCisAllTransitionDipoleMoments = "\t\tAll transition dipole moments: ";
   this->messageCisNumPrintCoefficients       = "\t\tNumber of printed coefficients of CIS-eigenvector: ";
   this->messageCisMulliken                   = "\t\tMulliken population of excited states: ";

   // memory
   this->messageMemoryConditions = "\tMemory conditions:\n";
   this->messageMemoryLimitHeap  = "\t\tHeap limit: ";

   // MD
   this->messageMdConditions = "\tMD conditions:\n";
   this->messageMdTotalSteps = "\t\tTotal steps: ";
   this->messageMdElecState  = "\t\tElectronic eigenstate: ";
   this->messageMdTimeWidth  = "\t\tTime width(dt): ";

   // MC
   this->messageMcConditions  = "\tMC conditions:\n";
   this->messageMcTotalSteps  = "\t\tTotal steps: ";
   this->messageMcElecState   = "\t\tElectronic eigenstate: ";
   this->messageMcStepWidth   = "\t\tStep width: ";
   this->messageMcTemperature = "\t\tTemperature: ";
   this->messageMcSeed        = "\t\tSeed: ";

   // RPMD
   this->messageRpmdConditions    = "\tRPMD conditions:\n";
   this->messageRpmdTotalSteps    = "\t\tTotal steps: ";
   this->messageRpmdElecState     = "\t\tElectronic eigenstate: ";
   this->messageRpmdNumElecStates = "\t\tNumber of the electronic eigenstates: ";
   this->messageRpmdTimeWidth     = "\t\tTime width: ";
   this->messageRpmdTemperature   = "\t\tTemperature: ";
   this->messageRpmdNumBeads      = "\t\tNumber of the beads in the Ring Polymer: ";
   this->messageRpmdSeed          = "\t\tSeed: ";

   // NASCO
   this->messageNascoConditions       = "\tNasco conditions:\n";
   this->messageNascoTotalSteps       = "\t\tTotal steps: ";
   this->messageNascoNumElecStates    = "\t\tNumber of the electronic eigenstates: ";
   this->messageNascoInitialElecState = "\t\tInitial electronic eigenstate: ";
   this->messageNascoTimeWidth        = "\t\tTime width: ";
   this->messageNascoSeed             = "\t\tSeed: ";

   // Optimization
   this->messageOptimizationConditions  = "\tOptimization conditions:\n";
   this->messageOptimizationMethod      = "\t\tMethod: ";
   this->messageOptimizationTotalSteps  = "\t\tTotal steps: ";
   this->messageOptimizationElecState   = "\t\tElectronic eigenstate: ";
   this->messageOptimizationMaxGradient = "\t\tMax gradient: ";
   this->messageOptimizationRmsGradient = "\t\tRms gradient: ";
   this->messageOptimizationTimeWidth   = "\t\tFictious time width: ";
   this->messageOptimizationInitialTrustRadius = "\t\tInitial trust radius: ";
   this->messageOptimizationMaxNormStep        = "\t\tMax size of the optimization step: ";

   // Frequencies (Normal modes)
   this->messageFrequenciesConditions    = "\tFrequencies (Normal modes) analysis conditions:\n";
   this->messageFrequenciesElecState     = "\t\tElectronic eigenstate: ";

   // MOPlot
   this->messageMOPlotConditions  = "\tMO plot conditions:\n";
   this->messageMOPlotIndex       = "\t\tMO index: ";
   this->messageMOPlotGridNumber  = "\t\tNumber of grid(x, y, z): ";
   this->messageMOPlotFrameLength = "\t\tFrame length[angst.](x, y, z): ";
   this->messageMOPlotFilePrefix  = "\t\tFile name prefix: ";

   // HolePlot
   this->messageHolePlotConditions  = "\tHole plot conditions:\n";
   this->messageHolePlotElecIndex   = "\t\tElectronic index: ";
   this->messageHolePlotGridNumber  = "\t\tNumber of grid(x, y, z): ";
   this->messageHolePlotFrameLength = "\t\tFrame length[angst.](x, y, z): ";
   this->messageHolePlotFilePrefix  = "\t\tFile name prefix: ";

   // ParticlePlot
   this->messageParticlePlotConditions  = "\tParticle plot conditions:\n";
   this->messageParticlePlotElecIndex   = "\t\tElectronic state: ";
   this->messageParticlePlotGridNumber  = "\t\tNumber of grid(x, y, z): ";
   this->messageParticlePlotFrameLength = "\t\tFrame length[angst.](x, y, z): ";
   this->messageParticlePlotFilePrefix  = "\t\tFile name prefix: ";

   // unit
   this->messageFs     = "[fs]";
   this->messageK      = "[K]";
   this->messageAngst  = "[Angst.]";
   this->messageMB     = "[MB]";

   // others
   this->stringYES   = "yes";
   this->stringNO    = "no";
   this->stringSpace = " ";
   this->stringTab   = "\t";

   // theory
   this->stringCommentOut    = "//";
   this->stringTheoryCNDO2   = "cndo/2";
   this->stringTheoryINDO    = "indo";
   this->stringTheoryZINDOS  = "zindo/s";
   this->stringTheoryMNDO    = "mndo";
   this->stringTheoryAM1     = "am1";
   this->stringTheoryAM1D    = "am1-d";
   this->stringTheoryPM3     = "pm3";
   this->stringTheoryPM3D    = "pm3-d";
   this->stringTheoryPM3PDDG = "pm3/pddg";
   this->stringTheory        = "theory";
   this->stringTheoryEnd     = "theory_end";

   // molecular configuratio
   this->stringGeometry    = "geometry";
   this->stringGeometryEnd = "geometry_end";

   // Ghost (ghost) atoms
   this->stringGhost       = "ghost";
   this->stringGhostEnd    = "ghost_end";

   // Environmental Point Charge
   this->stringEpc       = "epc";
   this->stringEpcEnd    = "epc_end";
   this->stringEpcCharge = "charge";

   // SCF
   this->stringScf                 = "scf";
   this->stringScfEnd              = "scf_end";
   this->stringScfMaxIter          = "max_iter";
   this->stringScfRmsDensity       = "rms_density";
   this->stringScfDampingThresh    = "damping_thresh";
   this->stringScfDampingWeight    = "damping_weight";
   this->stringScfDiisNumErrorVect = "diis_num_error_vect";
   this->stringScfDiisStartError   = "diis_start_error";
   this->stringScfDiisEndError     = "diis_end_error";
   this->stringScfVdW              = "vdw";
   this->stringScfVdWScalingFactor = "vdw_s6";
   this->stringScfVdWDampingFactor = "vdw_d";

   // MO plot
   this->stringMO                = "mo";
   this->stringMOPlot            = "moplot";
   this->stringMOPlotEnd         = "moplot_end";
   this->stringMOPlotGridNumber  = "grid_number";
   this->stringMOPlotFrameLength = "frame_length";
   this->stringMOPlotFilePrefix  = "file_prefix";

   // Hole plot
   this->stringHolePlot            = "holeplot";
   this->stringHolePlotEnd         = "holeplot_end";
   this->stringHolePlotElecIndex   = "electronic_state";
   this->stringHolePlotGridNumber  = "grid_number";
   this->stringHolePlotFrameLength = "frame_length";
   this->stringHolePlotFilePrefix  = "file_prefix";

   // MO plot
   this->stringParticlePlot            = "particleplot";
   this->stringParticlePlotEnd         = "particleplot_end";
   this->stringParticlePlotElecIndex   = "electronic_state";
   this->stringParticlePlotGridNumber  = "grid_number";
   this->stringParticlePlotFrameLength = "frame_length";
   this->stringParticlePlotFilePrefix  = "file_prefix";

   // Principal axes
   this->stringInertiaTensor       = "inertia";
   this->stringInertiaTensorEnd    = "inertia_end";
   this->stringInertiaTensorOrigin = "origin";

   // Rotate
   this->stringRotate                 = "rotate";
   this->stringRotateEnd              = "rotate_end";
   this->stringRotatingOrigin         = "origin";
   this->stringRotatingAxis           = "axis";
   this->stringRotatingAngle          = "angle";
   this->stringRotatingAngles         = "angles";
   this->stringRotatingType           = "type";
   this->stringRotatingTypeAxis       = "axis";
   this->stringRotatingTypeEularAngle = "eular_angle";

   // Translate
   this->stringTranslate             = "translate";
   this->stringTranslateEnd          = "translate_end";
   this->stringTranslatingDifference = "difference";

   // CIS
   this->stringCIS                           = "cis";
   this->stringCISEnd                        = "cis_end";
   this->stringCISActiveOcc                  = "active_occ";
   this->stringCISActiveVir                  = "active_vir";
   this->stringCISNStates                    = "nstates";
   this->stringCISDavidson                   = "davidson";
   this->stringCISMaxIter                    = "max_iter";
   this->stringCISMaxDimensions              = "max_dim";
   this->stringCISNormTolerance              = "norm_tol";
   this->stringCISExcitonEnergies            = "exciton_energies";
   this->stringCISAllTransitionDipoleMoments = "all_transition_dipole_moments";
   this->stringCISNumPrintCoefficients       = "num_print_coefficients";
   this->stringCISMulliken                   = "mulliken";
   this->stringCISUnpairedPop                = "unpaired_electron_population";

   // Memory
   this->stringMemory          = "memory";
   this->stringMemoryEnd       = "memory_end";
   this->stringMemoryLimitHeap = "limit_heap";

   // MD
   this->stringMD           = "md";
   this->stringMDEnd        = "md_end";
   this->stringMDTotalSteps = "total_steps";
   this->stringMDElecState  = "electronic_state";
   this->stringMDTimeWidth  = "dt";

   // MC
   this->stringMC            = "mc";
   this->stringMCEnd         = "mc_end";
   this->stringMCTotalSteps  = "total_steps";
   this->stringMCElecState   = "electronic_state";
   this->stringMCStepWidth   = "step_width";
   this->stringMCTemperature = "temperature";
   this->stringMCSeed        = "seed";

   // RPMD
   this->stringRPMD              = "rpmd";
   this->stringRPMDEnd           = "rpmd_end";
   this->stringRPMDTotalSteps    = "total_steps";
   this->stringRPMDElecState     = "electronic_state";
   this->stringRPMDNumElecStates = "num_electronic_states";
   this->stringRPMDTimeWidth     = "dt";
   this->stringRPMDTemperature   = "temperature";
   this->stringRPMDNumBeads      = "num_beads";
   this->stringRPMDSeed          = "seed";

   // NASCO
   this->stringNASCO                 = "nasco";
   this->stringNASCOEnd              = "nasco_end";
   this->stringNASCOTotalSteps       = "total_steps";
   this->stringNASCONumElecStates    = "num_electronic_states";
   this->stringNASCOInitialElecState = "initial_electronic_state";
   this->stringNASCOTimeWidth        = "dt";
   this->stringNASCOSeed             = "seed";

   // Opt
   this->stringOptimization                  = "optimization";
   this->stringOptimizationEnd               = "optimization_end";
   this->stringOptimizationMethod            = "method";
   this->stringOptimizationBFGS              = "bfgs";
   this->stringOptimizationGEDIIS            = "gediis";
   this->stringOptimizationConjugateGradient = "conjugate_gradient";
   this->stringOptimizationSteepestDescent   = "steepest_descent";
   this->stringOptimizationTotalSteps        = "total_steps";
   this->stringOptimizationElecState         = "electronic_state";
   this->stringOptimizationMaxGradient       = "max_gradient";
   this->stringOptimizationRmsGradient       = "rms_gradient";
   this->stringOptimizationTimeWidth         = "dt";
   this->stringOptimizationInitialTrustRadius = "initial_trust_radius";
   this->stringOptimizationMaxNormStep        = "max_norm_step";

   // Frequencies (Normal modes)
   this->stringFrequencies          = "frequencies";
   this->stringFrequenciesEnd       = "frequencies_end";
   this->stringFrequenciesElecState = "electronic_state";
}

vector<string> InputParser::GetInputTerms(int argc, char *argv[]) const{
   vector<string> inputTerms;
   if(argc==1){
      this->StoreInputTermsFromRedirect(inputTerms);
   }
   else{
      char* fileName = argv[1];
      this->StoreInputTermsFromFile(inputTerms,fileName);
   }
   if (inputTerms.size() == 0) {
      char* fileName = argv[1];
      stringstream ss; 
      ss << this->errorMessageInputFileEmpty;  
      if (argc > 1) {
         ss << this->errorMessageInputFile << fileName << endl; 
      }
      throw MolDSException(ss.str()); 
   }
   return inputTerms;
}

void InputParser::StoreInputTermsFromRedirect(vector<string>& inputTerms) const{
   string str;
   while(getline(cin, str)){
      this->AddInputTermsFromString(inputTerms, str);
   }
}

void InputParser::StoreInputTermsFromFile(vector<string>& inputTerms, char* fileName) const{
   string str;
   fstream ifs(fileName);
   if (ifs.fail()) { 
      stringstream ss;
      ss << this->errorMessageNotFoundInputFile; 
      ss << this->errorMessageInputFile << fileName << endl; 
      throw MolDSException(ss.str()); 
   }

   while(getline(ifs, str)){
      this->AddInputTermsFromString(inputTerms, str);
   }
}

void InputParser::AddInputTermsFromString(vector<string>& inputTerms, string str) const{
   //skip comment
   if(this->IsCommentOut(str)){return;}

   string inputTerm;
   bool isPreCharSpace = true;
   // get input terms
   for(int i=0; i<str.length(); i++){
      if(str.data()[i] != stringSpace.data()[0] && str.data()[i] != stringTab.data()[0]){
         // change to lower case.
         inputTerm += tolower(str.data()[i]);
         isPreCharSpace = false;
      }
      else{
         if(!isPreCharSpace){
            inputTerms.push_back(inputTerm);
            inputTerm = "";
            isPreCharSpace = true;
         }
      }
   }
   if(inputTerm.length()>0){
      inputTerms.push_back(inputTerm);
   }
}

int InputParser::ParseMolecularConfiguration(Molecule* molecule, vector<string>* inputTerms, int parseIndex) const{
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringGeometryEnd) != 0){
      double x = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
      double y = atof((*inputTerms)[parseIndex+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
      double z = atof((*inputTerms)[parseIndex+3].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
      AtomType atomType = H;
      if((*inputTerms)[parseIndex] == "h"){
        atomType = H;
      }
      else if((*inputTerms)[parseIndex] == "li"){
         atomType = Li;
      }
      else if((*inputTerms)[parseIndex] == "c"){
         atomType = C;
      }
      else if((*inputTerms)[parseIndex] == "n"){
         atomType = N;
      }
      else if((*inputTerms)[parseIndex] == "o"){
         atomType = O;
      }
      else if((*inputTerms)[parseIndex] == "s"){
         atomType = S;
      }
      else if((*inputTerms)[parseIndex] == "zn"){
         atomType = Zn;
      }
      int index = molecule->GetRealAtomVect().size() + molecule->GetGhostAtomVect().size();
      Atom* atom = AtomFactory::Create(atomType, index, x, y, z);
      molecule->AddRealAtom(atom);
      parseIndex += 4;
   }
   return parseIndex;
}

int InputParser::ParseGhostsConfiguration(Molecule* molecule, vector<string>* inputTerms, int parseIndex) const{
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringGhostEnd) != 0){
      double x = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
      double y = atof((*inputTerms)[parseIndex+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
      double z = atof((*inputTerms)[parseIndex+3].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
      AtomType atomType = H;
      if((*inputTerms)[parseIndex] == "h"){
        atomType = ghostH;
      }
      else if((*inputTerms)[parseIndex] == "li"){
         atomType = ghostLi;
      }
      else if((*inputTerms)[parseIndex] == "c"){
         atomType = ghostC;
      }
      else if((*inputTerms)[parseIndex] == "n"){
         atomType = ghostN;
      }
      else if((*inputTerms)[parseIndex] == "o"){
         atomType = ghostO;
      }
      else if((*inputTerms)[parseIndex] == "s"){
         atomType = ghostS;
      }
      else if((*inputTerms)[parseIndex] == "zn"){
         atomType = ghostZn;
      }
      int index = molecule->GetRealAtomVect().size() + molecule->GetGhostAtomVect().size();
      Atom* atom = AtomFactory::Create(atomType, index, x, y, z);
      molecule->AddGhostAtom(atom);
      parseIndex += 4;
   }
   return parseIndex;
}

int InputParser::ParseEpcsConfiguration(Molecule* molecule, vector<string>* inputTerms, int parseIndex) const{
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringEpcEnd) != 0){
      double x      = atof((*inputTerms)[parseIndex+0].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
      double y      = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
      double z      = atof((*inputTerms)[parseIndex+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
      parseIndex += 3;
      double charge = 0.0;
      if((*inputTerms)[parseIndex].compare(this->stringEpcCharge) == 0){
         charge = atof((*inputTerms)[parseIndex+1].c_str());
         parseIndex += 2;
      }
      AtomType atomType = EPC;
      int index = molecule->GetEpcVect().size();
      Atom* atom = AtomFactory::Create(atomType, index, x, y, z, charge);
      molecule->AddEpc(atom);
   }
   return parseIndex;
}

int InputParser::ParseConditionsSCF(vector<string>* inputTerms, int parseIndex) const{
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringScfEnd) != 0){
      // max iterations
      if((*inputTerms)[parseIndex].compare(this->stringScfMaxIter) == 0){
         Parameters::GetInstance()->SetMaxIterationsSCF(atoi((*inputTerms)[parseIndex+1].c_str()));
         parseIndex++;
      }
      // RMS density 
      if((*inputTerms)[parseIndex].compare(this->stringScfRmsDensity) == 0){
         Parameters::GetInstance()->SetThresholdSCF(atof((*inputTerms)[parseIndex+1].c_str()));
         parseIndex++;
      }
      // Damping Threshold 
      if((*inputTerms)[parseIndex].compare(this->stringScfDampingThresh) == 0){
         Parameters::GetInstance()->SetDampingThreshSCF(atof((*inputTerms)[parseIndex+1].c_str()));
         parseIndex++;
      }
      // Damping Weight
      if((*inputTerms)[parseIndex].compare(this->stringScfDampingWeight) == 0){
         Parameters::GetInstance()->SetDampingWeightSCF(atof((*inputTerms)[parseIndex+1].c_str()));
         parseIndex++;
      }
      // DIIS number of stored error vectors
      if((*inputTerms)[parseIndex].compare(this->stringScfDiisNumErrorVect) == 0){
         Parameters::GetInstance()->SetDiisNumErrorVectSCF(atoi((*inputTerms)[parseIndex+1].c_str()));
         parseIndex++;
      }
      // DIIS starting error
      if((*inputTerms)[parseIndex].compare(this->stringScfDiisStartError) == 0){
         Parameters::GetInstance()->SetDiisStartErrorSCF(atof((*inputTerms)[parseIndex+1].c_str()));
         parseIndex++;
      }
      // DIIS ending error
      if((*inputTerms)[parseIndex].compare(this->stringScfDiisEndError) == 0){
         Parameters::GetInstance()->SetDiisEndErrorSCF(atof((*inputTerms)[parseIndex+1].c_str()));
         parseIndex++;
      }
      // van der Waals correction 
      if((*inputTerms)[parseIndex].compare(this->stringScfVdW) == 0){
         if((*inputTerms)[parseIndex+1].compare(this->stringYES) == 0){
            Parameters::GetInstance()->SetRequiresVdWSCF(true);
         }
         else{
            Parameters::GetInstance()->SetRequiresVdWSCF(false);
         }
         parseIndex++;
      }
      // van der Waals (scaling factor) 
      if((*inputTerms)[parseIndex].compare(this->stringScfVdWScalingFactor) == 0){
         Parameters::GetInstance()->SetVdWScalingFactorSCF(atof((*inputTerms)[parseIndex+1].c_str()));
         parseIndex++;
      }
      // van der Waals (damping factor) 
      if((*inputTerms)[parseIndex].compare(this->stringScfVdWDampingFactor) == 0){
         Parameters::GetInstance()->SetVdWDampingFactorSCF(atof((*inputTerms)[parseIndex+1].c_str()));
         parseIndex++;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsPrincipalAxes(vector<string>* inputTerms, int parseIndex) const{
   Parameters::GetInstance()->SetCurrentSimulation(PrincipalAxes);
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringInertiaTensorEnd) != 0){
      // origin
      if((*inputTerms)[parseIndex].compare(this->stringInertiaTensorOrigin) == 0){
         double x = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         double y = atof((*inputTerms)[parseIndex+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         double z = atof((*inputTerms)[parseIndex+3].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         Parameters::GetInstance()->SetInertiaTensorOrigin(x, y, z);
         parseIndex+=3;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsTranslation(vector<string>* inputTerms, int parseIndex) const{
   Parameters::GetInstance()->SetCurrentSimulation(Translate);
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringTranslateEnd) != 0){
      // origin
      if((*inputTerms)[parseIndex].compare(this->stringTranslatingDifference) == 0){
         double x = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         double y = atof((*inputTerms)[parseIndex+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         double z = atof((*inputTerms)[parseIndex+3].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         Parameters::GetInstance()->SetTranslatingDifference(x, y, z);
         parseIndex+=3;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsRotation(vector<string>* inputTerms, int parseIndex) const{
   Parameters::GetInstance()->SetCurrentSimulation(Rotate);
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringRotateEnd) != 0){
      // origin
      if((*inputTerms)[parseIndex].compare(this->stringRotatingOrigin) == 0){
         double x = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         double y = atof((*inputTerms)[parseIndex+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         double z = atof((*inputTerms)[parseIndex+3].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         Parameters::GetInstance()->SetRotatingOrigin(x, y, z);
         parseIndex+=3;
      }
      // axis
      else if((*inputTerms)[parseIndex].compare(this->stringRotatingAxis) == 0){
         double x = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         double y = atof((*inputTerms)[parseIndex+2].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         double z = atof((*inputTerms)[parseIndex+3].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         Parameters::GetInstance()->SetRotatingAxis(x, y, z);
         parseIndex+=3;
      }
      // angle 
      else if((*inputTerms)[parseIndex].compare(this->stringRotatingAngle) == 0){
         double angle = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetDegree2Radian();
         Parameters::GetInstance()->SetRotatingAngle(angle);
         parseIndex++;
      }
      // angles (EularAngle)
      else if((*inputTerms)[parseIndex].compare(this->stringRotatingAngles) == 0){
         double alpha = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetDegree2Radian();
         double beta  = atof((*inputTerms)[parseIndex+2].c_str()) * Parameters::GetInstance()->GetDegree2Radian();
         double gamma = atof((*inputTerms)[parseIndex+3].c_str()) * Parameters::GetInstance()->GetDegree2Radian();
         Parameters::GetInstance()->SetRotatingEularAngles(alpha, beta, gamma);
         parseIndex += 3;
      }
      // type
      else if((*inputTerms)[parseIndex].compare(this->stringRotatingType) == 0){
         if((*inputTerms)[parseIndex+1].compare(this->stringRotatingTypeAxis) == 0){
            Parameters::GetInstance()->SetRotatingType(Axis);
         }
         else if((*inputTerms)[parseIndex+1].compare(this->stringRotatingTypeEularAngle) == 0){
            Parameters::GetInstance()->SetRotatingType(Eular);
         }
         parseIndex++;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsMOPlot(vector<string>* inputTerms, int parseIndex) const{
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringMOPlotEnd) != 0){
      // Frame length
      if((*inputTerms)[parseIndex].compare(this->stringMOPlotFrameLength) == 0){
         double lx = atof((*inputTerms)[parseIndex+1].c_str()) 
                    *Parameters::GetInstance()->GetAngstrom2AU();
         double ly = atof((*inputTerms)[parseIndex+2].c_str()) 
                    *Parameters::GetInstance()->GetAngstrom2AU();
         double lz = atof((*inputTerms)[parseIndex+3].c_str()) 
                    *Parameters::GetInstance()->GetAngstrom2AU();
         Parameters::GetInstance()->SetFrameLengthMOPlot(lx, ly, lz);
         parseIndex += 3;
      }
      // Grid number
      if((*inputTerms)[parseIndex].compare(this->stringMOPlotGridNumber) == 0){
         int nx = atof((*inputTerms)[parseIndex+1].c_str());
         int ny = atof((*inputTerms)[parseIndex+2].c_str());
         int nz = atof((*inputTerms)[parseIndex+3].c_str());
         Parameters::GetInstance()->SetGridNumberMOPlot(nx, ny, nz);
         parseIndex += 3;
      }
      // mo index
      if((*inputTerms)[parseIndex].compare(this->stringMO) == 0){
         int moIndex = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->AddIndexMOPlot(moIndex);
         parseIndex++;
      }
      // file prefix
      if((*inputTerms)[parseIndex].compare(this->stringMOPlotFilePrefix) == 0){
         string filePrefix((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetFileNamePrefixMOPlot(filePrefix);
         parseIndex++;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsHolePlot(vector<string>* inputTerms, int parseIndex) const{
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringHolePlotEnd) != 0){
      // Frame length
      if((*inputTerms)[parseIndex].compare(this->stringHolePlotFrameLength) == 0){
         double lx = atof((*inputTerms)[parseIndex+1].c_str()) 
                    *Parameters::GetInstance()->GetAngstrom2AU();
         double ly = atof((*inputTerms)[parseIndex+2].c_str()) 
                    *Parameters::GetInstance()->GetAngstrom2AU();
         double lz = atof((*inputTerms)[parseIndex+3].c_str()) 
                    *Parameters::GetInstance()->GetAngstrom2AU();
         Parameters::GetInstance()->SetFrameLengthHolePlot(lx, ly, lz);
         parseIndex += 3;
      }
      // Grid number
      if((*inputTerms)[parseIndex].compare(this->stringHolePlotGridNumber) == 0){
         int nx = atof((*inputTerms)[parseIndex+1].c_str());
         int ny = atof((*inputTerms)[parseIndex+2].c_str());
         int nz = atof((*inputTerms)[parseIndex+3].c_str());
         Parameters::GetInstance()->SetGridNumberHolePlot(nx, ny, nz);
         parseIndex += 3;
      }
      // hole index
      if((*inputTerms)[parseIndex].compare(this->stringHolePlotElecIndex) == 0){
         int elecIndex = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->AddElecIndexHolePlot(elecIndex);
         parseIndex++;
      }
      // file prefix
      if((*inputTerms)[parseIndex].compare(this->stringHolePlotFilePrefix) == 0){
         string filePrefix((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetFileNamePrefixHolePlot(filePrefix);
         parseIndex++;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsParticlePlot(vector<string>* inputTerms, int parseIndex) const{
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringParticlePlotEnd) != 0){
      // Frame length
      if((*inputTerms)[parseIndex].compare(this->stringParticlePlotFrameLength) == 0){
         double lx = atof((*inputTerms)[parseIndex+1].c_str()) 
                    *Parameters::GetInstance()->GetAngstrom2AU();
         double ly = atof((*inputTerms)[parseIndex+2].c_str()) 
                    *Parameters::GetInstance()->GetAngstrom2AU();
         double lz = atof((*inputTerms)[parseIndex+3].c_str()) 
                    *Parameters::GetInstance()->GetAngstrom2AU();
         Parameters::GetInstance()->SetFrameLengthParticlePlot(lx, ly, lz);
         parseIndex += 3;
      }
      // Grid number
      if((*inputTerms)[parseIndex].compare(this->stringParticlePlotGridNumber) == 0){
         int nx = atof((*inputTerms)[parseIndex+1].c_str());
         int ny = atof((*inputTerms)[parseIndex+2].c_str());
         int nz = atof((*inputTerms)[parseIndex+3].c_str());
         Parameters::GetInstance()->SetGridNumberParticlePlot(nx, ny, nz);
         parseIndex += 3;
      }
      // particle index
      if((*inputTerms)[parseIndex].compare(this->stringParticlePlotElecIndex) == 0){
         int elecIndex = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->AddElecIndexParticlePlot(elecIndex);
         parseIndex++;
      }
      // file prefix
      if((*inputTerms)[parseIndex].compare(this->stringParticlePlotFilePrefix) == 0){
         string filePrefix((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetFileNamePrefixParticlePlot(filePrefix);
         parseIndex++;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsCIS(vector<string>* inputTerms, int parseIndex) const{
   Parameters::GetInstance()->SetRequiresCIS(true);
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringCISEnd) != 0){
      // number of active occupied orbitals
      if((*inputTerms)[parseIndex].compare(this->stringCISActiveOcc) == 0){
         int activeOccCIS = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetActiveOccCIS(activeOccCIS);
         parseIndex++;
      }
      // number of active virtual orbitals
      if((*inputTerms)[parseIndex].compare(this->stringCISActiveVir) == 0){
         int activeVirCIS = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetActiveVirCIS(activeVirCIS);
         parseIndex++;
      }
      // number of excited states
      if((*inputTerms)[parseIndex].compare(this->stringCISNStates) == 0){
         int nStates = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetNumberExcitedStatesCIS(nStates);
         parseIndex++;
      }
      // Davidson is used or not
      if((*inputTerms)[parseIndex].compare(this->stringCISDavidson) == 0){
         if((*inputTerms)[parseIndex+1].compare(this->stringYES) == 0){
            Parameters::GetInstance()->SetIsDavidsonCIS(true);
         }
         else{
            Parameters::GetInstance()->SetIsDavidsonCIS(false);
         }
         parseIndex++;
      }
      // max iterations for the Davidson loop
      if((*inputTerms)[parseIndex].compare(this->stringCISMaxIter) == 0){
         int maxIter = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetMaxIterationsCIS(maxIter);
         parseIndex++;
      }
      // max dimensions for the Davidson expansion
      if((*inputTerms)[parseIndex].compare(this->stringCISMaxDimensions) == 0){
         int maxDim = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetMaxDimensionsCIS(maxDim);
         parseIndex++;
      }
      // nolm tolerance for the norm of the resiudal vectors of the Davidson.
      if((*inputTerms)[parseIndex].compare(this->stringCISNormTolerance) == 0){
         double normTol = atof((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetNormToleranceCIS(normTol);
         parseIndex++;
      }
      // max dimensions for the Davidson expansion
      if((*inputTerms)[parseIndex].compare(this->stringCISNumPrintCoefficients) == 0){
         int numPrintCoeff = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetNumberPrintCoefficientsCIS(numPrintCoeff);
         parseIndex++;
      }
      // exciton energies are calculated or not
      if((*inputTerms)[parseIndex].compare(this->stringCISExcitonEnergies) == 0){
         if((*inputTerms)[parseIndex+1].compare(this->stringYES) == 0){
            Parameters::GetInstance()->SetRequiresExcitonEnergiesCIS(true);
         }
         else{
            Parameters::GetInstance()->SetRequiresExcitonEnergiesCIS(false);
         }
         parseIndex++;
      }
      // all transition dipole moments are calculated or not
      if((*inputTerms)[parseIndex].compare(this->stringCISAllTransitionDipoleMoments) == 0){
         if((*inputTerms)[parseIndex+1].compare(this->stringYES) == 0){
            Parameters::GetInstance()->SetRequiresAllTransitionDipoleMomentsCIS(true);
         }
         else{
            Parameters::GetInstance()->SetRequiresAllTransitionDipoleMomentsCIS(false);
         }
         parseIndex++;
      }
      // mulliken
      if((*inputTerms)[parseIndex].compare(this->stringCISMulliken) == 0){
         int groundStateIndex = 0;
         int elecStateIndex = atoi((*inputTerms)[parseIndex+1].c_str());
         if(groundStateIndex<elecStateIndex){
            Parameters::GetInstance()->AddElectronicStateIndexMullikenCIS(elecStateIndex);
         }
         parseIndex++;
      }
      // unpaired electron population
      if((*inputTerms)[parseIndex].compare(this->stringCISUnpairedPop) == 0){
          if((*inputTerms)[parseIndex+1].compare(this->stringYES) == 0){
             Parameters::GetInstance()->SetRequiresUnpairedPopCIS(true);
          }
          else {
             Parameters::GetInstance()->SetRequiresUnpairedPopCIS(false);
          }
          parseIndex++;   
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsMC(vector<string>* inputTerms, int parseIndex) const{
   Parameters::GetInstance()->SetCurrentSimulation(MC);
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringMCEnd) != 0){
      // number of total steps 
      if((*inputTerms)[parseIndex].compare(this->stringMCTotalSteps) == 0){
         int totalSteps = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetTotalStepsMC(totalSteps);
         parseIndex++;
      }
      // index of electronic eigen state on whichi MC runs. 
      if((*inputTerms)[parseIndex].compare(this->stringMCElecState) == 0){
         int elecStateIndex = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetElectronicStateIndexMC(elecStateIndex);
         parseIndex++;
      }
      // temperature for MC.
      if((*inputTerms)[parseIndex].compare(this->stringMCTemperature) == 0){
         double temperature = atof((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetTemperatureMC(temperature);
         parseIndex++;
      }
      // step width for MC.
      if((*inputTerms)[parseIndex].compare(this->stringMCStepWidth) == 0){
         double stepWidth = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetAngstrom2AU();
         Parameters::GetInstance()->SetStepWidthMC(stepWidth);
         parseIndex++;
      }
      // seed for MC.
      if((*inputTerms)[parseIndex].compare(this->stringMCSeed) == 0){
         unsigned long seed = atol((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetSeedMC(seed);
         parseIndex++;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsMD(vector<string>* inputTerms, int parseIndex) const{
   Parameters::GetInstance()->SetCurrentSimulation(MD);
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringMDEnd) != 0){
      // number of total steps 
      if((*inputTerms)[parseIndex].compare(this->stringMDTotalSteps) == 0){
         int totalSteps = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetTotalStepsMD(totalSteps);
         parseIndex++;
      }
      // index of electronic eigen state on whichi MD runs. 
      if((*inputTerms)[parseIndex].compare(this->stringMDElecState) == 0){
         int elecStateIndex = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetElectronicStateIndexMD(elecStateIndex);
         parseIndex++;
      }
      // time width for MD.
      if((*inputTerms)[parseIndex].compare(this->stringMDTimeWidth) == 0){
         double dt = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetFs2AU();
         Parameters::GetInstance()->SetTimeWidthMD(dt);
         parseIndex++;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsRPMD(vector<string>* inputTerms, int parseIndex) const{
   Parameters::GetInstance()->SetCurrentSimulation(RPMD);
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringRPMDEnd) != 0){
      // number of total steps 
      if((*inputTerms)[parseIndex].compare(this->stringRPMDTotalSteps) == 0){
         int totalSteps = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetTotalStepsRPMD(totalSteps);
         parseIndex++;
      }
      // index of electronic eigen state on which RPMD runs. 
      if((*inputTerms)[parseIndex].compare(this->stringRPMDElecState) == 0){
         int elecStateIndex = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetElectronicStateIndexRPMD(elecStateIndex);
         parseIndex++;
      }
      // number of the electronic eigenstates for nonadiabatic RPMD.
      if((*inputTerms)[parseIndex].compare(this->stringRPMDNumElecStates) == 0){
         int numElecStates = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetNumberElectronicStatesRPMD(numElecStates);
         parseIndex++;
      }
      // temperature for RPMD.
      if((*inputTerms)[parseIndex].compare(this->stringRPMDTemperature) == 0){
         double temperature = atof((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetTemperatureRPMD(temperature);
         parseIndex++;
      }
      // time width for RPMD.
      if((*inputTerms)[parseIndex].compare(this->stringRPMDTimeWidth) == 0){
         double timeWidth = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetFs2AU();
         Parameters::GetInstance()->SetTimeWidthRPMD(timeWidth);
         parseIndex++;
      }
      // number of the beads in Ring Polymer.
      if((*inputTerms)[parseIndex].compare(this->stringRPMDNumBeads) == 0){
         int numBeads = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetNumberBeadsRPMD(numBeads);
         parseIndex++;
      }
      // seed for RPMD.
      if((*inputTerms)[parseIndex].compare(this->stringRPMDSeed) == 0){
         unsigned long seed = atol((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetSeedRPMD(seed);
         parseIndex++;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsNASCO(vector<string>* inputTerms, int parseIndex) const{
   Parameters::GetInstance()->SetCurrentSimulation(NASCO);
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringNASCOEnd) != 0){
      // number of total steps 
      if((*inputTerms)[parseIndex].compare(this->stringNASCOTotalSteps) == 0){
         int totalSteps = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetTotalStepsNASCO(totalSteps);
         parseIndex++;
      }
      // number of the electronic eigenstates for NASCO.
      if((*inputTerms)[parseIndex].compare(this->stringNASCONumElecStates) == 0){
         int numElecStates = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetNumberElectronicStatesNASCO(numElecStates);
         parseIndex++;
      }
      // initial electronic eigenstates for NASCO.
      if((*inputTerms)[parseIndex].compare(this->stringNASCOInitialElecState) == 0){
         int initElecState = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetInitialElectronicStateNASCO(initElecState);
         parseIndex++;
      }
      // time width for NASCO.
      if((*inputTerms)[parseIndex].compare(this->stringNASCOTimeWidth) == 0){
         double timeWidth = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetFs2AU();
         Parameters::GetInstance()->SetTimeWidthNASCO(timeWidth);
         parseIndex++;
      }
      // seed for NASCO.
      if((*inputTerms)[parseIndex].compare(this->stringNASCOSeed) == 0){
         unsigned long seed = atol((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetSeedNASCO(seed);
         parseIndex++;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsOptimization(vector<string>* inputTerms, int parseIndex) const{
   Parameters::GetInstance()->SetCurrentSimulation(Optimization);
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringOptimizationEnd) != 0){
      // method 
      if((*inputTerms)[parseIndex].compare(this->stringOptimizationMethod) == 0){
         if((*inputTerms)[parseIndex+1].compare(this->stringOptimizationConjugateGradient) == 0){
            Parameters::GetInstance()->SetMethodOptimization(ConjugateGradientMethod);
         }
         else if((*inputTerms)[parseIndex+1].compare(this->stringOptimizationSteepestDescent) == 0){
            Parameters::GetInstance()->SetMethodOptimization(SteepestDescentMethod);
         }
         else if((*inputTerms)[parseIndex+1].compare(this->stringOptimizationBFGS) == 0){
            Parameters::GetInstance()->SetMethodOptimization(BFGSMethod);
         }
         else if((*inputTerms)[parseIndex+1].compare(this->stringOptimizationGEDIIS) == 0){
            Parameters::GetInstance()->SetMethodOptimization(GEDIISMethod);
         }
         else{
         }
         parseIndex++;
      }
      // number of steps of the optimization
      if((*inputTerms)[parseIndex].compare(this->stringOptimizationTotalSteps) == 0){
         int totalSteps = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetTotalStepsOptimization(totalSteps);
         parseIndex++;
      }
      // index of electronic eigen state on which optimization is carried out. 
      if((*inputTerms)[parseIndex].compare(this->stringOptimizationElecState) == 0){
         int elecStateIndex = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetElectronicStateIndexOptimization(elecStateIndex);
         parseIndex++;
      }
      // time width for the optimization.
      if((*inputTerms)[parseIndex].compare(this->stringOptimizationTimeWidth) == 0){
         double timeWidth = atof((*inputTerms)[parseIndex+1].c_str()) * Parameters::GetInstance()->GetFs2AU();
         Parameters::GetInstance()->SetTimeWidthOptimization(timeWidth);
         parseIndex++;
      }
      // max gradient for the optimization.
      if((*inputTerms)[parseIndex].compare(this->stringOptimizationMaxGradient) == 0){
         double maxGradient = atof((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetMaxGradientOptimization(maxGradient);
         parseIndex++;
      }
      // rms gradient for the optimization.
      if((*inputTerms)[parseIndex].compare(this->stringOptimizationRmsGradient) == 0){
         double rmsGradient = atof((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetRmsGradientOptimization(rmsGradient);
         parseIndex++;
      }
      // Initial trust radius.
      if((*inputTerms)[parseIndex].compare(this->stringOptimizationInitialTrustRadius) == 0){
         double initialTrustRadius = atof((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetInitialTrustRadiusOptimization(initialTrustRadius);
         parseIndex++;
      }
      // Max size of optimization step.
      if((*inputTerms)[parseIndex].compare(this->stringOptimizationMaxNormStep) == 0){
         double maxNormStep = atof((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetMaxNormStepOptimization(maxNormStep);
         parseIndex++;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseConditionsFrequencies(vector<string>* inputTerms, int parseIndex) const{
   Parameters::GetInstance()->SetRequiresFrequencies(true);
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringFrequenciesEnd) != 0){
      // electronic state on which the frequencies (normal modes) are calculated
      if((*inputTerms)[parseIndex].compare(this->stringFrequenciesElecState) == 0){
         int elecIndex = atoi((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetElectronicStateIndexFrequencies(elecIndex);
         parseIndex++;
      }
      parseIndex++;   
   }
   return parseIndex;
}

int InputParser::ParseTheory(vector<string>* inputTerms, int parseIndex) const{
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringTheoryEnd) != 0){
      // CNDO/2
      if((*inputTerms)[parseIndex].compare(this->stringTheoryCNDO2) == 0){
         Parameters::GetInstance()->SetCurrentTheory(CNDO2);
      }
      // INDO
      else if((*inputTerms)[parseIndex].compare(this->stringTheoryINDO) == 0){
         Parameters::GetInstance()->SetCurrentTheory(INDO);
      }
      // ZINDO/S
      else if((*inputTerms)[parseIndex].compare(this->stringTheoryZINDOS) == 0){
         Parameters::GetInstance()->SetCurrentTheory(ZINDOS);
      }
      // MNDO
      else if((*inputTerms)[parseIndex].compare(this->stringTheoryMNDO) == 0){
         Parameters::GetInstance()->SetCurrentTheory(MNDO);
      }
      // AM1
      else if((*inputTerms)[parseIndex].compare(this->stringTheoryAM1) == 0){
         Parameters::GetInstance()->SetCurrentTheory(AM1);
      }
      // AM1-D
      else if((*inputTerms)[parseIndex].compare(this->stringTheoryAM1D) == 0){
         Parameters::GetInstance()->SetCurrentTheory(AM1D);
      }
      // PM3
      else if((*inputTerms)[parseIndex].compare(this->stringTheoryPM3) == 0){
         Parameters::GetInstance()->SetCurrentTheory(PM3);
      }
      // PM3-D
      else if((*inputTerms)[parseIndex].compare(this->stringTheoryPM3D) == 0){
         Parameters::GetInstance()->SetCurrentTheory(PM3D);
      }
      // PM3/PDG
      else if((*inputTerms)[parseIndex].compare(this->stringTheoryPM3PDDG) == 0){
         Parameters::GetInstance()->SetCurrentTheory(PM3PDDG);
      }
      parseIndex++;
   }
   return parseIndex;
}

int InputParser::ParseConditionsMemory(vector<string>* inputTerms, int parseIndex) const{
   parseIndex++;
   while((*inputTerms)[parseIndex].compare(this->stringMemoryEnd) != 0){
      // max of heap
      if((*inputTerms)[parseIndex].compare(this->stringMemoryLimitHeap) == 0){
         double limitHeap = atof((*inputTerms)[parseIndex+1].c_str());
         Parameters::GetInstance()->SetLimitHeapMemory(limitHeap);
         parseIndex++;
      }
      parseIndex++;   
   }
   return parseIndex;
}

void InputParser::Parse(Molecule* molecule, int argc, char *argv[]) const{

   this->OutputLog(messageStartParseInput);

   // read input
   vector<string> inputTerms = this->GetInputTerms(argc, argv);

   // parse input
   for(int i=0; i<inputTerms.size();i++){

      // theory
      if(inputTerms[i].compare(this->stringTheory) == 0){
         i = this->ParseTheory(&inputTerms, i);
      }

      // molecular configuration
      if(inputTerms[i].compare(this->stringGeometry) == 0){
         i = this->ParseMolecularConfiguration(molecule, &inputTerms, i);
      }

      // ghost atom configuration
      if(inputTerms[i].compare(this->stringGhost) == 0){
         i = this->ParseGhostsConfiguration(molecule, &inputTerms, i);
      }

      // Environmental Point Charges Configuration(EPC)
      if(inputTerms[i].compare(this->stringEpc) == 0){
         i = this->ParseEpcsConfiguration(molecule, &inputTerms, i);
      }

      // scf condition
      if(inputTerms[i].compare(this->stringScf) == 0){
         i = this->ParseConditionsSCF(&inputTerms, i);
      }
      
      // inertia tensor condition
      if(inputTerms[i].compare(this->stringInertiaTensor) == 0){
         i = this->ParseConditionsPrincipalAxes(&inputTerms, i);
      }
      
      // translating condition
      if(inputTerms[i].compare(this->stringTranslate) == 0){
         i = this->ParseConditionsTranslation(&inputTerms, i);
      }
      
      // rotating condition
      if(inputTerms[i].compare(this->stringRotate) == 0){
         i = this->ParseConditionsRotation(&inputTerms, i);
      }
      
      // mo plot condition
      if(inputTerms[i].compare(this->stringMOPlot) == 0){
         i = this->ParseConditionsMOPlot(&inputTerms, i);
      }

      // hole plot condition
      if(inputTerms[i].compare(this->stringHolePlot) == 0){
         i = this->ParseConditionsHolePlot(&inputTerms, i);
      }

      // particle plot condition
      if(inputTerms[i].compare(this->stringParticlePlot) == 0){
         i = this->ParseConditionsParticlePlot(&inputTerms, i);
      }

      // cis condition
      if(inputTerms[i].compare(this->stringCIS) == 0){
         i = this->ParseConditionsCIS(&inputTerms, i);
      }

      // Memory
      if(inputTerms[i].compare(this->stringMemory) == 0){
         i = this->ParseConditionsMemory(&inputTerms, i);
      }

      // MD condition
      if(inputTerms[i].compare(this->stringMD) == 0){
         i = this->ParseConditionsMD(&inputTerms, i);
      }

      // MC condition
      if(inputTerms[i].compare(this->stringMC) == 0){
         i = this->ParseConditionsMC(&inputTerms, i);
      }

      // RPMD condition
      if(inputTerms[i].compare(this->stringRPMD) == 0){
         i = this->ParseConditionsRPMD(&inputTerms, i);
      }

      // NASCO condition
      if(inputTerms[i].compare(this->stringNASCO) == 0){
         i = this->ParseConditionsNASCO(&inputTerms, i);
      }

      // Optimization condition
      if(inputTerms[i].compare(this->stringOptimization) == 0){
         i = this->ParseConditionsOptimization(&inputTerms, i);
      }

      // Frequencies
      if(inputTerms[i].compare(this->stringFrequencies) == 0){
         i = this->ParseConditionsFrequencies(&inputTerms, i);
      }

   }

   // calculate basics and validate conditions
   this->CalcMolecularBasics(molecule);
   this->ValidateVdWConditions();
   this->ValidateEpcConditions(*molecule);
   if(Parameters::GetInstance()->RequiresCIS()){
      this->ValidateCisConditions(*molecule);
   }
   if(Parameters::GetInstance()->RequiresFrequencies()){
      this->ValidateFrequenciesConditions();
   }
   if(Parameters::GetInstance()->GetCurrentSimulation()==MD){
      this->ValidateMdConditions(*molecule);
   }
   else if(Parameters::GetInstance()->GetCurrentSimulation()==MC){
      this->ValidateMcConditions(*molecule);
   }
   else if(Parameters::GetInstance()->GetCurrentSimulation()==RPMD){
      this->ValidateRpmdConditions(*molecule);
   }
   else if(Parameters::GetInstance()->GetCurrentSimulation()==NASCO){
      this->ValidateNascoConditions(*molecule);
   }
   else if(Parameters::GetInstance()->GetCurrentSimulation()==Optimization){
      this->ValidateOptimizationConditions(*molecule);
   }

   // output conditions
   this->OutputMolecularBasics(molecule);
   this->OutputScfConditions();
   this->OutputMemoryConditions();
   if(Parameters::GetInstance()->RequiresCIS()){
      this->OutputCisConditions();
   }
   if(Parameters::GetInstance()->RequiresFrequencies()){
      this->OutputFrequenciesConditions();
   }
   if(Parameters::GetInstance()->RequiresMOPlot()){
      this->OutputMOPlotConditions();
   }
   if(Parameters::GetInstance()->RequiresHolePlot()){
      this->OutputHolePlotConditions();
   }
   if(Parameters::GetInstance()->RequiresParticlePlot()){
      this->OutputParticlePlotConditions();
   }
   if(Parameters::GetInstance()->GetCurrentSimulation()==MD){
      this->OutputMdConditions();
   }
   else if(Parameters::GetInstance()->GetCurrentSimulation()==MC){
      this->OutputMcConditions();
   }
   else if(Parameters::GetInstance()->GetCurrentSimulation()==RPMD){
      this->OutputRpmdConditions();
   }
   else if(Parameters::GetInstance()->GetCurrentSimulation()==NASCO){
      this->OutputNascoConditions();
   }
   else if(Parameters::GetInstance()->GetCurrentSimulation()==Optimization){
      this->OutputOptimizationConditions();
   }

   // output inputs
   this->OutputInputTerms(inputTerms);
   this->OutputLog(messageDoneParseInput);

}

void InputParser::CalcMolecularBasics(Molecule* molecule) const{
   molecule->CalcBasics();
}

void InputParser::ValidateVdWConditions() const{
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   // Validate theory
   if(theory == PM3D || theory == AM1D){
      Parameters::GetInstance()->SetRequiresVdWSCF(true);
      Parameters::GetInstance()->SetVdWScalingFactorSCF();
      Parameters::GetInstance()->SetVdWDampingFactorSCF();
   }
}

void InputParser::ValidateEpcConditions(const Molecule& molecule) const{
   if(molecule.GetEpcVect().empty()){return;}
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   // Validate theory
   if(theory == MNDO || 
      theory == AM1  || 
      theory == AM1D || 
      theory == PM3  || 
      theory == PM3D || 
      theory == PM3PDDG ){
   }
   else{
      stringstream ss;
      ss << this->errorMessageNonValidTheoriesEpc;
      ss << this->errorMessageTheory << TheoryTypeStr(theory) << endl;
      throw MolDSException(ss.str());
   }
}

void InputParser::ValidateCisConditions(const Molecule& molecule) const{

   // direct CIS
   int numberOcc = molecule.GetTotalNumberValenceElectrons()/2;
   int numberVir = molecule.GetTotalNumberAOs() - numberOcc;

   // Validate the number of active occupied orbitals.
   if(numberOcc < Parameters::GetInstance()->GetActiveOccCIS()){
      Parameters::GetInstance()->SetActiveOccCIS(numberOcc);
   }   

   // Validate the number of active virtual orbitals.
   if(numberVir < Parameters::GetInstance()->GetActiveVirCIS()){
      Parameters::GetInstance()->SetActiveVirCIS(numberVir);
   }   

   // Validate the number of calculated excited states.
   int numberSlaterDeterminants = Parameters::GetInstance()->GetActiveOccCIS() 
                                 *Parameters::GetInstance()->GetActiveVirCIS();
   if(!Parameters::GetInstance()->IsDavidsonCIS()){
      Parameters::GetInstance()->SetNumberExcitedStatesCIS(numberSlaterDeterminants);
   }
   else{
      if(numberSlaterDeterminants < Parameters::GetInstance()->GetNumberExcitedStatesCIS()){
         Parameters::GetInstance()->SetNumberExcitedStatesCIS(numberSlaterDeterminants);
      }
      if(numberSlaterDeterminants < Parameters::GetInstance()->GetMaxDimensionsCIS()){
         Parameters::GetInstance()->SetMaxDimensionsCIS(numberSlaterDeterminants);
      }
   }
   
   // Validate the number of printing coefficients of CIS-eigenvector.
   int numPrintCoefficients = Parameters::GetInstance()->GetNumberPrintCoefficientsCIS();
   if(numberSlaterDeterminants < numPrintCoefficients){
      Parameters::GetInstance()->SetNumberPrintCoefficientsCIS(numberSlaterDeterminants);
   }   

   // Validate electronic state for Mulliken population analysis
   if(Parameters::GetInstance()->RequiresMullikenCIS()){
      vector<int>* indecesMulliken = Parameters::GetInstance()->GetElectronicStateIndecesMullikenCIS();
      int numExcitedStates = Parameters::GetInstance()->GetNumberExcitedStatesCIS();
      vector<int>::iterator it=(*indecesMulliken).begin();
      vector<int>::iterator end=(*indecesMulliken).end();
      while(it<end){
         if(numExcitedStates<*it){
            it  = (*indecesMulliken).erase(it);
            end = (*indecesMulliken).end();
         }
         ++it;
      }
   }
}

void InputParser::ValidateMdConditions(const Molecule& molecule) const{
   int groundStateIndex = 0;
   int targetStateIndex = Parameters::GetInstance()->GetElectronicStateIndexMD();
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   // Validate theory
   if(theory == CNDO2 || theory == INDO ){
      stringstream ss;
      ss << this->errorMessageNonValidTheoriesMD;
      ss << this->errorMessageTheory << TheoryTypeStr(theory) << endl;
      throw MolDSException(ss.str());
   }
   // Validate for the excited states dynamics
   if(groundStateIndex < targetStateIndex && !Parameters::GetInstance()->RequiresCIS()){
      Parameters::GetInstance()->SetNumberExcitedStatesCIS(targetStateIndex);
      Parameters::GetInstance()->SetRequiresCIS(true);
      this->ValidateCisConditions(molecule);
   }
   int numberExcitedStatesCIS = Parameters::GetInstance()->GetNumberExcitedStatesCIS();
   if(groundStateIndex < targetStateIndex && numberExcitedStatesCIS < targetStateIndex){
      stringstream ss;
      ss << this->errorMessageNonValidExcitedStatesMD;
      ss << this->errorMessageElecState << targetStateIndex << endl;
      ss << this->errorMessageNumberExcitedStateCIS << numberExcitedStatesCIS << endl;
      throw MolDSException(ss.str());
   } 
}

void InputParser::ValidateMcConditions(const Molecule& molecule) const{
   int groundStateIndex = 0;
   int targetStateIndex = Parameters::GetInstance()->GetElectronicStateIndexMC();
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   // CNDO2 and INDO do not support excited states.
   if((theory == CNDO2 || theory == INDO) && groundStateIndex < targetStateIndex){
      stringstream ss;
      ss << this->errorMessageNonValidExcitedStatesMC;
      ss << this->errorMessageElecState << targetStateIndex << endl;
      ss << this->errorMessageTheory << TheoryTypeStr(theory) << endl;
      throw MolDSException(ss.str());
   }
   // Validate for the excited states dynamics
   if(groundStateIndex < targetStateIndex && !Parameters::GetInstance()->RequiresCIS()){
      Parameters::GetInstance()->SetNumberExcitedStatesCIS(targetStateIndex);
      Parameters::GetInstance()->SetRequiresCIS(true);
      this->ValidateCisConditions(molecule);
   }
   int numberExcitedStatesCIS = Parameters::GetInstance()->GetNumberExcitedStatesCIS();
   if(numberExcitedStatesCIS < targetStateIndex){
      stringstream ss;
      ss << this->errorMessageNonValidExcitedStatesMC;
      ss << this->errorMessageElecState << targetStateIndex << endl;
      ss << this->errorMessageNumberExcitedStateCIS << numberExcitedStatesCIS << endl;
      throw MolDSException(ss.str());
   }
}

void InputParser::ValidateRpmdConditions(const Molecule& molecule) const{
   int groundStateIndex = 0;
   int targetStateIndex = Parameters::GetInstance()->GetElectronicStateIndexRPMD();
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   // Validate theory
   if(theory == CNDO2 || theory == INDO ){
      stringstream ss;
      ss << this->errorMessageNonValidTheoriesRPMD;
      ss << this->errorMessageTheory << TheoryTypeStr(theory) << endl;
      throw MolDSException(ss.str());
   }
   // Validate for the excited states dynamics
   if(groundStateIndex < targetStateIndex && !Parameters::GetInstance()->RequiresCIS()){
      Parameters::GetInstance()->SetNumberExcitedStatesCIS(targetStateIndex);
      Parameters::GetInstance()->SetRequiresCIS(true);
      this->ValidateCisConditions(molecule);
   }
   int numberExcitedStatesCIS = Parameters::GetInstance()->GetNumberExcitedStatesCIS();
   if(groundStateIndex < targetStateIndex && numberExcitedStatesCIS < targetStateIndex){
      stringstream ss;
      ss << this->errorMessageNonValidExcitedStatesRPMD;
      ss << this->errorMessageElecState << targetStateIndex << endl;
      ss << this->errorMessageNumberExcitedStateCIS << numberExcitedStatesCIS << endl;
      throw MolDSException(ss.str());
   } 
}

void InputParser::ValidateNascoConditions(const Molecule& molecule) const{
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   // Validate theory
   if(theory == CNDO2 || theory == INDO || theory == ZINDOS){
      stringstream ss;
      ss << this->errorMessageNonValidTheoriesNASCO;
      ss << this->errorMessageTheory << TheoryTypeStr(theory) << endl;
      throw MolDSException(ss.str());
   }
   // Validate requirement of CIS
   if(!Parameters::GetInstance()->RequiresCIS()){
      Parameters::GetInstance()->SetRequiresCIS(true);
      this->ValidateCisConditions(molecule);
   }
   // Validate number of excited states
   int numberExcitedStatesCIS   = Parameters::GetInstance()->GetNumberExcitedStatesCIS();
   int numberExcitedStatesNASCO = Parameters::GetInstance()->GetNumberElectronicStatesNASCO() - 1;
   if(numberExcitedStatesCIS < numberExcitedStatesNASCO){
      stringstream ss;
      ss << this->errorMessageNonValidNumberExcitedStatesNASCO;
      ss << this->errorMessageNumberElectronicStatesNASCO << numberExcitedStatesNASCO+1 << endl;
      ss << this->errorMessageNumberExcitedStateCIS << numberExcitedStatesCIS << endl;
      throw MolDSException(ss.str());
   } 
   // Validate initial electronic eigenstate
   int initialElectronicStateNASCO = Parameters::GetInstance()->GetInitialElectronicStateNASCO();
   int numberElectronicStatesNASCO = Parameters::GetInstance()->GetNumberElectronicStatesNASCO();
   if(numberElectronicStatesNASCO<=initialElectronicStateNASCO){
      stringstream ss;
      ss << this->errorMessageNonValidInitialElectronicStateNASCO;
      ss << this->errorMessageNumberElectronicStatesNASCO << numberElectronicStatesNASCO << endl;
      ss << this->errorMessageInitialElectronicStateNASCO << initialElectronicStateNASCO << endl;
      throw MolDSException(ss.str());
   }
}

void InputParser::ValidateOptimizationConditions(const Molecule& molecule) const{
   int groundStateIndex = 0;
   int targetStateIndex = Parameters::GetInstance()->GetElectronicStateIndexOptimization();
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   // Validate theory
   if(theory == CNDO2 || theory == INDO ){
      stringstream ss;
      ss << this->errorMessageNonValidTheoriesOptimization;
      ss << this->errorMessageElecState << targetStateIndex << endl;
      ss << this->errorMessageTheory << TheoryTypeStr(theory) << endl;
      throw MolDSException(ss.str());
   }
   // Validate for the excited states dynamics
   if(groundStateIndex < targetStateIndex && !Parameters::GetInstance()->RequiresCIS()){
      Parameters::GetInstance()->SetNumberExcitedStatesCIS(targetStateIndex);
      Parameters::GetInstance()->SetRequiresCIS(true);
      this->ValidateCisConditions(molecule);
   }
   int numberExcitedStatesCIS = Parameters::GetInstance()->GetNumberExcitedStatesCIS();
   if(groundStateIndex < targetStateIndex && numberExcitedStatesCIS < targetStateIndex){
      stringstream ss;
      ss << this->errorMessageNonValidExcitedStatesOptimization;
      ss << this->errorMessageElecState << targetStateIndex << endl;
      ss << this->errorMessageNumberExcitedStateCIS << numberExcitedStatesCIS << endl;
      throw MolDSException(ss.str());
   } 
}

void InputParser::ValidateFrequenciesConditions() const{
   // validate theory
   TheoryType theory = Parameters::GetInstance()->GetCurrentTheory();
   if(theory == CNDO2 || theory == INDO || theory == ZINDOS){
      stringstream ss;
      ss << this->errorMessageNonValidTheoryFrequencies;
      ss << this->errorMessageTheory << TheoryTypeStr(theory) << endl;
      throw MolDSException(ss.str());
   }
   // validate electronic state
   int groundStateIndex = 0;
   int targetStateIndex = Parameters::GetInstance()->GetElectronicStateIndexFrequencies();
   if(groundStateIndex < targetStateIndex){
      stringstream ss;
      ss << this->errorMessageNonValidElectronicStateFrequencies;
      ss << this->errorMessageElecState << targetStateIndex << endl;
      throw MolDSException(ss.str());
   } 
}

void InputParser::OutputMolecularBasics(Molecule* molecule) const{
   molecule->OutputTotalNumberAtomsAOsValenceelectrons();
   molecule->OutputConfiguration();
   molecule->OutputXyzCOM();
   molecule->OutputXyzCOC();
   molecule->OutputEpcs();
}

void InputParser::OutputScfConditions() const{
   this->OutputLog(this->messageScfConditions);
   this->OutputLog(boost::format("%s%d\n") % this->messageScfMaxIterations.c_str() 
                                           % Parameters::GetInstance()->GetMaxIterationsSCF());
   this->OutputLog(boost::format("%s%e\n") % this->messageScfRmsDensity.c_str() 
                                           % Parameters::GetInstance()->GetThresholdSCF());
   this->OutputLog(boost::format("%s%e\n") % this->messageScfDampingThresh.c_str()
                                           % Parameters::GetInstance()->GetDampingThreshSCF());
   this->OutputLog(boost::format("%s%e\n") % this->messageScfDampingWeight.c_str()
                                           % Parameters::GetInstance()->GetDampingWeightSCF());
   this->OutputLog(boost::format("%s%d\n") % this->messageScfDiisNumErrorVect.c_str() 
                                           % Parameters::GetInstance()->GetDiisNumErrorVectSCF());
   this->OutputLog(boost::format("%s%e\n") % this->messageScfDiisStartError.c_str() 
                                           % Parameters::GetInstance()->GetDiisStartErrorSCF());
   this->OutputLog(boost::format("%s%e\n") % this->messageScfDiisEndError.c_str() 
                                           % Parameters::GetInstance()->GetDiisEndErrorSCF());
   this->OutputLog(this->messageScfVdW);
   if(Parameters::GetInstance()->RequiresVdWSCF()){
      this->OutputLog(boost::format("%s\n") % this->stringYES.c_str());
      this->OutputLog(boost::format("%s%lf\n") % this->messageScfVdWScalingFactor.c_str() 
                                               % Parameters::GetInstance()->GetVdWScalingFactorSCF());
      this->OutputLog(boost::format("%s%lf\n") % this->messageScfVdWDampingFactor.c_str() 
                                               % Parameters::GetInstance()->GetVdWDampingFactorSCF());
   }
   else{
      this->OutputLog(boost::format("%s\n") % this->stringNO.c_str());
   }
   this->OutputLog("\n");
}

void InputParser::OutputMemoryConditions() const{
   this->OutputLog(this->messageMemoryConditions);
   this->OutputLog(boost::format("%s%e%s\n") % this->messageMemoryLimitHeap.c_str() 
                                           % Parameters::GetInstance()->GetLimitHeapMemory()
                                           % this->messageMB.c_str());
   this->OutputLog("\n");
}

void InputParser::OutputCisConditions() const{
   this->OutputLog(this->messageCisConditions);

   this->OutputLog(boost::format("%s%d\n") % this->messageCisNumberActiveOcc.c_str() 
                                           % Parameters::GetInstance()->GetActiveOccCIS());
   this->OutputLog(boost::format("%s%d\n") % this->messageCisNumberActiveVir.c_str() 
                                           % Parameters::GetInstance()->GetActiveVirCIS());
   this->OutputLog(boost::format("%s%d\n") % this->messageCisNumberExcitedStates.c_str() 
                                           % Parameters::GetInstance()->GetNumberExcitedStatesCIS());
   this->OutputLog(boost::format("%s%d\n") % this->messageCisNumPrintCoefficients.c_str() 
                                           % Parameters::GetInstance()->GetNumberPrintCoefficientsCIS());

   this->OutputLog(this->messageCisDavidson);
   if(Parameters::GetInstance()->IsDavidsonCIS()){
      this->OutputLog(boost::format("%s\n") % this->stringYES.c_str());
      this->OutputLog(boost::format("%s%d\n") % this->messageCisMaxIterations.c_str() 
                                              % Parameters::GetInstance()->GetMaxIterationsCIS());
      this->OutputLog(boost::format("%s%d\n") % this->messageCisMaxDimensions.c_str() 
                                              % Parameters::GetInstance()->GetMaxDimensionsCIS());
      this->OutputLog(boost::format("%s%e\n") % this->messageCisNormTolerance.c_str() 
                                              % Parameters::GetInstance()->GetNormToleranceCIS());
   }
   else{
      this->OutputLog(boost::format("%s\n") % this->stringNO.c_str());
   }

   this->OutputLog(this->messageCisExcitonEnergies);
   if(Parameters::GetInstance()->RequiresExcitonEnergiesCIS()){
      this->OutputLog(this->stringYES);
   }
   else{
      this->OutputLog(this->stringNO);
   }
   this->OutputLog("\n");

   this->OutputLog(this->messageCisAllTransitionDipoleMoments);
   if(Parameters::GetInstance()->RequiresAllTransitionDipoleMomentsCIS()){
      this->OutputLog(this->stringYES);
   }
   else{
      this->OutputLog(this->stringNO);
   }
   this->OutputLog("\n");

   if(Parameters::GetInstance()->RequiresMullikenCIS()){
      vector<int>* indeces = Parameters::GetInstance()->GetElectronicStateIndecesMullikenCIS();
      for(int i=0; i<indeces->size(); i++){
         this->OutputLog(boost::format("%s%d\n") % this->messageCisMulliken.c_str()
                                                 % (*indeces)[i]);
      }
      this->OutputLog("\n");
   }
   
   this->OutputLog("\n");
}

void InputParser::OutputMdConditions() const{
   this->OutputLog(this->messageMdConditions);

   this->OutputLog(boost::format("%s%d\n") % this->messageMdElecState.c_str() 
                                           % Parameters::GetInstance()->GetElectronicStateIndexMD());
   this->OutputLog(boost::format("%s%d\n") % this->messageMdTotalSteps.c_str() 
                                           % Parameters::GetInstance()->GetTotalStepsMD());
   this->OutputLog(boost::format("%s%lf%s\n") % this->messageMdTimeWidth.c_str() 
                                              % (Parameters::GetInstance()->GetTimeWidthMD()/Parameters::GetInstance()->GetFs2AU())
                                              % this->messageFs.c_str());

   this->OutputLog("\n");
}

void InputParser::OutputMcConditions() const{
   this->OutputLog(this->messageMcConditions);

   this->OutputLog(boost::format("%s%d\n") % this->messageMcElecState.c_str() 
                                           % Parameters::GetInstance()->GetElectronicStateIndexMC());
   this->OutputLog(boost::format("%s%d\n") % this->messageMcTotalSteps.c_str() 
                                           % Parameters::GetInstance()->GetTotalStepsMC());
   this->OutputLog(boost::format("%s%lf%s\n") % this->messageMcTemperature.c_str() 
                                              % Parameters::GetInstance()->GetTemperatureMC()
                                              % this->messageK.c_str());
   this->OutputLog(boost::format("%s%lf%s\n") % this->messageMcStepWidth.c_str() 
                                              % (Parameters::GetInstance()->GetStepWidthMC()/Parameters::GetInstance()->GetAngstrom2AU())
                                              % this->messageAngst.c_str());
   this->OutputLog(boost::format("%s%lu\n") % this->messageMcSeed.c_str() 
                                            % Parameters::GetInstance()->GetSeedMC());

   this->OutputLog("\n");
}

void InputParser::OutputRpmdConditions() const{
   this->OutputLog(this->messageRpmdConditions);

   this->OutputLog(boost::format("%s%d\n") % this->messageRpmdElecState.c_str() 
                                           % Parameters::GetInstance()->GetElectronicStateIndexRPMD());
   this->OutputLog(boost::format("%s%d\n") % this->messageRpmdNumElecStates.c_str() 
                                           % Parameters::GetInstance()->GetNumberElectronicStatesRPMD());
   this->OutputLog(boost::format("%s%d\n") % this->messageRpmdTotalSteps.c_str() 
                                           % Parameters::GetInstance()->GetTotalStepsRPMD());
   this->OutputLog(boost::format("%s%lf%s\n") % this->messageRpmdTemperature.c_str() 
                                              % Parameters::GetInstance()->GetTemperatureRPMD()
                                              % this->messageK.c_str());
   this->OutputLog(boost::format("%s%lf%s\n") % this->messageRpmdTimeWidth.c_str() 
                                              % (Parameters::GetInstance()->GetTimeWidthRPMD()/Parameters::GetInstance()->GetFs2AU()) 
                                              % this->messageFs.c_str());
   this->OutputLog(boost::format("%s%d\n") % this->messageRpmdNumBeads.c_str() 
                                           % Parameters::GetInstance()->GetNumberBeadsRPMD());
   this->OutputLog(boost::format("%s%lu\n") % this->messageRpmdSeed.c_str() 
                                            % Parameters::GetInstance()->GetSeedRPMD());

   this->OutputLog("\n");
}

void InputParser::OutputNascoConditions() const{
   this->OutputLog(this->messageNascoConditions);

   this->OutputLog(boost::format("%s%d\n") % this->messageNascoTotalSteps.c_str() 
                                           % Parameters::GetInstance()->GetTotalStepsNASCO());
   this->OutputLog(boost::format("%s%d\n") % this->messageNascoNumElecStates.c_str() 
                                           % Parameters::GetInstance()->GetNumberElectronicStatesNASCO());
   this->OutputLog(boost::format("%s%d\n") % this->messageNascoInitialElecState.c_str() 
                                           % Parameters::GetInstance()->GetInitialElectronicStateNASCO());
   this->OutputLog(boost::format("%s%lf%s\n") % this->messageNascoTimeWidth.c_str() 
                                              % (Parameters::GetInstance()->GetTimeWidthNASCO()/Parameters::GetInstance()->GetFs2AU()) 
                                              % this->messageFs.c_str());
   this->OutputLog(boost::format("%s%lu\n") % this->messageNascoSeed.c_str() 
                                            % Parameters::GetInstance()->GetSeedNASCO());

   this->OutputLog("\n");
}

void InputParser::OutputOptimizationConditions() const{
   this->OutputLog(this->messageOptimizationConditions);

   this->OutputLog(boost::format("%s%s\n") % this->messageOptimizationMethod.c_str() 
                                           % OptimizationMethodTypeStr(Parameters::GetInstance()->
                                                                        GetMethodOptimization()));
   this->OutputLog(boost::format("%s%d\n") % this->messageOptimizationTotalSteps.c_str() 
                                           % Parameters::GetInstance()->GetTotalStepsOptimization());
   this->OutputLog(boost::format("%s%d\n") % this->messageOptimizationElecState.c_str() 
                                           % Parameters::GetInstance()->GetElectronicStateIndexOptimization());
   this->OutputLog(boost::format("%s%lf\n") % this->messageOptimizationMaxGradient.c_str() 
                                            % Parameters::GetInstance()->GetMaxGradientOptimization());
   this->OutputLog(boost::format("%s%lf\n") % this->messageOptimizationRmsGradient.c_str() 
                                            % Parameters::GetInstance()->GetRmsGradientOptimization());

   switch(Parameters::GetInstance()->GetMethodOptimization()){
      case ConjugateGradientMethod:
      case SteepestDescentMethod:
         this->OutputLog(boost::format("%s%lf%s\n") % this->messageOptimizationTimeWidth.c_str() 
                                                    % (Parameters::GetInstance()->GetTimeWidthOptimization()/Parameters::GetInstance()->GetFs2AU())
                                                    % this->messageFs.c_str());
         break;
			case BFGSMethod:
         this->OutputLog(boost::format("%s%lf\n") % this->messageOptimizationInitialTrustRadius.c_str()
                                                  % Parameters::GetInstance()->GetInitialTrustRadiusOptimization());
         this->OutputLog(boost::format("%s%lf\n") % this->messageOptimizationMaxNormStep.c_str()
                                                  % Parameters::GetInstance()->GetMaxNormStepOptimization());
			case GEDIISMethod:
         this->OutputLog(boost::format("%s%lf\n") % this->messageOptimizationInitialTrustRadius.c_str()
                                                  % Parameters::GetInstance()->GetInitialTrustRadiusOptimization());
         this->OutputLog(boost::format("%s%lf\n") % this->messageOptimizationMaxNormStep.c_str()
                                                  % Parameters::GetInstance()->GetMaxNormStepOptimization());
      default:
         break;
   }

   this->OutputLog("\n");
}

void InputParser::OutputFrequenciesConditions() const{
   this->OutputLog(this->messageFrequenciesConditions);
   this->OutputLog(boost::format("%s%d\n") % this->messageFrequenciesElecState.c_str() 
                                           % Parameters::GetInstance()->GetElectronicStateIndexFrequencies());
   this->OutputLog("\n");
}

void InputParser::OutputMOPlotConditions() const{
   this->OutputLog(this->messageMOPlotConditions);
   vector<int>* moIndeces = Parameters::GetInstance()->GetIndecesMOPlot();
   for(int i=0; i<moIndeces->size(); i++){
      this->OutputLog(boost::format("%s%d\n") % this->messageMOPlotIndex.c_str() 
                                              % (*moIndeces)[i]);
   }
   const int* gridNum = Parameters::GetInstance()->GetGridNumberMOPlot();
   this->OutputLog(boost::format("%s%d %d %d\n") % this->messageMOPlotGridNumber.c_str() 
                                                 % gridNum[XAxis] 
                                                 % gridNum[YAxis]
                                                 % gridNum[ZAxis]);
   const double* frameLength = Parameters::GetInstance()->GetFrameLengthMOPlot();
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   this->OutputLog(boost::format("%s%e %e %e\n") % this->messageMOPlotFrameLength.c_str() 
                                                 % (frameLength[XAxis]/ang2AU) 
                                                 % (frameLength[YAxis]/ang2AU)
                                                 % (frameLength[ZAxis]/ang2AU));
   this->OutputLog(boost::format("%s%s\n") % this->messageMOPlotFilePrefix.c_str() 
                                           % Parameters::GetInstance()->GetFileNamePrefixMOPlot().c_str());

   this->OutputLog("\n");
}

void InputParser::OutputHolePlotConditions() const{
   this->OutputLog(this->messageHolePlotConditions);
   vector<int>* moIndeces = Parameters::GetInstance()->GetElecIndecesHolePlot();
   for(int i=0; i<moIndeces->size(); i++){
      this->OutputLog(boost::format("%s%d\n") % this->messageHolePlotElecIndex.c_str() 
                                              % (*moIndeces)[i]);
   }
   const int* gridNum = Parameters::GetInstance()->GetGridNumberHolePlot();
   this->OutputLog(boost::format("%s%d %d %d\n") % this->messageHolePlotGridNumber.c_str() 
                                                 % gridNum[XAxis] 
                                                 % gridNum[YAxis]
                                                 % gridNum[ZAxis]);
   const double* frameLength = Parameters::GetInstance()->GetFrameLengthHolePlot();
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   this->OutputLog(boost::format("%s%e %e %e\n") % this->messageHolePlotFrameLength.c_str() 
                                                 % (frameLength[XAxis]/ang2AU) 
                                                 % (frameLength[YAxis]/ang2AU)
                                                 % (frameLength[ZAxis]/ang2AU));
   this->OutputLog(boost::format("%s%s\n") % this->messageHolePlotFilePrefix.c_str() 
                                           % Parameters::GetInstance()->GetFileNamePrefixHolePlot().c_str());

   this->OutputLog("\n");
}

void InputParser::OutputParticlePlotConditions() const{
   this->OutputLog(this->messageParticlePlotConditions);
   const vector<int>* moIndeces = Parameters::GetInstance()->GetElecIndecesParticlePlot();
   for(int i=0; i<moIndeces->size(); i++){
      this->OutputLog(boost::format("%s%d\n") % this->messageParticlePlotElecIndex.c_str() 
                                              % (*moIndeces)[i]);
   }
   const int* gridNum = Parameters::GetInstance()->GetGridNumberParticlePlot();
   this->OutputLog(boost::format("%s%d %d %d\n") % this->messageParticlePlotGridNumber.c_str() 
                                                 % gridNum[XAxis] 
                                                 % gridNum[YAxis]
                                                 % gridNum[ZAxis]);
   const double* frameLength = Parameters::GetInstance()->GetFrameLengthParticlePlot();
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   this->OutputLog(boost::format("%s%e %e %e\n") % this->messageParticlePlotFrameLength.c_str() 
                                                 % (frameLength[XAxis]/ang2AU) 
                                                 % (frameLength[YAxis]/ang2AU)
                                                 % (frameLength[ZAxis]/ang2AU));
   this->OutputLog(boost::format("%s%s\n") % this->messageParticlePlotFilePrefix.c_str() 
                                           % Parameters::GetInstance()->GetFileNamePrefixParticlePlot().c_str());

   this->OutputLog("\n");
}

void InputParser::OutputInputTerms(vector<string> inputTerms) const{
   // output input terms
   this->OutputLog(this->messageInputTerms);
   for(int i=0; i<inputTerms.size();i++){
      this->OutputLog((inputTerms[i] + " | "));
      if(i%10 == 9){
         this->OutputLog("\n");
      }
   }
   this->OutputLog("\n\n");
}

/****
 *
 *  # or // are treated as comment out 
 *
 ****/
bool InputParser::IsCommentOut(string tempStr) const{
   string str = Utilities::TrimString(tempStr);

   string commentPrefix1 = "#";
   string prefix1;
   if(str.length()>=1){
      prefix1 += str.data()[0];
   }

   string commentPrefix2 = "//";
   string prefix2;
   if(str.length()>=2){
      prefix2 += str.data()[0];
      prefix2 += str.data()[1];
   }

   return 0==prefix1.compare(commentPrefix1) || 0==prefix2.compare(commentPrefix2) ;
}

}





