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
#ifndef INCLUDED_MC
#define INCLUDED_MC
namespace MolDS_mc{

/***
 *  Canonical MC is only implemented
 */
class MC : public MolDS_base::PrintController{
public:
   MC();
   ~MC();
   void SetMolecule(MolDS_base::Molecule* molecule);
   void DoMC();
   void DoMC(int totalSteps, int elecState, double temperature, double stepWidth, unsigned long seed);
private:
   std::string messageinitialConditionMC;
   std::string messageStartMC;
   std::string messageEndMC;
   std::string messageStartStepMC;
   std::string messageEndStepMC;
   std::string messageEnergies;
   std::string messageEnergiesTitle;
   std::string messageCoreRepulsionEnergy;
   std::string messageVdWCorrectionEnergy;
   std::string messageElectronicEnergy;
   std::string messageElectronicEnergyVdW;
   std::string messageTotalEnergy;
   std::string messageTransitionRate;
   std::string errorMessageNotEnebleExcitedTheoryType;
   std::string errorMessageTheoryType;
   MolDS_base::Molecule* molecule;
   void SetMessages();
   void CreateTrialConfiguration(MolDS_base::Molecule* trial,
                                 const MolDS_base::Molecule& current,
                                 boost::variate_generator<
                                    boost::mt19937&,
                                    boost::uniform_real<>
                                 > (*realRand),
                                 double dr) const;
   bool UsesTrial(const MolDS_base::ElectronicStructure& currentES, 
                  const MolDS_base::ElectronicStructure& trialES,
                  int elecState,
                  boost::variate_generator<
                     boost::mt19937&,
                     boost::uniform_real<>
                  > (*realRand),
                  double temperature) const;
   void OutputMolecule(const MolDS_base::ElectronicStructure& electronicStructure,
                       const MolDS_base::Molecule& molecule,
                       int elecState) const;
   void OutputEnergies(const MolDS_base::ElectronicStructure& electronicStructure,
                       int elecState) const;
};

}
#endif



