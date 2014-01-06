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
#ifndef INCLUDED_NASCO
#define INCLUDED_NASCO
namespace MolDS_nasco{

/***
 *  Non-Adiabatic SemiClassical kernel based on Overlap integrals
 *  see [F_2011]
 */
class NASCO : public MolDS_base::PrintController{
public:
   NASCO();
   ~NASCO();
   void DoNASCO(MolDS_base::Molecule& molecule);
private:
   std::string messageInitialConditionNASCO;
   std::string messageStartNASCO;
   std::string messageEndNASCO;
   std::string messageStartStepNASCO;
   std::string messageEndStepNASCO;
   std::string messageEnergies;
   std::string messageEnergiesTitle;
   std::string messageCoreKineticEnergy;
   std::string messageCoreRepulsionEnergy;
   std::string messageVdWCorrectionEnergy;
   std::string messageElectronicEnergy;
   std::string messageElectronicEnergyVdW;
   std::string messageTotalEnergy;
   std::string messageErrorEnergy;
   std::string messageElectronicState;
   std::string messageTime;
   std::string errorMessageNotEnebleTheoryType;
   std::string errorMessageTheoryType;
   std::vector<MolDS_base::TheoryType> enableTheoryTypes;
   void CheckEnableTheoryType(MolDS_base::TheoryType theoryType);
   void SetMessages();
   void SetEnableTheoryTypes();
   void UpdateMomenta(MolDS_base::Molecule& molecule, 
                      double const* const* matrixForce, 
                      const double dt) const;
   void UpdateCoordinates(MolDS_base::Molecule& tmpMolecule, 
                          const MolDS_base::Molecule& molecule,
                          const double dt) const;
   void DecideNextElecState(int* elecState, 
                            int* nonAdiabaticPhaseIndex,
                            int numElecStates,
                            double const* const* overlapESs,
                            boost::variate_generator<
                               boost::mt19937&,
                               boost::uniform_real<>
                            > (*realRand)) const;
   void OutputEnergies(const MolDS_base::ElectronicStructure& electronicStructure, 
                       const MolDS_base::Molecule& molecule, 
                       const double initialEnergy, 
                       const int elecState) const;
   double OutputEnergies(const MolDS_base::ElectronicStructure& electronicStructure, 
                         const MolDS_base::Molecule& molecule, 
                         const int elecState) const;
   void MallocOverlapsDifferentMolecules(double*** overlapAOs,
                                         double*** overlapMOs, 
                                         double*** overlapSingleSDs,
                                         double*** overlapESs, 
                                         const MolDS_base::Molecule& molecule) const;
   void FreeOverlapsDifferentMolecules(double*** overlapAOs,
                                       double*** overlapMOs, 
                                       double*** overlapSingleSDs,
                                       double*** overlapESs, 
                                       const MolDS_base::Molecule& molecule) const;
};

}
#endif



