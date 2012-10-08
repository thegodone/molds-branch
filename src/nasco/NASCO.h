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
   std::string messageinitialConditionNASCO;
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
   std::string messageTime;
   std::string errorMessageNotEnebleTheoryType;
   std::string errorMessageTheoryType;
   std::vector<MolDS_base::TheoryType> enableTheoryTypes;
   void CheckEnableTheoryType(MolDS_base::TheoryType theoryType);
   void SetMessages();
   void SetEnableTheoryTypes();
   void UpdateMomenta(MolDS_base::Molecule& molecule, double const* const* matrixForce, double dt) const;
   void SynchronousMolecularConfiguration(MolDS_base::Molecule& target, 
                                          const MolDS_base::Molecule& refference) const;
   void OutputEnergies(const MolDS_base::ElectronicStructure& electronicStructure, double initialEnergy, const MolDS_base::Molecule& molecule);
   double OutputEnergies(const MolDS_base::ElectronicStructure& electronicStructure, const MolDS_base::Molecule& molecule);
   void MallocOverlapsDifferentMolecules(double*** overlapAOs,
                                         double*** overlapMOs, 
                                         double*** overlapESs, 
                                         const MolDS_base::Molecule& molecule) const;
   void FreeOverlapsDifferentMolecules(double*** overlapAOs,
                                       double*** overlapMOs, 
                                       double*** overlapESs, 
                                       const MolDS_base::Molecule& molecule) const;
};

}
#endif



