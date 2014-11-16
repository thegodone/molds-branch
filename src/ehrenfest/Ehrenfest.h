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
#ifndef INCLUDED_EHRENFEST
#define INCLUDED_EHRENFEST
namespace MolDS_ehrenfest{

/***
 *  Velocty Verlet is used for nuclear motion.
 *  Unpublished algorythm by MF is used for electronic motion.
 */
class Ehrenfest : public MolDS_base::PrintController{
public:
   Ehrenfest();
   ~Ehrenfest();
   void SetMolecule(MolDS_base::Molecule* molecule);
   void DoEhrenfest();
private:
   std::string messageinitialConditionEhrenfest;
   std::string messageStartEhrenfest;
   std::string messageEndEhrenfest;
   std::string messageStartStepEhrenfest;
   std::string messageEndStepEhrenfest;
   std::string messageEnergies;
   std::string messageEnergiesVdW;
   std::string messageElecPopuTitle;
   std::string messageElecPopulations;
   std::string messageElecEigenEnergyTitle;
   std::string messageElecEigenEnergyAU;
   std::string messageElecEigenEnergyEV;
   std::string messageEnergiesTitle;
   std::string messageCoreKineticEnergy;
   std::string messageCoreRepulsionEnergy;
   std::string messageVdWCorrectionEnergy;
   std::string messageElectronicEnergy;
   std::string messageTotalEnergy;
   std::string messageErrorEnergy;
   std::string messageTime;
   std::string errorMessageNotEnebleTheoryType;
   std::string errorMessageTheoryType;
   MolDS_base::Molecule*               molecule;
   std::vector<int>*                   elecStates;
   std::complex<double>*               superpositionCoeff;
   std::vector<MolDS_base::TheoryType> enableTheoryTypes;
   void   CheckEnableTheoryType(MolDS_base::TheoryType theoryType);
   void   SetMessages();
   void   SetEnableTheoryTypes();
   double GetElecNorm();
   void   CalcMeanForce(double** matrixMeanForce, double const* const* const* matrixForce, double elecNorm);
   void   UpdateMomenta    (const MolDS_base::Molecule& molecule, double const* const* matrixForce, double dt) const;
   void   UpdateCoordinates(      MolDS_base::Molecule& molecule, double dt) const;
   void   OutputEnergies(const MolDS_base::ElectronicStructure&, double initialEnergy, double elecNorm);
   double OutputEnergies(const MolDS_base::ElectronicStructure& electronicStructure, double elecNorm);
   void   MallocOverlapsDifferentMolecules(double*** overlapAOs,
                                           double*** overlapMOs, 
                                           double*** overlapSingleSDs,
                                           double*** overlapESs, 
                                           const MolDS_base::Molecule& molecule) const;
   void   FreeOverlapsDifferentMolecules(double*** overlapAOs,
                                         double*** overlapMOs, 
                                         double*** overlapSingleSDs,
                                         double*** overlapESs, 
                                         const MolDS_base::Molecule& molecule) const;
};

}
#endif



