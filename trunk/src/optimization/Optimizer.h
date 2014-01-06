//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   // 
// Copyright (C) 2012-2014 Katsuhiko Nishimra                             // 
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
#ifndef INCLUDED_OPTIMIZER
#define INCLUDED_OPTIMIZER
namespace MolDS_optimization{

class Optimizer : public MolDS_base::PrintController{
public:
   Optimizer();
   virtual ~Optimizer();
   void Optimize(MolDS_base::Molecule& molecule);
protected:
   std::string errorMessageNotEnebleTheoryType;
   std::string errorMessageGeometyrOptimizationNotConverged;
   std::string messageLineSearchSteps;
   virtual void SetMessages();
   void UpdateMolecularCoordinates(MolDS_base::Molecule& molecule, double const* const* matrixForce, double dt) const;
   void UpdateMolecularCoordinates(MolDS_base::Molecule& molecule, double const* const* matrixForce) const;
   void UpdateElectronicStructure(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure, 
                                  MolDS_base::Molecule& molecule,
                                  bool requireGuess, 
                                  bool printsLogs) const;
   bool SatisfiesConvergenceCriterion(double const* const* matrixForce, 
                                      const MolDS_base::Molecule& molecule,
                                      double oldEnergy,
                                      double currentEnergy,
                                      double maxGradientThreshold,
                                      double rmsGradientThreshold) const;
   void OutputMoleculeElectronicStructure(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure, 
                                          MolDS_base::Molecule& molecule,
                                          bool printsLogs) const;
   void LineSearch(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                   MolDS_base::Molecule& molecule,
                   double &lineSearchCurrentEnergy,
                   double const* const* matrixForce,
                   int elecState,
                   double dt) const;
private:
   std::string errorMessageTheoryType;
   std::string errorMessageTotalSteps;
   std::string messageGeometyrOptimizationMetConvergence;
   std::string messageStartGeometryOptimization;
   std::string messageEndGeometryOptimization;
   std::string messageReducedTimeWidth;
   std::string messageOptimizationLog;
   std::string messageEnergyDifference;
   std::string messageMaxGradient;
   std::string messageRmsGradient;
   std::string messageAu;
   std::vector<MolDS_base::TheoryType> enableTheoryTypes;
   void SetEnableTheoryTypes();
   void CheckEnableTheoryType(MolDS_base::TheoryType theoryType) const;
   void ClearMolecularMomenta(MolDS_base::Molecule& molecule) const;
   virtual void SearchMinimum(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                              MolDS_base::Molecule& molecule,
                              double* lineSearchedEnergy,
                              bool* obainesOptimizedStructure) const = 0;
};

}
#endif



