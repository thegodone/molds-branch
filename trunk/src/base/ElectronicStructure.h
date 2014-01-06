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
#ifndef INCLUDED_ELECTRONICSTRUCTURE
#define INCLUDED_ELECTRONICSTRUCTURE
namespace MolDS_base{

class ElectronicStructure : public MolDS_base::PrintController{
public:
   virtual ~ElectronicStructure(){};
   virtual MolDS_base::TheoryType GetTheoryType() const = 0;
   virtual void SetMolecule(MolDS_base::Molecule* molecule) = 0;
   virtual void DoSCF(bool requiresGuess=true) = 0;
   virtual void OutputSCFResults() const = 0;
   virtual double const* const* GetFockMatrix() const = 0;
   virtual double const*        GetEnergiesMO() const = 0;
   virtual void DoCIS() = 0;
   virtual void OutputCISResults() const = 0;
   virtual double const* const* GetMatrixCIS() const = 0;
   virtual double const*        GetExcitedEnergies() const = 0;
   virtual double const* const*        GetForce(int elecState) = 0;
   virtual double const* const* const* GetForce(const std::vector<int>& elecStates) = 0;
   virtual double GetElectronicEnergy(int elecState) const = 0;
   virtual double GetCoreRepulsionEnergy() const = 0;
   virtual double GetVdWCorrectionEnergy() const = 0;
   virtual void CalcOverlapAOsWithAnotherConfiguration(double** overlapAOs, 
                                                       const MolDS_base::Molecule& lhsMoledule) const = 0;
   virtual void CalcOverlapMOsWithAnotherElectronicStructure(double** overlapMOs, 
                                                             double const* const* overlapAOs,
                                                             const MolDS_base::ElectronicStructure& lhsElectronicStructure) const = 0;
   virtual void CalcOverlapSingletSDsWithAnotherElectronicStructure(double** overlapSingletSDs, 
                                                                    double const* const* overlapMOs) const = 0;
   virtual void CalcOverlapESsWithAnotherElectronicStructure(double** overlapESs, 
                                                             double const* const* overlapSingletSDs,
                                                             const MolDS_base::ElectronicStructure& lhsElectronicStructure) const = 0;

};

}
#endif



