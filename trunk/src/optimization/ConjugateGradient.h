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
#ifndef INCLUDED_CONJUGATE_GRADIENT
#define INCLUDED_CONJUGATE_GRADIENT
namespace MolDS_optimization{

class ConjugateGradient : public MolDS_optimization::Optimizer{
private:
   class ConjugateGradientState: public OptimizerState{
   private:
      double** oldMatrixForce;
      double** matrixSearchDirection;
      size_t numAtoms;
      ConjugateGradientState(const ConjugateGradientState&); // delete default copy constructor
   public:
      ConjugateGradientState(MolDS_base::Molecule& molecule,
                             const boost::shared_ptr<MolDS_base::ElectronicStructure>& electronicStructure,
                             const boost::shared_ptr<MolDS_base_constraints::Constraint>& constraint);
      virtual ~ConjugateGradientState();
      double** GetOldMatrixForce(){return this->oldMatrixForce;}
      double** GetMatrixSearchDirection(){return this->matrixSearchDirection;}
   };
public:
   ConjugateGradient();
   ~ConjugateGradient();
protected:
   void SetMessages();
private:
   std::string messageStartConjugateGradientStep;
   const std::string& OptimizationStepMessage() const{
      return this->messageStartConjugateGradientStep;
   }
   OptimizerState* CreateState(MolDS_base::Molecule& molecule,
                               const boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                               const boost::shared_ptr<MolDS_base_constraints::Constraint> constraint) const{
      return new ConjugateGradientState(molecule, electronicStructure, constraint);
   }
   void InitializeState(OptimizerState &state, const MolDS_base::Molecule& molecule) const;
   virtual void PrepareState(OptimizerState& state,
                             const MolDS_base::Molecule& molecule,
                             const boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                             const int elecState) const;
   void CalcNextStepGeometry(MolDS_base::Molecule &molecule,
                             OptimizerState& state,
                             boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                             const int elecState,
                             const double dt) const;
   void UpdateState(OptimizerState& state) const;
   void UpdateSearchDirection(OptimizerState& state,
                              boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                              const MolDS_base::Molecule& molecule,
                              boost::shared_ptr<MolDS_base_constraints::Constraint> constraint,
                              int elecState) const;
};

}
#endif



