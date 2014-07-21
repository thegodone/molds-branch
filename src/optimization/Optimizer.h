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
protected:
   class OptimizerState{
   protected:
      MolDS_base::Molecule& molecule;
      const boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure;
      const boost::shared_ptr<MolDS_base_constraints::Constraint> constraint;
      const int elecState;
      const double dt;
      const int totalSteps;
      const double maxGradientThreshold;
      const double rmsGradientThreshold;
      double currentEnergy;
      double initialEnergy;
      double const* const* matrixForce;
      std::string errorMessageFailedToDowncastState;
      virtual void SetMessages();
   private:
      OptimizerState(const OptimizerState&); // delete default copy constructor
   public:
      OptimizerState(MolDS_base::Molecule& molecule,
                     const boost::shared_ptr<MolDS_base::ElectronicStructure>& electronicStructure,
                     const boost::shared_ptr<MolDS_base_constraints::Constraint>& constraint);
      virtual ~OptimizerState(){}
      double& GetCurrentEnergyRef(){return this->currentEnergy;}
      double GetCurrentEnergy()const{return this->currentEnergy;}
      MolDS_base::Molecule& GetMolecule(){return this->molecule;}
      const boost::shared_ptr<MolDS_base_constraints::Constraint> GetConstraint(){return this->constraint;}
      const boost::shared_ptr<MolDS_base::ElectronicStructure> GetElectronicStructure(){
         return this->electronicStructure;
      }
      int GetElecState()const{return this->elecState;}
      double GetDeltaT()const{return this->dt;}
      int GetTotalSteps()const{return this->totalSteps;}
      double GetMaxGradientThreshold()const{return this->maxGradientThreshold;}
      double GetRmsGradientThreshold()const{return this->rmsGradientThreshold;}
      double GetInitialEnergy()const{return this->initialEnergy;}
      double const* const*  GetMatrixForce()const{return this->matrixForce;}
      double const* const** GetMatrixForcePtr(){return &this->matrixForce;}
      void SetCurrentEnergy(double currentEnergy){this->currentEnergy = currentEnergy;}
      void SetInitialEnergy(double initialEnergy){this->initialEnergy = initialEnergy;}
      void SetMatrixForce(double const* const* matrixForce){this->matrixForce = matrixForce;}
      template<class State>
      State& CastRef() throw(MolDS_base::MolDSException){
         try{
            return dynamic_cast<State&>(*this);
         }
         catch(std::bad_cast& ex){
            throw MolDS_base::MolDSException(this->errorMessageFailedToDowncastState);
         }
      }
   };
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
   void OutputOptimizationStepMessage(int nthStep) const;
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
   void SearchMinimum(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                      MolDS_base::Molecule& molecule,
                      boost::shared_ptr<MolDS_base_constraints::Constraint> constraint,
                      double* lineSearchedEnergy,
                      bool* obainesOptimizedStructure) const;
   virtual OptimizerState* CreateState(MolDS_base::Molecule& molecule,
                                       const boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                                       const boost::shared_ptr<MolDS_base_constraints::Constraint> constraint) const{
      return new OptimizerState(molecule, electronicStructure, constraint);
   }
   virtual void InitializeState(OptimizerState &state, const MolDS_base::Molecule& molecule) const = 0;
   virtual void PrepareState(OptimizerState& state,
                             const MolDS_base::Molecule& molecule,
                             const boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                             const int elecState) const = 0;
   virtual void CalcNextStepGeometry(MolDS_base::Molecule& molecule,
                                     OptimizerState& state,
                                     boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                                     const int elecState,
                                     const double dt) const= 0;
   virtual void UpdateState(OptimizerState& state) const = 0;
   virtual const std::string& OptimizationStepMessage() const = 0;
};

}
#endif



