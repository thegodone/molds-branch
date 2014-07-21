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
#ifndef INCLUDED_GEDIIS
#define INCLUDED_GEDIIS

#include <list>

namespace MolDS_optimization{

class GEDIIS : public MolDS_optimization::BFGS{
private:
   class GEDIISHistory{
   public:
      GEDIISHistory();
      ~GEDIISHistory();
      void AddEntry(double energy,
                    const MolDS_base::Molecule& molecule,
                    double const* const* matrixForce);
      void DiscardEntries();
      void SolveGEDIISEquation(double* gediisEnergy, double** matrixCoordinate, double** matrixForce);
   private:
      class Entry{
      public:
         Entry(double energy,
               const MolDS_base::Molecule& molecule,
               double const* const* matrixForce);
         ~Entry();
         int GetNumberAtoms()                 const { return this->numAtoms;         }
         double GetEnergy()                   const { return this->energy;           }
         double const* const* GetCoordinate() const { return this->matrixCoordinate; }
         double const* const* GetForce()      const { return this->matrixForce;      }
      private:
         int      numAtoms;
         double   energy;
         double** matrixCoordinate;
         double** matrixForce;
      };
      const int maxEntryCount;
      typedef std::list< const Entry* > entryList_t;
      entryList_t entryList;
      void SetMessages();

      std::string errorMessageNegativeGEDIISCoefficient;
      std::string errorMessageNotSufficientHistory;
   };

   class GEDIISState: public BFGSState{
   protected:
      GEDIISHistory history;
      double   preRFOEnergy;
      double** matrixGEDIISCoordinates;
      double** matrixGEDIISForce;
      bool     isGEDIISStepTaken;
   public:
      GEDIISState(MolDS_base::Molecule& molecule,
                  const boost::shared_ptr<MolDS_base::ElectronicStructure>& electronicStructure,
                  const boost::shared_ptr<MolDS_base_constraints::Constraint>& constraint);
      virtual ~GEDIISState();
      GEDIISHistory& GetHistory(){return this->history;}
      double** GetMatrixGEDIISCoordinates(){return this->matrixGEDIISCoordinates;}
      double** GetMatrixGEDIISForce()      {return this->matrixGEDIISForce;}
      bool     GetIsGEDIISStepTaken()const {return this->isGEDIISStepTaken;}
      void SetIsGEDIISStepTaken(bool isGEDIISStepTaken){
         this->isGEDIISStepTaken = isGEDIISStepTaken;
      }
      void SetPreRFOEnergy(double preRFOEnergy){this->preRFOEnergy = preRFOEnergy;}
      double GetPreRFOEnergy() const{
         return GetIsGEDIISStepTaken() ? this->preRFOEnergy : this->GetInitialEnergy();
      }
      double const* GetVectorForceForRFO()const{
         return GetIsGEDIISStepTaken() ? &this->matrixGEDIISForce[0][0] : this->GetVectorForce();
      }
   };
public:
   GEDIIS();
   ~GEDIIS();
protected:
   void SetMessages();

   std::string messageStartGEDIISStep;
   std::string messageTakingGEDIISStep;
   std::string messageTakingRFOStep;
   std::string messageDiscardHistory;

private:
   const std::string& OptimizationStepMessage() const{
      return this->messageStartGEDIISStep;
   }
   virtual OptimizerState* CreateState(MolDS_base::Molecule& molecule,
                                       const boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                                       const boost::shared_ptr<MolDS_base_constraints::Constraint> constraint) const{
      return new GEDIISState(molecule, electronicStructure, constraint);
   }
protected:
   void InitializeState(OptimizerState &state, const MolDS_base::Molecule& molecule) const;
   void CalcNextStepGeometry(MolDS_base::Molecule &molecule,
                             OptimizerState& state,
                             boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                             const int elecState,
                             const double dt) const;
   void UpdateState(OptimizerState& state) const;
};

}
#endif
