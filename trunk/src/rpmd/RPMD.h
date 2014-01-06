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
#ifndef INCLUDED_RPMD
#define INCLUDED_RPMD
namespace MolDS_rpmd{

/***
 *  Velocty Verlet is used here.
 */
class RPMD : public MolDS_base::PrintController{
public:
   RPMD();
   ~RPMD();
   void DoRPMD(const MolDS_base::Molecule& refferenceMolecule);
private:
   std::string messageStartInitialRPMD;
   std::string messageEndInitialRPMD;
   std::string messageinitialConditionRPMD;
   std::string messageStartRPMD;
   std::string messageEndRPMD;
   std::string messageStartStepRPMD;
   std::string messageEndStepRPMD;
   std::string messageBeadsNum;
   std::string messageEnergies;
   std::string messageEnergiesTitle;
   std::string messageBeadsKineticEnergy;
   std::string messageBeadsHarmonicEnergy;
   std::string messageElecStateEnergy;
   std::string messageTotalEnergy;
   std::string messageErrorEnergy;
   std::string messageTime;
   std::string errorMessageNotEnebleTheoryType;
   std::string errorMessageTheoryType;
   std::string errorMessageElecState;
   std::vector<MolDS_base::TheoryType> enableGroundStateTheoryTypes;
   std::vector<MolDS_base::TheoryType> enableExcitedStateTheoryTypes;
   void SetMessages();
   void SetEnableTheoryTypes();
   void CheckEnableTheoryType(MolDS_base::TheoryType theoryType, int elecState);
   void CreateBeads(std::vector<boost::shared_ptr<MolDS_base::Molecule> >& molecularBeads, 
                    std::vector<boost::shared_ptr<MolDS_base::ElectronicStructure> >& electronicStructureBeads,
                    const MolDS_base::Molecule& refferenceMolecule,
                    int numBeads);
   void UpdateMomenta(const std::vector<boost::shared_ptr<MolDS_base::Molecule> >& molecularBeads, 
                      const std::vector<boost::shared_ptr<MolDS_base::ElectronicStructure> >& electronicStructureBeads,
                      int elecState,
                      double dt,
                      double templerature);
   void UpdateCoordinates(const std::vector<boost::shared_ptr<MolDS_base::Molecule> >& molecularBeads,
                          double dt);
   void UpdateElectronicStructure(const std::vector<boost::shared_ptr<MolDS_base::ElectronicStructure> >& electronicStructureBeads);
   void BroadcastPhaseSpacepointsToAllProcesses(std::vector<boost::shared_ptr<MolDS_base::Molecule> >& molecularBeads, int root) const;
   void FluctuateBeads(const std::vector<boost::shared_ptr<MolDS_base::Molecule> >& molecularBeads,
                       int elecState,
                       double temperature,
                       unsigned long seed);
   //void OutputEnergies(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure, double initialEnergy);
   //double OutputEnergies(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure);
   double OutputEnergies(const std::vector<boost::shared_ptr<MolDS_base::Molecule> >& molecularBeads, 
                         const std::vector<boost::shared_ptr<MolDS_base::ElectronicStructure> >& electronicStructureBeads,
                         int elecState,
                         double temperature);
   void OutputEnergies(const std::vector<boost::shared_ptr<MolDS_base::Molecule> >& molecularBeads, 
                       const std::vector<boost::shared_ptr<MolDS_base::ElectronicStructure> >& electronicStructureBeads,
                       int elecState,
                       double temperature,
                       double initialEnergy);
};

}
#endif



