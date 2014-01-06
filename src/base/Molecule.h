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
#ifndef INCLUDED_MOLECULE
#define INCLUDED_MOLECULE
namespace MolDS_base{

class Molecule : public MolDS_base::PrintController{
public:
   Molecule();
   explicit Molecule(const Molecule& rhs);
   Molecule& operator=(const Molecule& rhs);
   ~Molecule();
   inline const std::vector<MolDS_base_atoms::Atom*>& GetAtomVect() const{
#ifdef MOLDS_DBG
      if(this->atomVect==NULL) throw MolDS_base::MolDSException(this->errorMessageGetAtomVectNull);
#endif
      return *this->atomVect;
   }
   inline const std::vector<MolDS_base_atoms::Atom*>& GetRealAtomVect() const{
#ifdef MOLDS_DBG
      if(this->realAtomVect==NULL) throw MolDS_base::MolDSException(this->errorMessageGetRealAtomVectNull);
#endif
      return *this->realAtomVect;
   }
   inline const std::vector<MolDS_base_atoms::Atom*>& GetGhostAtomVect() const{
#ifdef MOLDS_DBG
      if(this->ghostAtomVect==NULL) throw MolDS_base::MolDSException(this->errorMessageGetGhostAtomVectNull);
#endif
      return *this->ghostAtomVect;
   }
   inline const std::vector<MolDS_base_atoms::Atom*>& GetEpcVect() const{
#ifdef MOLDS_DBG
      if(this->epcVect==NULL) throw MolDS_base::MolDSException(this->errorMessageGetEPCVectNull);
#endif
      return *this->epcVect;
   }
   void AddAtom(MolDS_base_atoms::Atom* atom);
   void AddRealAtom(MolDS_base_atoms::Atom* atom);
   void AddGhostAtom(MolDS_base_atoms::Atom* atom);
   void AddEpc(MolDS_base_atoms::Atom* epc);
   double const* GetXyzCOM() const; // Get the Cartesian coordinates of the center of atom's mass
   double const* GetXyzCOC() const; // Get the Cartesian coordinates of the cneter of core's mass
   double const* GetXyzDipoleCenter() const{return this->GetXyzCOC();}
   void CalcBasics();
   void CalcBasicsConfiguration();
   int GetTotalNumberAOs() const{return this->totalNumberAOs;}
   inline int GetTotalNumberValenceElectrons() const{return this->totalNumberValenceElectrons;}
   double GetTotalCoreMass() const{return this->totalCoreMass;};
   void OutputXyzCOM() const;
   void OutputXyzCOC() const;
   void OutputTotalNumberAtomsAOsValenceelectrons() const;
   void OutputConfiguration() const;
   void OutputMomenta() const;
   void OutputEpcs() const;
   void CalcPrincipalAxes();
   void Rotate();
   void Translate();
   inline double GetDistanceAtoms(int indexAtomA, int indexAtomB) const{return this->distanceAtoms[indexAtomA][indexAtomB];};
   inline double GetDistanceAtoms(const MolDS_base_atoms::Atom& atomA, 
                                  const MolDS_base_atoms::Atom& atomB) const{return this->GetDistanceAtoms(atomA.GetIndex(), atomB.GetIndex());};
   inline double GetDistanceEpcs(int indexEpcA, int indexEpcB) const{return this->distanceEpcs[indexEpcA][indexEpcB];};
   inline double GetDistanceEpcs(const MolDS_base_atoms::Atom& epcA, 
                                 const MolDS_base_atoms::Atom& epcB) const{return this->GetDistanceEpcs(epcA.GetIndex(), epcB.GetIndex());};
   double GetDistanceAtomEpc(int indexAtom, int indexEpc) const{return this->distanceAtomsEpcs[indexAtom][indexEpc];};
   double GetDistanceAtomEpc(const MolDS_base_atoms::Atom& atom, 
                             const MolDS_base_atoms::Atom& epc) const{return this->GetDistanceAtomEpc(atom.GetIndex(), epc.GetIndex());};
   void SynchronizeConfigurationTo  (const Molecule& ref);
   void SynchronizeMomentaTo        (const Molecule& ref);
   void SynchronizePhaseSpacePointTo(const Molecule& ref);
   void BroadcastConfigurationToAllProcesses(int root) const;
   void BroadcastMomentaToAllProcesses(int root) const;
   void BroadcastPhaseSpacePointToAllProcesses(int root) const;
private:
   std::vector<MolDS_base_atoms::Atom*>* atomVect;
   std::vector<MolDS_base_atoms::Atom*>* realAtomVect; // Vector of real (=not ghost) atoms
   std::vector<MolDS_base_atoms::Atom*>* ghostAtomVect;   // Vector of ghost atoms
   std::vector<MolDS_base_atoms::Atom*>* epcVect;      // Vector of Environmental Point Charges
   double*  xyzCOM; // x, y, z coordinates of the center of atomic mass;
   double*  xyzCOC; // x, y, z coordinates of the center of core's mass;
   double** distanceAtoms;    // distance between each atom;
   double** distanceEpcs;     // distance between each environmental point charge;
   double** distanceAtomsEpcs;// distance between each atom and environmental point charge;
   int totalNumberAOs;
   int totalNumberValenceElectrons;
   double totalCoreMass;
   void Initialize();
   void CopyInitialize(const Molecule& rhs);
   void Finalize(std::vector<MolDS_base_atoms::Atom*>** atomVect, 
                 std::vector<MolDS_base_atoms::Atom*>** realAtomVect,
                 std::vector<MolDS_base_atoms::Atom*>** ghostAtomVect,
                 std::vector<MolDS_base_atoms::Atom*>** epcVect,
                 double** xyzCOM, 
                 double** xyzCOC, 
                 double*** distanceAtoms,
                 double*** distanceEpcs,
                 double*** distanceAtomsEpcs);
   void SetMessages();
   void CopyRealGhostAtom2Atom();
   void CalcTotalNumberAOs();
   void CalcTotalNumberValenceElectrons();
   void CalcTotalCoreMass();
   void CalcXyzCOM();
   void CalcXyzCOC();
   void CalcDistanceAtoms();
   void CalcDistanceEpcs();
   void CalcDistanceAtomsEpcs();
   void CalcInertiaTensor(double** inertiaTensor, 
                          double const* inertiaTensorOrigin);
   void FreeInertiaTensorMoments(double*** inertiaTensor, 
                                 double** inertiaMoments);
   void Rotate(MolDS_base::EularAngle eularAngle, 
               const double* rotatingOrigin, 
               RotatedObjectType rotatedObj);
   void OutputPrincipalAxes(double const* const* inertiaTensor, 
                            double const* inertiaMoments) const;
   void OutputInertiaTensorOrigin(double* inertiaTensorOrigin) const;
   void OutputRotatingConditions(RotatingType rotatingType, 
                                 double const* rotatingOrigin, 
                                 double const* rotatingAxis, 
                                 double rotatingAngle, 
                                 MolDS_base::EularAngle rotatingEularAngles)const;
   void OutputTranslatingConditions(double const* translatingDifference) const;
   std::string errorMessageGetAtomVectNull;
   std::string errorMessageGetRealAtomVectNull;
   std::string errorMessageGetGhostAtomVectNull;
   std::string errorMessageGetEPCVectNull;
   std::string errorMessageAddAtomNull;
   std::string errorMessageAddRealAtomNull;
   std::string errorMessageAddGhostAtomNull;
   std::string errorMessageAddEPCNull;
   std::string errorMessageCopyRealGhostAtom2AtomNotEmpty;
   std::string errorMessageGetXyzCOMNull;
   std::string errorMessageGetXyzCOCNull;
   std::string errorMessageCalcXyzCOMNull;
   std::string errorMessageCalcXyzCOCNull;
   std::string messageTotalNumberAOs;
   std::string messageTotalNumberAtoms;
   std::string messageTotalNumberValenceElectrons;
   std::string messageAtomCoordinates;
   std::string messageAtomCoordinatesTitle;
   std::string messageAtomMomenta;
   std::string messageAtomMomentaTitle;
   std::string messageEpcConfiguration;
   std::string messageEpcCoordinates;
   std::string messageEpcCoordinatesTitle;
   std::string messageCOM;
   std::string messageCOC;
   std::string messageCOMTitle;
   std::string messageStartPrincipalAxes;
   std::string messageDonePrincipalAxes;
   std::string messagePrincipalAxes;
   std::string messagePrincipalAxesNote;
   std::string messagePrincipalAxesTitle;
   std::string messageInertiaTensorOrigin;
   std::string messageInertiaTensorOriginTitle;
   std::string messageStartRotate;
   std::string messageDoneRotate;
   std::string messageRotatingOrigin;
   std::string messageRotatingOriginTitle;
   std::string messageRotatingAxis;
   std::string messageRotatingAxisTitle;
   std::string messageRotatingAngle;
   std::string messageRotatingType;
   std::string messageRotatingEularAngles;
   std::string messageRotatingEularAnglesTitle;
   std::string messageStartTranslate;
   std::string messageDoneTranslate;
   std::string messageTranslatingDifference;
   std::string messageTranslatingDifferenceTitle;
};
}
#endif





