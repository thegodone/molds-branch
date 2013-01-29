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
#ifndef INCLUDED_MOLECULE
#define INCLUDED_MOLECULE
namespace MolDS_base{

class Molecule : public MolDS_base::PrintController{
public:
   Molecule();
   explicit Molecule(const Molecule& rhs);
   Molecule& operator=(const Molecule& rhs);
   ~Molecule();
   inline int GetNumberAtoms() const{
#ifdef MOLDS_DBG
      if(this->atomVect==NULL) throw MolDS_base::MolDSException(this->errorMessageGetNumberAtomsNull);
#endif
      return this->atomVect->size();
   }
   inline MolDS_base_atoms::Atom* GetAtom(int atomIndex) const{
#ifdef MOLDS_DBG
      if(this->atomVect==NULL) throw MolDS_base::MolDSException(this->errorMessageGetAtomNull);
#endif
      return (*this->atomVect)[atomIndex];
   }
   void AddAtom(MolDS_base_atoms::Atom* atom);
   double* GetXyzCOM() const;
   double* GetXyzCOM();
   double* GetXyzCOC() const;
   double* GetXyzCOC();
   void CalcXyzCOM();
   void CalcXyzCOC();
   void CalcBasics();
   int GetTotalNumberAOs() const;
   int GetTotalNumberValenceElectrons() const;
   double GetTotalCoreMass() const;
   void OutputXyzCOM() const;
   void OutputXyzCOC() const;
   void OutputTotalNumberAtomsAOsValenceelectrons() const;
   void OutputConfiguration() const;
   void OutputMomenta() const;
   void CalcPrincipalAxes();
   void Rotate();
   void Translate();
   double GetDistanceAtoms(int indexAtomA, int indexAtomB) const;
   double GetDistanceAtoms(const MolDS_base_atoms::Atom& atomA, 
                           const MolDS_base_atoms::Atom& atomB) const;
   void SynchronizeConfigurationTo  (const Molecule& ref);
   void SynchronizeMomentaTo        (const Molecule& ref);
   void SynchronizePhaseSpacePointTo(const Molecule& ref);
private:
   std::vector<MolDS_base_atoms::Atom*>* atomVect;
   double* xyzCOM; // x, y, z coordinates of Center of Mass;
   double* xyzCOC; // x, y, z coordinates of Center of Core;
   bool wasCalculatedXyzCOM;
   bool wasCalculatedXyzCOC;
   int totalNumberAOs;
   int totalNumberValenceElectrons;
   double totalCoreMass;
   void Initialize();
   void CopyInitialize(const Molecule& rhs);
   void Finalize(std::vector<MolDS_base_atoms::Atom*>** atomVect, double** xyzCOM, double**xyzCOC);
   void SetMessages();
   void CalcTotalNumberValenceElectrons();
   void CalcTotalNumberAOs();
   void CalcTotalCoreMass();
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
   std::string errorMessageGetAtomNull;
   std::string errorMessageAddAtomNull;
   std::string errorMessageGetNumberAtomsNull;
   std::string errorMessageGetXyzCOCNull;
   std::string errorMessageGetXyzCOMNull;
   std::string messageTotalNumberAOs;
   std::string messageTotalNumberAtoms;
   std::string messageTotalNumberValenceElectrons;
   std::string messageAtomCoordinates;
   std::string messageAtomCoordinatesTitle;
   std::string messageAtomMomenta;
   std::string messageAtomMomentaTitle;
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





