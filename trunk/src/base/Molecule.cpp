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
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<vector>
#include<stdexcept>
#include<boost/format.hpp>
#include"PrintController.h"
#include"MolDSException.h"
#include"Uncopyable.h"
#include"../wrappers/Lapack.h"
#include"Enums.h"
#include"MathUtilities.h"
#include"MallocerFreer.h"
#include"EularAngle.h"
#include"Parameters.h"
#include"atoms/Atom.h"
#include"factories/AtomFactory.h"
#include"Molecule.h"
using namespace std;
using namespace MolDS_base_atoms;
using namespace MolDS_base_factories;
namespace MolDS_base{

Molecule::Molecule(){
   this->Initialize();
}

Molecule::Molecule(const Molecule& rhs){
   this->CopyInitialize(rhs);
}

Molecule& Molecule::operator=(const Molecule& rhs){
   double* oldXyzCOM = this->xyzCOM;
   double* oldXyzCOC = this->xyzCOC;
   vector<Atom*>* oldAtomVect = this->atomVect;
   this->CopyInitialize(rhs);
   this->Finalize(&oldAtomVect, &oldXyzCOM, &oldXyzCOC);
   return *this;
}

Molecule::~Molecule(){
   this->Finalize(&this->atomVect, &this->xyzCOM, &this->xyzCOC);
   //this->OutputLog("molecule deleted\n");
}

void Molecule::CopyInitialize(const Molecule& rhs){
   this->Initialize();
   for(int i=0; i<CartesianType_end; i++){
      this->xyzCOM[i] = rhs.xyzCOM[i];
      this->xyzCOC[i] = rhs.xyzCOC[i];
   }
   this->wasCalculatedXyzCOM = rhs.wasCalculatedXyzCOM;
   this->wasCalculatedXyzCOC = rhs.wasCalculatedXyzCOC;
   this->totalNumberAOs = rhs.totalNumberAOs;
   this->totalNumberValenceElectrons = rhs.totalNumberValenceElectrons;
   this->totalCoreMass = rhs.totalCoreMass;
   for(int i=0; i<rhs.atomVect->size(); i++){
      Atom* atom = (*rhs.atomVect)[i];
      this->atomVect->push_back(AtomFactory::Create(atom->GetAtomType(),
                                                    atom->GetXyz()[XAxis],
                                                    atom->GetXyz()[YAxis],
                                                    atom->GetXyz()[ZAxis],
                                                    atom->GetPxyz()[XAxis],
                                                    atom->GetPxyz()[YAxis],
                                                    atom->GetPxyz()[ZAxis]));
      (*this->atomVect)[i]->SetFirstAOIndex(atom->GetFirstAOIndex());
   }                                                                     
}

void Molecule::Initialize(){
   this->SetMessages();
   this->wasCalculatedXyzCOM = false;
   this->wasCalculatedXyzCOC = false;
   this->xyzCOM = NULL;
   this->xyzCOC = NULL;
   try{
      this->atomVect = new vector<Atom*>;
      MallocerFreer::GetInstance()->Malloc<double>(&this->xyzCOM, CartesianType_end);
      MallocerFreer::GetInstance()->Malloc<double>(&this->xyzCOC, CartesianType_end);
   }
   catch(exception ex){
      this->Finalize(&this->atomVect, &this->xyzCOM, &this->xyzCOC);
      throw MolDSException(ex.what());
   }
}

void Molecule::Finalize(vector<Atom*>** atomVect, double** xyzCOM, double**xyzCOC){
   if(*atomVect != NULL){
      for(int i=0; i<(*atomVect)->size(); i++){
         if((**atomVect)[i] != NULL){
            delete (**atomVect)[i];
            (**atomVect)[i] = NULL;
         }
      }
      (*atomVect)->clear();
      delete *atomVect;
      *atomVect = NULL;
      //this->OutputLog("atomVect deleted\n");
   }
   MallocerFreer::GetInstance()->Free<double>(xyzCOM, CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(xyzCOC, CartesianType_end);
}

void Molecule::SetMessages(){
   this->errorMessageGetAtomVectNull = "Error in base::Molecule::GetAtomVect: atomVect is NULL.\n";
   this->errorMessageGetXyzCOCNull = "Error in base::Molecule::GetXyzCOC: xyzCOC is NULL.\n";
   this->errorMessageGetXyzCOMNull = "Error in base::Molecule::GetXyzCOM: xyzCOM is NULL.\n";
   this->messageTotalNumberAOs = "\tTotal number of valence AOs: ";
   this->messageTotalNumberAtoms = "\tTotal number of atoms: ";
   this->messageTotalNumberValenceElectrons = "\tTotal number of valence electrons: ";
   this->messageAtomCoordinatesTitle = "\t\t\t\t| i-th | atom type |   x[a.u.]   |   y[a.u.]   |   z[a.u.]   |\t\t|  x[angst.]  |  y[angst.]  |  z[angst.]  |\n";
   this->messageAtomCoordinates = "\tAtom coordinates:";
   this->messageAtomMomenta = "\tAtom momenta:";
   this->messageAtomMomentaTitle = "\t\t\t| i-th | atom type |   px[a.u.]   |   py[a.u.]   |   pz[a.u.]   |\t\t|   px[u]   |   py[u]   |   pz[u]   | [u] = [(g/Mol)*(angst/fs)]\n";
   this->messageCOM = "\tCenter of Mass:";
   this->messageCOC = "\tCenter of Core:";
   this->messageCOMTitle = "\t\t\t|   x[a.u.]   |   y[a.u.]   |   z[a.u.]   |\t\t|  x[angst.]  |  y[angst.]  |  z[angst.]  |\n";
   this->messageStartPrincipalAxes = "**********  START: Principal Axes of Inertia  **********\n";
   this->messageDonePrincipalAxes =  "**********  DONE: Principal Axes of Inertia  ***********\n\n\n";
   this->messagePrincipalAxes = "\tPrincipal Axis:";
   this->messagePrincipalAxesNote = "\tThe principal Axes in [a.u.] is normalized while the one in [angst.] is not normalized.\n";
   this->messagePrincipalAxesTitle = "\t\t\t| inertia moments [a.u.] | x[a.u.] | y[a.u.] | z[a.u.] |\t\t | inertia moments [g*angust**2/mol] | x[angst.] | y[angst.] | z[angst.] |\n";
   this->messageInertiaTensorOrigin = "\tInertia Tensor Origin:";
   this->messageInertiaTensorOriginTitle = "\t\t\t\t|   x[a.u.]   |   y[a.u.]   |   z[a.u.]   |\t\t| x[angst.] | y[angst.] | z[angst.] |\n";
   this->messageStartRotate = "**********  START: Rotate molecule  **********\n";
   this->messageDoneRotate =  "**********  DONE: Rotate molecule  ***********\n\n\n";
   this->messageRotatingOrigin = "\tRotating Origin:";
   this->messageRotatingOriginTitle = "\t\t\t\t|   x[a.u.]   |   y[a.u.]   |   z[a.u.]   |\t\t| x[angst.] | y[angst.] | z[angst.] |\n";
   this->messageRotatingAxis = "\tRotating Axis:";
   this->messageRotatingAxisTitle = "\t\t\t|   x[a.u.]   |   y[a.u.] |   z[a.u.]   |\t\t| x[angst.] | y[angst.] | z[angst.] |\n";
   this->messageRotatingAngle = "\tRotating Angle [degree]: ";
   this->messageRotatingType = "\tRotating Type: ";
   this->messageRotatingEularAngles = "\tRotating Eular Angles:";
   this->messageRotatingEularAnglesTitle = "\t\t\t\t| alpha[degree] | beta[degree] | gamma[degree] |\n";
   this->messageStartTranslate = "**********  START: Translate molecule  **********\n";
   this->messageDoneTranslate =  "**********  DONE: Translate molecule  ***********\n\n\n";
   this->messageTranslatingDifference = "\tTranslating Difference:";
   this->messageTranslatingDifferenceTitle = "\t\t\t\t|   x[a.u.]   |   y[a.u.]   |   z[a.u.]   |\t\t| x[angst.] | y[angst.] | z[angst.] |\n";
}

void Molecule::AddAtom(Atom* atom){
   this->atomVect->push_back(atom);
}

double* Molecule::GetXyzCOM() const{
   if(this->xyzCOM==NULL){
      stringstream ss;
      ss << this->errorMessageGetXyzCOMNull;
      throw MolDSException(ss.str());
   }
   return this->xyzCOM;
}

double* Molecule::GetXyzCOM(){
   if(!this->wasCalculatedXyzCOM){
      this->CalcXyzCOM();
   }
   return this->xyzCOM;
}

double* Molecule::GetXyzCOC() const{
   if(this->xyzCOC==NULL){
      stringstream ss;
      ss << this->errorMessageGetXyzCOCNull;
      throw MolDSException(ss.str());
   }
   return this->xyzCOC;
}

double* Molecule::GetXyzCOC(){
   if(!this->wasCalculatedXyzCOC){
      this->CalcXyzCOC();
   }
   return this->xyzCOC;
}

void Molecule::CalcXyzCOM(){
   MallocerFreer::GetInstance()->Malloc<double>(&this->xyzCOM, CartesianType_end);
   double totalAtomicMass = 0.0;
   double* atomicXyz;
   double atomicMass = 0.0;

   for(int j=0; j<3; j++){
      this->xyzCOM[j] = 0.0;
   }
      
   for(int i=0; i<this->atomVect->size(); i++){
      const Atom& atom = *(*this->atomVect)[i]; 
      atomicXyz = atom.GetXyz();
      atomicMass = atom.GetAtomicMass();
      totalAtomicMass += atomicMass;
      for(int j=0; j<3; j++){
         this->xyzCOM[j] += atomicXyz[j] * atomicMass;
      }
   }
   for(int i=0; i<3; i++){
      this->xyzCOM[i]/=totalAtomicMass;
   }
   this->wasCalculatedXyzCOM = true;
}

void Molecule::CalcXyzCOC(){
   MallocerFreer::GetInstance()->Malloc<double>(&this->xyzCOC, CartesianType_end);
   double totalCoreMass = 0.0;
   double* atomicXyz;
   double coreMass = 0.0;

   for(int j=0; j<3; j++){
      this->xyzCOC[j] = 0.0;
   }
      
   for(int i=0; i<this->atomVect->size(); i++){
      const Atom& atom = *(*this->atomVect)[i]; 
      atomicXyz = atom.GetXyz();
      coreMass = atom.GetCoreMass();
      totalCoreMass += coreMass;
      for(int j=0; j<3; j++){
         this->xyzCOC[j] += atomicXyz[j] * coreMass;
      }
   }
   for(int i=0; i<3; i++){
      this->xyzCOC[i]/=totalCoreMass;
   }
   this->wasCalculatedXyzCOC = true;
}

int Molecule::GetTotalNumberAOs() const{
   return this->totalNumberAOs;
}

void Molecule::CalcBasics(){
   this->CalcTotalNumberAOs();
   this->CalcTotalNumberValenceElectrons();
   this->CalcTotalCoreMass();
   this->CalcXyzCOM();
   this->CalcXyzCOC();
}

void Molecule::CalcTotalNumberAOs(){
   this->totalNumberAOs = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      (*this->atomVect)[i]->SetFirstAOIndex(totalNumberAOs);
      this->totalNumberAOs += (*this->atomVect)[i]->GetValenceSize();
   }
}

int Molecule::GetTotalNumberValenceElectrons() const{
   return this->totalNumberValenceElectrons;
}

void Molecule::CalcTotalNumberValenceElectrons(){
   this->totalNumberValenceElectrons = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      this->totalNumberValenceElectrons += (*this->atomVect)[i]->GetNumberValenceElectrons();
   }
}

double Molecule::GetTotalCoreMass() const{
   return this->totalCoreMass;
}

void Molecule::CalcTotalCoreMass(){
   this->totalCoreMass = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      const Atom& atom = *(*this->atomVect)[i]; 
      double coreMass = atom.GetAtomicMass() - static_cast<double>(atom.GetNumberValenceElectrons());
      this->totalCoreMass += coreMass;
   }
}

void Molecule::OutputConfiguration() const{
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   this->OutputLog(this->messageAtomCoordinatesTitle);
   for(int a=0; a<this->atomVect->size(); a++){
      const Atom& atom = *(*this->atomVect)[a];
      this->OutputLog(boost::format("%s\t%d\t%s\t%e\t%e\t%e\t\t%e\t%e\t%e\n") 
         % this->messageAtomCoordinates
         % a
         % AtomTypeStr(atom.GetAtomType()) 
         % atom.GetXyz()[0]
         % atom.GetXyz()[1]
         % atom.GetXyz()[2]
         % (atom.GetXyz()[0]/ang2AU)
         % (atom.GetXyz()[1]/ang2AU)
         % (atom.GetXyz()[2]/ang2AU));
   }
   this->OutputLog("\n");
}

void Molecule::OutputMomenta() const{
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double fs2AU = Parameters::GetInstance()->GetFs2AU();
   double gMolin2AU = Parameters::GetInstance()->GetGMolin2AU();
   double momentumUnit2AU = ang2AU*gMolin2AU/fs2AU;
   this->OutputLog(this->messageAtomMomentaTitle);
   for(int a=0; a<this->atomVect->size(); a++){
      const Atom& atom = *(*this->atomVect)[a];
      this->OutputLog(boost::format("%s\t%d\t%s\t%e\t%e\t%e\t\t%e\t%e\t%e\n") 
         % this->messageAtomMomenta
         % a
         % AtomTypeStr(atom.GetAtomType())
         % (atom.GetPxyz()[0]/momentumUnit2AU)
         % (atom.GetPxyz()[1]/momentumUnit2AU)
         % (atom.GetPxyz()[2]/momentumUnit2AU)
         % atom.GetPxyz()[0]
         % atom.GetPxyz()[1]
         % atom.GetPxyz()[2]);
   }
   this->OutputLog("\n");
}

void Molecule::OutputXyzCOM() const{
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   this->OutputLog(this->messageCOMTitle);
   this->OutputLog(boost::format("%s\t%e\t%e\t%e\t\t%e\t%e\t%e\n") 
      % this->messageCOM
      % this->xyzCOM[0]
      % this->xyzCOM[1]
      % this->xyzCOM[2]
      % (this->xyzCOM[0]/ang2AU)
      % (this->xyzCOM[1]/ang2AU)
      % (this->xyzCOM[2]/ang2AU));
   this->OutputLog("\n");
}

void Molecule::OutputXyzCOC() const{
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   this->OutputLog(this->messageCOMTitle);
   this->OutputLog(boost::format("%s\t%e\t%e\t%e\t\t%e\t%e\t%e\n") 
      % this->messageCOC
      % this->xyzCOC[0]
      % this->xyzCOC[1]
      % this->xyzCOC[2]
      % (this->xyzCOC[0]/ang2AU)
      % (this->xyzCOC[1]/ang2AU)
      % (this->xyzCOC[2]/ang2AU));
   this->OutputLog("\n");
}

void Molecule::OutputTotalNumberAtomsAOsValenceelectrons() const{
   this->OutputLog(boost::format("%s%d\n") % this->messageTotalNumberAtoms
                                           % this->atomVect->size());
   this->OutputLog(boost::format("%s%d\n") % this->messageTotalNumberAOs
                                           % this->totalNumberAOs);
   this->OutputLog(boost::format("%s%d\n\n") % this->messageTotalNumberValenceElectrons
                                             % this->totalNumberValenceElectrons);
}

void Molecule::OutputPrincipalAxes(double const* const* inertiaTensor, 
                                   double const* inertiaMoments) const{
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double gMolin2AU = Parameters::GetInstance()->GetGMolin2AU();
   this->OutputLog(this->messagePrincipalAxesTitle);
   for(int i=0; i<3; i++){
      this->OutputLog(boost::format("%s\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\t%e\n") 
         % this->messagePrincipalAxes
         % inertiaMoments[i]
         % inertiaTensor[i][0]
         % inertiaTensor[i][1]
         % inertiaTensor[i][2]
         % (inertiaMoments[i]/gMolin2AU)
         % (inertiaTensor[i][0]/ang2AU)
         % (inertiaTensor[i][1]/ang2AU)
         % (inertiaTensor[i][2]/ang2AU));
   }
   this->OutputLog(this->messagePrincipalAxesNote);
   this->OutputLog("\n");
}

void Molecule::OutputInertiaTensorOrigin(double* inertiaTensorOrigin) const{
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   this->OutputLog(this->messageInertiaTensorOriginTitle);
   this->OutputLog(boost::format("%s\t%e\t%e\t%e\t\t%e\t%e\t%e\n") 
      % this->messageInertiaTensorOrigin
      % inertiaTensorOrigin[0]
      % inertiaTensorOrigin[1]
      % inertiaTensorOrigin[2]
      % (inertiaTensorOrigin[0]/ang2AU)
      % (inertiaTensorOrigin[1]/ang2AU)
      % (inertiaTensorOrigin[2]/ang2AU));

   this->OutputLog("\n");
}

void Molecule::CalcPrincipalAxes(){
   this->OutputLog(this->messageStartPrincipalAxes);
   if(!this->wasCalculatedXyzCOM){
      this->CalcXyzCOM();
   }
   double inertiaTensorOrigin[3] = {this->xyzCOM[0], this->xyzCOM[1], this->xyzCOM[2]};
   if(Parameters::GetInstance()->GetInertiaTensorOrigin() != NULL){
      inertiaTensorOrigin[0] = Parameters::GetInstance()->GetInertiaTensorOrigin()[0];
      inertiaTensorOrigin[1] = Parameters::GetInstance()->GetInertiaTensorOrigin()[1];
      inertiaTensorOrigin[2] = Parameters::GetInstance()->GetInertiaTensorOrigin()[2];
   }

   double** inertiaTensor = NULL;
   double*  inertiaMoments = NULL;
   try{
      MallocerFreer::GetInstance()->Malloc<double>(&inertiaTensor, CartesianType_end, CartesianType_end);
      MallocerFreer::GetInstance()->Malloc<double>(&inertiaMoments, CartesianType_end);
      this->CalcInertiaTensor(inertiaTensor, inertiaTensorOrigin);
      
      bool calcEigenVectors = true;
      MolDS_wrappers::Lapack::GetInstance()->Dsyevd(inertiaTensor,
                                                    inertiaMoments,
                                                    3,
                                                    calcEigenVectors);
      this->OutputPrincipalAxes(inertiaTensor, inertiaMoments);
      this->OutputInertiaTensorOrigin(inertiaTensorOrigin);
   }
   catch(MolDSException ex){
      this->FreeInertiaTensorMoments(&inertiaTensor, &inertiaMoments);
      throw ex;
   }
   this->FreeInertiaTensorMoments(&inertiaTensor, &inertiaMoments);
   this->OutputLog(this->messageDonePrincipalAxes);
}

void Molecule::CalcInertiaTensor(double** inertiaTensor, const double* inertiaTensorOrigin){

   double x;
   double y;
   double z;
   double atomicMass;
   for(int a=0; a<this->atomVect->size(); a++){
      const Atom& atom = *(*this->atomVect)[a];
      atomicMass = atom.GetAtomicMass();
      x = atom.GetXyz()[0] - inertiaTensorOrigin[0];
      y = atom.GetXyz()[1] - inertiaTensorOrigin[1];
      z = atom.GetXyz()[2] - inertiaTensorOrigin[2];

      inertiaTensor[0][0] += atomicMass*(y*y + z*z);
      inertiaTensor[0][1] -= atomicMass*x*y;
      inertiaTensor[0][2] -= atomicMass*x*z;

      inertiaTensor[1][0] -= atomicMass*y*x;
      inertiaTensor[1][1] += atomicMass*(x*x + z*z);
      inertiaTensor[1][2] -= atomicMass*y*z;

      inertiaTensor[2][0] -= atomicMass*z*x;
      inertiaTensor[2][1] -= atomicMass*z*y;
      inertiaTensor[2][2] += atomicMass*(x*x + y*y);

   }
}

void Molecule::FreeInertiaTensorMoments(double*** inertiaTensor, double** inertiaMoments){
      MallocerFreer::GetInstance()->Free<double>(inertiaTensor, CartesianType_end, CartesianType_end);
      MallocerFreer::GetInstance()->Free<double>(inertiaMoments, CartesianType_end);
}

void Molecule::Rotate(){

   this->OutputLog(this->messageStartRotate);

   // Default values are set if some conditions are not specified.
   if(!this->wasCalculatedXyzCOM){
      this->CalcXyzCOM();
   }
   double rotatingOrigin[3] = {this->xyzCOM[0], this->xyzCOM[1], this->xyzCOM[2]};
   if(Parameters::GetInstance()->GetRotatingOrigin() != NULL){
      rotatingOrigin[0] = Parameters::GetInstance()->GetRotatingOrigin()[0];
      rotatingOrigin[1] = Parameters::GetInstance()->GetRotatingOrigin()[1];
      rotatingOrigin[2] = Parameters::GetInstance()->GetRotatingOrigin()[2];
   }

   RotatingType rotatingType = Parameters::GetInstance()->GetRotatingType();
   double* rotatingAxis = Parameters::GetInstance()->GetRotatingAxis();
   EularAngle rotatingEularAngles = Parameters::GetInstance()->GetRotatingEularAngles();
   double rotatingAngle = Parameters::GetInstance()->GetRotatingAngle();

   this->OutputRotatingConditions(rotatingType, rotatingOrigin, 
                                  rotatingAxis, rotatingAngle, 
                                  rotatingEularAngles);

   // rotate
   if(rotatingType == Axis){
      EularAngle setZAxisEularAngles(rotatingAxis[0], rotatingAxis[1], rotatingAxis[2]);
      EularAngle angleAroundAxis;
      angleAroundAxis.SetAlpha(rotatingAngle);

      this->Rotate(setZAxisEularAngles, rotatingOrigin, Frame);
      this->Rotate(angleAroundAxis, rotatingOrigin, System);
      this->Rotate(setZAxisEularAngles, rotatingOrigin, System);
   }
   else if(rotatingType == Eular){
      this->Rotate(rotatingEularAngles, rotatingOrigin, System);
   }
   
   this->OutputConfiguration();
   this->OutputLog(this->messageDoneRotate);
}

/***
 * rotatedObj == System: Molecule is rotated.
 * rotatedObj == Frame: De Cartesian is rotated.
 */
void Molecule::Rotate(EularAngle eularAngle, const double* rotatingOrigin, RotatedObjectType rotatedObj){

   double rotatingMatrixAlpha[3][3];
   double rotatingMatrixBeta[3][3];
   double rotatingMatrixGamma[3][3];
   double inv = 1.0;
   if(rotatedObj == System){
      inv = -1.0;
   }

   CalcRotatingMatrix(rotatingMatrixAlpha, inv*eularAngle.GetAlpha(), ZAxis);
   CalcRotatingMatrix(rotatingMatrixBeta, inv*eularAngle.GetBeta(), YAxis);
   CalcRotatingMatrix(rotatingMatrixGamma, inv*eularAngle.GetGamma(), ZAxis);

   double temp1[3][3];
   for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
         temp1[i][j] = 0.0;
         for(int k=0; k<3; k++){
            if(rotatedObj == System){
               temp1[i][j] += rotatingMatrixBeta[i][k] * rotatingMatrixGamma[k][j];
            }
            else if(rotatedObj == Frame){
               temp1[i][j] += rotatingMatrixBeta[i][k] * rotatingMatrixAlpha[k][j];
            }
         }
      }
   }

   double temp2[3][3];
   for(int i=0; i<3; i++){
      for(int j=0; j<3; j++){
         temp2[i][j] = 0.0;
         for(int k=0; k<3; k++){
            if(rotatedObj == System){
               temp2[i][j] += rotatingMatrixAlpha[i][k] * temp1[k][j];
            }
            else if(rotatedObj == Frame){
               temp2[i][j] += rotatingMatrixGamma[i][k] * temp1[k][j];
            }
         }
      }
   }

   double rotatedXyz[3];
   for(int i=0; i<this->atomVect->size(); i++){
         const Atom& atom = *(*this->atomVect)[i]; 
         for(int j=0; j<3; j++){
            rotatedXyz[j] = 0.0;
            for(int k=0; k<3; k++){
               rotatedXyz[j] += temp2[j][k] * (atom.GetXyz()[k] - rotatingOrigin[k]);
            }
         }
         for(int j=0; j<3; j++){
            atom.GetXyz()[j] = rotatedXyz[j] + rotatingOrigin[j];
         }
   }
}

void Molecule::OutputRotatingConditions(RotatingType rotatingType, 
                                        double const* rotatingOrigin, 
                                        double const* rotatingAxis, 
                                        double rotatingAngle, 
                                        EularAngle rotatingEularAngles) const{

   double angst2AU = Parameters::GetInstance()->GetAngstrom2AU();
   double degree2Radian = Parameters::GetInstance()->GetDegree2Radian();

   // type
   this->OutputLog(boost::format("%s%s\n\n") % this->messageRotatingType.c_str() 
                                              % RotatingTypeStr(rotatingType));

   // rotating origin
   this->OutputLog(this->messageRotatingOriginTitle);
   this->OutputLog(boost::format("%s\t%e\t%e\t%e\t\t%e\t%e\t%e\n\n") 
      % this->messageRotatingOrigin
      % rotatingOrigin[0]
      % rotatingOrigin[1]
      % rotatingOrigin[2]
      % (rotatingOrigin[0]/angst2AU)
      % (rotatingOrigin[1]/angst2AU)
      % (rotatingOrigin[2]/angst2AU));

   if(rotatingType == Axis){
      // rotating axis
      this->OutputLog(this->messageRotatingAxisTitle);
      this->OutputLog(boost::format("%s\t%e\t%e\t%e\t\t%e\t%e\t%e\n\n") 
         % this->messageRotatingAxis
         % rotatingAxis[0]
         % rotatingAxis[1]
         % rotatingAxis[2]
         % (rotatingAxis[0]/angst2AU)
         % (rotatingAxis[1]/angst2AU)
         % (rotatingAxis[2]/angst2AU));

      // angle
      this->OutputLog(boost::format("%s%e\n\n") % this->messageRotatingAngle.c_str() 
                                                % (rotatingAngle/degree2Radian));
   }
   else if (rotatingType == Eular){
      // Eular angles
      this->OutputLog(this->messageRotatingEularAnglesTitle);
      this->OutputLog(boost::format("%s\t%e\t%e\t%e\n\n") % this->messageRotatingEularAngles
                                                          % (rotatingEularAngles.GetAlpha()/degree2Radian)
                                                          % (rotatingEularAngles.GetBeta()/degree2Radian)
                                                          % (rotatingEularAngles.GetGamma()/degree2Radian));
   }
}


void Molecule::Translate(){
   this->OutputLog(this->messageStartTranslate);
   double x = Parameters::GetInstance()->GetTranslatingDifference()[0];
   double y = Parameters::GetInstance()->GetTranslatingDifference()[1];
   double z = Parameters::GetInstance()->GetTranslatingDifference()[2];

   this->OutputTranslatingConditions(Parameters::GetInstance()->GetTranslatingDifference()); 

   for(int i=0; i<this->atomVect->size(); i++){
         const Atom& atom = *(*this->atomVect)[i]; 
         atom.GetXyz()[0] += x;
         atom.GetXyz()[1] += y;
         atom.GetXyz()[2] += z;
   }
   
   this->wasCalculatedXyzCOM = false;
   this->CalcXyzCOM();
   this->wasCalculatedXyzCOC = false;
   this->CalcXyzCOC();

   this->OutputConfiguration();
   this->OutputXyzCOM();
   this->OutputXyzCOC();
   this->OutputLog(this->messageDoneTranslate);
}

void Molecule::OutputTranslatingConditions(double const* translatingDifference) const{
   double angst2AU = Parameters::GetInstance()->GetAngstrom2AU();
   this->OutputLog(this->messageTranslatingDifferenceTitle);
   this->OutputLog(boost::format("%s\t%e\t%e\t%e\t\t%e\t%e\t%e\n\n") 
      % this->messageTranslatingDifference
      % translatingDifference[0]
      % translatingDifference[1]
      % translatingDifference[2]
      % (translatingDifference[0]/angst2AU)
      % (translatingDifference[1]/angst2AU)
      % (translatingDifference[2]/angst2AU));
}

double Molecule::GetDistanceAtoms(int atomAIndex, int atomBIndex) const{
   const Atom& atomA = *(*this->atomVect)[atomAIndex];
   const Atom& atomB = *(*this->atomVect)[atomBIndex];
   return this->GetDistanceAtoms(atomA, atomB);
}

double Molecule::GetDistanceAtoms(const Atom& atomA, const Atom& atomB) const{

   double distance=0.0;
   distance = sqrt( pow(atomA.GetXyz()[0] - atomB.GetXyz()[0], 2.0)
                   +pow(atomA.GetXyz()[1] - atomB.GetXyz()[1], 2.0)
                   +pow(atomA.GetXyz()[2] - atomB.GetXyz()[2], 2.0) );
   return distance;

}

}





