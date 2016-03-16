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
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<sstream>
#include<math.h>
#include<string>
#include<vector>
#include<stdexcept>
#include<algorithm>
#include<omp.h>
#include<boost/format.hpp>
#include"Enums.h"
#include"Uncopyable.h"
#include"PrintController.h"
#include"MolDSException.h"
#include"MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
#include"../wrappers/Lapack.h"
#include"MathUtilities.h"
#include"EularAngle.h"
#include"Parameters.h"
#include"RealSphericalHarmonicsIndex.h"
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
   double*  oldXyzCOM = this->xyzCOM;
   double*  oldXyzCOC = this->xyzCOC;
   double** oldDistanceAtoms     = this->distanceAtoms;
   double** oldDistanceEpcs      = this->distanceEpcs;
   double** oldDistanceAtomsEpcs = this->distanceAtomsEpcs;
   vector<Atom*>* oldAtomVect     = this->atomVect;
   vector<Atom*>* oldRealAtomVect = this->realAtomVect;
   vector<Atom*>* oldGhostAtomVect= this->ghostAtomVect;
   vector<Atom*>* oldEpcVect      = this->epcVect;
   this->CopyInitialize(rhs);
   this->Finalize(&oldAtomVect, 
                  &oldRealAtomVect, 
                  &oldGhostAtomVect, 
                  &oldEpcVect, 
                  &oldXyzCOM, 
                  &oldXyzCOC, 
                  &oldDistanceAtoms,
                  &oldDistanceEpcs,
                  &oldDistanceAtomsEpcs);
   return *this;
}

Molecule::~Molecule(){
   this->Finalize(&this->atomVect, 
                  &this->realAtomVect, 
                  &this->ghostAtomVect, 
                  &this->epcVect, 
                  &this->xyzCOM, 
                  &this->xyzCOC, 
                  &this->distanceAtoms,
                  &this->distanceEpcs,
                  &this->distanceAtomsEpcs);
   //this->OutputLog("molecule deleted\n");
}

void Molecule::CopyInitialize(const Molecule& rhs){
   this->Initialize();
   for(int i=0; i<CartesianType_end; i++){
      this->xyzCOM[i] = rhs.xyzCOM[i];
      this->xyzCOC[i] = rhs.xyzCOC[i];
   }
   this->totalNumberAOs = rhs.totalNumberAOs;
   this->totalNumberValenceElectrons = rhs.totalNumberValenceElectrons;
   this->totalCoreMass = rhs.totalCoreMass;
   if(rhs.realAtomVect != NULL){
      int realAtomNum = rhs.realAtomVect->size();
      for(int i=0; i<realAtomNum; i++){
         Atom* atom = (*rhs.realAtomVect)[i];
         this->realAtomVect->push_back(AtomFactory::Create(atom->GetAtomType(),
                                                           atom->GetIndex(),
                                                           atom->GetXyz()[XAxis],
                                                           atom->GetXyz()[YAxis],
                                                           atom->GetXyz()[ZAxis],
                                                           atom->GetPxyz()[XAxis],
                                                           atom->GetPxyz()[YAxis],
                                                           atom->GetPxyz()[ZAxis]));
         (*this->realAtomVect)[i]->SetFirstAOIndex(atom->GetFirstAOIndex());
      }
   }
   if(rhs.ghostAtomVect != NULL){
      int ghostAtomNum = rhs.ghostAtomVect->size();
      for(int i=0; i<ghostAtomNum; i++){
         Atom* atom = (*rhs.ghostAtomVect)[i];
         this->ghostAtomVect->push_back(AtomFactory::Create(atom->GetAtomType(),
                                                         atom->GetIndex(),
                                                         atom->GetXyz()[XAxis],
                                                         atom->GetXyz()[YAxis],
                                                         atom->GetXyz()[ZAxis],
                                                         atom->GetPxyz()[XAxis],
                                                         atom->GetPxyz()[YAxis],
                                                         atom->GetPxyz()[ZAxis]));
         (*this->realAtomVect)[i]->SetFirstAOIndex(atom->GetFirstAOIndex());
      }
   }
   if(rhs.atomVect != NULL){
      this->CopyRealGhostAtom2Atom();
      int atomNum = this->atomVect->size();
      MallocerFreer::GetInstance()->Malloc<double>(&this->distanceAtoms, atomNum, atomNum);
      for(int i=0; i<atomNum; i++){
         for(int j=0; j<atomNum; j++){
            this->distanceAtoms[i][j] = rhs.distanceAtoms[i][j];
         }
      }
   }
   if(rhs.epcVect != NULL){
      int epcNum = rhs.epcVect->size();
      for(int i=0; i<epcNum; i++){
         Atom* epc = (*rhs.epcVect)[i];
         this->epcVect->push_back(
            AtomFactory::Create(EPC,
                                epc->GetIndex(),
                                epc->GetXyz()[XAxis],
                                epc->GetXyz()[YAxis],
                                epc->GetXyz()[ZAxis],
                                epc->GetPxyz()[XAxis],
                                epc->GetPxyz()[YAxis],
                                epc->GetPxyz()[ZAxis],
                                epc->GetCoreCharge()));
         (*this->epcVect)[i]->SetFirstAOIndex(epc->GetFirstAOIndex());
      }
      MallocerFreer::GetInstance()->Malloc<double>(&this->distanceEpcs, epcNum, epcNum);
      for(int i=0; i<epcNum; i++){
         for(int j=0; j<epcNum; j++){
            this->distanceEpcs[i][j] = rhs.distanceEpcs[i][j];
         }
      }
   }
   if(rhs.atomVect != NULL && rhs.epcVect != NULL){
      int atomNum = rhs.atomVect->size();
      int epcNum = rhs.epcVect->size();
      MallocerFreer::GetInstance()->Malloc<double>(&this->distanceAtomsEpcs, atomNum, epcNum);
      for(int i=0; i<atomNum; i++){
         for(int j=0; j<epcNum; j++){
            this->distanceAtomsEpcs[i][j] = rhs.distanceAtomsEpcs[i][j];
         }
      }
   }
}

void Molecule::Initialize(){
   this->SetMessages();
   this->xyzCOM             = NULL;
   this->xyzCOC             = NULL;
   this->distanceAtoms      = NULL;
   this->distanceEpcs       = NULL;
   this->distanceAtomsEpcs  = NULL;
   this->atomVect           = NULL;
   this->realAtomVect       = NULL;
   this->ghostAtomVect         = NULL;
   this->epcVect            = NULL;
   try{
      this->atomVect     = new vector<Atom*>;
      this->realAtomVect = new vector<Atom*>;
      this->ghostAtomVect= new vector<Atom*>;
      this->epcVect      = new vector<Atom*>;
      MallocerFreer::GetInstance()->Malloc<double>(&this->xyzCOM, CartesianType_end);
      MallocerFreer::GetInstance()->Malloc<double>(&this->xyzCOC, CartesianType_end);
   }
   catch(exception ex){
      this->Finalize(&this->atomVect, 
                     &this->realAtomVect, 
                     &this->ghostAtomVect, 
                     &this->epcVect, 
                     &this->xyzCOM, 
                     &this->xyzCOC, 
                     &this->distanceAtoms,
                     &this->distanceEpcs,
                     &this->distanceAtomsEpcs);
      throw MolDSException(ex.what());
   }
}

void Molecule::Finalize(vector<Atom*>** atomVect, 
                        vector<Atom*>** realAtomVect, 
                        vector<Atom*>** ghostAtomVect, 
                        vector<Atom*>** epcVect, 
                        double** xyzCOM, 
                        double** xyzCOC, 
                        double*** distanceAtoms,
                        double*** distanceEpcs,
                        double*** distanceAtomsEpcs){
   MallocerFreer::GetInstance()->Free<double>(xyzCOM, CartesianType_end);
   MallocerFreer::GetInstance()->Free<double>(xyzCOC, CartesianType_end);
   int atomNum=0;
   int epcNum = 0;
   if(*atomVect != NULL && *epcVect != NULL){
      atomNum = (*atomVect)->size();
      epcNum = (*epcVect)->size();
      MallocerFreer::GetInstance()->Free<double>(distanceAtomsEpcs, atomNum, epcNum);
   }
   if(*atomVect != NULL){
      (*atomVect)->clear();
      delete *atomVect;
      *atomVect = NULL;
      //this->OutputLog("atomVect deleted\n");
      MallocerFreer::GetInstance()->Free<double>(distanceAtoms, atomNum, atomNum);
   }
   if(*realAtomVect != NULL){
      int realAtomNum = (*realAtomVect)->size();
      for(int i=0; i<realAtomNum; i++){
         if((**realAtomVect)[i] != NULL){
            delete (**realAtomVect)[i];
            (**realAtomVect)[i] = NULL;
         }
      }
      (*realAtomVect)->clear();
      delete *realAtomVect;
      *realAtomVect = NULL;
      //this->OutputLog("realAtomVect deleted\n");
   }
   if(*ghostAtomVect != NULL){
      int ghostAtomNum = (*ghostAtomVect)->size();
      for(int i=0; i<ghostAtomNum; i++){
         if((**ghostAtomVect)[i] != NULL){
            delete (**ghostAtomVect)[i];
            (**ghostAtomVect)[i] = NULL;
         }
      }
      (*ghostAtomVect)->clear();
      delete *ghostAtomVect;
      *ghostAtomVect = NULL;
      //this->OutputLog("ghostAtomVect deleted\n");
   }
   if(*epcVect != NULL){
      epcNum = (*epcVect)->size();
      for(int i=0; i<epcNum; i++){
         if((**epcVect)[i] != NULL){
            delete (**epcVect)[i];
            (**epcVect)[i] = NULL;
         }
      }
      (*epcVect)->clear();
      delete *epcVect;
      *epcVect = NULL;
      //this->OutputLog("epcVect deleted\n");
      MallocerFreer::GetInstance()->Free<double>(distanceEpcs, epcNum, epcNum);
   }
}

void Molecule::SetMessages(){
   this->errorMessageGetAtomVectNull  = "Error in base::Molecule::GetAtomVect: atomVect is NULL.\n";
   this->errorMessageGetRealAtomVectNull  = "Error in base::Molecule::GetRealAtomVect: realAtomVect is NULL.\n";
   this->errorMessageGetGhostAtomVectNull  = "Error in base::Molecule::GetGhostAtomVect: ghostAtomVect is NULL.\n";
   this->errorMessageGetEPCVectNull  = "Error in base::Molecule::GetEpcVect: epcVect is NULL.\n";
   this->errorMessageAddAtomNull = "Error in base::Molecule::AddAtom: atomVect is NULL.\n";
   this->errorMessageAddRealAtomNull = "Error in base::Molecule::AddAtom: realAtomVect is NULL.\n";
   this->errorMessageAddGhostAtomNull = "Error in base::Molecule::AddAtom: ghostAtomVect is NULL.\n";
   this->errorMessageAddEPCNull  = "Error in base::Molecule::AddEnviromentalPointCharge: epcVect is NULL.\n";
   this->errorMessageCopyRealGhostAtom2AtomNotEmpty = "Error in base::Molecule::CopyRealGhostAtom2Atom: atomVect is not empty.\n";
   this->errorMessageGetXyzCOMNull = "Error in base::Molecule::GetXyzCOM: xyzCOM is NULL.\n";
   this->errorMessageGetXyzCOCNull = "Error in base::Molecule::GetXyzCOC: xyzCOC is NULL.\n";
   this->errorMessageCalcXyzCOMNull = "Error in base::Molecule::CalcXyzCOM: xyzCOM is NULL.\n";
   this->errorMessageCalcXyzCOCNull = "Error in base::Molecule::CalcXyzCOC: xyzCOC is NULL.\n";
   this->messageTotalNumberAOs = "\t\tTotal number of valence AOs: ";
   this->messageTotalNumberAtoms = "\t\tTotal number of atoms: ";
   this->messageTotalNumberValenceElectrons = "\t\tTotal number of valence electrons: ";
   this->messageAtomCoordinatesTitle = "\t\t\t\t\t| i-th | atom type |   x[a.u.]   |   y[a.u.]   |   z[a.u.]   |\t\t|  x[angst.]  |  y[angst.]  |  z[angst.]  |\n";
   this->messageAtomCoordinates = "\t\tAtom coordinates:";
   this->messageAtomMomenta = "\t\tAtom momenta:";
   this->messageAtomMomentaTitle = "\t\t\t\t| i-th | atom type |   px[a.u.]   |   py[a.u.]   |   pz[a.u.]   |\t\t|   px[u]   |   py[u]   |   pz[u]   | [u] = [(g/Mol)*(angst/fs)]\n";
   this->messageCOM = "\t\tCenter of Mass:";
   this->messageCOC = "\t\tCenter of Core:";
   this->messageCOMTitle = "\t\t\t\t|   x[a.u.]   |   y[a.u.]   |   z[a.u.]   |\t\t|  x[angst.]  |  y[angst.]  |  z[angst.]  |\n";
   this->messageEpcConfiguration = "\tEnvironmental Point Charge(EPC) configuration:\n";
   this->messageEpcCoordinates = "\t\tEPC coordinates:";
   this->messageEpcCoordinatesTitle = "\t\t\t\t\t| i-th |  charge[a.u.]  |   x[a.u.]   |   y[a.u.]   |   z[a.u.]   |\t\t|  x[angst.]  |  y[angst.]  |  z[angst.]  |\n";
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
#ifdef MOLDS_DBG
   if(this->atomVect==NULL) throw MolDSException(this->errorMessageAddAtomNull);
#endif
   this->atomVect->push_back(atom);
}

void Molecule::AddRealAtom(Atom* atom){
#ifdef MOLDS_DBG
   if(this->realAtomVect==NULL) throw MolDSException(this->errorMessageAddRealAtomNull);
#endif
   this->realAtomVect->push_back(atom);
}

void Molecule::AddGhostAtom(Atom* atom){
#ifdef MOLDS_DBG
   if(this->ghostAtomVect==NULL) throw MolDSException(this->errorMessageAddGhostAtomNull);
#endif
   this->ghostAtomVect->push_back(atom);
}

void Molecule::AddEpc(Atom* epc){
#ifdef MOLDS_DBG
   if(this->epcVect==NULL) throw MolDSException(this->errorMessageAddEPCNull);
#endif
   this->epcVect->push_back(epc);
}

double const* Molecule::GetXyzCOM() const{
#ifdef MOLDS_DBG
   if(this->xyzCOM==NULL) throw MolDSException(this->errorMessageGetXyzCOMNull);
#endif
   return this->xyzCOM;
}

double const* Molecule::GetXyzCOC() const{
#ifdef MOLDS_DBG
   if(this->xyzCOC==NULL) throw MolDSException(this->errorMessageGetXyzCOCNull);
#endif
   return this->xyzCOC;
}

void Molecule::CalcXyzCOM(){
#ifdef MOLDS_DBG
   if(this->xyzCOM==NULL) throw MolDSException(this->errorMessageCalcXyzCOMNull);
#endif
   double  totalAtomicMass = 0.0;
   double* atomicXyz;
   double  atomicMass = 0.0;

   for(int j=0; j<3; j++){
      this->xyzCOM[j] = 0.0;
   }
      
   for(int i=0; i<this->realAtomVect->size(); i++){
      const Atom& atom = *(*this->realAtomVect)[i]; 
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
}

void Molecule::CalcXyzCOC(){
#ifdef MOLDS_DBG
   if(this->xyzCOC==NULL) throw MolDSException(this->errorMessageCalcXyzCOCNull);
#endif
   double  totalCoreMass = 0.0;
   double* atomicXyz;
   double  coreMass = 0.0;

   for(int j=0; j<3; j++){
      this->xyzCOC[j] = 0.0;
   }
      
   for(int i=0; i<this->realAtomVect->size(); i++){
      const Atom& atom = *(*this->realAtomVect)[i]; 
      atomicXyz = atom.GetXyz();
      coreMass  = atom.GetCoreMass();
      totalCoreMass += coreMass;
      for(int j=0; j<3; j++){
         this->xyzCOC[j] += atomicXyz[j] * coreMass;
      }
   }
   for(int i=0; i<3; i++){
      this->xyzCOC[i]/=totalCoreMass;
   }
}

void Molecule::CalcDistanceAtoms(){
   if(this->distanceAtoms==NULL){
      MallocerFreer::GetInstance()->Malloc<double>(&this->distanceAtoms, this->atomVect->size(), this->atomVect->size());
   }
   for(int a=0; a<this->atomVect->size(); a++){
      const Atom& atomA = *(*this->atomVect)[a];
      for(int b=a; b<this->atomVect->size(); b++){
         const Atom& atomB = *(*this->atomVect)[b];
         double distance=0.0;
         distance = sqrt( pow(atomA.GetXyz()[0] - atomB.GetXyz()[0], 2.0)
                         +pow(atomA.GetXyz()[1] - atomB.GetXyz()[1], 2.0)
                         +pow(atomA.GetXyz()[2] - atomB.GetXyz()[2], 2.0) );
         this->distanceAtoms[a][b] = distance;
         this->distanceAtoms[b][a] = distance;
      }
   }
}

void Molecule::CalcDistanceEpcs(){
   if(this->epcVect == NULL){return;}
   if(this->distanceEpcs==NULL){
      MallocerFreer::GetInstance()->Malloc<double>(&this->distanceEpcs, this->epcVect->size(), this->epcVect->size());
   }
   for(int a=0; a<this->epcVect->size(); a++){
      const Atom& epcA = *(*this->epcVect)[a];
      for(int b=a; b<this->epcVect->size(); b++){
         const Atom& epcB = *(*this->epcVect)[b];
         double distance=0.0;
         distance = sqrt( pow(epcA.GetXyz()[0] - epcB.GetXyz()[0], 2.0)
                         +pow(epcA.GetXyz()[1] - epcB.GetXyz()[1], 2.0)
                         +pow(epcA.GetXyz()[2] - epcB.GetXyz()[2], 2.0) );
         this->distanceEpcs[a][b] = distance;
         this->distanceEpcs[b][a] = distance;
      }
   }
}

void Molecule::CalcDistanceAtomsEpcs(){
   if(this->epcVect == NULL){return;}
   if(this->distanceAtomsEpcs==NULL){
      MallocerFreer::GetInstance()->Malloc<double>(&this->distanceAtomsEpcs, this->atomVect->size(), this->epcVect->size());
   }
   for(int a=0; a<this->atomVect->size(); a++){
      const Atom& atom = *(*this->atomVect)[a];
      for(int b=0; b<this->epcVect->size(); b++){
         const Atom& epc = *(*this->epcVect)[b];
         double distance=0.0;
         distance = sqrt( pow(atom.GetXyz()[0] - epc.GetXyz()[0], 2.0)
                         +pow(atom.GetXyz()[1] - epc.GetXyz()[1], 2.0)
                         +pow(atom.GetXyz()[2] - epc.GetXyz()[2], 2.0) );
         this->distanceAtomsEpcs[a][b] = distance;
      }
   }
}

void Molecule::CalcBasics(){
   this->CopyRealGhostAtom2Atom();
   this->CalcTotalNumberAOs();
   this->CalcTotalNumberValenceElectrons();
   this->CalcTotalCoreMass();
   this->CalcBasicsConfiguration();
}

void Molecule::CalcBasicsConfiguration(){
   this->CalcXyzCOM();
   this->CalcXyzCOC();
   this->CalcDistanceAtoms();
   this->CalcDistanceEpcs();
   this->CalcDistanceAtomsEpcs();
}

void Molecule::CopyRealGhostAtom2Atom(){
   if(!this->atomVect->empty()){
      throw MolDSException(this->errorMessageCopyRealGhostAtom2AtomNotEmpty);
   }
   int realAtomNum = this->realAtomVect->size();
   int ghostAtomNum   = this->ghostAtomVect->size();
   for(int i=0; i<realAtomNum; i++){
      Atom* atom = (*this->realAtomVect)[i];
      this->atomVect->push_back(atom);
   }
   for(int i=0; i<ghostAtomNum; i++){
      Atom* atom = (*this->ghostAtomVect)[i];
      this->atomVect->push_back(atom);
   }
   std::sort(this->realAtomVect->begin(),  this->realAtomVect->end(),  MolDS_base_atoms::LessAtomIndex());
   std::sort(this->ghostAtomVect->begin(), this->ghostAtomVect->end(), MolDS_base_atoms::LessAtomIndex());
   std::sort(this->atomVect->begin(),      this->atomVect->end(),      MolDS_base_atoms::LessAtomIndex());
}

void Molecule::CalcTotalNumberAOs(){
   this->totalNumberAOs = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      (*this->atomVect)[i]->SetFirstAOIndex(totalNumberAOs);
      this->totalNumberAOs += (*this->atomVect)[i]->GetValenceSize();
   }
}

void Molecule::CalcTotalNumberValenceElectrons(){
   this->totalNumberValenceElectrons = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      this->totalNumberValenceElectrons += (*this->atomVect)[i]->GetNumberValenceElectrons();
   }
}

void Molecule::CalcTotalCoreMass(){
   this->totalCoreMass = 0; 
   for(int i=0; i<this->atomVect->size(); i++){
      const Atom& atom = *(*this->atomVect)[i]; 
      double coreMass = atom.GetCoreMass();
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

void Molecule::OutputEpcs() const{
   if(this->epcVect == NULL || this->epcVect->size() <= 0) {return;}
   double ang2AU = Parameters::GetInstance()->GetAngstrom2AU();
   this->OutputLog(this->messageEpcConfiguration);
   this->OutputLog(this->messageEpcCoordinatesTitle);
   for(int a=0; a<this->epcVect->size(); a++){
      const Atom& atom = *(*this->epcVect)[a];
      this->OutputLog(boost::format("%s\t%d\t%e\t%e\t%e\t%e\t\t%e\t%e\t%e\n") 
         % this->messageEpcCoordinates
         % a
         % atom.GetCoreCharge()
         % atom.GetXyz()[0]
         % atom.GetXyz()[1]
         % atom.GetXyz()[2]
         % (atom.GetXyz()[0]/ang2AU)
         % (atom.GetXyz()[1]/ang2AU)
         % (atom.GetXyz()[2]/ang2AU));
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
   this->CalcXyzCOM();
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
   for(int a=0; a<this->realAtomVect->size(); a++){
      const Atom& atom = *(*this->realAtomVect)[a];
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
   this->CalcXyzCOM();
   double rotatingOrigin[3] = {this->xyzCOM[0], this->xyzCOM[1], this->xyzCOM[2]};
   if(Parameters::GetInstance()->GetRotatingOrigin() != NULL){
      rotatingOrigin[0] = Parameters::GetInstance()->GetRotatingOrigin()[0];
      rotatingOrigin[1] = Parameters::GetInstance()->GetRotatingOrigin()[1];
      rotatingOrigin[2] = Parameters::GetInstance()->GetRotatingOrigin()[2];
   }

   RotatingType  rotatingType        = Parameters::GetInstance()->GetRotatingType();
   const double* rotatingAxis        = Parameters::GetInstance()->GetRotatingAxis();
   EularAngle    rotatingEularAngles = Parameters::GetInstance()->GetRotatingEularAngles();
   double        rotatingAngle       = Parameters::GetInstance()->GetRotatingAngle();

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
   this->CalcXyzCOM();
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

void Molecule::SynchronizeConfigurationTo(const Molecule& ref){
   for(int a=0; a<this->GetAtomVect().size(); a++){
      Atom& atom = *this->GetAtomVect()[a];
      const Atom& refAtom = *ref.GetAtomVect()[a];
      for(int i=0; i<CartesianType_end; i++){
         atom.GetXyz()[i] = refAtom.GetXyz()[i];
      }
   }
   this->CalcBasicsConfiguration();
}

void Molecule::SynchronizeMomentaTo(const Molecule& ref){
   for(int a=0; a<this->GetAtomVect().size(); a++){
      Atom& atom = *this->GetAtomVect()[a];
      const Atom& refAtom = *ref.GetAtomVect()[a];
      for(int i=0; i<CartesianType_end; i++){
         atom.GetPxyz()[i] = refAtom.GetPxyz()[i];
      }
   }
}

void Molecule::SynchronizePhaseSpacePointTo(const Molecule& ref){
   for(int a=0; a<this->GetAtomVect().size(); a++){
      Atom& atom = *this->GetAtomVect()[a];
      const Atom& refAtom = *ref.GetAtomVect()[a];
      for(int i=0; i<CartesianType_end; i++){
         atom.GetXyz() [i] = refAtom.GetXyz() [i];
         atom.GetPxyz()[i] = refAtom.GetPxyz()[i];
      }
   }
   this->CalcBasicsConfiguration();
}

void Molecule::BroadcastConfigurationToAllProcesses(int root) const{
   int numTransported = this->GetAtomVect().size()*CartesianType_end;
   double* tmp=NULL;
   try{
      MolDS_base::MallocerFreer::GetInstance()->Malloc<double>(&tmp, numTransported);
      for(int a=0; a<this->GetAtomVect().size(); a++){
         Atom& atom = *this->GetAtomVect()[a];
         for(int i=0; i<CartesianType_end; i++){
            tmp[a*CartesianType_end+i] = atom.GetXyz()[i];
         }
      }
      MolDS_mpi::MpiProcess::GetInstance()->Broadcast(tmp, numTransported, root);
      for(int a=0; a<this->GetAtomVect().size(); a++){
         Atom& atom = *this->GetAtomVect()[a];
         for(int i=0; i<CartesianType_end; i++){
            atom.GetXyz()[i] = tmp[a*CartesianType_end+i];
         }
      }
   }
   catch(MolDS_base::MolDSException ex){
      MolDS_base::MallocerFreer::GetInstance()->Free<double>(&tmp, numTransported);
      throw ex;
   }
   MolDS_base::MallocerFreer::GetInstance()->Free<double>(&tmp, numTransported);
}

void Molecule::BroadcastMomentaToAllProcesses(int root) const{
   int numTransported = this->GetAtomVect().size()*CartesianType_end;
   double* tmp=NULL;
   try{
      MolDS_base::MallocerFreer::GetInstance()->Malloc<double>(&tmp, numTransported);
      for(int a=0; a<this->GetAtomVect().size(); a++){
         Atom& atom = *this->GetAtomVect()[a];
         for(int i=0; i<CartesianType_end; i++){
            tmp[a*CartesianType_end+i] = atom.GetPxyz()[i];
         }
      }
      MolDS_mpi::MpiProcess::GetInstance()->Broadcast(tmp, numTransported, root);
      for(int a=0; a<this->GetAtomVect().size(); a++){
         Atom& atom = *this->GetAtomVect()[a];
         for(int i=0; i<CartesianType_end; i++){
            atom.GetPxyz()[i] = tmp[a*CartesianType_end+i];
         }
      }
   }
   catch(MolDS_base::MolDSException ex){
      MolDS_base::MallocerFreer::GetInstance()->Free<double>(&tmp, numTransported);
      throw ex;
   }
   MolDS_base::MallocerFreer::GetInstance()->Free<double>(&tmp, numTransported);
}

void Molecule::BroadcastPhaseSpacePointToAllProcesses(int root) const{
   int numTransported = 2*this->GetAtomVect().size()*CartesianType_end;
   double* tmp=NULL;
   try{
      MolDS_base::MallocerFreer::GetInstance()->Malloc<double>(&tmp, numTransported);
      for(int a=0; a<this->GetAtomVect().size(); a++){
         Atom& atom = *this->GetAtomVect()[a];
         for(int i=0; i<CartesianType_end; i++){
            int k = a*CartesianType_end+i;
            tmp[2*k  ] = atom.GetXyz() [i];
            tmp[2*k+1] = atom.GetPxyz()[i];
         }
      }
      MolDS_mpi::MpiProcess::GetInstance()->Broadcast(tmp, numTransported, root);
      for(int a=0; a<this->GetAtomVect().size(); a++){
         Atom& atom = *this->GetAtomVect()[a];
         for(int i=0; i<CartesianType_end; i++){
            int k = a*CartesianType_end+i;
            atom.GetXyz() [i] = tmp[2*k  ];
            atom.GetPxyz()[i] = tmp[2*k+1];
         }
      }
   }
   catch(MolDS_base::MolDSException ex){
      MolDS_base::MallocerFreer::GetInstance()->Free<double>(&tmp, numTransported);
      throw ex;
   }
   MolDS_base::MallocerFreer::GetInstance()->Free<double>(&tmp, numTransported);
}

}





