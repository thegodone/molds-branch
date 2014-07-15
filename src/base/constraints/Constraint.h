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
#ifndef INCLUDED_CONSTRAIN
#define INCLUDED_CONSTRAIN
namespace MolDS_base_constraints{

class Constraint : public MolDS_base::PrintController{
public:
   Constraint(const MolDS_base::Molecule* molecule,
              const boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure);
   virtual ~Constraint(){};
   virtual void                 SetConstraintCondition() = 0;
   virtual double const* const* GetForce(int elecState) = 0;
   inline MolDS_base::ConstraintType GetType() const{return this->type;}
protected:
   Constraint(){};
   MolDS_base::ConstraintType type;
   const MolDS_base::Molecule* molecule;
   const MolDS_base::Molecule* refMolecule;
   double** constrainedMatrixForce;
   const boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure;
private:
   
};

}
#endif



