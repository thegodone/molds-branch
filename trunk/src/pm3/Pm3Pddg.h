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
#ifndef INCLUDED_PM3PDDG
#define INCLUDED_PM3PDDG
namespace MolDS_pm3{

/***
 *  Main References for PM3/PDDG are [RCJ_2002, BGRJ_2003, and BGJ_2003]
 */
class Pm3Pddg : public MolDS_pm3::Pm3{
public:
   Pm3Pddg();
   virtual ~Pm3Pddg();
protected:
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual double GetDiatomCoreRepulsionEnergy(const MolDS_base_atoms::Atom& atomA,
                                               const MolDS_base_atoms::Atom& atomB) const;
   virtual double GetDiatomCoreRepulsion1stDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                      const MolDS_base_atoms::Atom& atomB, 
                                                      MolDS_base::CartesianType axisA) const;
   virtual double GetDiatomCoreRepulsion2ndDerivative(const MolDS_base_atoms::Atom& atomA, 
                                                      const MolDS_base_atoms::Atom& atomB, 
                                                      MolDS_base::CartesianType axisA1,
                                                      MolDS_base::CartesianType axisA2) const;
private:
   double GetPddgAdditonalDiatomCoreRepulsionTerm(int na, double pa, double da, 
                                                  int nb, double pb, double db,
                                                  double distance) const;
   double GetPddgAdditonalDiatomCoreRepulsionTerm1stDerivative(int na, double pa, double da, 
                                                               int nb, double pb, double db,
                                                               double distance) const;
   double GetPddgAdditonalDiatomCoreRepulsionTerm2ndDerivative(int na, double pa, double da, 
                                                               int nb, double pb, double db,
                                                               double distance) const;
                                                
};

}
#endif



