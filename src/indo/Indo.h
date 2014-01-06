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
#ifndef INCLUDED_INDO
#define INCLUDED_INDO
namespace MolDS_indo{

/***
 *  References for Indo are [PB_1970] and [PS_1966].
 */
class Indo : public MolDS_cndo::Cndo2{
public:
   Indo();
   virtual ~Indo();
protected:
   virtual void SetMessages();
   virtual void SetEnableAtomTypes();
   virtual double GetFockDiagElement(const MolDS_base_atoms::Atom& atomA, 
                                     int indexAtomA, 
                                     int mu, 
                                     const MolDS_base::Molecule& molecule, 
                                     double const* const* gammaAB,
                                     double const* const* orbitalElectronPopulation, 
                                     double const* atomicElectronPopulation,
                                     double const* const* const* const* const* const* twoElecsTwoAtomCores,
                                     bool isGuess) const;
   virtual double GetFockOffDiagElement(const MolDS_base_atoms::Atom& atomA, 
                                        const MolDS_base_atoms::Atom& atomB, 
                                        int indexAtomA, 
                                        int indexAtomB, 
                                        int mu, 
                                        int nu, 
                                        const MolDS_base::Molecule& molecule, 
                                        double const* const* gammaAB, 
                                        double const* const* overelap,
                                        double const* const* orbitalElectronPopulation,
                                        double const* const* const* const* const* const* twoElecsTwoAtomCores,
                                        bool isGuess) const;
   virtual double GetMolecularIntegralElement(int moI, int moJ, int moK, int moL, 
                                              const MolDS_base::Molecule& molecule, 
                                              double const* const* fockMatrix, 
                                              double const* const* gammaAB) const;
private:
   double GetCoulombInt(MolDS_base::OrbitalType orbital1, 
                        MolDS_base::OrbitalType orbital2, 
                        double gamma, 
                        const MolDS_base_atoms::Atom& atom) const; // Indo Coulomb Interaction, (3.87) - (3.91) in J. A. Pople book.
   double GetExchangeInt(MolDS_base::OrbitalType orbital1, 
                         MolDS_base::OrbitalType orbital2, 
                         double gamma, 
                         const MolDS_base_atoms::Atom& atom) const; // Indo Exchange Interaction, (3.87) - (3.91) in J. A. Pople book.
};

}
#endif



