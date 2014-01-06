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
#ifndef INCLUDED_MOLDS
#define INCLUDED_MOLDS
namespace MolDS_base{
class MolDS: public PrintController{
public:
   void Run(int argc, char *argv[]);
private:
   time_t  startTime;
   clock_t startTick;
   double  ompStartTime;
   void Initialize();
   void Finalize(bool runsNormally) const;
   void CalculateElectronicStructureOnce(Molecule* molecule, bool* runsNormally) const;
   void DoMC(Molecule* molecule, bool* runsNormally) const;
   void DoMD(Molecule* molecule, bool* runsNormally) const;
   void DoRPMD(Molecule* molecule, bool* runsNormally) const;
   void DoNASCO(Molecule* molecule, bool* runsNormally) const;
   void OptimizeGeometry(Molecule* molecule, bool* runsNormally) const;
   void DiagonalizePrincipalAxes(Molecule* molecule, bool* runsNormally) const;
   void TranslateMolecule(Molecule* molecule, bool* runsNormally) const;
   void RotateMolecule(Molecule* molecule, bool* runsNormally) const;
};
}
#endif



