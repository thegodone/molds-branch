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
#ifndef INCLUDED_MOLOGGER
#define INCLUDED_MOLOGGER
namespace MolDS_base_loggers{

class MOLogger: public MolDS_base::PrintController{
public:
   MOLogger(const MolDS_base::Molecule& molecule, double const* const* fockMatrix, MolDS_base::TheoryType theory);
   void DrawMO(int moIndex);
   void DrawMO(std::vector<int> moIndeces);
private:
   std::string errorMessageFockMatrixNULL;
   std::string messageCubeHeaderComment1;
   std::string messageCubeHeaderComment2;
   std::string messageStartMOPlot;
   std::string messageEndMOPlot;
   std::string messageSkippedMOIndex;
   std::string messageOmpElapsedTimeMOPlot;
   std::string messageUnitSec;
   std::string stringCubeExtension;
   MOLogger();
   MolDS_base::Molecule const* molecule;
   double const* const* fockMatrix;
   MolDS_base::TheoryType theory;
   void MatricesNullCheck() const;
   void SetMessages();
   void CalcGridDisplacement(double* dx, double* dy, double* dz) const;
   void CalcOrigin(double* origin) const;
   std::string GetFileName(int moIndex, int digit) const;
   void OutputHeaderToFile(std::ofstream& ofs,
                           double const* origin,
                           double dx,
                           double dy,
                           double dz) const;
   void OutputMoleculeToFile(std::ofstream& ofs, const MolDS_base::Molecule& molecule) const;
   double GetMoValue(int moIndex, 
                     const MolDS_base::Molecule& molecule, 
                     double const* const* fockMatrix, 
                     double x, 
                     double y, 
                     double z) const;
};

}
#endif
