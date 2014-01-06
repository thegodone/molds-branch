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
#ifndef INCLUDED_HOLE_DENSITY_LOGGER
#define INCLUDED_HOLE_DENSITY_LOGGER
namespace MolDS_base_loggers{

class HoleDensityLogger: public DensityLogger{
public:
   HoleDensityLogger(const MolDS_base::Molecule& molecule, 
                     double const* const* fockMatrix, 
                     double const* const* cisMatrix, 
                     MolDS_base::TheoryType theory);
   ~HoleDensityLogger();
protected:
   void SetMessages();
   double GetDensityValue(int elecStateIndex, 
                          double const* const* const* const* activeOccMOs,
                          double const* const* const* const* activeVirMOs,
                          double const* const* cisMatrix,
                          int ix, 
                          int iy, 
                          int iz) const;
   std::string GetFileName(int elecStateIndex, int digit) const;
   double const* GetFrameLength() const{return MolDS_base::Parameters::GetInstance()->GetFrameLengthHolePlot();}
   int const*    GetGridNumber()  const{return MolDS_base::Parameters::GetInstance()->GetGridNumberHolePlot();}
private:
   HoleDensityLogger();
};

}
#endif
