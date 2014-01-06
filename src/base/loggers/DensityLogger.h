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
#ifndef INCLUDED_DENSITY_LOGGER
#define INCLUDED_DENSITY_LOGGER
namespace MolDS_base_loggers{

class DensityLogger: public MolDS_base::PrintController{
public:
   DensityLogger(const MolDS_base::Molecule& molecule, 
                 double const* const* fockMatrix, 
                 double const* const* cisMatrix, 
                 MolDS_base::TheoryType theory);
   virtual ~DensityLogger();
   void DrawDensity(int elecStateIndex) const;
   void DrawDensity(std::vector<int> elecStateIndeces) const;
protected:
   std::string errorMessageCISMatrixNULL;
   std::string errorMessageFockMatrixNULL;
   std::string messageCubeHeaderComment1;
   std::string messageStartDensityPlot;
   std::string messageEndDensityPlot;
   std::string messageOmpElapsedTimeDensityPlot;
   std::string stringCubeExtension;
   DensityLogger();
   void SetMessages();
   double GetMOValue(int moIndex, 
                     const MolDS_base::Molecule& molecule, 
                     double const* const* fockmatrix,
                     double x, 
                     double y, 
                     double z) const;
   virtual std::string GetFileName(int elecStateIndex, int digit) const = 0;
   virtual double GetDensityValue(int elecStateIndex, 
                                  double const* const* const* const* activeOccMOs,
                                  double const* const* const* const* activeVirMOs,
                                  double const* const* cisMatrix,
                                  int ix, 
                                  int iy, 
                                  int iz) const =0;
   virtual double const* GetFrameLength() const =0;
   virtual int const* GetGridNumber() const =0;
private:
   std::string messageCubeHeaderComment2;
   std::string messageSkippedElecStateIndex;
   std::string messageUnitSec; 
   MolDS_base::Molecule const* molecule;
   double const* const* fockMatrix;
   double const* const* cisMatrix;
   MolDS_base::TheoryType theory;
   void MatricesNullCheck() const;
   void CalcGridDisplacement(double* dx, double* dy, double* dz) const;
   void CalcOrigin(double* origin) const;
   void OutputHeaderToFile(std::ofstream& ofs, 
                           double const* origin, 
                           double dx, 
                           double dy,
                           double dz) const;
   void OutputMoleculeToFile(std::ofstream& ofs, const MolDS_base::Molecule& molecule)const ;
   void CalcActiveMOs(double**** activeOccMOs, 
                      double**** activeVirMOs,
                      double dx, double dy, double dz,
                      double const* origin,
                      const MolDS_base::Molecule& molecule, 
                      double const* const* fockMatrix,
                      double const* const* cisMatrix) const;
   void MallocTemporaryActiveMOs(double***** activeOccMOs, double***** activeVirMOs) const;
   void FreeTemporaryActiveMOs(double***** activeOccMOs, double***** activeVirMOs) const;
};

}
#endif
