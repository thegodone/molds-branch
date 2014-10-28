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
#ifndef INCLUDED_SCALAPACK
#define INCLUDED_SCALAPACK
namespace MolDS_wrappers{
//typedef intptr_t molds_scalapack_int;
#ifdef __FCC_VERSION
typedef int molds_scalapack_int;
#else
typedef intptr_t molds_scalapack_int;
#endif
// ScaLapacke is singleton
class ScaLapack: public MolDS_base::PrintController, private MolDS_base::Uncopyable{
public:
   static ScaLapack* GetInstance();
   static void DeleteInstance();
   molds_scalapack_int Pdsyevd(double** matrix, double* eigenValues, molds_scalapack_int size, bool calcEigenVectors);
private:
   ScaLapack();
   ~ScaLapack();
   static ScaLapack* scaLapack;
   std::string errorMessagePdsyevdInfo;
   std::string errorMessagePdsyevdSize;
   std::string errorMessagePdsyevdNotSupported;
};
}
#endif
