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
#ifndef INCLUDED_LAPACK
#define INCLUDED_LAPACK
#include"mkl.h"
namespace MolDS_wrappers{
// Lapacke is singleton
class Lapack: public MolDS_base::PrintController, private MolDS_base::Uncopyable{
public:
   static Lapack* GetInstance();
   static void DeleteInstance();
   int Dsyevd(double** matrix, double* eigenValues, int size, bool calcEigenVectors);
   int Dsysv(double const* const* matrix, double* b, int size);
private:
   Lapack();
   ~Lapack();
   static Lapack* lapack;
   bool calculatedDsysvBlockSize;
   int dsysvBlockSize;
   std::string errorMessageDsyevdInfo;
   std::string errorMessageDsyevdSize;
   std::string errorMessageDsysvInfo;
   std::string errorMessageInfo;
   std::string errorMessageDsysvSize;
};
}
#endif
