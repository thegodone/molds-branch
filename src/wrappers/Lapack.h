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
#ifndef INCLUDED_LAPACK
#define INCLUDED_LAPACK
namespace MolDS_wrappers{
//typedef intptr_t molds_lapack_int;
#ifdef __FCC_VERSION
typedef int molds_lapack_int;
#else
typedef intptr_t molds_lapack_int;
#endif
// Lapacke is singleton
class Lapack: public MolDS_base::PrintController, private MolDS_base::Uncopyable{
public:
   static Lapack* GetInstance();
   static void DeleteInstance();
   molds_lapack_int Dsyevd(double** matrix, double* eigenValues, molds_lapack_int size, bool calcEigenVectors);
   molds_lapack_int Dsysv(double** matrix, double* b, molds_lapack_int size);
   molds_lapack_int Dgetrs(double** matrix, double** b, molds_lapack_int size, molds_lapack_int nrhs) const;
   molds_lapack_int Dgetrf(double** matrix, molds_lapack_int sizeM, molds_lapack_int sizeN) const;
   molds_lapack_int Dgetrf(double** matrix, molds_lapack_int* ipiv, molds_lapack_int sizeM, molds_lapack_int sizeN) const;
private:
   Lapack();
   ~Lapack();
   static Lapack* lapack;
   std::string errorMessageDsyevdInfo;
   std::string errorMessageDsyevdSize;
   std::string errorMessageDsysvInfo;
   std::string errorMessageDsysvSize;
   std::string errorMessageDgetrsInfo;
   std::string errorMessageDgetrsSize;
   std::string errorMessageDgetrfInfo;
   molds_lapack_int Dgetrf(double* matrix, molds_lapack_int* ipiv, molds_lapack_int sizeM, molds_lapack_int sizeN) const;
};
}
#endif
