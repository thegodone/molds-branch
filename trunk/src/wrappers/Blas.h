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
#ifndef INCLUDED_BLAS
#define INCLUDED_BLAS
namespace MolDS_wrappers{
// Blas is singleton
class Blas: public MolDS_base::PrintController, private MolDS_base::Uncopyable{
public:
   static Blas* GetInstance();
   static void DeleteInstance();
   void Dcopy(int n,
              double const* vectorX,
              double*       vectorY) const;
   void Dcopy(int n,
              double const* vectorX, int incrementX,
              double*       vectorY, int incrementY) const;
   void Daxpy(int n, double alpha,
              double const* vectorX,
              double*       vectorY) const;
   void Daxpy(int n, double alpha,
              double const* vectorX, int incrementX,
              double*       vectorY, int incrementY) const;
   double Ddot(int n,
               double const* vectorX,
               double const* vectorY) const;
   double Ddot(int n,
               double const* vectorX, int incrementX,
               double const* vectorY, int incrementY)const;
   void Dgemv(int m, int n,
              double const* const* matrixA,
              double const* vectorX,
              double* vectorY) const;
   void Dgemv(bool isColumnMajorMatrixA,
              int m, int n,
              double alpha,
              double const* const* matrixA,
              double const* vectorX,
              int incrementX,
              double beta,
              double* vectorY,
              int incrementY) const;
   void Dgemm(int m, int n, int k, 
              double const* const* matrixA, 
              double const* const* matrixB, 
              double**             matrixC) const;
   void Dgemm(bool isColumnMajorMatrixA, 
              bool isColumnMajorMatrixB, 
              int m, int n, int k, 
              double alpha,
              double const* const* matrixA,
              double const* const* matrixB,
              double beta,
              double**             matrixC) const;
private:
   Blas();
   ~Blas();
   static Blas* blas;
};
}
#endif
