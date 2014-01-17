//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   //
// Copyright (C) 2012-2014 Katsuhiko Nishimra                             //
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
//typedef intptr_t molds_blas_int;
typedef intptr_t molds_blas_int;
// Blas is singleton
class Blas: public MolDS_base::PrintController, private MolDS_base::Uncopyable{
public:
   static Blas* GetInstance();
   static void DeleteInstance();
   void Dcopy(molds_blas_int n,
              double const* vectorX,
              double*       vectorY) const;
   void Dcopy(molds_blas_int n,
              double const* vectorX, molds_blas_int incrementX,
              double*       vectorY, molds_blas_int incrementY) const;
   void Daxpy(molds_blas_int n, double alpha,
              double const* vectorX,
              double*       vectorY) const;
   void Daxpy(molds_blas_int n, double alpha,
              double const* vectorX, molds_blas_int incrementX,
              double*       vectorY, molds_blas_int incrementY) const;
   double Ddot(molds_blas_int n,
               double const* vectorX,
               double const* vectorY) const;
   double Ddot(molds_blas_int n,
               double const* vectorX, molds_blas_int incrementX,
               double const* vectorY, molds_blas_int incrementY)const;
   double Dnrm2(molds_blas_int n,
                double const* vectorX) const;
   double Dnrm2(molds_blas_int n,
                double const* vectorX,
                molds_blas_int incrementX) const;
   double Damax(molds_blas_int n,
                      double const* vectorX) const;
   double Damax(molds_blas_int n,
                      double const* vectorX, molds_blas_int incrementX)const;
   void Dgemv(molds_blas_int m, molds_blas_int n,
              double const* const* matrixA,
              double const* vectorX,
              double* vectorY) const;
   void Dgemv(bool isColumnMajorMatrixA,
              molds_blas_int m, molds_blas_int n,
              double alpha,
              double const* const* matrixA,
              double const* vectorX,
              molds_blas_int incrementX,
              double beta,
              double* vectorY,
              molds_blas_int incrementY) const;
   void Dsymv(molds_blas_int n,
              double const* const* matrixA,
              double const* vectorX,
              double* vectorY) const;
   void Dsymv(molds_blas_int n, double alpha,
              double const* const* matrixA,
              double const* vectorX, molds_blas_int incrementX,
              double beta,
              double*       vectorY, molds_blas_int incrementY) const;
   void Dsyr(molds_blas_int n, double alpha,
             double const* vectorX,
             double ** matrixA)const;
   void Dsyr(molds_blas_int n, double alpha,
             double const* vectorX, molds_blas_int incrementX,
             double ** matrixA)const;
   void Dgemm(molds_blas_int m, molds_blas_int n, molds_blas_int k, 
              double const* const* matrixA, 
              double const* const* matrixB, 
              double**             matrixC) const;
   void Dgemm(bool isColumnMajorMatrixA, 
              bool isColumnMajorMatrixB, 
              molds_blas_int m, molds_blas_int n, molds_blas_int k, 
              double alpha,
              double const* const* matrixA,
              double const* const* matrixB,
              double beta,
              double**             matrixC) const;
   void Dgemm(bool isColumnMajorMatrixA, 
              bool isColumnMajorMatrixB, 
              bool isColumnMajorMatrixC, 
              molds_blas_int m, molds_blas_int n, molds_blas_int k, 
              double alpha,
              double const* const* matrixA,
              double const* const* matrixB,
              double beta,
              double**             matrixC) const;
   void Dgemm(bool isColumnMajorMatrixA, 
              bool isColumnMajorMatrixB, 
              molds_blas_int m, molds_blas_int n, molds_blas_int k, 
              double alpha,
              double const* const* matrixA,
              double const* const* matrixB,
              double beta,
              double**             matrixC,
              double*              tmpC) const;
   void Dgemm(bool isColumnMajorMatrixA, 
              bool isColumnMajorMatrixB, 
              bool isColumnMajorMatrixC, 
              molds_blas_int m, molds_blas_int n, molds_blas_int k, 
              double alpha,
              double const* const* matrixA,
              double const* const* matrixB,
              double beta,
              double**             matrixC,
              double*              tmpC) const;
   void Dgemmm(molds_blas_int m, molds_blas_int n, molds_blas_int k, molds_blas_int l,
               double const* const* matrixA, 
               double const* const* matrixB, 
               double const* const* matrixC, 
               double**             matrixD) const;
   void Dgemmm(bool isColumnMajorMatrixA, 
               bool isColumnMajorMatrixB, 
               bool isColumnMajorMatrixC, 
               molds_blas_int m, molds_blas_int n, molds_blas_int k, molds_blas_int l,
               double alpha,
               double const* const* matrixA,
               double const* const* matrixB,
               double const* const* matrixC,
               double beta,
               double**             matrixD) const;
   void Dgemmm(bool isColumnMajorMatrixA, 
               bool isColumnMajorMatrixB, 
               bool isColumnMajorMatrixC, 
               molds_blas_int m, molds_blas_int n, molds_blas_int k, molds_blas_int l,
               double alpha,
               double const* const* matrixA,
               double const* const* matrixB,
               double const* const* matrixC,
               double beta,
               double**             matrixD,
               double**             tmpMatrixBC) const;
   void Dgemmm(bool isColumnMajorMatrixA, 
               bool isColumnMajorMatrixB, 
               bool isColumnMajorMatrixC, 
               molds_blas_int m, molds_blas_int n, molds_blas_int k, molds_blas_int l,
               double alpha,
               double const* const* matrixA,
               double const* const* matrixB,
               double const* const* matrixC,
               double beta,
               double**             matrixD,
               double*              tmpVectorD,
               double**             tmpMatrixBC,
               double*              tmpVectorBC) const;
  void Dsyrk(molds_blas_int n, molds_blas_int k,
             double const *const* matrixA,
             double**             matrixC)const;
  void Dsyrk(molds_blas_int n, molds_blas_int k,
             bool isMatrixAColumnMajor,
             bool isMatrixATransposed,
             bool isLowerTriangularPartMatrixCUsed,
             double alpha, double const* const* matrixA,
             double beta,  double**             matrixC)const;
private:
   Blas();
   ~Blas();
   static Blas* blas;
};
}
#endif
