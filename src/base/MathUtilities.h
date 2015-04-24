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
#ifndef INCLUDED_MATHUTILITIES
#define INCLUDED_MATHUTILITIES
namespace MolDS_base{
// n!
int Factorial(int n);
// nCk
int Conbination(int n, int k);
// rotating matrix
void CalcRotatingMatrix(double matrix[][3], double theta, CartesianType cartesianType);
// calculate determinant of the matrix. Note taht the matrix will be destroid
double GetDeterminant(double** matrix, int dim);
// Normalization of vector
void NormalizeVector(double* vector, int dim);
// Orthonormalization by Schmidt method
void OrthoNormalizeSchmidt(double** vectors, int dim, int numVector);
}
#endif
 
