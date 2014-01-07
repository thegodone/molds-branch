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

compile for 64bit OS with intel mkl:
mpicxx ../../base/Enums.cpp ../../base/PrintController.cpp ../../base/MolDSException.cpp  ../../base/MathUtilities.cpp ../../base/MallocerFreer.cpp ../../mpi/MpiProcess.cpp ../../wrappers/Blas.cpp ../../wrappers/Lapack.cpp ../../base/EularAngle.cpp ../../base/Parameters.cpp  paramYZ.cpp -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -O0 -openmp -DMKL_INT=intptr_t -DMKL_ILP64 -I/usr/local/boost/include -L/usr/local/boost/lib -Wl,-rpath=/usr/local/boost/lib -lboost_serialization -lboost_mpi -lboost_thread -lboost_serialization >& err.dat 
