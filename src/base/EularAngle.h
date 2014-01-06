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
#ifndef INCLUDED_EULARANGLE
#define INCLUDED_EULARANGLE
namespace MolDS_base{

class EularAngle{
public:
   EularAngle();
   EularAngle(double x, double y, double z);
   explicit EularAngle(double* angles);
   double GetAlpha() const{return this->alpha;}
   double GetBeta()  const{return this->beta;}
   double GetGamma() const{return this->gamma;}
   void SetAlpha(double alpha){this->alpha = alpha;}
   void SetBeta (double beta) {this->beta  = beta;}
   void SetGamma(double gamma){this->gamma = gamma;}
private:
   std::string errorMessageInvalidXYZ;
   double alpha;
   double beta;
   double gamma;
   void SetMessage();
};

}
#endif
