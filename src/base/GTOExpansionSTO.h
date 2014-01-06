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
#ifndef INCLUDED_GTOEXPANSIONSTO
#define INCLUDED_GTOEXPANSIONSTO
namespace MolDS_base{

// GTOExpansionSTO is singleton
class GTOExpansionSTO: private Uncopyable{
public:
   static GTOExpansionSTO* GetInstance();
   static void DeleteInstance();
   double GetExponent(MolDS_base::STOnGType stonG, 
                      MolDS_base::ShellType shellType, 
                      OrbitalType orbitalType, 
                      int index) const;
   double GetCoefficient(MolDS_base::STOnGType stonG, 
                         MolDS_base::ShellType shellType, 
                         OrbitalType orbitalType, 
                         int index) const;

private:
   static GTOExpansionSTO* gTOExpansionSTO;
   GTOExpansionSTO();
   ~GTOExpansionSTO();

   std::string errorMessageGetCoefficientNonValidOrbital;
   std::string errorMessageGetExponentNonValidOrbital;
   std::string errorMessageOrbitalType;
   std::string errorMessageSTOnGType;
   void SetCoefficientsExponents();
   double exponents[MolDS_base::STOnGType_end]
                   [MolDS_base::ShellType_end]
                   [MolDS_base::AzimuthalType_end]
                   [6];    
      //[N:expansion number][Shelltype][quasi orbital type:s, p, or d][expansion index]. 
      //This is alpha in (3) of [S_1970]. See Table I and II in [S_1970]
   double coefficients[MolDS_base::STOnGType_end]
                      [MolDS_base::ShellType_end]
                      [MolDS_base::AzimuthalType_end]
                      [6]; 
      //[N:expansion number][Shelltype][quasi orbital type:s, p, or d][expansion index]. 
      //This is d in (3) of [S_1970]. See Table I and II in [S_1970]
};

}
#endif





