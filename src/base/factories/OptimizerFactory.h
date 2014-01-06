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
#ifndef INCLUDED_OPTIMIZER_FACTORY
#define INCLUDED_OPTIMIZER_FACTORY
namespace MolDS_base_factories{

class OptimizerFactory{
public:
   static MolDS_optimization::Optimizer* Create(MolDS_base::OptimizationMethodType methodType);
   static MolDS_optimization::Optimizer* Create();
private:
   OptimizerFactory();
   ~OptimizerFactory();
   static std::string errorMessageNotEnableOptimizationMethod;
   static std::string errorMessageOptimizationMethodType;
};

}
#endif





