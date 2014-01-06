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
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<vector>
#include<stdexcept>
#include<boost/shared_ptr.hpp>
#include<boost/format.hpp>
#include"base/Enums.h"
#include"base/Uncopyable.h"
#include"base/PrintController.h"
#include"base/MolDSException.h"
#include"base/MallocerFreer.h"
#include"mpi/MpiInt.h"
#include"mpi/MpiProcess.h"
#include"base/EularAngle.h"
#include"base/RealSphericalHarmonicsIndex.h"
#include"base/atoms/Atom.h"
#include"base/Molecule.h"
#include"base/MolDS.h"
using namespace std;
using namespace MolDS_base;
int main(int argc, char *argv[]){
   string optionHelp="-h";
   string optionVersion="-v";
   string messageHelp="See README.txt: \"http://sourceforge.jp/projects/molds/scm/svn/tree/head/trunk/doc/\"\n";
   string messageVersion="MolDS 0.4.0(Under development)\n";
   for(int i=0; i<argc; i++){
      if(optionHelp.compare(argv[i])==0){
         std::cout << messageHelp;
         return 0;
      }
      if(optionVersion.compare(argv[i])==0){
         std::cout << messageVersion;
         return 0;
      }
   }
   try{
      MolDS_mpi::MpiProcess::CreateInstance(argc, argv);
      boost::shared_ptr<MolDS_base::MolDS> molds(new MolDS_base::MolDS());
      molds->Run(argc, argv);
      MolDS_mpi::MpiProcess::DeleteInstance();
   }
   catch(exception& ex){
      cout << ex.what();
   }
   return 0;
}

