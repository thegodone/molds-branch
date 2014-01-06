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
#include<sstream>
#include<string>
#include<math.h>
#include<time.h>
#include<omp.h>
#include<boost/format.hpp>
#include"Uncopyable.h"
#include"PrintController.h"
#include"MolDSException.h"
#include"MallocerFreer.h"
#include"../mpi/MpiInt.h"
#include"../mpi/MpiProcess.h"
#include"Utilities.h"
using namespace std;

namespace MolDS_base{
// output welcome message
string Utilities::GetWelcomeMessage(){
   return "\n\n     >>>>>  Welcome to the MolDS world at " + Utilities::GetDateString() + "  <<<<<\n\n\n";
}

// output farewell message
string Utilities::GetFarewellMessage(time_t startTime, clock_t startTick, double ompStartTime, bool runingNormally){
   time_t endTime;
   time(&endTime);
   clock_t endTick = clock();
   double consumedTime = static_cast<double>(endTick - startTick)/static_cast<double>(CLOCKS_PER_SEC);
   double ompEndTime = omp_get_wtime();
   stringstream ss;
   if(runingNormally){
      ss << "\n\n     >>>>>  The MolDS finished normally!  <<<<<\n";
   }
   else{
      ss << "\n\n     >>>>>  The MolDS finished abnormally..............  <<<<<\n";
   }
   ss << "     >>>>>  CPU time: " << consumedTime << "[s].  <<<<<\n";
   ss << "     >>>>>  Elapsed time: " << endTime - startTime << "[s].  <<<<<\n";
   ss << "     >>>>>  Elapsed time(OMP): " << ompEndTime - ompStartTime << "[s].  <<<<<\n";
   ss << "     >>>>>  See you.  <<<<<\n\n\n";
   return ss.str();
}

// string of today
string Utilities::GetDateString(){
   time_t current;
   struct tm *local;
   char  wday_name[][10] = {"Sun.", "Mon.", "Thu.", "Wed.", "Thu.", "Fri.", "Sat."};
   time(&current);
   local = localtime(&current);
   stringstream ss;
   ss << local->tm_year + 1900 
      << "/" 
      << local->tm_mon + 1 
      << "/" 
      << local->tm_mday 
      << "(" 
      << wday_name[local->tm_wday] 
      << ") ";
   ss << local->tm_hour << ":" << local->tm_min << ":" << local->tm_sec;
   return ss.str();
}

// trim the string
string Utilities::TrimString(const string str){
   int nStart = 0;
   int nEnd = str.length() - 1;
   // left trim 
   for(int n = 0; n < str.length(); n++ ){
      if( str.data()[n] != ' ' ){
         nStart = n;
         break;
      }
   }
   // right trim 
   for(int n = str.length() - 1; n >= 0; n-- ){
      if( str.data()[n] != ' ' ){
         nEnd = n;
         break;
      }
   }
   return(str.substr( nStart, nEnd - nStart + 1 ));
}

string Utilities::Num2String(int number, int digit){
   stringstream ss;
   int numberDigit = static_cast<int>(log10(static_cast<double>(number))) + 1;
   for(int i=0; i<digit-numberDigit; i++){
      ss << "0";
   }
   ss << number;
   return ss.str();
}
}

