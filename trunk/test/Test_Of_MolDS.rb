#!/usr/bin/env ruby
#//************************************************************************//
#// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
#// Copyright (C) 2012-2012 Katsuhiko Nishimra                             // 
#//                                                                        // 
#// This file is part of MolDS.                                            // 
#//                                                                        // 
#// MolDS is free software: you can redistribute it and/or modify          // 
#// it under the terms of the GNU General Public License as published by   // 
#// the Free Software Foundation, either version 3 of the License, or      // 
#// (at your option) any later version.                                    // 
#//                                                                        // 
#// MolDS is distributed in the hope that it will be useful,               // 
#// but WITHOUT ANY WARRANTY; without even the implied warranty of         // 
#// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          // 
#// GNU General Public License for more details.                           // 
#//                                                                        // 
#// You should have received a copy of the GNU General Public License      // 
#// along with MolDS.  If not, see <http://www.gnu.org/licenses/>.         // 
#//************************************************************************//
class TesterOmp
   @@surfixDat = ".dat"
   @@surfixInp = ".in"
   @@tempFile = "temp.dat"
   @@moldsBin = "../src/MolDS.out"
   @@command = "command: "
   @@deleteDiff = " | gawk '{if(($4!=\"RMS\")){print $0}}' | gawk '{if(($4!=\"time:\")){print $0}}' | gawk '{if(($3!=\"Elapsed\")){print $0}}' | gawk '{if(($2!=\"Elapsed\")){print $0}}' | gawk '{if(($3!=\"Welcome\")){print $0}}' | gawk '{if(($7!=\"residual\")){print $0}}'"
   def doesTestOmp(prefix, mklNumThreads, ompNumThreads)
      setPrefix(prefix)
      ENV["MKL_NUM_THREADS"] = mklNumThreads
      ENV["OMP_NUM_THREADS"] = ompNumThreads
      system("echo MPI:no")
      system("echo MKL_NUM_THREADS:$MKL_NUM_THREADS")
      system("echo OMP_NUM_THREADS:$OMP_NUM_THREADS")
      print @@command + @moldsCommand + "\n"
      system(@moldsCommand)
      print @@command + @diffCommand + @@deleteDiff + "\n"
      system(@diffCommand + @@deleteDiff)
      system("echo '\n\n'")
   end
   def setPrefix(prefix)
      @inputFile = prefix + @@surfixInp
      @outputFile = prefix + @@surfixDat
      @moldsCommand = @@moldsBin + " < " + @inputFile + " > " + @@tempFile 
      @diffCommand = "diff " + @outputFile + " " + @@tempFile
   end
   private :setPrefix
end

system("echo ")
system("echo '*****************************************'")
system("echo '***                                   ***'")
system("echo '***                                   ***'")
system("echo '***       Start Test for MolDS        ***'")
system("echo '***                                   ***'")
system("echo '***                    Power by Ruby  ***'")
system("echo '*****************************************\n\n'")

testerOmp = TesterOmp.new

system("echo '---------------------------------------------------'")
system("echo '-----------  Test of principal axes  --------------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> F8BT <<<\n'")
prefix = "FNC1_principal"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '-----------  Test of rotate  ----------------------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> F8BT <<<\n'")
prefix = "FNC1_rot120"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '-----------  Test of translate  -------------------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> F8BT <<<\n'")
prefix = "FNC1_translate"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '-------------------------------------------'")
system("echo '----------   Test of CNDO2/HF     ---------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_cndo2"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_cndo2"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> H2S <<<\n'")
prefix = "h2s_cndo2"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '-------------------------------------------'")
system("echo '----------   Test of INDO/HF    -----------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_indo"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_indo"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of ZINDO/CIS-singlet     ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_zindos_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_zindos_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> H2S <<<\n'")
prefix = "h2s_zindos_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> H2O <<<\n'")
prefix = "h2o_zindos_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of ZINDO/CIS-singlet  ---------'")
system("echo '----------  With Davidson for the CIS  ---------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_zindos_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_zindos_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> H2S <<<\n'")
prefix = "h2s_zindos_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of ZINDO/HF-Force  ------------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_zindos_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '-------------------------------------------'")
system("echo '----------   Test of MNDO/HF     ----------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_mndo"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of MNDO/CIS-singlet      ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_mndo_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of MNDO/CIS-singlet      ---------'")
system("echo '----------  With Davidson for the CIS     ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_mndo_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of MNDO/HF-Force  ------------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of MNDO/CIS-singlet-force --------'")
system("echo '----------  Without Davidson for the CIS   --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo_directCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of MNDO/CIS-singlet-force --------'")
system("echo '----------  With Davidson for the CIS      --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo_davidsonCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '-------------------------------------------'")
system("echo '----------   Test of AM1/HF    ------------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_am1"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_am1"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of AM1/CIS-singlet       ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_am1_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_am1_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of AM1/CIS-singlet       ---------'")
system("echo '----------  With Davidson for the CIS     ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_am1_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_am1_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of AM1/HF-Force  ------------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_am1_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of AM1/CIS-singlet-force  --------'")
system("echo '----------  Without Davidson for the CIS   --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_am1_directCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of AM1/CIS-singlet-force  --------'")
system("echo '----------  With Davidson for the CIS      --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_am1_davidsonCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '-------------------------------------------'")
system("echo '----------   Test of PM3/HF    ------------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_pm3"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of PM3/CIS-singlet       ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_pm3_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of PM3/CIS-singlet       ---------'")
system("echo '----------  With Davidson for the CIS     ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> CH4 <<<\n'")
prefix = "ch4_pm3_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of PM3/HF-Force  --------------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of PM3/CIS-singlet-force  --------'")
system("echo '----------  Without Davidson for the CIS   --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_directCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of PM3/CIS-singlet-force  --------'")
system("echo '----------  With Davidson for the CIS      --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_davidsonCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '-------------------------------------------'")
system("echo '---------- Test of PM3/PDDG/HF ------------'")
system("echo '-------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of PM3/PDDG/CIS-singlet  ---------'")
system("echo '----------  Without Davidson for the CIS  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_directCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '----------  Test of PM3/PDDG/CIS-singlet  ---------'")
system("echo '----------  With Davidson for the CIS     ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_davidsonCIS_singlet"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of PM3/PDDG/HF-Force  ---------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '---------  Test of PM3/PDDG/CIS-singlet-force  ----'")
system("echo '---------  Without Davidson for the CIS    --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_directCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '---------  Test of PM3/PDDG/CIS-singlet-force  ----'")
system("echo '---------  With Davidson for the CIS      --------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_davidsonCIS_singlet_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '------  Test of PM3/PDDG/Steepest Descent ------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_opt_steepest"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '---------  Test of PM3/HF-MC  ---------------------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_MC"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '---------  Test of PM3/CIS-singlet-MC    ----------'")
system("echo '---------  Without Davidson for the CIS  ----------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_directCIS_singlet_MC"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '------------------------------------------------'")
system("echo '----------  Test of PM3/PDDG/RPMD  -------------'")
system("echo '------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_rpmd"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '---------  Test of PM3/PDDG/CIS/RPMD      ---------'")
system("echo '---------  With Davidson for the CIS      ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3pddg_davidsonCIS_singlet_rpmd"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '------  Test of vdw correction in PM3/HF  ---------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_vdw"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
system("echo '---------------------------------------------------'")
system("echo '----  Test of vdw correction in PM3/HF-Force  -----'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_vdw_force"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
system("echo '---------------------------------------------------'")
system("echo '------  Test of vdw correction in PM3/HF-MC  ------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_pm3_vdw_MC"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)

system("echo '---------------------------------------------------'")
system("echo '-----------  Test of limitation of Heap  ----------'")
system("echo '---------------------------------------------------\n'")
system("echo '\t\t\t>>> C2H6 <<<\n'")
prefix = "c2h6_mndo_directCIS_singlet_force_heap_limit"
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(prefix, mklNumThreads,ompNumThreads)




system("rm -rf temp.dat")
