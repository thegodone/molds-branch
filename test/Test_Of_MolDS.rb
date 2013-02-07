#!/usr/bin/env ruby
#//************************************************************************//
#// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
#// Copyright (C) 2012-2013 Katsuhiko Nishimra                             // 
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

Dir.chdir(File.dirname(__FILE__))

MolDSBin = "../src/MolDS.out".freeze

module AllInclude
	def include? *arg
		true
	end
end

if ARGV.empty?
	Tests = Object.new.extend(AllInclude)
else
	Tests = ARGV.collect do |s|
		s.sub(/\.(in|dat)$/,'').freeze
	end.freeze
end

class TesterOmp
   @@surfixDat = ".dat"
   @@surfixInp = ".in"
   @@tempFile = "temp.dat"
   @@moldsBin = MolDSBin
   @@command = "command: "
   @@deleteDiff = " | gawk '{if(($2!=\"SCF\")&&($3!=\"iter\")){print $0}}' | gawk '{if(($4!=\"time:\")){print $0}}' | gawk '{if(($3!=\"Elapsed\")){print $0}}' | gawk '{if(($2!=\"Elapsed\")){print $0}}' | gawk '{if(($3!=\"Welcome\")){print $0}}' | gawk '{if(($7!=\"residual\")){print $0}}' | gawk '{if(($3!=\"mode(nmw):\") ){print $0}}' | gawk '{if( !(($3==\"mode(mw):\")&&($4<6)) ){print $0}}'" 
	 @@printed_section = []
   def doesTestOmp(mklNumThreads, ompNumThreads)
      return unless should_run?
      ENV["MKL_NUM_THREADS"] = mklNumThreads
      ENV["OMP_NUM_THREADS"] = ompNumThreads
			puts <<EOS % [ENV["MKL_NUM_THREADS"],ENV["OMP_NUM_THREADS"]]
MPI:no
MKL_NUM_THREADS:%s
OMP_NUM_THREADS:%s
EOS
      puts @@command + @moldsCommand
      system(@moldsCommand)
      puts @@command + @diffCommand + @@deleteDiff
      system(@diffCommand + @@deleteDiff)
			puts '','',''
   end
   #def initialize(prefix, section=nil, title)
   #Old ruby workaround.
   #Old versioned ruby accept default values for only last arguments.
   def initialize(prefix, section, title=nil)
      #So swap arguments if section is ommitted.
      if title.nil?
         section,title = title,section
      end
      @prefix = prefix
      @inputFile = prefix + @@surfixInp
      @outputFile = prefix + @@surfixDat
      @moldsCommand = @@moldsBin + " < " + @inputFile + " > " + @@tempFile 
      @diffCommand = "diff " + @outputFile + " " + @@tempFile
			@title = title
			# Update section title if given, otherwise reuse previous one.
			@@section = section unless section.nil?
			print_title
   end

	 private
	 def should_run?
		 @should ||= Tests.include?(@prefix)
	 end

	 def print_title
		 return unless should_run?
		 unless @@printed_section.include?(@@section)
			 @@printed_section << @@section
			 puts @@section,''
		 end
		 puts @title, ''
	 end
end

puts <<EOS

*****************************************
***                                   ***
***                                   ***
***       Start Test for MolDS        ***
***                                   ***
***                  Powered by Ruby  ***
*****************************************
EOS

puts 'MD5 sum of the MolDS.out to be tested:'
system "md5sum #{MolDSBin}"
puts '',''

prefix = "FNC1_principal"
testerOmp = TesterOmp.new(prefix, <<"SECTION",<<"TITLE")
---------------------------------------------------
-----------  Test of principal axes  --------------
---------------------------------------------------
SECTION
\t\t\t>>> F8BT <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "FNC1_rot120"
testerOmp = TesterOmp.new(prefix, <<"SECTION",<<"TITLE")
---------------------------------------------------
-----------  Test of rotate  ----------------------
---------------------------------------------------
SECTION
\t\t\t>>> F8BT <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "FNC1_translate"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
-----------  Test of translate  -------------------
---------------------------------------------------
SECTION
\t\t\t>>> F8BT <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_cndo2"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
-------------------------------------------
----------   Test of CNDO2/HF     ---------
-------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_cndo2"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "h2s_cndo2"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> H2S <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_indo"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
-------------------------------------------
----------   Test of INDO/HF    -----------
-------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_indo"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_zindos_directCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of ZINDO/CIS-singlet     ---------
----------  Without Davidson for the CIS  ---------
---------------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_zindos_directCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "h2s_zindos_directCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> H2S <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "h2o_zindos_directCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> H2O <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_zindos_davidsonCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
------------------------------------------------
----------  Test of ZINDO/CIS-singlet  ---------
----------  With Davidson for the CIS  ---------
------------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_zindos_davidsonCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "h2s_zindos_davidsonCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> H2S <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_zindos_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
------------------------------------------------
----------  Test of ZINDO/HF-Force  ------------
------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_mndo"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
-------------------------------------------
----------   Test of MNDO/HF     ----------
-------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_mndo"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_mndo_directCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of MNDO/CIS-singlet      ---------
----------  Without Davidson for the CIS  ---------
---------------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_mndo_directCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_mndo_davidsonCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of MNDO/CIS-singlet      ---------
----------  With Davidson for the CIS     ---------
---------------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_mndo_davidsonCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_mndo_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
------------------------------------------------
----------  Test of MNDO/HF-Force  ------------
------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_mndo_directCIS_singlet_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of MNDO/CIS-singlet-force --------
----------  Without Davidson for the CIS   --------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_mndo_davidsonCIS_singlet_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of MNDO/CIS-singlet-force --------
----------  With Davidson for the CIS      --------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_am1"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
-------------------------------------------
----------   Test of AM1/HF    ------------
-------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_am1"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_am1_directCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of AM1/CIS-singlet       ---------
----------  Without Davidson for the CIS  ---------
---------------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_am1_directCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<\n
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_am1_davidsonCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of AM1/CIS-singlet       ---------
----------  With Davidson for the CIS     ---------
---------------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_am1_davidsonCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_am1_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
------------------------------------------------
----------  Test of AM1/HF-Force  ------------
------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_am1_directCIS_singlet_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of AM1/CIS-singlet-force  --------
----------  Without Davidson for the CIS   --------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_am1_davidsonCIS_singlet_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of AM1/CIS-singlet-force  --------
----------  With Davidson for the CIS      --------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_pm3"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
-------------------------------------------
----------   Test of PM3/HF    ------------
-------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_pm3_directCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of PM3/CIS-singlet       ---------
----------  Without Davidson for the CIS  ---------
---------------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3_directCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "ch4_pm3_davidsonCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of PM3/CIS-singlet       ---------
----------  With Davidson for the CIS     ---------
---------------------------------------------------
SECTION
\t\t\t>>> CH4 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3_davidsonCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"TITLE")
\t\t\t>>> C2H6 <<<\n
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
------------------------------------------------
----------  Test of PM3/HF-Force  --------------
------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3_directCIS_singlet_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of PM3/CIS-singlet-force  --------
----------  Without Davidson for the CIS   --------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3_davidsonCIS_singlet_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of PM3/CIS-singlet-force  --------
----------  With Davidson for the CIS      --------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3pddg"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
-------------------------------------------
---------- Test of PM3/PDDG/HF ------------
-------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3pddg_directCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of PM3/PDDG/CIS-singlet  ---------
----------  Without Davidson for the CIS  ---------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3pddg_davidsonCIS_singlet"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----------  Test of PM3/PDDG/CIS-singlet  ---------
----------  With Davidson for the CIS     ---------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3pddg_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
------------------------------------------------
----------  Test of PM3/PDDG/HF-Force  ---------
------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3pddg_directCIS_singlet_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
---------  Test of PM3/PDDG/CIS-singlet-force  ----
---------  Without Davidson for the CIS    --------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3pddg_davidsonCIS_singlet_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
---------  Test of PM3/PDDG/CIS-singlet-force  ----
---------  With Davidson for the CIS      --------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3pddg_opt_steepest"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
------------------------------------------------
------  Test of PM3/PDDG/Steepest Descent ------
------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<\n
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3pddg_opt_conjugate"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
------------------------------------------------
----  Test of PM3/PDDG/Conjugate gradient ------
------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3pddg_opt_bfgs"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
------------------------------------------------
------  Test of PM3/PDDG/BFGS ------------------
------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3_MC"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
---------  Test of PM3/HF-MC  ---------------------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3_directCIS_singlet_MC"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
---------  Test of PM3/CIS-singlet-MC    ----------
---------  Without Davidson for the CIS  ----------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3pddg_rpmd"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
------------------------------------------------
----------  Test of PM3/PDDG/RPMD  -------------
------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3pddg_davidsonCIS_singlet_rpmd"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
---------  Test of PM3/PDDG/CIS/RPMD      ---------
---------  With Davidson for the CIS      ---------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3d"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
-------------------------------------------
------------ Test of PM3-D/HF -------------
-------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE

mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3_vdw"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
------  Test of vdw correction in PM3/HF  ---------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3_vdw_force"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
----  Test of vdw correction in PM3/HF-Force  -----
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_pm3_vdw_MC"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
------  Test of vdw correction in PM3/HF-MC  ------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6-h2o-cluster_pm3pddg_freq"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
-------------------------------------------
----   Test of PM3/PDDG/HF/FREQUENCIES  ---
-------------------------------------------
SECTION
\t\t\t>>> C2H6 H2O cluster <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6-nh3-cluster_pm3d_freq"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
-------------------------------------------
-----   Test of PM3-D/HF/FREQUENCIES   ----
-------------------------------------------
SECTION
\t\t\t>>> C2H6 NH3 cluster <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)

prefix = "c2h6_mndo_directCIS_singlet_force_heap_limit"
testerOmp = TesterOmp.new(prefix, <<"SECTION", <<"TITLE")
---------------------------------------------------
-----------  Test of limitation of Heap  ----------
---------------------------------------------------
SECTION
\t\t\t>>> C2H6 <<<
TITLE
mklNumThreads = "1"
ompNumThreads = "1"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)
mklNumThreads = "2"
ompNumThreads = "2"
testerOmp.doesTestOmp(mklNumThreads,ompNumThreads)




system("rm -rf temp.dat")
