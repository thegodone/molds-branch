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

The files in this directory were prepared for ghost atoms used in basis 
set super position error (BSSE). They can be specified in input files as:
         
   <GHOST>
      To set geometry of ghost atoms which provide only basis funcions,  
      write ghost-directive. The grammar to write the ghost-directive is 
      same with the grammar of the geometry-directive. The each ghost at-
      om has 0-mass, 0-core-charge, 0-electrons, 0-vdW-coefficient, and 
      0-pairwise-parameters in PDDG.

       E.g.
         GHOST
            O -0.3810  1.1411 5.0000
            H -0.2681 -0.5205 5.9016
            H -0.3681 -0.4725 5.8016
         GHOST_END
        
HOWEVER, treatment of ghost atoms with semiempirical methods has not been 
unfortunately established yet. In paticular, treatment of the core integ-
rals, U_{\mu\mu}, U_{\mu\nu}, and beta_{ij}, has not been established yet. 
As far as M. F. knows, threre are a trial for the BSSE with CNDO/INDO:
"Ghost orbitals in semiempirical methods. Estimation of basis set superp-
osition error", Józef Lipiński, Henryk Chojnacki, Int. J. of Quan. Chem.,
Vol. 19, 891(1981). M.F. could not find any paper for the BSSE with MNDOs.
Therefore, M. F. decided that these ghost atoms were not used anywhere a-
lthough these files for the ghost atoms are included in the source tree.

