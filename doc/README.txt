//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
// Copyright (C) 2012-2012 Katsuhiko Nishimra                             // 
// Copyright (C) 2013-2013 Michihiro Okuyama                              // 
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

==============================================================================


   MolDS ("Mol"ecular "D"ynamics "S"imulation package) ver. 0.2.0
      developed by Mikiya Fujii, Ph.D., Katsuhiko Nishimra, and Michihiro Okuyama, Ph.D..
      For Questions and bug reports: molds-dev@lists.sourceforge.jp


==============================================================================
REQUIREMENTS:
   MolDS requires c/c++ compiler of Intel (icpc) or GNU (g++) and boost-libraries. 
   Valid versions of these compiler are icpc 12.0.4(MkL 10.3 update 4), g++ 4.4, or later 
   because the MolDS is implemented with openMP 3.0. 
   To compile MolDS with g++, furthermore, openBLAS (version 0.2.5 or later) is also required. 
   The default compiler is the intel c++ compiler (icpc). 

   To get and install the boost-libraries, see the HP:<http://www.boost.org/>.
   The version of the boost would be no problem if 1.48.0 or later is used.

   To get and install the openBLAS-libraries, see the HP:<http://xianyi.github.com/OpenBLAS/>.
   Note that "USE_OPENMP = 1" should be set for the installation of the opneBLAS.
   Furthermore, "INTERFACE64 = 1" is also needed when you install the openBLAS into 64-bits machines

==============================================================================
COMPILE(using GNUmake): 
   In the "src" directory in the MolDS package.

   Case i) Using Intel c/c++ compiler (icpc)
      Change the "BOOST_TOP_DIR" in Makefile to the top directory of the 
      boost-libraries in your systems.

      To compile MolDS on 32 bits machine,
      $ make INTEL=32

      To compile MolDS on 64 bits machine,
      $ make INTEL=64

   Case ii) Using GNU c/c++ compiler (g++)
      Change the "BOOST_TOP_DIR" in "Makefile_GNU" to the top directory of the 
      boost-libraries in your systems.
      Change the "OPENBLAS_TOP_DIR" in "Makefile_GNU" to the top directory of the 
      boost-libraries in your systems.
      
      Then, just type: 
      $ make -f Makefile_GNU 

   For both case, the compile succeeded if you could fine "MolDS.out" in the "src" directory. 
   Type "$ make clean" when you wanna clean the compilation.
   If you want to compile MolDS in debug-mode, 
   -g, -rdynamic(for function names in backtrace) and -DMOLDS_DBG should be added to CFLAGS,
   that is, hit the following command:
   $make CFLAGS="-O0 -g -rdynamic -DMOLDS_DBG"

==============================================================================
CARRY OUT MolDS:
   After the compile, in the "src" directory,
   $ ./MolDS.out < input.in
   or
   $ ./MolDS.out input.in

==============================================================================
SAMPLE and TEST
   See files in "test" directories for sample files.
   In the "test" directory, *.in files are input files, then *.dat files are
   associated output files. To execute all test cases, carry out below ruby-script
   in the "test" directory. This script will finished in a few minutes with big
   output(a few thausands lines).

   $ ruby Test_Of_MolDS.rb

   To execute some specific test cases, carry out below ruby-script with the test names you want. 
   When the test names are specifed, ".in" and ".dat" in the arguments will be ignored.

   $ ruby Test_Of_MolDS.rb test1.in test2.dat test3 ...

==============================================================================
CAPABILITIES:

   Electronic state and molecular dynamics:
            | HF  | CIS | MD(gs) | MD(es) | MC(gs) | MC(es) | RPMD(gs) | RPMD(es) | Optimize(gs) | Optimize(es) | Frequencies(gs) |
   ---------|-----|-----|--------|--------|--------|--------|----------|----------|--------------|--------------|-----------------|
   CNDO2    | OK  | --  | --     | --     | OK     | --     | --       | --       | --           | --           | --              |
   ---------|-----|-----|--------|--------|--------|--------|----------|----------|--------------|--------------|-----------------|
   INDO     | OK  | --  | --     | --     | OK     | --     | --       | --       | --           | --           | --              |
   ---------|-----|-----|--------|--------|--------|--------|----------|----------|--------------|--------------|-----------------|
   ZINDO/S  | OK  | OK  | OK     | OK     | OK     | OK     | OK       | OK       | OK           | OK           | --              |
   ---------|-----|-----|--------|--------|--------|--------|----------|----------|--------------|--------------|-----------------|
   MNDO     | OK  | OK  | OK     | OK     | OK     | OK     | OK       | OK       | OK           | OK           | OK              |
   ---------|-----|-----|--------|--------|--------|--------|----------|----------|--------------|--------------|-----------------|
   AM1      | OK  | OK  | OK     | OK     | OK     | OK     | OK       | OK       | OK           | OK           | OK              |
   ---------|-----|-----|--------|--------|--------|--------|----------|----------|--------------|--------------|-----------------|
   AM1-D    | OK  | OK  | OK     | OK     | OK     | OK     | OK       | OK       | OK           | OK           | OK              |
   ---------|-----|-----|--------|--------|--------|--------|----------|----------|--------------|--------------|-----------------|
   PM3      | OK  | OK  | OK     | OK     | OK     | OK     | OK       | OK       | OK           | OK           | OK              |
   ---------|-----|-----|--------|--------|--------|--------|----------|----------|--------------|--------------|-----------------|
   PM3-D    | OK  | OK  | OK     | OK     | OK     | OK     | OK       | OK       | OK           | OK           | OK              |
   ---------|-----|-----|--------|--------|--------|--------|----------|----------|--------------|--------------|-----------------|
   PM3/PDDG | OK  | OK  | OK     | OK     | OK     | OK     | OK       | OK       | OK           | OK           | OK              |

      "OK", "Sch", and "--" mean available, shceduled, and non-scheduled methods, respectively.
      "gs" and "es" mean ground and excited states, respectively.
      i.e., MD(gs) and MD(es) mean Born-Oppenheimer Molecular Dynamics on ground and excited states, respectively. 

   Elements:
   CNDO2    | H, Li, C, N, O, and S
   INDO     | H, Li, C, N, and O
   ZINDO/S  | H, C, N, O, and S
   MNDO     | H, C, N, O, and S 
   AM1      | H, C, N, O, and S 
   AM1-D    | H, C, N, O, and S 
   PM3      | H, C, N, O, and S 
   PM3-D    | H, C, N, O, and S 
   PM3/PDDG | H, C, N, O, and S 

==============================================================================
HOW TO WRITE INPUT:

   <Comment Out>
      Lines starting with "//" or "#" in input-files are treated as comments.

   <SCF>
      Write "cndo/2", "indo", "zindo/s", "mndo", "am1", "am1-d",
      "pm3", "pm3-d", or "pm3/pddg" in theory-directive.
      This theory-directive indicate a electronic structure theory used in your simulations.
      MNDO only supports (can calculate) Heats of formation.
      SCF module outputs also the dipole moment arrond the center of cores of the molecule.
      To calculate the dipole moment, STO-6G [DY_1977] is used.

      E.g. 
         THEORY
            indo 
         THEORY_END
   
      -options
       "max_iter", "rms_density", "damping_thresh", "damping_weight", 
       "diis_num_error_vect", "diis_start_error", "diis_end_error",
       "vdW", "vdW_s6", and "vdW_d" are prepared as options.

       The default value of the "max_iter" is 100.
       The default value of the "rms_density" is 10**(-8.0).
       The default value of the "damping_thresh" is 1.
       The default value of the "damping_weight" is 0.8.
       The default value of the "diis_num_error_vect" is 5.
       The default value of the "diis_start_error" is 0.01.
       The default value of the "diis_end_error" is 10**(-8.0).

       "vdW" should be set as "yes" or "no". 
       When "yes" is set, Grimmes's empirical van der Waals correction([G_2004]) is applied.
       Note that this empirical van der Waals correction is applied to the semiempirical theories 
       of which semiempirical parameters are not modified.
       If user wants to use PM3-D or AM1-D of which semiempirical parameters are modified to be suite for vdW, 
       set theory-directive as "PM3-D" or "AM1-D". 
       When PM3-D or AM1-D is used, users do not need to set "vdW", "vdW_s6", and "vdW_d".
       Generally, PM3-D and AM1-D are recommended for noncovalent complexes 
       than naitive PM3 and AM1 with this empirical vdW, respectively.
       The default value of the "vdW" with the theories except for PM3-D and AM1-D is "no". 
       For PM3-D and AM1-D, this "vdW" is always "yes" whethere user sets or not.

       "vdW_s6" is a scaling factor in the Grimme's van der Waals correction([G_2004]).
       The default value of the "vdW_s6" is 1.4. 
       For PM3-D and AM1-D, this "vdW_s6" is forced to be set as 1.4.

       "vdW_d" is a damping factor in the Grimme's van der Waals correction([G_2004]).
       The default value of the "vdW_d" is 23.0.
       For PM3-D and AM1-D, this "vdW_s6" is forced to be set as 23.0.

       E.g.
         SCF
            max_iter 200
            rms_density 0.00000001
            damping_thresh 0.1
            damping_weight 0.7
            diis_num_error_vect 6
            diis_start_error 0.01
            diis_end_error 0.00000001
            vdW yes
            vdW_s6 0.75
            vdW_d 30
         SCF_END
   
   <MEMORY>
      For settings of memory usage, write options in memory-directive.

      E.g.
         MEMORY
            (options)
         MEMORY_END

      -options
       "limit_heap" is only prepared. Note that this limitation is not 
       the exact limitation of heap usage. Please consider this option as a rough limitation.
       The value of this option should be written with the MByte unit.
       The default value is 256[MB].

      E.g.
         MEMORY
            limit_heap 512
         MEMORY_END

   <Frequencies (Normal modes analysis)>
      write frequencies-directive. Note taht not only the frequencies but also the normal modes are calculated.

      E.g.
         FREQUENCIES
            (options)
         FREQUENCIES_END

      -options
       "electronic_state" is only prepared.
       "electronic_state" is index of the electronic state used for calculating the normal modes. 
       electronic_state=0 means the electronic ground state.
       electronic_state=1 means, then, first electornic excited state.
       The default value of the "electronic_state" is 0.

       E.g. 
         FREQUENCIES
            electronic_state 0
         FREQUENCIES_END

   <MO Plot>
      write MO plot directive.

      E.g.
         MOPLOT 
            (options)
         MOPLOT_END

      -options
       "mo", "grid_number", "frame_length", and "file_prefix" are prepared.

       "mo" is index of the molcular orbital. mo=0 means the lowest energy MO.
       The default value of the "mo" is not set.

       "grid_number" is the grid number of the frame in xyz-coordinates.
       The default values are 25, 25, and 25 for x, y, and z coordinates, respectively.

       "frame_length" is the length of the frame of each coordinate.
       The default values are 10, 10, and 10[angst.] for x, y, and z coordinates.

       "file_prefix" is a prefix of the file name to which the MO is written.
       The default values is "MO_".

      E.g.
         MOPLOT
            mo 5
            mo 8
            grid_number 30 30 30
            frame_length 10 10 10
            file_prefix MOPlot_
         MOPLOT_END

   <CIS>
      Write CIS-directive.

      E.g.
         CIS
            (options)
         CIS_END
   
      -options
       "davidson", "active_occ", "active_vir", "max_iter", "max_dim", "norm_tol", 
       "nstates", "exciton_energies", "all_transition_dipole_moments", 
       "mulliken", "unpaired_electron_population", and "num_print_coefficients" are prepared as options.

       "davidson" should be set as "yes" or "no". 
       The default value of the "davidson" is "yes".

       "active_occ" ("active_vir") is set to the number of occupied (virtual) orbitals
       if user set "active_occ" ("active_vir") to be greater than 
       the number of occupied (virtual) orbitals. 
       The default value of the "active_occ" is 10. The default value of the "active_vir" is 10.

       "nstates" means the number of the target electronic excited stetes.
       "nstates" is valid for the Davidson algorithm only, 
       hence "nstates" is set to "active_occ*active_vir" 
       in direct CIS algorithm (without the Davidson algorithem). 
       The default value of the "nstates" is 5 for the Davidson algorithem.

       "max_iter" is valid for the Davidson algorithm only. 
       This option means the number of times of the maximum Davidson roop. 
       The default value of the "max_iter" is 100.

       "max_dim" is valid for the Davidson algorithm only. 
       This option means the number of slater determinans used by expansion of the excited states. 
       Note that Hartree-Fock state (groudn state) is not included in the "max_dim".
       The default value of the "max_dim" is 100.

       "norm_tol" is valid for the Davidson algorithm only. 
       This option means the max tolerance for the norm of the residual vectors.
       The default value of the "norm_tol" is 10**(-6.0).

       "exciton_energies" should be set as "yes" or "no".
       In the case "yes" is set, free exciton and exciton biding energies are calculated 
       for each excited states.
       The default value of the "exciton_energies" is "no".

       "all_transition_dipole_moments" should be set as "yes" or "no". 
       The default value of the "all_transition_dipole" is "no".
       If user set this "all_transition_dipole" as "yes", all transition dipole moments 
       including between excited states would be calculated. 
       Otherwise "no", transition dipole moments from ground state to each excited state are calculated.

       "mulliken" is a option of mulliken popultaion analysis of the excited state.
       When "mulliken x" is included in CIS-directive, the mulliken popultaion of xth excited state is calculated.
       Multiple indication of these mulliken options is possible. 
       Note that "mulliken 0" is ignored because 0th excited state is the ground state.
       Default setting of this "mulliken" option is nothing.

       "unpaired_electron_population" is a option of unpaired electron population(UEP) analysis of the excited state.
       When "unpaired electron population yes" and "mulliken x" are included in CIS-directive, 
       the UEP of xth excited state is calculated.
       By multiple indication of these mulliken option, the UEP on multiple excited states are possible.
       Note that the UEP on ground state is ignored. 
       Default setting is "unpaired_electron_population" option is nothing.

       "num_print_coefficients" is a number of the coefficients of CIS-eigenvector shown in output.
       The default value of the "num_print_coefficients" is 1.

       E.g.
         CIS
            davidson no
            active_occ 2
            active_vir 2
            nstates 1000
            max_iter 100
            max_dim 100
            norm_tol 0.000001
            mulliken 1 
            mulliken 2
	         unpaired_electron_population yes
         CIS_END

   <Hole Plot>
      Write hole plot directive for the output of the density of the hole.
      Note that this hole plot is valid only when CIS is required.

      E.g.
         HOlEPLOT
            (options)
         HOLEPLOT_END

      -options
       "electronic_state", "grid_number", "frame_length", and "file_prefix" are prepared.

       "electronic_state" is index of the electronic state. 
       electronic_state=0 means the electronic ground state.
       electronic_state=1 means, then, first electornic excited state.
       The default value of the "hole" is not set.

       "grid_number" is the grid number of the frame in xyz-coordinates.
       The default values are 25, 25, and 25 for x, y, and z coordinates, respectively.

       "frame_length" is the length of the frame of each coordinate.
       The default values are 10, 10, and 10[angst.] for x, y, and z coordinates.

       "file_prefix" is a prefix of the file name to which the density of the hole is written.
       The default values is "hole_".

      E.g.
         HOLEPLOT
            electronic_state 5
            electronic_state 8
            grid_number 30 30 30
            frame_length 10 10 10
            file_prefix HOLEPlot_
         HOLEPLOT_END

   <Particle Plot>
      Write particle plot directive for the output of the density of 
      the particle which lives avobe Fermi's sea.
      Note that this plot is valid only when CIS is required.

      E.g.
         PARTICLEPLOT
            (options)
         PARTICLEPLOT_END

      -options
       "electronic_state", "grid_number", "frame_length", and "file_prefix" are prepared.

       "electronic_state" is index of the electronic state. 
       electronic_state=0 means the electronic ground state.
       electronic_state=1 means, then, first electornic excited state.
       The default value of the "hole" is not set.

       "grid_number" is the grid number of the frame in xyz-coordinates.
       The default values are 25, 25, and 25 for x, y, and z coordinates, respectively.

       "frame_length" is the length of the frame of each coordinate.
       The default values are 10, 10, and 10[angst.] for x, y, and z coordinates.

       "file_prefix" is a prefix of the file name to which the density of the particle is written.
       The default values is "particle_".

      E.g.
         HOLEPLOT
            electronic_state 5
            electronic_state 8
            grid_number 30 30 30
            frame_length 10 10 10
            file_prefix HOLEPlot_
         HOLEPLOT_END

   <OPT (geometry optimization)>
      Write OPT-directive. This module uses line search and steepest descent algorythms.
      In the early stage the line search algorythm is used, 
      then the algorythm used in this module is switched to steepest descent algorythm.
      Note that ZINDO/S is not suitable for geometry optimizations.

      E.g.
         OPTIMIZE
            (options)
         OPTIMIZE_END
  
      -options
       "method", "total_steps", "electronic_state", "max_gradient", "rms_gradient", 
        "dt", "initial_trust_radius" and "max_norm_step" are prepared as options.

       "method" should be set as "conjugate_gradient", "steepest_descent", or "bfgs". 
       The default of the "method" is conjugate gradient.

       "electronic_state" means the electronic eigenstate 
       on which the system runs.
       The default value of the "electronic_state" is 0. That is, 
       electronic ground state is default.

      "line_search_times" means the times of line-search trials.
      The default value of the "line_search_times" is 50.
      This parameter have no effect if method is "bfgs".

      "steep_step" means the number of steps of the steepest descent.
      The default value of the "steep_step" is 50.

      "max_gradient" and "rms_gradient" are threshold of the steepest descent.
      The "max_gradient" and "rms_gradient" means maximum and root-mean-squre of the gradient, respectively.
      The default value of the "max_gradient" is 0.00045.
      The default value of the "rms_gradient" is 0.00030.

      "dt" is initial fictious time steps for the steepest descent algorythms.
      The default value of the "dt" is 50[fs].
      This parameter have no effect if method is "bfgs".

      "initial_trust_radius" is an initial value for trust radius used by BFGS method.
      The default value of the "initial_trust_radius" is 0.3.
      This parameter have no effect if method is "steepest_descent" or "conjugate_gradient".

      "max_norm_step" is the maximum value for trust radius used by BFGS method.
      The default value of the "max_norm_step" is 0.3.
      This parameter have no effect if method is "steepest_descent" or "conjugate_gradient".

      E.g.
         OPTIMIZE
            method steepest_descent
            total_steps 50
            electronic_state 0
            max_gradient 0.00045
            rms_gradient 0.00030
            dt 50
         OPTIMIZE_END
  
      

   <MD (Molecular dynamics)>
      Write MD-directive. Note that ZINDO/S is not suitable for molcular dynamics simulations.

      E.g.
         MD 
            (options)
         MD_END
  
      -options
       "total_steps", "electronic_state", 
       and "dt" are prepared as options.

       "electronic_state" means the electronic eigenstate 
       on which the system runs.
       The default value of the "electronic_state" is 0. That is, 
       electronic ground state is default.

       The default value of the "total_steps" is 10. 

       "dt" means the time width of molecular dynamics.
       "dt" should be set in femto-second.
       The default value of the "dt" is 0.1[fs].

      E.g.
         MD
            total_steps 50
            electronic_state 0
            dt 0.05
         MD_END

    <Nonadiabatic semiclassical kernel based on ovelap integrals (NASCO)>
       Write NASCO-directive. This module can calculate only nonadiabatic trajectories,
       that is, semiclassical prefactor can not be calculated. 
       If user did not set CIS-directive, the default settings for CIS are used for NASCO, 
       see also <CIS> sections.
 
       E.g.
          NASCO
             (options)
          NASCO_END
   
       -options
        "total_steps", "num_electronic_states", "initial_electronic_state", "seed", 
        and "dt" are prepared as options.
 
        The default value of the "total_steps" is 10. 
 
        "num_electronic_states" means the number of the electronic 
        eigenstates used in NASCO simulations. 
        "num_electronic_states" minus 1 should be not over "nstates" in CIS-conditons.
        The default value of the "num_electronic_states" is 3, 
        that is, ground, 1st, and 2nd excited states.
 
        "initial_electronic_state" means the electronc eigenstates 
        form which trajectories start to run.
        "initial_electronic_state=0" means that the trajectories run from ground state.
        The "initial_electronic_state should be less than the "num_electronic_states".
        i.e., "initial_electronic_state=3" with "num_electronic_states=3" leads to error.
        The default value of the "initial_electronic_state" is 0.
 
        "seed" means the seed of the random-number-generator.
        The random numbers are used for trajectory-hopping.
        Default seed is generated by the time. 
        When you want to carry out many time jobs with same condition,
        "seed" should be set to an identic positive integer in each job.
 
        "dt" means the time width of molecular dynamics.
        "dt" should be set in femto-second.
        The default value of the "dt" is 0.1[fs].
 
       E.g.
          NASCO
             total_steps 50
             num_electronic_states 10
             seed 398
             dt 0.05
          NASCO_END
 
   <MC (Monte Carlo)>
      Write MC-directive. The canonical sampling is only implemented.

      E.g.
         MC 
            (options)
         MC_END
  
      -options
       "total_steps", "electronic_state", "temperature", "seed"
       and "step_width" are prepared as options.

       The default value of the "total_steps" is 10. 

       "electronic_state" means the electronic eigenstate 
       on which the system walks.
       The default value of the "electronic_state" is 0. That is, 
       electronic ground state is default.

       "temperature" means the temperature in the MC sampling.
       The default value of the "temeprture" is 300[K].

       "seed" means the seed of the random-number-generator.
       The random numbers are used during MC sampling.
       Default seed is generated by the time. 
       When you want to carry out many time jobs with same condition,
       "seed" should be set to an identic positive integer in each job.

       "step_width" means the max absolute displacement (step) width of each Cartesian coordinate.
       Namely, the actual displacement in the MC is in the range [-step_width, step_width).
       "step_width" should be set in angstrom unit.
       The default value of the "step_width" is 0.05[angstrom].

      E.g.
         MC
            total_steps 50
            electronic_state 0
            step_width 0.08
            seed 398
         MC_END

   <RPMD (Ring Polymer Molecular Dynamics)>
      Write RPMD-directive. 

      E.g.
         RPMD 
            (options)
         RPMD_END
  
      -options
       "total_steps", "electronic_state", "num_electonic_states", "temperature", 
       "num_beads", "seed", and "dt" are prepared as options.

       The default value of the "total_steps" is 10. 

       "electronic_state" means the electronic eigenstate 
       on which the system walks.
       The default value of the "electronic_state" is 0. 
       That is, electronic ground state is default.

       "num_electronic_states" means the number of the electronic eigenstate used in RPMD simulation.
       This option is used only for multisurface RPMD only.
       For the single surface RPMD, this "num_electronic_states" is automatically set to 1.
       The default value of the "num_electronic_states" is 1. 

       "temperature" means the temperature in the RPMD.
       The default value of the "temeprture" is 300[K].

       "num_beads" means the number of beads for the ring polymer.
       The default value of the "num_beads" is 10.

       "seed" means the seed of the random-number-generator.
       The random numbers are used during initial condition sampling.
       Default seed is generated by the time. 
       When you want to carry out many time jobs with same condition,
       "seed" should be set to an identic positive integer in each job.

       "dt" means the time width of molecular dynamics.
       "dt" should be set in femto-second.
       The default value of the "dt" is 0.1[fs].

      E.g.
         RPMD 
            total_steps 20
            electronic_state 0
            num_electronic_states 10
            temperature 100
            seed 398
            num_beads 20
            dt 0.5
         RPMD_END

   <Principal Axes (Diagonalizing the inertia tensor)>
      Write inertia-directive.

      E.g.
        INERTIA
           (options)
        INERTIA_END

      -options
       option is "origin" only for setting the origin of the inertia tensor.
       options are written in inertia-directive in angstrom unit.
       Center of mass is used as origin when the "origin" is not set.

      E.g.
        INERTIA
           origin 1.2 2.3 3.4
        INERTIA_END

   <Rotate Molecule>
      Write rotate-directive.

      E.g. 
         ROTATE
            (options) 
         ROTATE_END

      -options
       "type", "origin", "axis", "angle" and "angles" are prepared as options.
       These options are written in rotate-directive. Examples are shown below.

       "type" indicates whether the rotating is carring out around a axis or acording to Euler angles.
       The default value of the "type" is axis.

       "origin" indicates the origin of the rotation in angstrom unit.
       The default value of the "origin" is center of mass.

       "axis" indicates a axis around which the rotation is carried out in angstrom unit.
       The default value of the "axis" is z-axis.
       This option is valid only for "type" set as axis.

       "angle" indicates angle for the rotation around the "axis" in degree unit.
       The default value of the "angle" is 0.
       This option is valid only for "type" set as axis.

       "angles" indicates Euler angles for the rotation in degree unit.
       The default values of "angles" are 0, 0, and 0.
       This option is valid only for "type" set as Euler angles.
   

       E.g. for "type" set as axis
         ROTATE
            type axis
            origin 1.0 2.0 3.0 
            axis 3.0 4.0 5.0 
            angle 30
         ROTATE_END

       E.g. for "type" set as Euler angles
         ROTATE
            type eular_angle
            angles 15 25 35
         ROTATE_END


   <Translate Molecule>
      Write translate-directive.

      E.g. 
         TRANSLATE
            (options)
         TRANSLATE_END

      -options
       "difference" indicates difference for the translation in angstrom unit.
       This option is written in translate-directive.
       The default values are 0, 0, and 0.

       E.g. 
         TRANSLATE
            difference 12 30 45
         TRANSLATE_END



