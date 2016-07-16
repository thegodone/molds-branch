//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   // 
// Copyright (C) 2012-2014 Katsuhiko Nishimra                             // 
// Copyright (C) 2013-2014 Michihiro Okuyama                              // 
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

   MolDS ("Mol"ecular "D"ynamics simulation package with "S"emiempirical quantum chemistry) ver. 0.4.0 (Under developement)
      Developers: Mikiya Fujii, Ph.D.(project lead), Katsuhiko Nishimra, and Michihiro Okuyama, Ph.D..
      Other contributors: Michael Banck 
      Questions and bug reports: molds-dev@lists.sourceforge.jp


   The goal of MolDS is classical/semiclassical/quantum molecular dynamics simulations 
   based on semiempirical quantum chemistry. Besides, MolDS is implemented to be able 
   to run on huge cluster system by using OpenMP/MPI hybrid parallelization technique.

==============================================================================
PLATFORMS: 
   Linux clusters, K-computer, FX10, and Mac OS X.

==============================================================================
REQUIREMENTS:
   -Compilers:
    MolDS requires C++ mpi compiler (e.g. Intel MPI or Open MPI) 
    that is wrapping Intel (icpc with MKL) or GNU (g++) C++ compiler.
    Valid versions of the mpi compilers are Intel MPI 4.0.2, Open MPI 1.4.5, or later.
    Valid versions of the wrapped c++ compilers are icpc 12.0.4(MkL 10.3 update 4), 
    g++ 4.4, or later because the MolDS is implemented with openMP 3.0. 
    For the K-computer and FX10, Fujitsu C++ compiler, mpiFCCpx, is available.

   -Boost C++ Libraries
    Boost C++ Libraries  builded with MPI is needed.
    To get and install the Boost, see the HP:<http://www.boost.org/>.
    The version of the boost would be no problem if 1.46.0 or later is used.
    Especially, the Boost should be builded with MPI 
    because MolDS needs boost_mpi-library(i.e. -lboost_mpi).
    An example of manually building of the boost 1.48.0 by M.F. is shown in:
    http://d.hatena.ne.jp/futofuji/20120320/p2

   -Linear Algebra Packages (i.e. BLAS and LAPACK)
    MolDS needs a linear algebra package. In the current implementation of MolDS, 
    MKL (Intel's Math Kernel Library) or OpenBLAS is assumed as the linear algebra package 
    for the Intel or GNU compilers, respectively.
    See also the section of compilers about the version of the MKL.
    To get and install the OpenBLAS-libraries, see the HP:<http://xianyi.github.com/OpenBLAS/>.
    The version of the OpenBLAS would be no problem if 0.2.5 or later is used.
    Note that "USE_OPENMP = 1" should be set for the installation of the opneBLAS.
    Furthermore, "BINARY = 64" and "INTERFACE64 = 1" are also needed 
    when you install the OpenBLAS into 64-bits machines.
    An example of manually building of the openBLAS 0.2.5 by M.F. is shown in:
    http://d.hatena.ne.jp/futofuji/20130627/p1

    In the case that you compile molds on the K-computer or FX10, 
    the linear algebra package builted in K-computer or FX10 should be used.
    So, you do not need to install BLAS and LAPACK by yourself.

==============================================================================
COMPILE:
   GNUMake is used to compile the MolDS in the "src" directory of the MolDS package.
   MolDS officially suport the following cases.

   I) Linux clusters
      Case i) The Intel mpi compiler (mpiicpc) wrapping the Intel c++ compiler (icpc)
         Change the "BOOST_TOP_DIR" in Makefile to the top directory of the 
         Boost C++ Libraries in your systems.

         To compile MolDS on 32 bits machine,
         $ make INTEL=32

         To compile MolDS on 64 bits machine,
         $ make INTEL=64
      
      Case ii) The openMPI compiler (mpicxx) wrapping the Intel c++ compiler (icpc)
         Change the "BOOST_TOP_DIR" in Makefile to the top directory of the 
         Boost C++ Libraries in your systems.

         To compile MolDS on 32 bits machine,
         $ make INTEL=32 CC=mpicxx

         To compile MolDS on 64 bits machine,
         $ make INTEL=64 CC=mpicxx
      
      Case iii) The openMPI compiler (mpicxx) wrapping the GNU c++ compiler (g++):
         Change the "BOOST_TOP_DIR" in "Makefile_GNU" to the top directory of the 
         Boost C++ Libraries in your systems.
         Change the "OPENBLAS_TOP_DIR" in "Makefile_GNU" to the top directory of the 
         OpneBLAS in your systems.
         
         Then, just type: 
         $ make -f Makefile_GNU 

   II) K-computer or FX10
         Fujitsu C++ compiler on K-computer of FX10
         Change the "BOOST_TOP_DIR" in Makefile to the top directory of the 
         Boost C++ Libraries in the K-computer or FX10. 

         Then, just type: 
         $ make -f Makefile_FX10

   III) Mac OS X
         For mac OS X, the openMPI wrapping Intel compiler is only tested. 
         The dynamic library of boost should be added to your .bashrc,
         i.e., write following line to your .bashrc:
            export DYLD_LIBRARY_PATH=/usr/local/boost/lib:$DYLD_LIBRARY_PATH
         so, LDFLAGS in the makefiles should be empty.
         Besides, "-fvisibility=hidden" is needed to be added to CFLAGS in the makefiles.
         
         e.g., (for openMPI wrapping Intel compiler)
         $ make INTEL=64 CC=mpicxx CFLAGS="-O2 -fvisibility=hidden" LDFLAGS=""

   For above all cases, the compile succeeded if you could fine "molds" in the "src" directory. 
   If you want to clean the compilation, type 
   $ make clean

   If you want to compile MolDS in debug-mode, 
   -g, -rdynamic(for function names in backtrace) and -DMOLDS_DBG should be added to CFLAGS,
   namely, hit the following command:
   $ make CFLAGS="-O0 -g -rdynamic -DMOLDS_DBG"

==============================================================================
CARRY OUT MolDS:
   After the compile, in the "src" directory,

   For the calculations with single process:
   $ ./molds < input.in
   or
   $ ./molds input.in

   For the calculations with muliple threads, type
   $ export OMP_NUM_THREADS=n1
   $ ./molds input.in
   , where n1 is the number of threads.

   For the calculations with multiple processes by MPI:
   $ mpirun -np n2 molds input.in
   , where n2 after the "-np" is the number of process.

   For the calculations with muliple threads and muliple processes, type
   $ export OMP_NUM_THREADS=n1
   $ mpirun -np n2 molds input.in
   , where n1 is the number of cores of each node and n2 is the number of nodes.

   In the multiple processes calculations, process-0 can only output results.
   If you want to get all output from the all processes, 
   -DMOLDS_DBG should be added to CFLAGS at the compilation.
   Then, make only one process on each node and output results to 
   node unique file (e.g. local file system of each node.),
   namely, 
   $ make CFLAGS="-DMOLDS_DBG"
   $ export OMP_NUM_THREADS=n1
   $ mpirun -np n2 molds input.in > /localFileSyste/output.dat
   , where n1 is the number of cores of each node and n2 is the number of nodes.

==============================================================================
SAMPLE and TEST
   See sample files in "test" directory or 
   "http://sourceforge.jp/projects/molds/scm/svn/tree/head/trunk/test/"
   In the "test" directory, *.in files are input files, then *.dat files are
   associated output files. To execute all test cases, carry out below ruby-script
   in the "test" directory. This script will finished in a few minutes with big
   output(a few thausands lines).

   $ ruby Test_Of_MolDS.rb

   To execute some specific test cases, carry out below ruby-script with the test names you want. 
   When the test names are specifed, ".in" and ".dat" in the arguments will be ignored.

   $ ruby Test_Of_MolDS.rb test1.in test2.dat test3 ...

   Note that this test script needs at least 4 cores.

==============================================================================
CAPABILITIES:

   -Electronic state and molecular dynamics
             | HF  | CIS |   MD    |   MC    |  RPMD   | Nonadiabatic | Optimize | Frequencies      | Frequencies       |
             |     |     | (gs,es) | (gs,es) | (gs,es) | Ehrenfest    | (gs,es)  | analytic (gs,es) | numerical (gs,es) |
    ---------|-----|-----|---------|---------|---------|--------------|----------|------------------|-------------------|
    CNDO2    | OK  | --  | --,--   | OK,--   | --,--   |      --      | --,--    | --,--            | --,--             |
    ---------|-----|-----|---------|---------|---------|--------------|----------|------------------|-------------------|
    INDO     | OK  | --  | --,--   | OK,--   | --,--   |      --      | --,--    | --,--            | --,--             |
    ---------|-----|-----|---------|---------|---------|--------------|----------|------------------|-------------------|
    ZINDO/S  | OK  | OK  | OK,OK   | OK,OK   | OK,OK   |      OK      | OK,OK    | --,--            | OK,OK             |
    ---------|-----|-----|---------|---------|---------|--------------|----------|------------------|-------------------|
    MNDO     | OK  | OK  | OK,OK   | OK,OK   | OK,OK   |      OK      | OK,OK    | OK,--            | OK,OK             |
    ---------|-----|-----|---------|---------|---------|--------------|----------|------------------|-------------------|
    AM1      | OK  | OK  | OK,OK   | OK,OK   | OK,OK   |      OK      | OK,OK    | OK,--            | OK,OK             |
    ---------|-----|-----|---------|---------|---------|--------------|----------|------------------|-------------------|
    AM1-D    | OK  | OK  | OK,OK   | OK,OK   | OK,OK   |      OK      | OK,OK    | OK,--            | OK,OK             |
    ---------|-----|-----|---------|---------|---------|--------------|----------|------------------|-------------------|
    PM3      | OK  | OK  | OK,OK   | OK,OK   | OK,OK   |      OK      | OK,OK    | OK,--            | OK,OK             |
    ---------|-----|-----|---------|---------|---------|--------------|----------|------------------|-------------------|
    PM3-D    | OK  | OK  | OK,OK   | OK,OK   | OK,OK   |      OK      | OK,OK    | OK,--            | OK,OK             |
    ---------|-----|-----|---------|---------|---------|--------------|----------|------------------|-------------------|
    PM3/PDDG | OK  | OK  | OK,OK   | OK,OK   | OK,OK   |      OK      | OK,OK    | OK,--            | OK,OK             |
                                                                      
      "OK", "Sch", and "--" mean available, shceduled, and non-scheduled methods, respectively.
      "gs" and "es" mean ground and excited states, respectively.
      i.e., MD(gs) and MD(es) mean Born-Oppenheimer Molecular Dynamics on ground and excited states, respectively. 

   -Elements
    CNDO2    | H, Li, C, N, O, F, S, Cl
    INDO     | H, Li, C, N, O, F, 
    ZINDO/S  | H,     C, N, O, F, S, Cl
    MNDO     | H,     C, N, O, F, S, Cl, Zn
    AM1      | H,     C, N, O, F, S, Cl, Zn
    AM1-D    | H,     C, N, O, F, S, Cl       *
    PM3      | H,     C, N, O, F, S, Cl, Zn
    PM3-D    | H,     C, N, O, F, S, Cl       *
    PM3/PDDG | H,     C, N, O,    S 

    * For fluorine(F) and chlorine(Cl), semiempirical parameters used in SCF of AM1-D and PM3-D are same with 
      the original semiempirical parameters used in SCF of AM1 and PM3, respectively.

   -Parallelization
    Open MP parallelization: Everywhere in MolDS
    MPI parallelization: Everywhere in MolDS, basically. But some modules do not parallelized yet.

==============================================================================
HOW TO WRITE INPUT:

   <Terminiology>
      "hoge-directive" means line block stating "hoge" and ending "hoge_end" in input files.
      Uppercase and lowercase letters are treated as identical in input files.

   <Comment Out>
      Lines starting with "//" or "#" in input-files are treated as comments.

   <Theory>
      Write "cndo/2", "indo", "zindo/s", "mndo", "am1", "am1-d",
      "pm3", "pm3-d", or "pm3/pddg" in theory-directive.
      This theory-directive indicate a electronic structure theory used in your simulations.
      MNDO only supports (can calculate) Heats of formation.

      E.g. 
         THEORY
            indo 
         THEORY_END

   <SCF>
      Write SCF-directive. In the SCF-directive, settings(options) of SCF should be written.
      SCF-module in theories after MNDO requires (number of atoms)**2 * 4**4 * 8 byte for two electron integral.

      E.g.
         SCF
            (options)
         SCf_END
   
      -options
       Write below options in SCF-directive.
       "max_iter", "rms_density", "damping_thresh", "damping_weight", 
       "diis_num_error_vect", "diis_start_error", "diis_end_error", "two_elec_int", 
       "vdW", "vdW_s6", "vdW_d", "sum_charges", "mpi", and "scalapack" are prepared as options.
       SCF module outputs also the dipole moment arrond the center of core's mass
       To calculate the dipole moment, STO-6G [DY_1977] is used.

       The default value of the "max_iter" is 100.
       The default value of the "rms_density" is 10**(-8.0).
       The default value of the "damping_thresh" is 1.
       The default value of the "damping_weight" is 0.8. For systems consisting of several hundreds of atoms, the developers recommend setting "damping_weight" more than 0.95.
       The default value of the "diis_num_error_vect" is 5.
       The default value of the "diis_start_error" is 0.01.
       The default value of the "diis_end_error" is 10**(-8.0).

       "two_elec_int" should be set as "on_node" or "direct". When the "on_node" is specifed, 
       all two-core two-electron (2C2E) integrals are stored on every node. Namely, all node have all 2C2E integrals.
       So, the "on_node" consumes a large amount of memory, but the "on_node" is fast. 
       On the otherhand, when "direct" is specifed, no 2C2E element is stored on memory nor disk.
       Namely the 2C2E integrals are calulated every time when they are required. 
       The distributed storing of 2C2E integrals on multiple nodes is not implemented, unfortunately. 
       The default value of the "two_elec_int" is "on_node"

       "vdW" should be set as "yes" or "no". 
       When "yes" is set, Grimmes's empirical van der Waals correction(D1, [G_2004]) is applied.
       Note that this empirical van der Waals correction is applied to the semiempirical theories 
       of which semiempirical parameters are not modified.
       This "vdw" option can be used for H, C, N, O, F, S, and Cl.
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

       "sum_charges" is an option to calculate of summation of atomic charges in the ground state.
       To use this option, write "sum_charges first_atom_index last_atom_index"
       in the SCF-directive. The index starts from 0 for atoms written in geometry-directive.
       The atoms indicated with first_atom_index and last_atom_index are also included in the atoms which charge is summed.
       Multiple setting of this "sum_charges" option is approvable, of course.
       If you want to calculate summation, same "sum_charges" option is available in CIS-directive.

       "mpi" is an option to indicate wheter the MPI is used or not in SCF module
       even if MolDS is carried on multiple processors with command "mpirun(or mpiexec)".
       This "mpi" should be set as yes or no as. So to use this option, 
       write "mpi yes" or "mpi no" in SCF-directive.
       If "mpi" is set to "no", SCF module does not use the MPI.
       This option can change the setting of parallelization for SCF module only. Namely, 
       this option cannot change the setting of parallelization for 1st(2nd) gradients, CIS, etc.
       The default value of the "mpi" is "yes". 

       "scalapack" should be set as "yes" or "no". 
       When "scalapack yes", the diagonalizing SCF matrix is carried out with ScaLapack.
       The default setting of "scalapack" option is "NO".
       Note that "scalapack" is supported on FX10 only.

       E.g.
         SCF
            max_iter 200
            rms_density 1e-8
            two_elec_int on_node
            damping_thresh 0.1
            damping_weight 0.7
            diis_num_error_vect 6
            diis_start_error 0.01
            diis_end_error 1e-8
            sum_charges 0 4
            sum_charges 3 7
            vdW yes
            vdW_s6 0.75
            vdW_d 30
            mpi no
         SCF_END

   <GEOMETRY>
      To set geometry(configuration) of a system calculated by MolDS, 
      the configuration should be written in geometry-directive.
      Each line inside the geometry-directive indicates each atom of the system.  
      Namely, each line should containe one character and three doubles.
      The character indicates atomtype and three doubles indicate the cartesian coordinates of
      each atom in angstrom unit.

       E.g.
         GEOMETRY
            C -0.1000 0.1000 0.0000
            C 1.6938 0.0000 -0.1000
            H -0.381 1.1411 0.0000
            H -0.2681 -0.5205 -0.9016
            H -0.3681 -0.4725 0.8016
            H 1.9519 0.5200 -0.9007
            H 1.8519 0.5300 0.8007
            H 1.7519 -1.0401 -0.1000
         GEOMETRY_END
         
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

   <MO Plot>
      Write MO plot directive. Note that All threads draw different MOs. 
      In the hybrid parallelization with openMP/MPI, each MO is written 
      to the local disk of a node on which each process are running. 
      Therefore, staging-out should be carried out for the drawn MOs. 

      E.g.
         MOPLOT 
            (options)
         MOPLOT_END

      -options
       "mo", "mos", "grid_number", "frame_length", and "file_prefix" are prepared.

       "mo" is index of the molcular orbital. mo=0 means the lowest energy MO.
       The default value of the "mo" is not set.

       "mos" is a set of indeces of the molcular orbital. 
       To draw many and sequencial MOs, this "mos" is prepared as
       "mos first_mo_index last_mo_index".
       The default value of the "mos" is also not set

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
            mos 4 8
            grid_number 30 30 30
            frame_length 10 10 10
            file_prefix MOPlot_
         MOPLOT_END

   <Environmental Point Charge(EPC) method>
      Environmental point charge method is a simplified method of the QM/MM,
      namely the environmental point changes are treated as atoms in the MM region.
      The differences between the QM/MM and EPC are summarized below:
         - Electrostatic interaction between QM and MM region:
            QM/MM: Electrostatic interaction may be mutually added to QM and MM atoms.
            EPC  : Electrostatic field caused by the EPCs affects the QM region
                   although the each EPC is not affected by electrostatic field 
                   caused by the QM atoms and other EPCs. 
                   Namely, each EPC is fixed at point of space.
         - Van der Waals interaction between QM and MM region:
            QM/MM: Included.
            EPC  : Not included.
      In this EPC method, core-core replustion between QM and MM atoms is implemented with
      the method I (simple coulomb interaction: qq/r) of ref [LRCL_2000].
      This EPC method can be used with MNDO-series (MNDO, AM1, AM1-D, PM3, PM3-D, and PDDG/PM3) only.
      To use this environmental point charges method, write EPC-directive. 

      E.g.
         EPC
            (options)
         EPC_END
      
      -options
       "the cartesian coordinates and charge" is only prepared.
       Namely, each line should containe 4 doubles and a term. 
       The first three doubles are the cartesian coordinates of 
       each environmental point charge in angstrom unit.
       The term is "charge". The last double following the term, "charge", 
       is the charge in atomic unit, 
       e.g. -1 and 1 mean charge of an electron and a proton, respectively.
       Multiple setting of the environmental point charge is approvable, of course.

       E.g. 
         EPC 
            0.0 0.0 0.0 charge -1   
            2.2 1.5 3.0 charge -1.5
            0.0 2.0 5.0 charge  0.3
         EPC_END

   <Frequencies (Normal modes analysis)>
      write frequencies-directive. Note taht not only the frequencies but also the normal modes are calculated.

      E.g.
         FREQUENCIES
            (options)
         FREQUENCIES_END

      -options
       "electronic_state", "derivative", "numerical_dr", "projection", "projection_dphi" are prepared.

       "electronic_state" is index of the electronic state used for calculating the normal modes. 
       electronic_state=0 means the electronic ground state.
       electronic_state=1 means, then, first electornic excited state.
       The default value of the "electronic_state" is 0.

       "derivative" is how to calculate the hessian matrix, then
       "derivative" should be "analytic" or "numerical".
       Note that the analytic derivative can calcute 
       the hessian matrix for the grond state only.
       The default value of the "derivative" is "analytic".

       "numerical_dr" means numerical displacement to calculate the Hessian matrix.
       "numerical_dr" is valid in the numerical derivative only.
       "numerical_dr" should be set in atomic units (a.u.).
       The default value of the "numerical_dr" is 0.0001 a.u.(0.00053 angstrom).
       
       "projection" means projecting the translational and rotational 
       modes out from mass-weighted Hessian matrix.
       "projection" should be set as "yes" or "no". 
       The default value of the "projection" is "yes".

       "projection_dphi" means rotational angle to make projection vector 
       for the rotatinal modes.
       "projection_dphi" should be set in radian units.
       The default value of the "projection_dphi" is 0.001 radian.

       E.g. 
         FREQUENCIES
            electronic_state 0
            derivative numerical
         FREQUENCIES_END

   <CIS>
      Write CIS-directive. This module requires at least (active_occ*active_vir)**2 * 8 byte.

      E.g.
         CIS
            (options)
         CIS_END
   
      -options
       "davidson", "active_occ", "active_vir", "max_iter", "max_dim", "norm_tol", 
       "nstates", "exciton_energies", "all_transition_dipole_moments", 
       "mulliken", "unpaired_electron_population", "sum_charges", "scalapack", and "num_print_coefficients" are prepared as options.

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
       The center of the transition dipole moments is same with 
       the center of the dipole moment of the ground state. 
       The dipole moments and transition dipole moments are calculated 
       with GTO expansion of STO orbiltals(ST0-6G in [S_1970] is used).
       Besides, oscillator strength is also shown with each transition dipole moment.

       "mulliken" is an option of mulliken popultaion analysis of the excited state.
       When "mulliken x" is included in CIS-directive, the mulliken popultaion of xth excited state is calculated.
       Multiple indication of these mulliken options is possible. 
       Note that "mulliken 0" is ignored because 0th excited state is the ground state.
       The default setting of this "mulliken" option is nothing.

       "mullikens" is a set of indeces of the mulliken population analysis in the excites states.
       To calculate mulliken population for many and sequencial excited states, 
       this "mullikens" is prepared as
       "mullikens first_elec_index last_elec_index".
       The default value of the "mullikens" is also not set

       "unpaired_electron_population" is an option of unpaired electron population(UEP) analysis of the excited state.
       When "unpaired electron population yes" and "mulliken x" (amd/or "mullikens x y") are included in CIS-directive, 
       the UEP of xth excited state is calculated.
       By multiple indication of these mulliken option, the UEP on multiple excited states are possible.
       Note that the UEP on ground state is ignored. 
       The default setting is "unpaired_electron_population" option is nothing.

       "scalapack" should be set as "yes" or "no". 
       When "scalapack yes", the diagonalizing CIS matrix is carried out with ScaLapack.
       The Default setting of "scalapack" option is "NO".
       Note that "scalapack" is supported on FX10 only.

       "sum_charges" is an option to calculate of summation of the mulliken charges and unpaired electron population
       in each excited state which indicated with "mulliken" option.
       How to set this option is same way with the ground state.

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
            mullikens 1 3
	         unpaired_electron_population yes
            sum_charges 0 4
            sum_charges 3 7
         CIS_END

   <Hole Plot>
      Write hole plot directive for the output of the density of the hole.
      This module is valid only when CIS is carried out.
      Note that All threads draw different hole deinsities. 
      In the hybrid parallelization with openMP/MPI, each hole densitis is written 
      to the local disk of a node on which each process are running. 
      Therefore, staging-out should be carried out for the drawn hole deinsities. 
      This module requires at least (active_occ+active_vir)*(grid_number_x, grid_number_y, grid_number_z) * 8 byte.

      E.g.
         HOlEPLOT
            (options)
         HOLEPLOT_END

      -options
       "electronic_state", "electronic_states", "grid_number", "frame_length", and "file_prefix" are prepared.

       "electronic_state" is index of the electronic state. 
       electronic_state=0 means the electronic ground state.
       electronic_state=1 means, then, first electornic excited state.
       The default value of the "electronic_states" is not set.

       "electronic_states" is a set of indeces of the hole densities. 
       To draw many and sequencial densities, this "electronic_states" is prepared as
       "electronic_states first_electronic_index last_electronic_index".
       The default value of the "electronic_states" is also not set

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
            electronic_states 10 12
            grid_number 30 30 30
            frame_length 10 10 10
            file_prefix HOLEPlot_
         HOLEPLOT_END

   <Particle Plot>
      Write particle plot directive for the output of the density of 
      the particle which lives avobe Fermi's sea.
      This module is valid only when CIS is carried out.
      Note that All threads draw different particle deinsities. 
      In the hybrid parallelization with openMP/MPI, each particle densitis is written 
      to the local disk of a node on which each process are running. 
      Therefore, staging-out should be carried out for the drawn particle deinsities. 
      This module requires at least (active_occ+active_vir)*(grid_number_x, grid_number_y, grid_number_z) * 8 byte.

      E.g.
         PARTICLEPLOT
            (options)
         PARTICLEPLOT_END

      -options
       "electronic_state", "electronic_states", "grid_number", "frame_length", and "file_prefix" are prepared.

       "electronic_state" is index of the electronic state. 
       electronic_state=0 means the electronic ground state.
       electronic_state=1 means, then, first electornic excited state.
       The default value of the "electronic_state" is not set.

       "electronic_states" is a set of indeces of the particle densities. 
       To draw many and sequencial densities, this "electronic_states" is prepared as
       "electronic_states first_electronic_index last_electronic_index".
       The default value of the "electronic_states" is also not set

       "grid_number" is the grid number of the frame in xyz-coordinates.
       The default values are 25, 25, and 25 for x, y, and z coordinates, respectively.

       "frame_length" is the length of the frame of each coordinate.
       The default values are 10, 10, and 10[angst.] for x, y, and z coordinates.

       "file_prefix" is a prefix of the file name to which the density of the particle is written.
       The default values is "particle_".

      E.g.
         PARTICLEPLOT
            electronic_state 5
            electronic_state 8
            electronic_states 10 12
            grid_number 30 30 30
            frame_length 10 10 10
            file_prefix PARTPlot_
         PARTICLEPLOT_END

   <OPT (geometry optimization)>
      Write optimization-directive. This module uses line search and steepest descent algorythms.
      In the early stage the line search algorythm is used, 
      then the algorythm used in this module is switched to steepest descent algorythm.
      Note that ZINDO/S is not suitable for geometry optimizations.

      E.g.
         OPTIMIZATION
            (options)
         OPTIMIZATION_END
  
      -options
       "method", "total_steps", "electronic_state", "max_gradient", "rms_gradient", 
        "dt", "initial_trust_radius", "max_norm_step", and "space_fixed_atoms" are prepared as options.

       "method" should be set as "conjugate_gradient", "steepest_descent", "bfgs" or "gediis". 
       The default of the "method" is conjugate gradient.

       "electronic_state" means the electronic eigenstate 
       on which the system runs.
       The default value of the "electronic_state" is 0. That is, 
       electronic ground state is default.

      "line_search_times" means the times of line-search trials.
      The default value of the "line_search_times" is 50.
      This parameter have no effect if method is "bfgs" or "gediis".

      "steep_step" means the number of steps of the steepest descent.
      The default value of the "steep_step" is 50.

      "max_gradient" and "rms_gradient" are threshold of the steepest descent.
      The "max_gradient" and "rms_gradient" means maximum and root-mean-squre of the gradient, respectively.
      The default value of the "max_gradient" is 0.00045.
      The default value of the "rms_gradient" is 0.00030.

      "dt" is initial fictious time steps for the steepest descent algorythms.
      The default value of the "dt" is 50[fs].
      This parameter have no effect if method is "bfgs" or "gediis".

      "initial_trust_radius" is an initial value for trust radius used by RFO step method.
      The default value of the "initial_trust_radius" is 0.3.
      This parameter have no effect if method is "steepest_descent" or "conjugate_gradient".

      "max_norm_step" is the maximum value for trust radius used by RFO ssep method.
      The default value of the "max_norm_step" is 0.3.
      This parameter have no effect if method is "steepest_descent" or "conjugate_gradient".

      "space_fixed_atom" or "space_fixed_atoms" mean specification of 
      fixed atom or atoms in space during optimization.
      The implementation to fix the atom or atoms is ignoreing force of specified atoms.
      To use this option, write
      "space_fixed_atom atom_index"
      or
      "space_fixed_atoms first_atom_index last_atom_index"
      in the OPTIMIZATION-directive.
      Note that the index starts from 0 for atoms written in geometry-directive.
      The atoms indicated with first_atom_index and last_atom_index are also included in the atoms fixed to space.
      Multiple setting of this "space_fixed_atom" or "space_fixed_atoms" option is approvable, of course.

      E.g.
         OPTIMIZATION
            method steepest_descent
            total_steps 50
            electronic_state 0
            max_gradient 0.00045
            rms_gradient 0.00030
            dt 50
            space_fixed_atoms 2 4 // 2nd 3rd 4th atoms are fixed. Note that first atom is indexed with 0.
            space_fixed_atoms 5 5 // Only 5th atom is fixed.      Note that first atom is indexed with 0.
            space_fixed_atom  6   // Only 6th atom is fixed.      Note that first atom is indexed with 0.
         OPTIMIZATION_END

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
       "total_steps", "electronic_state", "num_electronic_states", "temperature", 
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

   <Ehrenfest MD>
      Write Ehrenfest-directive to run the Ehrenfest molecular dynamics. 
      The Ehrenfest MD is a semiclassical treatment known as mean-field dynamics for nonadiabatic systems.
      *Implementation of Ehrenfest in MolDS is based on MF's unpublished algorythm.

      E.g.
         Ehrenfest 
            (options)
         Ehrenfest_END
  
      -options
       "initial_electronic_state","highest_electronic_state",
       "lowest_electronic_state", "total_steps", and "dt" are prepared as options.

       "initial_electronic_state" means the electronic eigenstate 
       on which the system start to run.
       The default value of the "electronic_state" is 1. That is, 
       the electronic ground state is default.

       "highest_electronic_state" and "lowest_electronic_state" are
       the highest and lowest adiabatic surfaces, respectively, which are included calculations.
       Namely, a number of adiabatic surfaces used in the calculations is 
       "highest_electronic_state - lowest_electronic_state +1".
       Note that a state indicated by "highest(lowest)_electronic_state 0" means the ground states.
       The default values of "highest_electronic_state" and "lowest_electronic_state"
       are 2 and 0, respectively.

       The default value of the "total_steps" is 10. 

       "dt" means the time width of molecular dynamics.
       "dt" should be set in femto-second.
       The default value of the "dt" is 0.1[fs].

      E.g.
         EHRENFEST
            initial_electronic_state 1   // Ehrenfest start from the first electronic state
            highest_electronic_state 10
            lowest_electronic_state  0   // this is the ground state.
            total_steps 50
            dt 0.05
         EHRENFEST_END


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



