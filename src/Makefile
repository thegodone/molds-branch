#//************************************************************************//
#// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
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
CC = icpc
LIBSBASE = -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
LIBS32 = -lmkl_intel       $(LIBSBASE)
LIBS64 = -lmkl_intel_ilp64 $(LIBSBASE)
CFLAGSBASE = -O2 -openmp -openmp-report2 -DMKL_INT=intptr_t 
CFLAGS32   = $(CFLAGSBASE)
CFLAGS64   = $(CFLAGSBASE) -DMKL_INT64
BOOST_TOP_DIR = /usr/local/boost/
BOOST_INC_DIR = $(BOOST_TOP_DIR)include/
EXENAME = MolDS.out
DEPFILE = obj/objfile.dep

ALL_CPP_FILES = Main.cpp base/MolDSException.cpp wrappers/Blas.cpp wrappers/Lapack.cpp base/Utilities.cpp base/Enums.cpp base/MathUtilities.cpp base/MallocerFreer.cpp base/EularAngle.cpp base/Parameters.cpp base/atoms/Atom.cpp base/atoms/Hatom.cpp base/atoms/Liatom.cpp base/atoms/Catom.cpp base/atoms/Natom.cpp base/atoms/Oatom.cpp base/atoms/Satom.cpp base/factories/AtomFactory.cpp base/Molecule.cpp base/InputParser.cpp base/GTOExpansionSTO.cpp base/RealSphericalHarmonicsIndex.cpp base/loggers/MOLogger.cpp base/loggers/DensityLogger.cpp base/loggers/HoleDensityLogger.cpp base/loggers/ParticleDensityLogger.cpp cndo/Cndo2.cpp indo/Indo.cpp zindo/ZindoS.cpp mndo/Mndo.cpp am1/Am1.cpp am1/Am1D.cpp pm3/Pm3.cpp pm3/Pm3D.cpp pm3/Pm3Pddg.cpp base/factories/ElectronicStructureFactory.cpp md/MD.cpp mc/MC.cpp rpmd/RPMD.cpp nasco/NASCO.cpp optimization/Optimizer.cpp optimization/ConjugateGradient.cpp optimization/SteepestDescent.cpp optimization/BFGS.cpp base/factories/OptimizerFactory.cpp  base/MolDS.cpp
ALL_HEAD_FILES = base/PrintController.h base/MolDSException.h base/Uncopyable.h wrappers/Blas.h wrappers/Lapack.h base/Utilities.h base/Enums.h base/MathUtilities.h base/MallocerFreer.h base/EularAngle.h base/Parameters.h base/atoms/Atom.h base/atoms/Hatom.h base/atoms/Liatom.h base/atoms/Catom.h base/atoms/Natom.h base/atoms/Oatom.h base/atoms/Satom.h base/factories/AtomFactory.h base/Molecule.h base/InputParser.h base/GTOExpansionSTO.h base/RealSphericalHarmonicsIndex.h base/loggers/MOLogger.h base/loggers/DensityLogger.h base/loggers/HoleDensityLogger.h base/loggers/ParticleDensityLogger.h base/ElectronicStructure.h cndo/Cndo2.h cndo/ReducedOverlapAOsParameters.h indo/Indo.h zindo/ZindoS.h mndo/Mndo.h am1/Am1.h am1/Am1D.h pm3/Pm3.h pm3/Pm3D.h pm3/Pm3Pddg.h base/factories/ElectronicStructureFactory.h md/MD.h mc/MC.h rpmd/RPMD.h nasco/NASCO.h optimization/Optimizer.h optimization/ConjugateGradient.h optimization/SteepestDescent.h optimization/BFGS.h base/factories/OptimizerFactory.h base/MolDS.h
ALL_OBJ_FILES = obj/Main.o obj/MolDSException.o obj/Blas.o obj/Lapack.o obj/Utilities.o obj/Enums.o obj/MathUtilities.o obj/MallocerFreer.o obj/EularAngle.o obj/Parameters.o obj/Atom.o obj/Hatom.o obj/Liatom.o obj/Catom.o obj/Natom.o obj/Oatom.o obj/Satom.o obj/AtomFactory.o obj/Molecule.o obj/InputParser.o obj/GTOExpansionSTO.o obj/RealSphericalHarmonicsIndex.o obj/MOLogger.o obj/DensityLogger.o obj/HoleDensityLogger.o obj/ParticleDensityLogger.o obj/Cndo2.o obj/Indo.o obj/ZindoS.o obj/Mndo.o obj/Am1.o obj/Am1D.o obj/Pm3.o obj/Pm3D.o obj/Pm3Pddg.o obj/ElectronicStructureFactory.o obj/MD.o obj/MC.o obj/RPMD.o obj/NASCO.o obj/Optimizer.o obj/ConjugateGradient.o obj/SteepestDescent.o obj/BFGS.o obj/OptimizerFactory.o obj/MolDS.o

$(EXENAME): $(DEPFILE) $(ALL_OBJ_FILES)
ifeq ($(INTEL), 64)
	$(CC) -o $@  $(ALL_OBJ_FILES) $(LIBS64) 
else
	$(CC) -o $@  $(ALL_OBJ_FILES) $(LIBS32)  
endif

-include $(DEPFILE)

.PHONY: depend
depend: $(DEPFILE)

$(DEPFILE): $(ALL_HEAD_FILES) $(ALL_CPP_FILES)
ifeq ($(INTEL), 64)
	$(CC) -MM $(ALL_CPP_FILES) -I$(BOOST_INC_DIR) | sed 's/^\([^ ]\)/obj\/\1/g' | sed 's/\($$*\)\.o[ :]*/\1.o : /g' | sed 's/h$$/h\n\t$$(CC) -I$$(BOOST_INC_DIR) -o $$@ $$< $$(CFLAGS64) -c /g' > $(DEPFILE)
else
	$(CC) -MM $(ALL_CPP_FILES) -I$(BOOST_INC_DIR) | sed 's/^\([^ ]\)/obj\/\1/g' | sed 's/\($$*\)\.o[ :]*/\1.o : /g' | sed 's/h$$/h\n\t$$(CC) -I$$(BOOST_INC_DIR) -o $$@ $$< $$(CFLAGS32) -c /g' > $(DEPFILE)
endif

.PHONY: clean
clean:
	rm -f $(ALL_OBJ_FILES) $(EXENAME) $(DEPFILE)

all: clean $(EXENAME)
