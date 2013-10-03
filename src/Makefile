#//************************************************************************//
#// Copyright (C) 2011-2013 Mikiya Fujii                                   //
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
CC = mpiicpc
CFLAGS = -O0
override CFLAGS += -openmp -DMKL_INT=intptr_t
ifeq ($(INTEL), 64)
override CFLAGS += -DMKL_ILP64
endif
BOOST_TOP_DIR = /usr/local/boost/
BOOST_INC_DIR = $(BOOST_TOP_DIR)include/
BOOST_LIB_DIR = $(BOOST_TOP_DIR)lib/
BOOST_LIBS    = -lboost_serialization -lboost_mpi -lboost_thread
LIBSBASE = -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
ifeq ($(INTEL), 64)
LIBS = -lmkl_intel_ilp64 $(LIBSBASE) $(BOOST_LIBS)
else
LIBS = -lmkl_intel       $(LIBSBASE) $(BOOST_LIBS)
endif
EXENAME = MolDS.out
DEPFILE = obj/objfile.dep
LDFLAGS =

ALL_CPP_FILES  = base/Enums.cpp                   base/PrintController.cpp base/MolDSException.cpp                                     base/MallocerFreer.cpp              mpi/MpiProcess.cpp mpi/AsyncCommunicator.cpp wrappers/Blas.cpp wrappers/Lapack.cpp base/Utilities.cpp base/MathUtilities.cpp base/EularAngle.cpp base/Parameters.cpp base/RealSphericalHarmonicsIndex.cpp base/atoms/Atom.cpp base/atoms/Hatom.cpp base/atoms/Liatom.cpp base/atoms/Catom.cpp base/atoms/Natom.cpp base/atoms/Oatom.cpp base/atoms/Satom.cpp base/factories/AtomFactory.cpp base/Molecule.cpp base/InputParser.cpp base/GTOExpansionSTO.cpp base/loggers/MOLogger.cpp base/loggers/DensityLogger.cpp base/loggers/HoleDensityLogger.cpp base/loggers/ParticleDensityLogger.cpp                            cndo/Cndo2.cpp                                    indo/Indo.cpp zindo/ZindoS.cpp mndo/Mndo.cpp am1/Am1.cpp am1/Am1D.cpp pm3/Pm3.cpp pm3/Pm3D.cpp pm3/Pm3Pddg.cpp base/factories/ElectronicStructureFactory.cpp md/MD.cpp mc/MC.cpp rpmd/RPMD.cpp nasco/NASCO.cpp optimization/Optimizer.cpp optimization/ConjugateGradient.cpp optimization/SteepestDescent.cpp optimization/BFGS.cpp optimization/GEDIIS.cpp base/factories/OptimizerFactory.cpp  base/MolDS.cpp     Main.cpp
ALL_HEAD_FILES = base/Enums.h   base/Uncopyable.h base/PrintController.h   base/MolDSException.h   base/containers/ThreadSafeQueue.h   base/MallocerFreer.h   mpi/MpiInt.h mpi/MpiProcess.h   mpi/AsyncCommunicator.h   wrappers/Blas.h   wrappers/Lapack.h   base/Utilities.h   base/MathUtilities.h   base/EularAngle.h   base/Parameters.h   base/RealSphericalHarmonicsIndex.h   base/atoms/Atom.h   base/atoms/Hatom.h   base/atoms/Liatom.h   base/atoms/Catom.h   base/atoms/Natom.h   base/atoms/Oatom.h   base/atoms/Satom.h   base/factories/AtomFactory.h   base/Molecule.h   base/InputParser.h   base/GTOExpansionSTO.h   base/loggers/MOLogger.h   base/loggers/DensityLogger.h   base/loggers/HoleDensityLogger.h   base/loggers/ParticleDensityLogger.h   base/ElectronicStructure.h cndo/Cndo2.h   cndo/ReducedOverlapAOsParameters.h indo/Indo.h   zindo/ZindoS.h   mndo/Mndo.h   am1/Am1.h   am1/Am1D.h   pm3/Pm3.h   pm3/Pm3D.h   pm3/Pm3Pddg.h   base/factories/ElectronicStructureFactory.h   md/MD.h   mc/MC.h   rpmd/RPMD.h   nasco/NASCO.h   optimization/Optimizer.h   optimization/ConjugateGradient.h   optimization/SteepestDescent.h   optimization/BFGS.h   optimization/GEDIIS.h   base/factories/OptimizerFactory.h    base/MolDS.h
ALL_OBJ_FILES  =  obj/Enums.o                      obj/PrintController.o    obj/MolDSException.o                                        obj/MallocerFreer.o                obj/MpiProcess.o   obj/AsyncCommunicator.o        obj/Blas.o        obj/Lapack.o    obj/Utilities.o    obj/MathUtilities.o    obj/EularAngle.o    obj/Parameters.o    obj/RealSphericalHarmonicsIndex.o    obj/Atom.o                obj/Hatom.o          obj/Liatom.o          obj/Catom.o          obj/Natom.o          obj/Oatom.o          obj/Satom.o              obj/AtomFactory.o    obj/Molecule.o    obj/InputParser.o    obj/GTOExpansionSTO.o            obj/MOLogger.o            obj/DensityLogger.o            obj/HoleDensityLogger.o            obj/ParticleDensityLogger.o                               obj/Cndo2.o                                       obj/Indo.o     obj/ZindoS.o    obj/Mndo.o   obj/Am1.o   obj/Am1D.o   obj/Pm3.o   obj/Pm3D.o   obj/Pm3Pddg.o              obj/ElectronicStructureFactory.o  obj/MD.o  obj/MC.o    obj/RPMD.o     obj/NASCO.o            obj/Optimizer.o            obj/ConjugateGradient.o            obj/SteepestDescent.o            obj/BFGS.o            obj/GEDIIS.o              obj/OptimizerFactory.o     obj/MolDS.o   obj/Main.o

$(EXENAME): $(ALL_OBJ_FILES)
	$(CC) -o $@ -Wl,-rpath=$(BOOST_LIB_DIR) -L$(BOOST_LIB_DIR) $(LDFLAGS) $(ALL_OBJ_FILES) $(LIBS)

-include $(DEPFILE)

.PHONY: depend
depend:
	rm $(DEPFILE)
	+make $(DEPFILE)

$(DEPFILE): Makefile
	$(CC) -MM $(ALL_CPP_FILES) -I$(BOOST_INC_DIR) | sed 's/^\([^ ]\)/obj\/\1/g' | sed 's/\($$*\)\.o[ :]*/\1.o : /g' > $(DEPFILE)

$(ALL_OBJ_FILES):
	$(CC) -I$(BOOST_INC_DIR) -o $@ $< $(CFLAGS) -c

.PHONY: clean
clean:
	rm -f $(ALL_OBJ_FILES) $(EXENAME) $(DEPFILE)

all: clean $(EXENAME)
