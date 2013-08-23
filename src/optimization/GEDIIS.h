//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   //
// Copyright (C) 2012-2013 Katsuhiko Nishimra                             //
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
#ifndef INCLUDED_GEDIIS
#define INCLUDED_GEDIIS

#include <list>

namespace MolDS_optimization{

class GEDIIS : public MolDS_optimization::BFGS{
public:
   GEDIIS();
   ~GEDIIS();
protected:
   void SetMessages();

   std::string messageStartGEDIISStep;

   class GEDIISHistory{
   public:
      GEDIISHistory();
      ~GEDIISHistory();
      void AddEntry(double energy,
                    const MolDS_base::Molecule& molecule,
                    double const* const* matrixForce);
      void SolveGEDIISEquation()const;
   private:
      class Entry{
      public:
         Entry(double energy,
               const MolDS_base::Molecule& molecule,
               double const* const* matrixForce);
         ~Entry();
      private:
         int      numAtoms;
         double   energy;
         double** matrixCoordinate;
         double** matrixForce;
      };
      typedef std::list< const Entry* > entryList_t;
      entryList_t entryList;
   };

private:
   virtual void SearchMinimum(boost::shared_ptr<MolDS_base::ElectronicStructure> electronicStructure,
                              MolDS_base::Molecule& molecule,
                              double* lineSearchedEnergy,
                              bool* obainesOptimizedStructure) const;
};

}
#endif
