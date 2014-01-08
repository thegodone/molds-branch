//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   // 
// Copyright (C) 2012-2014 Katsuhiko Nishimra                             // 
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
#ifndef INCLUDED_ENUMS
#define INCLUDED_ENUMS
#ifndef RENUMSTR_BODY
#define RENUMSTR_BODY 0
#endif
#include "../third_parties/rEnumStr/rEnumStr.h"

namespace MolDS_base{

RENUMSTR_BEGIN( SimulationType, SimulationTypeStr )
   RENUMSTR( Once,  "Once" )
   RENUMSTR( MD,    "MD" )
   RENUMSTR( MC,    "MC" )
   RENUMSTR( RPMD,  "RPMD" )
   RENUMSTR( NASCO, "NASCO" )
   RENUMSTR( PrincipalAxes, "PrincipalAxes" )
   RENUMSTR( Translate, "Translate" )
   RENUMSTR( Rotate, "Rotate" )
   RENUMSTR( Optimization, "Optimization" )
   RENUMSTR( SimulationType_end,  "SimulationType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( TheoryType, TheoryTypeStr )
   RENUMSTR( CNDO2,  "CNDO/2" )
   RENUMSTR( INDO,   "INDO" )
   RENUMSTR( ZINDOS, "ZINDO/S" )
   RENUMSTR( MNDO,   "MNDO" )
   RENUMSTR( AM1,    "AM1" )
   RENUMSTR( AM1D,    "AM1-D" )
   RENUMSTR( PM3,    "PM3" )
   RENUMSTR( PM3D,    "PM3-D" )
   RENUMSTR( PM3PDDG,    "PM3/PDDG" )
   RENUMSTR( TheoryType_end,  "TheoryType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( OptimizationMethodType, OptimizationMethodTypeStr )
   RENUMSTR( ConjugateGradientMethod,  "Conjugate gradient" )
   RENUMSTR( SteepestDescentMethod,  "Steepest descent" )
   RENUMSTR( BFGSMethod,  "BFGS" )
   RENUMSTR( GEDIISMethod,  "GEDIIS" )
   RENUMSTR( OptimizationMethodType_end,  "OptimizationMethodType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( RotatingType, RotatingTypeStr )
   RENUMSTR( Axis,  "Axis" )
   RENUMSTR( Eular,  "EularAngle" )
   RENUMSTR( RotatingType_end,  "RotatingType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( RotatedObjectType, RotatedObjectTypeStr )
   RENUMSTR( System,  "System" )
   RENUMSTR( Frame,  "Frame" )
   RENUMSTR( RotatedObjectType_end,  "RotatedObjectType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( ShellType, ShellTypeStr )
   RENUMSTR( kShell,  "k-shell" )
   RENUMSTR( lShell,  "l-shell" )
   RENUMSTR( mShell,  "m-shell" )
   RENUMSTR( nShell,  "n-shell" )
   RENUMSTR( ShellType_end,  "ShellType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( CartesianType, CartesianTypeStr )
   RENUMSTR( XAxis,  "XAxis" )
   RENUMSTR( YAxis,  "YAxis" )
   RENUMSTR( ZAxis,  "ZAxis" )
   RENUMSTR( CartesianType_end,  "CartesianType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( AzimuthalType, AzimuthalTypeStr )
   RENUMSTR( sAzimuthal,  "s-azimuthal-quantum-number" )
   RENUMSTR( pAzimuthal,  "p-azimuthal-quantum-number" )
   RENUMSTR( dAzimuthal,  "d-azimuthal-quantum-number" )
   RENUMSTR( AzimuthalType_end,  "AzimuthalType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( OrbitalType, OrbitalTypeStr )
   RENUMSTR( s,  "s" )
   RENUMSTR( py,  "py" )
   RENUMSTR( pz,  "pz" )
   RENUMSTR( px,  "px" )
   RENUMSTR( dxy,  "dxy" )
   RENUMSTR( dyz,  "dyz" )
   RENUMSTR( dzz,  "dzz" )
   RENUMSTR( dzx,  "dzx" )
   RENUMSTR( dxxyy,  "dxxyy" )
   RENUMSTR( OrbitalType_end,  "OrbitalType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( AtomType, AtomTypeStr )
   RENUMSTR( H,    "H" )
   RENUMSTR( He,   "He" )
   RENUMSTR( Li,   "Li" )
   RENUMSTR( Be,   "Be" )
   RENUMSTR( B,    "B" )
   RENUMSTR( C,    "C" )
   RENUMSTR( N,    "N" )
   RENUMSTR( O,    "O" )
   RENUMSTR( F,    "F" )
   RENUMSTR( Ne,   "Ne" )
   RENUMSTR( Na,   "Na" )
   RENUMSTR( Mg,   "Mg" )
   RENUMSTR( Al,   "Al" )
   RENUMSTR( Si,   "Si" )
   RENUMSTR( P,    "P" )
   RENUMSTR( S,    "S" )
   RENUMSTR( Cl,   "Cl" )
   RENUMSTR( Ar,   "Ar" )
   RENUMSTR( Zn,   "Zn" )
   RENUMSTR( ghostH,  "ghost-H" )
   RENUMSTR( ghostHe, "ghost-He" )
   RENUMSTR( ghostLi, "ghost-Li" )
   RENUMSTR( ghostBe, "ghost-Be" )
   RENUMSTR( ghostB,  "ghost-B" )
   RENUMSTR( ghostC,  "ghost-C" )
   RENUMSTR( ghostN,  "ghost-N" )
   RENUMSTR( ghostO,  "ghost-O" )
   RENUMSTR( ghostF,  "ghost-F" )
   RENUMSTR( ghostNe, "ghost-Ne" )
   RENUMSTR( ghostNa, "ghost-Na" )
   RENUMSTR( ghostMg, "ghost-Mg" )
   RENUMSTR( ghostAl, "ghost-Al" )
   RENUMSTR( ghostSi, "ghost-Si" )
   RENUMSTR( ghostP,  "ghost-P" )
   RENUMSTR( ghostS,  "ghost-S" )
   RENUMSTR( ghostCl, "ghost-Cl" )
   RENUMSTR( ghostAr, "ghost-Ar" )
   RENUMSTR( ghostZn, "ghost-Zn" )
   RENUMSTR( EPC,  "Environmental Point Charge" )
   RENUMSTR( AtomType_end,  "AtomType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( STOnGType, STOnGTypeStr )
   RENUMSTR( STO1G,  "STO1G" )
   RENUMSTR( STO2G,  "STO2G" )
   RENUMSTR( STO3G,  "STO3G" )
   RENUMSTR( STO4G,  "STO4G" )
   RENUMSTR( STO5G,  "STO5G" )
   RENUMSTR( STO6G,  "STO6G" )
   RENUMSTR( STOnGType_end,  "STOnGType_end" )
RENUMSTR_END()

// For the definition of the MultipopleType, see appendix in [DT_1977].
RENUMSTR_BEGIN( MultipoleType, MultipoleTypeStr )
   RENUMSTR( sQ,  "q(small Q)" )
   RENUMSTR( Qxx, "Qxx" )
   RENUMSTR( Qyy, "Qyy" )
   RENUMSTR( Qzz, "Qzz" )
   RENUMSTR( Qxz, "Qxz" )
   RENUMSTR( Qyz, "Qyz" )
   RENUMSTR( Qxy, "Qxy" )
   RENUMSTR( mux, "mux" )
   RENUMSTR( muy, "muy" )
   RENUMSTR( muz, "muz" )
   RENUMSTR( MultipoleType_end,  "MultipoleType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( ExceptionKey, ExceptionKeyStr )
   RENUMSTR( LapackInfo, "LapackInfo" )
   RENUMSTR( EmptyQueue, "EmptyQueue" )
   RENUMSTR( GEDIISErrorID, "GEDIISErrorID" )
   RENUMSTR( ExceptionKey_end,  "ExceptionKey_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( MpiFunctionType, MpiFunctionTypeStr )
   RENUMSTR( Send,      "Send" )
   RENUMSTR( Recv,      "Recv" )
   RENUMSTR( Broadcast, "Broadcast" )
   RENUMSTR( MpiFunctionType_end,  "MpiFunctionType_end" )
RENUMSTR_END()

RENUMSTR_BEGIN( GEDIISErrorID, GEDIISErrorStr )
   RENUMSTR( GEDIISNotSufficientHistory, "GEDIISNotSufficientHistory" )
   RENUMSTR( GEDIISNegativeCoefficient, "GEDIISNegativeCoefficient" )
   RENUMSTR( GEDIISLapackInfo, "GEDIISLapackInfo" )
   RENUMSTR( GEDIISErrorID_end,  "GEDIISErrorID_end" )
RENUMSTR_END()

}
#endif

