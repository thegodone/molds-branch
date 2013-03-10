//************************************************************************//
// Copyright (C) 2011-2012 Mikiya Fujii                                   // 
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
#ifndef INCLUDED_ATOM
#define INCLUDED_ATOM
namespace MolDS_base_atoms{

class Atom : public MolDS_base::PrintController{
public:
   Atom();
   virtual ~Atom();
   MolDS_base::AtomType GetAtomType() const;
   double GetAtomicMass() const;
   double GetCoreMass() const;
   double* GetXyz() const;
   void SetXyz(double x, double y, double z) const;
   double* GetPxyz() const;
   void SetPxyz(double px, double py, double pz) const;
   inline int GetValenceSize() const{return this->valence.size();}
   inline MolDS_base::OrbitalType GetValence(int index) const{return his->valence[index];}
   double GetVdWCoefficient() const;
   double GetVdWRadii() const;
   double GetAtomicBasisValue(double x, 
                              double y, 
                              double z, 
                              int valenceIndex,
                              MolDS_base::TheoryType theory) const;
   double GetBondingParameter() const;
   double GetBondingParameter(MolDS_base::TheoryType theory, 
                              MolDS_base::OrbitalType orbital) const;
   double GetCoreCharge() const;
   int GetFirstAOIndex() const;
   void SetFirstAOIndex(int firstAOIndex);
   int GetLastAOIndex() const;
   MolDS_base::ShellType GetValenceShellType() const;
   int GetNumberValenceElectrons() const;
   double GetOrbitalExponent(MolDS_base::ShellType shellType, 
                             MolDS_base::OrbitalType orbitalType, 
                             MolDS_base::TheoryType theory) const;  // See (1.73) in J. A. Pople book for CNDO, INDO, and ZINDOS. See [BT_1977] for MNDO. See [DZHS_1985, DY_1990] for AM1. See [S_1989] for PM3.
   double GetCoreIntegral(MolDS_base::OrbitalType orbital, 
                          double gamma, 
                          bool isGuess, 
                          MolDS_base::TheoryType theory) const; // P82 - 83 in J. A. Pople book for INDO or Eq. (13) in [BZ_1979] for ZINDO/S. See [BT_1977] for MNDO. See [DZHS_1985, DY_1990] for AM1. See [S_1989] for PM3.
   double GetCoreIntegral(MolDS_base::OrbitalType orbital, 
                          bool isGuess, 
                          MolDS_base::TheoryType theory) const;
   double GetIndoF2() const;
   double GetIndoG1() const;
   double GetZindoF0ss() const;                // Table 1 in ref. [RZ_1976], Table 1 in [AEZ_1986], or Table 1 in [GD_1972]
   double GetZindoF0sd() const;                // Table 1 in [AEZ_1986]
   double GetZindoF0dd() const;                // Table 1 in [AEZ_1986]
   double GetZindoG1sp() const;                // Table 3 in ref. [BZ_1979]
   double GetZindoF2pp() const;                // Table 3 in ref. [BZ_1979]
   double GetZindoG2sd() const;                // Table 3 in ref. [BZ_1979]
   double GetZindoG1pd() const;                // Table 3 in ref. [BZ_1979]
   double GetZindoF2pd() const;                // Table 3 in ref. [BZ_1979]
   double GetZindoG3pd() const;                // Table 3 in ref. [BZ_1979]
   double GetZindoF2dd() const;                // Table 3 in ref. [BZ_1979]
   double GetZindoF4dd() const;                // Table 3 in ref. [BZ_1979]
   double GetZindoF0ssLower() const;           // Apendix in ref. [BZ_1979] 
   double GetZindoF0sdLower() const;           // Apendix in ref. [BZ_1979]
   double GetZindoF0ddLower() const;           // Apendix in ref. [BZ_1979]
   double GetZindoG1spLower() const;           // Apendix in ref. [BZ_1979]
   double GetZindoF2ppLower() const;           // Apendix in ref. [BZ_1979]
   double GetZindoG2sdLower() const;           // Apendix in ref. [BZ_1979]
   double GetZindoG1pdLower() const;           // Apendix in ref. [BZ_1979]
   double GetZindoF2pdLower() const;           // Apendix in ref. [BZ_1979]
   double GetZindoG3pdLower() const;           // Apendix in ref. [BZ_1979]
   double GetZindoF2ddLower() const;           // Apendix in ref. [BZ_1979]
   double GetZindoF4ddLower() const;           // Apendix in ref. [BZ_1979]
   double GetZindoIonPot(MolDS_base::OrbitalType orbital) const;
   double GetNddoAlpha(MolDS_base::TheoryType theory) const; // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S for MNDO. Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S for AM1. [S_1989] for PM3.
   double GetNddoDerivedParameterD(MolDS_base::TheoryType theory, 
                                   MolDS_base::MultipoleType multipole) const;    // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S. Or, calculated in tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double GetNddoDerivedParameterRho(MolDS_base::TheoryType theory, 
                                     MolDS_base::MultipoleType multipole) const;  // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S. Or, calculated in tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double GetMndoElecEnergyAtom() const;        // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double GetMndoHeatsFormAtom() const;         // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double GetNddoGss(MolDS_base::TheoryType theory) const;
   double GetNddoGpp(MolDS_base::TheoryType theory) const;
   double GetNddoGsp(MolDS_base::TheoryType theory) const;
   double GetNddoGpp2(MolDS_base::TheoryType theory) const;
   double GetNddoHsp(MolDS_base::TheoryType theory) const;
   double GetNddoHpp(MolDS_base::TheoryType theory) const;
   double GetNddoParameterK(MolDS_base::TheoryType theory, int kIndex) const;//Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S for AM1. [S_1989] for PM3.
   double GetNddoParameterL(MolDS_base::TheoryType theory, int lIndex) const;//Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S for AM1. [S_1989] for PM3.
   double GetNddoParameterM(MolDS_base::TheoryType theory, int mIndex) const;//Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S for AM1. [S_1989] for PM3.
   double GetPm3PddgParameterPa(int paIndex) const;
   double GetPm3PddgParameterDa(int daIndex) const;
protected:
   double* xyz; // coordinates
   double* pxyz; // momentum. Note that this is not velocity!! 
   MolDS_base::AtomType atomType;
   double atomicMass;  // Appendix 1 in [I_1998]
   std::vector<MolDS_base::OrbitalType> valence;
   MolDS_base::ShellType valenceShellType;
   int firstAOIndex;
   int numberValenceElectrons;
   double vdWCoefficient;               // Table 1 in [G_2004] and [G_2006]
   double vdWRadii;                     // Table 1 in [G_2004] and [G_2006]
   double imuAmuS;                      // Table 3.4 or 3.5 in J. A. Pople book
   double imuAmuP;                      // Table 3.4 or 3.5 in J. A. Pople book
   double imuAmuD;                      // Table 3.4 or 3.5 in J. A. Pople book
   double bondingParameter;             // Table 3.2 and 3.4 in J. A. Pople book
   double coreCharge;                   // = Z_A in J. A. Pople book.
   double effectiveNuclearChargeK;      // Table 1.5 in J. A. Pople book or table 1 in [HKLWNZ_1982]
   double effectiveNuclearChargeL;      // Table 1.5 in J. A. Pople book or table 1 in [HKLWNZ_1982]
   double effectiveNuclearChargeMsp;    // Table 1.5 in J. A. Pople book
   double effectiveNuclearChargeMd;     // Table 1.5 in J. A. Pople book
   double indoF2;                   // Table 3.6 in J. A. Pople book
   double indoG1;                   // Table 3.6 in J. A. Pople book
   double indoF0CoefficientS;       // (3.93-3.99) in J. A. Pople book
   double indoF0CoefficientP;       // (3.93-3.99) in J. A. Pople book
   double indoG1CoefficientS;       // (3.93-3.99) in J. A. Pople book
   double indoG1CoefficientP;       // (3.93-3.99) in J. A. Pople book
   double indoF2CoefficientS;       // (3.93-3.99) in J. A. Pople book
   double indoF2CoefficientP;       // (3.93-3.99) in J. A. Pople book
   double zindoBondingParameterS;        // Table 1 in [RZ_1976], table 1 in [HKLWNZ_1982], or table 3 in [AEZ_1986]
   double zindoBondingParameterD;        // Table 1 in [RZ_1976], table 1 in [HKLWNZ_1982], or table 3 in [AEZ_1986]
   double zindoF0ss;                // Table 1 in ref. [RZ_1976], Table 1 in [AEZ_1986], or Table 1 in [GD_1972]
   double zindoF0sd;        // Table 1 in [AEZ_1986]
   double zindoF0dd;        // Table 1 in [AEZ_1986]
   double zindoG1sp;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoF2pp;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoG2sd;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoG1pd;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoF2pd;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoG3pd;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoF2dd;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   double zindoF4dd;       // Table 3 in ref. [BZ_1979] or table 1 in [HKLWNZ_1982]
   int zindoL;              // see l of (13) in [BZ_1979]
   int zindoM;              // see m of (13) in [BZ_1979]
   int zindoN;              // see n (13) in [BZ_1979]
   double zindoIonPotS;   // Ionization potential, Table 4 in [BZ_1979]
   double zindoIonPotP;   // Ionization potential, Table 4 in [BZ_1979]
   double zindoIonPotD;   // Ionization potential, Table 4 in [BZ_1979]
   double mndoCoreintegralS;         // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoCoreintegralP;         // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S. 
   double mndoOrbitalExponentS;      // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoOrbitalExponentP;      // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoBondingParameterS;     // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoBondingParameterP;     // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoAlpha;        // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoDerivedParameterD[3];    // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S. Or, calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double mndoDerivedParameterRho[3];  // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S. Or, calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double mndoElecEnergyAtom;        // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoHeatsFormAtom;         // Table III in ref. [DT_1977-2] for H, B, C, N, O, and F. Table I & II in ref. [DMR_1978] and Table I in ref. [DR_1986] for S.
   double mndoGss;   //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double mndoGpp;   //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double mndoGsp;   //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double mndoGpp2;  //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double mndoHsp;   //Table I in ref. [BDL_1975] for H, B, C, N, O, F, Si, P, S, and Cl.
   double am1CoreintegralS; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1CoreintegralP; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1OrbitalExponentS;// Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1OrbitalExponentP;// Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1BondingParameterS; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1BondingParameterP; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1Alpha;// Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1Gss; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1Gpp; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1Gsp; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1Gpp2; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1Hsp; // Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1DerivedParameterD[3];    // Calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double am1DerivedParameterRho[3];  // Calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double am1ParameterK[4];// Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1ParameterL[4];// Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1ParameterM[4];// Table I in ref. [DZHS_1985] for H, C, N, O, and Table I in re. [DY_1990] for S.
   double am1DCoreintegralS; // Table II in ref. [MH_2007] for H, C, N, and O and Table IV in ref. [MMHBV_2007] for S.
   double am1DCoreintegralP; // Table II in ref. [MH_2007] for H, C, N, and O and Table IV in ref. [MMHBV_2007] for S.
   double am1DBondingParameterS; // Table II in ref. [MH_2007] for H, C, N, O, and Table IV in re. [MMHBV_2007] for S.
   double am1DBondingParameterP; // Table II in ref. [MH_2007] for H, C, N, O, and Table IV in re. [MMHBV_2007] for S.
   double am1DAlpha; // Table II in ref. [MH_2007] for H, C, N, O, and Table IV in re. [MMHBV_2007] for S.
   double pm3CoreintegralS; // Table II in ref. [S_1989].
   double pm3CoreintegralP; // Table II in ref. [S_1989].
   double pm3OrbitalExponentS;// Table II in ref. [S_1989].
   double pm3OrbitalExponentP;// Table II in ref. [S_1989].
   double pm3BondingParameterS; // Table II in ref. [S_1989].
   double pm3BondingParameterP; // Table II in ref. [S_1989].
   double pm3Alpha;// Table II in ref. [S_1989].
   double pm3DerivedParameterD[3];    // Calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double pm3DerivedParameterRho[3];  // Calculated by tools/deriveParametersNDDO/deriveParametersNDDO.cpp.
   double pm3ParameterK[4];// Table II in ref. [S_1989].
   double pm3ParameterL[4];// Table II in ref. [S_1989].
   double pm3ParameterM[4];// Table II in ref. [S_1989].
   double pm3Gss; // Table II in ref. [S_1989].
   double pm3Gpp; // Table II in ref. [S_1989].
   double pm3Gsp; // Table II in ref. [S_1989].
   double pm3Gpp2; // Table II in ref. [S_1989].
   double pm3Hsp; // Table II in ref. [S_1989].
   double pm3PddgCoreintegralS; // Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgCoreintegralP; // Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgOrbitalExponentS;// Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgOrbitalExponentP;// Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgBondingParameterS; // Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgBondingParameterP; // Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgAlpha;// Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgDerivedParameterD[3];    // Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgDerivedParameterRho[3];  // Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgParameterK[4];// Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgParameterL[4];// Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgParameterM[4];// Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgParameterPa[2];// Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3PddgParameterDa[2];// Table II in ref. [RCJ_2002] for H, C, N, O, and Table II in re. [BGJ_2003] for S.
   double pm3DCoreintegralS; // Table II in ref. [MH_2007] for H, C, N, and O and Table IV in ref. [MMHBV_2007] for S.
   double pm3DCoreintegralP; // Table II in ref. [MH_2007] for H, C, N, and O and Table IV in ref. [MMHBV_2007] for S.
   double pm3DBondingParameterS; // Table II in ref. [MH_2007] for H, C, N, O, and Table IV in re. [MMHBV_2007] for S.
   double pm3DBondingParameterP; // Table II in ref. [MH_2007] for H, C, N, O, and Table IV in re. [MMHBV_2007] for S.
   double pm3DAlpha; // Table II in ref. [MH_2007] for H, C, N, O, and Table IV in re. [MMHBV_2007] for S.
private:
   std::string errorMessageIonPot;
   std::string errorMessageAtomType;
   std::string errorMessageNumberValences;
   std::string errorMessageValenceIndex;
   std::string errorMessageOrbitalType;
   std::string errorMessageOrbitalExponent;
   std::string errorMessageShellType;
   std::string errorMessageEffectivPrincipalQuantumNumber;
   std::string errorMessageCndo2CoreIntegral;
   std::string errorMessageIndoCoreIntegral;
   std::string errorMessageZindoCoreIntegral;
   std::string errorMessageMndoCoreIntegral;
   std::string errorMessageAm1CoreIntegral;
   std::string errorMessageAm1DCoreIntegral;
   std::string errorMessagePm3CoreIntegral;
   std::string errorMessagePm3DCoreIntegral;
   std::string errorMessagePm3PddgCoreIntegral;
   std::string errorMessageGetAtomicBasisValueBadValenceIndex;
   std::string errorMessageGetRealAngularPartAOBadValence;
   std::string errorMessageGetOrbitalExponentBadTheory;
   std::string errorMessageTheoryType;
   std::string errorMessageGetBondingParameterBadTheoryBadOrbital;
   std::string errorMessageGetNddoAlphaBadTheory;
   std::string errorMessageGetNddoDerivedParameterDBadTheory;
   std::string errorMessageGetNddoDerivedParameterDBadMultipoleType;
   std::string errorMessageMultipoleType;
   std::string errorMessageGetNddoDerivedParameterRhoBadMultipoleType;
   std::string errorMessageGetNddoDerivedParameterRhoBadTheory;
   std::string errorMessageRhoIndex;
   std::string errorMessageGetNddoParameterKBadKIndex;
   std::string errorMessageGetNddoParameterKBadTheory;
   std::string errorMessageKIndex;
   std::string errorMessageGetNddoParameterLBadLIndex;
   std::string errorMessageGetNddoParameterLBadTheory;
   std::string errorMessageLIndex;
   std::string errorMessageGetNddoParameterMBadMIndex;
   std::string errorMessageGetNddoParameterMBadTheory;
   std::string errorMessageMIndex;
   std::string errorMessageGetPm3PddgParameterPaBadPaIndex;
   std::string errorMessagePaIndex;
   std::string errorMessageGetPm3PddgParameterDaBadDaIndex;
   std::string errorMessageDaIndex;
   std::string errorMessageGetNddoGssBadTheory;
   std::string errorMessageGetNddoGppBadTheory;
   std::string errorMessageGetNddoGspBadTheory;
   std::string errorMessageGetNddoGpp2BadTheory;
   std::string errorMessageGetNddoHspBadTheory;
   std::string errorMessageGetNddoHppBadTheory;
   std::string errorMessageGetXyzCoordinatesNull;
   std::string errorMessageSetXyzCoordinatesNull;
   std::string errorMessageGetPxyzMomentaNull;
   std::string errorMessageSetPxyzMomentaNull;
   void SetMessages();
   double GetRealAngularPartAO(double theta, 
                               double phi, 
                               MolDS_base::OrbitalType orbital) const;
   double GetRadialPartAO(double dr, 
                          double orbitalExponent, 
                          MolDS_base::ShellType shell) const;
   int GetEffectivePrincipalQuantumNumber(MolDS_base::ShellType shellType) const; // Table 1.4 in J. A. Pople book
   double GetZindoJss() const;  // Part of Eq. (13) in [BZ_1979]
   double GetZindoJsp() const;  // Part of Eq. (13) in [BZ_1979]
   double GetZindoJsd() const;  // Part of Eq. (13) in [BZ_1979]
   double GetZindoJpp() const;  // Part of Eq. (13) in [BZ_1979]
   double GetZindoJpd() const;  // Part of Eq. (13) in [BZ_1979]
   double GetZindoJdd() const;  // Part of Eq. (13) in [BZ_1979]
   double GetCndo2CoreIntegral(MolDS_base::OrbitalType orbital, double gamma, bool isGuess) const;
   double GetIndoCoreIntegral(MolDS_base::OrbitalType orbital, double gamma, bool isGuess) const;
   double GetZindoCoreIntegral(MolDS_base::OrbitalType orbital) const; // Eq. (13) in [BZ_1979]
   double GetMndoCoreIntegral(MolDS_base::OrbitalType orbital) const; 
   double GetAm1CoreIntegral(MolDS_base::OrbitalType orbital) const; 
   double GetAm1DCoreIntegral(MolDS_base::OrbitalType orbital) const; 
   double GetPm3CoreIntegral(MolDS_base::OrbitalType orbital) const; 
   double GetPm3DCoreIntegral(MolDS_base::OrbitalType orbital) const; 
   double GetPm3PddgCoreIntegral(MolDS_base::OrbitalType orbital) const; 
   virtual void SetAtomicParameters() = 0;
};
}
#endif

