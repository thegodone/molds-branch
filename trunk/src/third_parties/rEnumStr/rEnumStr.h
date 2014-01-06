//************************************************************************//
// Copyright (C) 2011-2014 Mikiya Fujii                                   // 
//                                                                        // 
// This file is part of MolDS.                                            // 
// This file is a modified version of the original rEnumStr.h.            // 
// The modification is carried in 2013/03/01 to delete Japanese.          // 
// The original rEnumStr.h is distributed at                              //
// http://www.geocities.jp/rage2050a/rEnumStr/rEnumStr.html               //
// The copyright of the original rEnumStr.h is shown below.               //
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

/*---------------------------------------------------------------------------
    Copyright (C) 2004  rage2050  rage2050@mail.goo.ne.jp
    (name MD5/UTF-8: ab93b9b4f0bec5870b90ceb98d08b7e2)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
---------------------------------------------------------------------------*/

//---------------------------------------------------------------------------
//#ifndef rEnumStrH
//#define rEnumStrH

//---------------------------------------------------------------------------
#ifndef RENUMSTR_BODY //{
#error RENUMSTR_BODY should be defined before #include "rEnumStr.h"
#define RENUMSTR_BEGIN(nameEnum, nameStr)
#define RENUMSTR(idEnum, idStr)
#define RENUMSTR1(idEnum)
#define RENUMSTR_LET(idEnum, val, idStr)
#define RENUMSTR_ENUMONLY(line)
#define RENUMSTR_END()
#define RENUMSTR_BODY
#define RENUMSTR_BODYFLAG

// doxygen }
#else //}{

//---------------------------------------------------------------------------
#if RENUMSTR_BODY //{

#undef RENUMSTR_BEGIN
#undef RENUMSTR
#undef RENUMSTR1
#undef RENUMSTR_LET
#undef RENUMSTR_ENUMONLY
#undef RENUMSTR_END

#ifndef RENUMSTR_BODYFLAG
#define RENUMSTR_BODYFLAG 1
#endif
#ifndef RENUMSTR_DEFAULTSTR
#define RENUMSTR_DEFAULTSTR "???"
#endif


#define RENUMSTR_BEGIN(nameEnum, nameStr)          const char *nameStr(int i) { \
                                                    switch (i) { \
                                                    default: return RENUMSTR_DEFAULTSTR;
#if RENUMSTR_BODYFLAG //{
#define RENUMSTR(idEnum, idStr)                case idEnum: return idStr;
#define RENUMSTR1(idEnum)                       case idEnum: return #idEnum;
#define RENUMSTR_LET(idEnum, val, idStr)        case idEnum: return idStr;
#else //}{
#define RENUMSTR(idEnum, idStr)
#define RENUMSTR1(idEnum)
#define RENUMSTR_LET(idEnum, val, idStr)
#endif //}
#define RENUMSTR_ENUMONLY(line)

#define RENUMSTR_END()                               } \
                                                  }

//---------------------------------------------------------------------------
#else //}{ !RENUMSTR_BODY

#undef RENUMSTR_BEGIN
#undef RENUMSTR
#undef RENUMSTR1
#undef RENUMSTR_LET
#undef RENUMSTR_ENUMONLY
#undef RENUMSTR_END

#ifdef RENUMSTR_BODYFLAG
#error RENUMSTR_BODYFLAG should be difeined in the source files
#endif

#define RENUMSTR_BEGIN(nameEnum, nameStr)          const char *nameStr(int); \
                                                  enum nameEnum {
#define RENUMSTR(idEnum, idStr)                idEnum,
#define RENUMSTR1(idEnum)                       idEnum,
#define RENUMSTR_LET(idEnum, val, idStr)        idEnum=val,
#define RENUMSTR_ENUMONLY(line)                      line,
#define RENUMSTR_END()                             };

//---------------------------------------------------------------------------
#endif //}

//#undef RENUMSTR_DEFAULTSTR
#undef RENUMSTR_BODYFLAG
#undef RENUMSTR_BODY
//---------------------------------------------------------------------------
#endif //}
//#endif
