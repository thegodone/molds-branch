/*---------------------------------------------------------------------------
    Copyright (C) 2004  rage2050  rage2050@mail.goo.ne.jp
    (�{��MD5/UTF-8: ab93b9b4f0bec5870b90ceb98d08b7e2)

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

// �� ���C�Z���X�� GPL �ɂ��Ă܂����A����Ȃ̂͂������č��x�ȃe�N���Ă킯�ł�
// �Ȃ��̂ŁA���̃t�@�C�����u���̂܂�܁v�g���̂łȂ���΁A���l�̃e�N��
// PROPRIETARY�ȃR�[�h�ɓ���邱�Ƃ͉\�Ƃ��܂��B
// ���Ȃ��Ƃ����͂���PROPRIETARY�ȃR�[�h�̏��L�҂�i�����肵�܂���B
// (���̐l�����̃e�N��Ɛ肵�悤�Ƃ��āArEnumStr.h ���g���Ă��鑼�̐l��
//  �i����c�ȂǂƂ��������Ȃ��Ƃ����Ȃ������!)

//---------------------------------------------------------------------------
/** @file
    @brief enum �ɑΉ����镶����𓯎��ɒ�`����}�N���Q

enum �ƁA����ɑΉ����镶������A�����ɒ�`���邱�Ƃ��ł��܂��B
���̂��߁Aenum �ƕ����񂪁A�Y������H��������葫��Ȃ�������A�Ƃ����~�X��h�����Ƃ��ł��܂��B

<hr>

@section whatis �����֗��Ȃ�?

@subsection old_problem �]���̕��@�̌��_

enum �萔�Ƃ���ɑΉ����镶������`�������Ƃ��A���ʂ͈ȉ��̂悤�ɏ����܂��B�F�̔ԍ��Ɩ��O���`�����ł��B

@code
  //////////////////// �w�b�_�t�@�C�� mycolor.h ////////////////////
  #ifndef MYCOLOR_H_INCLUDE
  #define MYCOLOR_H_INCLUDE
  enum MyColor { ///< �F��\�� enum
    C_BLACK,
    C_BLUE,
    C_RED,
    C_GREEN,
    C_WHITE,
  };
  extern char * gMyColor[C_MAX]; ///< �F��enum�ɑΉ�����������
  void DrawCircleWithColor(MyColor c, int x, int y, int r);///< �`��c�Ȃ�
  #endif //MYCOLOR_H_INCLUDE
@endcode

@code
  //////////////////// �\�[�X�t�@�C�� mycolor.cpp ////////////////////
  #include "mycolor.h"
  char * gMyColor[C_MAX] = {
    "��", // C_BLACK
    "��", // C_BLUE
    "��", // C_RED
    "��", // C_GREEN
    "��", // C_WHITE
  };
@endcode

�������A���̕��@�ɂ́Aenum��ǉ��E�폜�E�����ύX���邽�тɁA�Ή����镶����� <B>�Y�ꂸ�ɁE�ԈႦ��</B> �ɏC�����Ȃ���΂Ȃ�Ȃ��A�Ƃ������_������܂��B

rEnumStr.h ���g���ƁAenum �ƕ�����𓯎��ɒ�`���邱�Ƃ��ł��āAenum ��ǉ��E�폜�E�����ύX�����Ƃ��ɕ�������C�����Y�ꂽ��ԈႦ���肷�邱�Ƃ�h���܂��B

<hr>

@subsection improved rEnumStr.h ���g���Ƃǂ��Ȃ��?

��قǂ̐F�̔ԍ��Ɩ��O�̗�́A�ȉ��̂悤�ɏ����܂��B

@code
  //////////////////// �w�b�_�t�@�C�� mycolor.h ////////////////////
  //---- ��^�I�Ȃ��܂��Ȃ� ----
  #ifndef MYCOLOR_DEFINED
  #define MYCOLOR_DEFINED
  #ifndef RENUMSTR_BODY
  #define RENUMSTR_BODY 0
  #endif
  #include "rEnumStr.h"
  //---- �F�̔ԍ��ƕ�������` ----
  RENUMSTR_BEGIN( MyColor, MyColorStr )
    RENUMSTR( C_BLACK,  "��" )
    RENUMSTR( C_BLUE,   "��" )
    RENUMSTR( C_RED,    "��" )
    RENUMSTR( C_GREEN,  "��" )
    RENUMSTR( C_WHITE,  "��" )
  RENUMSTR_END()
  #endif //MYCOLOR_DEFINED
  //---- �������畁�ʂ̃w�b�_�t�@�C�� ----
  #ifndef MYCOLOR_H_INCLUDE
  #define MYCOLOR_H_INCLUDE
  void DrawCircleWithColor(MyColor c, int x, int y, int r);///< �`��c�Ȃ�
  #endif //MYCOLOR_H_INCLUDE
@endcode

@code
  //////////////////// �\�[�X�t�@�C�� mycolor.cpp ////////////////////
  #include "mycolor.h"
  #undef MYCOLOR_DEFINED        // ���������荞�ނ��߂�
  #define RENUMSTR_BODY 1       // �ݒ��ς���
  #include "mycolor.h"          // ������xinclude���܂��B
@endcode



<hr>

@section usage �g����

@subsection intro �g���n�߂̂��܂��Ȃ�

@code
  //////////////////// �w�b�_�t�@�C�� mycolor.h ////////////////////
  //---- ��^�I�Ȃ��܂��Ȃ� ----
  #ifndef MYCOLOR_DEFINED
  #define MYCOLOR_DEFINED
  #ifndef RENUMSTR_BODY
  #define RENUMSTR_BODY 0
  #endif
  #include "rEnumStr.h"
    :
@endcode

���G�Ɍ����܂����ArEnumStr.h ���g�������Ȃ����Ƃ��͂��̕������R�s�y���āA MYCOLOR_DEFINED �}�N���̖��O��K���ɒu�����邾���ő��v�ł��B


@subsection main ��`���Ă��镔��

@code
    :
  //---- �F�̔ԍ��ƕ�������` ----
  RENUMSTR_BEGIN( MyColor, MyColorStr )
    RENUMSTR( C_BLACK,  "��" )
    RENUMSTR( C_BLUE,   "��" )
    RENUMSTR( C_RED,    "��" )
    RENUMSTR( C_GREEN,  "��" )
    RENUMSTR( C_WHITE,  "��" )
  RENUMSTR_END()
  #endif //MYCOLOR_DEFINED
@endcode

�ŏ��� RENUMSTR_BEGIN( enum�^�� , ������֐��� ) �Œ�`���J�n���A�Ō�� RENUMSTR_END() �Œ�`���I���܂��B

���̊Ԃ� RENUMSTR( enum�� , ������ ) ����ׂāAenum �ƕ�����𓯎��ɒ�`���܂��B



@subsection cpp ��`�����̉����镔��

@code
  //////////////////// �\�[�X�t�@�C�� mycolor.cpp ////////////////////
  #include "mycolor.h"
  #undef MYCOLOR_DEFINED        // ���������荞�ނ��߂�
  #define RENUMSTR_BODY 1       // �ݒ��ς���
  #include "mycolor.h"          // ������xinclude���܂��B
@endcode

1��ڂ� include �� enum MyColor �^�Ƃ��̗v�f C_BLACK, C_BLUE, ... ����`����܂��B1��ڂ� include �ł� RENUMSTR_BODY �}�N���� 0 �ɂȂ邽�߁A mycolor.h ���ȉ��̂悤�ɓW�J����邩��ł��B
@code
  ////////// 1��ڂ͂��̂悤�ɓW�J����� //////////
  enum MyColor {
    C_BLUE,
    C_RED,
    C_GREEN,
    C_WHITE,
  };
  const char * MyColorStr(int);
@endcode

2��ڂ� include �ŁAenum�ɑΉ����镶�����Ԃ��֐� MyColorStr() ������܂��B2��ڂ� mycolor.cpp �Ń}�N�� RENUMSTR_BODY �� 1 �ɒ�`���Ă���̂ŁA mycolor.h ��W�J����ƁA�ȉ��̂悤�Ȋ֐� MyColorStr() ����������܂��B

@code
  ////////// 2��ڂ͂��̂悤�ɓW�J����� //////////
  const char * MyColorStr(int i)
  {
    switch (i) {
    default: return "???";
    case C_BLUE: return "��";
    case C_RED: return "��";
    case C_GREEN: return "��";
    case C_WHITE: return "��";
    }
  }
@endcode

@subsection invoke enum�╶��������o���ɂ�

mycolor.h ���}�N���W�J����� enum ��`�ƕ�����擾�֐� MyColorStr() �̐錾���ł��Ă���̂ŁA����� include ���邾���ŊȒP�� enum ����������g����悤�ɂȂ�܂��B

@code
  //////////////////// ���C�����[�`�� main.cpp ////////////////////
  #include <stdio.h>
  #include "mycolor.h"
  int main(int, char**)
  {
    MyColor c = C_BLUE;
    printf( "%d = %s\n", c, MyColorStr(c) );
    return 0;
  }
@endcode



<hr>

@attention

������̎��̒�` On/Off �� enum���Ƃɐݒ肷�邽�߂ɁA
RENUMSTR_BODYFLAG �}�N���́u���[�J���v�ɂȂ��Ă��܂��B
(rEnumStr.h��include���邽�т� \#undef RENUMSTR_BODYFLAG �����)
����� RENUMSTR_BODYFLAG �}�N���́A<B>�v���W�F�N�g�S�̂Őݒ肷��̂ł͂Ȃ�</B>�A
\#include "rEnumStr.h" ���g�����тɂ��̒��O�Œ�`���Ă��������B
�v���W�F�N�g�ݒ�� RENUMSTR_BODYFLAG ���`����Ɨ\�z�O�̋����ɂȂ�܂��B

RENUMSTR_BODYFLAG �}�N���͕K���u�\�[�X�t�@�C���v�� \#include "rEnumStr.h" ��
���O�Œ�`���Ă��������B�w�b�_�t�@�C����v���W�F�N�g�ݒ�Œ�`�����
�R���p�C���G���[�ɂȂ�܂��B

*/
//---------------------------------------------------------------------------
//#ifndef rEnumStrH
//#define rEnumStrH



//---------------------------------------------------------------------------
#ifndef RENUMSTR_BODY //{
#error #include "rEnumStr.h" �̑O�� RENUMSTR_BODY ����`����Ă��܂���

// doxygen�p�R�����g�͂����ɏ����Ȃ��ƃ_��? {

/// enum������\�̒�`�J�n�B
  /** @param nameEnum          enum�^�̖��O
      @param nameStr           ������擾�֐��̖��O

      �\�������I�������K�� RENUMSTR_END() ��u���Ă��������B
   */
#define RENUMSTR_BEGIN(nameEnum, nameStr)

  /** @brief enum�ƕ�������ЂƂ��B
      @param idEnum            enum�萔��
      @param idStr             enum�ɑΉ��Â��镶����
   */
#define RENUMSTR(idEnum, idStr)

  /** @brief enum���ЂƂ��A�����enum�����̂��̂̕���������蓖�Ă�B
      @param idEnum            enum�萔��
   */
#define RENUMSTR1(idEnum)

  /** @brief enum�ɒl��ݒ肵�A����enum�ƕ������Ή�������B
      @param idEnum            enum�萔��
      @param val               enum�ɐݒ肷��l
      @param idStr             enum�ɑΉ��Â��镶����
  */
#define RENUMSTR_LET(idEnum, val, idStr)

  /** @brief enum�̒��Ɏ��R�Ȉ�s�������B
      @param line              �s�̓��e

      enum�萔�ɕʖ�������ȂǂɎg���܂��B
  */
#define RENUMSTR_ENUMONLY(line)

  /** @brief enum������\���I������B

      �\�������I�������K�� RENUMSTR_END() ��u���Ă��������B
   */
#define RENUMSTR_END()

  /** @brief [���[�U��`] [�K�{] rEnumStr.h�̓��샂�[�h��ݒ�B

      rEnumStr.h ���g���w�b�_/�\�[�X�t�@�C���ł́A
      #include "rEnumStr.h" �̑O�ɕK�����̃}�N�����`���Ă��������B

      - 0 �ɂ���ƁAenum�錾�ƕ�����擾�֐��̃v���g�^�C�v�錾�𐶐����܂��B
        (�w�b�_���[�h)
      - 1 �ɂ���ƁA������擾�֐��̎��̂𐶐����܂��B
        (�\�[�X���[�h)
   */
#define RENUMSTR_BODY

  /** @brief [���[�U��`] ��������̂̒�`��On/Off����B

      �\�[�X�t�@�C������ #include "rEnumStr.h" �̒��O�ł��̃}�N����
      ��`����ƁA��������̂̒�`�� On/Off �ł��܂��B
      �����[�X�łł͕�������`���Ȃ��A�Ȃǂ̗p�r�Ɏg���܂��B

      - 0 �ɂ���ƁA������𐶐����܂���B
      - 1 �ɂ���ƁA������𐶐����܂��B
   */
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
/// ������̎��̂��`���邩? (���[�J��)
/** �����񂪃f�o�b�O�p�ł��胊���[�X�łł͒�`���Ȃ��Ȃǂ̏ꍇ�A
    �����[�X�ł̂Ƃ��� #define RENUMSTR_BODYFLAG 0 �ɂ��邱�Ƃ�
    ��������`���Ȃ��悤�ɂł��܂��B
 */
#define RENUMSTR_BODYFLAG 1
#endif

#ifndef RENUMSTR_DEFAULTSTR
/// enum�͈͊O�̂Ƃ���default������ (���[�J��)
/** ������擾�֐��ɔ͈͊O�̈������n���ꂽ�Ƃ��ɕԂ���������w�肵�܂��B
 */
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

// RENUMSTR_BODYFLAG ���Ԉ���ăw�b�_�Œ�`�����Ƃ��̓G���[�ɂ���
// �v���W�F�N�g�ݒ�Œ�`�����Ƃ����G���[�ɂȂ����Ⴄ���c
#ifdef RENUMSTR_BODYFLAG
#error RENUMSTR_BODYFLAG �́u�\�[�X�t�@�C���v�Œ�`���Ă�������
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
