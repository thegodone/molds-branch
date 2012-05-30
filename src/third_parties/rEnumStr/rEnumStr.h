/*---------------------------------------------------------------------------
    Copyright (C) 2004  rage2050  rage2050@mail.goo.ne.jp
    (本名MD5/UTF-8: ab93b9b4f0bec5870b90ceb98d08b7e2)

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

// ↑ ライセンスは GPL にしてますが、こんなのはたいして高度なテクってわけでも
// ないので、このファイルを「そのまんま」使うのでなければ、同様のテクを
// PROPRIETARYなコードに入れることは可能とします。
// 少なくとも私はそのPROPRIETARYなコードの所有者を訴えたりしません。
// (その人がこのテクを独占しようとして、rEnumStr.h を使っている他の人を
//  訴える…などという愚かなことをしない限りは!)

//---------------------------------------------------------------------------
/** @file
    @brief enum に対応する文字列を同時に定義するマクロ群

enum と、それに対応する文字列を、同時に定義することができます。
このため、enum と文字列が、ズレたり食い違ったり足りなかったり、というミスを防ぐことができます。

<hr>

@section whatis 何が便利なの?

@subsection old_problem 従来の方法の欠点

enum 定数とそれに対応する文字列を定義したいとき、普通は以下のように書きます。色の番号と名前を定義する例です。

@code
  //////////////////// ヘッダファイル mycolor.h ////////////////////
  #ifndef MYCOLOR_H_INCLUDE
  #define MYCOLOR_H_INCLUDE
  enum MyColor { ///< 色を表す enum
    C_BLACK,
    C_BLUE,
    C_RED,
    C_GREEN,
    C_WHITE,
  };
  extern char * gMyColor[C_MAX]; ///< 色のenumに対応した文字列
  void DrawCircleWithColor(MyColor c, int x, int y, int r);///< 描画…など
  #endif //MYCOLOR_H_INCLUDE
@endcode

@code
  //////////////////// ソースファイル mycolor.cpp ////////////////////
  #include "mycolor.h"
  char * gMyColor[C_MAX] = {
    "黒", // C_BLACK
    "青", // C_BLUE
    "赤", // C_RED
    "緑", // C_GREEN
    "白", // C_WHITE
  };
@endcode

しかし、この方法には、enumを追加・削除・順序変更するたびに、対応する文字列も <B>忘れずに・間違えず</B> に修正しなければならない、という欠点があります。

rEnumStr.h を使うと、enum と文字列を同時に定義することができて、enum を追加・削除・順序変更したときに文字列を修正し忘れたり間違えたりすることを防げます。

<hr>

@subsection improved rEnumStr.h を使うとどうなるの?

先ほどの色の番号と名前の例は、以下のように書けます。

@code
  //////////////////// ヘッダファイル mycolor.h ////////////////////
  //---- 定型的なおまじない ----
  #ifndef MYCOLOR_DEFINED
  #define MYCOLOR_DEFINED
  #ifndef RENUMSTR_BODY
  #define RENUMSTR_BODY 0
  #endif
  #include "rEnumStr.h"
  //---- 色の番号と文字列を定義 ----
  RENUMSTR_BEGIN( MyColor, MyColorStr )
    RENUMSTR( C_BLACK,  "黒" )
    RENUMSTR( C_BLUE,   "青" )
    RENUMSTR( C_RED,    "赤" )
    RENUMSTR( C_GREEN,  "緑" )
    RENUMSTR( C_WHITE,  "白" )
  RENUMSTR_END()
  #endif //MYCOLOR_DEFINED
  //---- ここから普通のヘッダファイル ----
  #ifndef MYCOLOR_H_INCLUDE
  #define MYCOLOR_H_INCLUDE
  void DrawCircleWithColor(MyColor c, int x, int y, int r);///< 描画…など
  #endif //MYCOLOR_H_INCLUDE
@endcode

@code
  //////////////////// ソースファイル mycolor.cpp ////////////////////
  #include "mycolor.h"
  #undef MYCOLOR_DEFINED        // 文字列を取り込むために
  #define RENUMSTR_BODY 1       // 設定を変えて
  #include "mycolor.h"          // もう一度includeします。
@endcode



<hr>

@section usage 使い方

@subsection intro 使い始めのおまじない

@code
  //////////////////// ヘッダファイル mycolor.h ////////////////////
  //---- 定型的なおまじない ----
  #ifndef MYCOLOR_DEFINED
  #define MYCOLOR_DEFINED
  #ifndef RENUMSTR_BODY
  #define RENUMSTR_BODY 0
  #endif
  #include "rEnumStr.h"
    :
@endcode

複雑に見えますが、rEnumStr.h を使いたくなったときはこの部分をコピペして、 MYCOLOR_DEFINED マクロの名前を適当に置換するだけで大丈夫です。


@subsection main 定義している部分

@code
    :
  //---- 色の番号と文字列を定義 ----
  RENUMSTR_BEGIN( MyColor, MyColorStr )
    RENUMSTR( C_BLACK,  "黒" )
    RENUMSTR( C_BLUE,   "青" )
    RENUMSTR( C_RED,    "赤" )
    RENUMSTR( C_GREEN,  "緑" )
    RENUMSTR( C_WHITE,  "白" )
  RENUMSTR_END()
  #endif //MYCOLOR_DEFINED
@endcode

最初に RENUMSTR_BEGIN( enum型名 , 文字列関数名 ) で定義を開始し、最後に RENUMSTR_END() で定義を終わります。

その間に RENUMSTR( enum名 , 文字列 ) を並べて、enum と文字列を同時に定義します。



@subsection cpp 定義を実体化する部分

@code
  //////////////////// ソースファイル mycolor.cpp ////////////////////
  #include "mycolor.h"
  #undef MYCOLOR_DEFINED        // 文字列を取り込むために
  #define RENUMSTR_BODY 1       // 設定を変えて
  #include "mycolor.h"          // もう一度includeします。
@endcode

1回目の include で enum MyColor 型とその要素 C_BLACK, C_BLUE, ... が定義されます。1回目の include では RENUMSTR_BODY マクロが 0 になるため、 mycolor.h が以下のように展開されるからです。
@code
  ////////// 1回目はこのように展開される //////////
  enum MyColor {
    C_BLUE,
    C_RED,
    C_GREEN,
    C_WHITE,
  };
  const char * MyColorStr(int);
@endcode

2回目の include で、enumに対応する文字列を返す関数 MyColorStr() が作られます。2回目は mycolor.cpp でマクロ RENUMSTR_BODY を 1 に定義しているので、 mycolor.h を展開すると、以下のような関数 MyColorStr() が生成されます。

@code
  ////////// 2回目はこのように展開される //////////
  const char * MyColorStr(int i)
  {
    switch (i) {
    default: return "???";
    case C_BLUE: return "青";
    case C_RED: return "赤";
    case C_GREEN: return "緑";
    case C_WHITE: return "白";
    }
  }
@endcode

@subsection invoke enumや文字列を取り出すには

mycolor.h がマクロ展開されて enum 定義と文字列取得関数 MyColorStr() の宣言ができているので、これを include するだけで簡単に enum も文字列も使えるようになります。

@code
  //////////////////// メインルーチン main.cpp ////////////////////
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

文字列の実体定義 On/Off を enumごとに設定するために、
RENUMSTR_BODYFLAG マクロは「ローカル」になっています。
(rEnumStr.hをincludeするたびに \#undef RENUMSTR_BODYFLAG される)
よって RENUMSTR_BODYFLAG マクロは、<B>プロジェクト全体で設定するのではなく</B>、
\#include "rEnumStr.h" を使うたびにその直前で定義してください。
プロジェクト設定で RENUMSTR_BODYFLAG を定義すると予想外の挙動になります。

RENUMSTR_BODYFLAG マクロは必ず「ソースファイル」の \#include "rEnumStr.h" の
直前で定義してください。ヘッダファイルやプロジェクト設定で定義すると
コンパイルエラーになります。

*/
//---------------------------------------------------------------------------
//#ifndef rEnumStrH
//#define rEnumStrH



//---------------------------------------------------------------------------
#ifndef RENUMSTR_BODY //{
#error #include "rEnumStr.h" の前に RENUMSTR_BODY が定義されていません

// doxygen用コメントはここに書かないとダメ? {

/// enum文字列表の定義開始。
  /** @param nameEnum          enum型の名前
      @param nameStr           文字列取得関数の名前

      表を書き終わったら必ず RENUMSTR_END() を置いてください。
   */
#define RENUMSTR_BEGIN(nameEnum, nameStr)

  /** @brief enumと文字列をひとつ作る。
      @param idEnum            enum定数名
      @param idStr             enumに対応づける文字列
   */
#define RENUMSTR(idEnum, idStr)

  /** @brief enumをひとつ作り、それにenum名そのものの文字列を割り当てる。
      @param idEnum            enum定数名
   */
#define RENUMSTR1(idEnum)

  /** @brief enumに値を設定し、そのenumと文字列を対応させる。
      @param idEnum            enum定数名
      @param val               enumに設定する値
      @param idStr             enumに対応づける文字列
  */
#define RENUMSTR_LET(idEnum, val, idStr)

  /** @brief enumの中に自由な一行を書く。
      @param line              行の内容

      enum定数に別名をつけるなどに使います。
  */
#define RENUMSTR_ENUMONLY(line)

  /** @brief enum文字列表を終了する。

      表を書き終わったら必ず RENUMSTR_END() を置いてください。
   */
#define RENUMSTR_END()

  /** @brief [ユーザ定義] [必須] rEnumStr.hの動作モードを設定。

      rEnumStr.h を使うヘッダ/ソースファイルでは、
      #include "rEnumStr.h" の前に必ずこのマクロを定義してください。

      - 0 にすると、enum宣言と文字列取得関数のプロトタイプ宣言を生成します。
        (ヘッダモード)
      - 1 にすると、文字列取得関数の実体を生成します。
        (ソースモード)
   */
#define RENUMSTR_BODY

  /** @brief [ユーザ定義] 文字列実体の定義をOn/Offする。

      ソースファイル中の #include "rEnumStr.h" の直前でこのマクロを
      定義すると、文字列実体の定義を On/Off できます。
      リリース版では文字列を定義しない、などの用途に使います。

      - 0 にすると、文字列を生成しません。
      - 1 にすると、文字列を生成します。
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
/// 文字列の実体を定義するか? (ローカル)
/** 文字列がデバッグ用でありリリース版では定義しないなどの場合、
    リリース版のときは #define RENUMSTR_BODYFLAG 0 にすることで
    文字列を定義しないようにできます。
 */
#define RENUMSTR_BODYFLAG 1
#endif

#ifndef RENUMSTR_DEFAULTSTR
/// enum範囲外のときのdefault文字列 (ローカル)
/** 文字列取得関数に範囲外の引数が渡されたときに返す文字列を指定します。
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

// RENUMSTR_BODYFLAG を間違ってヘッダで定義したときはエラーにする
// プロジェクト設定で定義したときもエラーになっちゃうが…
#ifdef RENUMSTR_BODYFLAG
#error RENUMSTR_BODYFLAG は「ソースファイル」で定義してください
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
