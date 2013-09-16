/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2010-2012                                                 |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                | 
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      Via Mesiano 77, I-38050 Trento, Italy                               |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
 |      version: 0.2                                                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

/*
 *    This file is part of GF2 Toolkit.
 *
 *    GF2 toolkit -- A small library for Factorization in Galois Field F2.
 *    Copyright (C) 2010-2012 by Enrico Bertolazzi.
 *    All rights reserved.
 *
 *    GF2 toolkit is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    GF2 toolkit is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADO Toolkit; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "GF2toolkit_4Russian.hh"

namespace GF2toolkit {

  /*
  //   __  __ __  __ _  _   ____                _           _____     _     _
  //  |  \/  |  \/  | || | |  _ \ _   _ ___ ___(_) __ _ _ _|_   _|_ _| |__ | | ___
  //  | |\/| | |\/| | || |_| |_) | | | / __/ __| |/ _` | '_ \| |/ _` | '_ \| |/ _ \
  //  | |  | | |  | |__   _|  _ <| |_| \__ \__ \ | (_| | | | | | (_| | |_) | |  __/
  //  |_|  |_|_|  |_|  |_| |_| \_\\__,_|___/___/_|\__,_|_| |_|_|\__,_|_.__/|_|\___|
  */

  #define MM4R_TABLE1  Table[kk^0x01] = Table[kk]^M[0x0] ; kk ^= 0x01 ;
  #define MM4R_TABLE2  MM4R_TABLE1  ; Table[kk^0x02]  = Table[kk]^M[0x1] ; kk ^= 0x02  ; MM4R_TABLE1
  #define MM4R_TABLE3  MM4R_TABLE2  ; Table[kk^0x04]  = Table[kk]^M[0x2] ; kk ^= 0x04  ; MM4R_TABLE2
  #define MM4R_TABLE4  MM4R_TABLE3  ; Table[kk^0x08]  = Table[kk]^M[0x3] ; kk ^= 0x08  ; MM4R_TABLE3
  #define MM4R_TABLE5  MM4R_TABLE4  ; Table[kk^0x10]  = Table[kk]^M[0x4] ; kk ^= 0x10  ; MM4R_TABLE4 
  #define MM4R_TABLE6  MM4R_TABLE5  ; Table[kk^0x20]  = Table[kk]^M[0x5] ; kk ^= 0x20  ; MM4R_TABLE5
  #define MM4R_TABLE7  MM4R_TABLE6  ; Table[kk^0x40]  = Table[kk]^M[0x6] ; kk ^= 0x40  ; MM4R_TABLE6
  #define MM4R_TABLE8  MM4R_TABLE7  ; Table[kk^0x80]  = Table[kk]^M[0x7] ; kk ^= 0x80  ; MM4R_TABLE7
  #define MM4R_TABLE9  MM4R_TABLE8  ; Table[kk^0x100] = Table[kk]^M[0x8] ; kk ^= 0x100 ; MM4R_TABLE8
  #define MM4R_TABLE10 MM4R_TABLE9  ; Table[kk^0x200] = Table[kk]^M[0x9] ; kk ^= 0x200 ; MM4R_TABLE9
  #define MM4R_TABLE11 MM4R_TABLE10 ; Table[kk^0x400] = Table[kk]^M[0xA] ; kk ^= 0x400 ; MM4R_TABLE10
  #define MM4R_TABLE12 MM4R_TABLE11 ; Table[kk^0x800] = Table[kk]^M[0xB] ; kk ^= 0x800 ; MM4R_TABLE11

  #define MAKETABLE_BODY(TYPE,N)                                    \
  template <> void MM4RmakeTable<N>( TYPE const M[], TYPE Table[] ) \
  { unsigned kk = 0 ; Table[0] = 0 ; MM4R_TABLE##N ; }

  MAKETABLE_BODY(uint16_t,1)
  MAKETABLE_BODY(uint16_t,2)
  MAKETABLE_BODY(uint16_t,3)
  MAKETABLE_BODY(uint16_t,4)
  MAKETABLE_BODY(uint16_t,5)
  MAKETABLE_BODY(uint16_t,6)
  MAKETABLE_BODY(uint16_t,7)
  MAKETABLE_BODY(uint16_t,8)
  MAKETABLE_BODY(uint16_t,9)
  MAKETABLE_BODY(uint16_t,10)
  MAKETABLE_BODY(uint16_t,11)
  MAKETABLE_BODY(uint16_t,12)

  MAKETABLE_BODY(uint32_t,1)
  MAKETABLE_BODY(uint32_t,2)
  MAKETABLE_BODY(uint32_t,3)
  MAKETABLE_BODY(uint32_t,4)
  MAKETABLE_BODY(uint32_t,5)
  MAKETABLE_BODY(uint32_t,6)
  MAKETABLE_BODY(uint32_t,7)
  MAKETABLE_BODY(uint32_t,8)
  MAKETABLE_BODY(uint32_t,9)
  MAKETABLE_BODY(uint32_t,10)
  MAKETABLE_BODY(uint32_t,11)
  MAKETABLE_BODY(uint32_t,12)

  MAKETABLE_BODY(uint64_t,1)
  MAKETABLE_BODY(uint64_t,2)
  MAKETABLE_BODY(uint64_t,3)
  MAKETABLE_BODY(uint64_t,4)
  MAKETABLE_BODY(uint64_t,5)
  MAKETABLE_BODY(uint64_t,6)
  MAKETABLE_BODY(uint64_t,7)
  MAKETABLE_BODY(uint64_t,8)
  MAKETABLE_BODY(uint64_t,9)
  MAKETABLE_BODY(uint64_t,10)
  MAKETABLE_BODY(uint64_t,11)
  MAKETABLE_BODY(uint64_t,12)

  #undef MM4R_TABLE1
  #undef MM4R_TABLE2
  #undef MM4R_TABLE3
  #undef MM4R_TABLE4
  #undef MM4R_TABLE5
  #undef MM4R_TABLE6
  #undef MM4R_TABLE7
  #undef MM4R_TABLE8
  #undef MM4R_TABLE9
  #undef MM4R_TABLE10
  #undef MM4R_TABLE11
  #undef MM4R_TABLE12
  #undef MAKETABLE_BODY

  #define MM4R_TABLE1  Table[kk^0x01].i128 = _mm_xor_si128( Table[kk].i128, M[0x0].i128 ) ; kk ^= 0x01 ;
  #define MM4R_TABLE2  MM4R_TABLE1  ; Table[kk^0x02].i128  = _mm_xor_si128( Table[kk].i128, M[0x1].i128 ) ; kk ^= 0x02  ; MM4R_TABLE1
  #define MM4R_TABLE3  MM4R_TABLE2  ; Table[kk^0x04].i128  = _mm_xor_si128( Table[kk].i128, M[0x2].i128 ) ; kk ^= 0x04  ; MM4R_TABLE2
  #define MM4R_TABLE4  MM4R_TABLE3  ; Table[kk^0x08].i128  = _mm_xor_si128( Table[kk].i128, M[0x3].i128 ) ; kk ^= 0x08  ; MM4R_TABLE3
  #define MM4R_TABLE5  MM4R_TABLE4  ; Table[kk^0x10].i128  = _mm_xor_si128( Table[kk].i128, M[0x4].i128 ) ; kk ^= 0x10  ; MM4R_TABLE4 
  #define MM4R_TABLE6  MM4R_TABLE5  ; Table[kk^0x20].i128  = _mm_xor_si128( Table[kk].i128, M[0x5].i128 ) ; kk ^= 0x20  ; MM4R_TABLE5
  #define MM4R_TABLE7  MM4R_TABLE6  ; Table[kk^0x40].i128  = _mm_xor_si128( Table[kk].i128, M[0x6].i128 ) ; kk ^= 0x40  ; MM4R_TABLE6
  #define MM4R_TABLE8  MM4R_TABLE7  ; Table[kk^0x80].i128  = _mm_xor_si128( Table[kk].i128, M[0x7].i128 ) ; kk ^= 0x80  ; MM4R_TABLE7
  #define MM4R_TABLE9  MM4R_TABLE8  ; Table[kk^0x100].i128 = _mm_xor_si128( Table[kk].i128, M[0x8].i128 ) ; kk ^= 0x100 ; MM4R_TABLE8
  #define MM4R_TABLE10 MM4R_TABLE9  ; Table[kk^0x200].i128 = _mm_xor_si128( Table[kk].i128, M[0x9].i128 ) ; kk ^= 0x200 ; MM4R_TABLE9
  #define MM4R_TABLE11 MM4R_TABLE10 ; Table[kk^0x400].i128 = _mm_xor_si128( Table[kk].i128, M[0xA].i128 ) ; kk ^= 0x400 ; MM4R_TABLE10
  #define MM4R_TABLE12 MM4R_TABLE11 ; Table[kk^0x800].i128 = _mm_xor_si128( Table[kk].i128, M[0xB].i128 ) ; kk ^= 0x800 ; MM4R_TABLE11

  #define MAKETABLE_BODY(TYPE,N)                                    \
  template <> void MM4RmakeTable<N>( TYPE const M[], TYPE Table[] ) \
  { unsigned kk = 0 ; Table[0].i128 = _mm_setzero_si128() ; MM4R_TABLE##N ; }

  MAKETABLE_BODY(uint128_t,1)
  MAKETABLE_BODY(uint128_t,2)
  MAKETABLE_BODY(uint128_t,3)
  MAKETABLE_BODY(uint128_t,4)
  MAKETABLE_BODY(uint128_t,5)
  MAKETABLE_BODY(uint128_t,6)
  MAKETABLE_BODY(uint128_t,7)
  MAKETABLE_BODY(uint128_t,8)
  MAKETABLE_BODY(uint128_t,9)
  MAKETABLE_BODY(uint128_t,10)
  MAKETABLE_BODY(uint128_t,11)
  MAKETABLE_BODY(uint128_t,12)
  
  #define MAKETABLE_REPEAT2( A ) A ; A
  #define MAKETABLE_REPEAT4( A ) A ; A ; A ; A
  #define MAKETABLE_REPEAT8( A ) A ; A ; A ; A ; A ; A ; A ; A
  #define MAKETABLE_REPEAT16( A ) MAKETABLE_REPEAT8( A ) ; MAKETABLE_REPEAT8( A )
  #define MAKETABLE_REPEAT32( A ) MAKETABLE_REPEAT16( A ) ; MAKETABLE_REPEAT16( A )

  /*
  //  ____    _     _ _   
  // |___ \  | |__ (_) |_ 
  //   __) | | '_ \| | __|
  //  / __/  | |_) | | |_ 
  // |_____| |_.__/|_|\__|
  */

  template <>
  void
  MM4RmakeTables<uint32_t,2>( uint32_t const M[], uint32_t * Table ) {
    MAKETABLE_REPEAT16( MM4RmakeTable<2>( M, Table ) ; M += 2 ; Table += 1<<2 ) ;
  }

  template <>
  void
  MM4RmakeTables<uint64_t,2>( uint64_t const M[], uint64_t * Table ) {
    MAKETABLE_REPEAT32( MM4RmakeTable<2>( M, Table ) ; M += 2 ; Table += 1<<2 ) ;
  }

  template <>
  void
  MM4RmakeTables<uint128_t,2>( uint128_t const M[], uint128_t * Table ) {
    MAKETABLE_REPEAT32( MM4RmakeTable<2>( M, Table ) ; M += 2 ; Table += 1<<2 ) ;
    MAKETABLE_REPEAT32( MM4RmakeTable<2>( M, Table ) ; M += 2 ; Table += 1<<2 ) ;
  }

  /*
  //   _____   _     _ _   
  //  |___ /  | |__ (_) |_ 
  //    |_ \  | '_ \| | __|
  //   ___) | | |_) | | |_ 
  //  |____/  |_.__/|_|\__|  
  */

  template <>
  void
  MM4RmakeTables<uint32_t,3>( uint32_t const M[], uint32_t * Table ) {
    MAKETABLE_REPEAT8( MM4RmakeTable<3>( M, Table ) ; M += 3 ; Table += 1<<3 ) ;
    MM4RmakeTable<4>( M, Table ) ; M += 4 ; Table += 1<<4 ;
    MM4RmakeTable<4>( M, Table ) ; // ultima due a 4 bit
  }

  template <>
  void
  MM4RmakeTables<uint64_t,3>( uint64_t const M[], uint64_t * Table ) {
    MAKETABLE_REPEAT16( MM4RmakeTable<3>( M, Table ) ; M += 3 ; Table += 1<<3 ) ;
    MAKETABLE_REPEAT4 ( MM4RmakeTable<3>( M, Table ) ; M += 3 ; Table += 1<<3 ) ;
    MM4RmakeTable<4>( M, Table ) ; // ultima a 4 bit
  }

  template <>
  void
  MM4RmakeTables<uint128_t,3>( uint128_t const M[], uint128_t * Table ) {
    MAKETABLE_REPEAT32( MM4RmakeTable<3>( M, Table ) ; M += 3 ; Table += 1<<3 ) ;
    MAKETABLE_REPEAT8 ( MM4RmakeTable<3>( M, Table ) ; M += 3 ; Table += 1<<3 ) ;
    MM4RmakeTable<4>( M, Table ) ; M += 4 ; Table += 1<<4 ;
    MM4RmakeTable<4>( M, Table ) ; // ultime due a 4 bit
  }

  /*
  //   _  _   _     _ _
  //  | || | | |__ (_) |_ 
  //  | || |_| '_ \| | __|
  //  |__   _| |_) | | |_ 
  //     |_| |_.__/|_|\__|
  */
  template <>
  void
  MM4RmakeTables<uint32_t,4>( uint32_t const M[], uint32_t * Table ) {
    MAKETABLE_REPEAT8( MM4RmakeTable<4>( M, Table ) ; M += 4 ; Table += 1<<4 ) ;
  }

  template <>
  void
  MM4RmakeTables<uint64_t,4>( uint64_t const M[], uint64_t * Table ) {
    MAKETABLE_REPEAT16( MM4RmakeTable<4>( M, Table ) ; M += 4 ; Table += 1<<4 ) ;
  }

  template <>
  void
  MM4RmakeTables<uint128_t,4>( uint128_t const M[], uint128_t * Table ) {
    MAKETABLE_REPEAT32( MM4RmakeTable<4>( M, Table ) ; M += 4 ; Table += 1<<4 ) ;
  }

  /*
  //  ____    _     _ _   
  // | ___|  | |__ (_) |_ 
  // |___ \  | '_ \| | __|
  //  ___) | | |_) | | |_ 
  // |____/  |_.__/|_|\__|
  */

  template <>
  void
  MM4RmakeTables<uint32_t,5>( uint32_t const M[], uint32_t * Table ) {
    MAKETABLE_REPEAT4( MM4RmakeTable<5>( M, Table ) ; M += 5 ; Table += 1<<5 ) ;
    MM4RmakeTable<5>( M, Table ) ; M += 5 ; Table += 1<<5 ;
    MM4RmakeTable<7>( M, Table ) ; // ultima tabella a 7 bit
  }

  template <>
  void
  MM4RmakeTables<uint64_t,5>( uint64_t const M[], uint64_t * Table ) {
    MAKETABLE_REPEAT8( MM4RmakeTable<5>( M, Table ) ; M += 5 ; Table += 1<<5 ) ;
    MAKETABLE_REPEAT4( MM4RmakeTable<5>( M, Table ) ; M += 5 ; Table += 1<<5 ) ;
    MM4RmakeTable<4>( M, Table ) ; // ultima tabella a 4 bit
  }

  template <>
  void
  MM4RmakeTables<uint128_t,5>( uint128_t const M[], uint128_t * Table ) {
    MAKETABLE_REPEAT16( MM4RmakeTable<5>( M, Table ) ; M += 5 ; Table += 1<<5 ) ;
    MAKETABLE_REPEAT8 ( MM4RmakeTable<5>( M, Table ) ; M += 5 ; Table += 1<<5 ) ;
    MM4RmakeTable<5>( M, Table ) ; M += 5 ; Table += 1<<5 ;
    MM4RmakeTable<3>( M, Table ) ; // ultima tabella a 3 bit
  }

  /*
  //    __     _     _ _   
  //   / /_   | |__ (_) |_ 
  //  | '_ \  | '_ \| | __|
  //  | (_) | | |_) | | |_ 
  //   \___/  |_.__/|_|\__|
  */
  template <>
  void
  MM4RmakeTables<uint32_t,6>( uint32_t const M[], uint32_t * Table ) {
    MM4RmakeTable<6>( M, Table ) ; M += 6 ; Table += 1<<6 ;
    MM4RmakeTable<6>( M, Table ) ; M += 6 ; Table += 1<<6 ;
    MM4RmakeTable<6>( M, Table ) ; M += 6 ; Table += 1<<6 ;
    MM4RmakeTable<7>( M, Table ) ; M += 7 ; Table += 1<<7 ;
    MM4RmakeTable<7>( M, Table ) ; // ultima 2 tabelle a 7 bit
  }

  template <>
  void
  MM4RmakeTables<uint64_t,6>( uint64_t const M[], uint64_t * Table ) {
    MAKETABLE_REPEAT8( MM4RmakeTable<6>( M, Table ) ; M += 6 ; Table += 1<<6 ) ;
    MM4RmakeTable<6>( M, Table ) ; M += 6 ; Table += 1<<6 ; // 1
    MM4RmakeTable<6>( M, Table ) ; M += 6 ; Table += 1<<6 ; // 2
    MM4RmakeTable<4>( M, Table ) ; // ultima tabella da 4 bit
  }

  template <>
  void
  MM4RmakeTables<uint128_t,6>( uint128_t const M[], uint128_t * Table ) {
    MAKETABLE_REPEAT16( MM4RmakeTable<6>( M, Table ) ; M += 6 ; Table += 1<<6 ) ;
    MM4RmakeTable<6>( M, Table ) ; M += 6 ; Table += 1<<6 ;
    MM4RmakeTable<6>( M, Table ) ; M += 6 ; Table += 1<<6 ;
    MM4RmakeTable<6>( M, Table ) ; M += 6 ; Table += 1<<6 ;
    MM4RmakeTable<7>( M, Table ) ; M += 7 ; Table += 1<<7 ;
    MM4RmakeTable<7>( M, Table ) ; // ultima 2 tabella da 7 bit
  }

  /*
  //   _____   _     _ _   
  //  |___  | | |__ (_) |_ 
  //     / /  | '_ \| | __|
  //    / /   | |_) | | |_ 
  //   /_/    |_.__/|_|\__|
  */
  template <>
  void
  MM4RmakeTables<uint32_t,7>( uint32_t const M[], uint32_t * Table ) {
    MAKETABLE_REPEAT4( MM4RmakeTable<7>( M, Table ) ; M += 7 ; Table += 1<<7 ) ;
    MM4RmakeTable<4>( M, Table ) ; // ultima tabella a 4 bit
  }

  template <>
  void
  MM4RmakeTables<uint64_t,7>( uint64_t const M[], uint64_t * Table ) {
    MAKETABLE_REPEAT8( MM4RmakeTable<7>( M, Table ) ; M += 7 ; Table += 1<<7 ) ;
    MM4RmakeTable<8>( M, Table ) ; // ultima tabella da 8 bit
  }

  template <>
  void
  MM4RmakeTables<uint128_t,7>( uint128_t const M[], uint128_t * Table ) {
    MAKETABLE_REPEAT8( MM4RmakeTable<7>( M, Table ) ; M += 7 ; Table += 1<<7 ) ;
    MM4RmakeTable<8>( M, Table ) ; M += 8 ; Table += 1<<8 ; // ultime tabella da 8 bit
    MAKETABLE_REPEAT8( MM4RmakeTable<7>( M, Table ) ; M += 7 ; Table += 1<<7 ) ;
    MM4RmakeTable<8>( M, Table ) ; // ultime tabella da 8 bit
  }

  /*
  //    ___  _     _ _   
  //   ( _ )| |__ (_) |_ 
  //   / _ \| '_ \| | __|
  //  | (_) | |_) | | |_ 
  //   \___/|_.__/|_|\__|
  */
  template <>
  void
  MM4RmakeTables<uint32_t,8>( uint32_t const M[], uint32_t * Table ) {
    MM4RmakeTable<8>( M, Table ) ; M += 8 ; Table += 1<<8 ;
    MM4RmakeTable<8>( M, Table ) ; M += 8 ; Table += 1<<8 ;
    MM4RmakeTable<8>( M, Table ) ; M += 8 ; Table += 1<<8 ;
    MM4RmakeTable<8>( M, Table ) ;
  }

  template <>
  void
  MM4RmakeTables<uint64_t,8>( uint64_t const M[], uint64_t * Table ) {
    MAKETABLE_REPEAT8( MM4RmakeTable<8>( M, Table ) ; M += 8 ; Table += 1<<8 ) ;
  }

  template <>
  void
  MM4RmakeTables<uint128_t,8>( uint128_t const M[], uint128_t * Table ) {
    MAKETABLE_REPEAT16( MM4RmakeTable<8>( M, Table ) ; M += 8 ; Table += 1<<8 ) ;
  }

  /*
  //    ___    _     _ _   
  //   / _ \  | |__ (_) |_ 
  //  | (_) | | '_ \| | __|
  //   \__, | | |_) | | |_ 
  //     /_/  |_.__/|_|\__|
  */

  template <>
  void
  MM4RmakeTables<uint32_t,9>( uint32_t const M[], uint32_t * Table ) {
    MM4RmakeTable<9>( M, Table ) ; M += 9 ; Table += 1<<9 ;
    MM4RmakeTable<9>( M, Table ) ; M += 9 ; Table += 1<<9 ;
    MM4RmakeTable<9>( M, Table ) ; M += 9 ; Table += 1<<9 ;
    MM4RmakeTable<5>( M, Table ) ; // ultima tabella 9 bit
  }

  template <>
  void
  MM4RmakeTables<uint64_t,9>( uint64_t const M[], uint64_t * Table ) {
    MAKETABLE_REPEAT4( MM4RmakeTable<9>( M, Table ) ; M += 9 ; Table += 1<<9 ) ;
    MM4RmakeTable<9>( M, Table ) ; M += 9 ; Table += 1<<9 ;
    MM4RmakeTable<9>( M, Table ) ; M += 9 ; Table += 1<<9 ;
    MM4RmakeTable<10>( M, Table ) ; // ultima tabella 10 bit
  }

  template <>
  void
  MM4RmakeTables<uint128_t,9>( uint128_t const M[], uint128_t * Table ) {
    MAKETABLE_REPEAT4( MM4RmakeTable<9>( M, Table ) ; M += 9 ; Table += 1<<9 ) ;
    MM4RmakeTable<9> ( M, Table ) ; M += 9  ; Table += 1<<9  ;
    MM4RmakeTable<9> ( M, Table ) ; M += 9  ; Table += 1<<9  ;
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ; // ultima tabella da 10 bit
    MAKETABLE_REPEAT4( MM4RmakeTable<9>( M, Table ) ; M += 9 ; Table += 1<<9 ) ;
    MM4RmakeTable<9> ( M, Table ) ; M += 9  ; Table += 1<<9  ;
    MM4RmakeTable<9> ( M, Table ) ; M += 9  ; Table += 1<<9  ;
    MM4RmakeTable<10>( M, Table ) ; // ultima tabella da 10 bit
  }

  /*
  //   _  ___    _     _ _   
  //  / |/ _ \  | |__ (_) |_ 
  //  | | | | | | '_ \| | __|
  //  | | |_| | | |_) | | |_ 
  //  |_|\___/  |_.__/|_|\__|
  */                   

  template <>
  void
  MM4RmakeTables<uint32_t,10>( uint32_t const M[], uint32_t * Table ) {
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ;
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ;
    MM4RmakeTable<12>( M, Table ) ;
  }

  template <>
  void
  MM4RmakeTables<uint64_t,10>( uint64_t const M[], uint64_t * Table ) {
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ;
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ;
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ;
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ;
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ;
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ;
    MM4RmakeTable<4>( M, Table ) ;
  }

  template <>
  void
  MM4RmakeTables<uint128_t,10>( uint128_t const M[], uint128_t * Table ) {
    MAKETABLE_REPEAT8( MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ) ;
    MAKETABLE_REPEAT4( MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ) ;
    MM4RmakeTable<8>( M, Table ) ;
  }

  /*
  //   _ _   _     _ _   
  //  / / | | |__ (_) |_ 
  //  | | | | '_ \| | __|
  //  | | | | |_) | | |_ 
  //  |_|_| |_.__/|_|\__|
  */

  template <>
  void
  MM4RmakeTables<uint32_t,11>( uint32_t const M[], uint32_t * Table ) {
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<10>( M, Table ) ;
  }

  template <>
  void
  MM4RmakeTables<uint64_t,11>( uint64_t const M[], uint64_t * Table ) {
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ;
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<10>( M, Table ) ;
  }

  template <>
  void
  MM4RmakeTables<uint128_t,11>( uint128_t const M[], uint128_t * Table ) {
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ; //-----
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ; //-----
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<10>( M, Table ) ; M += 10 ; Table += 1<<10 ; //-----
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<11>( M, Table ) ; M += 11 ; Table += 1<<11 ;
    MM4RmakeTable<10>( M, Table ) ; //-----
  }

}

// EOF GF2toolkit_4RussianTables.cc
