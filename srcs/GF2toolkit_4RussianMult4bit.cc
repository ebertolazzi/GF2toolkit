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
#include "GF2toolkit_4RussianMacros.hh"

namespace GF2toolkit {

  //////////////////////////////////////////////////////////////////////////////
  /*
  //  +---+   +---+ +---+    
  //  | C | = | A | | B |    
  //  |   |   |   | +---+
  //  |   |   |   |
  //  +---+   +---+
  */

  /*
  //   _  _   _     _ _   
  //  | || | | |__ (_) |_ 
  //  | || |_| '_ \| | __|
  //  |__   _| |_) | | |_ 
  //     |_| |_.__/|_|\__|
  */

  #define TT(N) B_Tables[N*(1<<4)+BITS4(tmp,N)]
  #define TTT   TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6)^B_Tables[7*16+(tmp>>28)]

  template <>
  void
  MM4RmultAssBlock<uint32_t,4>( uint32_t const * A, unsigned nRowsA,
                                uint32_t const * B_Tables,
                                uint32_t       * C ) {

    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; if ( tmp ) *C = TTT ; else *C = 0 ; ++C ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAddBlock<uint32_t,4>( uint32_t const * A, unsigned nRowsA,
                                uint32_t const * B_Tables,
                                uint32_t       * C ) {

    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; if ( tmp ) *C ^= TTT ; ++C ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultRightBlock<uint32_t,4>( uint32_t       * A, unsigned nRowsA,
                                  uint32_t const * B_Tables ) {

    #define MM4R_SINGLE_OP { uint32_t tmp = *A ; if ( tmp ) *A = TTT ; ++A ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////
  #undef  TTT
  #define TTT TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6)^TT(7)^ \
              TT(8)^TT(9)^TT(10)^TT(11)^TT(12)^TT(13)^TT(14)^ \
              B_Tables[15*(1<<4)+(tmp>>60)]
  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAssBlock<uint64_t,4>( uint64_t const * A, unsigned nRowsA,
                                uint64_t const * B_Tables,
                                uint64_t       * C ) {

    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; if ( tmp ) *C = TTT ; else *C = 0 ; ++C ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAddBlock<uint64_t,4>( uint64_t const * A, unsigned nRowsA,
                                uint64_t const * B_Tables,
                                uint64_t       * C ) {

    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; if  (tmp ) *C ^= TTT ; ++C ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultRightBlock<uint64_t,4>( uint64_t       * A, unsigned nRowsA,
                                  uint64_t const * B_Tables ) {

    #define MM4R_SINGLE_OP { uint64_t tmp = *A ; if ( tmp ) *A = TTT ; ++A ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }
  

  //////////////////////////////////////////////////////////////////////////////
  #undef  TT
  #undef  TTT
  #define TT(IDX1,IDX2) tmp = _mm_xor_si128( tmp, B_Tables[IDX1*(1<<4)+BITS4(tt,IDX2)].i128 )
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) {                                                    \
    TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ;                                    \
    TT(4,4) ; TT(5,5) ; TT(6,6) ; TT(7,7) ;                                    \
  }                                                                            \
  if ( (tt = A->i32.n1) ) {                                                    \
    TT(8,0)  ; TT(9,1)  ; TT(10,2) ; TT(11,3) ;                                \
    TT(12,4) ; TT(13,5) ; TT(14,6) ; TT(15,7) ;                                \
  }                                                                            \
  if ( (tt = A->i32.n2) ) {                                                    \
    TT(16,0) ; TT(17,1) ; TT(18,2) ; TT(19,3) ;                                \
    TT(20,4) ; TT(21,5) ; TT(22,6) ; TT(23,7) ;                                \
  }                                                                            \
  if ( (tt = A->i32.n3) ) {                                                    \
    TT(24,0) ; TT(25,1) ; TT(26,2) ; TT(27,3) ;                                \
    TT(28,4) ; TT(29,5) ; TT(30,6) ; TT(31,7) ;                                \
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAssBlock<uint128_t,4>( uint128_t const * A, unsigned nRowsA,
                                 uint128_t const * B_Tables,
                                 uint128_t       * C ) {

    uint64_t tt ;
    #define MM4R_SINGLE_OP { __m128i tmp = _mm_setzero_si128() ; MM4R_SINGLE_OP_INTERNAL ; C++->i128 = tmp ; ++A ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAddBlock<uint128_t,4>( uint128_t const * A, unsigned nRowsA,
                                 uint128_t const * B_Tables,
                                 uint128_t       * C ) {

    uint64_t tt ;
    #define MM4R_SINGLE_OP { __m128i tmp = C->i128 ; MM4R_SINGLE_OP_INTERNAL ; C++->i128 = tmp ; ++A ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultRightBlock<uint128_t,4>( uint128_t       * A, unsigned nRowsA,
                                   uint128_t const * B_Tables ) {

    uint64_t tt ;
    #define MM4R_SINGLE_OP { __m128i tmp = _mm_setzero_si128() ; MM4R_SINGLE_OP_INTERNAL ; A++->i128 = tmp ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  #undef TT
  #undef MM4R_SINGLE_OP_INTERNAL
}

// EOF GF2toolkit_4RussianMult4bit.cc
