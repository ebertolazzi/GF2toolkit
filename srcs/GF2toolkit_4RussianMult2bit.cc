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
  //  ____    _     _ _   
  // |___ \  | |__ (_) |_ 
  //   __) | | '_ \| | __|
  //  / __/  | |_) | | |_ 
  // |_____| |_.__/|_|\__|
  */

  #define TT(N) B_Tables[N*4+BITS2(tmp,N)]
  #define TTT   TT(0)^TT(1)^TT(2)^ TT(3)^ TT(4)^ TT(5)^ TT(6)^ TT(7)^ \
                TT(8)^TT(9)^TT(10)^TT(11)^TT(12)^TT(13)^TT(14)^TT(15)

  template <>
  void
  MM4RmultAssBlock<uint32_t,2>( uint32_t const * A, unsigned nRowsA,
                                uint32_t const * B_Tables,
                                uint32_t       * C ) {

    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; *C++ = TTT ; }    
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP

  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAddBlock<uint32_t,2>( uint32_t const * A, unsigned nRowsA,
                                uint32_t const * B_Tables,
                                uint32_t       * C ) {

    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; *C++ ^= TTT ; }    
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultRightBlock<uint32_t,2>( uint32_t       * A, unsigned nRowsA,
                                  uint32_t const * B_Tables ) {

    #define MM4R_SINGLE_OP { uint32_t tmp = *A ; *A++ = TTT ; }    
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////
  #undef  TTT
  #define TTT TT(0) ^TT(1) ^TT(2) ^TT(3) ^TT(4) ^TT(5) ^TT(6) ^TT(7) ^  \
              TT(8) ^TT(9) ^TT(10)^TT(11)^TT(12)^TT(13)^TT(14)^TT(15)^  \
              TT(16)^TT(17)^TT(18)^TT(19)^TT(20)^TT(21)^TT(22)^TT(23)^  \
              TT(24)^TT(25)^TT(26)^TT(27)^TT(28)^TT(29)^TT(30)^TT(31)
  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAssBlock<uint64_t,2>( uint64_t const * A, unsigned nRowsA,
                                uint64_t const * B_Tables,
                                uint64_t       * C ) {

    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ = TTT ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAddBlock<uint64_t,2>( uint64_t const * A, unsigned nRowsA,
                                uint64_t const * B_Tables,
                                uint64_t       * C ) {

    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ ^= TTT ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultRightBlock<uint64_t,2>( uint64_t       * A, unsigned nRowsA,
                                  uint64_t const * B_Tables ) {

    #define MM4R_SINGLE_OP { uint64_t tmp = *A ; *A++ = TTT ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////
  
  #undef TTT
  #undef TT

  #define TT(IDX) tmp = _mm_xor_si128( tmp, TB[IDX*4+BITS2(tt,IDX)] )
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  tt = A->i32.n0 ;                                                             \
  TB = (__m128i const *)B_Tables ;                                             \
  TT(0) ; TT(1)  ; TT(2)  ; TT(3)  ; TT(4)  ; TT(5)  ; TT(6)  ; TT(7)  ;       \
  TT(8) ; TT(9)  ; TT(10) ; TT(11) ; TT(12) ; TT(13) ; TT(14) ; TT(15) ;       \
  tt = A->i32.n1 ;                                                             \
  TB += 16*4 ;                                                                 \
  TT(0) ; TT(1)  ; TT(2)  ; TT(3)  ; TT(4)  ; TT(5)  ; TT(6)  ; TT(7)  ;       \
  TT(8) ; TT(9)  ; TT(10) ; TT(11) ; TT(12) ; TT(13) ; TT(14) ; TT(15) ;       \
  tt = A->i32.n2 ;                                                             \
  TB += 16*4 ;                                                                 \
  TT(0) ; TT(1)  ; TT(2)  ; TT(3)  ; TT(4)  ; TT(5)  ; TT(6)  ; TT(7)  ;       \
  TT(8) ; TT(9)  ; TT(10) ; TT(11) ; TT(12) ; TT(13) ; TT(14) ; TT(15) ;       \
  tt = A->i32.n3 ;                                                             \
  TB += 16*4 ;                                                                 \
  TT(0) ; TT(1)  ; TT(2)  ; TT(3)  ; TT(4)  ; TT(5)  ; TT(6)  ; TT(7)  ;       \
  TT(8) ; TT(9)  ; TT(10) ; TT(11) ; TT(12) ; TT(13) ; TT(14) ; TT(15) ; 

  template <>
  void
  MM4RmultAssBlock<uint128_t,2>( uint128_t const * A, unsigned nRowsA,
                                 uint128_t const * B_Tables,
                                 uint128_t       * C ) {

    uint32_t        tt ;
    __m128i const * TB ;

    #define MM4R_SINGLE_OP { __m128i tmp = _mm_setzero_si128() ; MM4R_SINGLE_OP_INTERNAL ; C++->i128 = tmp ; ++A ; }
    MM4R_LOOP_UNROLLED_BY_4 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAddBlock<uint128_t,2>( uint128_t const * A, unsigned nRowsA,
                                 uint128_t const * B_Tables,
                                 uint128_t       * C ) {

    uint32_t        tt ;
    __m128i const * TB ;

    #define MM4R_SINGLE_OP { __m128i tmp = C->i128 ; MM4R_SINGLE_OP_INTERNAL ; C++->i128 = tmp ; ++A ; }
    MM4R_LOOP_UNROLLED_BY_4 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultRightBlock<uint128_t,2>( uint128_t       * A, unsigned nRowsA,
                                   uint128_t const * B_Tables ) {

    uint32_t        tt ;
    __m128i const * TB ;

    #define MM4R_SINGLE_OP { __m128i tmp = _mm_setzero_si128() ; MM4R_SINGLE_OP_INTERNAL ; A++->i128 = tmp ; }
    MM4R_LOOP_UNROLLED_BY_4 ;
    #undef MM4R_SINGLE_OP
  }

  #undef TT
  #undef MM4R_SINGLE_OP_INTERNAL
}

// EOF GF2_4RussianMult2bit.cc
