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
 *    GF2 -- A small library for Factorization in Galois Field F2.
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
  //   _____   _     _ _   
  //  |___  | | |__ (_) |_ 
  //     / /  | '_ \| | __|
  //    / /   | |_) | | |_ 
  //   /_/    |_.__/|_|\__|
  */

  #define TT(N) B_Tables[N*(1<<7)+BITS7(tmp,N)]
  #define TTT   TT(0)^TT(1)^TT(2)^TT(3)^B_Tables[4*(1<<7)+(tmp>>28)]

  template <>
  void
  MM4RmultAssBlock<uint32_t,7>( uint32_t const * A, unsigned nRowsA,
                                uint32_t const * B_Tables,
                                uint32_t       * C ) {

    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; *C++ = TTT ; }    
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP

  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAddBlock<uint32_t,7>( uint32_t const * A, unsigned nRowsA,
                                uint32_t const * B_Tables,
                                uint32_t       * C ) {

    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; *C++ ^= TTT ; }    
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultRightBlock<uint32_t,7>( uint32_t       * A, unsigned nRowsA,
                                  uint32_t const * B_Tables ) {

    #define MM4R_SINGLE_OP { uint32_t tmp = *A ; *A++ = TTT ; }    
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }
  
  //////////////////////////////////////////////////////////////////////////////
  #undef  TTT
  #define TTT TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6)^TT(7)^B_Tables[8*(1<<7)+(tmp>>56)]
  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAssBlock<uint64_t,7>( uint64_t const * A, unsigned nRowsA,
                                uint64_t const * B_Tables,
                                uint64_t       * C ) {

    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ = TTT ; }
    MM4R_LOOP_UNROLLED_BY_8 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAddBlock<uint64_t,7>( uint64_t const * A, unsigned nRowsA,
                                uint64_t const * B_Tables,
                                uint64_t       * C ) {

    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ ^= TTT ; }
    MM4R_LOOP_UNROLLED_BY_8 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultRightBlock<uint64_t,7>( uint64_t       * A, unsigned nRowsA,
                                  uint64_t const * B_Tables ) {

    #define MM4R_SINGLE_OP { uint64_t tmp = *A ; *A++ = TTT ; }
    MM4R_LOOP_UNROLLED_BY_8 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////
  #undef  TTT
  #undef  TT
  #define TT(IDX) tmp = _mm_xor_si128( tmp, TB[IDX*(1<<7)+BITS7(tt,IDX)] )
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  tt = A->i64.lo ;                                                             \
  TB = (__m128i const *)B_Tables ;                                             \
  TT(0); TT(1); TT(2); TT(3); TT(4); TT(5); TT(6); TT(7);                      \
  TB += 8*(1<<7) ;                                                             \
  tmp = _mm_xor_si128( tmp, TB[tt>>56] ) ;                                     \
  tt = A->i64.hi ;                                                             \
  TB += (1<<8) ;                                                               \
  TT(0); TT(1); TT(2); TT(3); TT(4); TT(5); TT(6); TT(7);                      \
  TB += 8*(1<<7) ;                                                             \
  tmp = _mm_xor_si128( tmp, TB[tt>>56])

  template <>
  void
  MM4RmultAssBlock<uint128_t,7>( uint128_t const * A, unsigned nRowsA,
                                 uint128_t const * B_Tables,
                                 uint128_t       * C ) {

    uint64_t        tt ;
    __m128i const * TB ;

    #define MM4R_SINGLE_OP { __m128i tmp = _mm_setzero_si128() ; MM4R_SINGLE_OP_INTERNAL ; C++->i128 = tmp ; ++A ; }
    MM4R_LOOP_UNROLLED_BY_4 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultAddBlock<uint128_t,7>( uint128_t const * A, unsigned nRowsA,
                                 uint128_t const * B_Tables,
                                 uint128_t       * C ) {

    uint64_t        tt ;
    __m128i const * TB ;

    #define MM4R_SINGLE_OP { __m128i tmp = C->i128 ; MM4R_SINGLE_OP_INTERNAL ; C++->i128 = tmp ; ++A ; }
    MM4R_LOOP_UNROLLED_BY_4 ;
    #undef MM4R_SINGLE_OP
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4RmultRightBlock<uint128_t,7>( uint128_t       * A, unsigned nRowsA,
                                   uint128_t const * B_Tables ) {

    uint64_t        tt ;
    __m128i const * TB ;

    #define MM4R_SINGLE_OP { __m128i tmp = _mm_setzero_si128() ; MM4R_SINGLE_OP_INTERNAL ; A++->i128 = tmp ; }
    MM4R_LOOP_UNROLLED_BY_4 ;
    #undef MM4R_SINGLE_OP
  }

  #undef TT
  #undef MM4R_SINGLE_OP_INTERNAL
}

// EOF GF2toolkit_4RussianMult7bit.cc
