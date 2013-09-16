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
 |      Dipartimento di Ingegneria Industriale           e                  |
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

  #define TT(N) B_Tables[N*(1<<8)+BITS8(tmp,N)]

  static
  void
  MM4RmultRightBlock1( uint32_t * A, unsigned nRowsA, uint32_t const * B_Tables, unsigned rshift ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = (*A)>>rshift ; *A++ = TT(0) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }
  
  static
  void
  MM4RmultRightBlock2( uint32_t * A, unsigned nRowsA, uint32_t const * B_Tables, unsigned rshift ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = (*A)>>rshift ; *A++ = TT(0)^TT(1) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock3( uint32_t * A, unsigned nRowsA, uint32_t const * B_Tables, unsigned rshift ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = (*A)>>rshift ; *A++ = TT(0)^TT(1)^TT(2) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock4( uint32_t * A, unsigned nRowsA, uint32_t const * B_Tables, unsigned rshift ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = (*A)>>rshift ; *A++ = TT(0)^TT(1)^TT(2)^TT(3) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  template <>
  void
  MM4Rreduced<uint32_t>::multRightShift( uint32_t * A, unsigned nrows, unsigned rshift ) const {
    switch ( nblock ) {
      case 1: MM4RmultRightBlock1(A,nrows,Tables,rshift) ; break ;
      case 2: MM4RmultRightBlock2(A,nrows,Tables,rshift) ; break ;
      case 3: MM4RmultRightBlock3(A,nrows,Tables,rshift) ; break ;
      case 4: MM4RmultRightBlock4(A,nrows,Tables,rshift) ; break ;
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////

  static
  void
  MM4RmultAssBlock1( uint32_t const * A, unsigned nRowsA,
                     uint32_t const * B_Tables,
                     uint32_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = (*A++)>>rshift ; *C++ = TT(0) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }
  
  static
  void
  MM4RmultAssBlock2( uint32_t const * A, unsigned nRowsA,
                     uint32_t const * B_Tables,
                     uint32_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = (*A++)>>rshift ; *C++ = TT(0)^TT(1) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock3( uint32_t const * A, unsigned nRowsA,
                     uint32_t const * B_Tables,
                     uint32_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = (*A++)>>rshift ; *C++ = TT(0)^TT(1)^TT(2) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock4( uint32_t const * A, unsigned nRowsA,
                     uint32_t const * B_Tables,
                     uint32_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = (*A++)>>rshift ; *C++ = TT(0)^TT(1)^TT(2)^TT(3) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  template <>
  void
  MM4Rreduced<uint32_t>::multRightAssShift( uint32_t const * A, unsigned nrows, 
                                            uint32_t       * C,
                                            unsigned         rshift ) const {
    switch ( nblock ) {
      case 1: MM4RmultAssBlock1(A,nrows,Tables,C,rshift) ; break ;
      case 2: MM4RmultAssBlock2(A,nrows,Tables,C,rshift) ; break ;
      case 3: MM4RmultAssBlock3(A,nrows,Tables,C,rshift) ; break ;
      case 4: MM4RmultAssBlock4(A,nrows,Tables,C,rshift) ; break ;
    }  
  }
  
  //////////////////////////////////////////////////////////////////////////////

  static
  void
  MM4RmultAddBlock1( uint32_t const * A, unsigned nRowsA, 
                     uint32_t const * B_Tables,
                     uint32_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = (*A++)>>rshift ; *C++ ^= TT(0) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }
  
  static
  void
  MM4RmultAddBlock2( uint32_t const * A, unsigned nRowsA, 
                     uint32_t const * B_Tables,
                     uint32_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = (*A++)>>rshift ; *C++ ^= TT(0)^TT(1) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock3( uint32_t const * A, unsigned nRowsA, 
                     uint32_t const * B_Tables,
                     uint32_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = (*A++)>>rshift ; *C++ ^= TT(0)^TT(1)^TT(2) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock4( uint32_t const * A, unsigned nRowsA, 
                     uint32_t const * B_Tables,
                     uint32_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = (*A++)>>rshift ; *C++ ^= TT(0)^TT(1)^TT(2)^TT(3) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  template <>
  void
  MM4Rreduced<uint32_t>::multRightAddShift( uint32_t const * A, unsigned nrows, 
                                            uint32_t       * C,
                                            unsigned         rshift  ) const {
    switch ( nblock ) {
      case 1: MM4RmultAddBlock1(A,nrows,Tables,C,rshift) ; break ;
      case 2: MM4RmultAddBlock2(A,nrows,Tables,C,rshift) ; break ;
      case 3: MM4RmultAddBlock3(A,nrows,Tables,C,rshift) ; break ;
      case 4: MM4RmultAddBlock4(A,nrows,Tables,C,rshift) ; break ;
    }  
  }

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  static
  void
  MM4RmultRightBlock1( uint64_t       * A, unsigned nRowsA, 
                       uint64_t const * B_Tables,
                       unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A)>>rshift ; *A++ = TT(0) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock2( uint64_t       * A, unsigned nRowsA, 
                       uint64_t const * B_Tables,
                       unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A)>>rshift ; *A++ = TT(0)^TT(1) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock3( uint64_t       * A, unsigned nRowsA, 
                       uint64_t const * B_Tables,
                       unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A)>>rshift ; *A++ = TT(0)^TT(1)^TT(2) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock4( uint64_t       * A, unsigned nRowsA, 
                       uint64_t const * B_Tables,
                       unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A)>>rshift ; *A++ = TT(0)^TT(1)^TT(2)^TT(3) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock5( uint64_t       * A, unsigned nRowsA, 
                       uint64_t const * B_Tables,
                       unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A)>>rshift ; *A++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock6( uint64_t       * A, unsigned nRowsA, 
                       uint64_t const * B_Tables,
                       unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A)>>rshift ; *A++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock7( uint64_t       * A, unsigned nRowsA, 
                       uint64_t const * B_Tables,
                       unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A)>>rshift ; *A++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock8( uint64_t       * A, unsigned nRowsA, 
                       uint64_t const * B_Tables,
                       unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A)>>rshift ; *A++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6)^TT(7) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  template <>
  void
  MM4Rreduced<uint64_t>::multRightShift( uint64_t * A, unsigned nrows, unsigned rshift ) const {
    switch ( nblock ) {
      case 1: MM4RmultRightBlock1(A,nrows,Tables,rshift) ; break ;
      case 2: MM4RmultRightBlock2(A,nrows,Tables,rshift) ; break ;
      case 3: MM4RmultRightBlock3(A,nrows,Tables,rshift) ; break ;
      case 4: MM4RmultRightBlock4(A,nrows,Tables,rshift) ; break ;
      case 5: MM4RmultRightBlock5(A,nrows,Tables,rshift) ; break ;
      case 6: MM4RmultRightBlock6(A,nrows,Tables,rshift) ; break ;
      case 7: MM4RmultRightBlock7(A,nrows,Tables,rshift) ; break ;
      case 8: MM4RmultRightBlock8(A,nrows,Tables,rshift) ; break ;
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////

  static
  void
  MM4RmultAssBlock1( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ = TT(0) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }
  
  static
  void
  MM4RmultAssBlock2( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ = TT(0)^TT(1) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock3( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ = TT(0)^TT(1)^TT(2) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock4( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ = TT(0)^TT(1)^TT(2)^TT(3) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock5( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock6( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock7( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock8( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6)^TT(7) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  template <>
  void
  MM4Rreduced<uint64_t>::multRightAssShift( uint64_t const * A, unsigned nrows,
                                            uint64_t       * C,
                                            unsigned         rshift ) const {
    switch ( nblock ) {
      case 1: MM4RmultAssBlock1(A,nrows,Tables,C,rshift) ; break ;
      case 2: MM4RmultAssBlock2(A,nrows,Tables,C,rshift) ; break ;
      case 3: MM4RmultAssBlock3(A,nrows,Tables,C,rshift) ; break ;
      case 4: MM4RmultAssBlock4(A,nrows,Tables,C,rshift) ; break ;
      case 5: MM4RmultAssBlock5(A,nrows,Tables,C,rshift) ; break ;
      case 6: MM4RmultAssBlock6(A,nrows,Tables,C,rshift) ; break ;
      case 7: MM4RmultAssBlock7(A,nrows,Tables,C,rshift) ; break ;
      case 8: MM4RmultAssBlock8(A,nrows,Tables,C,rshift) ; break ;
    }  
  }
  
  //////////////////////////////////////////////////////////////////////////////

  static
  void
  MM4RmultAddBlock1( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ ^= TT(0) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }
  
  static
  void
  MM4RmultAddBlock2( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ ^= TT(0)^TT(1) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock3( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ ^= TT(0)^TT(1)^TT(2) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock4( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ ^= TT(0)^TT(1)^TT(2)^TT(3) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock5( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ ^= TT(0)^TT(1)^TT(2)^TT(3)^TT(4) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock6( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ ^= TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock7( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ ^= TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock8( uint64_t const * A, unsigned nRowsA,
                     uint64_t const * B_Tables, 
                     uint64_t       * C,
                     unsigned         rshift ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = (*A++)>>rshift ; *C++ ^= TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6)^TT(7) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  template <>
  void
  MM4Rreduced<uint64_t>::multRightAddShift( uint64_t const * A, unsigned nrows, 
                                            uint64_t       * C,
                                            unsigned         rshift ) const {
    switch ( nblock ) {
      case 1: MM4RmultAddBlock1(A,nrows,Tables,C,rshift) ; break ;
      case 2: MM4RmultAddBlock2(A,nrows,Tables,C,rshift) ; break ;
      case 3: MM4RmultAddBlock3(A,nrows,Tables,C,rshift) ; break ;
      case 4: MM4RmultAddBlock4(A,nrows,Tables,C,rshift) ; break ;
      case 5: MM4RmultAddBlock5(A,nrows,Tables,C,rshift) ; break ;
      case 6: MM4RmultAddBlock6(A,nrows,Tables,C,rshift) ; break ;
      case 7: MM4RmultAddBlock7(A,nrows,Tables,C,rshift) ; break ;
      case 8: MM4RmultAddBlock8(A,nrows,Tables,C,rshift) ; break ;
    }  
  }

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  #undef  TT
  #define TT(IDX1,IDX2) tmp = _mm_xor_si128( tmp, B_Tables[IDX1*(1<<8)+BITS8(tt,IDX2)].i128 )

  #define MULT_BLOCK(N)                                                        \
  static                                                                       \
  void                                                                         \
  MM4RmultRightBlock##N( uint128_t       * A, unsigned nRowsA,                 \
                         uint128_t const * B_Tables,                           \
                         unsigned          rshift ) {                          \
    uint32_t  tt ;                                                             \
    uint64_t  t64 ;                                                            \
    uint128_t t128 ;                                                           \
    if ( rshift == 64 ) {                                                      \
      GF2_UNROLL_BY_16( nRowsA, { __m128i tmp = _mm_setzero_si128() ;          \
                                  MM4R_SINGLE_OP_INTERNAL_EQ_64 ;              \
                                  A++->i128 = tmp ; } ) ;                      \
    } else if ( rshift > 64 ) {                                                \
      unsigned rs = rshift-64 ;                                                \
      GF2_UNROLL_BY_16( nRowsA, { __m128i tmp = _mm_setzero_si128() ;          \
                                  MM4R_SINGLE_OP_INTERNAL_GT_64 ;              \
                                  A++->i128 = tmp ; } ) ;                      \
    } else {                                                                   \
      unsigned lshift = 64-rshift ;                                            \
      GF2_UNROLL_BY_16( nRowsA, { __m128i tmp = _mm_setzero_si128() ;          \
                                  MM4R_SINGLE_OP_INTERNAL_LT_64 ;              \
                                  A++->i128 = tmp ; } ) ;                      \
    }                                                                          \
  }                                                                            \
                                                                               \
  static                                                                       \
  void                                                                         \
  MM4RmultAssBlock##N( uint128_t const * A, unsigned nRowsA,                   \
                       uint128_t const * B_Tables,                             \
                       uint128_t       * C,                                    \
                       unsigned          rshift ) {                            \
    uint32_t  tt ;                                                             \
    uint64_t  t64 ;                                                            \
    uint128_t t128 ;                                                           \
    if ( rshift == 64 ) {                                                      \
      GF2_UNROLL_BY_16( nRowsA, { __m128i tmp = _mm_setzero_si128() ;          \
                                  MM4R_SINGLE_OP_INTERNAL_EQ_64 ;              \
                                  C++->i128 = tmp ; ++A ; } ) ;                \
    } else if ( rshift > 64 ) {                                                \
      unsigned rs = rshift-64 ;                                                \
      GF2_UNROLL_BY_16( nRowsA, { __m128i tmp = _mm_setzero_si128() ;          \
                                  MM4R_SINGLE_OP_INTERNAL_GT_64 ;              \
                                  C++->i128 = tmp ; ++A ; } ) ;                \
    } else {                                                                   \
      unsigned lshift = 64-rshift ;                                            \
      GF2_UNROLL_BY_16( nRowsA, { __m128i tmp = _mm_setzero_si128() ;          \
                                  MM4R_SINGLE_OP_INTERNAL_LT_64 ;              \
                                  C++->i128 = tmp ; ++A ; } ) ;                \
    }                                                                          \
  }                                                                            \
                                                                               \
  static                                                                       \
  void                                                                         \
  MM4RmultAddBlock##N( uint128_t const * A, unsigned nRowsA,                   \
                       uint128_t const * B_Tables,                             \
                       uint128_t       * C,                                    \
                       unsigned          rshift ) {                            \
    uint32_t  tt ;                                                             \
    uint64_t  t64 ;                                                            \
    uint128_t t128 ;                                                           \
    if ( rshift == 64 ) {                                                      \
      GF2_UNROLL_BY_16( nRowsA, { __m128i tmp = C->i128 ;                      \
                                  MM4R_SINGLE_OP_INTERNAL_EQ_64 ;              \
                                  C++->i128 = tmp ; ++A ; } ) ;                \
    } else if ( rshift > 64 ) {                                                \
      unsigned rs = rshift-64 ;                                                \
      GF2_UNROLL_BY_16( nRowsA, { __m128i tmp = C->i128 ;                      \
                                  MM4R_SINGLE_OP_INTERNAL_GT_64 ;              \
                                  C++->i128 = tmp ; ++A ; } ) ;                \
    } else {                                                                   \
      unsigned lshift = 64-rshift ;                                            \
      GF2_UNROLL_BY_16( nRowsA, { __m128i tmp = C->i128 ;                      \
                                  MM4R_SINGLE_OP_INTERNAL_LT_64 ;              \
                                  C++->i128 = tmp ; ++A ; } ) ;                \
    }                                                                          \
  }

  #define MM4R_SINGLE_OP_INTERNAL_EQ_64                                        \
  if ( (tt = A->i32.n2) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }      \
  if ( (tt = A->i32.n3) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }

  #define MM4R_SINGLE_OP_INTERNAL_GT_64                                        \
  t64 = A->i64.hi>>rs ;                                                        \
  if ( (tt = (t64&0xFFFFFFFFU) ) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; } \
  if ( (tt = (t64>>32) ) )         { TT(4,0) ; TT(5,1) ; TT(6,2) ; TT(7,3) ; }

  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }    \
  if ( (tt = t128.i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }    \
  if ( (tt = t128.i32.n2) ) { TT(8,0)  ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }    \
  if ( (tt = t128.i32.n3) ) { TT(12,0) ; TT(13,1) ; TT(14,2) ; TT(15,3) ; }

  MULT_BLOCK(16)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }    \
  if ( (tt = t128.i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }    \
  if ( (tt = t128.i32.n2) ) { TT(8,0)  ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }    \
  if ( (tt = t128.i32.n3) ) { TT(12,0) ; TT(13,1) ; TT(14,2) ; }

  MULT_BLOCK(15)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }    \
  if ( (tt = t128.i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }    \
  if ( (tt = t128.i32.n2) ) { TT(8,0)  ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }    \
  tt = t128.i32.n3 ; TT(12,0) ; TT(13,1)

  MULT_BLOCK(14)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }    \
  if ( (tt = t128.i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }    \
  if ( (tt = t128.i32.n2) ) { TT(8,0)  ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }    \
  tt = t128.i32.n3 ; TT(12,0)

  MULT_BLOCK(13)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }    \
  if ( (tt = t128.i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }    \
  if ( (tt = t128.i32.n2) ) { TT(8,0)  ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }    \

  MULT_BLOCK(12)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }    \
  if ( (tt = t128.i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }    \
  if ( (tt = t128.i32.n2) ) { TT(8,0)  ; TT(9,1)  ; TT(10,2) ; }

  MULT_BLOCK(11)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }    \
  if ( (tt = t128.i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }    \
  tt = t128.i32.n2 ; TT(8,0) ; TT(9,1)

  MULT_BLOCK(10)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }    \
  if ( (tt = t128.i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }    \
  tt = t128.i32.n2 ; TT(8,0)

  MULT_BLOCK(9)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }    \
  if ( (tt = t128.i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }

  MULT_BLOCK(8)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }        \
  if ( (tt = t128.i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; }
  
  #undef  MM4R_SINGLE_OP_INTERNAL_EQ_64
  #define MM4R_SINGLE_OP_INTERNAL_EQ_64                                        \
  if ( (tt = A->i32.n2) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n3) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; }

  #undef  MM4R_SINGLE_OP_INTERNAL_GT_64
  #define MM4R_SINGLE_OP_INTERNAL_GT_64                                        \
  t64 = A->i64.hi>>rs ;                                                        \
  if ( (tt = (t64&0xFFFFFFFFU) ) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; } \
  if ( (tt = (t64>>32) ) )         { TT(4,0) ; TT(5,1) ; TT(6,2) ; }

  MULT_BLOCK(7)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }        \
  tt = t128.i32.n1 ; TT(4,0) ; TT(5,1)
  
  #undef  MM4R_SINGLE_OP_INTERNAL_EQ_64
  #define MM4R_SINGLE_OP_INTERNAL_EQ_64                                        \
  if ( (tt = A->i32.n2) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  tt = A->i32.n3 ; TT(4,0) ; TT(5,1)

  #undef  MM4R_SINGLE_OP_INTERNAL_GT_64
  #define MM4R_SINGLE_OP_INTERNAL_GT_64                                        \
  t64 = A->i64.hi>>rs ;                                                        \
  if ( (tt = (t64&0xFFFFFFFFU) ) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; } \
  tt = (t64>>32) ; TT(4,0) ; TT(5,1)

  MULT_BLOCK(6)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }        \
  tt = t128.i32.n1 ; TT(4,0)
  
  #undef  MM4R_SINGLE_OP_INTERNAL_EQ_64
  #define MM4R_SINGLE_OP_INTERNAL_EQ_64                                        \
  if ( (tt = A->i32.n2) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  tt = A->i32.n3 ; TT(4,0)

  #undef  MM4R_SINGLE_OP_INTERNAL_GT_64
  #define MM4R_SINGLE_OP_INTERNAL_GT_64                                        \
  t64 = A->i64.hi>>rs ;                                                        \
  if ( (tt = (t64&0xFFFFFFFFU) ) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; } \
  tt = (t64>>32) ; TT(4,0)

  MULT_BLOCK(5)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }
  
  #undef  MM4R_SINGLE_OP_INTERNAL_EQ_64
  #define MM4R_SINGLE_OP_INTERNAL_EQ_64                                        \
  if ( (tt = A->i32.n2) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }

  #undef  MM4R_SINGLE_OP_INTERNAL_GT_64
  #define MM4R_SINGLE_OP_INTERNAL_GT_64                                        \
  t64 = A->i64.hi>>rs ;                                                        \
  if ( (tt = (t64&0xFFFFFFFFU) ) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }

  MULT_BLOCK(4)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  if ( (tt = t128.i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; }
  
  #undef  MM4R_SINGLE_OP_INTERNAL_EQ_64
  #define MM4R_SINGLE_OP_INTERNAL_EQ_64                                        \
  if ( (tt = A->i32.n2) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; }

  #undef  MM4R_SINGLE_OP_INTERNAL_GT_64
  #define MM4R_SINGLE_OP_INTERNAL_GT_64                                        \
  t64 = A->i64.hi>>rs ;                                                        \
  if ( (tt = (t64&0xFFFFFFFFU) ) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; }

  MULT_BLOCK(3)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  tt = t128.i32.n0 ; TT(0,0) ; TT(1,1)
  
  #undef  MM4R_SINGLE_OP_INTERNAL_EQ_64
  #define MM4R_SINGLE_OP_INTERNAL_EQ_64                                        \
  tt = A->i32.n2 ; TT(0,0) ; TT(1,1)

  #undef  MM4R_SINGLE_OP_INTERNAL_GT_64
  #define MM4R_SINGLE_OP_INTERNAL_GT_64                                        \
  t64 = A->i64.hi>>rs ;                                                        \
  tt = (t64&0xFFFFFFFFU) ; TT(0,0) ; TT(1,1)

  MULT_BLOCK(2)

  #undef  MM4R_SINGLE_OP_INTERNAL_LT_64
  #define MM4R_SINGLE_OP_INTERNAL_LT_64                                        \
  t128.i64.lo = (A->i64.lo>>rshift) | (A->i64.hi<<lshift) ;                    \
  t128.i64.hi = A->i64.hi>>rshift ;                                            \
  tt = t128.i32.n0 ; TT(0,0)
  
  #undef  MM4R_SINGLE_OP_INTERNAL_EQ_64
  #define MM4R_SINGLE_OP_INTERNAL_EQ_64                                        \
  tt = A->i32.n2 ; TT(0,0)

  #undef  MM4R_SINGLE_OP_INTERNAL_GT_64
  #define MM4R_SINGLE_OP_INTERNAL_GT_64                                        \
  t64 = A->i64.hi>>rs ;                                                        \
  tt = (t64&0xFFFFFFFFU) ; TT(0,0)

  MULT_BLOCK(1)

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4Rreduced<uint128_t>::multRightShift( uint128_t * A, unsigned nrows, unsigned rshift ) const {
    switch ( nblock ) {
      case  1: MM4RmultRightBlock1 (A,nrows,Tables,rshift) ; break ;
      case  2: MM4RmultRightBlock2 (A,nrows,Tables,rshift) ; break ;
      case  3: MM4RmultRightBlock3 (A,nrows,Tables,rshift) ; break ;
      case  4: MM4RmultRightBlock4 (A,nrows,Tables,rshift) ; break ;
      case  5: MM4RmultRightBlock5 (A,nrows,Tables,rshift) ; break ;
      case  6: MM4RmultRightBlock6 (A,nrows,Tables,rshift) ; break ;
      case  7: MM4RmultRightBlock7 (A,nrows,Tables,rshift) ; break ;
      case  8: MM4RmultRightBlock8 (A,nrows,Tables,rshift) ; break ;
      case  9: MM4RmultRightBlock9 (A,nrows,Tables,rshift) ; break ;
      case 10: MM4RmultRightBlock10(A,nrows,Tables,rshift) ; break ;
      case 11: MM4RmultRightBlock11(A,nrows,Tables,rshift) ; break ;
      case 12: MM4RmultRightBlock12(A,nrows,Tables,rshift) ; break ;
      case 13: MM4RmultRightBlock13(A,nrows,Tables,rshift) ; break ;
      case 14: MM4RmultRightBlock14(A,nrows,Tables,rshift) ; break ;
      case 15: MM4RmultRightBlock15(A,nrows,Tables,rshift) ; break ;
      case 16: MM4RmultRightBlock16(A,nrows,Tables,rshift) ; break ;
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4Rreduced<uint128_t>::multRightAssShift( uint128_t const * A, unsigned nrows, 
                                             uint128_t       * C,
                                             unsigned          rshift ) const {
    switch ( nblock ) {
      case  1: MM4RmultAssBlock1 (A,nrows,Tables,C,rshift) ; break ;
      case  2: MM4RmultAssBlock2 (A,nrows,Tables,C,rshift) ; break ;
      case  3: MM4RmultAssBlock3 (A,nrows,Tables,C,rshift) ; break ;
      case  4: MM4RmultAssBlock4 (A,nrows,Tables,C,rshift) ; break ;
      case  5: MM4RmultAssBlock5 (A,nrows,Tables,C,rshift) ; break ;
      case  6: MM4RmultAssBlock6 (A,nrows,Tables,C,rshift) ; break ;
      case  7: MM4RmultAssBlock7 (A,nrows,Tables,C,rshift) ; break ;
      case  8: MM4RmultAssBlock8 (A,nrows,Tables,C,rshift) ; break ;
      case  9: MM4RmultAssBlock9 (A,nrows,Tables,C,rshift) ; break ;
      case 10: MM4RmultAssBlock10(A,nrows,Tables,C,rshift) ; break ;
      case 11: MM4RmultAssBlock11(A,nrows,Tables,C,rshift) ; break ;
      case 12: MM4RmultAssBlock12(A,nrows,Tables,C,rshift) ; break ;
      case 13: MM4RmultAssBlock13(A,nrows,Tables,C,rshift) ; break ;
      case 14: MM4RmultAssBlock14(A,nrows,Tables,C,rshift) ; break ;
      case 15: MM4RmultAssBlock15(A,nrows,Tables,C,rshift) ; break ;
      case 16: MM4RmultAssBlock16(A,nrows,Tables,C,rshift) ; break ;
    }  
  }
  
  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4Rreduced<uint128_t>::multRightAddShift( uint128_t const * A, unsigned nrows, 
                                             uint128_t       * C,
                                             unsigned          rshift ) const {
    switch ( nblock ) {
      case  1: MM4RmultAddBlock1 (A,nrows,Tables,C,rshift) ; break ;
      case  2: MM4RmultAddBlock2 (A,nrows,Tables,C,rshift) ; break ;
      case  3: MM4RmultAddBlock3 (A,nrows,Tables,C,rshift) ; break ;
      case  4: MM4RmultAddBlock4 (A,nrows,Tables,C,rshift) ; break ;
      case  5: MM4RmultAddBlock5 (A,nrows,Tables,C,rshift) ; break ;
      case  6: MM4RmultAddBlock6 (A,nrows,Tables,C,rshift) ; break ;
      case  7: MM4RmultAddBlock7 (A,nrows,Tables,C,rshift) ; break ;
      case  8: MM4RmultAddBlock8 (A,nrows,Tables,C,rshift) ; break ;
      case  9: MM4RmultAddBlock9 (A,nrows,Tables,C,rshift) ; break ;
      case 10: MM4RmultAddBlock10(A,nrows,Tables,C,rshift) ; break ;
      case 11: MM4RmultAddBlock11(A,nrows,Tables,C,rshift) ; break ;
      case 12: MM4RmultAddBlock12(A,nrows,Tables,C,rshift) ; break ;
      case 13: MM4RmultAddBlock13(A,nrows,Tables,C,rshift) ; break ;
      case 14: MM4RmultAddBlock14(A,nrows,Tables,C,rshift) ; break ;
      case 15: MM4RmultAddBlock15(A,nrows,Tables,C,rshift) ; break ;
      case 16: MM4RmultAddBlock16(A,nrows,Tables,C,rshift) ; break ;
    }  
  }

}

// EOF GF2toolkit_4RussianMultReducedShift.cc
