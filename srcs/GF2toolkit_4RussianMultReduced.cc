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

  template <>
  void
  MM4Rreduced<uint32_t>::makeTable( uint32_t const * M, unsigned row_used ) {
    nblock = (row_used+7)>>3 ;
    uint32_t * Tb = Tables ;

    while ( row_used > 7 )
      { MM4RmakeTable<8>( M, Tb ) ; M += 8 ; Tb += 1<<8 ; row_used -= 8 ; }
    
    if ( row_used > 0 ) {
      uint32_t MM[8] ;
      std::copy( M, M+row_used, MM ) ;
      while ( row_used < 8 ) MM[row_used++] = 0 ;
      MM4RmakeTable<8>( MM, Tb ) ;
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  #define TT(N) B_Tables[N*(1<<8)+BITS8(tmp,N)]

  static
  void
  MM4RmultRightBlock1( uint32_t * A, unsigned nRowsA, uint32_t const * B_Tables ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = *A ; *A++ = TT(0) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }
  
  static
  void
  MM4RmultRightBlock2( uint32_t * A, unsigned nRowsA, uint32_t const * B_Tables ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = *A ; *A++ = TT(0)^TT(1) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock3( uint32_t * A, unsigned nRowsA, uint32_t const * B_Tables ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = *A ; *A++ = TT(0)^TT(1)^TT(2) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock4( uint32_t * A, unsigned nRowsA, uint32_t const * B_Tables ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = *A ; *A++ = TT(0)^TT(1)^TT(2)^TT(3) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  template <>
  void
  MM4Rreduced<uint32_t>::multRight( uint32_t * A, unsigned nrows ) const {
    switch ( nblock ) {
      case 1: MM4RmultRightBlock1(A,nrows,Tables) ; break ;
      case 2: MM4RmultRightBlock2(A,nrows,Tables) ; break ;
      case 3: MM4RmultRightBlock3(A,nrows,Tables) ; break ;
      case 4: MM4RmultRightBlock4(A,nrows,Tables) ; break ;
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////

  static
  void
  MM4RmultAssBlock1( uint32_t const * A, unsigned nRowsA, uint32_t const * B_Tables, uint32_t * C ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; *C++ = TT(0) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }
  
  static
  void
  MM4RmultAssBlock2( uint32_t const * A, unsigned nRowsA, uint32_t const * B_Tables, uint32_t * C ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; *C++ = TT(0)^TT(1) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock3( uint32_t const * A, unsigned nRowsA, uint32_t const * B_Tables, uint32_t * C ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; *C++ = TT(0)^TT(1)^TT(2) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock4( uint32_t const * A, unsigned nRowsA, uint32_t const * B_Tables, uint32_t * C ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; *C++ = TT(0)^TT(1)^TT(2)^TT(3) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  template <>
  void
  MM4Rreduced<uint32_t>::multRightAss( uint32_t const * A, unsigned nrows, uint32_t * C ) const {
    switch ( nblock ) {
      case 1: MM4RmultAssBlock1(A,nrows,Tables,C) ; break ;
      case 2: MM4RmultAssBlock2(A,nrows,Tables,C) ; break ;
      case 3: MM4RmultAssBlock3(A,nrows,Tables,C) ; break ;
      case 4: MM4RmultAssBlock4(A,nrows,Tables,C) ; break ;
    }  
  }
  
  //////////////////////////////////////////////////////////////////////////////

  static
  void
  MM4RmultAddBlock1( uint32_t const * A, unsigned nRowsA, uint32_t const * B_Tables, uint32_t * C ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; *C++ ^= TT(0) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }
  
  static
  void
  MM4RmultAddBlock2( uint32_t const * A, unsigned nRowsA, uint32_t const * B_Tables, uint32_t * C ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; *C++ ^= TT(0)^TT(1) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock3( uint32_t const * A, unsigned nRowsA, uint32_t const * B_Tables, uint32_t * C ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; *C++ ^= TT(0)^TT(1)^TT(2) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock4( uint32_t const * A, unsigned nRowsA, uint32_t const * B_Tables, uint32_t * C ) {
    #define MM4R_SINGLE_OP { uint32_t tmp = *A++ ; *C++ ^= TT(0)^TT(1)^TT(2)^TT(3) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  template <>
  void
  MM4Rreduced<uint32_t>::multRightAdd( uint32_t const * A, unsigned nrows, uint32_t * C ) const {
    switch ( nblock ) {
      case 1: MM4RmultAddBlock1(A,nrows,Tables,C) ; break ;
      case 2: MM4RmultAddBlock2(A,nrows,Tables,C) ; break ;
      case 3: MM4RmultAddBlock3(A,nrows,Tables,C) ; break ;
      case 4: MM4RmultAddBlock4(A,nrows,Tables,C) ; break ;
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

  template <>
  void
  MM4Rreduced<uint64_t>::makeTable( uint64_t const * M, unsigned row_used ) {
    nblock = (row_used+7)>>3 ;
    uint64_t * Tb = Tables ;

    while ( row_used > 7 )
      { MM4RmakeTable<8>( M, Tb ) ; M += 8 ; Tb += 1<<8 ; row_used -= 8 ; }

    if ( row_used > 0 ) {
      uint64_t MM[8] ;
      std::copy( M, M+row_used, MM ) ;
      while ( row_used < 8 ) MM[row_used++] = 0 ;
      MM4RmakeTable<8>( MM, Tb ) ;
    }

  }
  
  //////////////////////////////////////////////////////////////////////////////

  static
  void
  MM4RmultRightBlock1( uint64_t * A, unsigned nRowsA, uint64_t const * B_Tables ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A ; *A++ = TT(0) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock2( uint64_t * A, unsigned nRowsA, uint64_t const * B_Tables ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A ; *A++ = TT(0)^TT(1) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock3( uint64_t * A, unsigned nRowsA, uint64_t const * B_Tables ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A ; *A++ = TT(0)^TT(1)^TT(2) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock4( uint64_t * A, unsigned nRowsA, uint64_t const * B_Tables ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A ; *A++ = TT(0)^TT(1)^TT(2)^TT(3) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock5( uint64_t * A, unsigned nRowsA, uint64_t const * B_Tables ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A ; *A++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock6( uint64_t * A, unsigned nRowsA, uint64_t const * B_Tables ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A ; *A++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock7( uint64_t * A, unsigned nRowsA, uint64_t const * B_Tables ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A ; *A++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultRightBlock8( uint64_t * A, unsigned nRowsA, uint64_t const * B_Tables ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A ; *A++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6)^TT(7) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  template <>
  void
  MM4Rreduced<uint64_t>::multRight( uint64_t * A, unsigned nrows ) const {
    switch ( nblock ) {
      case 1: MM4RmultRightBlock1(A,nrows,Tables) ; break ;
      case 2: MM4RmultRightBlock2(A,nrows,Tables) ; break ;
      case 3: MM4RmultRightBlock3(A,nrows,Tables) ; break ;
      case 4: MM4RmultRightBlock4(A,nrows,Tables) ; break ;
      case 5: MM4RmultRightBlock5(A,nrows,Tables) ; break ;
      case 6: MM4RmultRightBlock6(A,nrows,Tables) ; break ;
      case 7: MM4RmultRightBlock7(A,nrows,Tables) ; break ;
      case 8: MM4RmultRightBlock8(A,nrows,Tables) ; break ;
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////

  static
  void
  MM4RmultAssBlock1( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ = TT(0) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }
  
  static
  void
  MM4RmultAssBlock2( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ = TT(0)^TT(1) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock3( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ = TT(0)^TT(1)^TT(2) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock4( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ = TT(0)^TT(1)^TT(2)^TT(3) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock5( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock6( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock7( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAssBlock8( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ = TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6)^TT(7) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  template <>
  void
  MM4Rreduced<uint64_t>::multRightAss( uint64_t const * A, unsigned nrows, uint64_t * C ) const {
    switch ( nblock ) {
      case 1: MM4RmultAssBlock1(A,nrows,Tables,C) ; break ;
      case 2: MM4RmultAssBlock2(A,nrows,Tables,C) ; break ;
      case 3: MM4RmultAssBlock3(A,nrows,Tables,C) ; break ;
      case 4: MM4RmultAssBlock4(A,nrows,Tables,C) ; break ;
      case 5: MM4RmultAssBlock5(A,nrows,Tables,C) ; break ;
      case 6: MM4RmultAssBlock6(A,nrows,Tables,C) ; break ;
      case 7: MM4RmultAssBlock7(A,nrows,Tables,C) ; break ;
      case 8: MM4RmultAssBlock8(A,nrows,Tables,C) ; break ;
    }  
  }
  
  //////////////////////////////////////////////////////////////////////////////

  static
  void
  MM4RmultAddBlock1( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ ^= TT(0) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }
  
  static
  void
  MM4RmultAddBlock2( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ ^= TT(0)^TT(1) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock3( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ ^= TT(0)^TT(1)^TT(2) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock4( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ ^= TT(0)^TT(1)^TT(2)^TT(3) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock5( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ ^= TT(0)^TT(1)^TT(2)^TT(3)^TT(4) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock6( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ ^= TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock7( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ ^= TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  static
  void
  MM4RmultAddBlock8( uint64_t const * A, unsigned nRowsA, uint64_t const * B_Tables, uint64_t * C ) {
    #define MM4R_SINGLE_OP { uint64_t tmp = *A++ ; *C++ ^= TT(0)^TT(1)^TT(2)^TT(3)^TT(4)^TT(5)^TT(6)^TT(7) ; }
    MM4R_LOOP_UNROLLED_BY_16 ;
    #undef MM4R_SINGLE_OP
  }

  template <>
  void
  MM4Rreduced<uint64_t>::multRightAdd( uint64_t const * A, unsigned nrows, uint64_t * C ) const {
    switch ( nblock ) {
      case 1: MM4RmultAddBlock1(A,nrows,Tables,C) ; break ;
      case 2: MM4RmultAddBlock2(A,nrows,Tables,C) ; break ;
      case 3: MM4RmultAddBlock3(A,nrows,Tables,C) ; break ;
      case 4: MM4RmultAddBlock4(A,nrows,Tables,C) ; break ;
      case 5: MM4RmultAddBlock5(A,nrows,Tables,C) ; break ;
      case 6: MM4RmultAddBlock6(A,nrows,Tables,C) ; break ;
      case 7: MM4RmultAddBlock7(A,nrows,Tables,C) ; break ;
      case 8: MM4RmultAddBlock8(A,nrows,Tables,C) ; break ;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  MM4Rreduced<uint128_t>::makeTable( uint128_t const * M, unsigned row_used ) {
    nblock = (row_used+7)>>3 ;
    uint128_t * Tb = Tables ;

    while ( row_used > 7 )
      { MM4RmakeTable<8>( M, Tb ) ; M += 8 ; Tb += 1<<8 ; row_used -= 8 ; }

    if ( row_used > 0 ) {
      uint128_t MM[8] ;
      std::copy( M, M+row_used, MM ) ;
      while ( row_used < 8 ) MM[row_used++].i128 = _mm_setzero_si128() ;
      MM4RmakeTable<8>( MM, Tb ) ;
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////

  #undef  TT
  #define TT(IDX1,IDX2) tmp = _mm_xor_si128( tmp, B_Tables[IDX1*(1<<8)+BITS8(tt,IDX2)].i128 )

  #define MM4R_SINGLE_OP { __m128i tmp = _mm_setzero_si128() ; MM4R_SINGLE_OP_INTERNAL ; A++->i128 = tmp ; }

  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }      \
  if ( (tt = A->i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }      \
  if ( (tt = A->i32.n2) ) { TT(8,0)  ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }      \
  if ( (tt = A->i32.n3) ) { TT(12,0) ; TT(13,1) ; TT(14,2) ; TT(15,3) ; }

  static
  void
  MM4RmultRightBlock16( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }      \
  if ( (tt = A->i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }      \
  if ( (tt = A->i32.n2) ) { TT(8,0)  ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }      \
  if ( (tt = A->i32.n3) ) { TT(12,0) ; TT(13,1) ; TT(14,2) ; }

  static
  void
  MM4RmultRightBlock15( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }       \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }       \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }       \
  if ( (tt = A->i32.n3) ) { TT(12,0) ; TT(13,1) ; }

  static
  void
  MM4RmultRightBlock14( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2)  ; TT(3,3)  ; }        \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2)  ; TT(7,3)  ; }        \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1) ; TT(10,2) ; TT(11,3) ; }        \
  if ( (tt = A->i32.n3) ) { TT(12,0) ; }

  static
  void
  MM4RmultRightBlock13( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2)  ; TT(3,3)  ; }        \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2)  ; TT(7,3)  ; }        \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1) ; TT(10,2) ; TT(11,3) ; }

  static
  void
  MM4RmultRightBlock12( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2)  ; TT(3,3) ; }         \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2)  ; TT(7,3) ; }         \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1) ; TT(10,2) ; }

  static
  void
  MM4RmultRightBlock11( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; TT(7,3) ; }          \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1) ; }

  static
  void
  MM4RmultRightBlock10( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; TT(7,3) ; }          \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; }

  static
  void
  MM4RmultRightBlock9( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; TT(7,3) ; }

  static
  void
  MM4RmultRightBlock8( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; }

  static
  void
  MM4RmultRightBlock7( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; }

  static
  void
  MM4RmultRightBlock6( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) TT(4,0) ;

  static
  void
  MM4RmultRightBlock5( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }

  static
  void
  MM4RmultRightBlock4( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; }

  static
  void
  MM4RmultRightBlock3( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; }

  static
  void
  MM4RmultRightBlock2( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL if ( (tt = A->i32.n0) ) { TT(0,0) ; }

  static
  void
  MM4RmultRightBlock1( uint128_t * A, unsigned nRowsA, uint128_t const * B_Tables ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  template <>
  void
  MM4Rreduced<uint128_t>::multRight( uint128_t * A, unsigned nrows ) const {
    switch ( nblock ) {
      case  1: MM4RmultRightBlock1 (A,nrows,Tables) ; break ;
      case  2: MM4RmultRightBlock2 (A,nrows,Tables) ; break ;
      case  3: MM4RmultRightBlock3 (A,nrows,Tables) ; break ;
      case  4: MM4RmultRightBlock4 (A,nrows,Tables) ; break ;
      case  5: MM4RmultRightBlock5 (A,nrows,Tables) ; break ;
      case  6: MM4RmultRightBlock6 (A,nrows,Tables) ; break ;
      case  7: MM4RmultRightBlock7 (A,nrows,Tables) ; break ;
      case  8: MM4RmultRightBlock8 (A,nrows,Tables) ; break ;
      case  9: MM4RmultRightBlock9 (A,nrows,Tables) ; break ;
      case 10: MM4RmultRightBlock10(A,nrows,Tables) ; break ;
      case 11: MM4RmultRightBlock11(A,nrows,Tables) ; break ;
      case 12: MM4RmultRightBlock12(A,nrows,Tables) ; break ;
      case 13: MM4RmultRightBlock13(A,nrows,Tables) ; break ;
      case 14: MM4RmultRightBlock14(A,nrows,Tables) ; break ;
      case 15: MM4RmultRightBlock15(A,nrows,Tables) ; break ;
      case 16: MM4RmultRightBlock16(A,nrows,Tables) ; break ;
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////

  #undef MM4R_SINGLE_OP
  #undef MM4R_SINGLE_OP_INTERNAL

  #define MM4R_SINGLE_OP { __m128i tmp = _mm_setzero_si128() ; MM4R_SINGLE_OP_INTERNAL ; C++->i128 = tmp ; ++A ; }

  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }      \
  if ( (tt = A->i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }      \
  if ( (tt = A->i32.n2) ) { TT(8,0)  ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }      \
  if ( (tt = A->i32.n3) ) { TT(12,0) ; TT(13,1) ; TT(14,2) ; TT(15,3) ; }

  static
  void
  MM4RmultAssBlock16( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }      \
  if ( (tt = A->i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }      \
  if ( (tt = A->i32.n2) ) { TT(8,0)  ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }      \
  if ( (tt = A->i32.n3) ) { TT(12,0) ; TT(13,1) ; TT(14,2) ; }

  static
  void
  MM4RmultAssBlock15( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }       \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }       \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }       \
  tt = A->i32.n3 ; TT(12,0) ; TT(13,1)

  static
  void
  MM4RmultAssBlock14( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2)  ; TT(3,3)  ; }        \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2)  ; TT(7,3)  ; }        \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1) ; TT(10,2) ; TT(11,3) ; }        \
  tt = A->i32.n3 ; TT(12,0)

  static
  void
  MM4RmultAssBlock13( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2)  ; TT(3,3)  ; }        \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2)  ; TT(7,3)  ; }        \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1) ; TT(10,2) ; TT(11,3) ; }

  static
  void
  MM4RmultAssBlock12( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2)  ; TT(3,3) ; }         \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2)  ; TT(7,3) ; }         \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1) ; TT(10,2) ; }

  static
  void
  MM4RmultAssBlock11( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; TT(7,3) ; }          \
  tt = A->i32.n2 ; TT(8,0) ; TT(9,1)

  static
  void
  MM4RmultAssBlock10( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; TT(7,3) ; }          \
  tt = A->i32.n2 ; TT(8,0)

  static
  void
  MM4RmultAssBlock9( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; TT(7,3) ; }

  static
  void
  MM4RmultAssBlock8( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; }

  static
  void
  MM4RmultAssBlock7( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  tt = A->i32.n1 ; TT(4,0) ; TT(5,1)

  static
  void
  MM4RmultAssBlock6( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  tt = A->i32.n1 ; TT(4,0)

  static
  void
  MM4RmultAssBlock5( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }

  static
  void
  MM4RmultAssBlock4( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; }

  static
  void
  MM4RmultAssBlock3( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL tt = A->i32.n0 ; TT(0,0) ; TT(1,1)

  static
  void
  MM4RmultAssBlock2( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL tt = A->i32.n0 ; TT(0,0)

  static
  void
  MM4RmultAssBlock1( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  template <>
  void
  MM4Rreduced<uint128_t>::multRightAss( uint128_t const * A, unsigned nrows, uint128_t * C ) const {
    switch ( nblock ) {
      case  1: MM4RmultAssBlock1 (A,nrows,Tables,C) ; break ;
      case  2: MM4RmultAssBlock2 (A,nrows,Tables,C) ; break ;
      case  3: MM4RmultAssBlock3 (A,nrows,Tables,C) ; break ;
      case  4: MM4RmultAssBlock4 (A,nrows,Tables,C) ; break ;
      case  5: MM4RmultAssBlock5 (A,nrows,Tables,C) ; break ;
      case  6: MM4RmultAssBlock6 (A,nrows,Tables,C) ; break ;
      case  7: MM4RmultAssBlock7 (A,nrows,Tables,C) ; break ;
      case  8: MM4RmultAssBlock8 (A,nrows,Tables,C) ; break ;
      case  9: MM4RmultAssBlock9 (A,nrows,Tables,C) ; break ;
      case 10: MM4RmultAssBlock10(A,nrows,Tables,C) ; break ;
      case 11: MM4RmultAssBlock11(A,nrows,Tables,C) ; break ;
      case 12: MM4RmultAssBlock12(A,nrows,Tables,C) ; break ;
      case 13: MM4RmultAssBlock13(A,nrows,Tables,C) ; break ;
      case 14: MM4RmultAssBlock14(A,nrows,Tables,C) ; break ;
      case 15: MM4RmultAssBlock15(A,nrows,Tables,C) ; break ;
      case 16: MM4RmultAssBlock16(A,nrows,Tables,C) ; break ;
    }  
  }
  
  //////////////////////////////////////////////////////////////////////////////

  #undef MM4R_SINGLE_OP
  #undef MM4R_SINGLE_OP_INTERNAL

  #define MM4R_SINGLE_OP { __m128i tmp = C->i128 ; MM4R_SINGLE_OP_INTERNAL ; C++->i128 = tmp ; ++A ; }

  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }      \
  if ( (tt = A->i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }      \
  if ( (tt = A->i32.n2) ) { TT(8,0)  ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }      \
  if ( (tt = A->i32.n3) ) { TT(12,0) ; TT(13,1) ; TT(14,2) ; TT(15,3) ; }

  static
  void
  MM4RmultAddBlock16( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0)  ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }      \
  if ( (tt = A->i32.n1) ) { TT(4,0)  ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }      \
  if ( (tt = A->i32.n2) ) { TT(8,0)  ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }      \
  if ( (tt = A->i32.n3) ) { TT(12,0) ; TT(13,1) ; TT(14,2) ; }

  static
  void
  MM4RmultAddBlock15( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1)  ; TT(2,2)  ; TT(3,3)  ; }       \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1)  ; TT(6,2)  ; TT(7,3)  ; }       \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1)  ; TT(10,2) ; TT(11,3) ; }       \
  tt = A->i32.n3 ; TT(12,0) ; TT(13,1)

  static
  void
  MM4RmultAddBlock14( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2)  ; TT(3,3)  ; }        \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2)  ; TT(7,3)  ; }        \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1) ; TT(10,2) ; TT(11,3) ; }        \
  tt = A->i32.n3 ; TT(12,0)

  static
  void
  MM4RmultAddBlock13( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2)  ; TT(3,3)  ; }        \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2)  ; TT(7,3)  ; }        \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1) ; TT(10,2) ; TT(11,3) ; }

  static
  void
  MM4RmultAddBlock12( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2)  ; TT(3,3) ; }         \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2)  ; TT(7,3) ; }         \
  if ( (tt = A->i32.n2) ) { TT(8,0) ; TT(9,1) ; TT(10,2) ; }

  static
  void
  MM4RmultAddBlock11( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; TT(7,3) ; }          \
  tt = A->i32.n2 ; TT(8,0) ; TT(9,1)

  static
  void
  MM4RmultAddBlock10( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; TT(7,3) ; }          \
  tt = A->i32.n2 ; TT(8,0)

  static
  void
  MM4RmultAddBlock9( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; TT(7,3) ; }

  static
  void
  MM4RmultAddBlock8( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  if ( (tt = A->i32.n1) ) { TT(4,0) ; TT(5,1) ; TT(6,2) ; }

  static
  void
  MM4RmultAddBlock7( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  tt = A->i32.n1 ; TT(4,0) ; TT(5,1)

  static
  void
  MM4RmultAddBlock6( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }          \
  tt = A->i32.n1 ; TT(4,0)

  static
  void
  MM4RmultAddBlock5( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; TT(3,3) ; }

  static
  void
  MM4RmultAddBlock4( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL                                              \
  if ( (tt = A->i32.n0) ) { TT(0,0) ; TT(1,1) ; TT(2,2) ; }

  static
  void
  MM4RmultAddBlock3( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL tt = A->i32.n0 ; TT(0,0) ; TT(1,1)

  static
  void
  MM4RmultAddBlock2( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  #undef  MM4R_SINGLE_OP_INTERNAL
  #define MM4R_SINGLE_OP_INTERNAL tt = A->i32.n0 ; TT(0,0)

  static
  void
  MM4RmultAddBlock1( uint128_t const * A, unsigned nRowsA, uint128_t const * B_Tables, uint128_t * C ) {
    uint32_t tt ;
    MM4R_LOOP_UNROLLED_BY_16 ;
  }

  template <>
  void
  MM4Rreduced<uint128_t>::multRightAdd( uint128_t const * A, unsigned nrows, uint128_t * C ) const {
    switch ( nblock ) {
      case  1: MM4RmultAddBlock1 (A,nrows,Tables,C) ; break ;
      case  2: MM4RmultAddBlock2 (A,nrows,Tables,C) ; break ;
      case  3: MM4RmultAddBlock3 (A,nrows,Tables,C) ; break ;
      case  4: MM4RmultAddBlock4 (A,nrows,Tables,C) ; break ;
      case  5: MM4RmultAddBlock5 (A,nrows,Tables,C) ; break ;
      case  6: MM4RmultAddBlock6 (A,nrows,Tables,C) ; break ;
      case  7: MM4RmultAddBlock7 (A,nrows,Tables,C) ; break ;
      case  8: MM4RmultAddBlock8 (A,nrows,Tables,C) ; break ;
      case  9: MM4RmultAddBlock9 (A,nrows,Tables,C) ; break ;
      case 10: MM4RmultAddBlock10(A,nrows,Tables,C) ; break ;
      case 11: MM4RmultAddBlock11(A,nrows,Tables,C) ; break ;
      case 12: MM4RmultAddBlock12(A,nrows,Tables,C) ; break ;
      case 13: MM4RmultAddBlock13(A,nrows,Tables,C) ; break ;
      case 14: MM4RmultAddBlock14(A,nrows,Tables,C) ; break ;
      case 15: MM4RmultAddBlock15(A,nrows,Tables,C) ; break ;
      case 16: MM4RmultAddBlock16(A,nrows,Tables,C) ; break ;
    }  
  }

}

// EOF GF2toolkit_4RussianMultreduced.cc
