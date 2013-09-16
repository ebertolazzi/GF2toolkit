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

#include "GF2toolkit_Matrix.hh"
#include "GF2toolkit_MatrixMult.hh"
#include "GF2toolkit_4Russian.hh"
#include "GF2toolkit_Strassen.hh"

namespace GF2toolkit {

  //#define MM_ADD MMaddRecurr<UNSIGNED>
  #define MM_ADD MMadd<UNSIGNED>

  /*
  //   __  __ __  __     _
  //  |  \/  |  \/  |___| |_ _ __ __ _ ___ ___  ___ _ __
  //  | |\/| | |\/| / __| __| '__/ _` / __/ __|/ _ \ '_ \
  //  | |  | | |  | \__ \ |_| | | (_| \__ \__ \  __/ | | |
  //  |_|  |_|_|  |_|___/\__|_|  \__,_|___/___/\___|_| |_|
  */
  template <typename UNSIGNED>
  void
  MMaddStrassen( UNSIGNED * work, unsigned nwork,
                 unsigned nBlock,
                 UNSIGNED const * Amat, unsigned dimRowA,
                 UNSIGNED const * Bmat, unsigned dimRowB,
                 UNSIGNED       * Cmat, unsigned dimRowC ) {

    unsigned const numBits = sizeof(UNSIGNED) * CHAR_BIT ;
              
    unsigned rmem = (nBlock*nBlock*numBits*3)/4 ;

    if ( numBits*nBlock <= STRASSEN_MIN_SIZE || nwork < rmem ) {
      MM_ADD( Amat, dimRowA, nBlock * numBits, nBlock,
              Bmat, dimRowB, nBlock * numBits, nBlock,
              Cmat, dimRowC ) ;
      return ;
    }
    
    if ( nBlock & 0x01 ) {
      unsigned nBorder = nBlock & STRASSEN_MASK ;
      nBlock &= ~STRASSEN_MASK ;
      unsigned N = nBlock  * numBits ;

      UNSIGNED const * A11 = Amat ;
      UNSIGNED const * A21 = A11 + N ;
      UNSIGNED const * A12 = A11 + nBlock * dimRowA ;

      UNSIGNED const * B11 = Bmat ;
      UNSIGNED const * B21 = B11 + N ;
      UNSIGNED const * B12 = B11 + nBlock * dimRowB ;

      UNSIGNED       * C11 = Cmat ;
      UNSIGNED       * C21 = C11 + N ;
      UNSIGNED       * C12 = C11 + nBlock * dimRowC ;
    
      /*
      //  +---------+     +---------+   +---------+
      //  | A11 A12 |     | B11 B12 |   | C11 C12 |
      //  |         |  x  |         | = |         |
      //  | A21 A22 |     | B21 B22 |   | C21 C22 |
      //  +---------+     +---------+   +---------+
      //
      //                  +---------+
      //                  | B11 B12 |
      //  +---------+  x  |         | = +---------+
      //  | A21 A22 |     | B21 B22 |   | C21 C22 |
      //  +---------+     +---------+   +---------+
      //
      //  +=========================+
      //  | +-----+     +---------+ |   +-----+     +---------+   +-----+     +---------+   +---------+
      //  | | A11 |  x  | B11   0 | | + | A11 |  x  |  0  B12 | + | A12 |  x  | B21 B22 | = | C11 C12 |
      //  | +-----+     +---------+ |   +-----+     +---------+   +-----+     +---------+   +---------+
      //  +=========================+
      */
      MM_ADD( A21, dimRowA, nBorder*numBits,          nBlock+nBorder,
              B11, dimRowB, (nBlock+nBorder)*numBits, nBlock+nBorder,
              C21, dimRowC ) ;
    
      MM_ADD( A11, dimRowA, N, nBlock,
              B12, dimRowB, N, nBorder,
              C12, dimRowC ) ;

      MM_ADD( A12, dimRowA, N,               nBorder,
              B21, dimRowB, nBorder*numBits, nBlock+nBorder,
              C11, dimRowC ) ;
    }

    nBlock >>= 1 ;
    unsigned N = nBlock * numBits ;
    /*
    //  +         +     +         +   +         +
    //  | A11 A12 |     | B11 B12 |   | C11 C12 |
    //  |         |  x  |         | = |         |
    //  | A21 A22 |     | B21 B22 |   | C21 C22 |
    //  +         +     +         +   +         +
    */

    UNSIGNED const * A11 = Amat ;
    UNSIGNED const * A21 = A11 + N ;
    UNSIGNED const * A12 = A11 + nBlock * dimRowA ;
    UNSIGNED const * A22 = A12 + N ;

    UNSIGNED const * B11 = Bmat ;
    UNSIGNED const * B21 = B11 + N ;
    UNSIGNED const * B12 = B11 + nBlock * dimRowB ;
    UNSIGNED const * B22 = B12 + N ;

    UNSIGNED       * C11 = Cmat ;
    UNSIGNED       * C21 = C11 + N ;
    UNSIGNED       * C12 = C11 + nBlock * dimRowC ;
    UNSIGNED       * C22 = C12 + N ;
    
    UNSIGNED * BF1 = work ; work += N*nBlock ;
    UNSIGNED * BF2 = work ; work += N*nBlock ;
    UNSIGNED * PK  = work ; work += N*nBlock ;

    nwork -= rmem ;

    // P1 = (A11+A22)(B11+B22)
    add<UNSIGNED>( A11, dimRowA, A22, dimRowA, BF1, N, N, nBlock ) ; // BF1 = A11 + A22
    add<UNSIGNED>( B11, dimRowB, B22, dimRowB, BF2, N, N, nBlock ) ; // BF2 = B11 + B22
    zero<UNSIGNED>( PK, N*nBlock ) ;
    MMaddStrassen<UNSIGNED>( work, nwork, nBlock, BF1, N, BF2, N, PK, N ) ;
    addTo<UNSIGNED>( PK, N, C11, dimRowC, N, nBlock ) ;
    addTo<UNSIGNED>( PK, N, C22, dimRowC, N, nBlock ) ;

    // P2 = (A21+A22)B11
    add<UNSIGNED>( A21, dimRowA, A22, dimRowA, BF1, N, N, nBlock ) ; // BF1 = A21 + A22
    zero<UNSIGNED>( PK, N*nBlock ) ;
    MMaddStrassen<UNSIGNED>( work, nwork, nBlock, BF1, N, B11, dimRowB, PK, N ) ; 
    addTo<UNSIGNED>( PK, N, C21, dimRowC, N, nBlock ) ;
    addTo<UNSIGNED>( PK, N, C22, dimRowC, N, nBlock ) ;

    // P3 = A11(B12-B22)
    add<UNSIGNED>( B12, dimRowB, B22, dimRowB, BF1, N, N, nBlock ) ; // BF1 = B12 + B22
    zero<UNSIGNED>( PK, N*nBlock ) ;
    MMaddStrassen<UNSIGNED>( work, nwork, nBlock, A11, dimRowA, BF1, N, PK, N ) ; 
    addTo<UNSIGNED>( PK, N, C12, dimRowC, N, nBlock ) ;
    addTo<UNSIGNED>( PK, N, C22, dimRowC, N, nBlock ) ;

    // P4 = A22(B21-B11)
    add<UNSIGNED>( B21, dimRowB, B11, dimRowB, BF1, N, N, nBlock ) ; // BF1 = B21 + B11
    zero<UNSIGNED>( PK, N*nBlock ) ;
    MMaddStrassen<UNSIGNED>( work, nwork, nBlock, A22, dimRowA, BF1, N, PK, N ) ; 
    addTo<UNSIGNED>( PK, N, C11, dimRowC, N, nBlock ) ;
    addTo<UNSIGNED>( PK, N, C21, dimRowC, N, nBlock ) ;

    // P5 = (A11+A12)B22
    add<UNSIGNED>( A11, dimRowA, A12, dimRowA, BF1, N, N, nBlock ) ; // BF1 = A11 + A12
    zero<UNSIGNED>( PK, N*nBlock ) ;
    MMaddStrassen<UNSIGNED>( work, nwork, nBlock, BF1, N, B22, dimRowB, PK, N ) ; 
    addTo<UNSIGNED>( PK, N, C11, dimRowC, N, nBlock ) ;
    addTo<UNSIGNED>( PK, N, C12, dimRowC, N, nBlock ) ;

    // P6 = (A21-A11)(B11+B12)
    add<UNSIGNED>( A21, dimRowA, A11, dimRowA, BF1, N, N, nBlock ) ; // BF1 = A21 + A11
    add<UNSIGNED>( B11, dimRowB, B12, dimRowB, BF2, N, N, nBlock ) ; // BF2 = B11 + B12
    zero<UNSIGNED>( PK, N*nBlock ) ;
    MMaddStrassen<UNSIGNED>( work, nwork, nBlock, BF1, N, BF2, N, PK, N ) ; 
    addTo<UNSIGNED>( PK, N, C22, dimRowC, N, nBlock ) ;

    // P7 = (A12-A22)(B21+B22)
    add<UNSIGNED>( A12, dimRowA, A22, dimRowA, BF1, N, N, nBlock ) ; // BF1 = A12 + A22
    add<UNSIGNED>( B21, dimRowB, B22, dimRowB, BF2, N, N, nBlock ) ; // BF2 = B21 + B22
    zero<UNSIGNED>( PK, N*nBlock ) ;
    MMaddStrassen<UNSIGNED>( work, nwork, nBlock, BF1, N, BF2, N, PK, N ) ;
    addTo<UNSIGNED>( PK, N, C11, dimRowC, N, nBlock ) ;

    //////////////////
    // C11 += P1 +           P4 - P5 +     P7
    // C21 +=      P2 +      P4
    // C12 +=           P3 +      P5
    // C22 += P1 - P2 + P3 +           P6

  }

  //! C += A*B
  template <typename UNSIGNED>
  void
  MMaddStrassen( UNSIGNED * work, unsigned nwork,
                 UNSIGNED const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
                 UNSIGNED const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
                 UNSIGNED       * C, unsigned dimRowsC ) {

    unsigned const numBits = sizeof(UNSIGNED) * CHAR_BIT ;

    if ( nColsBlockA*numBits <= STRASSEN_MIN_SIZE ||
         nColsBlockB*numBits <= STRASSEN_MIN_SIZE ||
         nRowsA              <= STRASSEN_MIN_SIZE ||
         nRowsB              <= STRASSEN_MIN_SIZE ) {
      MM_ADD( A, dimRowsA, nRowsA, nColsBlockA,
              B, dimRowsB, nRowsB, nColsBlockB,
              C, dimRowsC ) ;
      return ;
    }

    unsigned nBlock = std::min( std::min( nColsBlockA, nColsBlockB),
                                std::min( nRowsA, nRowsB) /numBits ) & ~STRASSEN_MASK ;

    MMaddStrassen<UNSIGNED>( work, nwork, nBlock, A, dimRowsA, B, dimRowsB, C, dimRowsC ) ;

    unsigned N = nBlock * numBits ;

    UNSIGNED const * A11 = A ;
    UNSIGNED const * A21 = A11 + N ;
    UNSIGNED const * A12 = A11 + nBlock * dimRowsA ;

    UNSIGNED const * B11 = B ;
    UNSIGNED const * B21 = B11 + N ;
    UNSIGNED const * B12 = B11 + nBlock * dimRowsB ;

    UNSIGNED       * C11 = C ;
    UNSIGNED       * C21 = C11 + N ;
    UNSIGNED       * C12 = C11 + nBlock * dimRowsC ;

    /*
    //  +---------+     +---------+   +---------+
    //  | A11 A12 |     | B11 B12 |   | C11 C12 |
    //  |         |  x  |         | = |         |
    //  | A21 A22 |     | B21 B22 |   | C21 C22 |
    //  +---------+     +---------+   +---------+
    //
    //                  +---------+
    //                  | B11 B12 |
    //  +---------+  x  |         | = +---------+
    //  | A21 A22 |     | B21 B22 |   | C21 C22 |
    //  +---------+     +---------+   +---------+
    //
    //  +=========================+
    //  | +-----+     +---------+ |   +-----+     +---------+   +-----+     +---------+   +---------+
    //  | | A11 |  x  | B11   0 | | + | A11 |  x  |  0  B12 | + | A12 |  x  | B21 B22 | = | C11 C12 |
    //  | +-----+     +---------+ |   +-----+     +---------+   +-----+     +---------+   +---------+
    //  +=========================+
    */
    #if 1
    MM_ADD( A21, dimRowsA, nRowsA-N, nColsBlockA,
            B11, dimRowsB, nRowsB,   nColsBlockB,
            C21, dimRowsC ) ;
    
    MM_ADD( A11, dimRowsA, N, nBlock,
            B12, dimRowsB, N, nColsBlockB - nBlock,
            C12, dimRowsC ) ;

    MM_ADD( A12, dimRowsA, N,          nColsBlockA - nBlock,
            B21, dimRowsB, nRowsB - N, nColsBlockB,
            C11, dimRowsC ) ;
    #else
    MMaddStrassen<UNSIGNED>( work, nwork, 
                             A21, dimRowsA, nRowsA-N, nColsBlockA,
                             B11, dimRowsB, nRowsB,   nColsBlockB,
                             C21, dimRowsC ) ;
    
    MMaddStrassen<UNSIGNED>( work, nwork, 
                             A11, dimRowsA, N, nBlock,
                             B12, dimRowsB, N, nColsBlockB - nBlock,
                             C12, dimRowsC ) ;

    MMaddStrassen<UNSIGNED>( work, nwork,
                             A12, dimRowsA, N,          nColsBlockA - nBlock,
                             B21, dimRowsB, nRowsB - N, nColsBlockB,
                             C11, dimRowsC ) ;
    #endif
  }
  
  #define EXPLICT_INSTANTIATION(TYPE) \
    template void MMaddStrassen<TYPE>( TYPE * work, unsigned nwork, \
                                       unsigned nBlock, \
                                       TYPE const * Amat, unsigned dimRowA, \
                                       TYPE const * Bmat, unsigned dimRowB, \
                                       TYPE       * Cmat, unsigned dimRowC ) ; \
    \
    template void MMaddStrassen<TYPE>( TYPE * work, unsigned nwork, \
                                       TYPE const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA, \
                                       TYPE const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB, \
                                       TYPE       * C, unsigned dimRowsC )

  EXPLICT_INSTANTIATION(uint32_t) ;
  EXPLICT_INSTANTIATION(uint64_t) ;
  EXPLICT_INSTANTIATION(uint128_t) ;

}

// EOF GF2toolkit_Strassen.cc
