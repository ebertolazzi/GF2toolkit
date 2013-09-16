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

#include "GF2toolkit_MatrixMult.hh"
#include "GF2toolkit_4Russian.hh"

#define MM4R_32BIT_USE_10_BIT 6750u
#define MM4R_32BIT_USE_8_BIT  551u
#define MM4R_32BIT_USE_6_BIT  90u
#define MM4R_32BIT_USE_5_BIT  51u
#define MM4R_32BIT_USE_4_BIT  23u
#define MM4R_32BIT_USE_3_BIT  4u

#define MM4R_64BIT_USE_9_BIT  2012u
#define MM4R_64BIT_USE_8_BIT  942u
#define MM4R_64BIT_USE_7_BIT  212u
#define MM4R_64BIT_USE_6_BIT  83u
#define MM4R_64BIT_USE_5_BIT  24u
#define MM4R_64BIT_USE_3_BIT  4u

#define MM4R_128BIT_USE_8_BIT  1152u
#define MM4R_128BIT_USE_6_BIT  129u
#define MM4R_128BIT_USE_5_BIT  40u
#define MM4R_128BIT_USE_4_BIT  5u

namespace GF2toolkit {

  /*
  //  +---+ +---+    
  //  |   | | G |    
  //  |   | +---+
  //  |   |
  //  |   |
  //  |   |
  //  | A |
  //  |   |
  //  |   |
  //  |   |
  //  |   |
  //  |   |
  //  |   |
  //  +---+
  //
  //  moltiplico matrice colonna per matrice quadrata
  */
  
  template <>
  void
  MMassLeft( uint32_t * A, unsigned nRowsA, uint32_t const * G ) {
    if ( nRowsA > MM4R_32BIT_USE_10_BIT ) {
      MM4R<uint32_t,10> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_8_BIT ) {
      MM4R<uint32_t,8> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_6_BIT ) {
      MM4R<uint32_t,6> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_5_BIT ) {
      MM4R<uint32_t,5> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_4_BIT ) {
      MM4R<uint32_t,4> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_3_BIT ) {
      MM4R<uint32_t,3> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else {
      MM4R<uint32_t,2> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    }
  }

  template <>
  void
  MMassLeft( uint64_t * A, unsigned nRowsA, uint64_t const * G ) {
    if ( nRowsA > MM4R_64BIT_USE_9_BIT ) {
      MM4R<uint64_t,9> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_8_BIT ) {
      MM4R<uint64_t,8> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_7_BIT ) {
      MM4R<uint64_t,7> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_6_BIT ) {
      MM4R<uint64_t,6> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_5_BIT ) {
      MM4R<uint64_t,5> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_3_BIT ) {
      MM4R<uint64_t,3> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else {
      MM4R<uint64_t,2> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    }
  }

  template <>
  void
  MMassLeft( uint128_t * A, unsigned nRowsA, uint128_t const * G ) {
    if ( nRowsA > MM4R_128BIT_USE_8_BIT ) {
      MM4R<uint128_t,8> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_128BIT_USE_6_BIT ) {
      MM4R<uint128_t,6> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_128BIT_USE_5_BIT ) {
      MM4R<uint128_t,5> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else if ( nRowsA > MM4R_128BIT_USE_4_BIT ) {
      MM4R<uint128_t,4> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    } else {
      MM4R<uint128_t,2> mm4r( G ) ;
      mm4r.multRight( A, nRowsA ) ;
    }
  }

  // ---------------------------------------------------------------------------

  template <>
  void
  MMass( uint32_t const * A, unsigned nRowsA, uint32_t const * G, uint32_t * B ) {
    if ( nRowsA > MM4R_32BIT_USE_10_BIT ) {
      MM4R<uint32_t,10> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_8_BIT ) {
      MM4R<uint32_t,8> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_6_BIT ) {
      MM4R<uint32_t,6> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_5_BIT ) {
      MM4R<uint32_t,5> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_4_BIT ) {
      MM4R<uint32_t,4> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_3_BIT ) {
      MM4R<uint32_t,3> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else {
      MM4R<uint32_t,2> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    }
  }

  template <>
  void
  MMass( uint64_t const * A, unsigned nRowsA, uint64_t const * G, uint64_t * B ) {
    if ( nRowsA > MM4R_64BIT_USE_9_BIT ) {
      MM4R<uint64_t,9> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_8_BIT ) {
      MM4R<uint64_t,8> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_7_BIT ) {
      MM4R<uint64_t,7> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_6_BIT ) {
      MM4R<uint64_t,6> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_5_BIT ) {
      MM4R<uint64_t,5> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_3_BIT ) {
      MM4R<uint64_t,3> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else {
      MM4R<uint64_t,2> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    }
  }

  template <>
  void
  MMass( uint128_t const * A, unsigned nRowsA, uint128_t const * G, uint128_t * B ) {
    if ( nRowsA > MM4R_128BIT_USE_8_BIT ) {
      MM4R<uint128_t,8> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_128BIT_USE_6_BIT ) {
      MM4R<uint128_t,6> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_128BIT_USE_5_BIT ) {
      MM4R<uint128_t,5> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_128BIT_USE_4_BIT ) {
      MM4R<uint128_t,4> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    } else {
      MM4R<uint128_t,2> mm4r( G ) ;
      mm4r.multRightAss( A, nRowsA, B ) ;
    }
  }

  // ---------------------------------------------------------------------------

  template <>
  void
  MMadd( uint32_t const * A, unsigned nRowsA, uint32_t const * G, uint32_t * B ) {
    if ( nRowsA > MM4R_32BIT_USE_10_BIT ) {
      MM4R<uint32_t,10> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_8_BIT ) {
      MM4R<uint32_t,8> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_6_BIT ) {
      MM4R<uint32_t,6> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_5_BIT ) {
      MM4R<uint32_t,5> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_4_BIT ) {
      MM4R<uint32_t,4> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_32BIT_USE_3_BIT ) {
      MM4R<uint32_t,3> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else {
      MM4R<uint32_t,2> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    }
  }

  template <>
  void
  MMadd( uint64_t const * A, unsigned nRowsA, uint64_t const * G, uint64_t * B ) {
    if ( nRowsA > MM4R_64BIT_USE_9_BIT ) {
      MM4R<uint64_t,9> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_8_BIT ) {
      MM4R<uint64_t,8> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_7_BIT ) {
      MM4R<uint64_t,7> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_6_BIT ) {
      MM4R<uint64_t,6> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_5_BIT ) {
      MM4R<uint64_t,5> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_64BIT_USE_3_BIT ) {
      MM4R<uint64_t,3> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else {
      MM4R<uint64_t,2> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    }
  }

  template <>
  void
  MMadd( uint128_t const * A, unsigned nRowsA, uint128_t const * G, uint128_t * B ) {
    if ( nRowsA > MM4R_128BIT_USE_8_BIT ) {
      MM4R<uint128_t,8> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_128BIT_USE_6_BIT ) {
      MM4R<uint128_t,6> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_128BIT_USE_5_BIT ) {
      MM4R<uint128_t,5> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else if ( nRowsA > MM4R_128BIT_USE_4_BIT ) {
      MM4R<uint128_t,4> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    } else {
      MM4R<uint128_t,2> mm4r( G ) ;
      mm4r.multRightAdd( A, nRowsA, B ) ;
    }
  }
  
  /*
  //  +---+ +---+     +---+
  //  | G | | A | --> | A |
  //  +---+ +---+     +---+
  //
  //  moltiplico matrice quadrata per matrice riga
  */
  //! A = G*A
  template <>
  void
  MMassRight( uint32_t const * G, uint32_t * A ) {
    MM4R<uint32_t,4> mm4r ;
    mm4r.makeTable( A ) ;
    mm4r.multRightAss( G, 32, A ) ;      
  }

  template <>
  void
  MMassRight( uint64_t const * G, uint64_t * A ) {
    MM4R<uint64_t,5> mm4r ; // ottimo con 5
    //MM4R<uint64_t,4> mm4r ; // ottimo con 5
    mm4r.makeTable( A ) ;
    mm4r.multRightAss( G, 64, A ) ;      
  }

  template <>
  void
  MMassRight( uint128_t const * G, uint128_t * A ) {
    MM4R<uint128_t,6> mm4r ; // ottimo con 6
    mm4r.makeTable( A ) ;
    mm4r.multRightAss( G, 128, A ) ;      
  }

  /*
  //  +---+ +--------------------+     +--------------------+
  //  | G | |         A          | --> |         A          |
  //  +---+ +--------------------+     +--------------------+
  //
  //  moltiplico matrice quadrata per matrice riga
  */
  //! A = G*A
  template <>
  void
  MMassRight( uint32_t const * G, uint32_t * A, unsigned dimRowsA, unsigned nColsBlockA ) {
    MM4R<uint32_t,4> mm4r ;
    for ( unsigned colBlockA = 0 ; colBlockA < nColsBlockA ; ++colBlockA, A += dimRowsA ) {
      mm4r.makeTable( A ) ;
      mm4r.multRightAss( G, 32, A ) ;      
    }
  }

  template <>
  void
  MMassRight( uint64_t const * G, uint64_t * A, unsigned dimRowsA, unsigned nColsBlockA ) {
    MM4R<uint64_t,5> mm4r ; // ottimo con 5
    //MM4R<uint64_t,4> mm4r ; // ottimo con 5
    for ( unsigned colBlockA = 0 ; colBlockA < nColsBlockA ; ++colBlockA, A += dimRowsA ) {
      mm4r.makeTable( A ) ;
      mm4r.multRightAss( G, 64, A ) ;      
    }
  }

  template <>
  void
  MMassRight( uint128_t const * G, uint128_t * A, unsigned dimRowsA, unsigned nColsBlockA ) {
    //MM4R<uint128_t,5> mm4r ; // ottimo con 5
    MM4R<uint128_t,6> mm4r ; // ottimo con 6
    for ( unsigned colBlockA = 0 ; colBlockA < nColsBlockA ; ++colBlockA, A += dimRowsA ) {
      mm4r.makeTable( A ) ;
      mm4r.multRightAss( G, 128, A ) ;      
    }
  }

  /*
  //                   _        _                           _        _
  //   _ __ ___   __ _| |_ _ __(_)_  __     _ __ ___   __ _| |_ _ __(_)_  __
  //  | '_ ` _ \ / _` | __| '__| \ \/ /____| '_ ` _ \ / _` | __| '__| \ \/ /
  //  | | | | | | (_| | |_| |  | |>  <_____| | | | | | (_| | |_| |  | |>  < 
  //  |_| |_| |_|\__,_|\__|_|  |_/_/\_\    |_| |_| |_|\__,_|\__|_|  |_/_/\_\
  */

  //! C += A*B
  template <typename UNSIGNED>
  void
  MMass( UNSIGNED const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
         UNSIGNED const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
         UNSIGNED       * C, unsigned dimRowsC ) {

    unsigned const numBits = sizeof(UNSIGNED) * CHAR_BIT ;

    ASSERT( nColsBlockA * numBits == nRowsB,
            "Matrix incompatible A (" << nRowsA << " x " << nColsBlockA*numBits <<
                             " ) B (" << nRowsB << " x " << nColsBlockB*numBits << 
                             " ) " ) ;

    ASSERT( dimRowsA >= nRowsA, "GF2::MM4R nRowsA = " << nRowsA << " dimRowsA = " << dimRowsA ) ; 
    ASSERT( dimRowsB >= nRowsB, "GF2::MM4R nRowsB = " << nRowsB << " dimRowsB = " << dimRowsB ) ; 
    ASSERT( dimRowsC >= nRowsA, "GF2::MM4R nRowsA = " << nRowsA << " dimRowsC = " << dimRowsC ) ; 

    if ( nRowsA == 0 || nColsBlockB == 0 || nColsBlockA == 0 ) return ; // niente da moltiplicare

    // loop sui blocchi di B
    for ( unsigned colBlockB = 0 ; colBlockB < nColsBlockB ; ++colBlockB, B += dimRowsB, C += dimRowsC ) {
      MMass( A, nRowsA, B, C ) ;
      for ( unsigned colBlockA = 1 ; colBlockA < nColsBlockA ; ++colBlockA )
        MMadd( A + colBlockA * dimRowsA, nRowsA, B + colBlockA * numBits, C ) ;
    }
  }

  template
  void
  MMass( uint32_t const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
         uint32_t const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
         uint32_t       * C, unsigned dimRowsC ) ;

  template
  void
  MMass( uint64_t const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
         uint64_t const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
         uint64_t       * C, unsigned dimRowsC ) ;

  template
  void
  MMass( uint128_t const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
         uint128_t const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
         uint128_t       * C, unsigned dimRowsC ) ;

  //! C += A*B
  template <typename UNSIGNED>
  void
  MMadd( UNSIGNED const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
         UNSIGNED const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
         UNSIGNED       * C, unsigned dimRowsC ) {

    unsigned const numBits = sizeof(UNSIGNED) * CHAR_BIT ;

    ASSERT( nColsBlockA * numBits == nRowsB,
            "Matrix incompatible A (" << nRowsA << " x " << nColsBlockA*numBits <<
                             " ) B (" << nRowsB << " x " << nColsBlockB*numBits << 
                             " ) " ) ;

    ASSERT( dimRowsA >= nRowsA, "GF2::MM4R nRowsA = " << nRowsA << " dimRowsA = " << dimRowsA ) ; 
    ASSERT( dimRowsB >= nRowsB, "GF2::MM4R nRowsB = " << nRowsB << " dimRowsB = " << dimRowsB ) ; 
    ASSERT( dimRowsC >= nRowsA, "GF2::MM4R nRowsA = " << nRowsA << " dimRowsC = " << dimRowsC ) ; 

    if ( nRowsA == 0 || nColsBlockB == 0 || nColsBlockA == 0 ) return ; // niente da moltiplicare

    // loop sui blocchi di B
    if ( nColsBlockA < nColsBlockB ) {
      for ( unsigned colBlockA = 0 ; colBlockA < nColsBlockA ; ++colBlockA, A += dimRowsA, B += numBits )
        for ( unsigned colBlockB = 0 ; colBlockB < nColsBlockB ; ++colBlockB )
          MMadd( A, nRowsA, B + colBlockB * dimRowsB, C + colBlockB * dimRowsC ) ;
    } else {
      for ( unsigned colBlockB = 0 ; colBlockB < nColsBlockB ; ++colBlockB, B += dimRowsB, C += dimRowsC )
        for ( unsigned colBlockA = 0 ; colBlockA < nColsBlockA ; ++colBlockA )
          MMadd( A + colBlockA * dimRowsA, nRowsA, B + colBlockA * numBits, C ) ;
    }
  }

  template
  void
  MMadd( uint32_t const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
         uint32_t const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
         uint32_t       * C, unsigned dimRowsC ) ;

  template
  void
  MMadd( uint64_t const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
         uint64_t const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
         uint64_t       * C, unsigned dimRowsC ) ;

  template
  void
  MMadd( uint128_t const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
         uint128_t const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
         uint128_t       * C, unsigned dimRowsC ) ;

  //////////////////////////////////////////////////////////////////////////////

  //! C += A*B
  template <typename UNSIGNED>
  void
  MMaddRecurr( UNSIGNED const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
               UNSIGNED const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
               UNSIGNED       * C, unsigned dimRowsC ) {

    unsigned const numBits = sizeof(UNSIGNED) * CHAR_BIT ;

    ASSERT( nColsBlockA * numBits == nRowsB,
            "Matrix incompatible A (" << nRowsA << " x " << nColsBlockA*numBits <<
                             " ) B (" << nRowsB << " x " << nColsBlockB*numBits << 
                             " )\nnumBits = " << numBits ) ;

    ASSERT( dimRowsA >= nRowsA, "GF2::MM4R nRowsA = " << nRowsA << " dimRowsA = " << dimRowsA ) ; 
    ASSERT( dimRowsB >= nRowsB, "GF2::MM4R nRowsB = " << nRowsB << " dimRowsB = " << dimRowsB ) ; 
    ASSERT( dimRowsC >= nRowsA, "GF2::MM4R nRowsA = " << nRowsA << " dimRowsC = " << dimRowsC ) ; 

    if ( nRowsA == 0 || nColsBlockB == 0 || nColsBlockA == 0 ) return ; // niente da moltiplicare

    /*
    //  splitting tipo A
    //  +-------+-------+    +-------+ +-------+-------+
    //  |       |       |    |       | |       |       |
    //  |  C11  |  C12  | += |  A11  | |  B11  |  B12  |
    //  |       |       |    |       | |       |       |
    //  +-------+-------+    +-------+ +-------+-------+
    //
    //  splitting tipo B
    //  +-------+-------+    +-------+-------+ +-------+-------+
    //  |       |       |    |       |       | |       |       |
    //  |  C11  |  C12  | += |  A11  |  A12  | |  B11  |  B12  |
    //  |       |       |    |       |       | |       |       |
    //  +-------+-------+    +-------+-------+ +-------+-------+
    //                                         |       |       |
    //                                         |  B21  |  B22  |
    //                                         |       |       |
    //                                         +-------+-------+
    */
    
    // fine ricorsione
    if ( nColsBlockA < 2 || nColsBlockB < 2 ) {
      MMadd<UNSIGNED>( A, dimRowsA, nRowsA, nColsBlockA, 
                       B, dimRowsB, nRowsB, nColsBlockB,
                       C, dimRowsC ) ;
      return ;
    }

    // controllo il tipo di splitting
    unsigned const nColsBlockA1 = nColsBlockA>>1 ;
    unsigned const nColsBlockA2 = nColsBlockA - nColsBlockA1 ;

    unsigned const nColsBlockB1 = nColsBlockB>>1 ;
    unsigned const nColsBlockB2 = nColsBlockB - nColsBlockB1 ;

    unsigned const nRowsB1      = numBits * nColsBlockA1 ;
    unsigned const nRowsB2      = numBits * nColsBlockA2 ;

    UNSIGNED const * A11 = A ;
    UNSIGNED const * A12 = A + nColsBlockA1 * dimRowsA ;

    UNSIGNED const * B11 = B ;
    UNSIGNED const * B21 = B   + nRowsB1 ;
    UNSIGNED const * B12 = B   + nColsBlockB1 * dimRowsB ;
    UNSIGNED const * B22 = B12 + nRowsB1 ;

    UNSIGNED       * C11 = C ;
    UNSIGNED       * C12 = C + nColsBlockB1 * dimRowsC ;

    if ( nColsBlockA < nColsBlockB ) { // A
      // C11 += A11 * B11
      MMaddRecurr<UNSIGNED>( A11, dimRowsA, nRowsA, nColsBlockA,
                             B11, dimRowsB, nRowsB, nColsBlockB1,
                             C11, dimRowsC ) ;
      // C12 += A11 * B12
      MMaddRecurr<UNSIGNED>( A11, dimRowsA, nRowsA, nColsBlockA,
                             B12, dimRowsB, nRowsB, nColsBlockB2,
                             C12, dimRowsC ) ;
    } else { // tipo 2B
      // C11 += A11 * B11
      MMaddRecurr<UNSIGNED>( A11, dimRowsA, nRowsA,  nColsBlockA1,
                             B11, dimRowsB, nRowsB1, nColsBlockB1,
                             C11, dimRowsC ) ;
      // C11 += A12 * B21
      MMaddRecurr<UNSIGNED>( A12, dimRowsA, nRowsA,  nColsBlockA2,
                             B21, dimRowsB, nRowsB2, nColsBlockB1,
                             C11, dimRowsC ) ;
      // C12 += A11 * B12
      MMaddRecurr<UNSIGNED>( A11, dimRowsA, nRowsA,  nColsBlockA1,
                             B12, dimRowsB, nRowsB1, nColsBlockB2,
                             C12, dimRowsC ) ;
      // C12 += A12 * B22
      MMaddRecurr<UNSIGNED>( A12, dimRowsA, nRowsA,  nColsBlockA2,
                             B22, dimRowsB, nRowsB2, nColsBlockB2,
                             C12, dimRowsC ) ;
    }
  }

  template
  void
  MMaddRecurr( uint32_t const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
               uint32_t const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
               uint32_t       * C, unsigned dimRowsC ) ;

  template
  void
  MMaddRecurr( uint64_t const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
               uint64_t const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
               uint64_t       * C, unsigned dimRowsC ) ;

  template
  void
  MMaddRecurr( uint128_t const * A, unsigned dimRowsA, unsigned nRowsA, unsigned nColsBlockA,
               uint128_t const * B, unsigned dimRowsB, unsigned nRowsB, unsigned nColsBlockB,
               uint128_t       * C, unsigned dimRowsC ) ;

}

// EOF GF2toolkit_MatrixMult.cc
