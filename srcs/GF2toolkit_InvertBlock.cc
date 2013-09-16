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

#include "GF2toolkit_common.hh"
#include "GF2toolkit_popCounts.hh"
#include "GF2toolkit_InvertBlock.hh"

namespace GF2toolkit {

  #define SET(I)   uint32_t(1)<<(I)
  #define SET4(N)  SET(N), SET(N+1), SET(N+2), SET(N+3)
  #define SET16(N) SET4(N), SET4(N+4), SET4(N+8), SET4(N+12)

  uint32_t id32[32] = { SET16(0), SET16(16) } ;

  #undef SET
  #undef SET4
  #undef SET16

  #define SET(I)   uint64_t(1)<<(I)
  #define SET4(N)  SET(N), SET(N+1), SET(N+2), SET(N+3)
  #define SET16(N) SET4(N), SET4(N+4), SET4(N+8), SET4(N+12)

  uint64_t id64[64] = { SET16(0), SET16(16), SET16(32), SET16(48) } ;

  #undef SET
  #undef SET4
  #undef SET16

  static __m128i zero = _mm_setzero_si128() ;
  
  static
  inline
  uint128_t
  setBit128( unsigned i ) {
    uint128_t tmp ;
    if ( i < 64 ) {
      tmp.i64.lo = uint64_t(1)<<i ;
      tmp.i64.hi = 0 ;
    } else {
      tmp.i64.lo = 0 ;
      tmp.i64.hi = uint64_t(1)<<(i-64) ;
    }
    return tmp ;
  }

  #undef SET
  #undef SET4
  #undef SET16

  #define SET(I)   setBit128(I)
  #define SET4(N)  SET(N), SET(N+1), SET(N+2), SET(N+3)
  #define SET16(N) SET4(N), SET4(N+4), SET4(N+8), SET4(N+12)

  uint128_t id128[128] = {
    SET16(0),  SET16(16), SET16(32), SET16(48),
    SET16(64), SET16(80), SET16(96), SET16(112)
  } ;

  #undef SET
  #undef SET4
  #undef SET16
  
  /*
  //             _                       _____      _                  _   
  //    ___ ___ | |_   _ _ __ ___  _ __ | ____|_  _| |_ _ __ __ _  ___| |_ 
  //   / __/ _ \| | | | | '_ ` _ \| '_ \|  _| \ \/ / __| '__/ _` |/ __| __|
  //  | (_| (_) | | |_| | | | | | | | | | |___ >  <| |_| | | (_| | (__| |_ 
  //   \___\___/|_|\__,_|_| |_| |_|_| |_|_____/_/\_\\__|_|  \__,_|\___|\__|
  */
  template <typename UNSIGNED>
  void
  columnExtract( UNSIGNED const * Mmat, 
                 unsigned         kbit,
                 UNSIGNED       & xColumn ) ;

  template <>
  void
  columnExtract( uint32_t const * Mmat, 
                 unsigned         kbit,
                 uint32_t       & xColumn ) {
    // estrae primo pezzo di colonna
    xColumn = id32[kbit] ;
    for ( unsigned j = 0 ; j < kbit ; ++j )
      if ( Mmat[j] & id32[kbit] )
        xColumn |= id32[j] ;
  }

  template <>
  void
  columnExtract( uint64_t const * Mmat, 
                 unsigned         kbit,
                 uint64_t       & xColumn ) {
    // estrae primo pezzo di colonna
    xColumn = id64[kbit] ;
    for ( unsigned j = 0 ; j < kbit ; ++j )
      if ( Mmat[j] & id64[kbit] )
        xColumn |= id64[j] ;
  }

  template <>
  void
  columnExtract( uint128_t const * Mmat,
                 unsigned          kbit,
                 uint128_t       & xColumn ) {
    xColumn = id128[kbit] ;
    if ( kbit >= 64 ) {
      // estrae primo pezzo di colonna
      uint64_t const ik = id64[kbit-64] ;
      for ( unsigned j = 0 ; j < 64 ; ++j )
        if ( Mmat[j].i64.hi & ik )
          xColumn.i64.lo |= id64[j] ;
      for ( unsigned j = 64 ; j < kbit ; ++j )
        if ( Mmat[j].i64.hi & ik )
          xColumn.i64.hi |= id64[j-64] ;
    } else {
      // estrae primo pezzo di colonna
      uint64_t const ik = id64[kbit] ;
      for ( unsigned j = 0 ; j < kbit ; ++j )
        if ( Mmat[j].i64.lo & ik )
          xColumn.i64.lo |= id64[j] ;    
    }
  }

  /*
  //                           _     ____
  //   ___  ___  __ _ _ __ ___| |__ |  _ \ _____      __
  //  / __|/ _ \/ _` | '__/ __| '_ \| |_) / _ \ \ /\ / /
  //  \__ \  __/ (_| | | | (__| | | |  _ < (_) \ V  V /
  //  |___/\___|\__,_|_|  \___|_| |_|_| \_\___/ \_/\_/
  */
  template <typename UNSIGNED>
  inline
  int32_t
  searchRow( UNSIGNED const * Amat,
             unsigned         nRows, 
             unsigned         kbit,
             UNSIGNED const   xColumn ) ;

  template <>
  inline
  int32_t
  searchRow( uint32_t const * Amat, 
             unsigned         nRows, 
             unsigned         kbit,
             uint32_t const   xColumn ) {
    int32_t rj = 0 ;
    GF2_UNROLL_BY_16( nRows,
                      if ( *Amat != 0 && fastParity( *Amat & xColumn ) ) return rj ;
                      ++Amat ; 
                      ++rj ) ;
    return -1 ;
  }

  template <>
  inline
  int32_t
  searchRow( uint64_t const * Amat, 
             unsigned         nRows, 
             unsigned         kbit,
             uint64_t const   xColumn ) {
    int32_t rj = 0 ;
    GF2_UNROLL_BY_16( nRows,
                      if ( *Amat != 0 && fastParity( *Amat & xColumn ) ) return rj ;
                      ++Amat ; 
                      ++rj ) ;
    return -1 ;
  }

  template <>
  inline
  int32_t
  searchRow( uint128_t const * Amat, 
             unsigned          nRows, 
             unsigned          kbit,
             uint128_t const   xColumn ) {
    int32_t rj = 0 ;
    if ( kbit >= 64 ) {
      // cerco riga
      GF2_UNROLL_BY_8( nRows, 
                       if ( fastParity( _mm_and_si128( Amat->i128, xColumn.i128 ) ) ) return rj ;
                       ++rj ; ++Amat ) ; 
    } else {
      // cerco riga
      GF2_UNROLL_BY_8( nRows,
                       if ( //Amat->i64.lo != 0 && 
                                    fastParity( Amat->i64.lo & xColumn.i64.lo ) ) return rj ;
                       ++rj ; ++Amat ) ;
    }
    return -1 ;
  }

  /*
  //                   _       _       ____  _            _    
  //   _   _ _ __   __| | __ _| |_ ___| __ )| | ___   ___| | __
  //  | | | | '_ \ / _` |/ _` | __/ _ \  _ \| |/ _ \ / __| |/ /
  //  | |_| | |_) | (_| | (_| | ||  __/ |_) | | (_) | (__|   < 
  //   \__,_| .__/ \__,_|\__,_|\__\___|____/|_|\___/ \___|_|\_\
  //        |_|                                                
  //
  //  Multiply Q by a block k x numBits
  //
  //      +---+---+      +---+---+            +-------+-------+
  //      | I |   |      | I | x |            |I+x*y^T|   x   |
  //  L = +---+---+  U = +---+---+  T = U*L = +-------+-------+
  //      |y^T| 1 |      |   | 1 |            |  y^T  |1+y^T*x|
  //      +---+---+      +---+---+            +-------+-------+
  //
  //    +---+   +---+   +---+
  //    | A |   | A |   | x |
  //  T +---+ = +---+ + +---+ (y^T A + b^T)
  //    |b^T|   |0^T|   | 1 |
  //    +---+   +---+   +---+
  */
  template <typename UNSIGNED>
  void
  updateBlock( UNSIGNED *     Zmat,
               UNSIGNED *     Mmat,
               unsigned       kbit,
               UNSIGNED const xColumn ) {

    // y^T A + b^T
    UNSIGNED yRow = Mmat[kbit] ;
    for ( unsigned j = 0 ; j < kbit ; ++j )
      if ( (yRow>>j) & 0x01 ) {
        Zmat[kbit] ^= Zmat[j] ;
        Mmat[kbit] ^= Mmat[j]    ;
      }

    for ( unsigned j = 0 ; j < kbit ; ++j )
      if ( (xColumn>>j) & 0x01 ) {
        Zmat[j] ^= Zmat[kbit] ;
        Mmat[j] ^= Mmat[kbit]    ;
      }
  }

  template <>
  void
  updateBlock( uint128_t *     Zmat,
               uint128_t *     Mmat,
               unsigned        kbit,
               uint128_t const xColumn ) {

    // y^T A + b^T
    uint128_t yRow   = Mmat[kbit] ;
    unsigned  kbit64 = std::min( kbit, 64U ) ;
    __m128i & Zmatk  = Zmat[kbit].i128 ;
    __m128i & Mmatk  = Mmat[kbit].i128 ;

    // -------------------------------------------------------------------------
    for ( unsigned j = 0 ; j < kbit64 ; ++j )
      if ( (yRow.i64.lo>>j) & 0x01 ) {
        Zmatk = _mm_xor_si128( Zmatk, Zmat[j].i128 ) ;
        Mmatk = _mm_xor_si128( Mmatk, Mmat[j].i128 ) ;
      }
    for ( unsigned j = 64 ; j < kbit ; ++j )
      if ( (yRow.i64.hi>>(j-64)) & 0x01 ) {
        Zmatk = _mm_xor_si128( Zmatk, Zmat[j].i128 ) ;
        Mmatk = _mm_xor_si128( Mmatk, Mmat[j].i128 ) ;
      }

    // -------------------------------------------------------------------------
    for ( unsigned j = 0 ; j < kbit64 ; ++j )
      if ( (xColumn.i64.lo>>j) & 0x01 ) {
        Zmat[j].i128 = _mm_xor_si128( Zmatk, Zmat[j].i128 ) ;
        Mmat[j].i128 = _mm_xor_si128( Mmatk, Mmat[j].i128 ) ;
      }
    for ( unsigned j = 64 ; j < kbit ; ++j )
      if ( (xColumn.i64.hi>>(j-64)) & 0x01 ) {
        Zmat[j].i128 = _mm_xor_si128( Zmatk, Zmat[j].i128 ) ;
        Mmat[j].i128 = _mm_xor_si128( Mmatk, Mmat[j].i128 ) ;
      }
  }

  /*                                       _   _____                _   
  //   ___ ___  _ __ ___  _ __   __ _  ___| |_|__  /_ __ ___   __ _| |_ 
  //  / __/ _ \| '_ ` _ \| '_ \ / _` |/ __| __| / /| '_ ` _ \ / _` | __|
  // | (_| (_) | | | | | | |_) | (_| | (__| |_ / /_| | | | | | (_| | |_ 
  //  \___\___/|_| |_| |_| .__/ \__,_|\___|\__/____|_| |_| |_|\__,_|\__|
  //                     |_|
  */
  template <typename UNSIGNED>
  void
  compactZmat( UNSIGNED * Zmat, int32_t * Perm ) ;

  template <>
  void
  compactZmat( uint32_t * Zmat, int32_t * Perm ) {
    uint32_t m = 0 ;
    for ( unsigned i = 0 ; i < 32 ; ++i ) {
      if ( Perm[i] < 0 )
        for ( unsigned j = 0 ; j < 32 ; ++j )
          Zmat[j] = ( (Zmat[j]>>1) & (~m) ) | (Zmat[j]&m) ;
      else 
        m = (m<<1) | 0x1U ;
    }
  }

  template <>
  void
  compactZmat( uint64_t * Zmat, int32_t * Perm ) {
    uint64_t m = 0 ;
    for ( unsigned i = 0 ; i < 64 ; ++i ) {
      if ( Perm[i] < 0 )
        for ( unsigned j = 0 ; j < 64 ; ++j )
          Zmat[j] = ( (Zmat[j]>>1) & (~m) ) | (Zmat[j]&m) ;
      else 
        m = (m<<1) | 0x01UL ;
    }
  }

  template <>
  void
  compactZmat( uint128_t * Zmat, int32_t * Perm ) {
    uint128_t m ;
    m.i128 = _mm_setzero_si128() ;
    for ( unsigned i = 0 ; i < 128 ; ++i ) {
      if ( Perm[i] < 0 ) {
        for ( unsigned j = 0 ; j < 128 ; ++j ) {
          uint128_t tmp ;
          tmp.i64.lo = (Zmat[j].i64.lo>>1) | (Zmat[j].i64.hi<<63) ;
          tmp.i64.hi = Zmat[j].i64.hi>>1 ;
          __m128i t1 = _mm_andnot_si128( m.i128, tmp.i128 ) ;
          __m128i t2 = _mm_and_si128( m.i128, Zmat[j].i128 ) ;
          Zmat[j].i128 = _mm_or_si128( t1, t2 ) ;
        }
      } else {
        m.i64.hi = (m.i64.hi<<1) | (m.i64.lo>>63) ;
        m.i64.lo = (m.i64.lo<<1) | 0x01 ;
      }
    }
  }

  /*
  //   _                     _   ____  _            _    
  //  (_)_ ____   _____ _ __| |_| __ )| | ___   ___| | __
  //  | | '_ \ \ / / _ \ '__| __|  _ \| |/ _ \ / __| |/ /
  //  | | | | \ V /  __/ |  | |_| |_) | | (_) | (__|   < 
  //  |_|_| |_|\_/ \___|_|   \__|____/|_|\___/ \___|_|\_\
  */
  template <typename UNSIGNED>
  unsigned
  invertBlock( UNSIGNED * Amat, unsigned nRows,
               int32_t  * Perm,
               UNSIGNED * Zmat ) {

    uint32_t numUsedRows = 0 ;
    unsigned const numBits = CHAR_BIT*sizeof(UNSIGNED) ;
    UNSIGNED xColumn, Mmat[numBits] ;

    setIdentity( Zmat, numBits ) ;
    setIdentity( Mmat, numBits ) ;
    
    xColumn = Mmat[0] ; // = 0x01
    Perm[0] = searchRow<UNSIGNED>( Amat, nRows, 0, xColumn ) ;
    if ( Perm[0] >= 0 ) {
      std::swap( Amat[0], Amat[Perm[0]] ) ;
      Mmat[0] = Amat[0] ;
      ++numUsedRows ;
    }
    for ( unsigned kbit = 1 ; kbit < numBits ; ++kbit ) {
      columnExtract<UNSIGNED>( Mmat, kbit, xColumn ) ;
      Perm[kbit] = searchRow<UNSIGNED>( Amat+numUsedRows, nRows-numUsedRows, kbit, xColumn ) ;
      if ( Perm[kbit] >= 0 ) {
        Perm[kbit] += numUsedRows ; // shift della permutazione
        std::swap( Amat[numUsedRows], Amat[Perm[kbit]] ) ;
        Mmat[kbit] = Amat[numUsedRows] ;
        ++numUsedRows ;
      }
      updateBlock<UNSIGNED>( Zmat, Mmat, kbit, xColumn ) ;
    }
    compactZmat( Zmat, Perm ) ;
    return numUsedRows ; // numero righe usate dal blocco Amat
  }

  // instanziazione esplicita
  template
  unsigned
  invertBlock( uint32_t * Amat, unsigned nRows,
               int32_t  * Perm,
               uint32_t * AinvMmat ) ;

  template
  unsigned
  invertBlock( uint64_t * Amat, unsigned nRows,
               int32_t  * Perm,
               uint64_t * AinvMmat ) ;

  template
  unsigned
  invertBlock( uint128_t * Amat, unsigned nRows,
               int32_t   * Perm,
               uint128_t * AinvMmat ) ;

}

// EOF GF2toolkit_InvertBlock.cc
