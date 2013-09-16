/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2010                                                      |
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
 |      version: 0.2 29-03-2011                                             |
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

// http://graphics.stanford.edu/~seander/bithacks.html

#ifndef GF2TOOLKIT_POP_COUNT_HH
#define GF2TOOLKIT_POP_COUNT_HH

#include "GF2toolkit_common.hh"

///////////////////////////////////////////////////////////////////////////

// http://graphics.stanford.edu/~seander/bithacks.html
// http://www.azillionmonkeys.com/qed/
// http://www.dspguru.com/dsp/tricks/gray-code-conversion
// http://chessprogramming.wikispaces.com/
// http://home.pipeline.com/~hbaker1/hakmem/hakmem.html

// BIT_COUNT_PARALLEL  7.1 sec (12.4 sec 32 bit)
// BIT_COUNT_HAKMEM    8.1 sec (14.7 sec 32 bit)
// BIT_COUNT_STANFORD  5   sec (9.6  sec 32 bit)
// BIT_COUNT_GCC       7   sec (11.6 sec 32 bit)

// if undefined set default 64 bit and Stanford bit count

#if !defined(BIT_COUNT_PARALLEL) && \
    !defined(BIT_COUNT_HAKMEM)   && \
    !defined(BIT_COUNT_STANFORD) && \
    !defined(BIT_COUNT_GCC)
  #define BIT_COUNT_STANFORD
#endif

/*
//                      ____                  _       
//   _ __   ___  _ __  / ___|___  _   _ _ __ | |_ ___ 
//  | '_ \ / _ \| '_ \| |   / _ \| | | | '_ \| __/ __|
//  | |_) | (_) | |_) | |__| (_) | |_| | | | | |_\__ \
//  | .__/ \___/| .__/ \____\___/ \__,_|_| |_|\__|___/
//  |_|         |_|                                   
*/  

/*
//                      _ _     _ 
//  _ __  __ _ _ _ __ _| | |___| |
// | '_ \/ _` | '_/ _` | | / -_) |
// | .__/\__,_|_| \__,_|_|_\___|_|
// |_|                            
*/

extern unsigned popCountTable8bit[] ;
extern uint8_t  parityTable8bit[] ;

static
inline
unsigned
popCounts_parallel( uint16_t n ) {
  n = (n&UINT16(0x5555))+((n>>1)&UINT16(0x5555)) ;
  n = (n&UINT16(0x3333))+((n>>2)&UINT16(0x3333)) ;
  n = (n&UINT16(0x0F0F))+((n>>4)&UINT16(0x0F0F)) ;
  n = (n&UINT16(0x00FF))+((n>>8)&UINT16(0x00FF)) ;
  return unsigned(n) ;
}

static
inline
unsigned
popCounts_parallel( uint32_t n ) {
  n = (n&UINT32(0x55555555))+((n>>1)&UINT32(0x55555555)) ;
  n = (n&UINT32(0x33333333))+((n>>2)&UINT32(0x33333333)) ;
  n = (n&UINT32(0x0F0F0F0F))+((n>>4)&UINT32(0x0F0F0F0F)) ;
  n = (n&UINT32(0x00FF00FF))+((n>>8)&UINT32(0x00FF00FF)) ;
  n = (n&UINT32(0x0000FFFF))+((n>>16)&UINT32(0x0000FFFF)) ;
  return unsigned(n) ;
}

static
inline
unsigned
popCounts_parallel( uint64_t n ) {
  n = (n&UINT64(0x5555555555555555))+( (n>>1)&UINT64(0x5555555555555555)) ;
  n = (n&UINT64(0x3333333333333333))+( (n>>2)&UINT64(0x3333333333333333)) ;
  n = (n&UINT64(0x0F0F0F0F0F0F0F0F))+( (n>>4)&UINT64(0x0F0F0F0F0F0F0F0F)) ;
  n = (n&UINT64(0x00FF00FF00FF00FF))+( (n>>8)&UINT64(0x00FF00FF00FF00FF)) ;
  n = (n&UINT64(0x0000FFFF0000FFFF))+((n>>16)&UINT64(0x0000FFFF0000FFFF)) ;
  n = (n&UINT64(0x00000000FFFFFFFF))+((n>>32)&UINT64(0x00000000FFFFFFFF)) ;
  return unsigned(n) ;
}

static
inline
unsigned
popCounts_parallel( uint128_t n ) {
  return popCounts_parallel(n.i64.lo) +
         popCounts_parallel(n.i64.hi) ;
}

/*
//  _         _                   
// | |_  __ _| |___ __  ___ _ __  
// | ' \/ _` | / / '  \/ -_) '  \ 
// |_||_\__,_|_\_\_|_|_\___|_|_|_|
*/

static
inline
unsigned
popCounts_hakmem( uint16_t n ) {
  n -= ((n >> 1) & UINT32(0133333)) + 
       ((n >> 2) & UINT32(0111111)) ;
  n += n >> 3 ;
  return ( n & UINT32(070707) ) % 63 ;
}

static
inline
unsigned
popCounts_hakmem( uint32_t n ) {
  n -= ((n >> 1) & UINT32(033333333333)) + 
       ((n >> 2) & UINT32(011111111111)) ;
  n += n >> 3 ;
  return ( n & UINT32(030707070707) ) % 63 ;
}

static
inline
unsigned
popCounts_hakmem( uint64_t n ) {
  n -= ((n >> 1) & UINT64(0333333333333333333333)) + 
       ((n >> 2) & UINT64(0111111111111111111111)) ;
  n += n >> 3 ;
  return unsigned((n & UINT64(0707070707070707070700)) % 63) + unsigned(n & 7);
}

static
inline
unsigned
popCounts_hakmem( uint128_t n ) {
  return popCounts_hakmem(n.i64.lo) +
         popCounts_hakmem(n.i64.hi) ;
}

/*
//     _              __            _ 
//  __| |_ __ _ _ _  / _|___ _ _ __| |
// (_-<  _/ _` | ' \|  _/ _ \ '_/ _` |
// /__/\__\__,_|_||_|_| \___/_| \__,_|
*/

static
inline
unsigned
popCounts_stanford( uint16_t n ) {
  n = (n & UINT16(0x5555)) + ((n>>1) & UINT16(0x5555));
  n = (n & UINT16(0x3333)) + ((n>>2) & UINT16(0x3333));
  n = (n & UINT16(0x0F0F)) + ((n>>4) & UINT16(0x0F0F));
  return ((n * UINT16(0x0101)) >> 8) & 0xFF ;
}

static
inline
unsigned
popCounts_stanford( uint32_t n ) {
  n  = (n & UINT32(0x55555555)) + ((n>>1) & UINT32(0x55555555));
  n  = (n & UINT32(0x33333333)) + ((n>>2) & UINT32(0x33333333));
  n  = (n & UINT32(0x0F0F0F0F)) + ((n>>4) & UINT32(0x0F0F0F0F));
  n *= UINT32(0x01010101) ;
  return n>>24 ;
}

static
inline
unsigned
popCounts_stanford( uint64_t n ) {
  n -= (n>>1) & UINT64(0x5555555555555555) ;
  n  = (n & UINT64(0x3333333333333333)) + ((n>>2) & UINT64(0x3333333333333333));
  n  = ((n>>4)+n) & UINT64(0x0F0F0F0F0F0F0F0F) ;
  n *= UINT64(0x0101010101010101) ;
  return n>>56 ;
}

static
inline
unsigned
popCounts_stanford( uint128_t n ) {
  return popCounts_stanford(n.i64.lo) +
         popCounts_stanford(n.i64.hi) ;
}

/*
//   _        _    _     
//  | |_ __ _| |__| |___ 
//  |  _/ _` | '_ \ / -_)
//   \__\__,_|_.__/_\___|
*/

static
inline
unsigned
popCounts_table( uint16_t n ) {
  return popCountTable8bit[n&0xFF]+
         popCountTable8bit[n>>8] ;
}

static
inline
unsigned
popCounts_table( uint32_t n ) {
  return popCountTable8bit[n&0xFF]+
         popCountTable8bit[(n>>8)&0xFF]+
         popCountTable8bit[(n>>16)&0xFF]+
         popCountTable8bit[(n>>24)&0xFF] ;
}

static
inline
unsigned
popCounts_table( uint64_t n ) {
  return popCountTable8bit[n&0xFF]+
         popCountTable8bit[(n>>8)&0xFF]+
         popCountTable8bit[(n>>16)&0xFF]+
         popCountTable8bit[(n>>24)&0xFF]+
         popCountTable8bit[(n>>32)&0xFF]+
         popCountTable8bit[(n>>40)&0xFF]+
         popCountTable8bit[(n>>48)&0xFF]+
         popCountTable8bit[(n>>56)&0xFF] ;
}

static
inline
unsigned
popCounts_table( uint128_t n ) {
  return popCounts_table(n.i64.lo) +
         popCounts_table(n.i64.hi) ;
}

template <typename Type> static inline unsigned popCounts( Type n ) ;

////////////////////////////////////////////////////////////

template <>
inline unsigned popCounts( uint16_t n ) { return popCounts_table(n) ; }

template <>
inline unsigned popCounts( uint32_t n ) { return popCounts_table(n) ; }

template <>
inline unsigned popCounts( uint64_t n ) { return popCounts_stanford(n) ; }

template <>
inline
unsigned
popCounts( __m128i b ) {

  static const unsigned mu1 = 0x55555555 ;
  static const unsigned mu2 = 0x33333333 ;
  static const unsigned mu3 = 0x0F0F0F0F ;
  static const unsigned mu4 = 0x0000003F ;

  // Loading masks
  __m128i m1 = _mm_set_epi32 (mu1, mu1, mu1, mu1);
  __m128i m2 = _mm_set_epi32 (mu2, mu2, mu2, mu2);
  __m128i m3 = _mm_set_epi32 (mu3, mu3, mu3, mu3);
  __m128i m4 = _mm_set_epi32 (mu4, mu4, mu4, mu4);
  __m128i mcnt ;

  mcnt = _mm_setzero_si128(); // cnt = 0

  __m128i tmp1, tmp2 ;

  // b = (b & 0x55555555) + (b >> 1 & 0x55555555);
  tmp1 = _mm_srli_epi32(b, 1);                    // tmp1 = ((b>>1) & 0x55555555)
  tmp1 = _mm_and_si128(tmp1, m1); 
  tmp2 = _mm_and_si128(b, m1);                    // tmp2 = (b&0x55555555)
  b    = _mm_add_epi32(tmp1, tmp2);               // b    = tmp1 + tmp2

  // b = (b & 0x33333333) + (b >> 2 & 0x33333333);
  tmp1 = _mm_srli_epi32(b, 2);                    // (b >> 2 & 0x33333333)
  tmp1 = _mm_and_si128(tmp1, m2); 
  tmp2 = _mm_and_si128(b, m2);                    // (b & 0x33333333)
  b    = _mm_add_epi32(tmp1, tmp2);               // b = tmp1 + tmp2

  // b = (b + (b >> 4)) & 0x0F0F0F0F;
  tmp1 = _mm_srli_epi32(b, 4);                    // tmp1 = b >> 4
  b    = _mm_add_epi32(b, tmp1);                  // b = b + (b >> 4)
  b    = _mm_and_si128(b, m3);                    //     & 0x0F0F0F0F

  // b = b + (b >> 8);
  tmp1 = _mm_srli_epi32 (b, 8);                   // tmp1 = b >> 8
  b    = _mm_add_epi32(b, tmp1);                  // b = b + (b >> 8)

  // b = (b + (b >> 16)) & 0x0000003F;
  tmp1 = _mm_srli_epi32 (b, 16);                  // b >> 16
  b    = _mm_add_epi32(b, tmp1);                  // b + (b >> 16)
  b    = _mm_and_si128(b, m4);                    // (b >> 16) & 0x0000003F;

  uint128_t tcnt ;
  tcnt.i128 = _mm_add_epi32(mcnt, b);                  // mcnt += b

  return tcnt.i32.n0 + tcnt.i32.n1 + tcnt.i32.n2 + tcnt.i32.n3 ;
}

template <>
inline
unsigned
popCounts( uint128_t n ) 
{ return popCounts( n.i128 ) ; }

/*
//    __           _   ____            _ _         
//   / _| __ _ ___| |_|  _ \ __ _ _ __(_) |_ _   _ 
//  | |_ / _` / __| __| |_) / _` | '__| | __| | | |
//  |  _| (_| \__ \ |_|  __/ (_| | |  | | |_| |_| |
//  |_|  \__,_|___/\__|_|   \__,_|_|  |_|\__|\__, |
//                                           |___/ 
*/
static
inline
unsigned
fastParity( uint16_t n ) {
  return parityTable8bit[(n^(n>>8))&0xFF] ;
}

////////////////////////////////////////////////////////////

static
inline
unsigned
fastParity( uint32_t n ) {
  n ^= n>>16 ;
  return parityTable8bit[(n^(n>>8))&0xFF] ;
}

////////////////////////////////////////////////////////////

static
inline
unsigned
fastParity( uint64_t n ) {
  if ( n ) {
    n ^= n>>32 ;
    n ^= n>>16 ;
    return parityTable8bit[(n^(n>>8))&0xFF] ;
  } else {
    return 0 ;
  }
}

static
inline
unsigned
fastParity( uint128_t const & n ) {
  return fastParity(n.i64.lo^n.i64.hi) ;
}

static
inline
unsigned
fastParity( __m128i const & nn ) {
  uint64_t const * nv = (uint64_t const *) &nn ;
  return fastParity(nv[0]^nv[1]) ;
}

#endif

// EOF GF2_popCounts.hh
