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

#ifndef COMMON_HH
#define COMMON_HH

#include <cstdlib>
#include <climits>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

#ifdef __unix__
  // for memcpy on linux
  #include <string.h>
#endif

// blocco minimo ricorsione per inversione L
#define L_INV_APPLY_MIN_SIZE (1<<9)

#define STRASSEN_MIN_SIZE   (1<<9)
#define STRASSEN_MASK       0x0Fu

#ifdef MISSING_STDINT_H
  #include "pstdint.h"
#else
  #include <stdint.h>
#endif

#ifndef ASSERT
  #ifdef _WIN32
    #define ASSERT(COND,MSG) \
    if ( !(COND) ) { std::ostringstream ost ; ost << MSG << '\n' ; throw std::runtime_error(ost.str()) ; }
  #else
    #define ASSERT(COND,MSG) \
    if ( !(COND) ) { std::ostringstream ost ; ost << "File: " << __FILE__ << " line: " << __LINE__ << '\n' << MSG << '\n' ; throw std::runtime_error(ost.str()) ; }
  #endif
#endif

#ifdef UNSIGNED_IS_64BIT
  #define UINT16(A) A##u
  #define UINT32(A) A##ul
  #define UINT64(A) A##ull
#else
  #define UINT16(A) A##u
  #define UINT32(A) A##u
  #define UINT64(A) A##ull
#endif

#include <emmintrin.h>

typedef union uint128_t {
  __m128i i128 ;
  struct i64 { 
    uint64_t lo ;
    uint64_t hi ;
  } i64 ;
  struct i32 {
    uint32_t n0 ;
    uint32_t n1 ;
    uint32_t n2 ;
    uint32_t n3 ;
  } i32 ;
  struct i8 {
    uint8_t n0 ;
    uint8_t n1 ;
    uint8_t n2 ;
    uint8_t n3 ;
    uint8_t n4 ;
    uint8_t n5 ;
    uint8_t n6 ;
    uint8_t n7 ;
    uint8_t n8 ;
    uint8_t n9 ;
    uint8_t nA ;
    uint8_t nB ;
    uint8_t nC ;
    uint8_t nD ;
    uint8_t nE ;
    uint8_t nF ;
  } i8 ;
} uint128_t ;

/*
//                         _              _____     ____  _       
//   _ __  _   _ _ __ ___ | |__   ___ _ _|_   _|__ | __ )(_)_ __  
//  | '_ \| | | | '_ ` _ \| '_ \ / _ \ '__|| |/ _ \|  _ \| | '_ \ 
//  | | | | |_| | | | | | | |_) |  __/ |   | | (_) | |_) | | | | |
//  |_| |_|\__,_|_| |_| |_|_.__/ \___|_|   |_|\___/|____/|_|_| |_|
*/

template <typename T>
inline
std::string
numberToBin( T in ) {
  unsigned N = sizeof(T)*CHAR_BIT ;
  uint64_t v = in ;
  std::string res = "" ;
  for ( uint64_t mask = uint64_t(1)<<(N-1) ; mask != 0 ; mask >>= 1 )
    res += (v & mask) ? "1" : "0" ;
  return res ;
}

template <>
inline
std::string
numberToBin( uint128_t in ) {
  return numberToBin(in.i64.hi)+numberToBin(in.i64.lo) ;
}


template <typename T>
inline
std::string
numberToBinR( T in ) {
  unsigned N = sizeof(T)*CHAR_BIT ;
  uint64_t v = in ;
  std::string res = "" ;
  for ( uint64_t mask = uint64_t(1)<<(N-1) ; mask != 0 ; mask >>= 1 )
    res = ((v & mask) ? "1" : "0") + res ;
  return res ;
}

template <>
inline
std::string
numberToBinR( uint128_t in ) {
  return numberToBinR(in.i64.lo)+numberToBinR(in.i64.hi) ;
}

#define GF2_UNROLL_BY_32(N,OP) \
  { switch ( (N) & 0x1F ) {    \
    case 31: OP ;              \
    case 30: OP ;              \
    case 29: OP ;              \
    case 28: OP ;              \
    case 27: OP ;              \
    case 26: OP ;              \
    case 25: OP ;              \
    case 24: OP ;              \
    case 22: OP ;              \
    case 21: OP ;              \
    case 20: OP ;              \
    case 19: OP ;              \
    case 18: OP ;              \
    case 17: OP ;              \
    case 16: OP ;              \
    case 15: OP ;              \
    case 14: OP ;              \
    case 13: OP ;              \
    case 12: OP ;              \
    case 11: OP ;              \
    case 10: OP ;              \
    case 9:  OP ;              \
    case 8:  OP ;              \
    case 7:  OP ;              \
    case 6:  OP ;              \
    case 5:  OP ;              \
    case 4:  OP ;              \
    case 3:  OP ;              \
    case 2:  OP ;              \
    case 1:  OP ;              \
  }                            \
  unsigned idx = (N)>>5 ;      \
  while ( idx-- > 0 ) {        \
    OP ; OP ; OP ; OP ;        \
    OP ; OP ; OP ; OP ;        \
    OP ; OP ; OP ; OP ;        \
    OP ; OP ; OP ; OP ;        \
    OP ; OP ; OP ; OP ;        \
    OP ; OP ; OP ; OP ;        \
    OP ; OP ; OP ; OP ;        \
    OP ; OP ; OP ; OP ;        \
  }}


#define GF2_UNROLL_BY_16(N,OP) \
  { switch ( (N) & 0x0F ) {    \
    case 15: OP ;              \
    case 14: OP ;              \
    case 13: OP ;              \
    case 12: OP ;              \
    case 11: OP ;              \
    case 10: OP ;              \
    case 9:  OP ;              \
    case 8:  OP ;              \
    case 7:  OP ;              \
    case 6:  OP ;              \
    case 5:  OP ;              \
    case 4:  OP ;              \
    case 3:  OP ;              \
    case 2:  OP ;              \
    case 1:  OP ;              \
  }                            \
  unsigned idx = (N)>>4 ;      \
  while ( idx-- > 0 ) {        \
    OP ; OP ; OP ; OP ;        \
    OP ; OP ; OP ; OP ;        \
    OP ; OP ; OP ; OP ;        \
    OP ; OP ; OP ; OP ;        \
  }}

#define GF2_UNROLL_BY_8(N,OP) \
  { switch ( (N) & 0x07 ) {   \
    case 7:  OP ;             \
    case 6:  OP ;             \
    case 5:  OP ;             \
    case 4:  OP ;             \
    case 3:  OP ;             \
    case 2:  OP ;             \
    case 1:  OP ;             \
  }                           \
  unsigned idx = (N)>>3 ;     \
  while ( idx-- > 0 ) {       \
    OP ; OP ; OP ; OP ;       \
    OP ; OP ; OP ; OP ;       \
  }}

#define GF2_UNROLL_BY_4(N,OP) \
  { switch ( (N) & 0x03 ) {   \
    case 3:  OP ;             \
    case 2:  OP ;             \
    case 1:  OP ;             \
  }                           \
  unsigned idx = (N)>>2 ;     \
  while ( idx-- > 0 ) {       \
    OP ; OP ; OP ; OP ;       \
  }}

#if 1
#define SSE_ZERO(A)      (A)->i128 = _mm_setzero_si128()
#define SSE_XOR(A,B,C)   (A)->i128 = _mm_xor_si128( (B)->i128, (C)->i128 )
#define SSE_ASS_XOR(A,B) (A)->i128 = _mm_xor_si128( (A)->i128, (B)->i128 )
#else
#define SSE_ZERO(A)      A->i64.lo = 0 ; A->i64.hi = 0
#define SSE_XOR(A,B,C)   A->i64.lo = B->i64.lo^C->i64.lo ; A->i64.hi = B->i64.hi^C->i64.hi
#define SSE_ASS_XOR(A,B) A->i64.lo ^= B->i64.lo ; A->i64.hi ^= B->i64.hi
#endif

#endif

// EOF common.hh
