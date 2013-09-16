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
#include "GF2toolkit_MatrixLower.hh"
#include "GF2toolkit_InvertBlock.hh"
#include "GF2toolkit_4Russian.hh"

namespace GF2toolkit {

  using namespace std ;
  
  /*
  //   _     _                              
  //  | |   (_)_ ____   _____ _ __ ___  ___ 
  //  | |   | | '_ \ \ / / _ \ '__/ __|/ _ \
  //  | |___| | | | \ V /  __/ |  \__ \  __/
  //  |_____|_|_| |_|\_/ \___|_|  |___/\___|
  */
  
  /*
  //
  //  +---+---+ +---+---+   +---+---+
  //  | L | 0 | | iL| 0 |   | I | 0 |
  //  +---+---+ +---+---+ = +---+---+  -> z = -w iL
  //  | w | 1 | | z | 1 |   | 0 | 1 |
  //  +---+---+ +---+---+   +---+---+
  */
  
  /*
  //   _________    _     _ _   
  //  |___ /___ \  | |__ (_) |_ 
  //    |_ \ __) | | '_ \| | __|
  //   ___) / __/  | |_) | | |_ 
  //  |____/_____| |_.__/|_|\__|
  */
  template <>
  void
  Linverse ( uint32_t const * L, uint32_t * Linv ) {
    setIdentity( Linv, 32 ) ;
    for ( unsigned k = 1 ; k < 32 ; ++k ) {
      uint32_t   tmp   = L[k] ;
      uint32_t & Linvk = Linv[k] ;
      uint32_t * pLinv = Linv ;
      GF2_UNROLL_BY_8( k,
                       if ( tmp & 0x01 ) Linvk ^= *pLinv ;
                       ++pLinv ;
                       tmp >>= 1 ) ;
      //for ( unsigned j = 0 ; j < k ; ++j, tmp >>= 1 )
      //  if ( tmp & 0x01 )
      //    Linvk ^= Linv[j] ;
    }
  }

  template <>
  void
  Linverse ( uint32_t const * L, uint32_t * Linv, unsigned nrow ) {
    setIdentity( Linv, nrow ) ;
    for ( unsigned k = 1 ; k < nrow ; ++k ) {
      uint32_t   tmp   = L[k] ;
      uint32_t & Linvk = Linv[k] ;
      uint32_t * pLinv = Linv ;
      GF2_UNROLL_BY_8( k,
                       if ( tmp & 0x01 ) Linvk ^= *pLinv ;
                       ++pLinv ;
                       tmp >>= 1 ) ;
      //for ( unsigned j = 0 ; j < k ; ++j, tmp >>= 1 )
      //  if ( tmp & 0x01 )
      //    Linvk ^= Linv[j] ;
    }
  }
  
  /*
  //    __   _  _     _     _ _   
  //   / /_ | || |   | |__ (_) |_ 
  //  | '_ \| || |_  | '_ \| | __|
  //  | (_) |__   _| | |_) | | |_ 
  //   \___/   |_|   |_.__/|_|\__|
  */
  template <>
  void
  Linverse ( uint64_t const * L, uint64_t * Linv ) {
    setIdentity( Linv, 64 ) ;
    for ( unsigned k = 1 ; k < 64 ; ++k ) {
      uint64_t   tmp   = L[k] ;
      uint64_t & Linvk = Linv[k] ;
      uint64_t * pLinv = Linv ;
      GF2_UNROLL_BY_16( k,
                        if ( tmp & 0x01 ) Linvk ^= *pLinv ;
                        ++pLinv ;
                        tmp >>= 1 ) ;
      //for ( unsigned j = 0 ; j < k ; ++j, tmp >>= 1 )
      //  if ( tmp & 0x01 )
      //    Linvk ^= Linv[j] ;
    }
  }

  template <>
  void
  Linverse ( uint64_t const * L, uint64_t * Linv, unsigned nrow ) {
    setIdentity( Linv, nrow ) ;
    for ( unsigned k = 1 ; k < nrow ; ++k ) {
      uint64_t   tmp   = L[k] ;
      uint64_t & Linvk = Linv[k] ;
      uint64_t * pLinv = Linv ;
      GF2_UNROLL_BY_16( k,
                        if ( tmp & 0x01 ) Linvk ^= *pLinv ;
                        ++pLinv ;
                        tmp >>= 1 ) ;
      //for ( unsigned j = 0 ; j < k ; ++j, tmp >>= 1 )
      //  if ( tmp & 0x01 )
      //    Linvk ^= Linv[j] ;
    }
  }
  
  /*
  //   _ ____  ___    _     _ _   
  //  / |___ \( _ )  | |__ (_) |_ 
  //  | | __) / _ \  | '_ \| | __|
  //  | |/ __/ (_) | | |_) | | |_ 
  //  |_|_____\___/  |_.__/|_|\__|
  */

  template <>
  void
  Linverse ( uint128_t const * L, uint128_t * Linv ) {
    setIdentity( Linv, 128 ) ;

    for ( unsigned k = 1 ; k < 64 ; ++k ) {
      uint64_t    tmp   = L[k].i64.lo ;
      uint64_t  & Linvk = Linv[k].i64.lo ;
      uint128_t * pLinv = Linv ;
      GF2_UNROLL_BY_16( k,
                        if ( tmp & 0x01 ) Linvk ^= pLinv->i64.lo ;
                        ++pLinv ;
                        tmp >>= 1 ) ;
      //for ( unsigned j = 0 ; j < k ; ++j, tmp >>= 1 )
      //  if ( tmp & 0x01 )
      //    Linvk ^= Linv[j].i64.lo ;
    }

    for ( unsigned k = 64 ; k < 128 ; ++k ) {
      uint64_t    tmp   = L[k].i64.lo ;
      uint64_t  & Linvk = Linv[k].i64.lo ;
      uint128_t * pLinv = Linv ;
      for ( unsigned i = 0 ; i < 8 ; ++i ) {
        if ( tmp & 0x01 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x02 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x04 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x08 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x10 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x20 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x40 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x80 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        tmp >>= 8 ;
      }        
      //for ( unsigned j = 0 ; j < 64 ; ++j, tmp >>= 1 )
      //  if ( tmp & 0x01 )
      //    Linvk ^= Linv[j].i64.lo ;
      tmp = L[k].i64.hi ;
      __m128i & Lk = Linv[k].i128 ;
      GF2_UNROLL_BY_16( k-64, if ( tmp & 0x01 ) Lk = _mm_xor_si128( Lk, pLinv->i128 ) ; 
                              ++pLinv ;
                              tmp >>= 1 ) ;
      //for ( unsigned j = 64 ; j < k ; ++j, tmp >>= 1 )
      //  if ( tmp & 0x01 )
      //    Lk = _mm_xor_si128( Lk, Linv[j].i128  ) ;
    }
  }

  template <>
  void
  Linverse ( uint128_t const * L, uint128_t * Linv, unsigned nrow ) {
    setIdentity( Linv, nrow ) ;

    unsigned NROW64 = nrow < 64 ? nrow : 64 ;
    for ( unsigned k = 1 ; k < NROW64 ; ++k ) {
      uint64_t    tmp   = L[k].i64.lo ;
      uint64_t  & Linvk = Linv[k].i64.lo ;
      uint128_t * pLinv = Linv ;
      GF2_UNROLL_BY_16( k,
                        if ( tmp & 0x01 ) Linvk ^= pLinv->i64.lo ;
                        ++pLinv ;
                        tmp >>= 1 ) ;
      //for ( unsigned j = 0 ; j < k ; ++j, tmp >>= 1 )
      //  if ( tmp & 0x01 )
      //    Linvk ^= Linv[j].i64.lo ;
    }

    for ( unsigned k = 64 ; k < nrow ; ++k ) {
      uint64_t    tmp   = L[k].i64.lo ;
      uint64_t  & Linvk = Linv[k].i64.lo ;
      uint128_t * pLinv = Linv ;
      for ( unsigned i = 0 ; i < 8 ; ++i ) {
        if ( tmp & 0x01 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x02 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x04 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x08 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x10 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x20 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x40 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        if ( tmp & 0x80 ) Linvk ^= pLinv->i64.lo ; ++pLinv ;
        tmp >>= 8 ;
      }        
      //for ( unsigned j = 0 ; j < 64 ; ++j, tmp >>= 1 )
      //  if ( tmp & 0x01 )
      //   Linvk ^= Linv[j].i64.lo ;

      tmp = L[k].i64.hi ;
      __m128i & Lk = Linv[k].i128 ;
      GF2_UNROLL_BY_16( k-64,
                        if ( tmp & 0x01 ) Lk = _mm_xor_si128( Lk, pLinv->i128 ) ; 
                        ++pLinv ;
                        tmp >>= 1 ) ;
      //for ( unsigned j = 64 ; j < k ; ++j, tmp >>= 1 )
      //  if ( tmp & 0x01 )
      //    Lk = _mm_xor_si128( Lk, Linv[j].i128  ) ;
    }
  }
}

// EOF GF2toolkit_MatrixLower.cc
