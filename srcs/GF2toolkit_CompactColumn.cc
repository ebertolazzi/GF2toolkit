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

#include "GF2toolkit_common.hh"                                                 
#include "GF2toolkit_InvertBlock.hh" 
#include "GF2toolkit_Matrix.hh" 

namespace GF2toolkit {

  /*
  //                                        _    ____      _                       
  //    ___ ___  _ __ ___  _ __   __ _  ___| |_ / ___|___ | |_   _ _ __ ___  _ __  
  //   / __/ _ \| '_ ` _ \| '_ \ / _` |/ __| __| |   / _ \| | | | | '_ ` _ \| '_ \ 
  //  | (_| (_) | | | | | | |_) | (_| | (__| |_| |__| (_) | | |_| | | | | | | | | |
  //   \___\___/|_| |_| |_| .__/ \__,_|\___|\__|\____\___/|_|\__,_|_| |_| |_|_| |_|
  //                      |_|                                                      
  */

  template <>
  void
  compactColumn( uint32_t * Lmat,     // colonna da compattare
                 unsigned   ncol,     // quante colonne contiene
                 uint32_t * Amat,     // colonna di arrivo
                 uint32_t * Amat1,    // colonna successiva
                 unsigned   startcol, // quante colonne sono gia occupate in Amat
                 unsigned   nRows ) {
    if ( startcol > 0 ) {
      unsigned nAvailable = 32-startcol ;
      if ( nAvailable >= ncol ) { // abbastanza spazio per una sola colonna
        GF2_UNROLL_BY_16( nRows, *Amat++ |= (*Lmat)<<startcol ; *Lmat++ = 0 ) ;    
      } else {
        GF2_UNROLL_BY_16( nRows, *Amat++ |= (*Lmat)<<startcol ;
                                 *Amat1++ = (*Lmat)>>nAvailable ;
                                 *Lmat++  = 0 ) ;
      }
    } else {
      copy( Lmat, nRows, Amat ) ;
      zero( Lmat, nRows ) ;
    }
  }
  
  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  compactColumn( uint32_t * Lmat,     // colonna da compattare
                 unsigned   ncol,     // quante colonne contiene
                 uint32_t * Amat,     // colonna di arrivo
                 unsigned   startcol, // quante colonne sono gia occupate in Amat
                 unsigned   nRows ) {
    if ( startcol > 0 ) {
      unsigned nAvailable = 32-startcol ;
      if ( nAvailable >= ncol ) { // abbastanza spazio per una sola colonna
        GF2_UNROLL_BY_16( nRows, *Amat++ |= (*Lmat)<<startcol ; *Lmat++ = 0 ) ;    
      } else {
        GF2_UNROLL_BY_16( nRows, *Amat++  |= (*Lmat)<<startcol ;
                                  *Lmat++ >>= nAvailable ) ;
      }
    } else {
      copy ( Lmat, nRows, Amat ) ;
      zero ( Lmat, nRows ) ;
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  compactColumn( uint64_t * Lmat,     // colonna da compattare
                 unsigned   ncol,     // quante colonne contiene
                 uint64_t * Amat,     // colonna di arrivo
                 uint64_t * Amat1,    // colonna successiva
                 unsigned   startcol, // quante colonne sono gia occupate in Amat
                 unsigned   nRows ) {
    if ( startcol > 0 ) {
      unsigned nAvailable = 64-startcol ;
      if ( nAvailable >= ncol ) { // abbastanza spazio per una sola colonna
        GF2_UNROLL_BY_16( nRows, *Amat++ |= (*Lmat)<<startcol ; *Lmat++ = 0 ) ;    
      } else {
        GF2_UNROLL_BY_16( nRows, *Amat++ |= (*Lmat)<<startcol ;
                                 *Amat1++ = (*Lmat)>>nAvailable ;
                                 *Lmat++  = 0 ) ;
      }
    } else {
      copy ( Lmat, nRows, Amat ) ;
      zero ( Lmat, nRows ) ;
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  compactColumn( uint64_t * Lmat,     // colonna da compattare
                 unsigned   ncol,     // quante colonne contiene
                 uint64_t * Amat,     // colonna di arrivo
                 unsigned   startcol, // quante colonne sono gia occupate in Amat
                 unsigned   nRows ) {
    if ( startcol > 0 ) {
      unsigned nAvailable = 64-startcol ;
      if ( nAvailable >= ncol ) { // abbastanza spazio per una sola colonna
        GF2_UNROLL_BY_16( nRows, *Amat++ |= (*Lmat)<<startcol ; *Lmat++ = 0 ) ;
      } else {
        GF2_UNROLL_BY_16( nRows, *Amat++  |= (*Lmat)<<startcol ; 
                                 *Lmat++ >>= nAvailable ) ;
      }
    } else {
      copy ( Lmat, nRows, Amat ) ;
      zero ( Lmat, nRows ) ;
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  compactColumn( uint128_t * Lmat,     // colonna da compattare
                 unsigned    ncol,     // quante colonne contiene
                 uint128_t * Amat,     // colonna di arrivo
                 uint128_t * Amat1,    // colonna successiva
                 unsigned    startcol, // quante colonne sono gia occupate in Amat
                 unsigned    nRows ) {

    if ( startcol == 0 ) {

      copy ( Lmat, nRows, Amat ) ;
      zero ( Lmat, nRows ) ;

    } else if ( startcol+ncol <= 64 ) {

      GF2_UNROLL_BY_16( nRows, Amat++->i64.lo |= Lmat->i64.lo<<startcol ;
                               Lmat++->i128    = _mm_setzero_si128() ) ;

    } else if ( startcol < 64 ) { // startcol+ncol >= 64

      unsigned const lshift = startcol ;
      unsigned const rshift = 64-startcol ;

      if ( ncol <= 64 ) {

        GF2_UNROLL_BY_16( nRows, Amat->i64.lo  |= Lmat->i64.lo<<lshift ;
                                 Amat++->i64.hi = Lmat->i64.lo>>rshift ;
                                 Lmat++->i64.lo = 0 ) ;

      } else if ( startcol+ncol <= 128 ) { // ncol > 64

        GF2_UNROLL_BY_16( nRows, 
                          Amat->i64.lo  |= Lmat->i64.lo<<lshift ;
                          Amat++->i64.hi = (Lmat->i64.lo>>rshift) | (Lmat->i64.hi<<lshift) ;
                          Lmat++->i128   = _mm_setzero_si128() ) ;
      
      } else { // ncol > 64 & startcol+ncol > 128

        GF2_UNROLL_BY_16( nRows,
                          Amat->i64.lo   |= Lmat->i64.lo<<lshift ;
                          Amat++->i64.hi  = (Lmat->i64.lo>>rshift) | (Lmat->i64.hi<<lshift) ;
                          Amat1->i64.lo   = Lmat->i64.hi>>rshift ;
                          Amat1++->i64.hi = 0 ;
                          Lmat++->i128    = _mm_setzero_si128() ) ;
        
      }

    } else if ( startcol == 64 ) { // caso particolare
      
      if ( ncol <= 64 ) { // sta tutto nella parte alta
        GF2_UNROLL_BY_16( nRows, Amat++->i64.hi = Lmat->i64.lo ; Lmat++->i64.lo = 0 ) ;
      } else {
        GF2_UNROLL_BY_16( nRows, Amat++->i64.hi  = Lmat->i64.lo ;
                                 Amat1++->i64.lo = Lmat->i64.hi ;
                                 Lmat++->i64.hi  = 0 ) ;
      }

    } else { // startcol > 64

      unsigned const rshift = 128-startcol ;
      unsigned const lshift = startcol-64 ;

      if ( ncol <= 64 ) {

        if ( startcol+ncol <= 128 ) { // sta tutto nella parte alta

          GF2_UNROLL_BY_16( nRows, Amat++->i64.hi |= Lmat->i64.lo<<lshift ;
                                   Lmat++->i64.lo  = 0 ) ;

        } else {

          GF2_UNROLL_BY_16( nRows, Amat->i64.hi  |= Lmat->i64.lo<<lshift ;
                                   Amat++->i64.lo = Lmat->i64.lo>>rshift ;
                                   Lmat++->i64.lo = 0 ) ;

        }

      } else { // ncol > 64 e startcol >= 64

        GF2_UNROLL_BY_16( nRows,
                          Amat++->i64.hi |= Lmat->i64.lo<<lshift ;
                          Amat1->i64.lo   = (Lmat->i64.lo>>rshift) | (Lmat->i64.hi<<lshift) ;
                          Amat1++->i64.hi = Lmat->i64.hi>>rshift ;
                          Lmat++->i128    = _mm_setzero_si128() ) ;
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  template <>
  void
  compactColumn( uint128_t * Lmat,     // colonna da compattare
                 unsigned    ncol,     // quante colonne contiene
                 uint128_t * Amat,     // colonna di arrivo
                 unsigned    startcol, // quante colonne sono gia occupate in Amat
                 unsigned    nRows ) {
    if ( startcol == 0 ) {

      copy ( Lmat, nRows, Amat ) ;
      zero ( Lmat, nRows ) ;

    } else if ( startcol+ncol <= 64 ) {

      GF2_UNROLL_BY_16( nRows, Amat++->i64.lo |= Lmat->i64.lo<<startcol ;
                               Lmat++->i128    = _mm_setzero_si128() ) ;

    } else if ( startcol < 64 ) { // startcol+ncol >= 64

      unsigned const lshift = startcol ;
      unsigned const rshift = 64-startcol ;

      if ( ncol <= 64 ) {

        GF2_UNROLL_BY_16( nRows, Amat->i64.lo  |= Lmat->i64.lo<<lshift ;
                                 Amat++->i64.hi = Lmat->i64.lo>>rshift ;
                                 Lmat++->i64.lo = 0 ) ;

      } else if ( startcol+ncol <= 128 ) { // ncol > 64

        GF2_UNROLL_BY_16( nRows, 
                          Amat->i64.lo  |= Lmat->i64.lo<<lshift ;
                          Amat++->i64.hi = (Lmat->i64.lo>>rshift) | (Lmat->i64.hi<<lshift) ;
                          Lmat++->i128   = _mm_setzero_si128() ) ;
      
      } else { // ncol > 64 & startcol+ncol > 128

        GF2_UNROLL_BY_16( nRows,
                          Amat->i64.lo   |= Lmat->i64.lo<<lshift ;
                          Amat++->i64.hi  = (Lmat->i64.lo>>rshift) | (Lmat->i64.hi<<lshift) ;
                          Lmat->i64.lo    = Lmat->i64.hi>>rshift ;
                          Lmat++->i64.hi  = 0 ) ;
        
      }

    } else if ( startcol == 64 ) { // caso particolare
      
      if ( ncol <= 64 ) { // sta tutto nella parte alta
        GF2_UNROLL_BY_16( nRows, Amat++->i64.hi = Lmat->i64.lo ; Lmat++->i64.lo = 0 ) ;
      } else {
        GF2_UNROLL_BY_16( nRows, Amat++->i64.hi = Lmat->i64.lo ;
                                 Lmat->i64.lo   = Lmat->i64.hi ;
                                 Lmat++->i64.hi = 0 ) ;
      }

    } else { // startcol > 64

      unsigned const rshift = 128-startcol ;
      unsigned const lshift = startcol-64 ;

      if ( ncol <= 64 ) {

        if ( startcol+ncol <= 128 ) { // sta tutto nella parte alta

          GF2_UNROLL_BY_16( nRows, Amat++->i64.hi |= Lmat->i64.lo<<lshift ;
                                   Lmat++->i64.lo  = 0 ) ;

        } else {

          GF2_UNROLL_BY_16( nRows, Amat->i64.hi  |= Lmat->i64.lo<<lshift ;
                                   Amat++->i64.lo = Lmat->i64.lo>>rshift ;
                                   Lmat++->i64.lo = 0 ) ;

        }

      } else { // ncol > 64 e startcol >= 64

        GF2_UNROLL_BY_16( nRows,
                          Amat++->i64.hi  |= Lmat->i64.lo<<lshift ;
                          Lmat->i64.lo     = (Lmat->i64.lo>>rshift) | (Lmat->i64.hi<<lshift) ;
                          Lmat++->i64.hi >>= rshift ) ;
      }
    }
  }
}

// EOF GF2toolkit_CompactColumn.cc
