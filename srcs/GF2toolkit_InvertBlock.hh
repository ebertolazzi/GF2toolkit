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
 |      version: 0.2 21-01-2010                                             |
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

#ifndef GF2TOOLKIT_INVERT_BLOCK_HH
#define GF2TOOLKIT_INVERT_BLOCK_HH

#include "GF2toolkit_common.hh"

// for memcpy
#include <string.h>

namespace GF2toolkit {

  extern uint32_t  id32[32]   ;
  extern uint64_t  id64[64]   ;
  extern uint128_t id128[128] ;

  /*
  //            _   ___    _            _   _ _         
  //   ___  ___| |_|_ _|__| | ___ _ __ | |_(_) |_ _   _ 
  //  / __|/ _ \ __|| |/ _` |/ _ \ '_ \| __| | __| | | |
  //  \__ \  __/ |_ | | (_| |  __/ | | | |_| | |_| |_| |
  //  |___/\___|\__|___\__,_|\___|_| |_|\__|_|\__|\__, |
  //                                               |___/
  */
  template <typename UNSIGNED>
  inline
  void
  setIdentity( UNSIGNED * AinvMmat, unsigned nr ) ;

  template <>
  inline
  void
  setIdentity( uint32_t * AinvMmat, unsigned nr ) 
  { memcpy( AinvMmat, id32, sizeof(uint32_t) * nr ) ; }

  template <>
  inline
  void
  setIdentity( uint64_t * AinvMmat, unsigned nr )
  { memcpy( AinvMmat, id64, sizeof(uint64_t) * nr ) ; }

  template <>
  inline
  void
  setIdentity( uint128_t * AinvMmat, unsigned nr ) 
  { memcpy( AinvMmat, id128, sizeof(uint128_t) * nr ) ; }

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
               UNSIGNED * AinvMmat ) ;

  /*
  //                                        _    ____      _                       
  //    ___ ___  _ __ ___  _ __   __ _  ___| |_ / ___|___ | |_   _ _ __ ___  _ __  
  //   / __/ _ \| '_ ` _ \| '_ \ / _` |/ __| __| |   / _ \| | | | | '_ ` _ \| '_ \ 
  //  | (_| (_) | | | | | | |_) | (_| | (__| |_| |__| (_) | | |_| | | | | | | | | |
  //   \___\___/|_| |_| |_| .__/ \__,_|\___|\__|\____\___/|_|\__,_|_| |_| |_|_| |_|
  //                      |_|                                                      
  */
  template <typename UNSIGNED>
  void
  compactColumn( UNSIGNED * Lmat,     // colonna da compattare
                 unsigned   ncol,     // quante colonne contiene
                 UNSIGNED * Amat,     // colonna di arrivo
                 UNSIGNED * Amat1,    // colonna successiva
                 unsigned   startcol, // quante colonne sono gia occupate in Amat
                 unsigned   nRows ) ;

  template <typename UNSIGNED>
  void
  compactColumn( UNSIGNED * Lmat,     // colonna da compattare
                 unsigned   ncol,     // quante colonne contiene
                 UNSIGNED * Amat,     // colonna di arrivo (Amat e Lmat adiacenti)
                 unsigned   startcol, // quante colonne sono gia occupate in Amat
                 unsigned   nRows ) ;
}

#endif

// EOF GF2toolkit_InvertBlock.hh
