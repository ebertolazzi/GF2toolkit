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

#define GF2MM4R_BODY_OF_METHODS(TYPE,SIZE)                                        \
  MM4R<TYPE,SIZE>::MM4R() {}                                                      \
  MM4R<TYPE,SIZE>::MM4R( TYPE const * M ) { makeTable(M) ; }                      \
  MM4R<TYPE,SIZE>::MM4R( TYPE const * M, int32_t const * Perm )                   \
  { makeTable(M,Perm) ; }                                                         \
                                                                                  \
  void                                                                            \
  MM4R<TYPE,SIZE>::makeTable( TYPE const * M )                                    \
  { MM4RmakeTables<TYPE,SIZE>( M, Tables ) ; }                                    \
                                                                                  \
  void                                                                            \
  MM4R<TYPE,SIZE>::multRight( TYPE * A, unsigned nrows ) const                    \
  { MM4RmultRightBlock<TYPE,SIZE>( A, nrows, Tables ) ; }                         \
                                                                                  \
  void                                                                            \
  MM4R<TYPE,SIZE>::multRightAss( TYPE const * A, unsigned nrows, TYPE * B ) const \
  { MM4RmultAssBlock<TYPE,SIZE>( A, nrows, Tables, B ) ; }                        \
                                                                                  \
  void                                                                            \
  MM4R<TYPE,SIZE>::multRightAdd( TYPE const * A, unsigned nrows, TYPE * B ) const \
  { MM4RmultAddBlock<TYPE,SIZE>( A, nrows, Tables, B ) ; }

/*
//  __  __ __  __ _  _   ____                _           _____     _     _
// |  \/  |  \/  | || | |  _ \ _   _ ___ ___(_) __ _ _ _|_   _|_ _| |__ | | ___
// | |\/| | |\/| | || |_| |_) | | | / __/ __| |/ _` | '_ \| |/ _` | '_ \| |/ _ \
// | |  | | |  | |__   _|  _ <| |_| \__ \__ \ | (_| | | | | | (_| | |_) | |  __/
// |_|  |_|_|  |_|  |_| |_| \_\\__,_|___/___/_|\__,_|_| |_|_|\__,_|_.__/|_|\___|
*/

namespace GF2toolkit {

  #define GF2MM4R_BODY_OF_MAKETABLE(NBIT)                                      \
  void                                                                         \
  MM4R<uint32_t,NBIT>::makeTable( uint32_t const * M, int32_t const * Perm ) { \
    uint32_t MM[32] ;                                                          \
    MM4RexpandRows<uint32_t>( M, Perm, MM ) ;                                  \
    MM4RmakeTables<uint32_t,NBIT>( MM, Tables ) ;                              \
  }                                                                            \
                                                                               \
  void                                                                         \
  MM4R<uint64_t,NBIT>::makeTable( uint64_t const * M, int32_t const * Perm ) { \
    uint64_t MM[64] ;                                                          \
    MM4RexpandRows<uint64_t>( M, Perm, MM ) ;                                  \
    MM4RmakeTables<uint64_t,NBIT>( MM, Tables ) ;                              \
  }                                                                            \
                                                                               \
  void                                                                         \
  MM4R<uint128_t,NBIT>::makeTable( uint128_t const * M, int32_t const * Perm ){\
    uint128_t MM[128] ;                                                        \
    MM4RexpandRows<uint128_t>( M, Perm, MM ) ;                                 \
    MM4RmakeTables<uint128_t,NBIT>( MM, Tables ) ;                             \
  }
  
  /*                                   _ ____                   
  //    _____  ___ __   __ _ _ __   __| |  _ \ _____      _____ 
  //   / _ \ \/ / '_ \ / _` | '_ \ / _` | |_) / _ \ \ /\ / / __|
  //  |  __/>  <| |_) | (_| | | | | (_| |  _ < (_) \ V  V /\__ \
  //   \___/_/\_\ .__/ \__,_|_| |_|\__,_|_| \_\___/ \_/\_/ |___/
  //            |_|                                             
  */
  template <typename UNSIGNED> 
  void
  MM4RexpandRows( UNSIGNED const M[], int32_t const Perm[], UNSIGNED MM[] ) ;

  template <>
  inline 
  void
  MM4RexpandRows( uint32_t const M[], int32_t const Perm[], uint32_t MM[] ) {
    unsigned k = 0 ;
    for ( unsigned i = 0 ; i < 32 ; ++i )
      if ( Perm[i] < 0 ) MM[i] = 0 ;
      else               MM[i] = M[k++] ;
  }

  template <>
  inline 
  void
  MM4RexpandRows( uint64_t const M[], int32_t const Perm[], uint64_t MM[] ) {
    unsigned k = 0 ;
    for ( unsigned i = 0 ; i < 64 ; ++i )
      if ( Perm[i] < 0 ) MM[i] = 0 ;
      else               MM[i] = M[k++] ;
  }

  template <>
  inline 
  void
  MM4RexpandRows( uint128_t const M[], int32_t const Perm[], uint128_t MM[] ) {
    unsigned k = 0 ;
    for ( unsigned i = 0 ; i < 128 ; ++i )
      if ( Perm[i] < 0 ) MM[i].i128 = _mm_setzero_si128() ;
      else               MM[i].i128 = M[k++].i128 ;
  }

  /*
  //  ____    _     _ _   
  // |___ \  | |__ (_) |_ 
  //   __) | | '_ \| | __|
  //  / __/  | |_) | | |_ 
  // |_____| |_.__/|_|\__|
  */

  GF2MM4R_BODY_OF_METHODS(uint32_t,2)
  GF2MM4R_BODY_OF_METHODS(uint64_t,2)
  GF2MM4R_BODY_OF_METHODS(uint128_t,2)
  GF2MM4R_BODY_OF_MAKETABLE(2)

  /*
  //   _____   _     _ _   
  //  |___ /  | |__ (_) |_ 
  //    |_ \  | '_ \| | __|
  //   ___) | | |_) | | |_ 
  //  |____/  |_.__/|_|\__|  
  */

  GF2MM4R_BODY_OF_METHODS(uint32_t,3)
  GF2MM4R_BODY_OF_METHODS(uint64_t,3)
  GF2MM4R_BODY_OF_METHODS(uint128_t,3)
  GF2MM4R_BODY_OF_MAKETABLE(3)

  /*
  //   _  _   _     _ _
  //  | || | | |__ (_) |_ 
  //  | || |_| '_ \| | __|
  //  |__   _| |_) | | |_ 
  //     |_| |_.__/|_|\__|
  */

  GF2MM4R_BODY_OF_METHODS(uint32_t,4)
  GF2MM4R_BODY_OF_METHODS(uint64_t,4)
  GF2MM4R_BODY_OF_METHODS(uint128_t,4)
  GF2MM4R_BODY_OF_MAKETABLE(4)

  /*
  //  ____    _     _ _   
  // | ___|  | |__ (_) |_ 
  // |___ \  | '_ \| | __|
  //  ___) | | |_) | | |_ 
  // |____/  |_.__/|_|\__|
  */

  GF2MM4R_BODY_OF_METHODS(uint32_t,5)
  GF2MM4R_BODY_OF_METHODS(uint64_t,5)
  GF2MM4R_BODY_OF_METHODS(uint128_t,5)
  GF2MM4R_BODY_OF_MAKETABLE(5)

  /*
  //    __     _     _ _   
  //   / /_   | |__ (_) |_ 
  //  | '_ \  | '_ \| | __|
  //  | (_) | | |_) | | |_ 
  //   \___/  |_.__/|_|\__|
  */

  GF2MM4R_BODY_OF_METHODS(uint32_t,6)
  GF2MM4R_BODY_OF_METHODS(uint64_t,6)
  GF2MM4R_BODY_OF_METHODS(uint128_t,6)
  GF2MM4R_BODY_OF_MAKETABLE(6)

  /*
  //   _____   _     _ _   
  //  |___  | | |__ (_) |_ 
  //     / /  | '_ \| | __|
  //    / /   | |_) | | |_ 
  //   /_/    |_.__/|_|\__|
  */

  GF2MM4R_BODY_OF_METHODS(uint32_t,7)
  GF2MM4R_BODY_OF_METHODS(uint64_t,7)
  GF2MM4R_BODY_OF_METHODS(uint128_t,7)
  GF2MM4R_BODY_OF_MAKETABLE(7)

  /*
  //    ___  _     _ _   
  //   ( _ )| |__ (_) |_ 
  //   / _ \| '_ \| | __|
  //  | (_) | |_) | | |_ 
  //   \___/|_.__/|_|\__|
  */

  GF2MM4R_BODY_OF_METHODS(uint32_t,8)
  GF2MM4R_BODY_OF_METHODS(uint64_t,8)
  GF2MM4R_BODY_OF_METHODS(uint128_t,8)
  GF2MM4R_BODY_OF_MAKETABLE(8)

  /*
  //    ___    _     _ _   
  //   / _ \  | |__ (_) |_ 
  //  | (_) | | '_ \| | __|
  //   \__, | | |_) | | |_ 
  //     /_/  |_.__/|_|\__|
  */

  GF2MM4R_BODY_OF_METHODS(uint32_t,9)
  GF2MM4R_BODY_OF_METHODS(uint64_t,9)
  GF2MM4R_BODY_OF_METHODS(uint128_t,9)
  GF2MM4R_BODY_OF_MAKETABLE(9)

  /*
  //   _  ___    _     _ _   
  //  / |/ _ \  | |__ (_) |_ 
  //  | | | | | | '_ \| | __|
  //  | | |_| | | |_) | | |_ 
  //  |_|\___/  |_.__/|_|\__|
  */                   

  GF2MM4R_BODY_OF_METHODS(uint32_t,10)
  GF2MM4R_BODY_OF_METHODS(uint64_t,10)
  GF2MM4R_BODY_OF_METHODS(uint128_t,10)
  GF2MM4R_BODY_OF_MAKETABLE(10)

  /*
  //   _ _   _     _ _   
  //  / / | | |__ (_) |_ 
  //  | | | | '_ \| | __|
  //  | | | | |_) | | |_ 
  //  |_|_| |_.__/|_|\__|
  */

  GF2MM4R_BODY_OF_METHODS(uint32_t,11)
  GF2MM4R_BODY_OF_METHODS(uint64_t,11)
  GF2MM4R_BODY_OF_METHODS(uint128_t,11)
  GF2MM4R_BODY_OF_MAKETABLE(11)

}

// EOF GF2toolkit_4RussianMultClasses.hh
