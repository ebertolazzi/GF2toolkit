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

#ifndef GF2TOOLKIT_FOUR_RUSSIAN_HH
#define GF2TOOLKIT_FOUR_RUSSIAN_HH

#include "GF2toolkit_common.hh"

/*
//   __  __ __  __ _  _   ____                _           _____     _     _
//  |  \/  |  \/  | || | |  _ \ _   _ ___ ___(_) __ _ _ _|_   _|_ _| |__ | | ___
//  | |\/| | |\/| | || |_| |_) | | | / __/ __| |/ _` | '_ \| |/ _` | '_ \| |/ _ \
//  | |  | | |  | |__   _|  _ <| |_| \__ \__ \ | (_| | | | | | (_| | |_) | |  __/
//  |_|  |_|_|  |_|  |_| |_| \_\\__,_|___/___/_|\__,_|_| |_|_|\__,_|_.__/|_|\___|
*/

#define GF2MM4R_DEFINE_METHODS(TYPE,SIZE)                                      \
  MM4R();                                                                      \
  MM4R( TYPE const * M );                                                      \
  MM4R( TYPE const * M, int32_t const * Perm );                                \
  void makeTable( TYPE const * M );                                            \
  void makeTable( TYPE const * M, int32_t const * Perm ) ;                     \
  void multRight( TYPE * A, unsigned nrows ) const ;                           \
  void multRightAss( TYPE const * A, unsigned nrows, TYPE * B ) const ;        \
  void multRightAdd( TYPE const * A, unsigned nrows, TYPE * B ) const

namespace GF2toolkit {

  template <unsigned nbits> void MM4RmakeTable( uint16_t  const M[], uint16_t  Table[] ) ;
  template <unsigned nbits> void MM4RmakeTable( uint32_t  const M[], uint32_t  Table[] ) ;
  template <unsigned nbits> void MM4RmakeTable( uint64_t  const M[], uint64_t  Table[] ) ;
  template <unsigned nbits> void MM4RmakeTable( uint128_t const M[], uint128_t Table[] ) ;

  template <typename UNSIGNED, unsigned nbits> 
  void
  MM4RmakeTables( UNSIGNED const M[], UNSIGNED * Table ) ;

  //////////////////////////////////////////////////////////////////////////////
  /*
  //  +---+   +---+ +---+    
  //  | C | = | A | | B |    
  //  |   |   |   | +---+
  //  |   |   |   |
  //  +---+   +---+
  */
  template <typename UNSIGNED, unsigned NBIT>
  void
  MM4RmultAssBlock( UNSIGNED const * A, unsigned nRowsA,
                    UNSIGNED const * B_Tables,
                    UNSIGNED       * C ) ;

  /*
  //  +---+    +---+ +---+    
  //  | C | ^= | A | | B |    
  //  |   |    |   | +---+
  //  |   |    |   |
  //  +---+    +---+
  */
  template <typename UNSIGNED, unsigned NBIT>
  void
  MM4RmultAddBlock( UNSIGNED const * A, unsigned nRowsA,
                    UNSIGNED const * B_Tables,
                    UNSIGNED       * C ) ;

  /*
  //  +---+    +---+ +---+    
  //  | A | <= | A | | B |    
  //  |   |    |   | +---+
  //  |   |    |   |
  //  +---+    +---+
  */
  template <typename UNSIGNED, unsigned NBIT>
  void
  MM4RmultRightBlock( UNSIGNED * A, unsigned nRowsA, UNSIGNED const * B_Tables ) ;

  //////////////////////////////////////////////////////////////////////////////

  template <typename UNSIGNED, unsigned NBIT>
  class MM4R {
  public:
    MM4R( ) ;
    MM4R( UNSIGNED const * M ) ;
    MM4R( UNSIGNED const * M, int32_t const * Perm ) ;
    void makeTable( UNSIGNED const * M ) ;
    void makeTable( UNSIGNED const * M, int32_t const * Perm ) ;
    void multRight( UNSIGNED * A, unsigned nrows ) const ;
    void multRightAss( UNSIGNED const * A, unsigned nrows, UNSIGNED * B ) const ;
    void multRightAdd( UNSIGNED const * A, unsigned nrows, UNSIGNED * B ) const ;
  } ;

  //////////////////////////////////////////////////////////////////////////////

  template <typename UNSIGNED>
  class MM4Rreduced {
    unsigned nblock ;
    UNSIGNED Tables[sizeof(UNSIGNED)*(1<<8)] ;
  public:
    MM4Rreduced( ) {}
    MM4Rreduced( UNSIGNED const * M, unsigned row_used ) { makeTable( M, row_used ) ; }

    void makeTable( UNSIGNED const * M, unsigned row_used ) ;

    void multRight( UNSIGNED * A, unsigned nrows ) const ;
    void multRightAss( UNSIGNED const * A, unsigned nrows, UNSIGNED * B ) const ;
    void multRightAdd( UNSIGNED const * A, unsigned nrows, UNSIGNED * B ) const ;

    void multRightShift( UNSIGNED * A, unsigned nrows, unsigned rshift ) const ;
    void multRightAssShift( UNSIGNED const * A, unsigned nrows, UNSIGNED * B, unsigned rshift ) const ;
    void multRightAddShift( UNSIGNED const * A, unsigned nrows, UNSIGNED * B, unsigned rshift ) const ;
  } ;

  /*
  //  ____    _     _ _   
  // |___ \  | |__ (_) |_ 
  //   __) | | '_ \| | __|
  //  / __/  | |_) | | |_ 
  // |_____| |_.__/|_|\__|
  */

  template <>
  class MM4R<uint32_t,2> {
    uint32_t Tables[16*4] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint32_t,2) ;
  } ;

  template <>
  class MM4R<uint64_t,2> {
    uint64_t Tables[32*4] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint64_t,2) ;
  } ;

  template <>
  class MM4R<uint128_t,2> {
    uint128_t Tables[64*4] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint128_t,2) ;
  } ;

  /*
  //   _____   _     _ _   
  //  |___ /  | |__ (_) |_ 
  //    |_ \  | '_ \| | __|
  //   ___) | | |_) | | |_ 
  //  |____/  |_.__/|_|\__|  
  */
  
  template <>
  class MM4R<uint32_t,3> {
    uint32_t Tables[8*(1<<3)+2*(1<<4)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint32_t,3) ;
  } ;

  template <>
  class MM4R<uint64_t,3> {
    uint64_t Tables[20*(1<<3)+(1<<4)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint64_t,3) ;
  } ;

  template <>
  class MM4R<uint128_t,3> {
    uint128_t Tables[40*(1<<3)+2*(1<<4)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint128_t,3) ;
  } ;
  
  /*
  //   _  _   _     _ _
  //  | || | | |__ (_) |_ 
  //  | || |_| '_ \| | __|
  //  |__   _| |_) | | |_ 
  //     |_| |_.__/|_|\__|
  */

  template <>
  class MM4R<uint32_t,4> {
    uint32_t Tables[8*(1<<4)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint32_t,4) ;
  } ;

  template <>
  class MM4R<uint64_t,4> {
    uint64_t Tables[16*(1<<4)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint64_t,4) ;
  } ;

  template <>
  class MM4R<uint128_t,4> {
    uint128_t Tables[32*(1<<4)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint128_t,4) ;
  } ;

  /*
  //  ____    _     _ _   
  // | ___|  | |__ (_) |_ 
  // |___ \  | '_ \| | __|
  //  ___) | | |_) | | |_ 
  // |____/  |_.__/|_|\__|
  */

  template <>
  class MM4R<uint32_t,5> {
    uint32_t Tables[5*(1<<5)+(1<<7)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint32_t,5) ;
  } ;

  template <>
  class MM4R<uint64_t,5> {
    uint64_t Tables[12*(1<<5)+(1<<4)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint64_t,5) ;
  } ;

  template <>
  class MM4R<uint128_t,5> {
    uint128_t Tables[25*(1<<5)+(1<<3)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint128_t,5) ;
  } ;

  /*
  //    __     _     _ _   
  //   / /_   | |__ (_) |_ 
  //  | '_ \  | '_ \| | __|
  //  | (_) | | |_) | | |_ 
  //   \___/  |_.__/|_|\__|
  */

  template <>
  class MM4R<uint32_t,6> {
    uint32_t Tables[3*(1<<6)+2*(1<<7)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint32_t,6) ;
  } ;

  template <>
  class MM4R<uint64_t,6> {
    uint64_t Tables[10*(1<<6)+(1<<4)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint64_t,6) ;
  } ;

  template <>
  class MM4R<uint128_t,6> {
    uint128_t Tables[19*(1<<6)+2*(1<<7)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint128_t,6) ;
  } ;

  /*
  //   _____   _     _ _   
  //  |___  | | |__ (_) |_ 
  //     / /  | '_ \| | __|
  //    / /   | |_) | | |_ 
  //   /_/    |_.__/|_|\__|
  */

  template <>
  class MM4R<uint32_t,7> {
    uint32_t Tables[4*(1<<7)+(1<<4)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint32_t,7) ;
  } ;

  template <>
  class MM4R<uint64_t,7> {
    uint64_t Tables[8*(1<<7)+(1<<8)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint64_t,7) ;
  } ;

  template <>
  class MM4R<uint128_t,7> {
    uint128_t Tables[16*(1<<7)+2*(1<<8)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint128_t,7) ;
  } ;

  /*
  //    ___  _     _ _   
  //   ( _ )| |__ (_) |_ 
  //   / _ \| '_ \| | __|
  //  | (_) | |_) | | |_ 
  //   \___/|_.__/|_|\__|
  */

  template <>
  class MM4R<uint32_t,8> {
    uint32_t Tables[4*(1<<8)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint32_t,8) ;
  } ;

  template <>
  class MM4R<uint64_t,8> {
    uint64_t Tables[8*(1<<8)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint64_t,8) ;
  } ;

  template <>
  class MM4R<uint128_t,8> {
    uint128_t Tables[16*(1<<8)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint128_t,8) ;
  } ;

  /*
  //    ___    _     _ _   
  //   / _ \  | |__ (_) |_ 
  //  | (_) | | '_ \| | __|
  //   \__, | | |_) | | |_ 
  //     /_/  |_.__/|_|\__|
  */

  template <>
  class MM4R<uint32_t,9> {
    uint32_t Tables[3*(1<<9)+(1<<5)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint32_t,9) ;
  } ;

  template <>
  class MM4R<uint64_t,9> {
    uint64_t Tables[6*(1<<9)+(1<<10)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint64_t,9) ;
  } ;

  template <>
  class MM4R<uint128_t,9> {
    uint128_t Tables[12*(1<<9)+2*(1<<10)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint128_t,9) ;
  } ;

  /*
  //   _  ___    _     _ _   
  //  / |/ _ \  | |__ (_) |_ 
  //  | | | | | | '_ \| | __|
  //  | | |_| | | |_) | | |_ 
  //  |_|\___/  |_.__/|_|\__|
  */
  template <>
  class MM4R<uint32_t,10> {
    uint32_t Tables[2*(1<<10)+(1<<12)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint32_t,10) ;
  } ;

  template <>
  class MM4R<uint64_t,10> {
    uint64_t Tables[6*(1<<10)+(1<<4)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint64_t,10) ;
  } ;

  template <>
  class MM4R<uint128_t,10> {
    uint128_t Tables[12*(1<<10)+(1<<8)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint128_t,10) ;
  } ;

  /*
  //   _ _   _     _ _   
  //  / / | | |__ (_) |_ 
  //  | | | | '_ \| | __|
  //  | | | | |_) | | |_ 
  //  |_|_| |_.__/|_|\__|
  */

  template <>
  class MM4R<uint32_t,11> {
    uint32_t Tables[2*(1<<11)+(1<<10)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint32_t,11) ;
  } ;

  template <>
  class MM4R<uint64_t,11> {
    uint64_t Tables[4*(1<<11)+2*(1<<10)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint64_t,11) ;
  } ;

  template <>
  class MM4R<uint128_t,11> {
    uint128_t Tables[8*(1<<11)+4*(1<<10)] ;
  public:
    GF2MM4R_DEFINE_METHODS(uint128_t,11) ;
  } ;

}

#endif

// EOF GF2toolkit_4Russian.hh
