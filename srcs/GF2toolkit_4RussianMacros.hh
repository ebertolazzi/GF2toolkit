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

#ifndef GF2_FOUR_RUSSIAN_MACROS_HH
#define GF2_FOUR_RUSSIAN_MACROS_HH

#define BITS2(A,N)  (((A)>>((N)*2))&0x03)
#define BITS3(A,N)  (((A)>>((N)*3))&0x07)
#define BITS4(A,N)  (((A)>>((N)*4))&0x0F)
#define BITS5(A,N)  (((A)>>((N)*5))&0x1F)
#define BITS6(A,N)  (((A)>>((N)*6))&0x3F)
#define BITS7(A,N)  (((A)>>((N)*7))&0x7F)
#define BITS8(A,N)  (((A)>>((N)*8))&0xFF)
#define BITS9(A,N)  (((A)>>((N)*9))&0x1FF)
#define BITS10(A,N) (((A)>>((N)*10))&0x3FF)
#define BITS11(A,N) (((A)>>((N)*11))&0x7FF)
#define BITS12(A,N) (((A)>>((N)*12))&0xFFF)

#define MM4R_LOOP_UNROLLED_BY_16                                               \
  switch ( nRowsA & 0x0F ) {                                                   \
    case 15: MM4R_SINGLE_OP ;                                                  \
    case 14: MM4R_SINGLE_OP ;                                                  \
    case 13: MM4R_SINGLE_OP ;                                                  \
    case 12: MM4R_SINGLE_OP ;                                                  \
    case 11: MM4R_SINGLE_OP ;                                                  \
    case 10: MM4R_SINGLE_OP ;                                                  \
    case 9:  MM4R_SINGLE_OP ;                                                  \
    case 8:  MM4R_SINGLE_OP ;                                                  \
    case 7:  MM4R_SINGLE_OP ;                                                  \
    case 6:  MM4R_SINGLE_OP ;                                                  \
    case 5:  MM4R_SINGLE_OP ;                                                  \
    case 4:  MM4R_SINGLE_OP ;                                                  \
    case 3:  MM4R_SINGLE_OP ;                                                  \
    case 2:  MM4R_SINGLE_OP ;                                                  \
    case 1:  MM4R_SINGLE_OP ;                                                  \
  }                                                                            \
  unsigned idx = nRowsA>>4 ;                                                   \
  while ( idx-- > 0 ) {                                                        \
    MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ;        \
    MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ;        \
    MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ;        \
    MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ;        \
  }

#define MM4R_LOOP_UNROLLED_BY_8                                                \
  switch ( nRowsA & 0x07 ) {                                                   \
    case 7:  MM4R_SINGLE_OP ;                                                  \
    case 6:  MM4R_SINGLE_OP ;                                                  \
    case 5:  MM4R_SINGLE_OP ;                                                  \
    case 4:  MM4R_SINGLE_OP ;                                                  \
    case 3:  MM4R_SINGLE_OP ;                                                  \
    case 2:  MM4R_SINGLE_OP ;                                                  \
    case 1:  MM4R_SINGLE_OP ;                                                  \
  }                                                                            \
  unsigned idx = nRowsA>>3 ;                                                   \
  while ( idx-- > 0 ) {                                                        \
    MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ;        \
    MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ;        \
  }

#define MM4R_LOOP_UNROLLED_BY_4                                                \
  switch ( nRowsA & 0x03 ) {                                                      \
    case 3:  MM4R_SINGLE_OP ;                                                  \
    case 2:  MM4R_SINGLE_OP ;                                                  \
    case 1:  MM4R_SINGLE_OP ;                                                  \
  }                                                                            \
  unsigned idx = nRowsA>>2 ;                                                   \
  while ( idx-- > 0 ) {                                                        \
    MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ; MM4R_SINGLE_OP ;        \
  }

#define MM4R_LOOP_UNROLLED_BY_2                                                \
  if ( nRowsA & 0x01 ) MM4R_SINGLE_OP ;                                        \
  unsigned idx = nRowsA>>1 ;                                                   \
  while ( idx-- > 0 ) {                                                        \
    MM4R_SINGLE_OP ; MM4R_SINGLE_OP ;                                          \
  }


#endif

// EOF GF2_4RussianMacros.hh
