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
 *    GF2 is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    GF2 is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADO Toolkit; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

/* This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 3 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.  You should have received
   a copy of the GNU General Public License along with this program;
   if not, write to the Free Foundation, Inc., 59 Temple Place, Suite
   330, Boston, MA 02111-1307 USA
   
   This implementation copyright (C) 2010 Enrico Bertolazzi.
   I adapt the random number generator of Makoto Matsumoto and Takuji Nishimura,
   implemented in the GSL - GNU Scientific Library (http://www.gnu.org/software/gsl/).

   The original code of Matsumoto and Nishimura included the comment:
   "When you use this, send an email to: matumoto@math.keio.ac.jp with
   an appropriate reference to your work".

   Makoto Matsumoto has a web page with more information about the
   generator, http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html 

   The paper below has details of the algorithm.

   From: Makoto Matsumoto and Takuji Nishimura, "Mersenne Twister: A
   623-dimensionally equidistributerd uniform pseudorandom number
   generator". ACM Transactions on Modeling and Computer Simulation,
   Vol. 8, No. 1 (Jan. 1998), Pages 3-30

   You can obtain the paper directly from Makoto Matsumoto's web page.

   The period of this generator is 2^{19937} - 1.

*/


/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/


#include "Random.hh"

#include <cmath>
#include <stdlib.h>
#include <string>

namespace RNG {

  using namespace std ;

  /*   ____              
  //  |  _ \ _ __   __ _ 
  //  | |_) | '_ \ / _` |
  //  |  _ <| | | | (_| |
  //  |_| \_\_| |_|\__, |
  //               |___/ 
  */

  /* These real versions are due to Isaku Wada, 2002/01/09 added */

  /* generates a random number on [0,1) with 53-bit resolution */
  /*
  double
  Rng::rand53() {
    uint32_t a = this -> rand() >> 5 ;
    uint32_t b = this -> rand() >> 6 ; 
    return (a*67108864.0+b)*(1.0/9007199254740992.0); 
  }
  */

  /*   __  __                                    _____          _     _            
  //  |  \/  | ___ _ __ ___  ___ _ __  _ __   __|_   _|_      _(_)___| |_ ___ _ __ 
  //  | |\/| |/ _ \ '__/ __|/ _ \ '_ \| '_ \ / _ \| | \ \ /\ / / / __| __/ _ \ '__|
  //  | |  | |  __/ |  \__ \  __/ | | | | | |  __/| |  \ V  V /| \__ \ ||  __/ |   
  //  |_|  |_|\___|_|  |___/\___|_| |_|_| |_|\___||_|   \_/\_/ |_|___/\__\___|_|   
  */                                                                             

  // Period parameters
  #define N          624
  #define M          397
  #define MATRIX_A   uint32_t(0x9908b0df) // constant vector a
  #define UPPER_MASK uint32_t(0x80000000) // most significant w-r bits
  #define LOWER_MASK uint32_t(0x7fffffff) // least significant r bits

  MersenneTwister::MersenneTwister( uint32_t s )
  : Rng(two_pow_32_m_1,"MersenneTwister")
  , mti(N+1) // mti==N+1 means mt[N] is not initialized
  { seed(s) ; }
  
  MersenneTwister::~MersenneTwister()
  {}

  void
  MersenneTwister::seed( uint32_t s ) {
    mt[0]= s ;
    for ( mti = 1 ; mti < N ; ++mti ) {
      mt[mti] = (uint32_t(1812433253) * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti) ;
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      mt[mti] &= uint32_t(0xffffffff) ; /* for >32 bit machines */
    }
  }
  
  /* initialize by an array with array-length */
  /* init_key is the array for initializing keys */
  /* key_length is its length */
  /* slight change for C++, 2004/2/26 */
  void
  MersenneTwister::seed( uint32_t init_key[], int key_length ) {
    seed(uint32_t(19650218)) ;
    int i = 1 ;
    int j = 0 ;
    int k = key_length > N ? N : key_length ;
    for ( ; k ; --k ) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * uint32_t(1664525))) + init_key[j] + j; /* non linear */
      mt[i] &= uint32_t(0xffffffff); /* for WORDSIZE > 32 machines */
      ++i ; ++j ;
      if ( i >= N ) { mt[0] = mt[N-1] ; i=1 ; }
      if ( j >= key_length ) j=0 ;
    }
    for ( k = N-1 ; k ; --k ) {
      mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * uint32_t(1566083941))) - i ; /* non linear */
      mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      if ( ++i >= N ) { mt[0] = mt[N-1] ; i=1 ; }
    }
    mt[0] = UINT32(0x80000000); /* MSB is 1; assuring non-zero initial array */ 
  }

  /* generates a random number on [0,0xffffffff]-interval */
  uint32_t
  MersenneTwister::rand() {
    uint32_t y ;
    static uint32_t const mag01[2] = { UINT32(0x0), MATRIX_A } ;
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if ( mti >= N ) { /* generate N words at one time */
      int kk;
      if ( mti == N+1 ) seed(UINT32(5489)) ; // if seed() has not been called, a default initial seed is used
      for ( kk = 0 ; kk < N-M ; ++kk ) {
        y      = (mt[kk]&UPPER_MASK) | (mt[kk+1]&LOWER_MASK) ;
        mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y&0x1] ;
      }
      for ( ; kk < N-1 ; ++kk ) {
        y = (mt[kk]&UPPER_MASK) | (mt[kk+1]&LOWER_MASK) ;
        mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & UINT32(0x1)] ;
      }
      y       = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK) ;
      mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];
      mti     = 0 ;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11) ;
    y ^= (y << 7)  & UINT32(0x9d2c5680) ;
    y ^= (y << 15) & UINT32(0xefc60000) ;
    y ^= (y >> 18) ;
    return y ;
  }
  
  #undef N
  #undef M
  #undef MATRIX_A
  #undef UPPER_MASK
  #undef LOWER_MASK

  /*    ____ __  __ ____   ____ 
  //   / ___|  \/  |  _ \ / ___|
  //  | |   | |\/| | |_) | |  _ 
  //  | |___| |  | |  _ <| |_| |
  //   \____|_|  |_|_| \_\\____|
  */     
  
  /*
   
   This is a combined multiple recursive generator. The sequence is,

   z_n = (x_n - y_n) mod m1

   where the two underlying generators x and y are,

   x_n = (a_{1} x_{n-1} + a_{2} x_{n-2} + a_{3} x_{n-3}) mod m1
   y_n = (b_{1} y_{n-1} + b_{2} y_{n-2} + b_{3} y_{n-3}) mod m2

   with coefficients a11 ... a23,

   a_{1} = 0,     a_{2} = 63308, a_{3} = -183326
   b_{1} = 86098, b_{2} = 0,     b_{3} = -539608

   and moduli m1, m2,

   m1 = 2^31 - 1 = 2147483647
   m2 = 2^31 - 2000169 = 2145483479

   We initialize the generator with 

   x_1 = s_1 MOD m1, x_2 = s_2 MOD m1, x_3 = s_3 MOD m1
   y_1 = s_4 MOD m2, y_2 = s_5 MOD m2, y_3 = s_6 MOD m2

   where s_n = (69069 * s_{n-1}) mod 2^32 and s_0 = s is the
   user-supplied seed.

   NOTE: According to the paper the initial values for x_n must lie in
   the range 0 <= x_n <= (m1 - 1) and the initial values for y_n must
   lie in the range 0 <= y_n <= (m2 - 1), with at least one non-zero
   value -- our seeding procedure satisfies these constraints.

   We then use 7 iterations of the generator to "warm up" the internal
   state.

   The theoretical value of z_{10008} is 719452880. The subscript 10008
   means (1) seed the generator with s=1, (2) do the seven warm-up
   iterations that are part of the seeding process, (3) then do 10000
   actual iterations.

   The period of this generator is about 2^205.

   From: P. L'Ecuyer, "Combined Multiple Recursive Random Number
   Generators," Operations Research, 44, 5 (1996), 816--822.

   This is available on the net from L'Ecuyer's home page,

   http://www.iro.umontreal.ca/~lecuyer/myftp/papers/combmrg.ps
   ftp://ftp.iro.umontreal.ca/pub/simulation/lecuyer/papers/combmrg.ps
  */

  void
  CMRG::seed( uint32_t s_in ) {

    // An entirely adhoc way of seeding! This does **not** come from L'Ecuyer et al
    if ( s_in == 0 ) s_in = 1 ; // default seed is 1
    
    int64_t s = s_in ;

    #define LCG(n) ((69069 * n) & UINT32(0xffffffff))
    s = LCG(s) ; x1 = s % m1 ;
    s = LCG(s) ; x2 = s % m1 ;
    s = LCG(s) ; x3 = s % m1 ;

    s = LCG(s) ; y1 = s % m2 ;
    s = LCG(s) ; y2 = s % m2 ;
    s = LCG(s) ; y3 = s % m2 ;
    #undef LCG

    // "warm it up"
    rand(); rand(); rand(); rand(); rand(); rand(); rand();
  }

  uint32_t
  CMRG::rand() {

    // Component 1
    int64_t xNew = (1403580 * x2 -  810728 * x3) % m1 ;
    int64_t yNew = (527612  * y1 - 1370589 * y3) % m2 ;

    x3 = x2 ; x2 = x1 ; x1 = xNew ;
    y3 = y2 ; y2 = y1 ; y1 = yNew ;
 
    return static_cast<uint32_t>((xNew - yNew) % m1) ;

  }
  
  /*   _____ _     _                           ____       
  //  |  ___(_)___| |__  _ __ ___   __ _ _ __ |___ \__  __
  //  | |_  | / __| '_ \| '_ ` _ \ / _` | '_ \  __) \ \/ /
  //  |  _| | \__ \ | | | | | | | | (_| | | | |/ __/ >  < 
  //  |_|   |_|___/_| |_|_| |_| |_|\__,_|_| |_|_____/_/\_\
  */
  /* Fishman */
  #define AAA_F UINT32(48271)
  #define MMM_F UINT32(0x7fffffff)      /* 2 ^ 31 - 1 */
  #define QQQ_F UINT32(44488)
  #define RRR_F UINT32(3399)

  /* L'Ecuyer */
  #define AAA_L UINT32(40692)
  #define MMM_L UINT32(0x7fffff07)      /* 2 ^ 31 - 249 */
  #define QQQ_L UINT32(52774)
  #define RRR_L UINT32(3791)

  void
  Fishman2x::seed( uint32_t s ) {
    if ((s % MMM_F) == 0 || (s % MMM_L) == 0) s = 1 ; // default seed is 1
    sx = s % MMM_F;
    sy = s % MMM_L;
    sz = sx > sy ? sx - sy : MMM_F + sx - sy ;
  }

  uint32_t
  Fishman2x::rand() {
    int32_t r = RRR_F * (sx / QQQ_F);
    int32_t y = AAA_F * (sx % QQQ_F) - r;
    if (y < 0) y += MMM_F;
    sx = y;

    r = RRR_L * (sy / QQQ_L);
    y = AAA_L * (sy % QQQ_L) - r;
    if (y < 0) y += MMM_L;
    sy = y;

    sz = sx > sy ? sx - sy : MMM_F + sx - sy ;
    return sz ;
  }
  
  #undef AAA_F
  #undef MMM_F
  #undef QQQ_F
  #undef RRR_F
  #undef AAA_L
  #undef MMM_L
  #undef QQQ_L
  #undef RRR_L

  /*   _  __            _   _     
  //  | |/ /_ __  _   _| |_| |__  
  //  | ' /| '_ \| | | | __| '_ \ 
  //  | . \| | | | |_| | |_| | | |
  //  |_|\_\_| |_|\__,_|\__|_| |_|
  */
  #define BUFLEN 2009        // [Brian]: length of the buffer aa[]
  #define KK     100         // the long lag
  #define LL     37          // the short lag
  #define MM     (1L << 30)  // the modulus
  #define TT     70          // guaranteed separation between streams

  #define evenize(x)     ((x) & (MM - 2))     /* make x even */
  #define is_odd(x)      ((x) & 1)     /* the units bit of x */
  #define mod_diff(x, y) (((x) - (y)) & (MM - 1)) /* (x - y) mod MM */

  void
  Knuth::seed( uint32_t s ) {
    long x[KK + KK - 1]; // the preparation buffer

    register int j ;
    register int t ;
    register long ss = evenize(s + 2) ;

    for ( j = 0 ; j < KK ; ++j ) {
      x[j] = ss ; // bootstrap the buffer
      ss <<= 1 ;
      if (ss >= MM) ss -= MM - 2 ; // cyclic shift 29 bits
    }
    for (      ; j < KK + KK - 1; ++j )
      x[j] = 0 ;
    x[1]++ ; // make x[1] (and only x[1]) odd
    ss = s & (MM - 1) ;
    t  = TT - 1 ;
    while ( t ) {
      for (j = KK - 1      ; j > 0       ; --j    ) x[j + j] = x[j] ; // square
      for (j = KK + KK - 2 ; j > KK - LL ; j -= 2 ) x[KK + KK - 1 - j] = evenize(x[j]);
      for (j = KK + KK - 2; j >= KK; j--)
        if (is_odd (x[j])) {
          x[j - (KK - LL)] = mod_diff(x[j - (KK - LL)], x[j]);
          x[j - KK]        = mod_diff(x[j - KK], x[j]);
        }
      if ( is_odd(ss) ) { // multiply by "z"
        for ( j = KK ; j > 0 ; --j ) x[j] = x[j - 1] ;
        x[0] = x[KK] ; // shift the buffer cyclically
        if ( is_odd(x[KK]) ) x[LL] = mod_diff (x[LL], x[KK]);
      }
      if (ss) ss >>= 1 ;
      else    --t ;
    }
    idx = 0 ;
    for ( j = 0 ; j < LL ; ++j ) ran_x[j + KK - LL] = x[j] ;
    for (       ; j < KK ; ++j ) ran_x[j - LL] = x[j] ;
  }
  
  uint32_t
  Knuth::rand() {
    if (idx == 0) {
      /* fill buffer with new random numbers */
      unsigned i, j ;
      for ( j = 0 ; j < KK     ; ++j      ) aa[j]    = ran_x[j];
      for (       ; j < BUFLEN ; ++j      ) aa[j]    = mod_diff(aa[j - KK], aa[j - LL]);
      for ( i = 0 ; i < LL     ; ++i, ++j ) ran_x[i] = mod_diff(aa[j - KK], aa[j - LL]);
      for (       ; i < KK     ; ++i, ++j ) ran_x[i] = mod_diff(aa[j - KK], ran_x[i - LL]);
    }
    idx = (idx + 1) % BUFLEN;
    return aa[idx];
  }

  /*    ____ _____ ____  ____  _  _   
  //   / ___|  ___/ ___||  _ \| || |  
  //  | |  _| |_  \___ \| |_) | || |_ 
  //  | |_| |  _|  ___) |  _ <|__   _|
  //   \____|_|   |____/|_| \_\  |_|  
  */
  /* Magic numbers */
  #define A 471
  #define B 1586
  #define C 6988
  #define D 9689
  #define M 16383 /* = 2^14-1 */
  /* #define M 0x0003fff */
  #define LCG(n) ((69069 * n) & UINT32(0xffffffff))

  void
  GFSR4::seed( uint32_t s ) {
    // Masks for turning on the diagonal bit and turning off the leftmost bits
    uint32_t msb  = UINT32(0x80000000) ;
    uint32_t mask = UINT32(0xffffffff) ;

    if ( s == 0 ) s = 4357 ; // the default seed is 4357

    /* We use the congruence s_{n+1} = (69069*s_n) mod 2^32 to
       initialize the state. This works because ANSI-C unsigned long
       integer arithmetic is automatically modulo 2^32 (or a higher
       power of two), so we can safely ignore overflow. */

    // Brian Gough suggests this to avoid low-order bit correlations
    for ( int i = 0 ; i <= M ; ++i ) {
      uint32_t t   = 0 ;
      uint32_t bit = msb ;
      for ( int j = 0 ; j < 32 ; ++j ) {
        s = LCG(s) ;
        if ( s & msb ) t |= bit ;
        bit >>= 1 ;
      }
      ra[i] = t ;
    }

    /* Perform the "orthogonalization" of the matrix */
    /* Based on the orthogonalization used in r250, as suggested initially
     * by Kirkpatrick and Stoll, and pointed out to me by Brian Gough
     */

    /* BJG: note that this orthogonalisation doesn't have any effect
       here because the the initial 6695 elements do not participate in
       the calculation.  For practical purposes this orthogonalisation
       is somewhat irrelevant, because the probability of the original
       sequence being degenerate should be exponentially small. */

    for (int i = 0 ; i < 32 ; ++i ) {
      int k=7+i*3 ;
      ra[k] &= mask ; // Turn off bits left of the diagonal
      ra[k] |= msb  ; // Turn on the diagonal bit
      mask >>= 1 ;
      msb  >>= 1 ;
    }

    nd = 32 ;
  }
  
  uint32_t
  GFSR4::rand() {
    nd = (nd+1) & M ;
    ra[nd] = ra[(nd+(M+1-A))&M]^
             ra[(nd+(M+1-B))&M]^
             ra[(nd+(M+1-C))&M]^
             ra[(nd+(M+1-D))&M];
    return ra[nd] ;
  }
  
  /*   __  __ ____   ____ 
  //  |  \/  |  _ \ / ___|
  //  | |\/| | |_) | |  _ 
  //  | |  | |  _ <| |_| |
  //  |_|  |_|_| \_\\____|
  */
  #define LCG(n) ((69069 * n) & UINT32(0xffffffff))
  void
  MRG::seed( uint32_t s ) {
    int32_t const m  = 2147483647 ;

    // An entirely adhoc way of seeding! This does **not** come from L'Ecuyer et al
    if (s == 0) s = 1 ; // default seed is 1

    s = LCG(s) ; x1 = s % m;
    s = LCG(s) ; x2 = s % m;
    s = LCG(s) ; x3 = s % m;
    s = LCG(s) ; x4 = s % m;
    s = LCG(s) ; x5 = s % m;

    // "warm it up" with at least 5 calls to go through all the x values
    rand() ; rand() ; rand() ; rand() ; rand() ; rand() ;
  }
  
  uint32_t
  MRG::rand() {
    int32_t const m  = 2147483647 ;
    int32_t const a1 = 107374182 ;
    int32_t const q1 = 20 ;
    int32_t const r1 = 7 ;
    int32_t const a5 = 104480 ;
    int32_t const q5 = 20554 ;
    int32_t const r5 = 1727 ;

    int32_t h5 = x5 / q5;
    int32_t p5 = a5 * (x5 - h5 * q5) - h5 * r5;
    if ( p5 > 0 ) p5 -= m ;

    int32_t h1 = x1 / q1;
    int32_t p1 = a1 * (x1 - h1 * q1) - h1 * r1;
    if ( p1 < 0 ) p1 += m;

    x5 = x4 ;
    x4 = x3 ;
    x3 = x2 ;
    x2 = x1 ;

    x1 = p1 + p5;

    if ( x1 < 0 ) x1 += m ;

    return x1 ;
  }
  
  /*   _                         _               
  //  | |   _   _  ___  ___  ___| |__   ___ _ __ 
  //  | |  | | | |/ _ \/ __|/ __| '_ \ / _ \ '__|
  //  | |__| |_| |  __/\__ \ (__| | | |  __/ |   
  //  |_____\__,_|\___||___/\___|_| |_|\___|_|
  */
  void     
  Luescher::seed( uint32_t s ) {
    uint32_t const mask_lo = UINT32(0x00ffffff) ; // 2^24 - 1
    uint32_t const mask_hi = ~mask_lo ;

    if ( s == 0 ) s = 314159265 ;      /* default seed is 314159265 */

    int32_t seed = s ;
    // This is the initialization algorithm of F. James, widely in use for RANLUX.

    for ( int i = 0 ; i < 24 ; ++i ) {
      int32_t k = seed / 53668 ;
      seed = 40014 * (seed - k * 53668) - k * 12211 ;
      if ( seed < 0 ) seed += 2147483563 ;
      u[i] = seed % (uint32_t(1)<<24) ;
    }

    si   = 23 ;
    sj   = 9 ;
    sn   = 0 ;
    skip = luxury - 24 ;

    if ( u[23] & mask_hi ) carry = 1 ;
    else                   carry = 0 ;
  }

  uint32_t
  Luescher::incrementState() {
    uint32_t const mask_lo = UINT32(0x00ffffff) ; // 2^24 - 1
    uint32_t const mask_hi = ~mask_lo ;

    int32_t delta = u[sj] - u[si] - carry ;
    if ( delta & mask_hi ) { carry = 1 ; delta &= mask_lo ; }
    else                     carry = 0 ;

    u[si] = delta ;
    if ( si == 0 ) si = 23 ; else --si ;
    if ( sj == 0 ) sj = 23 ; else --sj ;
    return delta;
  }

  uint32_t
  Luescher::rand() {
    uint32_t skip = this -> skip ;
    uint32_t r    = incrementState();

    if ( ++sn == 24 ) {
      sn = 0 ;
      for ( uint32_t i = 0 ; i < skip ; ++i )
        incrementState() ;
    }
    return r;
  }
  
  /*   _                         _              ____  
  //  | |   _   _  ___  ___  ___| |__   ___ _ _|___ \ 
  //  | |  | | | |/ _ \/ __|/ __| '_ \ / _ \ '__|__) |
  //  | |__| |_| |  __/\__ \ (__| | | |  __/ |  / __/ 
  //  |_____\__,_|\___||___/\___|_| |_|\___|_| |_____|
  */
  #define RANLUX_STEP(x1,x2,i1,i2,i3)    \
          x1=xdbl[i1] - xdbl[i2];        \
          if (x2 < 0)                    \
          {                              \
            x1-=one_bit;                 \
            x2+=1;                       \
          }                              \
          xdbl[i3]=x2

  void
  Luescher2::seed( uint32_t s ) {
    double const one_bit = 1.0 / 281474976710656.0 ;  /* 1/2^48 */
    int xbit[31];

    if ( s == 0 ) s = 1 ; /* default seed is 1 */

    for ( int k = 0 ; k < 31 ; ++k ) { xbit[k] = s % 2 ; s /= 2 ; }

    int ibit = 0  ;
    int jbit = 18 ;
    for ( int k = 0 ; k < 12 ; ++k ) {
      double x = 0 ;
      for ( int m = 1 ; m <= 48 ; ++m ) {
        double y = (double) xbit[ibit] ;
        x += x + y ;
        xbit[ibit] = (xbit[ibit] + xbit[jbit]) % 2 ;
        ibit = (ibit + 1) % 31 ;
        jbit = (jbit + 1) % 31 ;
      }
      xdbl[k] = one_bit * x ;
    }

    carry  = 0 ;
    ir     = 0 ;
    jr     = 7 ;
    is     = 23 ;
    is_old = 0 ;
    pr     = luxury ;
  }
  
  void
  Luescher2::incrementState() {
  
    double const sbase    = 16777216.0 ;              /* 2^24 */
    double const sone_bit = 1.0 / 16777216.0 ;        /* 1/2^24 */
    double const one_bit  = 1.0 / 281474976710656.0 ; /* 1/2^48 */
    double const shift    = 268435456.0 ;             /* 2^28 */

    double x, y1, y2, y3;
    uint32_t k, m ;
    for ( k = 0 ; ir > 0 ; ++k ) {
      y1 = xdbl[jr] - xdbl[ir] ;
      y2 = y1 - carry ;
      if ( y2 < 0 ) { carry = one_bit ; y2 += 1 ; }
      else          { carry = 0 ; }
      xdbl[ir] = y2 ;
      ir = (ir+1) % 12 ;
      jr = (jr+1) % 12 ;
    }

    uint32_t kmax = pr - 12;

    for ( ; k <= kmax ; k += 12 ) {
      y1 = xdbl[7] - xdbl[0] ;
      y1 -= carry ;

      RANLUX_STEP (y2, y1, 8,  1,  0);
      RANLUX_STEP (y3, y2, 9,  2,  1);
      RANLUX_STEP (y1, y3, 10, 3,  2);
      RANLUX_STEP (y2, y1, 11, 4,  3);
      RANLUX_STEP (y3, y2, 0,  5,  4);
      RANLUX_STEP (y1, y3, 1,  6,  5);
      RANLUX_STEP (y2, y1, 2,  7,  6);
      RANLUX_STEP (y3, y2, 3,  8,  7);
      RANLUX_STEP (y1, y3, 4,  9,  8);
      RANLUX_STEP (y2, y1, 5,  10, 9);
      RANLUX_STEP (y3, y2, 6,  11, 10);

      if ( y3 < 0 ) { carry = one_bit ; y3 += 1 ; }
      else          { carry = 0 ; }
      xdbl[11] = y3;
    }

    kmax = pr;

    for ( ; k < kmax ; ++k ) {
      y1 = xdbl[jr] - xdbl[ir] ;
      y2 = y1 - carry ;
      if ( y2 < 0 ) { carry = one_bit ; y2 += 1 ; }
      else          { carry = 0 ; }
      xdbl[ir] = y2 ;
      ydbl[ir] = y2 + shift ;
      ir = (ir+1) % 12 ;
      jr = (jr+1) % 12 ;
    }

    ydbl[ir] = xdbl[ir] + shift ;

    for ( k = (ir+1) % 12 ; k > 0 ; k = (k+1) % 12 ) ydbl[k] = xdbl[k] + shift ;

    for ( k = 0, m = 0 ; k < 12 ; ++k ) {
      x  = xdbl[k];
      y2 = ydbl[k] - shift;
      if ( y2 > x ) y2 -= sone_bit;
      y1 = (x - y2) * sbase;

      xflt[m++] = (float) y1;
      xflt[m++] = (float) y2;
    }
    is = is_old = 2 * ir ;
  }

  uint32_t
  Luescher2::rand() {
    is = (is+1) % 24 ;
    if ( is == is_old ) incrementState() ;
    return uint32_t(xflt[is] * 16777216.0) ;
  }
  
  #undef RANLUX_STEP
  
}

// EOF Random.cc
