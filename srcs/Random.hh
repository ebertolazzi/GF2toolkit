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

#ifndef RANDOM_HH
#define RANDOM_HH

#include <vector>
#include <cmath>
#include <string>

#ifdef MISSING_STDINT_H
  #include "pstdint.h"
#else
  #include <stdint.h>
#endif

// for srandom on unix machine
#include <stdlib.h>

#ifndef UINT16
  #define UINT16(A) A##u
#endif

#ifndef UINT32
  #ifdef UNSIGNED_IS_64BIT
    #define UINT32(A) A##ul
  #else
    #define UINT32(A) A##u
  #endif
#endif

#ifndef UINT64
  #define UINT64(A) A##ull
#endif

namespace RNG {

  typedef double valueType ;

  using namespace std ;

  //using std::string ;
  //using std::vector ;
  //using std::ostream ;

  /*            _                          
  //   ___  ___| |__  _ __ __ _  __ _  ___ 
  //  / __|/ __| '_ \| '__/ _` |/ _` |/ _ \
  //  \__ \ (__| | | | | | (_| | (_| |  __/
  //  |___/\___|_| |_|_|  \__,_|\__, |\___|
  //                            |___/      
  */
  static
  inline
  uint32_t
  schrage( uint32_t a, uint32_t b, uint32_t m) {
    /* This is a modified version of Schrage's method. It ensures that no
     * overflow or underflow occurs even if a=ceil(sqrt(m)). Usual 
     * Schrage's method works only until a=floor(sqrt(m)).
     */
    if ( a == UINT32(0) ) return UINT32(0) ;
    uint32_t q = m / a;
    uint32_t t = 2 * m - (m % a) * (b / q) ;
    if (t >= m) t -= m ;
    t += a * (b % q) ;
    return (t >= m) ? (t - m) : t ;
  }

  static
  inline
  uint32_t
  schrage_mult( uint32_t a, uint32_t b, uint32_t m, uint32_t sqrtm ) {
    /* To multiply a and b use Schrage's method 3 times.
     * represent a in base ceil(sqrt(m))  a = a1*ceil(sqrt(m)) + a0  
     * a*b = (a1*ceil(sqrt(m)) + a0)*b = a1*ceil(sqrt(m))*b + a0*b   
     */
    uint32_t t0 = schrage (sqrtm, b, m);
    uint32_t t1 = schrage (a / sqrtm, t0, m);
    uint32_t t2 = schrage (a % sqrtm, b, m);
    uint32_t t = t1 + t2;
    return (t >= m) ? (t - m) : t;
  }

  static uint32_t const two_pow_32_m_1 = ~(uint32_t(0)) ;
  static uint32_t const two_pow_31_m_1 = two_pow_32_m_1>>1 ;

  #define RNG_INHERITH \
  inline valueType rand01()   { return Rng::rand01()   ; } \
  inline valueType rand01o()  { return Rng::rand01o()  ; } \
  inline valueType rando01()  { return Rng::rando01()  ; } \
  inline valueType rando01o() { return Rng::rando01o() ; } \
  inline uint32_t rand( uint32_t a, uint32_t b ) { return Rng::rand( a, b ) ; }

  /*
  //double   rand53()   { return Rng::rand53()   ; } \
  //uint32_t rand31()   { return Rng::rand31()   ; } \
  */

  /*   ____              
  //  |  _ \ _ __   __ _ 
  //  | |_) | '_ \ / _` |
  //  |  _ <| | | | (_| |
  //  |_| \_\_| |_|\__, |
  //               |___/ 
  */
  class Rng {
  
    string const _name ;

    Rng( Rng const & ) ;
    Rng const & operator = ( Rng const & ) const ;

  protected:
    uint32_t  rndMax ;
    valueType to01 ;
    valueType to01open ;

  public:
  
    Rng( uint32_t rmax, string const & _n )
    : _name(_n), rndMax(rmax) {
      to01     = 1.0/rndMax ;
      to01open = 1.0/(1.0+rndMax) ;
    }

    virtual ~Rng() {}

    virtual void     seed( uint32_t s ) = 0 ;
    virtual uint32_t rand() = 0 ;

    /* generates a random number on [0,1]-real-interval */
    inline valueType rand01() { return this -> rand()*to01 ; }

    /* generates a random number on [0,1)-real-interval */
    inline valueType rand01o() { return this -> rand()*to01open ; }

    /* generates a random number on (0,1]-real-interval */
    inline valueType rando01() { return 1.0 - this -> rand()*to01open ; }

    /* generates a random number on (0,1)-real-interval */
    inline valueType rando01o() { return (this->rand() + 0.5)*to01open ; }

    //double rand53() ;

    inline 
    uint32_t
    rand( uint32_t a, uint32_t b )
    { return uint32_t(0.5 + a + this -> rand01()*(b-a)) ; }
    
    string const & name() const { return _name ; }

  };

  /*
  //   _   _       _      ____            _ 
  //  | | | |_ __ (_)_  _|  _ \ _ __   __| |
  //  | | | | '_ \| \ \/ / |_) | '_ \ / _` |
  //  | |_| | | | | |>  <|  _ <| | | | (_| |
  //   \___/|_| |_|_/_/\_\_| \_\_| |_|\__,_|
  */                                      
  class UnixRnd : public Rng {

    UnixRnd( UnixRnd const & ) ;
    UnixRnd const & operator = ( UnixRnd const & ) const ;
  
  public:
  
    UnixRnd( uint32_t s ) : Rng(two_pow_31_m_1,"Unix Random Number") {}

    virtual ~UnixRnd() {}

    virtual void seed( uint32_t s ) { srandom(s) ; }
    virtual uint32_t rand() { return static_cast<uint32_t>(random()) ; }
    
    RNG_INHERITH ;
  };
 
  /*   __  __                                    _____          _     _            
  //  |  \/  | ___ _ __ ___  ___ _ __  _ __   __|_   _|_      _(_)___| |_ ___ _ __ 
  //  | |\/| |/ _ \ '__/ __|/ _ \ '_ \| '_ \ / _ \| | \ \ /\ / / / __| __/ _ \ '__|
  //  | |  | |  __/ |  \__ \  __/ | | | | | |  __/| |  \ V  V /| \__ \ ||  __/ |   
  //  |_|  |_|\___|_|  |___/\___|_| |_|_| |_|\___||_|   \_/\_/ |_|___/\__\___|_|   
  */
  class MersenneTwister : public Rng {

    uint32_t mt[624] ;
    int      mti ;

    MersenneTwister( MersenneTwister const & ) ;
    MersenneTwister const & operator = ( MersenneTwister const & ) const ;
  
  public:
  
    MersenneTwister( uint32_t s ) ;
    virtual ~MersenneTwister() ;

    virtual void seed( uint32_t s ) ;
            void seed( uint32_t init_key[], int key_length ) ;

    virtual uint32_t rand() ;
    
    RNG_INHERITH ;
  };
  
  /*   ____                      _     
  //  | __ )  ___  _ __ ___  ___| |__  
  //  |  _ \ / _ \| '__/ _ \/ __| '_ \ 
  //  | |_) | (_) | | | (_) \__ \ | | |
  //  |____/ \___/|_|  \___/|___/_| |_|
  //   _   _ _          _                    _ _            
  //  | \ | (_) ___  __| | ___ _ __ _ __ ___(_) |_ ___ _ __ 
  //  |  \| | |/ _ \/ _` |/ _ \ '__| '__/ _ \ | __/ _ \ '__|
  //  | |\  | |  __/ (_| |  __/ |  | | |  __/ | ||  __/ |   
  //  |_| \_|_|\___|\__,_|\___|_|  |_|  \___|_|\__\___|_|   
  */
  /*
   * This generator is taken from
   *
   * Donald E. Knuth
   * The Art of Computer Programming
   * Volume 2
   * Third Edition
   * Addison-Wesley
   * Page 106-108
   *
   * It is called "Borosh - Niederreiter"
   *
   * This implementation copyright (C) 2001 Carlo Perassi and
   * (C) 2003 Heiko Bauke.
   */
  class BoroshNiederreiter : public Rng {
  
    uint32_t x ;

    BoroshNiederreiter( BoroshNiederreiter const & ) ;
    BoroshNiederreiter const & operator = ( BoroshNiederreiter const & ) const ;
  
  public:
  
    BoroshNiederreiter( uint32_t s ) : Rng(two_pow_32_m_1,"BoroshNiederreiter") { seed( s ) ; }
    virtual ~BoroshNiederreiter() {} ;

    virtual void seed( uint32_t s ) { x = s == 0 ? 1 : s ; }
    virtual uint32_t rand() {
      x = static_cast<uint32_t>((UINT64(1812433253) * x) & UINT64(0xFFFFFFFF)) ;
      return x ;
    }
    
    RNG_INHERITH ;
  };
  
  /*    ____ __  __ ____   ____ 
  //   / ___|  \/  |  _ \ / ___|
  //  | |   | |\/| | |_) | |  _ 
  //  | |___| |  | |  _ <| |_| |
  //   \____|_|  |_|_| \_\\____|
  */                        
  class CMRG : public Rng {
  
    int64_t const m1 ;
    int64_t const m2 ;

    uint64_t x1, x2, x3 ; // first component
    uint64_t y1, y2, y3 ; // second component

    CMRG( CMRG const & ) ;
    CMRG const & operator = ( CMRG const & ) const ;
  
  public:
  
    CMRG( uint32_t s )
    : Rng(4294967086u,"CMRG")
    , m1(4294967087u) // 2^32 - 309
    , m2(4294944443u) // 2^32 − 22853
    { seed( s ) ; }

    virtual ~CMRG() {} ;

    virtual void seed( uint32_t s ) ;
    virtual uint32_t rand() ;
    
    RNG_INHERITH ;
  };
  
  /*   _____ _     _                           ____       
  //  |  ___(_)___| |__  _ __ ___   __ _ _ __ |___ \__  __
  //  | |_  | / __| '_ \| '_ ` _ \ / _` | '_ \  __) \ \/ /
  //  |  _| | \__ \ | | | | | | | | (_| | | | |/ __/ >  < 
  //  |_|   |_|___/_| |_|_| |_| |_|\__,_|_| |_|_____/_/\_\
  */
  /*
   * This generator is taken from
   *
   * Donald E. Knuth
   * The Art of Computer Programming
   * Volume 2
   * Third Edition
   * Addison-Wesley
   * Page 108
   *
   * It is called "Fishman - L'Ecuyer"
   *
   * This implementation copyright (C) 2001 Carlo Perassi
   * and (C) 2003 Heiko Bauke.
   */
  class Fishman2x : public Rng {
  
    int32_t sx, sy, sz ;

    Fishman2x( Fishman2x const & ) ;
    Fishman2x const & operator = ( Fishman2x const & ) const ;
  
  public:
  
    Fishman2x( uint32_t s ) : Rng(UINT32(0x7FFFFFFE),"Fishman2x") { seed( s ) ; } // 2^31-2
    virtual ~Fishman2x() {} ;

    virtual void seed( uint32_t s ) ;
    virtual uint32_t rand() ;

    RNG_INHERITH ;
  } ;

  /*   _____ _     _                           _  ___  
  //  |  ___(_)___| |__  _ __ ___   __ _ _ __ / |( _ ) 
  //  | |_  | / __| '_ \| '_ ` _ \ / _` | '_ \| |/ _ \ 
  //  |  _| | \__ \ | | | | | | | | (_| | | | | | (_) |
  //  |_|   |_|___/_| |_|_| |_| |_|\__,_|_| |_|_|\___/ 
  */
  /*
   * This generator is taken from
   *
   * Donald E. Knuth
   * The Art of Computer Programming
   * Volume 2
   * Third Edition
   * Addison-Wesley
   * Page 106-108
   *
   * It is called "Fishman - Moore III".
   *
   * This implementation copyright (C) 2001 Carlo Perassi
   * and (C) 2003 Heiko Bauke.
   */
  class Fishman18 : public Rng {
  
    uint32_t x ;

    Fishman18( Fishman2x const & ) ;
    Fishman18 const & operator = ( Fishman18 const & ) const ;
  
  public:
  
    Fishman18( uint32_t s ) : Rng(UINT32(0x7FFFFFFD),"Fishman18") { seed( s ) ; } // 2^31-3 = 2147483645
    virtual ~Fishman18() {} ;

    virtual void seed( uint32_t s ) {
      x = s % 0x7FFFFFFFU ; // 2 ^ 31 - 1
      if ( x == 0 ) x = 1 ; // default seed is 1
    }

    virtual uint32_t rand() {
      uint32_t const AA           = UINT32(62089911) ;
      uint32_t const MM           = UINT32(0x7FFFFFFF) ; // 2 ^ 31 - 1
      uint32_t const CEIL_SQRT_MM = UINT32(46341) ; // ceil(sqrt(2 ^ 31 - 1))
      x = schrage_mult(AA, x, MM, CEIL_SQRT_MM) ;
      return x-1 ;
    }

    RNG_INHERITH ;
  } ;

  /*
  //   _____ _     _                           ____   ___  
  //  |  ___(_)___| |__  _ __ ___   __ _ _ __ |___ \ / _ \ 
  //  | |_  | / __| '_ \| '_ ` _ \ / _` | '_ \  __) | | | |
  //  |  _| | \__ \ | | | | | | | | (_| | | | |/ __/| |_| |
  //  |_|   |_|___/_| |_|_| |_| |_|\__,_|_| |_|_____|\___/ 
  */                                                     
  /*
   * This generator is taken from
   *
   * Donald E. Knuth
   * The Art of Computer Programming
   * Volume 2
   * Third Edition
   * Addison-Wesley
   * Page 108
   *
   * It is called "Fishman"
   *
   * This implementation copyright (C) 2001 Carlo Perassi
   * and (C) 2003 Heiko Bauke.
   */
  class Fishman20 : public Rng {
  
    uint32_t x ;

    Fishman20( Fishman2x const & ) ;
    Fishman20 const & operator = ( Fishman20 const & ) const ;
  
  public:
  
    Fishman20( uint32_t s ) : Rng(0x7FFFFFFD,"Fishman20") { seed( s ) ; } // 2^31-3 = 2147483645
    virtual ~Fishman20() {} ;

    virtual void seed( uint32_t s ) {
      x = s % 2147483647 ; // 2 ^ 31 - 1
      if ( x == 0 ) x = 1 ; // default seed is 1
    }

    virtual uint32_t rand() {
      uint32_t const m = 2147483647 ; 
      uint32_t const a = 48271 ;
      uint32_t const q = 44488 ;
      uint32_t const r = 3399 ;

      int32_t h = x / q;
      int32_t t = a * (x - h * q) - h * r ;

      if ( t < 0 ) x = t + m ;
      else         x = t ;
      return x - 1 ;
    }

    RNG_INHERITH ;
  } ;
  
  /*   _  __            _   _     
  //  | |/ /_ __  _   _| |_| |__  
  //  | ' /| '_ \| | | | __| '_ \ 
  //  | . \| | | | |_| | |_| | | |
  //  |_|\_\_| |_|\__,_|\__|_| |_|
  */
  /*
   * This generator is taken from
   *
   * Donald E. Knuth
   * The Art of Computer Programming
   * Volume 2
   * Third Edition
   * Addison-Wesley
   * Section 3.6
   *
   * The comments are taken from the book
   * Our comments are signed
   */
  class Knuth : public Rng {
  
    uint32_t idx ;
    uint32_t aa[2009]   ; // [Carlo]: I can't pass n to ran_array like Knuth does
    uint32_t ran_x[100] ; // the generator state

    Knuth( Knuth const & ) ;
    Knuth const & operator = ( Knuth const & ) const ;
  
  public:
  
    Knuth( uint32_t s ) : Rng(UINT32(0x3FFFFFFF),"Knuth") { seed( s ) ; } // 2^30-1
    virtual ~Knuth() {} ;

    virtual void     seed( uint32_t s ) ;
    virtual uint32_t rand() ;

    RNG_INHERITH ;
  } ;

  /*   _  __            _   _     ____  
  //  | |/ /_ __  _   _| |_| |__ |___ \ 
  //  | ' /| '_ \| | | | __| '_ \  __) |
  //  | . \| | | | |_| | |_| | | |/ __/ 
  //  |_|\_\_| |_|\__,_|\__|_| |_|_____|
  */
  /*
   * This generator is taken from
   *
   * Donald E. Knuth
   * The Art of Computer Programming
   * Volume 2
   * Third Edition
   * Addison-Wesley
   * Page 108
   *
   * This implementation  copyright (C) 2001 Carlo Perassi
   * and (C) 2003 Heiko Bauke.
   */
  class Knuth2 : public Rng {
  
    uint32_t x0, x1 ;

    Knuth2( Knuth2 const & ) ;
    Knuth2 const & operator = ( Knuth2 const & ) const ;
  
  public:
  
    Knuth2( uint32_t s ) : Rng(UINT32(0x7FFFFFFE),"Knuth2") { seed( s ) ; } // 2^31-1
    virtual ~Knuth2() {} ;

    virtual void seed( uint32_t s ) {
      uint32_t const MM = UINT32(0x7fffffff) ; // 2 ^ 31 - 1
      if ((s % MM) == 0) s = 1 ; // default seed is 1
      x0 = s % MM ;
      x1 = s % MM ;
    }

    virtual uint32_t rand() {
      uint32_t const AA1          = UINT32(271828183) ;
      uint32_t const AA2          = UINT32(1833324378) ; // = -314159269 mod (2 ^ 31 -1)
      uint32_t const MM           = UINT32(0x7fffffff) ; // 2 ^ 31 - 1
      uint32_t const CEIL_SQRT_MM = UINT32(46341) ;      // sqrt(2 ^ 31 - 1)
      uint32_t const xtmp = x1 ;
      x1 = schrage_mult(AA1, x1, MM, CEIL_SQRT_MM) +
           schrage_mult(AA2, x0, MM, CEIL_SQRT_MM) ;
      if ( x1 >= MM ) x1 -= MM ;
      x0 = xtmp ;
      return x1;
    }

    RNG_INHERITH ;
  } ;

  /*    ____ _____ ____  ____  _  _   
  //   / ___|  ___/ ___||  _ \| || |  
  //  | |  _| |_  \___ \| |_) | || |_ 
  //  | |_| |  _|  ___) |  _ <|__   _|
  //   \____|_|   |____/|_| \_\  |_|  
  */
  /*
   From Robert M. Ziff, "Four-tap shift-register-sequence
   random-number generators," Computers in Physics 12(4), Jul/Aug
   1998, pp 385-392.  A generalized feedback shift-register (GFSR)
   is basically an xor-sum of particular past lagged values.  A 
   four-tap register looks like:
      ra[nd] = ra[nd-A] ^ ra[nd-B] ^ ra[nd-C] ^ ra[nd-D]
   
   Ziff notes that "it is now widely known" that two-tap registers
   have serious flaws, the most obvious one being the three-point
   correlation that comes from the defn of the generator.  Nice
   mathematical properties can be derived for GFSR's, and numerics
   bears out the claim that 4-tap GFSR's with appropriately chosen
   offsets are as random as can be measured, using the author's test.

   This implementation uses the values suggested the the author's
   example on p392, but altered to fit the GSL framework.  The "state"
   is 2^14 longs, or 64Kbytes; 2^14 is the smallest power of two that
   is larger than D, the largest offset.  We really only need a state
   with the last D values, but by going to a power of two, we can do a
   masking operation instead of a modulo, and this is presumably
   faster, though I haven't actually tried it.  The article actually
   suggested a short/fast hack:

   #define RandomInteger (++nd, ra[nd&M]=ra[(nd-A)&M]\
                          ^ra[(nd-B)&M]^ra[(nd-C)&M]^ra[(nd-D)&M])

   so that (as long as you've defined nd,ra[M+1]), then you ca do things
   like: 'if (RandomInteger < p) {...}'.

   Note that n&M varies from 0 to M, *including* M, so that the
   array has to be of size M+1.  Since M+1 is a power of two, n&M
   is a potentially quicker implementation of the equivalent n%(M+1).

   This implementation copyright (C) 1998 James Theiler, based on
   the example mt.c in the GSL, as implemented by Brian Gough.
  */
  class GFSR4 : public Rng {

    int      nd ;
    uint32_t ra[1<<14] ; // 2^14

    GFSR4( GFSR4 const & ) ;
    GFSR4 const & operator = ( GFSR4 const & ) const ;
  
  public:
  
    GFSR4( uint32_t s ) : Rng(UINT32(0xFFFFFFFF),"GFSR4") { seed( s ) ; } // 2^32-1
    virtual ~GFSR4() {} ;

    virtual void     seed( uint32_t s ) ;
    virtual uint32_t rand() ;

    RNG_INHERITH ;
  } ;

  /*   __  __ ____   ____ 
  //  |  \/  |  _ \ / ___|
  //  | |\/| | |_) | |  _ 
  //  | |  | |  _ <| |_| |
  //  |_|  |_|_| \_\\____|
  */
  /*
   This is a fifth-order multiple recursive generator. The sequence is,

   x_n = (a_1 x_{n-1} + a_5 x_{n-5}) mod m

   with a_1 = 107374182, a_2 = a_3 = a_4 = 0, a_5 = 104480 and m = 2^31-1.

   We initialize the generator with x_n = s_n MOD m for n = 1..5,
   where s_n = (69069 * s_{n-1}) mod 2^32, and s_0 = s is the
   user-supplied seed.

   NOTE: According to the paper the seeds must lie in the range [0,
   2^31 - 2] with at least one non-zero value -- our seeding procedure
   satisfies these constraints.

   We then use 6 iterations of the generator to "warm up" the internal
   state.

   With this initialization procedure the theoretical value of
   z_{10006} is 2064828650 for s = 1. The subscript 10006 means (1)
   seed the generator with s = 1, (2) do the 6 warm-up iterations
   that are part of the seeding process, (3) then do 10000 actual
   iterations.

   The period of this generator is about 2^155.

   From: P. L'Ecuyer, F. Blouin, and R. Coutre, "A search for good
   multiple recursive random number generators", ACM Transactions on
   Modeling and Computer Simulation 3, 87-98 (1993).
   */
   class MRG : public Rng {

    int32_t x1, x2, x3, x4, x5 ;

    MRG( MRG const & ) ;
    MRG const & operator = ( MRG const & ) const ;
  
  public:
  
    MRG( uint32_t s ) : Rng(0x7FFFFFFE,"MRG") { seed( s ) ; } // 2^31-2 = 2147483646
    virtual ~MRG() {} ;

    virtual void     seed( uint32_t s ) ;
    virtual uint32_t rand() ;

    RNG_INHERITH ;
  } ;
  
  /*   _                         _               
  //  | |   _   _  ___  ___  ___| |__   ___ _ __ 
  //  | |  | | | |/ _ \/ __|/ __| '_ \ / _ \ '__|
  //  | |__| |_| |  __/\__ \ (__| | | |  __/ |   
  //  |_____\__,_|\___||___/\___|_| |_|\___|_|
  */
  /*
   This is a lagged fibonacci generator with skipping developed by Luescher.
   The sequence is a series of 24-bit integers, x_n, 

   x_n = d_n + b_n

   where d_n = x_{n-10} - x_{n-24} - c_{n-1}, b_n = 0 if d_n >= 0 and
   b_n = 2^24 if d_n < 0, c_n = 0 if d_n >= 0 and c_n = 1 if d_n < 0,
   where after 24 samples a group of p integers are "skipped", to
   reduce correlations. By default p = 199, but can be increased to
   365.

   The period of the generator is around 10^171. 

   From: M. Luescher, "A portable high-quality random number generator
   for lattice field theory calculations", Computer Physics
   Communications, 79 (1994) 100-110.

   Available on the net as hep-lat/9309020 at http://xxx.lanl.gov/

   See also,

   F. James, "RANLUX: A Fortran implementation of the high-quality
   pseudo-random number generator of Luscher", Computer Physics
   Communications, 79 (1994) 111-114

   Kenneth G. Hamilton, F. James, "Acceleration of RANLUX", Computer
   Physics Communications, 101 (1997) 241-248

   Kenneth G. Hamilton, "Assembler RANLUX for PCs", Computer Physics
   Communications, 101 (1997) 249-253 
  */
  class Luescher : public Rng {
  
    uint32_t luxury ;
    uint32_t si ;
    uint32_t sj ;
    uint32_t sn ;
    uint32_t skip ;
    uint32_t carry ;
    uint32_t u[24] ;

    Luescher( Luescher const & ) ;
    Luescher const & operator = ( Luescher const & ) const ;
    
    uint32_t incrementState() ;

  public:
  
    Luescher( uint32_t s, uint32_t l = 223 ) 
    : Rng(UINT32(0x00ffffff),"Luescher")
    , luxury(l)
    { seed( s ) ; } // 2^24-1, 389
    
    virtual ~Luescher() {} ;

    virtual void     seed( uint32_t s ) ;
    virtual uint32_t rand() ;

    RNG_INHERITH ;
  } ;
  /*   _                         _              ____  
  //  | |   _   _  ___  ___  ___| |__   ___ _ _|___ \ 
  //  | |  | | | |/ _ \/ __|/ __| '_ \ / _ \ '__|__) |
  //  | |__| |_| |  __/\__ \ (__| | | |  __/ |  / __/ 
  //  |_____\__,_|\___||___/\___|_| |_|\___|_| |_____|
  */                                                
  /* This is an implementation of M. Luescher's second generation
     version of the RANLUX generator.

     Thanks to Martin Luescher for providing information on this
     generator.

   */
  class Luescher2 : public Rng {
  
    uint32_t luxury ;
    double   xdbl[12] ;
    double   ydbl[12] ;  /* doubles first so they are 8-byte aligned */
    float    xflt[24] ;
    double   carry    ;
    uint32_t ir ;
    uint32_t jr ;
    uint32_t is ;
    uint32_t is_old ;
    uint32_t pr ;

    Luescher2( Luescher2 const & ) ;
    Luescher2 const & operator = ( Luescher2 const & ) const ;
    
    void incrementState() ;

  public:
  
    Luescher2( uint32_t s, uint32_t l = 109 ) //202, 397
    : Rng(UINT32(0x00ffffff),"Luescher2")
    , luxury(l)
    { seed( s ) ; } // 2^24-1, 389
    
    virtual ~Luescher2() {} ;

    virtual void     seed( uint32_t s ) ;
    virtual uint32_t rand() ;

    RNG_INHERITH ;
  } ;

}

#endif

// EOF Random.hh
