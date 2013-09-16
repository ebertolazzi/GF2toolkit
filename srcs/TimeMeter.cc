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
 |      version: 0.1 04-01-2010                                             |
 |                                                                          |
 \*--------------------------------------------------------------------------*/

#include "TimeMeter.hh"
#include <stddef.h>

/*
//              _  _____ _                
//    __ _  ___| ||_   _(_)_ __ ___   ___ 
//   / _` |/ _ \ __|| | | | '_ ` _ \ / _ \
//  | (_| |  __/ |_ | | | | | | | | |  __/
//   \__, |\___|\__||_| |_|_| |_| |_|\___|
//   |___/                                
*/

#ifdef WIN32
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
  bool
  getTime( long & sec, long & usec ) {
    FILETIME   ft ;
    GetSystemTimeAsFileTime( &ft ) ;

    // FILETIME To UNIX Time
    unsigned __int64 u_sec = ft . dwHighDateTime ;
    u_sec <<= 32 ;
    u_sec |= ft . dwLowDateTime ;
    u_sec -= 116444736000000000ULL ;
    u_sec /= 10ULL ;

    sec  = long( u_sec / 1000000L ) ;
    usec = long( u_sec % 1000000L ) ;
    return true ;
  }
#else
  #include <sys/time.h>
  bool
  getTime( long & sec, long & usec ) {
    struct timeval now ;
    bool ok = gettimeofday(&now, NULL) == 0 ;
    if ( ok ) {
      sec  = now . tv_sec;
      usec = now . tv_usec;
    } else {
      sec = usec = 0 ;
    }
    return ok ;
  }
#endif

// EOF TimeMeter.cc
