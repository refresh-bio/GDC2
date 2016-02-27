/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#ifndef port_h
#define port_h
#include <limits.h>

#ifdef GCC
#define Inline inline
#else
#define Inline __inline
#endif

#if INT_MAX > 0x7FFF
typedef unsigned short uint2;  /* two-byte integer (large arrays)      */
typedef unsigned int   uint4;  /* four-byte integers (range needed)    */
#else
typedef unsigned int   uint2;
typedef unsigned long  uint4;
#endif

typedef unsigned int uint;     /* fast unsigned integer, 2 or 4 bytes  */

#endif