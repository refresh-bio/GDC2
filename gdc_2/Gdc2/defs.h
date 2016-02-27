/*
  This file is a part of GDC software distributed under GNU GPL 2 licence.
  The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
  
  Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
  
  Version: 2.0
  Date   : 2015-March-05
*/

#ifndef _DEFS_H
#define _DEFS_H

#define COMPRESS_DESC

#ifdef COMPRESS_DESC
#define my_getc(x)					gzgetc(x)
#define my_putc(x, y)				gzputc(y, x)
#define my_fread(a, b, c, d)		gzread(d, a, (b)*(c))
#define my_fwrite(a, b, c, d)		gzwrite(d, a, (b)*(c))
#define my_fclose(x)				gzclose(x)
#define my_feof(x)					gzeof(x)
#else
#define my_getc(x)					getc(x)
#define my_putc(x, y)				putc(x, y)
#define my_fread(a, b, c, d)		fread(a, b, c, d)
#define my_fwrite(a, b, c, d)		fwrite(a, b, c, d)
#define my_fclose(x)				fclose(x)
#define my_feof(x)					feof(x)
#endif


#ifdef WIN32
#define my_fopen	fopen
#define my_fseek	_fseeki64
#define my_ftell	_ftelli64
typedef unsigned short int uint16;
typedef int int32;
typedef short int int16;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;
typedef unsigned char uchar;
typedef unsigned char uchar_t;
#else
#define my_fopen	fopen
#define my_fseek	fseek
#define my_ftell	ftell
#define _TCHAR	char
#define _tmain	main

typedef unsigned short int uint16;
typedef int int32;
typedef short int int16;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;
typedef unsigned char uchar;
typedef unsigned char uchar_t;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#endif 

#endif
