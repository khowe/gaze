/*  Last edited: Jul 13 14:19 2002 (klh) */
/**********************************************************************
 ** FILE: util.h
 ** NOTES:
 **   This file contains general utiliy functions used throughout the 
 **   application, such as those for memory management and error
 **   messaging. I have used it as a place to put other general stuff
 **   until I have somewhere better to put it.
 **********************************************************************/

#ifndef _GAZE_UTIL
#define _GAZE_UTIL

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>


#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#define malloc_util(x) malloc(x)
#define free_util(x) free(x)
#define realloc_util(x,y) realloc(x,y)

/*
void *malloc_util( size_t );
void *realloc_util( void *, size_t );
void *free_util( void * );
*/

void *calloc_util( size_t, size_t );
void *malloc0_util(size_t );
void fatal_util( char *, ... );
void warning_util( char *, ... );


#endif
