/*  Last edited: Jul 15 11:55 2002 (klh) */
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

typedef unsigned char boolean;

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef	MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

#ifndef	MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif

#ifndef	ABS
#define ABS(a)	   (((a) < 0) ? -(a) : (a))
#endif

/*
#define malloc_util(x) malloc(x)
#define free_util(x) free(x)
#define realloc_util(x,y) realloc(x,y)
*/

void *malloc_util( size_t );
void *realloc_util( void *, size_t );
void *free_util( void * );


void *calloc_util( size_t, size_t );
void *malloc0_util(size_t );
void fatal_util( char *, ... );
void warning_util( char *, ... );

/* stuff for growable arrays */

/******* arrays **********************************************/

typedef struct _GRealArray  GRealArray;

struct _GRealArray
{
  char *data;
  int   len;
  int   alloc;
  int   elt_size;
  int   clear;
};

typedef struct _Array	Array;

struct _Array
{
  char *data;
  int len;
};

#define append_val_Array(a,v)	append_vals_Array (a, &v, 1)
#define prepend_val_Array(a,v)  prepend_vals_Array (a, &v, 1)
#define insert_val_Array(a,i,v) insert_vals_Array (a, i, &v, 1)
#define index_Array(a,t,i)      (((t*) (a)->data) [(i)])

Array* new_Array(int,
		 boolean );

void	free_Array (Array *,
		    boolean );

Array* append_vals_Array (Array *,
			  const void *,
			  int len);

Array* prepend_vals_Array (Array *,
			   const void *,
			   int len);

Array* insert_vals_Array (Array *array,
			  int,
			  const void *,
			  int);

Array* set_size_Array (Array *array,
		       int length);


Array *remove_index_Array (Array *,
			   int);

/******* strings ****************************/

char *strdup_util (const char *str);

long int how_many_bytes( void );


#endif
