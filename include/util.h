/*  Last edited: Jul 22 11:31 2002 (klh) */
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

/*********************************************************************/
/********************** Error reporting ******************************/
/*********************************************************************/

void fatal_util( char *, ... );
void warning_util( char *, ...);


/*********************************************************************/
/********************** String functions *****************************/
/*********************************************************************/

typedef struct {
  char *buf;
  int buf_size;
} Line;

void free_Line( Line *);
Line *new_Line( void );
int read_Line( FILE *, Line * );


char *strdup_util (const char *str);


/**********************************************************************/
/*************** Core memory allocation wrappers **********************/
/**********************************************************************/

void free_util( void * );
void *malloc_util( size_t );
void *realloc_util( void *, size_t );
void *calloc_util( size_t, size_t );
void *malloc0_util(size_t );


/**********************************************************************/
/*************** dynamically growable arrays **************************/
/**********************************************************************/

typedef struct _Array	Array;

struct _Array
{
  char *data;
  int len;
};

#define append_val_Array(a,v)	append_vals_Array (a, &v, 1)
#define prepend_val_Array(a,v)  prepend_vals_Array (a, &v, 1)
#define insert_val_Array(a,i,v) insert_vals_Array (a, i, &v, 1)
#define index_Array(a,t,i)      (((t*) (a)->data) [(int)(i)])

void  free_Array (Array *, boolean );
Array *new_Array( int, boolean );
Array *append_vals_Array (Array *, const void *, int len);
Array *prepend_vals_Array (Array *,  const void *,  int len);
Array *insert_vals_Array (Array *array, int, const void *, int);
Array *set_size_Array (Array *array, int length);
Array *remove_index_Array (Array *, int);


/*********************************************************************/
/********************** Dictionary functions *************************/
/*********************************************************************/

typedef Array Dict;

signed char dict_lookup( Dict *, const char *);


/*********************************************************************/
/********************** Debug functions ******************************/
/*********************************************************************/
long int how_many_bytes (void);



#endif
