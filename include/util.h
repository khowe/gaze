/**********************************************************************
 * File: util.h
 * Author: Kevin Howe
 * Copyright (C) Genome Research Limited, 2002-
 *-------------------------------------------------------------------
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *-------------------------------------------------------------------
 * NOTES:
 * This file contains general utiliy functions used throughout the 
 * application, such as those for memory management and error
 * messaging. I have used it as a place to put other general stuff
 * until I have somewhere better to put it.
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

short int dict_lookup( Dict *, const char *);


/*********************************************************************/
/********************** Debug functions ******************************/
/*********************************************************************/
long int how_many_bytes (void);



#endif
