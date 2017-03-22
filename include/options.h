/**********************************************************************
 * File: options.h
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
 * DESCRIPTION:
 * Rudimentary provision of command-line option processing
 **********************************************************************/
#ifndef _GETOPTIONS
#define _GETOPTIONS 

#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include "util.h"

/************************** constants *********************************/

#define MAX_DEF_LINE_SIZE 200

#define NO_ARGS    0 
#define INT_ARG    1
#define FLOAT_ARG  2
#define CHAR_ARG   3
#define STRING_ARG 4

/******************* structure definitions ****************************/

typedef struct {
  char *name;         /* name of option, e.g. "-option" */
  unsigned int type;  /* for typechecking, e.g. INT_ARG     */
} Option;


/******************* function prototypes ****************************/

int get_option(int, char **, Option *, int, int *, char **, char **, boolean *);
int process_default_Options( FILE *defs, 
			     boolean (*)( char *, char *) );



#endif

