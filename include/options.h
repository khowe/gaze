/*  Last edited: Jul 22 19:02 2002 (klh) */
/**********************************************************************
 ** FILE: options.h
 ** DESCRIPTION:
 **  Rudimentary provision of command-line option processing
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

