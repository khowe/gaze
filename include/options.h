/*  Last edited: Sep 14 14:19 2000 (klh) */
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
#include <glib.h>

/************************** constants *********************************/

#define MAX_DEF_LINE_SIZE 100

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

int get_option(int, char **, Option *, int, int *, char **, char **, gboolean *);
gboolean process_default_Options( FILE *defs, 
				  gboolean (*)( char *, char *) );



#endif

