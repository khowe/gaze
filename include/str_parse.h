/*  Last edited: Jul 13 12:53 2002 (klh) */
/**********************************************************************
 ** File: str_parse.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description :
     The parsing of the Gaze XML structure file
 **********************************************************************/


#ifndef _GAZE_STR_PARSE
#define _GAZE_STR_PARSE


#include <glib.h>
#include <stdlib.h>

#include "util.h"
#include "structure.h"
#include "info.h"
#include "features.h"
#include "expat.h"



Gaze_Structure *parse_Gaze_Structure( FILE * );


#endif
