/*  Last edited: Aug 30 15:15 2000 (klh) */
/**********************************************************************
 ** File: str_parse.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description :
     The parsing of the Gaze XML structure file
 **********************************************************************/


#ifndef _GAZE_STR_PARSE
#define _GAZE_STR_PARSE


#include <stdlib.h>
#include <glib.h>

#include "structure.h"
#include "info.h"
#include "features.h"
#include "xmlparse.h"



Gaze_Structure *parse_Gaze_Structure( FILE * );


#endif
