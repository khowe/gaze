/*  Last edited: Oct  2 14:04 2001 (klh) */
/**********************************************************************
 ** File: output.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#ifndef _GAZE_OUTPUT
#define _GAZE_OUTPUT

#include <glib.h>
#include <math.h>
#include "features.h"
#include "structure.h"

void print_GFF_Gaze_Features( FILE *, GArray *, Gaze_Structure *, char *);
void print_GFF_path( FILE *, GArray *, Gaze_Structure *, char *);
void print_post_probs( FILE *,  GArray *, double, Gaze_Structure *, char *);



#endif
