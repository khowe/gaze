/*  Last edited: Apr 23 12:45 2002 (klh) */
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
#include "engine.h"

/*******
  This object is a shopping bag of things necessary for
  GAZE to produce an output. It contains things that are
  general to all GAZE outputs 
********/

typedef struct {
  FILE *fh;
  char *seq_name;
  gboolean posterior;
  gboolean use_threshold;
  double threshold;
} Gaze_Output;

Gaze_Output *new_Gaze_Output( void );
void free_Gaze_Output( Gaze_Output * );

void print_GFF_Gaze_Features( Gaze_Output *, GArray *, Gaze_Structure * );
void print_GFF_path( Gaze_Output *, GArray *, Gaze_Structure * );

void write_GFF_header( Gaze_Output *, int, int );


#endif
