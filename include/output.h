/*  Last edited: Jul 22 14:58 2002 (klh) */
/**********************************************************************
 ** File: output.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#ifndef _GAZE_OUTPUT
#define _GAZE_OUTPUT

#include <math.h>
#include "util.h"
#include "gff.h"
#include "sequence.h"
#include "structure.h"
#include "engine.h"

/*******
  This object is a shopping bag of things necessary for
  GAZE to produce an output. It contains things that are
  general to all GAZE outputs 
********/

typedef struct {
  FILE *fh;
  boolean posterior;
  boolean use_threshold;
  double threshold;
} Gaze_Output;

Gaze_Output *new_Gaze_Output( FILE *,
			      boolean,
			      boolean,
			      double );
void free_Gaze_Output( Gaze_Output * );

void write_Gaze_Features( Gaze_Output *, Gaze_Sequence *, Gaze_Structure * );
void write_Gaze_path( Gaze_Output *, Gaze_Sequence *, Gaze_Structure * );

void write_Gaze_header( Gaze_Output *, Gaze_Sequence * );


#endif
