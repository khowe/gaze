/*  Last edited: Jul 13 12:54 2002 (klh) */

/**********************************************************************
 ** File: p_engine.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description :
       Function and structures associated with the general dynamic
       programming engine for the GAZE system
 **********************************************************************/

#ifndef _ENGINE
#define _ENGINE

#include <glib.h>
#include <math.h>
#include "util.h"
#include "features.h"
#include "structure.h"

#define NEG_INFINITY -9999999.0

enum DP_Calc_Mode {
  NO_SUM,
  STANDARD_SUM,
  PRUNED_SUM
};

enum DP_Traceback_Mode {
  NO_TRACEBACK,
  MAX_TRACEBACK,
  SAMPLE_TRACEBACK
};


/****************** Seg_Results *************************/

typedef struct {
  GArray *raw_scores;
  GArray *has_score;
} Seg_Results;

void free_Seg_Results( Seg_Results * );
Seg_Results *new_Seg_Results( int );




double calculate_segment_score( Feature *, Feature *,
				GArray *, 
				Gaze_Structure *,
				Seg_Results * );
	
GArray *calculate_post_accuracies( GArray *, int, double);			

gboolean is_legal_path( GArray *, Gaze_Structure *, FILE *);

#endif
