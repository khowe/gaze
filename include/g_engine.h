/*  Last edited: Oct  4 15:08 2001 (klh) */
/**********************************************************************
 ** File: engine.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#ifndef _G_ENGINE
#define _G_ENGINE

#include "engine.h"

typedef struct {
  double pth_score;
  int    pth_trace;
  double score;
  int last_selected;
  GArray *fringe_idx;
} Gaze_Results;

void free_Gaze_Results( Gaze_Results * );
Gaze_Results *new_Gaze_Results( int );


void calculate_path_score(GArray *, GArray *, Gaze_Structure *);

void forwards_calc(GArray *, 
		   GArray *,
		   Gaze_Structure *,
		   enum DP_Calc_Mode,
		   int,
		   FILE *);

void backwards_calc(GArray *, 
		    GArray *,
		    Gaze_Structure *,
		    enum DP_Calc_Mode,
		    int,
		    FILE *);

void scan_through_sources_dp(GArray *,
			     GArray *,
			     int,
			     Gaze_Results *,
			     Gaze_Structure *,
			     enum DP_Calc_Mode,
			     enum DP_Traceback_Mode,
			     int,
			     FILE *);

void scan_through_targets_dp(GArray *,
			     GArray *,
			     int,
			     Gaze_Results *,
			     Gaze_Structure *,
			     enum DP_Calc_Mode,
			     int,
			     FILE *);

GArray *trace_back_general(GArray *,
			   GArray *,
			   Gaze_Structure *,
			   enum DP_Traceback_Mode);


#endif


