/*  Last edited: Aug  3 15:24 2002 (klh) */
/**********************************************************************
 ** File: engine.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#ifndef _G_ENGINE
#define _G_ENGINE

#include "engine.h"
#include "output.h"

typedef struct {
  double pth_score;
  int    pth_trace;
  double score;
  int last_selected;

  Array ***feats;   /* indices of features processed so far, organised by type */
  int ***fringes;    /* indices of last "significant" feature, organised first by			
		        target type, then by source type, then by frame */
  Seg_Results *seg_res;
  
} Gaze_DP_struct;    

void free_Gaze_DP_struct( Gaze_DP_struct *, int );
Gaze_DP_struct *new_Gaze_DP_struct( int, int, int );


double calculate_path_score(Gaze_Sequence *, Gaze_Structure *);

void forwards_calc(Gaze_Sequence *,
		   Gaze_Structure *,
		   enum DP_Calc_Mode,
		   enum DP_Traceback_Mode,
		   Gaze_Output *);

void backwards_calc(Gaze_Sequence *,
		    Gaze_Structure *,
		    enum DP_Calc_Mode);

void scan_through_sources_dp(Gaze_Sequence *,
			     Gaze_Structure *,
			     int,
			     Gaze_DP_struct *,
			     enum DP_Calc_Mode,
			     enum DP_Traceback_Mode,
			     Gaze_Output *);

void scan_through_targets_dp(Gaze_Sequence *,
			     Gaze_Structure *,
			     int,
			     Gaze_DP_struct *,
			     enum DP_Calc_Mode);

void trace_back_general(Gaze_Sequence *,
			Gaze_Structure * );

#endif


