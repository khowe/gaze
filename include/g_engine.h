/*  Last edited: Apr  2 15:13 2002 (klh) */
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

  GArray ***feats;   /* indices of features processed so far, organised by type */
  int ***fringes;    /* indices of last "significant" feature, organised first by
		        target type, then by source type, then by frame */
} Gaze_DP_struct;    

void free_Gaze_DP_struct( Gaze_DP_struct *, int );
Gaze_DP_struct *new_Gaze_DP_struct( int, int );


double calculate_path_score(GArray *, GArray *, Gaze_Structure *);

void forwards_calc(GArray *, 
		   GArray *,
		   Gaze_Structure *,
		   enum DP_Calc_Mode,
		   FILE *,
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
			     Gaze_Structure *,
			     int,
			     Gaze_DP_struct *,
			     enum DP_Calc_Mode,
			     enum DP_Traceback_Mode,
			     FILE *,
			     int,
			     FILE *);

void scan_through_targets_dp(GArray *,
			     GArray *,
			     Gaze_Structure *,
			     int,
			     Gaze_DP_struct *,
			     enum DP_Calc_Mode,
			     int,
			     FILE *);

void scan_through_sources_dp_old(GArray *,
				 GArray *,
				 int,
				 Gaze_DP_struct *,
				 Gaze_Structure *,
				 enum DP_Calc_Mode,
				 enum DP_Traceback_Mode,
				 int,
				 FILE *);

void scan_through_targets_dp_old(GArray *,
				 GArray *,
				 int,
				 Gaze_DP_struct *,
				 Gaze_Structure *,
				 enum DP_Calc_Mode,
				 int,
				 FILE *);



GArray *trace_back_general(GArray *,
			   GArray *,
			   Gaze_Structure *,
			   enum DP_Traceback_Mode);

#endif


