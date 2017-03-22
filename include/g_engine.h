/**********************************************************************
 ** File: engine.h
 * Author: Kevin Howe
 * Copyright (C) Genome Research Limited, 2002-
 *-------------------------------------------------------------------
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *-------------------------------------------------------------------
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
		   boolean, 
		   Gaze_Output *);

void backwards_calc(Gaze_Sequence *,
		    Gaze_Structure *,
		    boolean);

void scan_through_sources_dp(Gaze_Sequence *,
			     Gaze_Structure *,
			     int,
			     Gaze_DP_struct *,
			     boolean,
			     Gaze_Output *);

void scan_through_targets_dp(Gaze_Sequence *,
			     Gaze_Structure *,
			     int,
			     Gaze_DP_struct *,
			     boolean);

void scan_through_sources_for_max_only( Gaze_Sequence *,
					Gaze_Structure *,
					int,
					Gaze_DP_struct *,
					boolean);

void trace_back_general(Gaze_Sequence * );


#endif


