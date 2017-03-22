
/**********************************************************************
 * File: p_engine.h
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
#ifndef _ENGINE
#define _ENGINE

#include <math.h>
#include "util.h"
#include "sequence.h"
#include "structure.h"

#define NEG_INFINITY -9999999.0
#define MOD3(x) ((x >= 0)?((x) % 3):((3 - (-(x) % 3))%3)) 

enum DP_Traceback_Mode {
  NO_TRACEBACK,
  MAX_TRACEBACK,
  SAMPLE_TRACEBACK
};

/****************** Index_list **************************/

struct _Index_list {
  boolean need_to_keep;
  int idx;
  struct _Index_list *next;
};

typedef struct _Index_list Index_list;


Index_list *add_to_Index_list( Index_list *, int, boolean );
Index_list *free_Index_list( Index_list *l, boolean);
Index_list *new_Index_list(void);
void traverse_Index_list( Index_list * );

/****************** Seg_Results *************************/

typedef struct {
  Array *raw_scores;
  Array *has_score;
  boolean has_exact_at_src;
  boolean has_exact_at_tgt;
  boolean exact_extends_beyond_tgt;
  boolean exact_extends_beyond_src;
} Seg_Results;

void free_Seg_Results( Seg_Results * );
Seg_Results *new_Seg_Results( int );




double calculate_segment_score(Gaze_Sequence *, 
			       Feature *, 
			       Feature *,
			       Gaze_Structure *,
			       Seg_Results * );
	
Array *calculate_post_accuracies( Array *, int, double);			

boolean is_legal_path( Array *, Gaze_Structure * );

#endif
