/**********************************************************************
 ** File: structure.h
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
 * Description :
 
     The Gaze structure file encodes three classes of information:

     1. Info about Features, and how they relate to each other
     2. Info about making features from DNA sequence
     3. Info about making features from a GFF file

 **********************************************************************/

#ifndef _GAZE_STRUCTURE
#define _GAZE_STRUCTURE

#include <stdlib.h>

#include "util.h"
#include "info.h"
#include "g_features.h"


/* For convenience, certain information about each feature 
   (namely the names of the features (feat_dict), and the dna
   to be taken for the features (take_dna) ) is situated 
   centrally here, rather than being distrubuted to the 
   appropriate Feature_Info object */

typedef struct {
  Dict *feat_dict;
  Dict *seg_dict;
  Dict *len_fun_dict;
  Dict *motif_dict;
  Array *feat_info;       /* of Feature_Info */
  Array *seg_info;        /* of Segment_Info */
  Array *length_funcs;    /* of Length_Function */
  Array *take_dna;        /* of StartEnd */
  Array *dna_to_feats;    /* of DNA_to_Features */
  Array *gff_to_feats;    /* of GFF_to_Features */

} Gaze_Structure;



/* Allocation and deallocation */

void free_Gaze_Structure( Gaze_Structure * );
Gaze_Structure *new_Gaze_Structure( void );
void write_Gaze_Structure( Gaze_Structure *, FILE *);
void fill_in_Gaze_Structure( Gaze_Structure *);

#endif
