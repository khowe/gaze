/**********************************************************************
 ** File: output.h
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
  boolean probability;
  boolean sample_gene;
  boolean regions;
  boolean features;
  boolean use_threshold;
  double threshold;
} Gaze_Output;

Gaze_Output *new_Gaze_Output( FILE *,
			      boolean,
			      boolean,
			      boolean,
			      boolean,
			      boolean,
			      double );
void free_Gaze_Output( Gaze_Output * );

void write_Gaze_Features( Gaze_Output *, Gaze_Sequence *, Gaze_Structure * );
void write_Gaze_path( Gaze_Output *, Gaze_Sequence *, Gaze_Structure * );

void write_Gaze_header( Gaze_Output *, Gaze_Sequence * );


#endif
