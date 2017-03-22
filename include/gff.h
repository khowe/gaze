/**********************************************************************
 ** File: gff.h
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
#ifndef _GAZE_GFF
#define _GAZE_GFF

#include "util.h"

typedef struct {
  Line *ln;
  char *seqname;
  char *source;
  char *type;
  int start;
  int end;
  double score;
  char *strand;
  char *frame;
  char *group;
} GFF_line;

void free_GFF_line( GFF_line * );
GFF_line *new_GFF_line( void );
int read_GFF_line(FILE *, GFF_line *);

void write_GFF_line( FILE *, char *, char *, char *, int, int, double, char *, char *, char *);
void write_GFF_header( FILE *, char *name, int start, int end );
void write_GFF_comment( FILE *, char *, ... );

#endif
