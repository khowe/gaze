/*  Last edited: Jul 22 11:43 2002 (klh) */
/**********************************************************************
 ** File: gff.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#ifndef _GAZE_GFF
#define _GAZE_GFF

#include "util.h"

typedef struct {
  Line *ln;
  char *seqname;
  char *source;
  char *type;
  char *strand;
  double score;
  int start;
  int end;
  char *group;
} GFF_line;

void free_GFF_line( GFF_line * );
GFF_line *new_GFF_line( void );
int read_GFF_line(FILE *, GFF_line *);

void write_GFF_line( FILE *, char *, char *, char *, int, int, double, char *, char *, char *);
void write_GFF_header( FILE *, char *name, int start, int end );
void write_GFF_comment( FILE *, char *, ... );

#endif
