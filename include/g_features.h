/*  Last edited: Aug  3 13:01 2002 (klh) */
/**********************************************************************
 ** File: g_features.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#ifndef _GAZE_FEATURES
#define _GAZE_FEATURES

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"
#include "gff.h"

/*************************************************************/
/************************ StartEnd ***************************/
/*************************************************************/

typedef struct {
  int s;
  int e;
} StartEnd;

void free_StartEnd(StartEnd *);
StartEnd *new_StartEnd(int, int);


/*************************************************************/
/********************** Feature ******************************/
/*************************************************************/

typedef struct {
  StartEnd real_pos;
  StartEnd adj_pos;
  short int feat_idx;
  short int dna;
  boolean is_selected;
  boolean is_antiselected;
  boolean is_correct;
  boolean invalid;
  int trace_pointer;
  double score;
  double path_score;
  double forward_score;
  double backward_score;

#ifdef PARAM  
  void *forward_uval;
  void *backward_uval;
#endif
} Feature;

Feature *clone_Feature(Feature *);
void free_Feature(Feature *);
Feature *new_Feature(void);
void write_Feature(FILE *, Feature *, Array *, Array *);


/*************************************************************/
/************************ Segment ****************************/
/*************************************************************/

typedef struct {
  char seg_idx;
  StartEnd pos;
  double score;
  int max_end_up;
  int max_end_up_idx;
} Segment;

Segment *clone_Segment(Segment *);
void free_Segment(Segment *);
Segment *new_Segment(void);
void write_Segment(Segment *, FILE *, Array *);
void index_Segments( Array * ); 
Array *project_Segments( Array * );


/*************************************************************/
/** Extraction of features and segments from GFF / DNA *******/
/*************************************************************/

typedef struct {
  int entity_idx;
  boolean has_score;
  double score;
} Gaze_entity;

void free_Gaze_entity(Gaze_entity *);
Gaze_entity *new_Gaze_entity(void);


/********************** From DNA *************************/

typedef struct {
  char *dna_motif;
  Array *features;     /* of DNA_to_entity * */
  Array *segments;     /* of DNA_to_entity * */
} DNA_to_Gaze_entities;

void free_DNA_to_Gaze_entities(DNA_to_Gaze_entities *);
DNA_to_Gaze_entities *new_DNA_to_Gaze_entities(void);

/********************** From GFF ************************/

typedef struct {
  char *gff_source;
  char *gff_feature;
  char *gff_strand;
  char *gff_frame;
  Array *features;
  Array *segments;
} GFF_to_Gaze_entities;

void free_GFF_to_Gaze_entities(GFF_to_Gaze_entities *);
GFF_to_Gaze_entities *new_GFF_to_Gaze_entities(void);

/***************** GFF parsing ***************************/

/**************** Sort comparison routines ************/
int order_features_for_dp(const void *, const void *);
int order_features(const void *, const void *);
int order_segments(const void *, const void *);


#endif


