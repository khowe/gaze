/*  Last edited: Jul 24 14:24 2002 (klh) */
/**********************************************************************
 ** File: info.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#ifndef _GAZE_INFO
#define _GAZE_INFO

#include <stdio.h>
#include "util.h"

/*******************************************************
 The following three object represent a specific
 relationship between Target and either a specific single
 upstream feature, or all legal upstream features
*******************************************************/

/*******************************************************
		Killer_Feature_Qualifier
 These objects represent interruption constraints
 that legal source-target pairs must adhere to 
*******************************************************/

typedef struct {
  int feat_idx;
  boolean has_src_phase;
  boolean has_tgt_phase;
  int phase;
} Killer_Feature_Qualifier;

Killer_Feature_Qualifier *clone_Killer_Feature_Qualifier( Killer_Feature_Qualifier * );
void free_Killer_Feature_Qualifier( Killer_Feature_Qualifier *);
Killer_Feature_Qualifier *new_Killer_Feature_Qualifier( void );


/*******************************************************
		Killer_Feature_Qualifier
 These objects represent DNA constraints
 that legal source-target pairs must adhere to 
*******************************************************/

typedef struct {
  char src_dna;
  char tgt_dna;
} Killer_DNA_Qualifier;

Killer_DNA_Qualifier *clone_Killer_DNA_Qualifier( Killer_DNA_Qualifier *);
void free_Killer_DNA_Qualifier( Killer_DNA_Qualifier *);
Killer_DNA_Qualifier *new_Killer_DNA_Qualifier( void );


/*******************************************************
		Output_Qualifier 
 This object represents the output of a Gaze "region",
 i.e. the annotation of a piece of DNA lying between
 a pair of features.
*******************************************************/

typedef struct {
  char *feature;
  char *strand;
  char *frame;
  boolean need_to_print;
} Output_Qualifier;

Output_Qualifier *clone_Output_Qualifier( Output_Qualifier * );
void free_Output_Qualifier( Output_Qualifier * );
Output_Qualifier *new_Output_Qualifier( void );




/********************************************************
                   Features
*********************************************************/

typedef struct {
  double multiplier;
  int start_offset;
  int end_offset;
  boolean is_killer_feat;
  Array *targets;         /* of Feature_Relation*, idx by feat id */
  Array *sources;         /* of Feature_Relation*, idx by feat id */
  Array *kill_feat_quals; /* of Killer_Feature_Qualifier*         */
  Array *seg_quals;       /* of Segment_Qualifier*                */
  Output_Qualifier *out_qual;
} Feature_Info;                       

Feature_Info *empty_Feature_Info( void );
void free_Feature_Info(Feature_Info *);
Feature_Info *new_Feature_Info( int, int, double );


/*********************************************************
		   Segments
 General information about a segments and how they are
 used in speciic source-target instances
**********************************************************/

typedef struct {
  double multiplier;
  boolean use_projected;   /* default: use standard segs */
  boolean score_sum;       /* default: score max */
  boolean partial;         /* default: all overlapping segs */
} Segment_Info;                       

Segment_Info *empty_Segment_Info(void);
void free_Segment_Info(Segment_Info *seg_info);
Segment_Info *new_Segment_Info(double);



typedef struct {
  int seg_idx;
  boolean use_projected;   /* default: use standard segs */
  boolean score_sum;       /* default: score max */
  boolean is_exact_src;
  boolean is_exact_tgt;
  /* the following two attributes are mutually exclusive */
  boolean has_tgt_phase;
  boolean has_src_phase;
  int phase;
  boolean partial;       /* only use segments falling completely 
			      within regions of interest */
} Segment_Qualifier;

Segment_Qualifier *clone_Segment_Qualifier( Segment_Qualifier * );
void free_Segment_Qualifier(Segment_Qualifier *);
Segment_Qualifier *new_Segment_Qualifier( void );



/********************************************************
		Length Functions
*********************************************************/

typedef struct {
  double multiplier;
  Array *value_map;
  Array *raw_x_vals;
  Array *raw_y_vals;
  boolean becomes_monotonic;
  int monotonic_point;
} Length_Function;                       


#define apply_Length_Function(l,n) ((n>=l->value_map->len) ? \
    (index_Array(l->value_map,double,l->value_map->len-1)*(2+n-(l->value_map->len))) \
      - (index_Array(l->value_map,double,l->value_map->len-2)*(n-(l->value_map->len-1))) : \
    index_Array(l->value_map,double,n)) 

void calc_Length_Function(Length_Function *);
void free_Length_Function(Length_Function *);
Length_Function *new_Length_Function( double );
void scale_Length_Function(Length_Function *, double );

/**********************************************************
                Feature_Relation							   
 and finally, this object ties everthing together defines
 the relationship beween a source and target, from
 constraints, to scoring, to output
**********************************************************/

typedef struct {      
  int target;                   
  int source;
  int *min_dist;
  int *max_dist;
  int *phase;
  int *len_fun;
  Array *seg_quals;         /* of Segment_Qualifier        */
  Array *kill_feat_quals;   /* of Killer_Feature_Qualifier */
  Array *kill_dna_quals;    /* of Killer_DNA_Qualifier     */
  Output_Qualifier *out_qual;
} Feature_Relation;                     

Feature_Relation *clone_Feature_Relation(Feature_Relation *);
void free_Feature_Relation(Feature_Relation *ft_src);
Feature_Relation *new_Feature_Relation(void);



#endif
