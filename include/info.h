/*  Last edited: Apr 23 15:37 2002 (klh) */
/**********************************************************************
 ** File: info.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#ifndef _GAZE_INFO
#define _GAZE_INFO

#include <glib.h>
#include <stdio.h>


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
  gboolean has_src_phase;
  gboolean has_tgt_phase;
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
  int src_dna;
  int tgt_dna;
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
  gboolean need_to_print;
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
  gboolean is_killer_feat;
  GArray *targets;         /* of Feature_Relation*, idx by feat id */
  GArray *sources;         /* of Feature_Relation*, idx by feat id */
  GArray *kill_feat_quals; /* of Killer_Feature_Qualifier*         */
  GArray *seg_quals;       /* of Segment_Qualifier*                */
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
  gboolean use_projected;   /* default: use standard segs */
  gboolean score_sum;       /* default: score max */
} Segment_Info;                       

Segment_Info *empty_Segment_Info(void);
void free_Segment_Info(Segment_Info *seg_info);
Segment_Info *new_Segment_Info(double);



typedef struct {
  int seg_idx;
  gboolean use_projected;   /* default: use standard segs */
  gboolean score_sum;       /* default: score max */
  gboolean is_exact_src;
  gboolean is_exact_tgt;
  /* the following two attributes are mutually exclusive */
  gboolean has_tgt_phase;
  gboolean has_src_phase;
  int phase;
} Segment_Qualifier;

Segment_Qualifier *clone_Segment_Qualifier( Segment_Qualifier * );
void free_Segment_Qualifier(Segment_Qualifier *);
Segment_Qualifier *new_Segment_Qualifier( void );



/********************************************************
		Length Functions
*********************************************************/

typedef struct {
  double multiplier;
  GArray *value_map;
  GArray *raw_x_vals;
  GArray *raw_y_vals;
} Length_Function;                       


#define apply_Length_Function(l,n) ((n>=l->value_map->len) ? \
    (g_array_index(l->value_map,double,l->value_map->len-1)*(2+n-(l->value_map->len))) \
      - (g_array_index(l->value_map,double,l->value_map->len-2)*(n-(l->value_map->len-1))) : \
    g_array_index(l->value_map,double,n)) 

void calc_Length_Function(Length_Function *);
void free_Length_Function(Length_Function *);
Length_Function *new_Length_Function( double );


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
  GArray *seg_quals;         /* of Segment_Qualifier        */
  GArray *kill_feat_quals;   /* of Killer_Feature_Qualifier */
  GArray *kill_dna_quals;    /* of Killer_DNA_Qualifier     */
  Output_Qualifier *out_qual;
} Feature_Relation;                     

Feature_Relation *clone_Feature_Relation(Feature_Relation *);
void free_Feature_Relation(Feature_Relation *ft_src);
Feature_Relation *new_Feature_Relation(void);



#endif
