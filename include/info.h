/*  Last edited: Oct  3 11:58 2001 (klh) */
/**********************************************************************
 ** File: info.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#ifndef _GAZE_INFO
#define _GAZE_INFO

#include <stdio.h>
#include <glib.h>


/*************************** data types  **********************/

typedef struct {
  gboolean has_phase;
  int phase;
} Killer_Feature_Qualifier;

Killer_Feature_Qualifier *clone_Killer_Feature_Qualifier( Killer_Feature_Qualifier * );
void free_Killer_Feature_Qualifier( Killer_Feature_Qualifier *);
Killer_Feature_Qualifier *new_Killer_Feature_Qualifier( void );


typedef struct {
  int src_dna;
  int tgt_dna;
} Killer_DNA_Qualifier;

Killer_DNA_Qualifier *clone_Killer_DNA_Qualifier( Killer_DNA_Qualifier *);
void free_Killer_DNA_Qualifier( Killer_DNA_Qualifier *);
Killer_DNA_Qualifier *new_Killer_DNA_Qualifier( void );



typedef struct {
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

/* General information about a particular feature type */

typedef struct {
  double multiplier;
  int start_offset;
  int end_offset;
  gboolean is_killer_feat;
  GArray *sources;                 /* of Feature_Relation *, idx by feat id */
  GArray *targets;                 /* of Feature_Relation *, idx by feat id */
  GArray *kill_feat_quals_up;             /* of int *, indexed by feature id */
  GArray *kill_feat_quals_down;           /* of int *, indexed by feature id */
  GArray *seg_quals;               /* of Segment_Qualifier *, indexed by seg id */
} Feature_Info;                       

Feature_Info *empty_Feature_Info( void );
void free_Feature_Info(Feature_Info *);
Feature_Info *new_Feature_Info( int, int, double );


/* Represents a length function */

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

/* the relationship of a feature with respect to a target */

typedef struct {      
  int target;                   
  int source;
  int *min_dist;
  int *max_dist;
  int *phase;
  int *len_fun;
  GArray *seg_quals;
  GArray *kill_feat_quals;
  GArray *kill_dna_quals;
  char *out_feature;
  char *out_strand;
  char *out_frame;
} Feature_Relation;                     

Feature_Relation *clone_Feature_Relation(Feature_Relation *);
void free_Feature_Relation(Feature_Relation *ft_src);
Feature_Relation *new_Feature_Relation(void);



/* General information about a particular segment type */

typedef struct {
  double multiplier;
} Segment_Info;                       

Segment_Info *empty_Segment_Info(void);
void free_Segment_Info(Segment_Info *seg_info);
Segment_Info *new_Segment_Info(double);

#endif
