/*  Last edited: Aug  3 15:54 2002 (klh) */
/**********************************************************************
 ** File: sequence.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#ifndef _GAZE_SEQUENCE
#define _GAZE_SEQUENCE

#include "features.h"
#include "structure.h"

/********************************************************************/
/**************** Gaze_Sequence *************************************/
/********************************************************************/

typedef struct Gaze_Sequence{
  char *seq_name;
  char *dna_seq;
  StartEnd seq_region;

  Array *features;
  Array *segment_lists;
  Array *path;

  int offset_dna;
  
  Feature *beg_ft;
  Feature *end_ft;

  Array *min_scores;

} Gaze_Sequence;


void free_Gaze_Sequence( Gaze_Sequence * );
void initialise_Gaze_Sequence( Gaze_Sequence *, Gaze_Structure * );
Gaze_Sequence *new_Gaze_Sequence( char *, 
				  int, 
				  int ); 
				  

void get_features_from_dna( Gaze_Sequence *,
			    Array * );

void get_dna_for_features( Gaze_Sequence *,
			   Array *,
			   Dict * );

void remove_duplicate_features( Gaze_Sequence *);


/********************************************************************/
/**************** Gaze_Sequence_list ********************************/
/********************************************************************/

typedef struct {
  Gaze_Sequence **seq_list;
  int num_seqs;
  Dict *seq_id_dict;

} Gaze_Sequence_list;


void free_Gaze_Sequence_list ( Gaze_Sequence_list * );
Gaze_Sequence_list *new_Gaze_Sequence_list ( Array * );

void read_dna_seqs( Gaze_Sequence_list *, Array * );

void get_features_from_gff( Gaze_Sequence_list *,
			    Array *,
			    Array * ); 

boolean get_correct_features( Gaze_Sequence_list *,
			      Array *, 
			      Dict *,
			      boolean);

/********************************************************************/
/********************** Segment_list ********************************/
/********************************************************************/

typedef struct {
  Array *orig;    /* The orginal segment lists */
  Array *proj;    /* The projected segment lists */

  int reg_len;
  double *per_base[3];
} Segment_list;

void index_Segment_list (Segment_list * );
void project_Segment_list( Segment_list * );
void scale_Segment_list( Segment_list *, double );
void sort_Segment_list ( Segment_list *);
void free_Segment_list( Segment_list * );
Segment_list *new_Segment_list( int, int );


#endif
