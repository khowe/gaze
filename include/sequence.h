/*  Last edited: Aug  3 15:54 2002 (klh) */
/**********************************************************************
 ** File: sequence.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#ifndef _GAZE_SEQUENCE
#define _GAZE_SEQUENCE

#include "g_features.h"
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

  char *dna_file_name;
  char *gene_file_name;
  Array *gff_file_names;
  Array *selected_file_names;

} Gaze_Sequence;


void free_Gaze_Sequence( Gaze_Sequence *, boolean );
void initialise_Gaze_Sequence( Gaze_Sequence *, Gaze_Structure * );
Gaze_Sequence *new_Gaze_Sequence( char *, 
				  int, 
				  int ); 
				  
void convert_dna_Gaze_Sequence ( Gaze_Sequence *,
				 Array *,
				 Array *, 
				 Dict *);

void convert_gff_Gaze_Sequence( Gaze_Sequence *,
				Array *,
				Array * ); 


boolean get_correct_feats_Gaze_Sequence( Gaze_Sequence *,
					 Array *, 
					 Dict *,
					 boolean);

void read_dna_Gaze_Sequence( Gaze_Sequence *, Array * );

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

void read_dna_Gaze_Sequence_list( Gaze_Sequence_list *, Array * );

void convert_gff_Gaze_Sequence_list( Gaze_Sequence_list *,
				     Array *,
				     Array * ); 

boolean get_correct_feats_Gaze_Sequence_list( Gaze_Sequence_list *,
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
