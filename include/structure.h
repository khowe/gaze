/*  Last edited: Jul 13 12:53 2002 (klh) */
/**********************************************************************
 ** File: structure.h
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description :
     This file takes care of the representation and parsing of the 
     XML encoded gaze structure file. 

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
#include "features.h"


/* For convenience, certain information about each feature 
   (namely the names of the features (feat_dict), and the dna
   to be taken for the features (take_dna) ) is situated 
   centrally here, rather than being distrubuted to the 
   appropriate Feature_Info object */

typedef struct {
  Array *feat_dict;
  Array *seg_dict;
  Array *len_fun_dict;
  Array *motif_dict;
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
void print_Gaze_Structure( Gaze_Structure *, FILE *);
void fill_in_Gaze_Structure( Gaze_Structure *);

#endif
