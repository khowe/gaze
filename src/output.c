/**********************************************************************
 ** File: output.c
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
 * Author : Kevin Howe
 * E-mail : klh@sanger.ac.uk
 * Description :
 **********************************************************************/

#include "output.h"


/*********************************************************************
 FUNCTION: Gaze_Output
 DESCRIPTION:
   Thisd 
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Gaze_Output *new_Gaze_Output( FILE *fh,
			      boolean probability,
			      boolean out_sample_gene,
			      boolean out_features,
			      boolean out_regions,
			      boolean use_threshold,
			      double threshold ) {
  Gaze_Output *out = (Gaze_Output *) malloc_util( sizeof( Gaze_Output ) );

  out->fh = fh;
  out->probability = probability;
  out->use_threshold = use_threshold;
  out->threshold = threshold;
  out->sample_gene = out_sample_gene;
  out->regions = out_regions;
  out->features = out_features;

  return out;
}


/*********************************************************************
 FUNCTION: free_Gaze_Output
 DESCRIPTION:
   An Output_line is a shopping bag of things needed to produce
   a line of GAZE output
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Gaze_Output( Gaze_Output *out ) {
  if (out != NULL) {
    free_util( out );
  }
}



/*********************************************************************
 FUNCTION: write_Gaze_path
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void write_Gaze_path( Gaze_Output *out,
		      Gaze_Sequence *g_seq,
		      Gaze_Structure *gs) {
  
  Feature *f1, *f2 = NULL;
  Feature_Info *f2_info;
  Feature_Relation *src;
  double feature_score, region_score;
  int i;

  write_GFF_comment( out->fh, 
		     "  Score of path : %.6f",
		     g_seq->end_ft->path_score );
  if (out->probability) {
    write_GFF_comment( out->fh, 
		       "  Forward score : %.6f", 
		       g_seq->end_ft->forward_score );
    write_GFF_comment( out->fh, 
		       "  Probability of path : %.6f", 
		       exp ( g_seq->end_ft->path_score -
			     g_seq->end_ft->forward_score) );
  }

  for( i=0; i < g_seq->path->len - 1; i++) {
    f1 = index_Array( g_seq->path, Feature *, i);
    f2 = index_Array( g_seq->path, Feature *, i+1);

    f2_info = index_Array( gs->feat_info, Feature_Info *, f2->feat_idx);

    src = index_Array( f2_info->sources, Feature_Relation *, f1->feat_idx );

    if (src == NULL) {
      fprintf( stderr, "Error: Illegal feature list\n");
      return;
    }

    /* first print feature itself, then print the region */
    feature_score = f1->score;
    if (out->probability)
      feature_score = exp( f1->forward_score + 
			   f1->backward_score - 
			   g_seq->beg_ft->backward_score);
    
    write_GFF_line( out->fh,
		    g_seq->seq_name,
		    "GAZE",
		    index_Array( gs->feat_dict, char *, f1->feat_idx ),
		    f1->real_pos.s, 
		    f1->real_pos.e,
		    feature_score,
		    NULL, NULL, NULL );

    region_score = f2->path_score - f1->path_score - f2->score; 
    
    if (out->probability) 
      region_score = exp( f1->forward_score +
			  f2->backward_score +
			  region_score + 
			  f2->score -
			  g_seq->beg_ft->backward_score );

    write_GFF_line( out->fh,
		    g_seq->seq_name,
		    "GAZE",
		    src->out_qual->feature,
		    f1->adj_pos.s, 
		    f2->adj_pos.e,
		    region_score,
		    src->out_qual->strand,
		    src->out_qual->frame,
		    NULL );
  }

  /* finally, the last feature */

  feature_score = f2->score;
  if (out->probability)
    feature_score = exp( f2->forward_score + 
			 f2->backward_score - 
			 g_seq->beg_ft->backward_score );
  if (f2 != NULL) {
    write_GFF_line( out->fh,
		    g_seq->seq_name, 
		    "GAZE",
		    index_Array( gs->feat_dict, char *, f2->feat_idx ),
		    f2->real_pos.s, 
		    f2->real_pos.e, 
		    feature_score,
		    NULL, NULL, NULL );
  }
}



/*********************************************************************
 FUNCTION: write_Gaze_Features
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void write_Gaze_Features( Gaze_Output *out,
			  Gaze_Sequence *g_seq,
			  Gaze_Structure *gs) {
  

  int i, discarded = 0;

  if (out->probability) {
    write_GFF_comment( out->fh, 
		       " GAZE feature scored by posterior probability");
    write_GFF_comment( out->fh, 
		       "   Fend = %.10f,     Bbegin = %.10f", 
		       g_seq->end_ft->forward_score,
		       g_seq->beg_ft->backward_score );
  }

  if (out->use_threshold)
    write_GFF_comment( out->fh, 
		       "  (only features scoring above %.2f are shown)", out->threshold );

  for(i=0; i < g_seq->features->len; i++) {
    Feature *f = index_Array( g_seq->features, Feature *, i);
    double score = f->score;

    if (out->probability) {
      score = exp( f->forward_score + 
		   f->backward_score - 
		   g_seq->beg_ft->backward_score );
    }

    if (! out->use_threshold || ! score < out->threshold) 
      write_GFF_line( out->fh,
		      g_seq->seq_name,
		      "GAZE",		      
		      index_Array( gs->feat_dict, char *,f->feat_idx ),
		      f->real_pos.s, 
		      f->real_pos.e, 
		      score,
		      NULL, NULL,
		      f->is_correct?"TRUE":NULL );
    else
      discarded++;
  }

  if (out->use_threshold)
    write_GFF_comment( out->fh,
		       " Discarded %d features with post. probs. %.2f or below", 
		       discarded, 
		       out->threshold );
}



/*********************************************************************
 FUNCTION: write_Gaze_header
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void write_Gaze_header( Gaze_Output *out, Gaze_Sequence *g_seq ) {

  write_GFF_header( out->fh, g_seq->seq_name, g_seq->seq_region.s, g_seq->seq_region.e );
  write_GFF_comment( out->fh, "source-version GAZE 1.0");

}

