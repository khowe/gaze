/*  Last edited: Oct  5 14:27 2001 (klh) */
/**********************************************************************
 ** File: params.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : Routines for the general dynamic programming engine 
                  for the GAZE system
 **********************************************************************/

#include "engine.h"



/*********************************************************************
 FUNCTION: calculate_segment_score
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
double calculate_segment_score( Feature *src, Feature *tgt, 
				GArray *segments, 
				Gaze_Structure *gs,
				Seg_Results *s_res) {
  GArray *seg_quals;
  int i,j;
  double score;

  /*
  Feature_Info *src_info = g_array_index( gs->feat_info, 
					  Feature_Info *,
					  src->feat_idx );
  */
  Feature_Info *tgt_info = g_array_index( gs->feat_info, 
					  Feature_Info *,
					  tgt->feat_idx );
  Feature_Relation *tgt_rel = g_array_index( tgt_info->sources, Feature_Relation *, src->feat_idx);

  /*
  int src_pos = src->real_pos.s + src_info->start_offset;
  int tgt_pos = tgt->real_pos.e - tgt_info->end_offset;
  */
  int src_pos = src->adj_pos.s; 
  int tgt_pos = tgt->adj_pos.e; 

  /* It is expensive to do this for every segment calc, but it has to be done. 
     Need to find way of making segment calc a whole lot more efficient. How
     about a sliding window that takes account of where we were last time?
  */

  for(i=0; i < s_res->has_score->len; i++) 
    g_array_index( s_res->has_score, gboolean, i ) = FALSE;

  if ( (seg_quals = tgt_rel->seg_quals) != NULL) {
    for (i=0; i < seg_quals->len; i++) {
      Segment_Qualifier *qual = g_array_index( seg_quals, Segment_Qualifier *, i);

      /* to do: each segment should have a type, either normal or projected.
	 for normal segs, the maximum original segment in the region being
	 considered is taken. For projected segments, the sum of all projected
	 segments in the region of interest is taken. Exact segments, by their
	 nature, are necessarily standard (not projected) segments. */

      if (qual != NULL) {
	Segment_lists *sl = g_array_index( segments, Segment_lists *, i);

	GArray *list = (qual->use_projected) ? sl->proj : sl->orig;
	GArray *segs;

	if (qual->has_tgt_phase) 
	  segs = g_array_index( list, GArray *, (tgt_pos - qual->phase + 1) % 3 );
	else if (qual->has_src_phase)
	  segs = g_array_index( list, GArray *, (src_pos + qual->phase) % 3 );
	else 
	  segs = g_array_index( list, GArray *, 3 );

	for (j=0; j < segs->len; j++) {
	  Segment *seg = g_array_index( segs, Segment *, j ); 
	  if ( seg->pos.s > tgt_pos )
	    break;
	  else if ( seg->pos.e < src_pos )
	    continue;
	  else {
	    int low = (seg->pos.s < src_pos)?src_pos:seg->pos.s;
	    int high = (seg->pos.e > tgt_pos)?tgt_pos:seg->pos.e;
	    
	    if ((! qual->is_exact_src || seg->pos.s == src_pos ) &&
		(! qual->is_exact_tgt || seg->pos.e == tgt_pos)) {
	      
	      score = seg->score * ( (double)(high - low + 1) / (double) (seg->pos.e - seg->pos.s + 1));
	      
	      if (! g_array_index( s_res->has_score, gboolean, i)) {
		g_array_index( s_res->raw_scores, double, i) = score;
		g_array_index( s_res->has_score, gboolean, i) = TRUE;
	      }
	      else {
		if (qual->score_sum)
		  /* now sum projected segment scores in a region rather than take max */
		  g_array_index( s_res->raw_scores, double, i) += score;		  
		else 
		  if (score > g_array_index( s_res->raw_scores, double, i ))
		    g_array_index( s_res->raw_scores, double, i) = score;
		
	      }
	    }
	  }
	}
      }
    }
  }

  score = 0.0;
  for(i=0; i < s_res->has_score->len; i++)
    if (g_array_index( s_res->has_score, gboolean, i))
	score += g_array_index( s_res->raw_scores, double, i );
     
  return score;
}


/*********************** Seg_Results *********************************/

/*********************************************************************
 FUNCTION: free_Seg_Results
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Seg_Results( Seg_Results *s_res ) {
  if (s_res != NULL) { 
    g_array_free( s_res->raw_scores, TRUE );
    g_array_free( s_res->has_score, TRUE );
    g_free( s_res );
  }
}


/*********************************************************************
 FUNCTION: new_Seg_Results
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Seg_Results *new_Seg_Results( int seg_dict_size ) {
  Seg_Results *s_res;

  s_res = (Seg_Results *) g_malloc (sizeof(Seg_Results));
  s_res->raw_scores = g_array_new( FALSE, TRUE, sizeof( double ) ); 
  s_res->has_score = g_array_new( FALSE, TRUE, sizeof( gboolean ) ); 
  g_array_set_size( s_res->raw_scores, seg_dict_size );
  g_array_set_size( s_res->has_score, seg_dict_size );

  return s_res;
}
