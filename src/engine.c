/*  Last edited: Nov  5 13:06 2001 (klh) */
/**********************************************************************
 ** File: params.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : Routines for the general dynamic programming engine 
                  for the GAZE system
 **********************************************************************/

#include "engine.h"


/*********************************************************************
 FUNCTION: calculate_post_accuracies
 DESCRIPTION:
    This function takes the given list of function, and measures the
    accuracies of the posterior probabilities of each feature. It
    does this by plotting a histogram (with the number of bins given)
    of the proportion of features that are correct for various ranges
    of posterior probability (idea: 70% of features that have a post.
    prob. of 0.7 should be correct).
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
GArray *calculate_post_accuracies( GArray *feats, int bins, double sigma ) {
  int i;
  GArray *res = g_array_new( FALSE, TRUE, sizeof( double ) );
  GArray *trues = g_array_new( FALSE, TRUE, sizeof( int ) );
  GArray *totals = g_array_new( FALSE, TRUE, sizeof( int  ) );


  g_array_set_size( res, bins );
  g_array_set_size( trues, bins );
  g_array_set_size( totals, bins );

  for(i=0; i < feats->len; i++) {
    Feature *ft = g_array_index( feats, Feature *, i);

    double post_prob = exp( ft->forward_score +
			    ft->backward_score -
			    g_array_index( feats, Feature *, 0)->backward_score);
    int index = (int) (post_prob * (double) bins);

    /* the following is to allow for the special case of a post prob of 1 */
    
    if (index == bins)
      index = bins - 1;

    g_array_index( totals, int, index ) = g_array_index( totals, int, index ) + 1;			    
    if (ft->is_correct) 
      g_array_index( trues, int, index ) = g_array_index( trues, int, index ) + 1;

  }
    
  for(i=0; i < res->len; i++) {
    if (g_array_index( totals, int, i ) > 0)
      g_array_index( res, double, i ) = 
	(double) g_array_index( trues, int, i ) /
	(double) g_array_index( totals, int, i );
  }

  g_array_free( trues, TRUE );
  g_array_free( totals, TRUE );

  /* put output in here for now. */

  printf( "Posterior probability accuracy plot - sigma = %.3f\n", sigma );
  printf( "Post prob     Prop. correct\n" );
  for (i=0; i < res->len; i++) {
    printf( "%4.3f:%4.3f\t%.3f\n", 
	    1.0 / (double) bins * i, 
	    1.0 / (double) bins * (i+1),
	    g_array_index( res, double, i ));
  }

  return res;
}



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

  Feature_Info *tgt_info = g_array_index( gs->feat_info, 
					  Feature_Info *,
					  tgt->feat_idx );
  Feature_Relation *tgt_rel = g_array_index( tgt_info->sources, Feature_Relation *, src->feat_idx);

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
