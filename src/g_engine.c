/*  Last edited: Oct 10 13:17 2001 (klh) */
/**********************************************************************
 ** File: engine.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#include "g_engine.h"


/*********************************************************************
 FUNCTION: free_Gaze_Results
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Gaze_Results( Gaze_Results *g_res ) {
  int i;

  if (g_res != NULL) { 
    for(i=0; i < g_res->fringe_idx->len; i++) 
      g_free( g_array_index( g_res->fringe_idx, int *, i) );
    g_array_free( g_res->fringe_idx, TRUE );
    g_free( g_res );
  }
}


/*********************************************************************
 FUNCTION: new_Gaze_Results
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Gaze_Results *new_Gaze_Results( int feat_dict_size ) {
  Gaze_Results *g_res;
  int i, j;

  g_res = (Gaze_Results *) g_malloc (sizeof(Gaze_Results));

  g_res->pth_score = g_res->score = 0.0;
  g_res->last_selected = -1;
  g_res->fringe_idx = g_array_new( FALSE, TRUE, sizeof( int *) ); 
  g_array_set_size( g_res->fringe_idx, feat_dict_size );
  for (i=0; i < g_res->fringe_idx->len; i++) {
    g_array_index( g_res->fringe_idx, int *, i) = 
      (int *) g_malloc( 3 * sizeof(int) );
    for(j=0; j < 3; j++)
      g_array_index( g_res->fringe_idx, int *, i)[j] = 0;
  }

  return g_res;
}



/*********************************************************************
 FUNCTION: calculate_path_score
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void calculate_path_score(GArray *path, 
			  GArray *segs,
			  Gaze_Structure *gs) {
  
  int idx, left_pos, right_pos, distance; 
  Feature_Info *tgt_info;
  Feature_Relation *reg_info;
  Feature *src, *tgt;
  double  trans_score, len_pen, seg_score;
  Seg_Results *s_res = new_Seg_Results( gs->seg_dict->len );

  trans_score = len_pen = seg_score = 0.0;

  for( idx=0; idx < path->len - 1; idx++) { 
    /* for the score to mean anything, all paths must begin with "BEGIN"
       and end with "END." Therefore ignoring the local score of the first 
       feature (done here) has no effect, since the score of "BEGIN" is 0
       (and has to be for the DP to work */

    src = g_array_index( path, Feature *, idx );
    
    tgt = g_array_index( path, Feature *, idx + 1 );
    tgt_info = g_array_index( gs->feat_info, Feature_Info *, tgt->feat_idx );
    
    if (! src->invalid && ! tgt->invalid) {
      /*
      left_pos = src->real_pos.s + src_info->start_offset;
      right_pos = tgt->real_pos.e - tgt_info->end_offset; */
      left_pos = src->adj_pos.s;
      right_pos = tgt->adj_pos.e;

      distance = right_pos - left_pos + 1;

      if ((reg_info = g_array_index(tgt_info->sources, Feature_Relation *, src->feat_idx)) != NULL) {

	if ((reg_info->phase == NULL) || (*(reg_info->phase) == distance % 3)) {
	  if ((reg_info->min_dist == NULL) || (*(reg_info->min_dist)) <= distance) {
	    if ((reg_info->max_dist == NULL) || (*(reg_info->max_dist)) >= distance) {

	      trans_score = 0.0;
	      len_pen = 0.0;
	      
	      seg_score = calculate_segment_score( src, tgt, segs, gs, s_res );
	      trans_score += seg_score;
	      
	      if (reg_info->len_fun != NULL) {
		Length_Function *lf = 
		  g_array_index(gs->length_funcs, Length_Function *, *(reg_info->len_fun));
		len_pen = apply_Length_Function( lf, distance );
	      }
	      trans_score -= len_pen;
	      /* tgt->path_score = src->path_score + trans_score + tgt->score; */
	      
	      tgt->path_score = src->path_score + trans_score + tgt->score;
	      
	    } 
	  }
	}
      } 
    }
  }

  free_Seg_Results( s_res );
}






/*********************************************************************
 FUNCTION: forwards_calc
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void forwards_calc( GArray *features,
		    GArray *segments,
		    Gaze_Structure *gs,
		    enum DP_Calc_Mode sum_mode,
		    int trace,
		    FILE *trace_fh) {
  
  int ft_idx;
  Gaze_Results *g_res = new_Gaze_Results( gs->feat_dict->len );

  if (trace > 1)
    fprintf(trace_fh, "\nForward calculation:\n\n");
  for (ft_idx = 1; ft_idx < features->len; ft_idx++) {
    scan_through_sources_dp( features, 
			     segments, 
			     ft_idx,  
			     g_res, 
			     gs, 
			     sum_mode,
			     MAX_TRACEBACK,
			     trace,
			     trace_fh);

    g_array_index( features, Feature *, ft_idx )->forward_score = g_res->score;
    g_array_index( features, Feature *, ft_idx )->path_score = g_res->pth_score;
#ifndef LOW_MEM
    g_array_index( features, Feature *, ft_idx )->trace_pointer = g_res->pth_trace;
#endif
  }

  free_Gaze_Results( g_res );
}



/*********************************************************************
 FUNCTION: backwards_calc
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
   This will run a whole lot faster if forwards_calc has been called
   first, invalidating lots of features. 
 *********************************************************************/
void backwards_calc( GArray *features,
		     GArray *segments,
		     Gaze_Structure *gs,
		     enum DP_Calc_Mode sum_mode,
		     int trace,
		     FILE *trace_fh) {
  int ft_idx, frame;
  Gaze_Results *g_res = new_Gaze_Results( gs->feat_dict->len );

  /* need to set the fringe indices, because they are not 0 for the
     backward calculation */
  for(ft_idx=0; ft_idx < g_res->fringe_idx->len; ft_idx++)
    for(frame=0; frame < 3; frame++)
      g_array_index( g_res->fringe_idx, int *, ft_idx)[frame] = features->len - 1;

  g_res->last_selected = features->len + 1;

  if (trace > 1)
    fprintf(trace_fh, "\nBackward calculation:\n\n");
  for (ft_idx = features->len-2; ft_idx >= 0; ft_idx--) {
    scan_through_targets_dp( features, 
			     segments, 
			     ft_idx,  
			     g_res,
			     gs,
			     sum_mode,
			     trace, 
			     trace_fh);

    g_array_index( features, Feature *, ft_idx )->backward_score = g_res->score;

  }

  free_Gaze_Results( g_res );
}



/*********************************************************************
 FUNCTION: scan_through_sources_dp
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void scan_through_sources_dp(GArray *features,
			     GArray *segments,
			     int tgt_idx,
			     Gaze_Results *g_res,
			     Gaze_Structure *gs,
			     enum DP_Calc_Mode sum_mode,
			     enum DP_Traceback_Mode trace_mode,
			     int trace,
			     FILE *trace_fh) {

  int src_idx, max_index, left_pos, right_pos, distance, k;
  Killer_Feature_Qualifier *kq;
  int local_fringe, last_fringe;
  double trans_score, viterbi_temp, forward_temp, len_pen, seg_score;

  gboolean touched_score;
  GArray *all_scores, *all_indices;
  GArray *poss_killer_feats;
  Feature_Info *src_info, *tgt_info;
  Feature_Relation *reg_info;
  Feature *src, *tgt;
  Seg_Results *seg_res;


  gboolean passed_killer_feat = FALSE;
  double max_score = NEG_INFINITY;
  double max_forward = NEG_INFINITY;
  double max_forpluslen = NEG_INFINITY;

  g_res->pth_score = 0.0;
  g_res->pth_trace = 0,0;
  g_res->score = 0.0;

  trans_score = viterbi_temp = forward_temp = len_pen = seg_score = 0.0;

  tgt = g_array_index( features, Feature *, tgt_idx );
  tgt_info = g_array_index( gs->feat_info, Feature_Info *, tgt->feat_idx );

  if (trace) {
    /*    fprintf( trace_fh, "Target %d %s %d %d %.3f (fringe %d)", tgt_idx,
	     g_array_index(gs->feat_dict, char *, tgt->feat_idx ), 
	     tgt->real_pos.s, tgt->real_pos.e, tgt->score,
	     g_array_index( g_res->fringe_idx, int *, tgt->feat_idx)[tgt->real_pos.s % 3]); */
    fprintf( trace_fh, "Target %d %s %d %d %.3f", tgt_idx,
	     g_array_index(gs->feat_dict, char *, tgt->feat_idx ), 
	     tgt->real_pos.s, tgt->real_pos.e, tgt->score); 
    if (trace > 1)
      fprintf( trace_fh, "\n" );
  }

  poss_killer_feats = g_array_new( FALSE, TRUE, sizeof(Feature *) );
  seg_res = new_Seg_Results( gs->seg_dict->len );

  if (sum_mode == STANDARD_SUM || sum_mode == PRUNED_SUM || trace_mode == SAMPLE_TRACEBACK) {
    all_scores = g_array_new( FALSE, TRUE, sizeof(double) );
    all_indices = g_array_new( FALSE, TRUE, sizeof(int) );
  }

  touched_score = FALSE;
  last_fringe = g_array_index( g_res->fringe_idx, int *, tgt->feat_idx)[tgt->real_pos.s % 3];
  for( src_idx = tgt_idx - 1;
       src_idx >= last_fringe
	 && (src_idx >= g_res->last_selected) 
	 && ! passed_killer_feat; 
       src_idx--) {

    src = g_array_index( features, Feature *, src_idx );
    src_info = g_array_index( gs->feat_info, Feature_Info *, src->feat_idx );

    if (trace > 1)
      fprintf( trace_fh, "  Source %d %s %d %d  ", src_idx,
	       g_array_index(gs->feat_dict, char *, src->feat_idx ), 
	       src->real_pos.s, src->real_pos.e );

    /* left_pos = src->real_pos.s + src_info->start_offset;
       right_pos = tgt->real_pos.e - tgt_info->end_offset; */
    left_pos = src->adj_pos.s;
    right_pos = tgt->adj_pos.e;

    distance = right_pos - left_pos + 1;
    
    if (! src->invalid ) {
      if (trace > 1)
	fprintf( trace_fh, "dist=%d  ", distance );

      if ((reg_info = g_array_index(tgt_info->sources, Feature_Relation *, src->feat_idx)) != NULL) {

	if (trace > 1)
	  fprintf( trace_fh, "TYPE_MATCH ");

	if ((reg_info->phase == NULL) || (*(reg_info->phase) == distance % 3)) {

	  if (trace > 1)
	    fprintf( trace_fh, "PHASE_MATCH ");

	  if ((reg_info->min_dist == NULL) || (*(reg_info->min_dist)) <= distance) {

	    if ((reg_info->max_dist == NULL) || (*(reg_info->max_dist)) >= distance) {
	      gboolean invalid_pair = FALSE;
	      
	      if (tgt_info->kill_feat_quals_up == NULL) {
		/* 
		   Looking at global upstream killers. If there are none, then
		   look for killer features specific to this pair of features
		*/
		
		if (reg_info->kill_feat_quals != NULL) {
		  
		  for(k=0; k < poss_killer_feats->len; k++) {
		    Feature *pos_kil = g_array_index( poss_killer_feats, Feature *, k );
		    if ( (kq = g_array_index( reg_info->kill_feat_quals,
					      Killer_Feature_Qualifier *,
					      pos_kil->feat_idx )) != NULL) {
		      /* Feature_Info *kill_info = g_array_index( gs->feat_info,
			 Feature_Info *,
			 pos_kil->feat_idx ); */
		      int kill_dist = pos_kil->adj_pos.e - left_pos + 1;
		      
		      if ((kill_dist >=0 && (kq->phase == (kill_dist %3))) || ! kq->has_phase) {
			invalid_pair = TRUE;
			break;
		      }
		    }
		  }
		}
	      }
	      
	      if (invalid_pair) {
		/* This src will always be killed by the same offending killer for all
		   targets of this type in the same frame. Therefore no need to update
		   fringe index */
		
		if (trace > 1)
		  fprintf( trace_fh, "KILLED_BY_STOP\n" );
	      }
	      else {
		/* check for DNA killers now */

		if (src->dna >= 0 && tgt->dna >= 0) {
		  if (reg_info->kill_dna_quals != NULL) {
		    for(k=0; k < reg_info->kill_dna_quals->len; k++) {
		      Killer_DNA_Qualifier *kdq = g_array_index( reg_info->kill_dna_quals, 
								 Killer_DNA_Qualifier *,
								 k );
		      if (src->dna == kdq->src_dna && tgt->dna == kdq->tgt_dna) {
			invalid_pair = TRUE;
			break;
		      }
		    }
		  }
		}
		
		/* This src may not be killed for future targets of this type in the same
		   frame - therefore need to update fringe index
		*/
		
		if (invalid_pair) {
		  if (sum_mode == PRUNED_SUM)
		    local_fringe = src_idx;
		  
		  if (trace > 1)
		    fprintf( trace_fh, "KILLED_BY_DNA\n" );
		  
		}
	      }
	      
	      if (! invalid_pair) {
		/* a legal pair of features, so calculate the score */
		trans_score = 0.0;
		len_pen = 0.0;
		
		seg_score = calculate_segment_score( src, tgt, segments, gs, seg_res );
		trans_score += seg_score;
		
		if (reg_info->len_fun != NULL) {
		  Length_Function *lf = 
		    g_array_index(gs->length_funcs, Length_Function *, *(reg_info->len_fun));
		  len_pen = apply_Length_Function( lf, distance );
		}
		trans_score -= len_pen;
		
		viterbi_temp = src->path_score +
		  + trans_score
		  + tgt->score;
		
		if (! touched_score || (viterbi_temp > max_score) ) {
		  max_score = viterbi_temp;
		  max_index = src_idx;
		}
		
		if (sum_mode == STANDARD_SUM || sum_mode == PRUNED_SUM || trace_mode == SAMPLE_TRACEBACK) {
		  forward_temp = src->forward_score 
		    + trans_score
		    + tgt->score;
		  
		  /************ DEBUG: ************/
		  /* forward_temp = src->forward_score;  */
		  /************ DEBUG: ************/

		  g_array_append_val( all_scores, forward_temp);
		  g_array_append_val( all_indices, src_idx);
		  
		  if (! touched_score || (forward_temp > max_forward))
		    max_forward = forward_temp;
		  
		  if (sum_mode == PRUNED_SUM) {
		    
		    if (! touched_score ) {
		      /* add back in the length penalty, because when judging for dominance, 
			 the length penalty will be different for future features */
		      max_forpluslen = forward_temp + len_pen;
		      local_fringe = src_idx;
		    }
		    else {
		      /* compare this one to max_forward, to see if it is dominated */
		      if (forward_temp + len_pen > max_forpluslen) 
			max_forpluslen = forward_temp + len_pen;
		      
		      if ( (max_forpluslen - (forward_temp + len_pen) < 20.0)
			   || (max_index == src_idx) )
			local_fringe = src_idx;
		      
		      /* 
			 this feature is not dominated if it has a "significant" contribution
			 to the forward score, OR it is the maximum source (sometimes the maximum
			 source contributes insignificantly to the forward score!)
		      */
		    }
		  }
		}
		
		touched_score = TRUE;
		
		if (trace > 1) 
		  fprintf( trace_fh, "scre: v=%.3f, f=%.8f (seg:%.3f len:%.3f)\n",
			   viterbi_temp, forward_temp, seg_score, len_pen );
		
	      }
 
	    } /* if max dist */
	    else {
	      /* This src will always be too distant for future targets of this type in the same
		 frame - therefore NO need to update fringe index
	      */
	      if (trace > 1)
		fprintf( trace_fh, "TOO DISTANT\n" );
	    }
	  } /* if min dist */
	  else {
	    /* This src may not be too close for future targets of this type in the same
	       frame - therefore need to update fringe index
	    */
	    if (sum_mode == PRUNED_SUM)
	      local_fringe = src_idx;
	    
	    if (trace > 1)
	      fprintf( trace_fh, "TOO CLOSE\n" );
	  }
	  
	} /* if mod3 */
	else {
	  /* This src will always be out of phase with all targets of this type in the same 
	     frame. Therefore no need to update fringe index */
	  
	  if (trace > 1)
	    fprintf( trace_fh, "PHASE-MISMATCH\n" );
	}

      } /* if match */
      else {
	  /* This src will always be out of type with all targets of this type in the same 
	     frame. Therefore no need to update fringe index */

	if (trace > 1)
	  fprintf( trace_fh, "TYPE-MISMATCH\n" );
      }

    } /* if valid */
    else {
      if (trace > 1)
	fprintf( trace_fh, "INVALID\n" );
    }

    if ( tgt_info->kill_feat_quals_up != NULL) {
      if ((kq = g_array_index( tgt_info->kill_feat_quals_up, 
			       Killer_Feature_Qualifier *, 
			       src->feat_idx)) != NULL) {
	if ((distance >= 0 && (kq->phase == (distance % 3))) || ! kq->has_phase)
	  passed_killer_feat = TRUE;
      }
    }
    
    if (src_info->is_killer_feat)
      g_array_append_val( poss_killer_feats, src );

  } /* for each source */         


  /* 
   update the position of the last forced feature. However, this
   needs to be the position of the first feature in the last
   forced 'block' (where 5'0, 5'1, 5'2 form a force block, for
   example). This is to get around the bother of the fact that
   a single splice site in the source data ends up as three
   splice sites in the feature list
  */

  if (tgt->is_selected) {
    if (g_res->last_selected < 0)
      g_res->last_selected = tgt_idx;
    else {
      Feature *last_feat = g_array_index( features, 
					  Feature *, 
					  g_res->last_selected );
      if (last_feat->real_pos.s != tgt->real_pos.s || last_feat->real_pos.e != tgt->real_pos.e) 
	g_res->last_selected = tgt_idx;
    }
  }

  if (touched_score) {
    g_res->pth_score = max_score;

    if (sum_mode == STANDARD_SUM || sum_mode == PRUNED_SUM) {
      /* the trick of subtracting the max before exponentiating avoids
	 overflow errors. Just need to add it back when logging back down */
      for (src_idx=0; src_idx < all_scores->len; src_idx++) {
	g_res->score += exp( g_array_index(all_scores, double, src_idx) 
				- max_forward );  

      }
      g_res->score = log( g_res->score ) + max_forward;

      /* if we passed killer, then passes for future instances of this target type
	 in the same absolute frame will be terminated by the same killer. Therefore,
	 we do not have to update the fringe index */
      if (sum_mode == PRUNED_SUM && local_fringe > 0 && ! passed_killer_feat)
	g_array_index( g_res->fringe_idx, int *, tgt->feat_idx )[tgt->real_pos.s % 3] = local_fringe;
    }

    if (trace_mode == SAMPLE_TRACEBACK) {
      double random_number = (double) rand() / (double) RAND_MAX;
      double sum = 0.0;
      for(src_idx=0; src_idx < all_indices->len; src_idx++) {
	double ft_prob = exp( g_array_index( all_scores, double, src_idx ) -
			      tgt->forward_score ); 
	sum += ft_prob;

	if (sum >= random_number) {
	  Feature *tmp;

	  g_res->pth_trace = g_array_index( all_indices, int, src_idx);

	  tmp = g_array_index(features, Feature *, g_res->pth_trace ); 

	  /* note that we are just returning the transition + local score 
	     here for the source-target pair. This can be accumulated
	     by the caller to get the score of the sample path. Alternatively,
	     the score of the path can be computed with calculate_path_score */

	  g_res->pth_score = g_array_index( all_scores, double, src_idx) -
	    tmp->forward_score;

	  break;
	}
      }
    }
    else if (trace_mode == MAX_TRACEBACK)
      g_res->pth_trace = max_index;

    if (trace)
      fprintf(trace_fh, "  RESULT: v=%.3f, max=%d, f=%.8f\n", 
	      g_res->pth_score,
	      max_index,
	      g_res->score );
  }
  else {
    if (trace)
      fprintf( trace_fh, "  *** Invalidating\n");

   /* We don't need to set the score to negative infinity for the 
       calculation of Fend - however we need to do it for the sake
       of individual feature posterior probabilities */

    g_res->score = NEG_INFINITY; 
    tgt->invalid = TRUE;
  }

  if (sum_mode == STANDARD_SUM || sum_mode == PRUNED_SUM ||  trace_mode == SAMPLE_TRACEBACK) {
    g_array_free( all_indices, TRUE );
    g_array_free( all_scores, TRUE );
  }
  free_Seg_Results( seg_res );
  g_array_free( poss_killer_feats, TRUE );


}      



/*********************************************************************
 FUNCTION: scan_through_targets_dp
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void scan_through_targets_dp(GArray *features,
			     GArray *segments,
			     int src_idx,
			     Gaze_Results *g_res,
			     Gaze_Structure *gs,
			     enum DP_Calc_Mode sum_mode,
			     int trace,
			     FILE *trace_fh) {

  int tgt_idx, left_pos, right_pos, distance, k; 
  Killer_Feature_Qualifier *kq;
  int local_fringe, last_fringe;
  double trans_score, len_pen, seg_score, backward_temp;

  gboolean touched_score;
  GArray *poss_killer_feats;
  GArray *all_scores;
  Feature_Info *tgt_info, *src_info;
  Feature_Relation *reg_info;
  Feature *src, *tgt;
  Seg_Results *seg_res;

  gboolean passed_killer_feat = FALSE;
  double max_backward = 0.0;
  double max_backpluslen = 0.0;

  g_res->score = 0.0;
  trans_score = len_pen = seg_score = backward_temp = 0.0;

  src = g_array_index( features, Feature *, src_idx );
  src_info = g_array_index( gs->feat_info, Feature_Info *, src->feat_idx );

  if (trace) {
    /*    fprintf( trace_fh, "Source %d %s %d %d %.3f (fringe %d)", src_idx,
	     g_array_index(gs->feat_dict, char *, src->feat_idx ), 
	     src->real_pos.s, src->real_pos.e, src->score,
	     g_array_index( g_res->fringe_idx, int *, src->feat_idx)[src->real_pos.s % 3]); */
    fprintf( trace_fh, "Source %d %s %d %d %.3f", src_idx,
	     g_array_index(gs->feat_dict, char *, src->feat_idx ), 
	     src->real_pos.s, src->real_pos.e, src->score ); 
    
    if (trace > 1)
      fprintf( trace_fh, "\n" );
  }

  poss_killer_feats = g_array_new( FALSE, TRUE, sizeof(Feature *));
  seg_res = new_Seg_Results( gs->seg_dict->len );
  all_scores = g_array_new( FALSE, TRUE, sizeof(double));

  if (! src->invalid ) { 
    touched_score = FALSE;
    last_fringe = g_array_index( g_res->fringe_idx, int *, src->feat_idx)[src->real_pos.s % 3];
    for( tgt_idx = src_idx + 1;
	 tgt_idx <= last_fringe
	   && (tgt_idx < g_res->last_selected) 
	   && ! passed_killer_feat; 
	 tgt_idx++) {

      tgt = g_array_index( features, Feature *, tgt_idx );
      tgt_info = g_array_index( gs->feat_info, Feature_Info *, tgt->feat_idx );
      
      if (trace > 1)
	fprintf( trace_fh, "  Target %d %s %d %d  ", tgt_idx,
		 g_array_index(gs->feat_dict, char *, tgt->feat_idx ), 
		 tgt->real_pos.s, tgt->real_pos.e );

      /* left_pos = src->real_pos.s + src_info->start_offset;
	 right_pos = tgt->real_pos.e - tgt_info->end_offset; */
      left_pos = src->adj_pos.s;
      right_pos = tgt->adj_pos.e;
      
      distance = right_pos - left_pos + 1;

      if (! tgt->invalid ) {
	if (trace > 1)
	  fprintf( trace_fh, "dist=%d  ", distance );
	
	if ((reg_info = g_array_index(tgt_info->sources, Feature_Relation *, src->feat_idx)) != NULL) {

	  if (trace > 1)
	    fprintf( trace_fh, "TYPE_MATCH ");
	  
	  if ((reg_info->phase == NULL) || (*(reg_info->phase) == distance % 3)) {
	    if (trace > 1)
	      fprintf( trace_fh, "PHASE_MATCH ");
	    
	    if ((reg_info->min_dist == NULL) || (*(reg_info->min_dist)) <= distance) {

	      if ((reg_info->max_dist == NULL) || (*(reg_info->max_dist)) >= distance) {
		gboolean invalid_pair = FALSE;

		if (src_info->kill_feat_quals_down == NULL) {
		  /* 
		     Looking at global downstream killers. If there are none, then
		     look for killer features specific to this pair of features
		  */
		  
		  if (reg_info->kill_feat_quals != NULL) {

		    for(k=0; k < poss_killer_feats->len; k++) {
		      Feature *pos_kil = g_array_index( poss_killer_feats, Feature *, k );
		      if ((kq = g_array_index( reg_info->kill_feat_quals,
					       Killer_Feature_Qualifier *,
					       pos_kil->feat_idx )) != NULL) {
			/* Feature_Info *kill_info = g_array_index( gs->feat_info,
			   Feature_Info *,
			   pos_kil->feat_idx ); */
			int kill_dist = right_pos - pos_kil->adj_pos.s + 1;
			
			if ((kill_dist >=0 && (kq->phase == (kill_dist %3))) || ! kq->has_phase) {
			  invalid_pair = TRUE;
			  break;
			}
		      }
		    }
		  }
		}
		
		if (invalid_pair) {
		  /* this source will always be killed by the offending killer for
		     all future instances of this target. Therefore do NOT update
		     the fringe index */
		  if (trace > 1)
		    fprintf( trace_fh, "KILLED BY STOP\n" );
		}
		else {
		  /* check for DNA killers now */
		  
		  if (src->dna >= 0 && tgt->dna >= 0) {
		    if (reg_info->kill_dna_quals != NULL) {
		      for(k=0; k < reg_info->kill_dna_quals->len; k++) {
			Killer_DNA_Qualifier *kdq = g_array_index( reg_info->kill_dna_quals, 
								   Killer_DNA_Qualifier *,
								   k );
			if (src->dna == kdq->src_dna && tgt->dna == kdq->tgt_dna) {
			  invalid_pair = TRUE;
			  break;
			}
		      }
		    }
		  }

		  /* This src may not be killed for future targets of this type in the same
		     frame - therefore need to update fringe index
		  */
		  
		  if (invalid_pair) {
		    if (sum_mode == PRUNED_SUM)
		      local_fringe = tgt_idx;
		    
		    if (trace > 1)
		      fprintf( trace_fh, "KILLED_BY_DNA\n" );
		  }
		}
		
		if (! invalid_pair) {
		  /* a legal pair of features, so calculate the score */
		  trans_score = 0.0;
		  len_pen = 0.0;
		  
		  seg_score = calculate_segment_score( src, tgt, segments, gs, seg_res );
		  trans_score += seg_score;
		  
		  if (reg_info->len_fun != NULL) {
		    Length_Function *lf = 
		      g_array_index(gs->length_funcs, Length_Function *, *(reg_info->len_fun));
		    len_pen = apply_Length_Function( lf, distance );
		  }
		  trans_score -= len_pen;
		  
		  backward_temp = tgt->backward_score
		    + trans_score
		    + tgt->score;
		  
		  /************ DEBUG: ************/
		  /* backward_temp = tgt->backward_score; */
		  /************ DEBUG: ************/

		  if (! touched_score || backward_temp > max_backward)
		    max_backward = backward_temp;
		  
		  g_array_append_val( all_scores, backward_temp);
		    
		  if (sum_mode == PRUNED_SUM) {
		      
		    if (! touched_score ) {
		      /* add back in the length penalty, because when judging for dominance, 
			 the length penalty will be different for future features */
		      max_backpluslen = backward_temp + len_pen;
		      local_fringe = tgt_idx;
		    }
		    else {
		      /* compare this one to min_backward, to see if it is dominated */
		      
		      if (backward_temp + len_pen > max_backpluslen) 
			max_backpluslen = backward_temp + len_pen;		  
		      
		      if (max_backpluslen - (backward_temp + len_pen) < 20.0)
			local_fringe = tgt_idx;
		      
		      /* this feature is not dominated. */
		    }
		  }
		  
		  touched_score = TRUE;
		  
		  if (trace > 1) 
		    fprintf( trace_fh, "Score: b=%.3f, (seg:%.3f len:%.3f)\n",
			     backward_temp, seg_score, len_pen );
		  
		} 
	      } /* if max dist */
	      else {
		/* This tgt will always be too distant for future targets of this type in the same
		   frame - therefore NO need to update fringe index
		*/
		if (trace > 1)
		  fprintf( trace_fh, "TOO DISTANT\n" );
	      }
	    } /* if min dist */
	    else {
	      /* This tgt may not be too close for future sources of this type in the same
		 frame - therefore need to update fringe index
	      */
	      if (sum_mode == PRUNED_SUM)
		local_fringe = tgt_idx;
	      
	      if (trace > 1)
		fprintf( trace_fh, "TOO CLOSE\n" );
	    }
	  } /* if mod3 */
	  else {
	    if (trace > 1)
	      fprintf( trace_fh, "PHASE-MISMATCH\n" );
	  }
	  
	} /* if match */
	else {
	  if (trace > 1)
	    fprintf( trace_fh, "TYPE-MISMATCH\n" );
	}
	
      } /* if valid */
      else {
	if (trace > 1)
	  fprintf( trace_fh, "INVALID\n" );
      }
      
      if ( src_info->kill_feat_quals_down != NULL) {
	if ((kq = g_array_index( src_info->kill_feat_quals_down, 
				 Killer_Feature_Qualifier *, 
				 tgt->feat_idx)) != NULL) {
	  if ((distance >= 0 && (kq->phase == (distance % 3))) || ! kq->has_phase)
	    passed_killer_feat = TRUE;
	}
      }
      
      if (tgt_info->is_killer_feat)
	g_array_append_val( poss_killer_feats, tgt );
      
    } /* for each target */         

    /* 
       update the position of the last forced feature. However, this
       needs to be the position of the first feature in the last
       forced 'block' (where 5'0, 5'1, 5'2 form a force block, for
       example). This is to get around the bother of the fact that
       a single splice site in the source data ends up as three
       splice sites in the feature list
    */
    
    if (src->is_selected) {
      if (g_res->last_selected > features->len)
	g_res->last_selected = src_idx;
      else {
	Feature *last_feat = g_array_index( features, 
					    Feature *, 
					    g_res->last_selected );
	if (last_feat->real_pos.s != src->real_pos.s || last_feat->real_pos.e != src->real_pos.e) 
	  g_res->last_selected = src_idx;
      }
    }
    
    if (touched_score) {
      for (src_idx=0; src_idx < all_scores->len; src_idx++) {
	g_res->score += exp( g_array_index(all_scores, double, src_idx) 
				- max_backward );  
      }
      g_res->score = log( g_res->score ) + max_backward;
      
      /* if we passed killer, then passes for future instances of this target type
	 in the same absolute frame will be terminated by the same killer. Therefore,
	 we do not have to update the fringe index */
      if (sum_mode == PRUNED_SUM && local_fringe > 0 && ! passed_killer_feat)
	g_array_index( g_res->fringe_idx, int *, src->feat_idx )[src->real_pos.s % 3] = local_fringe;
      
      if (trace)
	fprintf(trace_fh, "  RESULT: b=%.3f\n", 
		g_res->score );
    }
    else {
      if (trace) 
	fprintf( trace_fh, "  *** Invalidating\n");
      
      src->invalid = TRUE;
    }

  }
  else
    if (trace)
      fprintf( trace_fh, "  *** Invalid\n" );

  /* We don't need to set the score to negative infinity for the 
     calculation of Bbegin - however we need to do it for the sake
     of individual feature posterior probabilities */
  
  if (src->invalid)
    g_res->score = NEG_INFINITY;

  free_Seg_Results( seg_res );
  g_array_free( all_scores, TRUE );
  g_array_free( poss_killer_feats, TRUE );
   
}      



/*********************************************************************
 FUNCTION: trace_back_general
 DESCRIPTION:
   This function performs the dp tracback without using traceback
   pointers. It assumes however that the calculation of the 
   dp cumulative measures has already been calculated (via 
   'forward_calc')
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
GArray *trace_back_general ( GArray *feats,
			     GArray *segs,
			     Gaze_Structure *gs,
			     enum DP_Traceback_Mode tracemode) {
  
  int i;
  Feature *temp;

  GArray *stack = g_array_new( FALSE, TRUE, sizeof(Feature *));
  GArray *feat_path = g_array_new( FALSE, TRUE, sizeof(Feature *));
  Gaze_Results *g_res = new_Gaze_Results( gs->feat_dict->len );
  int pos = feats->len - 1;  

  temp = g_array_index( feats, Feature *, pos );
  g_array_append_val( stack, temp );


  while (pos > 0) {
#ifndef LOW_MEM
    if (tracemode != SAMPLE_TRACEBACK)  
      pos = temp->trace_pointer;    
    else {
      scan_through_sources_dp( feats, segs, pos,  g_res, gs, NO_SUM, tracemode, FALSE, NULL );
      pos = g_res->pth_trace;
    }
#else
    scan_through_sources_dp( feats, segs, pos,  g_res, gs, NO_SUM, tracemode, FALSE, NULL );
    pos = g_res->pth_trace;
#endif

    temp = g_array_index( feats, Feature *, pos );
    g_array_append_val( stack, temp );
  }


  if (pos == 0) {
    for (i=stack->len-1; i>=0; i--) {
      temp =  g_array_index( stack, Feature *, i );
      g_array_append_val( feat_path, temp);
    }
    /* for standard tracebacks, the path score in END will already be correct.
       However, we need to recalcuate the path score for sampled trances. Might 
       as well do it anyway, since it's cheap */
    calculate_path_score( feat_path, segs, gs );
  }
  else {
    g_array_free( feat_path, TRUE );
    feat_path = NULL;
  }

  free_Gaze_Results( g_res );
  g_array_free( stack, TRUE );

  return feat_path;
}
