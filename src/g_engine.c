/*  Last edited: Jul 13 14:02 2002 (klh) */
/**********************************************************************
 ** File: engine.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#include "g_engine.h"


/*********************************************************************
 FUNCTION: free_Gaze_DP_struct
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Gaze_DP_struct( Gaze_DP_struct *g_res, int feat_types ) {
  int i,j;

  if (g_res != NULL) { 

    /* The lists of features that were kept during the dp */
    if (g_res->feats != NULL) {
      for (i=0; i < feat_types; i++) {
	if ( g_res->feats[i] != NULL) {
	  for (j=0; j < 3; j++)
	    g_array_free( g_res->feats[i][j], TRUE );
	  free_util( g_res->feats[i] );
	}
      }

      free_util( g_res->feats);
    }

    /* The fringe indices that were kept during the dp */
    if (g_res->fringes != NULL) {
      for (i=0; i < feat_types; i++) {
	for (j=0; j < feat_types; j++)
	  free_util( g_res->fringes[i][j] );		
	free_util( g_res->fringes[i] );
      }
      free_util( g_res->fringes );
    }

    free_util( g_res );
  }  
}


/*********************************************************************
 FUNCTION: new_Gaze_DP_struct
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Gaze_DP_struct *new_Gaze_DP_struct( int feat_dict_size, int fringe_init ) {
  Gaze_DP_struct *g_res;
  int i, j, k;

  g_res = (Gaze_DP_struct *) malloc_util (sizeof(Gaze_DP_struct));

  g_res->pth_score = g_res->score = 0.0;
  g_res->last_selected = -1;
 
  /* The lists of features that will be kept during the dp */
  g_res->feats = (GArray ***) malloc_util( feat_dict_size * sizeof( GArray ** ));
  for(i=0; i < feat_dict_size; i++) {
    g_res->feats[i] = (GArray **) malloc_util( 3 * sizeof( GArray * ) );
    for(j=0; j < 3; j++) 
      g_res->feats[i][j] = g_array_new( FALSE, TRUE, sizeof(int) );
  }

  /* The indices of the fringes */
  g_res->fringes = (int ***) malloc_util( feat_dict_size * sizeof( int **) );
  for (i=0; i < feat_dict_size; i++) {
    g_res->fringes[i] = (int **) malloc_util( feat_dict_size * sizeof( int *) );
    for( j=0; j < feat_dict_size; j++) {
      g_res->fringes[i][j] = (int *) malloc_util( 3 * sizeof( int ) );
      for( k=0; k < 3; k++ )
	g_res->fringes[i][j][k] = fringe_init;
    }
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
double calculate_path_score(GArray *path, 
			    GArray *segs,
			    Gaze_Structure *gs) {
  
  int idx; 
  Feature_Relation *reg_info;
  Feature *src, *tgt;
  double total_score, trans_score, len_pen, seg_score;
  Seg_Results *s_res = new_Seg_Results( gs->seg_dict->len );

  total_score = 0.0;

  for( idx=0; idx < path->len - 1; idx++) { 
    /* for the score to mean anything, all paths must begin with "BEGIN"
       and end with "END." Therefore ignoring the local score of the first 
       feature (done here) has no effect, since the score of "BEGIN" is 0
       (and has to be for the DP to work */

    src = g_array_index( path, Feature *, idx );
    tgt = g_array_index( path, Feature *, idx + 1 );
    reg_info = g_array_index( g_array_index( gs->feat_info, Feature_Info *, tgt->feat_idx )->sources, 
			      Feature_Relation *, 
			      src->feat_idx);

    trans_score = len_pen = 0.0;
    
    seg_score = calculate_segment_score( src, tgt, segs, gs, s_res );
    trans_score += seg_score;
    
    if (reg_info->len_fun != NULL) {
      Length_Function *lf = 
	g_array_index(gs->length_funcs, Length_Function *, *(reg_info->len_fun));
      len_pen = apply_Length_Function( lf,  tgt->adj_pos.e - src->adj_pos.s + 1 );
    }
    trans_score -= len_pen;
    
    tgt->path_score = src->path_score + trans_score + tgt->score;
    total_score += trans_score + tgt->score;
  }

  free_Seg_Results( s_res );

  return total_score;
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
		    Gaze_Output *g_out) {
  
  int ft_idx, prev_idx;
  GArray *temp;
  Feature *prev_feat;

  Gaze_DP_struct *g_res = new_Gaze_DP_struct( gs->feat_dict->len, 0 );

#ifdef TRACE
    fprintf(stderr, "\nForward calculation:\n\n");
#endif

  for (ft_idx = 1; ft_idx < features->len; ft_idx++) {
    /* push the index of the last feature onto the list of
       sorted indices */
    prev_idx = ft_idx - 1;
    prev_feat = g_array_index( features, Feature *, prev_idx );
    temp = g_res->feats[prev_feat->feat_idx][prev_feat->adj_pos.s % 3];
    g_array_append_val( temp, prev_idx );

    scan_through_sources_dp( features, 
			     segments, 
			     gs, 
			     ft_idx,  
			     g_res, 
			     sum_mode,
			     MAX_TRACEBACK,
			     g_out);

    g_array_index( features, Feature *, ft_idx )->forward_score = g_res->score;
    g_array_index( features, Feature *, ft_idx )->path_score = g_res->pth_score;
    g_array_index( features, Feature *, ft_idx )->trace_pointer = g_res->pth_trace;
  }

  free_Gaze_DP_struct( g_res, gs->feat_dict->len );
}



/*********************************************************************
 FUNCTION: backwards_calc
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void backwards_calc( GArray *features,
		     GArray *segments,
		     Gaze_Structure *gs,
		     enum DP_Calc_Mode sum_mode) {

  int ft_idx, prev_idx;
  Feature *prev_feat;
  GArray *temp;

  Gaze_DP_struct *g_res = new_Gaze_DP_struct( gs->feat_dict->len, features->len - 1 );

  g_res->last_selected = features->len + 1;

#ifdef TRACE
  fprintf(stderr, "\nBackward calculation:\n\n");
#endif
  
  for (ft_idx = features->len-2; ft_idx >= 0; ft_idx--) {
    /* push the index if the last feature onto the list of
       sorted indices */
    prev_idx = ft_idx + 1;
    prev_feat = g_array_index( features, Feature *, prev_idx );
    temp = g_res->feats[prev_feat->feat_idx][prev_feat->adj_pos.e % 3];
    g_array_append_val( temp, prev_idx );

    scan_through_targets_dp( features, 
			     segments, 
			     gs,
			     ft_idx,  
			     g_res,
			     sum_mode);

    g_array_index( features, Feature *, ft_idx )->backward_score = g_res->score;

  }

  free_Gaze_DP_struct( g_res,  gs->feat_dict->len );
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
			     Gaze_Structure *gs,
			     int tgt_idx,
			     Gaze_DP_struct *g_res,
			     enum DP_Calc_Mode sum_mode,
			     enum DP_Traceback_Mode trace_mode,
			     Gaze_Output *g_out) {
  
  int src_type, src_idx, kill_idx, max_index = 0; /* Initialsied to get arounc gcc warnings */
  int frame, k, index_count[3];
  int last_necessary_idx, local_fringe, last_idx_for_frame[3];
  int left_pos, right_pos, distance;
  Killer_Feature_Qualifier *kq;

  gboolean touched_score, touched_score_local;
  GArray *all_scores = NULL;
  GArray *all_indices = NULL;  /* Initialsied to get around gcc warnings */
  Feature_Info *tgt_info;
  Feature_Relation *reg_info;
  Feature *src, *tgt;
  Seg_Results *seg_res;

  double max_forpluslen = 0.0;
  double max_score = NEG_INFINITY;
  double max_forward = NEG_INFINITY;
  int *killer_source_dna = NULL;
  int *danger_source_dna = NULL;

  g_res->pth_score = 0.0;
  g_res->pth_trace = 0,0;
  g_res->score = 0.0;

  tgt = g_array_index( features, Feature *, tgt_idx );
  tgt_info = g_array_index( gs->feat_info, Feature_Info *, tgt->feat_idx );
  right_pos = tgt->adj_pos.e;

#ifdef TRACE 
  fprintf( stderr, "Target %d %s %d %d %.3f", tgt_idx,
	   g_array_index(gs->feat_dict, char *, tgt->feat_idx ), 
	   tgt->real_pos.s, tgt->real_pos.e, tgt->score );
  
  if (TRACE > 1)
    fprintf( stderr, "\n" );
#endif
  seg_res = new_Seg_Results( gs->seg_dict->len );

  if (sum_mode == STANDARD_SUM || sum_mode == PRUNED_SUM || trace_mode == SAMPLE_TRACEBACK) {
    all_scores = g_array_new( FALSE, TRUE, sizeof(double) );
    all_indices = g_array_new( FALSE, TRUE, sizeof(int) );
  }

  if (! tgt->invalid) {
    touched_score = FALSE;
    
    /* set up the boundaries for the scan. We do not want to go past:
       1. The last forced feature */

    last_necessary_idx = 0; 

    if (g_res->last_selected > last_necessary_idx)
      last_necessary_idx = g_res->last_selected;
    
    /* Look through the sources themselves */
    
    for ( src_type = 0; src_type < tgt_info->sources->len; src_type++) {
      if ((reg_info = g_array_index(tgt_info->sources, Feature_Relation *, src_type)) != NULL) {
	GArray **feats = g_res->feats[src_type];
	
	for(frame = 0; frame < 3; frame++) {
	  
	  last_idx_for_frame[frame] = last_necessary_idx;
	  
	  /* first, identify the killers local to this feature pair, and make
	     sure that our search back through the sources does not go past a killer */
	  
	  if (reg_info->kill_feat_quals != NULL) {
	    
	    for(kill_idx=0; kill_idx < reg_info->kill_feat_quals->len; kill_idx++) {
	      if ( (kq = g_array_index( reg_info->kill_feat_quals,
					Killer_Feature_Qualifier *,
					kill_idx )) != NULL) {
		/* The following checks to see if last killer in the correct
		   frame has an index greater than the source; if not, none
		   of them can have. However, this does not allow for weird 
		   edge effects that occur when killers overlap with other
		   features. This may or may not make a difference with the
		   new feature ordering */
		
		if (kq->has_tgt_phase) {
		  GArray *apt_list = (GArray *) g_res->feats[kq->feat_idx][(right_pos - kq->phase + 1) % 3];
		  /* rationale: (right_pos - left_pos + 1) % 3 == phase -->
		     (right_pos - {left_pos % 3} + 1) % 3 == phase -->
		     (right_pos - {left_pos % 3} + 1) % 3 - phase == 0 -->
		     (right_pos - {left_pos % 3} + 1 - phase) % 3 == 0 -->
		     (right_pos - phase + 1) % 3 - {left_pos % 3} == 0 -->
		     (right_pos - phase + 1) % 3 == {left_pos % 3} */
		  if (apt_list->len > 0 
		      && g_array_index( apt_list, int, apt_list->len - 1) > last_idx_for_frame[frame] )
		    last_idx_for_frame[frame] = g_array_index( apt_list, int, apt_list->len - 1);
		}		  
		else if (kq->has_src_phase) {
		    /* the frame calculation here is slightly hacky; we have to allow for the
		       fact that we need the distance from the source forward to the apt. killer.
		       BUT the the killers are stored by the frame of their adjusted START, 
		       rather than their end. So, we rely on the fact that all killers that have
		       a phase are width 3, which might not be so unreasonable */
		  GArray *apt_list = (GArray *) g_res->feats[kq->feat_idx][(frame + kq->phase) % 3];
		  if (apt_list->len > 0 
		      && g_array_index( apt_list, int, apt_list->len - 1) > last_idx_for_frame[frame] )
		    last_idx_for_frame[frame] = g_array_index( apt_list, int, apt_list->len - 1);
		}
		else {
		  /* phaseless killer - need to check all frames */
		  for (k=0; k < 3; k++) {
		    GArray *apt_list = (GArray *) g_res->feats[kq->feat_idx][k];
		    if (apt_list->len > 0 
			&& g_array_index( apt_list, int, apt_list->len - 1) > last_idx_for_frame[frame]) {
		      last_idx_for_frame[frame] = g_array_index( apt_list, int, apt_list->len - 1);
		    }			
		  }
		}

	      }
	    }
	  }
	  
	  /* finally, make sure that we do not proceed past the fringe for
	     this feature pair */
	  if (g_res->fringes[tgt->feat_idx][src_type][tgt->real_pos.s % 3] > last_idx_for_frame[frame])
	    last_idx_for_frame[frame] = g_res->fringes[tgt->feat_idx][src_type][tgt->real_pos.s % 3];
	}
	
	/* Before actually scanning through the features themselves, we need to check
	   if there are potential dna killers. If so, set flags for the source
	   dna entries that will cause problems */
	
	if (reg_info->kill_dna_quals != NULL) {
	  danger_source_dna = (int *) malloc0_util( gs->motif_dict->len * sizeof(int) );
	  
	  for(k=0; k < reg_info->kill_dna_quals->len; k++) {
	    Killer_DNA_Qualifier *kdq = g_array_index( reg_info->kill_dna_quals, 
						       Killer_DNA_Qualifier *,
						       k );
	    
	    danger_source_dna[kdq->src_dna] = 1; 
	    
	    if (tgt->dna >= 0 && tgt->dna == kdq->tgt_dna) {
	      if (killer_source_dna == NULL) 
		killer_source_dna = (int *) malloc0_util( gs->motif_dict->len * sizeof(int) );
	      killer_source_dna[kdq->src_dna] = 1;
	    }
	  }
	}
	
	/* at this point, we have the list of features that need to be processed (feats),
	   and the index that we must not proceed past in each frame. We can now process
	   the features themselves, in a frame-dependent or frame-independent way */
	
#ifdef TRACE
	if (TRACE > 1)
	  fprintf( stderr, "  %s (fringes: %d %d %d)\n",
		   g_array_index(gs->feat_dict, char *, src_type ), 
		   last_idx_for_frame[0], last_idx_for_frame[1], last_idx_for_frame[2] );
#endif

	for(k=0; k < 3; k++) 
	  index_count[k] = feats[k]->len - 1; 
	
	frame = reg_info->phase != NULL ? (right_pos - *(reg_info->phase) + 1) % 3 : 0;
	
	max_forpluslen = NEG_INFINITY;
	touched_score_local = FALSE;
	/* the following aggressively assumes that if this target has no 
	   potential sources for this source type, then we need go no 
	   further back than the index of the target itself when consdering
	   furture instances of the target */
	local_fringe = tgt_idx;

	while (1) {
	  /* This loop is terminated by break statements scattered throughout. 
	     I know, I know... */
	  
	  if (reg_info->phase == NULL) {
	    /* For frameless feature pairs, we need to examine all frames, but 
	       for the pruning to work effectively, the featured need to be examined 
	       in order. Therefore, we have to re-create the original list
	       by effectively a merge sort. */
	    gboolean gotone = FALSE;
	    
	    for (k=0; k < 3; k++) {
	      if (index_count[k] >= 0) {
		if (gotone) {
		  if (g_array_index( feats[k], int, index_count[k] ) > 
		      g_array_index( feats[frame], int, index_count[frame] ) ) {
		    frame = k;
		  }
		}
		else {
		  frame = k;
		  gotone = TRUE;
		}
	      }
	    }
	  }
	  
	  /* The following tests if there was anything in the list at all */
	  if (index_count[frame] < 0)
	    break;
	  
	  src_idx = g_array_index( feats[frame], int, index_count[frame]-- );
	  
	  if (src_idx < last_idx_for_frame[frame]) {
	    /* we must be careful not to simply break out of the loop at this 
	       point, because for phaseless sources we are flipping between frames,
	       so there may be others sources in different frames still to consider.
	       However, we do know that we need consider no more features
	       if this type in THIS frame, which can be achieved by the 
	       following trick: */
	    index_count[frame] = -1;
	    continue;
	  }
	  
	  src = g_array_index( features, Feature *, src_idx );
	  
#ifdef TRACE
	  if (TRACE > 1)
	    fprintf( stderr, "     Source %d %s %d %d ", src_idx,
		     g_array_index(gs->feat_dict, char *, src_type ),
		     src->real_pos.s, src->real_pos.e );
#endif
	  
	  if (! src->invalid) {
	    
	    left_pos = src->adj_pos.s;
	    distance = right_pos - left_pos + 1;
	    
#ifdef TRACE
	    if (TRACE > 1)
	      fprintf( stderr, "dist=%d  ", distance );
#endif	    

	    if ((reg_info->max_dist == NULL) || (*(reg_info->max_dist)) >= distance) {
	      
	      if ((reg_info->min_dist == NULL) || (*(reg_info->min_dist)) <= distance) {
		/* Finally, if this source does not result in a DNA kill, we can calc the score */
		if (killer_source_dna == NULL || src->dna < 0 || ! killer_source_dna[src->dna]) {
		  double trans_score, len_pen, seg_score, forward_temp, viterbi_temp;
		  Length_Function *lf = NULL;
		  trans_score = len_pen = forward_temp = viterbi_temp = 0.0;

		  seg_score = calculate_segment_score( src, tgt, segments, gs, seg_res );
		  trans_score += seg_score;
		  
		  if (reg_info->len_fun != NULL) {
		    lf = g_array_index(gs->length_funcs, Length_Function *, *(reg_info->len_fun));
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
		    
		    g_array_append_val( all_scores, forward_temp);
		    g_array_append_val( all_indices, src_idx);
		  
		    if (! touched_score || (forward_temp > max_forward))
		      max_forward = forward_temp;
		    
		    if (sum_mode == PRUNED_SUM) {
		      /* There are two assumptions for my pruning method:
			 1. Because the scores are log scores, if two scores differ
			    by 25 (say) or more, then the first score is e^25 times bigger
			    than the other; the smaller score will not register given 
			    machine precision, so can be ignored. 
			 2. If all features to the left of a given source are "dominated" in 
			    this way, they will be dominated for all subsequence occurrences
			    of the current target, so can be pruned away
			 
			    However, assumption 2 does not quite hold when the dominant source
			    for a given target is not valid for a future target of this type,
			    due to DNA killers. Therefore, we ensure that sources
			    do not dominate if they might be illegal with respect to future
			    targets of this type */
		      
		      if (! touched_score_local ) {
			
			if (danger_source_dna == NULL 
			    || src->dna < 0  
			    || ! danger_source_dna[src->dna]) {

			  /* strictly speaking, it is only sound to register this source
			     as "dominant" if we are into the monotonic part of the length
			     function (i.e. the point past which the penalty never decreases 
			     with increasing distance). Therefore, we update the fringe, but
			     don't flag touched_local_score. The efect of this is that the
			     fringe will be updated for all scoring sources until we get past
			     the point of monotonicity, at which point the pruning kicks in */

			  if (lf == NULL || (lf->becomes_monotonic && lf->monotonic_point <= distance)) { 
			    /* add back in the length penalty, because when judging for dominance, 
			       the length penalty will be different for future features */
			    
			    max_forpluslen = forward_temp + len_pen;
			    touched_score_local = TRUE;
			  } 
			}
			
			local_fringe = src_idx;
		      }
		      else {
			/* compare this one to max_forward, to see if it is dominated */
			if (forward_temp + len_pen > max_forpluslen 
			    && (danger_source_dna == NULL 
				|| src->dna < 0  
				|| ! danger_source_dna[src->dna]))
			  max_forpluslen = forward_temp + len_pen;
			
			if ( max_forpluslen - (forward_temp + len_pen) < 25.0)
			  local_fringe = src_idx;
		      }
		    }		    
		  }
		  
		  touched_score = TRUE;
		  
#ifdef TRACE
		  if (TRACE > 1) 
		    fprintf( stderr, "scre: v=%.3f, f=%.8f (seg:%.3f len:%.3f)\n",
			     viterbi_temp, forward_temp, seg_score, len_pen );
#endif
		  if (g_out != NULL) {
		    /* the following should be changed to a structure file lookup */
		    if (reg_info->out_qual != NULL && reg_info->out_qual->need_to_print) {
		      
		      double reg_score = trans_score;
		      
		      if (g_out->posterior) 
			reg_score = exp( src->forward_score + 
					 trans_score + tgt->score +				   
					 tgt->backward_score - 
					 g_array_index( features, Feature *, 0)->backward_score );
		      
		      if (! g_out->use_threshold || reg_score >= g_out->threshold)
			fprintf(g_out->fh, "%s\tGAZE\t%s\t%d\t%d\t%.5f\t%s\t%s\t\n",
				g_out->seq_name, 
				reg_info->out_qual->feature != NULL ? reg_info->out_qual->feature : "Anonymous",
				left_pos, 
				right_pos, 
				reg_score,
				reg_info->out_qual->strand != NULL ? reg_info->out_qual->strand : ".", 
				reg_info->out_qual->frame != NULL ?  reg_info->out_qual->frame : ".");
		      
		    }
		  }
		} /* if killed by DNA */
		else {
		  /* source might not be killed for future incidences, so update fringe index */
		  if (sum_mode == PRUNED_SUM)
		    local_fringe = src_idx;
		  
#ifdef TRACE
		  if (TRACE > 1)
		    fprintf( stderr, "KILLED_BY_DNA\n" );
#endif
		}
	      } /* if min dist */
	      else {
		/* source might not be too close for future incidences, so update fringe index */
		if (sum_mode == PRUNED_SUM)
		  local_fringe = src_idx;
		
#ifdef TRACE
		if (TRACE > 1)
		  fprintf( stderr, "TOO CLOSE\n" );
#endif
	      }
	    } /* if max dist */
	    else {
#ifdef TRACE
	      if (TRACE > 1)
		fprintf( stderr, "TOO DISTANT\n" );
#endif
	      /* we can break out of the loop here; all other sources will be too distant */
	      break;
	    }
	  }
#ifdef TRACE
	  else 
	    if (TRACE > 1)
	      fprintf( stderr, "INVALID\n" );
#endif
	}
      
	if (sum_mode == PRUNED_SUM) {
	  /* We conservatively only prune in the frame of the target if this
	     feature pair has a phase constraint, or if there are potential
	     killers (which also might have a phase constraint). Otherwise, 
	     prune in all frames */
	  if (reg_info->phase != NULL || reg_info->kill_feat_quals != NULL) {
	    g_res->fringes[tgt->feat_idx][src_type][tgt->real_pos.s % 3] = local_fringe;
	  }
	  else {
	    for(k=0; k < 3; k++)
	      g_res->fringes[tgt->feat_idx][src_type][k] = local_fringe;
	  }
	}
	
	if (danger_source_dna != NULL) {
	  free_util( danger_source_dna );
	  danger_source_dna = NULL;
	}
	if (killer_source_dna != NULL) {
	  free_util( killer_source_dna );
	  killer_source_dna = NULL;
	}
      }  
    }
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
      
      if (sum_mode == STANDARD_SUM || sum_mode == PRUNED_SUM) {
	/* the trick of subtracting the max before exponentiating avoids
	   overflow errors. Just need to add it back when logging back down */
	for (src_idx=0; src_idx < all_scores->len; src_idx++) {
	  g_res->score += exp( g_array_index(all_scores, double, src_idx) 
			       - max_forward );  
	  
	}
	g_res->score = log( g_res->score ) + max_forward;	
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
      else if (trace_mode == MAX_TRACEBACK) {
	g_res->pth_trace = max_index;
	g_res->pth_score = max_score;	
      }
      
#ifdef TRACE
      fprintf(stderr, "  RESULT: v=%.3f, max=%d, f=%.8f\n", 
	      g_res->pth_score,
	      max_index,
	      g_res->score );
#endif
    }
    else {
      tgt->invalid = TRUE;

#ifdef TRACE
      fprintf( stderr, "  *** Invalidating\n");
#endif      

    }
  }
#ifdef TRACE
  else
    fprintf( stderr, "  *** Invalid\n" );
#endif  

  /* We don't need to set the score to negative infinity for the 
     calculation of Fend - however we need to do it for the sake
     of individual feature posterior probabilities */
  
  if (tgt->invalid)
    g_res->score = NEG_INFINITY;

  if (sum_mode == STANDARD_SUM || sum_mode == PRUNED_SUM ||  trace_mode == SAMPLE_TRACEBACK) {
    g_array_free( all_indices, TRUE );
    g_array_free( all_scores, TRUE );
  }
  free_Seg_Results( seg_res );
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
			     Gaze_Structure *gs,
			     int src_idx,
			     Gaze_DP_struct *g_res,
			     enum DP_Calc_Mode sum_mode) {

  int tgt_type, tgt_idx, kill_idx;
  int frame, k, index_count[3]; 
  int left_pos, right_pos, distance;
  int last_necessary_idx, local_fringe, last_idx_for_frame[3];
  Killer_Feature_Qualifier *kq;

  gboolean touched_score, touched_score_local;
  GArray *all_scores;
  Feature_Info *src_info, *tgt_info;
  Feature_Relation *reg_info;
  Feature *src, *tgt;
  Seg_Results *seg_res;

  double max_backward = 0.0;
  double max_backpluslen = 0.0;
  int *killer_target_dna = NULL;
  int *danger_target_dna = NULL;

  g_res->score = 0.0;

  src = g_array_index( features, Feature *, src_idx );
  src_info = g_array_index( gs->feat_info, Feature_Info *, src->feat_idx );
  left_pos = src->adj_pos.s;

#ifdef TRACE
    fprintf( stderr, "Source %d %s %d %d %.3f", src_idx,
	     g_array_index(gs->feat_dict, char *, src->feat_idx ), 
	     src->real_pos.s, src->real_pos.e, src->score ); 
    
    if (TRACE > 1)
      fprintf( stderr, "\n" );
#endif
  
  seg_res = new_Seg_Results( gs->seg_dict->len );
  all_scores = g_array_new( FALSE, TRUE, sizeof(double));

  if (! src->invalid ) { 
    touched_score = FALSE;

    /* set up the boundaries for the scan. We do not want to go past 
       1. The last forced feature */

    last_necessary_idx = features->len - 1; 

    if (g_res->last_selected < last_necessary_idx)
      last_necessary_idx = g_res->last_selected;
    
    /* Look through the targets themselves */
    
    for ( tgt_type = 0; tgt_type < gs->feat_dict->len; tgt_type++) {
      tgt_info = g_array_index( gs->feat_info, Feature_Info *, tgt_type );

      /* this guards against looking for sources of BEGIN, of which there are none */
      if (tgt_info->sources == NULL)
	continue;

      if ((reg_info = g_array_index(tgt_info->sources, Feature_Relation *, src->feat_idx)) != NULL) {
	GArray **feats = g_res->feats[tgt_type];

	for(frame = 0; frame < 3; frame ++) {
	  
	  last_idx_for_frame[frame] = last_necessary_idx;

	  if (reg_info->kill_feat_quals != NULL) {
	    /* the following can be speeded up quite easiy, now that we
	       have full indexing into feature types. Requires the killers
	       to be stored as a simple list in the structure. Ill do it 
	       the slow way for now though */

	    for(kill_idx=0; kill_idx < reg_info->kill_feat_quals->len; kill_idx++) {
	      if ( (kq = g_array_index( reg_info->kill_feat_quals,
					Killer_Feature_Qualifier *,
					kill_idx )) != NULL) {
		/* The following checks to see if last killer in the correct
		   frame has an index less than the target; if not, none
		   of them can have. However, this does not allow for weird 
		   edge effects that occur when killers overlap with other
		   features. This may or may not make a difference with the
		   new feature ordering */
		
		if (kq->has_src_phase) {
		  GArray *apt_list = (GArray *) g_res->feats[kq->feat_idx][(left_pos + kq->phase - 1) % 3];
		  /* rationale: (right_pos - left_pos + 1) % 3 == phase -->
	                  (right_pos - left_pos + 1) % 3 - phase == 0 -->
			  (right_pos - left_pos + 1 - phase) % 3 == 0 -->
			  (left_pos - right_pos - 1 + phase) % 3 == 0 -->
			  (left_pos - {right_pos % 3} + phase - 1) % 3 == 0 -->
			  (left_pos + phase - 1) % 3 - {right_pos % 3} == 0 -->
			  (left_pos + phase - 1) % 3 == {right_pos % 3} */
		  if (apt_list->len > 0 
		      && g_array_index( apt_list, int, apt_list->len - 1) < last_idx_for_frame[frame] ) {
		    last_idx_for_frame[frame] = g_array_index( apt_list, int, apt_list->len - 1);
		  }
		}
		if (kq->has_tgt_phase) {
		  /* the frame calcualtion here is slightly hacky; we have to allow for the
		     fact that we need the distance from the target back to the apt. killer.
		     BUT, the the killers are stored by the frame of their adjusted END, 
		     rather than their start. So, we rely on the fact that all killers that have 
		     a phase are width 3, which might not be so unreasonable */
		  GArray *apt_list = (GArray *) g_res->feats[kq->feat_idx][(frame + 3 - kq->phase) % 3];
		  if (apt_list->len > 0 
		      && g_array_index( apt_list, int, apt_list->len - 1) < last_idx_for_frame[frame] ) {
		    last_idx_for_frame[frame] = g_array_index( apt_list, int, apt_list->len - 1);
		  }
		}	
		else {
		  /* phaseless killer - need to check all frames */
		  for (k=0; k < 3; k++) {
		    GArray *apt_list = (GArray *) g_res->feats[kill_idx][k];
		    if (apt_list->len > 0 
			&& g_array_index( apt_list, int, apt_list->len - 1) < last_idx_for_frame[frame]) {
		      last_idx_for_frame[frame] = g_array_index( apt_list, int, apt_list->len - 1);
		    }			
		  }
		}
	      }
	    }
	  }

	  /* finally, make sure that we do not proceed past the fringe for
	     this feature pair */
	  if (g_res->fringes[src->feat_idx][tgt_type][src->real_pos.s % 3] < last_idx_for_frame[frame])
	    last_idx_for_frame[frame] = g_res->fringes[src->feat_idx][tgt_type][src->real_pos.s % 3];
	}


	/* Before actually scanning through the features themselves, we need to check
	   if there are potential dna killers. If so, set flags for the source
	   dna entries that will cause problems */
	
	if (reg_info->kill_dna_quals != NULL) {
	  danger_target_dna = (int *) malloc0_util( gs->motif_dict->len * sizeof(int) );
	  
	  for(k=0; k < reg_info->kill_dna_quals->len; k++) {
	    Killer_DNA_Qualifier *kdq = g_array_index( reg_info->kill_dna_quals, 
						       Killer_DNA_Qualifier *,
						       k );
	    
	    danger_target_dna[kdq->tgt_dna] = 1;
	    
	    if (src->dna >= 0 && src->dna == kdq->src_dna) {
	      if (killer_target_dna == NULL)
		killer_target_dna = (int *) malloc0_util( gs->motif_dict->len * sizeof(int) );
	      killer_target_dna[kdq->tgt_dna] = 1;
	    }
	  }
	}
	
	/* at this point, we have the list of features that need to be processed (feats),
	   and the index that we must not proceed past in each frame. We can now process
	   the features themselves, in a frame-dependent or frame-independent way */

#ifdef TRACE
	if (TRACE > 1)
	  fprintf( stderr, "  %s (fringes: %d %d %d)\n",
		   g_array_index(gs->feat_dict, char *, tgt_type ), 
		   last_idx_for_frame[0], last_idx_for_frame[1], last_idx_for_frame[2] );
#endif

	for(k=0; k < 3; k++) 
	  index_count[k] = feats[k]->len - 1; 
	
	frame = reg_info->phase != NULL ? (left_pos + *(reg_info->phase) - 1) % 3 : 0;
	
	max_backpluslen = NEG_INFINITY;
	touched_score_local = FALSE;
	/* the following aggressively assumes that if this target has no 
	   potential sources for this source type, then we need go no 
	   further back than the index of the target itself when consdering
	   furture instances of the target */
	local_fringe = src_idx;

	while(1) {
	  /* This loop is terminated by break statements scattered throughout. 
	     I know, I know... */
	  
	  if (reg_info->phase == NULL) {
	    /* For frameless feature pairs, we need to examine all frames, but 
	       for the pruning to work effectively, the featured need to be examined 
	       in order. Therefore, we have to re-create the original list
	       by effectively a merge sort. */
	    gboolean gotone = FALSE;
	    
	    for (k=0; k < 3; k++) {
	      if (index_count[k] >= 0) {
		if (gotone) {
		  if (g_array_index( feats[k], int, index_count[k] ) < 
		      g_array_index( feats[frame], int, index_count[frame] ) ) {
		    frame = k;
		  }
		}
		else {
		  frame = k;
		  gotone = TRUE;
		}
	      }
	    }
	  }
	  
	  /* The following tests if there was anything in the list at all */
	  if (index_count[frame] < 0)
	    break;
	  
	  tgt_idx = g_array_index( feats[frame], int, index_count[frame]-- );

	  if (tgt_idx > last_idx_for_frame[frame]) {
	    /* we must be careful not to simply break out of the loop at this 
	       point, because for phaseless sources we are mixing the frame,
	       so there may be others in different frames still to consider.
	       However, we do know that we need consider no more features
	       if this type in this frame, which can be achieved by the 
	       following trick: */
	    index_count[frame] = -1;
	    continue;
	  }

	  tgt = g_array_index( features, Feature *, tgt_idx );
	    
#ifdef TRACE
	  if (TRACE > 1)
	    fprintf( stderr, "  Target %d %s %d %d  ", tgt_idx,
		     g_array_index(gs->feat_dict, char *, tgt->feat_idx ), 
		     tgt->real_pos.s, tgt->real_pos.e );
#endif	  
	  if (! tgt->invalid) {
	      
	    right_pos = tgt->adj_pos.e;
	    distance = right_pos - left_pos + 1;
	    
#ifdef TRACE
	    if (TRACE > 1)
	      fprintf( stderr, "dist=%d  ", distance );
#endif	    
	    if ((reg_info->max_dist == NULL) || (*(reg_info->max_dist)) >= distance) {
	      
	      if ((reg_info->min_dist == NULL) || (*(reg_info->min_dist)) <= distance) {

		/* Finally, if this source does not result in a DNA kill, we can calc the score */
		if (killer_target_dna == NULL || tgt->dna < 0 || ! killer_target_dna[tgt->dna]) {
		  double trans_score, len_pen, seg_score, backward_temp;
		  Length_Function *lf = NULL;
		  trans_score = len_pen = 0.0;
		  
		  seg_score = calculate_segment_score( src, tgt, segments, gs, seg_res );
		  trans_score += seg_score;
		  
		  if (reg_info->len_fun != NULL) {
		      lf = g_array_index(gs->length_funcs, Length_Function *, *(reg_info->len_fun));
		      len_pen = apply_Length_Function( lf, distance );
		  }
		  trans_score -= len_pen;
		  
		  backward_temp = tgt->backward_score
		    + trans_score
		    + tgt->score;
		  
		  if (! touched_score || backward_temp > max_backward)
		    max_backward = backward_temp;
		  
		  g_array_append_val( all_scores, backward_temp);

		  if (sum_mode == PRUNED_SUM) {
		    if (! touched_score_local ) {
		      /* add back in the length penalty, because when judging for dominance, 
			 the length penalty will be different for future features */

		      if (danger_target_dna == NULL 
			  || tgt->dna < 0  
			  || ! danger_target_dna[src->dna]) {

			/* strictly speaking, it is only sound to register this source
			   as "dominant" if we are into the monotonic part of the length
			   function (i.e. the point past which the penalty never decreases 
			   with increasing distance). Therefore, we update the fringe, but
			   don't flag touched_local_score. The efect of this is that the
			   fringe will be updated for all scoring sources until we get past
			   the point of monotonicity, at which point the pruning kicks in */

			if (lf == NULL || (lf->becomes_monotonic && lf->monotonic_point <= distance)) {
			  touched_score_local = TRUE;
			  max_backpluslen = backward_temp + len_pen;
			}
		      }

		      local_fringe = tgt_idx;

		    }
		    else {
		      /* compare this one to max_forward, to see if it is dominated */
		      if (backward_temp + len_pen > max_backpluslen
			  && (danger_target_dna == NULL 
			      || tgt->dna < 0  
			      || ! danger_target_dna[tgt->dna])) 
			max_backpluslen = backward_temp + len_pen;
		      
		      if ( max_backpluslen - (backward_temp + len_pen) < 25.0)
			local_fringe = tgt_idx;
		    }
		  }
		  
		  touched_score = TRUE;
		  
#ifdef TRACE
		  if (TRACE > 1) 
		    fprintf( stderr, "Score: b=%.3f, (seg:%.3f len:%.3f)\n",
			     backward_temp, seg_score, len_pen );
#endif
		  
		} /* if killed by DNA */
		else {
		  /* tgt might not be killed for future incidences, so update fringe index */
		  if (sum_mode == PRUNED_SUM)
		    local_fringe = tgt_idx;

#ifdef TRACE
		  if (TRACE > 1)
		    fprintf( stderr, "KILLED_BY_DNA\n" );
#endif
		}  
	      } /* if min dist */
	      else {
		/* target might not be too close for future incidences, so update fringe index */
		if (sum_mode == PRUNED_SUM)
		  local_fringe = tgt_idx;

#ifdef TRACE		
		if (TRACE > 1)
		  fprintf( stderr, "TOO CLOSE\n" );
#endif
	      }
	    } /* if max dist */
	    else {
#ifdef TRACE
	      if (TRACE > 1)
		fprintf( stderr, "TOO DISTANT\n" );
#endif
	      /* we can break out of the loop here; they will all be too distant */
	      break;
	    }
	  } /* if valid */
#ifdef TRACE
	  else {
	    if (TRACE > 1)
	      fprintf( stderr, "INVALID\n" );
	  }
#endif
	}

	if (sum_mode == PRUNED_SUM) {
	  /* We conservatively only prune in the frame of the target if this
	     feature pair has a phase constraint, or if there are potential
	     killers (which also might have a phase constraint). Otherwise, 
	     prune in all frames */
	  if (reg_info->phase != NULL || reg_info->kill_feat_quals != NULL) {
	    g_res->fringes[src->feat_idx][tgt_type][src->real_pos.s % 3] = local_fringe;
	  }
	  else {
	    for(k=0; k < 3; k++)
	      g_res->fringes[src->feat_idx][tgt_type][k] = local_fringe;
	  }
	}

	if (danger_target_dna != NULL) {
	  free_util( danger_target_dna );
	  danger_target_dna = NULL;
	}
	if (killer_target_dna != NULL) {
	  free_util( killer_target_dna );
	  killer_target_dna = NULL;
	}
      }      
    }

    /* update the position of the last forced feature. However, this
       needs to be the position of the first feature in the last
       forced 'block' (where 5'0, 5'1, 5'2 form a force block, for
       example). This is to get around the bother of the fact that
       a single splice site in the source data ends up as three
       splice sites in the feature list */
    
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
      for (tgt_idx=0; tgt_idx < all_scores->len; tgt_idx++) {
	g_res->score += exp( g_array_index(all_scores, double, tgt_idx) 
				- max_backward );  
      }
      g_res->score = log( g_res->score ) + max_backward;
      
#ifdef TRACE
      fprintf(stderr, "  RESULT: b=%.8f\n", 
		g_res->score );
#endif      
    }
    else {
      src->invalid = TRUE;

#ifdef TRACE
      fprintf( stderr, "  *** Invalidating\n");
#endif      
    }

  }
#ifdef TRACE
  else
    fprintf( stderr, "  *** Invalid\n" );
#endif

  /* We don't need to set the score to negative infinity for the 
     calculation of Bbegin - however we need to do it for the sake
     of individual feature posterior probabilities */
  
  if (src->invalid)
    g_res->score = NEG_INFINITY;

  free_Seg_Results( seg_res );
  g_array_free( all_scores, TRUE );
   
}      




/*********************************************************************
 FUNCTION: trace_back_general
 DESCRIPTION:
   This function performs the dp tracback. If a sample trace back
   is required, the traceback pointers are not used
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
GArray *trace_back_general ( GArray *feats,
			     GArray *segs,
			     Gaze_Structure *gs,
			     enum DP_Traceback_Mode tb_mode) {
  
  int i;
  Feature *temp;

  GArray *stack = g_array_new( FALSE, TRUE, sizeof(Feature *));
  GArray *feat_path = g_array_new( FALSE, TRUE, sizeof(Feature *));
  int pos = feats->len - 1;  

  temp = g_array_index( feats, Feature *, pos );
  g_array_append_val( stack, temp );

  if (tb_mode == SAMPLE_TRACEBACK) {
    Gaze_DP_struct *g_res = new_Gaze_DP_struct( gs->feat_dict->len, 0 );

    while (pos > 0) {
      scan_through_sources_dp( feats, segs, gs, pos, g_res,
			       NO_SUM, SAMPLE_TRACEBACK, NULL );
      pos = g_res->pth_trace;
      
      temp = g_array_index( feats, Feature *, pos );
      g_array_append_val( stack, temp );
    }
    
    free_Gaze_DP_struct( g_res, gs->feat_dict->len );
  }
  else {
    while (pos > 0) {
      pos = temp->trace_pointer;    
      
      temp = g_array_index( feats, Feature *, pos );
      g_array_append_val( stack, temp );
    }
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

  g_array_free( stack, TRUE );

  return feat_path;
}
