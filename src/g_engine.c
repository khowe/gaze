/*  Last edited: Aug  3 15:26 2002 (klh) */
/**********************************************************************
 ** File: engine.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#include "g_engine.h"
#include "time.h"

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
	    free_Array( g_res->feats[i][j], TRUE );
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


    if (g_res->seg_res != NULL)
      free_Seg_Results( g_res->seg_res );

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
Gaze_DP_struct *new_Gaze_DP_struct( int feat_dict_size, 
				    int seg_dict_size,
				    int fringe_init ) {
  Gaze_DP_struct *g_res;
  int i, j, k;

  g_res = (Gaze_DP_struct *) malloc_util (sizeof(Gaze_DP_struct));

  g_res->pth_score = g_res->score = 0.0;
  g_res->last_selected = -1;
 
  /* The lists of features that will be kept during the dp */
  g_res->feats = (Array ***) malloc_util( feat_dict_size * sizeof( Array ** ));
  for(i=0; i < feat_dict_size; i++) {
    g_res->feats[i] = (Array **) malloc_util( 3 * sizeof( Array * ) );
    for(j=0; j < 3; j++) 
      g_res->feats[i][j] = new_Array( sizeof(int), TRUE );
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

  g_res->seg_res = new_Seg_Results( seg_dict_size );

  return g_res;
}



/*********************************************************************
 FUNCTION: calculate_path_score
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
double calculate_path_score(Gaze_Sequence *g_seq,
			    Gaze_Structure *gs) {
  
  int idx; 
  Feature_Relation *reg_info;
  Feature *src, *tgt;
  double total_score, trans_score, len_pen, seg_score;
  Seg_Results *s_res = new_Seg_Results( gs->seg_dict->len );

  total_score = 0.0;

  for( idx=0; idx < g_seq->path->len - 1; idx++) { 
    /* for the score to mean anything, all paths must begin with "BEGIN"
       and end with "END." Therefore ignoring the local score of the first 
       feature (done here) has no effect, since the score of "BEGIN" is 0
       (and has to be for the DP to work */

    src = index_Array( g_seq->path, Feature *, idx );
    tgt = index_Array( g_seq->path, Feature *, idx + 1 );
    reg_info = index_Array( index_Array( gs->feat_info, Feature_Info *, tgt->feat_idx )->sources, 
			      Feature_Relation *, 
			      src->feat_idx);

    trans_score = len_pen = 0.0;
    
    seg_score = calculate_segment_score( g_seq, src, tgt, gs, s_res );
    trans_score += seg_score;
    
    if (reg_info->len_fun != NULL) {
      Length_Function *lf = 
	index_Array(gs->length_funcs, Length_Function *, *(reg_info->len_fun));

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
void forwards_calc( Gaze_Sequence *g_seq,
		    Gaze_Structure *gs,
		    boolean use_pruning,
		    Gaze_Output *g_out) {
  
  int ft_idx, prev_idx;
  Array *temp;
  Feature *prev_feat;
  
  Gaze_DP_struct *g_res = new_Gaze_DP_struct( gs->feat_dict->len, 
					      gs->seg_dict->len,
					      0 );
  
#ifdef TRACE
  fprintf(stderr, "\nForward calculation:\n\n");
#endif
  
  if ( g_out->sample_gene)
    srand( time(NULL) );

  for (ft_idx = 1; ft_idx < g_seq->features->len; ft_idx++) {

    /* only do the dp if a path is not already present for the sequence,
       or if a path is present by the user wanted to see probabilities */

    if (g_seq->path == NULL || g_out->probability) {
      prev_idx = ft_idx - 1;
      prev_feat = index_Array( g_seq->features, Feature *, prev_idx );
      temp = g_res->feats[prev_feat->feat_idx][MOD3(prev_feat->adj_pos.s)];
      append_val_Array( temp, prev_idx );

      if (g_out->sample_gene || g_out->regions || g_out->probability)
	scan_through_sources_dp( g_seq,
				 gs, 
				 ft_idx,  
				 g_res, 
				 use_pruning,
				 g_out);
      else
	scan_through_sources_for_max_only( g_seq,
					   gs, 
					   ft_idx,  
					   g_res,
					   use_pruning);
    }

    index_Array( g_seq->features, Feature *, ft_idx )->forward_score = g_res->score;
    index_Array( g_seq->features, Feature *, ft_idx )->path_score = g_res->pth_score;
    index_Array( g_seq->features, Feature *, ft_idx )->trace_pointer = g_res->pth_trace;
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
void backwards_calc( Gaze_Sequence *g_seq,
		     Gaze_Structure *gs,
		     boolean use_pruning) {

  int ft_idx, prev_idx;
  Feature *prev_feat;
  Array *temp;

  Gaze_DP_struct *g_res = new_Gaze_DP_struct( gs->feat_dict->len, 
					      gs->seg_dict->len,
					      g_seq->features->len - 1 );

  g_res->last_selected = g_seq->features->len + 1;

#ifdef TRACE
  fprintf(stderr, "\nBackward calculation:\n\n");
#endif
  
  for (ft_idx = g_seq->features->len-2; ft_idx >= 0; ft_idx--) {
    /* push the index if the last feature onto the list of
       sorted indices */
    prev_idx = ft_idx + 1;
    prev_feat = index_Array( g_seq->features, Feature *, prev_idx );
    temp = g_res->feats[prev_feat->feat_idx][MOD3(prev_feat->adj_pos.e)];
    append_val_Array( temp, prev_idx );

    scan_through_targets_dp( g_seq,
			     gs,
			     ft_idx,  
			     g_res,
			     use_pruning);

    index_Array( g_seq->features, Feature *, ft_idx )->backward_score = g_res->score;

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
void scan_through_sources_dp( Gaze_Sequence *g_seq,
			      Gaze_Structure *gs,
			      int tgt_idx,
			      Gaze_DP_struct *g_res,
			      boolean use_pruning,
			      Gaze_Output *g_out) {
  
  int src_type, src_idx, kill_idx, max_index = 0; /* Initialsied to get arounc gcc warnings */
  int frame, k, index_count[3];
  int last_necessary_idx, local_fringe, last_idx_for_frame[3];
  int left_pos, right_pos, distance;
  Killer_Feature_Qualifier *kq;

  boolean touched_score, touched_score_local;
  Array *all_scores = NULL;
  Array *all_indices = NULL;  /* Initialised to get around gcc warnings */
  Feature_Info *tgt_info;
  Feature_Relation *reg_info;
  Feature *src, *tgt;

  boolean gone_far_enough = FALSE;
  double max_forpluslen = 0.0;
  double max_score = NEG_INFINITY;
  double max_forward = NEG_INFINITY;
  int *killer_source_dna = NULL;
  int *danger_source_dna = NULL;

  g_res->pth_score = 0.0;
  g_res->pth_trace = 0,0;
  g_res->score = 0.0;

  tgt = index_Array( g_seq->features, Feature *, tgt_idx );
  tgt_info = index_Array( gs->feat_info, Feature_Info *, tgt->feat_idx );
  right_pos = tgt->adj_pos.e;

#ifdef TRACE 
  fprintf( stderr, "Target %d %s %d %d %.3f", tgt_idx,
	   index_Array(gs->feat_dict, char *, tgt->feat_idx ), 
	   tgt->real_pos.s, tgt->real_pos.e, tgt->score );
  if (TRACE > 1)
    fprintf( stderr, "\n" );
#endif

    all_scores = new_Array( sizeof(double), TRUE );
    all_indices = new_Array( sizeof(int), TRUE );

  if (! tgt->invalid) {
    touched_score = FALSE;

    /* set up the boundaries for the scan. We do not want to go past:
       1. The last forced feature */

    last_necessary_idx = 0; 

    if (g_res->last_selected > last_necessary_idx)
      last_necessary_idx = g_res->last_selected;
    
    /* Look through the sources themselves */
    
    for ( src_type = 0; src_type < tgt_info->sources->len; src_type++) {
      if ((reg_info = index_Array(tgt_info->sources, Feature_Relation *, src_type)) != NULL) {
	Array **feats = g_res->feats[src_type];
	
	for(frame = 0; frame < 3; frame++) {
	  
	  last_idx_for_frame[frame] = last_necessary_idx;
	  
	  /* first, identify the killers local to this feature pair, and make
	     sure that our search back through the sources does not go past a killer */
	  
	  if (reg_info->kill_feat_quals != NULL) {
	    
	    for(kill_idx=0; kill_idx < reg_info->kill_feat_quals->len; kill_idx++) {
	      if ( (kq = index_Array( reg_info->kill_feat_quals,
					Killer_Feature_Qualifier *,
					kill_idx )) != NULL) {

		Array *apt_list;
		int k = 0;
		boolean more_frames = TRUE;

		while (more_frames) {
		  if (kq->has_tgt_phase) {
		    apt_list = (Array *) g_res->feats[(int)kq->feat_idx][MOD3(right_pos - kq->phase + 1)];
		    /* rationale: (right_pos - left_pos + 1) % 3 == phase -->
		       (right_pos - {left_pos % 3} + 1) % 3 == phase -->
		       (right_pos - {left_pos % 3} + 1) % 3 - phase == 0 -->
		       (right_pos - {left_pos % 3} + 1 - phase) % 3 == 0 -->
		       (right_pos - phase + 1) % 3 - {left_pos % 3} == 0 -->
		       (right_pos - phase + 1) % 3 == {left_pos % 3} */
		    more_frames = FALSE;
		  }		  
		  else if (kq->has_src_phase) {
		    /* the frame calculation here is slightly hacky; we have to allow for the
		       fact that we need the distance from the source forward to the apt. killer.
		       BUT the the killers are stored by the frame of their adjusted START, 
		       rather than their end. So, we rely on the fact that all killers that have
		       a phase have width that is 3-mutlple, which is sensible */
		    apt_list = (Array *) g_res->feats[(int)kq->feat_idx][MOD3(frame + kq->phase)];
		    more_frames = FALSE;
		  }
		  else {
		    /* phaseless killer - need to check all frames */
		    apt_list = (Array *) g_res->feats[(int)kq->feat_idx][k++];
		    if (k > 2)
		      more_frames = FALSE;
		  }
		  
		  if (apt_list->len > 0) {
		    /* first search back for the first occurrence that does not overlap with target */
		    int this_kill_idx = apt_list->len - 1;
		    int boundary_index = index_Array( apt_list, int, this_kill_idx );
		    int local_idx;
		    Feature *killer_feat = index_Array(g_seq->features, Feature *, boundary_index );
		    while ( killer_feat != NULL && killer_feat->real_pos.e > tgt->adj_pos.e) {
		      if (--this_kill_idx >= 0) {
			boundary_index = index_Array( apt_list, int, this_kill_idx );
			killer_feat = index_Array(g_seq->features, Feature *, boundary_index );
		      }
		      else
			killer_feat = NULL;
		    }
		    /* now search back for sources beyond the killer of this type that overlap the killer */

		    if (killer_feat != NULL) {
		      local_idx = boundary_index - 1;
		      while (  local_idx >= 0 ) {
			Feature *candidate = index_Array(g_seq->features, Feature *, local_idx );
			if (candidate->adj_pos.s <= killer_feat->real_pos.s)
			break;
			else if (candidate->feat_idx == src_type)
			  boundary_index = local_idx;
			local_idx--;
		      }

		      if (boundary_index > last_idx_for_frame[frame]) 
			last_idx_for_frame[frame] = boundary_index;
		    }
		  }
		}
	      }
	    }
	  }
	  
	  /* finally, make sure that we do not proceed past the fringe for
	     this feature pair */
	  if (g_res->fringes[(int)tgt->feat_idx][src_type][MOD3(tgt->real_pos.s)] > last_idx_for_frame[frame])
	    last_idx_for_frame[frame] = g_res->fringes[(int)tgt->feat_idx][src_type][MOD3(tgt->real_pos.s)];
	}
	
	/* Before actually scanning through the features themselves, we need to check
	   if there are potential dna killers. If so, set flags for the source
	   dna entries that will cause problems */
	
	if (reg_info->kill_dna_quals != NULL) {
	  danger_source_dna = (int *) malloc0_util( gs->motif_dict->len * sizeof(int) );
	  
	  for(k=0; k < reg_info->kill_dna_quals->len; k++) {
	    Killer_DNA_Qualifier *kdq = index_Array( reg_info->kill_dna_quals, 
						       Killer_DNA_Qualifier *,
						       k );
	    
	    danger_source_dna[(int)kdq->src_dna] = 1; 
	    
	    if (tgt->dna >= 0 && tgt->dna == kdq->tgt_dna) {
	      if (killer_source_dna == NULL) 
		killer_source_dna = (int *) malloc0_util( gs->motif_dict->len * sizeof(int) );
	      killer_source_dna[(int)kdq->src_dna] = 1;
	    }
	  }
	}
	
	/* at this point, we have the list of features that need to be processed (feats),
	   and the index that we must not proceed past in each frame. We can now process
	   the features themselves, in a frame-dependent or frame-independent way */
	
#ifdef TRACE
	if (TRACE > 1)
	  fprintf( stderr, "  %s (fringes: %d %d %d)\n",
		   index_Array(gs->feat_dict, char *, src_type ), 
		   last_idx_for_frame[0], last_idx_for_frame[1], last_idx_for_frame[2] );
#endif

	for(k=0; k < 3; k++) 
	  index_count[k] = feats[k]->len - 1; 
	
	frame = reg_info->phase != NULL ? MOD3(right_pos - *(reg_info->phase) + 1) : 0;


	max_forpluslen = NEG_INFINITY;
	touched_score_local = FALSE;
	/* the following aggressively assumes that if this target has no 
	   potential sources for this source type, then we need go no 
	   further back than the index of the target itself when consdering
	   furture instances of the target */
	local_fringe = tgt_idx;
	gone_far_enough = FALSE;

	while( ! gone_far_enough ) {

	  if (reg_info->phase == NULL) {
	    /* For frameless feature pairs, we need to examine all frames, but 
	       for the pruning to work effectively, the featured need to be examined 
	       in order. Therefore, we have to re-create the original list
	       by effectively a merge sort. */
	    boolean gotone = FALSE;

	    for (k=0; k < 3; k++) {
	      if (index_count[k] >= 0) {
		if (gotone) {
		  if (index_Array( feats[k], int, index_count[k] ) > 
		      index_Array( feats[frame], int, index_count[frame] ) ) {
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
	  
	  /* The following tests if there was anything in the list at all
	     For targets very close to the start of the sequence, frame
	     can be < 0, but there will be no sources of this type for
	     such targets anyway */
	  if (frame < 0 || index_count[frame] < 0) {
	    gone_far_enough = TRUE;
	    continue;
	  }
	  
	  src_idx = index_Array( feats[frame], int, index_count[frame]-- );
	  
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

	  src = index_Array( g_seq->features, Feature *, src_idx );
	  
#ifdef TRACE
	  if (TRACE > 1)
	    fprintf( stderr, "     Source %d %s %d %d ", src_idx,
		     index_Array(gs->feat_dict, char *, src_type ),
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
		if (killer_source_dna == NULL || src->dna < 0 || ! killer_source_dna[(int)src->dna]) {
		  double trans_score, len_pen, seg_score, forward_temp, viterbi_temp;
		  Length_Function *lf = NULL;
		  trans_score = len_pen = forward_temp = viterbi_temp = 0.0;

		  seg_score = calculate_segment_score( g_seq, src, tgt, gs, g_res->seg_res );
		  trans_score += seg_score;
		  
		  if (reg_info->len_fun != NULL) {
		    lf = index_Array(gs->length_funcs, Length_Function *, *(reg_info->len_fun));
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
		  
		  forward_temp = src->forward_score 
		    + trans_score
		    + tgt->score;
		    
		  append_val_Array( all_scores, forward_temp);
		  append_val_Array( all_indices, src_idx);
		  
		  if (! touched_score || (forward_temp > max_forward))
		    max_forward = forward_temp;
		  
		  if (use_pruning) {
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
			  || ! danger_source_dna[(int)src->dna]) {
			
			/* strictly speaking, it is only sound to register this source
			   as "dominant" if we are into the monotonic part of the length
			   function (i.e. the point past which the penalty never decreases 
			   with increasing distance). Therefore, we update the fringe, but
			   don't flag touched_local_score. The efect of this is that the
			   fringe will be updated for all scoring sources until we get past
			   the point of monotonicity, at which point the pruning kicks in */
			
			if (lf == NULL || (lf->becomes_monotonic && lf->monotonic_point <= distance)) { 
			  /* finally, check that the source is not likely to be involved in
			     an exact segment either here or at some point down the line. If so,
			     it's unfair to consider the source as omnipotent
			  */
			  
			  if (! g_res->seg_res->has_exact_at_src) { 
			    /* add back in the length penalty, because when judging for dominance, 
			       the length penalty will be different for future features */
			    
			    max_forpluslen = forward_temp + len_pen;
			    touched_score_local = TRUE;
			  }
			} 
		      }
		      
		      local_fringe = src_idx;
		    }
		    else {
		      /* compare this one to max_forward, to see if it is dominated */
		      if (forward_temp + len_pen > max_forpluslen 
			  && (danger_source_dna == NULL 
			      || src->dna < 0  
			      || ! danger_source_dna[(int)src->dna])
			  && ! g_res->seg_res->has_exact_at_src)
			max_forpluslen = forward_temp + len_pen;
		      
		      if ( max_forpluslen - (forward_temp + len_pen) < 25.0 
			   || (g_res->seg_res->has_exact_at_src 
			       && g_res->seg_res->exact_extends_beyond_tgt))
			local_fringe = src_idx;
		    }
		  }
		
		  touched_score = TRUE;		  
#ifdef TRACE
		  if (TRACE > 1) 
		    fprintf( stderr, "scre: v=%.3f, f=%.8f (seg:%.5f len:%.3f)\n",
			     viterbi_temp, forward_temp, seg_score, len_pen );
#endif
		  if (g_out->regions) {
		    if (reg_info->out_qual != NULL && reg_info->out_qual->need_to_print) {
		      
		      double reg_score = trans_score;
		      
		      if (g_out->probability) 
			reg_score = exp( src->forward_score + 
					 trans_score + tgt->score +				   
					 tgt->backward_score - 
					 g_seq->beg_ft->backward_score );
		      
		      if (! g_out->use_threshold || reg_score >= g_out->threshold)
			fprintf(g_out->fh, "%s\tGAZE\t%s\t%d\t%d\t%.5f\t%s\t%s\t\n",
				g_seq->seq_name, 
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
		  if (use_pruning)
		    local_fringe = src_idx;
		  
#ifdef TRACE
		  if (TRACE > 1)
		    fprintf( stderr, "KILLED_BY_DNA\n" );
#endif
		}
	      } /* if min dist */
	      else {
		/* source might not be too close for future incidences, so update fringe index */
		if (use_pruning)
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
	      gone_far_enough = TRUE;
	    }
	  }
#ifdef TRACE
	  else 
	    if (TRACE > 1)
	      fprintf( stderr, "INVALID\n" );
#endif
	} /* while !gone_far_enough */
      
	if (use_pruning) {
	  /* We conservatively only prune in the frame of the target if this
	     feature pair has a phase constraint, or if there are potential
	     killers (which also might have a phase constraint). Otherwise, 
	     prune in all frames */
	  if (reg_info->phase != NULL || reg_info->kill_feat_quals != NULL) {
	    g_res->fringes[(int)tgt->feat_idx][src_type][MOD3(tgt->real_pos.s)] = local_fringe;
	  }
	  else {
	    for(k=0; k < 3; k++)
	      g_res->fringes[(int)tgt->feat_idx][src_type][k] = local_fringe;
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

    /* update the position of the last forced feature. */
    
    if (tgt->is_selected)
      g_res->last_selected = tgt_idx;
    
    if (touched_score) {      
      /* the trick of subtracting the max before exponentiating avoids
	 overflow errors. Just need to add it back when logging back down */
      for (src_idx=0; src_idx < all_scores->len; src_idx++) {
	g_res->score += exp( index_Array(all_scores, double, src_idx) 
			     - max_forward );  
	  
      }
      g_res->score = log( g_res->score ) + max_forward;	
            
      if (g_out->sample_gene) {
	double random_number = (double) rand() / (double) RAND_MAX;
	double sum = 0.0;

	for(src_idx=0; src_idx < all_indices->len; src_idx++) {
	  double ft_prob = exp( index_Array( all_scores, double, src_idx ) -
				g_res->score ); 
	  sum += ft_prob;
	  
	  if (sum >= random_number) {
	    Feature *tmp;
	    
	    g_res->pth_trace = index_Array( all_indices, int, src_idx);

	    tmp = index_Array( g_seq->features, Feature *, g_res->pth_trace ); 

	    /* note that we are just returning the transition + local score 
	       here for the source-target pair. This can be accumulated
	       by the caller to get the score of the sample path. Alternatively,
	       the score of the path can be computed with calculate_path_score */
	    
	    g_res->pth_score = index_Array( all_scores, double, src_idx) -
	      tmp->forward_score;
	    
	    break;
	  }
	}
      }
      else {
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

  free_Array( all_indices, TRUE );
  free_Array( all_scores, TRUE );
}



/*********************************************************************
 FUNCTION: scan_through_targets_dp
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void scan_through_targets_dp( Gaze_Sequence *g_seq,
			      Gaze_Structure *gs,
			      int src_idx,
			      Gaze_DP_struct *g_res,
			      boolean use_pruning) {

  int tgt_type, tgt_idx, kill_idx;
  int frame, k, index_count[3]; 
  int left_pos, right_pos, distance;
  int last_necessary_idx, local_fringe, last_idx_for_frame[3];
  Killer_Feature_Qualifier *kq;

  boolean touched_score, touched_score_local;
  Array *all_scores;
  Feature_Info *tgt_info;
  Feature_Relation *reg_info;
  Feature *src, *tgt;

  boolean gone_far_enough = FALSE;
  double max_backward = 0.0;
  double max_backpluslen = 0.0;
  int *killer_target_dna = NULL;
  int *danger_target_dna = NULL;

  g_res->score = 0.0;

  src = index_Array( g_seq->features, Feature *, src_idx );
  left_pos = src->adj_pos.s;

  /* if the user specified unusual offsets, it may be that this feature is
     "off the end of the sequence" when viewed as a target. */
  /*
  if (left_pos > g_seq->end_ft->real_pos.e ||
      src->is_antiselected ) {
    src->invalid = TRUE;
  }
  */
#ifdef TRACE
  fprintf( stderr, "Source %d %s %d %d %.3f", src_idx,
	   index_Array(gs->feat_dict, char *, src->feat_idx ), 
	   src->real_pos.s, src->real_pos.e, src->score ); 
  
  if (TRACE > 1)
    fprintf( stderr, "\n" );
#endif
  
  all_scores = new_Array( sizeof(double), TRUE);

  if (! src->invalid ) { 
    touched_score = FALSE;

    /* set up the boundaries for the scan. We do not want to go past 
       1. The last forced feature */

    last_necessary_idx = g_seq->features->len - 1; 

    if (g_res->last_selected < last_necessary_idx)
      last_necessary_idx = g_res->last_selected;
    
    /* Look through the targets themselves */
    
    for ( tgt_type = 0; tgt_type < gs->feat_dict->len; tgt_type++) {
      tgt_info = index_Array( gs->feat_info, Feature_Info *, tgt_type );

      /* this guards against looking for sources of BEGIN, of which there are none */
      if (tgt_info->sources == NULL)
	continue;

      if ((reg_info = index_Array(tgt_info->sources, Feature_Relation *, src->feat_idx)) != NULL) {
	Array **feats = g_res->feats[tgt_type];

	for(frame = 0; frame < 3; frame ++) {
	  
	  last_idx_for_frame[frame] = last_necessary_idx;

	  /* first, identify the killers local to this feature pair, and make
	     sure that our search forward through the targets does not go past 
	     a killer */

	  if (reg_info->kill_feat_quals != NULL) {

	    for(kill_idx=0; kill_idx < reg_info->kill_feat_quals->len; kill_idx++) {
	      if ( (kq = index_Array( reg_info->kill_feat_quals,
				      Killer_Feature_Qualifier *,
				      kill_idx )) != NULL) {
		Array *apt_list;
		int k = 0;
		boolean more_frames = TRUE;
		
		while( more_frames) {
		  if (kq->has_src_phase) {
		    apt_list = (Array *) g_res->feats[(int)kq->feat_idx][MOD3(left_pos + kq->phase - 1)];
		    /* rationale: (right_pos - left_pos + 1) % 3 == phase -->
		       (right_pos - left_pos + 1) % 3 - phase == 0 -->
		       (right_pos - left_pos + 1 - phase) % 3 == 0 -->
		       (left_pos - right_pos - 1 + phase) % 3 == 0 -->
		       (left_pos - {right_pos % 3} + phase - 1) % 3 == 0 -->
		       (left_pos + phase - 1) % 3 - {right_pos % 3} == 0 -->
		       (left_pos + phase - 1) % 3 == {right_pos % 3} */
		    more_frames = FALSE;
		  }
		  else if (kq->has_tgt_phase) {
		    /* the frame calcualtion here is slightly hacky; we have to allow for the
		       fact that we need the distance from the target back to the apt. killer.
		       BUT, the the killers are stored by the frame of their adjusted END, 
		       rather than their start. So, we rely on the fact that all killers that have 
		       a phase are width 3, which might not be so unreasonable */
		    apt_list = (Array *) g_res->feats[(int)kq->feat_idx][MOD3(frame + 3 - kq->phase)];
		    more_frames = FALSE;
		  }	
		  else {
		    /* phaseless killer - need to check all frames */
		    apt_list = (Array *) g_res->feats[kill_idx][k++];
		    if (k > 2)
		     more_frames = FALSE;
		  }
		  
		  if (apt_list->len > 0) {
		    /* first search forward for the first occurrence that does not overlap with the src */
		    int this_kill_idx = apt_list->len - 1;
		    int boundary_index = index_Array( apt_list, int, this_kill_idx );
		    int local_idx;
		    Feature *killer_feat = index_Array(g_seq->features, Feature *, boundary_index );

		    while ( killer_feat != NULL && killer_feat->real_pos.s < src->adj_pos.s ) {
		      if (--this_kill_idx >= 0) {
			boundary_index = index_Array( apt_list, int, this_kill_idx );
			killer_feat = index_Array(g_seq->features, Feature *, boundary_index );
		      }
		      else
			killer_feat = NULL;
		    }
		    /* now search forward for targets beyond the killer of this type that overlap the killer */
		    
		    if (killer_feat != NULL) {
		      local_idx = boundary_index + 1;
		      while (  local_idx < g_seq->features->len ) {
			Feature *candidate = index_Array(g_seq->features, Feature *, local_idx );
			if (candidate->adj_pos.e >= killer_feat->real_pos.e)
			  break;
			else if (candidate->feat_idx == tgt_type)
			  boundary_index = local_idx;
			local_idx++;
		      }
		      
		      if (boundary_index < last_idx_for_frame[frame]) 
		      last_idx_for_frame[frame] = boundary_index;
		    }
		  }
		}
	      }
	    }
	  }
	  
	  /* finally, make sure that we do not proceed past the fringe for
	     this feature pair */
	  if (g_res->fringes[(int)src->feat_idx][tgt_type][MOD3(src->real_pos.s)] < last_idx_for_frame[frame])
	    last_idx_for_frame[frame] = g_res->fringes[(int)src->feat_idx][tgt_type][MOD3(src->real_pos.s)];
	}


	/* Before actually scanning through the features themselves, we need to check
	   if there are potential dna killers. If so, set flags for the source
	   dna entries that will cause problems */
	
	if (reg_info->kill_dna_quals != NULL) {
	  danger_target_dna = (int *) malloc0_util( gs->motif_dict->len * sizeof(int) );
	  
	  for(k=0; k < reg_info->kill_dna_quals->len; k++) {
	    Killer_DNA_Qualifier *kdq = index_Array( reg_info->kill_dna_quals, 
						       Killer_DNA_Qualifier *,
						       k );
	    
	    danger_target_dna[(int)kdq->tgt_dna] = 1;
	    
	    if (src->dna >= 0 && src->dna == kdq->src_dna) {
	      if (killer_target_dna == NULL)
		killer_target_dna = (int *) malloc0_util( gs->motif_dict->len * sizeof(int) );
	      killer_target_dna[(int)kdq->tgt_dna] = 1;
	    }
	  }
	}
	
	/* at this point, we have the list of features that need to be processed (feats),
	   and the index that we must not proceed past in each frame. We can now process
	   the features themselves, in a frame-dependent or frame-independent way */

#ifdef TRACE
	if (TRACE > 1)
	  fprintf( stderr, "  %s (fringes: %d %d %d)\n",
		   index_Array(gs->feat_dict, char *, tgt_type ), 
		   last_idx_for_frame[0], last_idx_for_frame[1], last_idx_for_frame[2] );
#endif

	for(k=0; k < 3; k++) 
	  index_count[k] = feats[k]->len - 1; 
	
	frame = reg_info->phase != NULL ? MOD3(left_pos + *(reg_info->phase) - 1) : 0;

	max_backpluslen = NEG_INFINITY;
	touched_score_local = FALSE;
	/* the following aggressively assumes that if this target has no 
	   potential sources for this source type, then we need go no 
	   further back than the index of the target itself when consdering
	   furture instances of the target */
	local_fringe = src_idx;
	gone_far_enough = FALSE;

	while(! gone_far_enough) {
	  
	  if (reg_info->phase == NULL) {
	    /* For frameless feature pairs, we need to examine all frames, but 
	       for the pruning to work effectively, the featured need to be examined 
	       in order. Therefore, we have to re-create the original list
	       by effectively a merge sort. */
	    boolean gotone = FALSE;
	    
	    for (k=0; k < 3; k++) {
	      if (index_count[k] >= 0) {
		if (gotone) {
		  if (index_Array( feats[k], int, index_count[k] ) < 
		      index_Array( feats[frame], int, index_count[frame] ) ) {
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
	  
	  /* The following tests if there was anything in the list at all
	     For targets very close to the start of the sequence, frame
	     can be < 0, but there will be no sources of this type for
	     such targets anyway */
	  if (frame < 0 || index_count[frame] < 0) {
	    gone_far_enough = TRUE;
	    continue;
	  }
	  
	  tgt_idx = index_Array( feats[frame], int, index_count[frame]-- );

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

	  tgt = index_Array( g_seq->features, Feature *, tgt_idx );
	    
#ifdef TRACE
	  if (TRACE > 1)
	    fprintf( stderr, "  Target %d %s %d %d  ", tgt_idx,
		     index_Array(gs->feat_dict, char *, tgt->feat_idx ), 
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
		if (killer_target_dna == NULL || tgt->dna < 0 || ! killer_target_dna[(int)tgt->dna]) {
		  double trans_score, len_pen, seg_score, backward_temp;
		  Length_Function *lf = NULL;
		  trans_score = len_pen = 0.0;
		  
		  seg_score = calculate_segment_score( g_seq, src, tgt, gs, g_res->seg_res );
		  trans_score += seg_score;
		  
		  if (reg_info->len_fun != NULL) {
		      lf = index_Array(gs->length_funcs, Length_Function *, *(reg_info->len_fun));
		      len_pen = apply_Length_Function( lf, distance );
		  }
		  trans_score -= len_pen;
		  
		  backward_temp = tgt->backward_score
		    + trans_score
		    + tgt->score;
		  
		  if (! touched_score || backward_temp > max_backward)
		    max_backward = backward_temp;
		  
		  append_val_Array( all_scores, backward_temp);

		  if (use_pruning) {
		    if (! touched_score_local ) {

		      if (danger_target_dna == NULL 
			  || tgt->dna < 0  
			  || ! danger_target_dna[(int)src->dna]) {

			if (lf == NULL || (lf->becomes_monotonic && lf->monotonic_point <= distance)) {
			  if (! g_res->seg_res->has_exact_at_tgt) {

			    touched_score_local = TRUE;
			    max_backpluslen = backward_temp + len_pen;
			  }
			}
		      }

		      local_fringe = tgt_idx;

		    }
		    else {
		      /* compare this one to max_forward, to see if it is dominated */
		      if (backward_temp + len_pen > max_backpluslen
			  && (danger_target_dna == NULL 
			      || tgt->dna < 0  
			      || ! danger_target_dna[(int)tgt->dna])
			  && !g_res->seg_res->has_exact_at_tgt) 
			max_backpluslen = backward_temp + len_pen;
		      
		      if ( max_backpluslen - (backward_temp + len_pen) < 25.0
			   || (g_res->seg_res->has_exact_at_tgt 
			       && g_res->seg_res->exact_extends_beyond_src) )
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
		  if (use_pruning)
		    local_fringe = tgt_idx;

#ifdef TRACE
		  if (TRACE > 1)
		    fprintf( stderr, "KILLED_BY_DNA\n" );
#endif
		}  
	      } /* if min dist */
	      else {
		/* target might not be too close for future incidences, so update fringe index */
		if (use_pruning)
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
	      gone_far_enough = TRUE;
	    }
	  } /* if valid */
#ifdef TRACE
	  else {
	    if (TRACE > 1)
	      fprintf( stderr, "INVALID\n" );
	  }
#endif
	} /* while ! gone_far_enough */

	if (use_pruning) {
	  /* We conservatively only prune in the frame of the target if this
	     feature pair has a phase constraint, or if there are potential
	     killers (which also might have a phase constraint). Otherwise, 
	     prune in all frames */
	  if (reg_info->phase != NULL || reg_info->kill_feat_quals != NULL) {
	    g_res->fringes[(int)src->feat_idx][tgt_type][MOD3(src->real_pos.s)] = local_fringe;
	  }
	  else {
	    for(k=0; k < 3; k++)
	      g_res->fringes[(int)src->feat_idx][tgt_type][k] = local_fringe;
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

    /* update the position of the last forced feature. */
    
    if (src->is_selected)
      g_res->last_selected = src_idx;
    
    if (touched_score) {
      for (tgt_idx=0; tgt_idx < all_scores->len; tgt_idx++) {
	g_res->score += exp( index_Array(all_scores, double, tgt_idx) 
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

  free_Array( all_scores, TRUE );
   
}      



/*********************************************************************
 FUNCTION: scan_through_source_for_max_only
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void scan_through_sources_for_max_only( Gaze_Sequence *g_seq,
					Gaze_Structure *gs,
					int tgt_idx,
					Gaze_DP_struct *g_res,
					boolean use_pruning) {
  
  int src_type, src_idx, kill_idx, max_index = 0; /* Initialsied to get arounc gcc warnings */
  int frame, k, index_count[3];
  int last_necessary_idx, local_fringe, last_idx_for_frame[3];
  int left_pos, right_pos, distance;
  Killer_Feature_Qualifier *kq;

  boolean touched_score, touched_score_local;
  Feature_Info *tgt_info;
  Feature_Relation *reg_info;
  Feature *src, *tgt;

  boolean gone_far_enough = FALSE;
  double max_vit_plus_len = 0.0;
  double max_score = NEG_INFINITY;
  int *killer_source_dna = NULL;
  int *danger_source_dna = NULL;

  g_res->pth_score = 0.0;
  g_res->pth_trace = 0,0;
  g_res->score = 0.0;

  tgt = index_Array( g_seq->features, Feature *, tgt_idx );
  tgt_info = index_Array( gs->feat_info, Feature_Info *, tgt->feat_idx );
  right_pos = tgt->adj_pos.e;

#ifdef TRACE 
  fprintf( stderr, "Target %d %s %d %d %.3f", tgt_idx,
	   index_Array(gs->feat_dict, char *, tgt->feat_idx ), 
	   tgt->real_pos.s, tgt->real_pos.e, tgt->score );
  if (TRACE > 1)
    fprintf( stderr, "\n" );
#endif

  if (! tgt->invalid) {
    touched_score = FALSE;

    /* set up the boundaries for the scan. We do not want to go past:
       1. The last forced feature */

    last_necessary_idx = 0; 

    if (g_res->last_selected > last_necessary_idx)
      last_necessary_idx = g_res->last_selected;
    
    /* Look through the sources themselves */
    
    for ( src_type = 0; src_type < tgt_info->sources->len; src_type++) {
      if ((reg_info = index_Array(tgt_info->sources, Feature_Relation *, src_type)) != NULL) {
	Array **feats = g_res->feats[src_type];
	
	for(frame = 0; frame < 3; frame++) {
	  last_idx_for_frame[frame] = last_necessary_idx;
	  
	  /* first, identify the killers local to this feature pair, and make
	     sure that our search back through the sources does not go past a killer */
	  	  
	  if (reg_info->kill_feat_quals != NULL) {
	    
	    for(kill_idx=0; kill_idx < reg_info->kill_feat_quals->len; kill_idx++) {
	      if ( (kq = index_Array( reg_info->kill_feat_quals,
					Killer_Feature_Qualifier *,
					kill_idx )) != NULL) {

		Array *apt_list;
		int k = 0;
		boolean more_frames = TRUE;

		while (more_frames) {
		  if (kq->has_tgt_phase) {
		    apt_list = (Array *) g_res->feats[(int)kq->feat_idx][MOD3(right_pos - kq->phase + 1)];
		    /* rationale: (right_pos - left_pos + 1) % 3 == phase -->
		       (right_pos - {left_pos % 3} + 1) % 3 == phase -->
		       (right_pos - {left_pos % 3} + 1) % 3 - phase == 0 -->
		       (right_pos - {left_pos % 3} + 1 - phase) % 3 == 0 -->
		       (right_pos - phase + 1) % 3 - {left_pos % 3} == 0 -->
		       (right_pos - phase + 1) % 3 == {left_pos % 3} */
		    more_frames = FALSE;
		  }		  
		  else if (kq->has_src_phase) {
		    /* the frame calculation here is slightly hacky; we have to allow for the
		       fact that we need the distance from the source forward to the apt. killer.
		       BUT the the killers are stored by the frame of their adjusted START, 
		       rather than their end. So, we rely on the fact that all killers that have
		       a phase have width that is 3-mutlple, which is sensible */
		    apt_list = (Array *) g_res->feats[(int)kq->feat_idx][MOD3(frame + kq->phase)];
		    more_frames = FALSE;
		  }
		  else {
		    /* phaseless killer - need to check all frames */
		    apt_list = (Array *) g_res->feats[(int)kq->feat_idx][k++];
		    if (k > 2)
		      more_frames = FALSE;
		  }
		  
		  if (apt_list->len > 0) {
		    /* first search back for the first occurrence that does not overlap with target */
		    int this_kill_idx = apt_list->len - 1;
		    int boundary_index = index_Array( apt_list, int, this_kill_idx );
		    int local_idx;
		    Feature *killer_feat = index_Array(g_seq->features, Feature *, boundary_index );
		    while ( killer_feat != NULL && killer_feat->real_pos.e > tgt->adj_pos.e) {
		      if (--this_kill_idx >= 0) {
			boundary_index = index_Array( apt_list, int, this_kill_idx );
			killer_feat = index_Array(g_seq->features, Feature *, boundary_index );
		      }
		      else
			killer_feat = NULL;
		    }
		    /* now search back for sources beyond the killer of this type that overlap the killer */

		    if (killer_feat != NULL) {
		      local_idx = boundary_index - 1;
		      while (  local_idx >= 0 ) {
			Feature *candidate = index_Array(g_seq->features, Feature *, local_idx );
			if (candidate->adj_pos.s <= killer_feat->real_pos.s)
			break;
			else if (candidate->feat_idx == src_type)
			  boundary_index = local_idx;
			local_idx--;
		      }

		      if (boundary_index > last_idx_for_frame[frame]) 
			last_idx_for_frame[frame] = boundary_index;
		    }
		  }
		}
	      }
	    }
	  }

	  /* finally, make sure that we do not proceed past the fringe for
	     this feature pair */
	  if (g_res->fringes[(int)tgt->feat_idx][src_type][MOD3(tgt->real_pos.s)] > last_idx_for_frame[frame])
	    last_idx_for_frame[frame] = g_res->fringes[(int)tgt->feat_idx][src_type][MOD3(tgt->real_pos.s)];
	}
	
	/* Before actually scanning through the features themselves, we need to check
	   if there are potential dna killers. If so, set flags for the source
	   dna entries that will cause problems */
	
	if (reg_info->kill_dna_quals != NULL) {
	  danger_source_dna = (int *) malloc0_util( gs->motif_dict->len * sizeof(int) );
	  
	  for(k=0; k < reg_info->kill_dna_quals->len; k++) {
	    Killer_DNA_Qualifier *kdq = index_Array( reg_info->kill_dna_quals, 
						       Killer_DNA_Qualifier *,
						       k );
	    
	    danger_source_dna[(int)kdq->src_dna] = 1; 
	    
	    if (tgt->dna >= 0 && tgt->dna == kdq->tgt_dna) {
	      if (killer_source_dna == NULL) 
		killer_source_dna = (int *) malloc0_util( gs->motif_dict->len * sizeof(int) );
	      killer_source_dna[(int)kdq->src_dna] = 1;
	    }
	  }
	}
	
	/* at this point, we have the list of features that need to be processed (feats),
	   and the index that we must not proceed past in each frame. We can now process
	   the features themselves, in a frame-dependent or frame-independent way */
	
#ifdef TRACE
	if (TRACE > 1)
	  fprintf( stderr, "  %s (fringes: %d %d %d)\n",
		   index_Array(gs->feat_dict, char *, src_type ), 
		   last_idx_for_frame[0], last_idx_for_frame[1], last_idx_for_frame[2] );
#endif

	for(k=0; k < 3; k++) 
	  index_count[k] = feats[k]->len - 1; 
	
	frame = reg_info->phase != NULL ? MOD3(right_pos - *(reg_info->phase) + 1) : 0;

	max_vit_plus_len = NEG_INFINITY;
	touched_score_local = FALSE;
	/* the following aggressively assumes that if this target has no 
	   potential sources for this source type, then we need go no 
	   further back than the index of the target itself when consdering
	   furture instances of the target */
	local_fringe = tgt_idx;
	gone_far_enough = FALSE;

	while( ! gone_far_enough ) {

	  if (reg_info->phase == NULL) {
	    /* For frameless feature pairs, we need to examine all frames, but 
	       for the pruning to work effectively, the featured need to be examined 
	       in order. Therefore, we have to re-create the original list
	       by effectively a merge sort. */
	    boolean gotone = FALSE;
	    
	    for (k=0; k < 3; k++) {
	      if (index_count[k] >= 0) {
		if (gotone) {
		  if (index_Array( feats[k], int, index_count[k] ) > 
		      index_Array( feats[frame], int, index_count[frame] ) ) {
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
	  if (frame < 0 || index_count[frame] < 0) {
	    gone_far_enough = TRUE;
	    continue;
	  }
	  
	  src_idx = index_Array( feats[frame], int, index_count[frame]-- );
	  
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
	  
	  src = index_Array( g_seq->features, Feature *, src_idx );
	  
#ifdef TRACE
	  if (TRACE > 1)
	    fprintf( stderr, "     Source %d %s %d %d ", src_idx,
		     index_Array(gs->feat_dict, char *, src_type ),
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
		if (killer_source_dna == NULL || src->dna < 0 || ! killer_source_dna[(int)src->dna]) {
		  double trans_score, len_pen, seg_score, viterbi_temp;
		  Length_Function *lf = NULL;
		  trans_score = len_pen = viterbi_temp = 0.0;

		  seg_score = calculate_segment_score( g_seq, src, tgt, gs, g_res->seg_res );
		  trans_score += seg_score;
		  
		  if (reg_info->len_fun != NULL) {
		    lf = index_Array(gs->length_funcs, Length_Function *, *(reg_info->len_fun));
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
		  		  
		  if (! touched_score_local ) {
			
		    if (danger_source_dna == NULL 
			|| src->dna < 0  
			|| ! danger_source_dna[(int)src->dna]) {

		      if (lf == NULL || (lf->becomes_monotonic && lf->monotonic_point <= distance)) { 

			if (! g_res->seg_res->has_exact_at_src) {

			  max_vit_plus_len = viterbi_temp + len_pen;
			  touched_score_local = TRUE;
			}
		      } 
		    }

		    local_fringe = src_idx;
		  }
		  else {
		    /* compare this one to max_viterbi, to see if it is dominated */
		    if (viterbi_temp + len_pen > max_vit_plus_len) {
 
		      if ( (danger_source_dna == NULL 
			    || src->dna < 0  
			    || ! danger_source_dna[(int)src->dna])
			   && ! g_res->seg_res->has_exact_at_src )
			max_vit_plus_len = viterbi_temp + len_pen;		      		      
		      
			local_fringe = src_idx;
		    }
		    else if (g_res->seg_res->has_exact_at_src 
			     && g_res->seg_res->exact_extends_beyond_tgt)
		      local_fringe = src_idx;
		  }		  
  
		  touched_score = TRUE;
		  
#ifdef TRACE
		  if (TRACE > 1) 
		    fprintf( stderr, "scre: v=%.3f (seg:%.5f len:%.3f)\n",
			     viterbi_temp, seg_score, len_pen );
#endif
		} /* if killed by DNA */
		else {
		  /* source might not be killed for future incidences, so update fringe index */
		  local_fringe = src_idx;
		  
#ifdef TRACE
		  if (TRACE > 1)
		    fprintf( stderr, "KILLED_BY_DNA\n" );
#endif
		}
	      } /* if min dist */
	      else {
		/* source might not be too close for future incidences, so update fringe index */
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
	      gone_far_enough = TRUE;
	    }
	  }
#ifdef TRACE
	  else 
	    if (TRACE > 1)
	      fprintf( stderr, "INVALID\n" );
#endif
	} /* while !gone_far_enough */
      
	/* We conservatively only prune in the frame of the target if this
	   feature pair has a phase constraint, or if there are potential
	   killers (which also might have a phase constraint). Otherwise, 
	   prune in all frames */
	if (use_pruning) {
	  if (reg_info->phase != NULL || reg_info->kill_feat_quals != NULL) {
	    g_res->fringes[(int)tgt->feat_idx][src_type][MOD3(tgt->real_pos.s)] = local_fringe;
	  }
	  else {
	    for(k=0; k < 3; k++)
	      g_res->fringes[(int)tgt->feat_idx][src_type][k] = local_fringe;
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

    /* update the position of the last forced feature. */
    
    if (tgt->is_selected)
      g_res->last_selected = tgt_idx;
    
    if (touched_score) {

      g_res->pth_trace = max_index;
      g_res->pth_score = max_score;	
      
      
#ifdef TRACE
      fprintf(stderr, "  RESULT: v=%.3f, max=%d\n", 
	      g_res->pth_score,
	      max_index);
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
void trace_back_general ( Gaze_Sequence *g_seq ) {
  
  int i;
  Feature *temp;

  Array *stack = new_Array( sizeof(Feature *), TRUE);
  Array *feat_path = new_Array( sizeof(Feature *), TRUE);
  int pos = g_seq->features->len - 1;  

  temp = index_Array( g_seq->features, Feature *, pos );
  append_val_Array( stack, temp );

  while (pos > 0) {
    pos = temp->trace_pointer;    
    
    temp = index_Array( g_seq->features, Feature *, pos );
    append_val_Array( stack, temp );
  }

  if (pos == 0) {
    for (i=stack->len-1; i>=0; i--) {
      temp =  index_Array( stack, Feature *, i );
      append_val_Array( feat_path, temp);
    }
  }
  else {
    free_Array( feat_path, TRUE );
    feat_path = NULL;
  }

  free_Array( stack, TRUE );

  g_seq->path = feat_path;
}
