/*  Last edited: Aug  3 13:03 2002 (klh) */
/**********************************************************************
 ** File: features.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#include "g_features.h"


/********** Dictionary lookup ***********************/

/*************** StartEnd **************************/

/*********************************************************************
 FUNCTION: free_StartEnd
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_StartEnd( StartEnd *se) {
  if (se != NULL) 
    free_util( se );
}

/*********************************************************************
 FUNCTION: new_StartEnd
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
StartEnd *new_StartEnd( int st, int en) {
  StartEnd *se = (StartEnd *) malloc_util( sizeof( StartEnd ) );
  se->s = st;
  se->e = en;
  return se;
}



/**************** Feature ***************************/

/*********************************************************************
 FUNCTION: clone_Feature
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Feature *clone_Feature(Feature *source) {
  Feature *temp;

  temp = new_Feature();
  temp->feat_idx = source->feat_idx;

  temp->real_pos.s = source->real_pos.s;
  temp->real_pos.e = source->real_pos.e;
  temp->adj_pos.s = source->adj_pos.s;
  temp->adj_pos.e = source->adj_pos.e;

  temp->dna = source->dna;
  temp->score = source->score;
  temp->is_selected = source->is_selected;
  temp->is_antiselected = source->is_antiselected;
  temp->is_correct = source->is_correct;

  return temp;
}


/*********************************************************************
 FUNCTION: free_Feature
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Feature(Feature *f) {
  if (f != NULL) {
    free_util( f );
  }

}



/*********************************************************************
 FUNCTION: new_Feature
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Feature *new_Feature(void) {
  Feature *temp;

  temp = (Feature *) malloc_util( sizeof(Feature) );
  temp->feat_idx = 0;
  temp->real_pos.s = 0;
  temp->real_pos.e = 0;
  temp->adj_pos.s = 0;
  temp->adj_pos.e = 0;
  temp->dna = -1;
  temp->score = 0.0;
  temp->is_selected = FALSE;
  temp->is_antiselected = FALSE;
  temp->is_correct = FALSE;
  temp->path_score = 0.0;
  temp->forward_score = 0.0;
  temp->backward_score = 0.0;
  temp->trace_pointer = 0;
  temp->invalid = FALSE;

#ifdef PARAM
  temp->forward_uval = NULL;
  temp->backward_uval = NULL;
#endif

  return temp;
}



/*********************************************************************
 FUNCTION: write_Feature
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void write_Feature( FILE *fh, Feature *ft, Array *ft_dict, Array *motif_dict) {
  if (ft_dict != NULL) 
    fprintf( fh, "%-10s ", index_Array( ft_dict, char *, ft->feat_idx ));
  else
    fprintf( fh, "%-4d ", ft->feat_idx );    

  fprintf( fh, "%6d %6d %f",	     
	   ft->real_pos.s,
	   ft->real_pos.e,
	   ft->score);
  if (ft->is_selected) 
    fprintf( fh, " Selected");
  if (ft->is_antiselected) 
    fprintf( fh, " Anti-Selected");
  if (ft->is_correct) 
    fprintf( fh, " TRUE");
  if (ft->dna >= 0 && motif_dict != NULL)
    fprintf( fh, " %s", index_Array( motif_dict, char *, ft->dna ) );
  fprintf( fh, "\n");
}



/*********** Segments ***************************/

/*********************************************************************
 FUNCTION: clone_Segment
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Segment *clone_Segment(Segment *source) {
  Segment *temp;

  temp = new_Segment();
  temp->seg_idx = source->seg_idx;

  temp->pos.s = source->pos.s;
  temp->pos.e = source->pos.e;
  temp->score = source->score;
  temp->max_end_up = source->max_end_up;
  temp->max_end_up_idx = source->max_end_up_idx;

  return temp;
}


/*********************************************************************
 FUNCTION: free_Segment
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Segment(Segment *seg) {
  if (seg != NULL)
    free_util( seg );
}


/*********************************************************************
 FUNCTION: new_Segment
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Segment *new_Segment(void) {
  Segment *temp;

  temp = (Segment *) malloc_util( sizeof(Segment) );
  temp->seg_idx = 0;

  temp->pos.s = 0;
  temp->pos.e = 0;

  temp->score = 0.0;
  temp->max_end_up = 0;
  temp->max_end_up_idx = 0;

  return temp;
}



/*********************************************************************
 FUNCTION: write_Segment
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void write_Segment(Segment *seg, FILE *fh, Array *seg_dict) {
  fprintf( fh, "%-10s %6d %6d %f\n",
	   index_Array( seg_dict, char *, seg->seg_idx ),
	   seg->pos.s,
	   seg->pos.e,
	   seg->score * (seg->pos.e - seg->pos.s + 1));
}


/*********************************************************************
 FUNCTION: index_Segments
 DESCRIPTION:
   This function stores, for each segment, the maximum right-hand position
   of all segments to the left (upstream), along with their indices. 
   It is the basis for many list indexing strategies, allowing fast
   look-up of the segments in a given range.  
 RETURNS:
 ARGS: 
 NOTES:
   This procdure is perfectly defined if the segments in the given 
   list are not sorted by their start point. However, there would
   be little point in performing it if that is the case. Hint:
   Sort the segment list by start-point before calling this function.
 *********************************************************************/
void index_Segments(Array *segments) {
  int i;
  int max_end_upstream = -1;
  int idx_of_max_end_upstream = -1;

  for (i=0; i < segments->len; i++) {
    Segment *this_seg = index_Array( segments, Segment *, i);

    if (this_seg->pos.e > max_end_upstream) {
      max_end_upstream = this_seg->pos.e;
      idx_of_max_end_upstream = i;
    }

    this_seg->max_end_up = max_end_upstream;
    this_seg->max_end_up_idx = idx_of_max_end_upstream;
  }
}


/*********************************************************************
 FUNCTION: project_Segments
 DESCRIPTION:
   This function takes a list of possible overlapping scored segments,
   and returns a non-overlapping list of segment that results from
   the projection of the original list onto the sequence.
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Array *project_Segments(Array *segs) {
  boolean found_match;
  Segment *last_seg;
  int i,j;
  Array *ret;
  Array *proj = new_Array( sizeof( Segment * ), TRUE );


  for(i=0; i < segs->len; i++) {
    Segment *current = clone_Segment( index_Array( segs, Segment *, i ) );
    Segment *this = NULL;

    /* conjecture: current either needs to be appended to the list, 
       or its start must lie strictly within one of the segments 
       in the list so far. Method assumes this */
    
    /* all scores are per-residue scores */

    /* 1. locate position in list where the start of current lies */

    found_match = FALSE;
    for(j=proj->len-1; j >= 0; j--) {
      
      this = index_Array( proj, Segment *, j );

      if (current->pos.s <= this->pos.e && current->pos.s >= this->pos.s) {
	/* we've reached the place of interest;  */
	found_match = TRUE;
	break;
      }
      else if (current->pos.s > this->pos.e)
	 break; 
    }
	   
    /* At this point, j is the index of the segment in proj that overlaps with
       current->pos.s. */

    if (! found_match)
      /*  list was empty, or current->pos.s was > proj[last]->pos.e, just append seg */
      append_val_Array( proj, current );
    else {
      /* "this" is overlapping seg and j is the index of that seg */
      if (this->pos.s != current->pos.s) {
	Segment *insert = clone_Segment( current );
	insert->pos.e = this->pos.e;
	this->pos.e = current->pos.s - 1;
	insert->score = this->score;

	insert_val_Array( proj, ++j, insert );
      }
      /* If the current->pos.s matches this->pos.s, We may need to adjust 
	 the score of "this", but that will be done when we go back up the 
	 list looking for current->pos.e */

      /* we need to walk to the end of the segment and do the same there */

      found_match = FALSE;
      for( ; j < proj->len; j++) {
	this = index_Array( proj, Segment *, j );
	if (current->pos.e > this->pos.e) 
	  this->score = MAX( this->score, current->score );
	else if (current->pos.e >= this->pos.s) {
	  if (current->pos.e == this->pos.e) {
	    /* nothing to do except decide on the score */
	    this->score = MAX( current->score, this->score );
	    found_match = TRUE;
	    continue;
	  }
	  else {
	    /* overlap; do the same as last time but the other way around */

	    Segment *insert = clone_Segment( current );
	    insert->pos.s = this->pos.s;
	    this->pos.s = current->pos.e + 1;
	    insert->score = MAX( current->score, this->score );
	    insert_val_Array( proj, j, insert );
	    found_match = TRUE;
	    break;
	  }
	}
	/* the case current->pos.e < this->pos.s should never happen */
      }

      /* if we didn't find a match above, then current extends past the
	 end of the last seg in the list, so we've on more thing to add */

      if (! found_match) {
	current->pos.s = this->pos.e + 1;
	append_val_Array( proj, current );
      }
      else
	free_Segment( current );
    }     
  }

 
  /* last stage: merge adjacent segments with same score (within limits) */

  ret = new_Array( sizeof( Segment * ), TRUE );

  last_seg = NULL;
  for (i=0; i < proj->len; i++) {
    Segment *this_seg = index_Array( proj, Segment *, i );

    if (last_seg == NULL || 
	this_seg->pos.s - last_seg->pos.e != 1 || 
	ABS(this_seg->score - last_seg->score) > 1.0e-10) {

      /* add the seg */
      append_val_Array( ret, this_seg );
      last_seg = this_seg;
    }
    else {
	/* adjust the end of last seg and don't push this one */
	last_seg->pos.e = this_seg->pos.e;
	free_Segment( this_seg );
    }
  }
  free_Array( proj, TRUE );

  return ret;
}



/*************** DNA to features and segments *********************/

/*********************************************************************
 FUNCTION: free_DNA_to_features
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_DNA_to_features(DNA_to_features *dna2fts) {
  int i;

  if (dna2fts != NULL) {
    if (dna2fts->dna_motif != NULL) {
      free_util( dna2fts->dna_motif );
    }
    
    if (dna2fts->features) {
      for(i=0; i < dna2fts->features->len; i++) {
	free_util( index_Array( dna2fts->features, Feature *, i));
      }
      free_Array( dna2fts->features, TRUE );
    }
    if (dna2fts->segments) {
      for(i=0; i < dna2fts->segments->len; i++) {
	  free_util( index_Array( dna2fts->segments, Segment *, i));
      }
      free_Array( dna2fts->segments, TRUE );
    }

    free_util( dna2fts );
  }
}


/*********************************************************************
 FUNCTION: new_DNA_to_features
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
DNA_to_features *new_DNA_to_features(void) {
  DNA_to_features *temp;

  temp = (DNA_to_features *) malloc_util( sizeof(DNA_to_features) );
  temp->dna_motif = NULL;
  temp->has_score = FALSE;
  temp->score = 0.0;
  temp->features = new_Array( sizeof(Feature *), TRUE);
  temp->segments = new_Array( sizeof(Segment *), TRUE);

  return temp;
}



/************* GFF to features and segments ***********************/

/*********************************************************************
 FUNCTION: free_GFF_to_features
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_GFF_to_features(GFF_to_features *gff2fts) {
  int i;

  if (gff2fts != NULL) {
    if (gff2fts->features) {
      for(i=0; i < gff2fts->features->len; i++) {
	free_util( index_Array( gff2fts->features, Feature*, i));
      }
      free_Array( gff2fts->features, TRUE );
    }
    if (gff2fts->segments) {
      for(i=0; i < gff2fts->segments->len; i++) {
	  free_util( index_Array( gff2fts->segments, Segment*, i));
      }
      free_Array( gff2fts->segments, TRUE );
    }
    if (gff2fts->gff_source != NULL) {
      free_util( gff2fts->gff_source );
    }
    if (gff2fts->gff_feature != NULL) {
      free_util( gff2fts->gff_feature );
    }
    if (gff2fts->gff_strand != NULL) {
      free_util( gff2fts->gff_strand );
    }
    if (gff2fts->gff_frame != NULL) {
      free_util( gff2fts->gff_frame );
    }

    free_util( gff2fts );
  }
}


/*********************************************************************
 FUNCTION: new_GFF_to_features
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
GFF_to_features *new_GFF_to_features(void) {
  GFF_to_features *temp;

  temp = (GFF_to_features *) malloc_util( sizeof(GFF_to_features) );
  temp->gff_source = NULL;
  temp->gff_feature = NULL;
  temp->gff_strand = NULL;
  temp->gff_frame = NULL;
  temp->features = new_Array( sizeof(Feature *), TRUE);
  temp->segments = new_Array( sizeof(Segment *), TRUE);

  return temp;
}



/**************** Sorting features and segments **********************/

/*********************************************************************
 FUNCTION: order_features_for_dp
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
   This function is deprecated. It was used in an older version of
   GAZE which tried to form a definitive total order over features
   such that killers only fall in an interval if they really do 
   kill that region. For various reasons, this was unsatisfactory,
   so this problem is now dealt with in the D.P. itself. 
 *********************************************************************/
int order_features_for_dp(const void *a, const void *b) {
  Feature *fa = (Feature *) *((Feature **)a);
  Feature *fb = (Feature *) *((Feature **)b);
  int ret = 0;

  if (fa->real_pos.s > fb->real_pos.e) {
    ret = 1;
  }
  else if (fb->real_pos.s > fa->real_pos.e) {
    ret = -1;
  }
  else {
    /* It is critical what we do with overlapping features here.
       Killer features cut out illegal paths, but this relies on them
       being in the correct position with respect to other features.
       The way I have handled this is slightly hacky. If there is what
       I judge to be a "dangerous" overlap between a pair of features
       (if either their starts, ends match, or they lie at the same
       "position"), then they are ordered according to their type.
    */

    if (fa->real_pos.s == fb->real_pos.s 
	|| fa->real_pos.e == fb->real_pos.e) {
      
      ret = fa->feat_idx - fb->feat_idx;
    }
    else {
      /* both the ends and starts are different */
      
      if (fa->real_pos.s < fb->real_pos.s) {
	ret = -1;

	if ( (fa->adj_pos.s == fb->adj_pos.e) && (fa->feat_idx != fb->feat_idx))
	  ret = fa->feat_idx - fb->feat_idx;
      }
      else {
	ret = 1;

	if ( (fb->adj_pos.s == fa->adj_pos.e) && (fa->feat_idx != fb->feat_idx) )
	  ret = fa->feat_idx - fb->feat_idx;
      }

    }
  }

  return ret;
}


/*********************************************************************
 FUNCTION: order_features
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
int order_features(const void *a, const void *b) {
  Feature *fa = (Feature *) *((Feature **)a);
  Feature *fb = (Feature *) *((Feature **)b);
  int ret = 0;
 
  if (! (ret = fa->real_pos.s - fb->real_pos.s))
    if (! (ret = fa->real_pos.e - fb->real_pos.e))
      ret = fa->feat_idx - fb->feat_idx;
    
  return ret;
}


/*********************************************************************
 FUNCTION: order_segments
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
int order_segments(const void *a, const void *b) {
  Segment *sa = (Segment *) *((Segment **)a);
  Segment *sb = (Segment *) *((Segment **)b);
  int ret = 0;

  if (! (ret = sa->pos.s - sb->pos.s))
    if (! (ret = sa->pos.e - sb->pos.e))
      ret = sa->seg_idx - sb->seg_idx;

  return ret;
}

