/*  Last edited: Jul 24 16:28 2002 (klh) */
/**********************************************************************
 ** File: sequence.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#include "sequence.h"

#define ALLOC_STEP 100


/********************************************************************/
/**************** Gaze_Sequence *************************************/
/********************************************************************/


/*********************************************************************
 FUNCTION: free_Gaze_Sequence
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Gaze_Sequence( Gaze_Sequence *g_seq ) {
  int i;
  
  if (g_seq != NULL) {
    if (g_seq->seq_name != NULL)
      free_util( g_seq->seq_name );
    
    if (g_seq->features != NULL) {
      for(i=0; i < g_seq->features->len; i++)
	free_Feature( index_Array( g_seq->features, Feature *, i));
      free_Array( g_seq->features, TRUE);
    }
    if (g_seq->segment_lists != NULL) {
      for(i=0; i < g_seq->segment_lists->len; i++)
	free_Segment_lists( index_Array( g_seq->segment_lists, Segment_lists *, i ));
      free_Array( g_seq->segment_lists, TRUE);
    }

    if (g_seq->min_scores != NULL)
      free_Array( g_seq->min_scores, TRUE );

    if (g_seq->path != NULL) {
      /* the features themselved belong to the main feature
	 list, so we don't free them */
      free_Array( g_seq->path, TRUE );
    }

    /* freed elsewhere: dna_seq */ 
    /* freed elsewhere: beg_dt */    
    /* freed elsewhere: end_ft */    

    free_util( g_seq );
  }
}


/*********************************************************************
 FUNCTION: new_Gaze_Sequence
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Gaze_Sequence *new_Gaze_Sequence( char *seq_name, 
				  int sta, 
				  int end,
				  Gaze_Structure *gs ) {
  int i;

  Gaze_Sequence *g_seq = (Gaze_Sequence *) malloc_util( sizeof( Gaze_Sequence ) );

  g_seq->seq_name = strdup_util( seq_name );
  g_seq->dna_seq = NULL;

  g_seq->features = new_Array( sizeof(Feature *), TRUE);

  g_seq->beg_ft = new_Feature();
  g_seq->beg_ft->feat_idx = dict_lookup( gs->feat_dict, "BEGIN" );
  g_seq->beg_ft->real_pos.s = sta;  /* may be overridden when reading DNA */
  g_seq->beg_ft->real_pos.e = sta;  /* may be overridden when reading DNA */
  append_val_Array( g_seq->features, g_seq->beg_ft );
  
  g_seq->end_ft = new_Feature();
  g_seq->end_ft->feat_idx = dict_lookup( gs->feat_dict, "END" );
  g_seq->end_ft->real_pos.s = end;  /* may be overridden when reading DNA */
  g_seq->end_ft->real_pos.e = end;  /* may be overridden when reading DNA */
  append_val_Array( g_seq->features, g_seq->end_ft );

  g_seq->segment_lists = new_Array( sizeof(Segment_lists *), TRUE);
  set_size_Array( g_seq->segment_lists, gs->seg_dict->len );

  for(i=0; i < g_seq->segment_lists->len; i++) 
    index_Array( g_seq->segment_lists, Segment_lists *, i) = new_Segment_lists(); 

  g_seq->path = NULL;

  g_seq->offset_dna = 1;
  g_seq->seq_region.s = sta;
  g_seq->seq_region.e = end;

  g_seq->min_scores = new_Array( sizeof( double ), TRUE );
  set_size_Array( g_seq->min_scores, gs->feat_dict->len );

  return g_seq;
}



/*********************************************************************
 FUNCTION: get_dna_for_features
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void get_dna_for_features( Gaze_Sequence *g_seq,			  
			   Array *offsets, 
			   Dict *motif_dict ) {

  int i,j,k; 
  StartEnd *off;

  for (i=0; i < g_seq->features->len; i++) {
    Feature *feat = index_Array( g_seq->features, Feature *, i); 

    if ((off = index_Array( offsets, StartEnd *, feat->feat_idx)) != NULL) {
      int dna_start = feat->real_pos.s + off->s;
      int dna_end = feat->real_pos.e - off->e;
      if ( (dna_start >= g_seq->seq_region.s && dna_end <= g_seq->seq_region.e) &&
	   (dna_end - dna_start + 1 > 0) ) {
	char *temp = (char *) malloc_util( (dna_end - dna_start + 2) * sizeof( char ) + 1);
	for(j = dna_start, k=0; j <= dna_end; j++, k++)
	  temp[k] = g_seq->dna_seq[j - g_seq->seq_region.s];
	temp[k] = '\0';

	feat->dna = dict_lookup( motif_dict, temp );
	free_util( temp );
      }
    }
  }
}




/*********************************************************************
 FUNCTION: get_features_from_dna
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void get_features_from_dna( Gaze_Sequence *g_seq,
			    Array *dna2fts ) {
  int i,j,offset;

  /* dna_str[i] = residue (dna_off + i) */
  
  if (g_seq->dna_seq == NULL)
    return;

  for (i=0; i < dna2fts->len; i++) {
    DNA_to_features *con = index_Array(dna2fts, DNA_to_features *, i);
    char *pattern = con->dna_motif;
    int pattern_len = strlen( pattern );
    char *match, *ptr;

    ptr = g_seq->dna_seq;
    offset = 0;
    while ((match = strstr( ptr, pattern )) != NULL) {
      char save = *match;
      int start_match, end_match;

      *match = '\0';
      offset += strlen( ptr );

      start_match = g_seq->seq_region.s + offset;
      end_match = start_match + pattern_len - 1;

      for(j=0; j < con->features->len; j++) {
	Feature *ft = clone_Feature( index_Array( con->features,
						  Feature *,
						  j ));
	ft->real_pos.s = start_match;
	ft->real_pos.e = end_match;

	ft->score = con->has_score ? con->score : index_Array( g_seq->min_scores, double, ft->feat_idx );
	append_val_Array( g_seq->features, ft );
      }

      for(j=0; j < con->segments->len; j++) {
	Segment *sg1, *sg2;
	Segment_lists *correct_list;

	sg1 = clone_Segment( index_Array( con->segments, Segment *, j ));
	sg1->pos.s = start_match;
	sg1->pos.e = end_match;
	sg2 = clone_Segment( sg1 );

	correct_list = index_Array( g_seq->segment_lists, Segment_lists *, sg1->seg_idx );
	append_val_Array( index_Array( correct_list->orig, Array *, start_match % 3 ), sg1 );
	append_val_Array( index_Array( correct_list->orig, Array *, 3 ), sg2 );
      }

      /* and restore */

      *match = save;
      ptr = match + 1; offset++;
    }
  }
}




/*********************************************************************
 FUNCTION: remove_duplicate_features
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
    This function makes the feature list non-redundant by keeping
    only the best-scoring feature of each type at each location
 *********************************************************************/
void remove_duplicate_features( Gaze_Sequence *g_seq) {
  int i,j;

  Feature *last_one = NULL;

  for(i=0; i < g_seq->features->len; i++ ) {
    Feature *this_one = index_Array( g_seq->features, Feature *, i );

    if (last_one != NULL) {
      if (last_one->real_pos.s == this_one->real_pos.s 
	  && last_one->real_pos.e == this_one->real_pos.e
	  && last_one->feat_idx == this_one->feat_idx) {
	
	if (this_one->score > last_one->score) {
	  free_Feature( last_one );
	  index_Array( g_seq->features, Feature *, i-1) = NULL;
	}
	else {
	  free_Feature( this_one );
	  index_Array( g_seq->features, Feature *, i ) = NULL; 
	  this_one = last_one;
	}
      }
    }

    last_one = this_one;
  }

  /* now consense the array */

  for (i=0, j=-1; i < g_seq->features->len && j < g_seq->features->len; i++, j++) {
    if (index_Array( g_seq->features, Feature *, i) == NULL) {
      /* move from the next non-null location to here */
      if (j < i) j = i+1;
      while ( index_Array( g_seq->features, Feature *, j) == NULL ) { j++; }
      index_Array( g_seq->features, Feature *, i ) = index_Array( g_seq->features, Feature *, j );
      index_Array( g_seq->features, Feature *, j ) = NULL;
    }
  }

  g_seq->features->len = i;

}


/********************************************************************/
/**************** Gaze_Sequence_list ********************************/
/********************************************************************/

/*********************************************************************
 FUNCTION: free_Gaze_Sequence_list
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/

void free_Gaze_Sequence_list ( Gaze_Sequence_list *list ) {
  int i;

  if (list != NULL) {
    for (i=0; i < list->num_seqs; i++) {
      free_Gaze_Sequence( list->seq_list[i] );
      free_util( index_Array( list->seq_id_dict, char *, i ) );
    }

    free_Array( list->seq_id_dict, TRUE );
  }
}

/*********************************************************************
 FUNCTION: new_Gaze_Sequence_list
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Gaze_Sequence_list *new_Gaze_Sequence_list ( Array *names,
					     Array *starts,
					     Array *ends,
					     Gaze_Structure *gs ) {
  int i;

  Gaze_Sequence_list *list = (Gaze_Sequence_list *) malloc_util( sizeof(Gaze_Sequence_list));

  list->num_seqs = names->len;
  list->seq_list = (Gaze_Sequence **) malloc_util( list->num_seqs * sizeof( Gaze_Sequence *) );
  
  list->seq_id_dict = new_Array( sizeof( char *), TRUE );
  set_size_Array( list->seq_id_dict, list->num_seqs );
    
  for(i=0; i < names->len; i++) {

    index_Array( list->seq_id_dict, char *, i) = index_Array( names, char *, i);

    list->seq_list[i] =  new_Gaze_Sequence( index_Array( names, char *, i),
					    index_Array( starts, int, i),
					    index_Array( ends, int, i),
					    gs );
    
  }

  return list;
}



/*********************************************************************
 FUNCTION: read_dna_seqs
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
   This routine changes the value of offset_dna for each Gaze_Sequence,
   but since thisis the only routine that makes use of it, this is
   not harmful
 *********************************************************************/
void read_dna_seqs( Gaze_Sequence_list *glist,
		    Array *file_list ) {
		   
  char *name, c;
  Line *ln = new_Line(); 
  boolean interested = FALSE;
  Gaze_Sequence *g_seq;
  int line_len, g_seq_idx, i, f, num_bases;

  for (f=0; f < file_list->len; f++) {
    FILE *dna_file = index_Array( file_list, FILE *, f );

    while( (line_len = read_Line( dna_file, ln )) != 0) {
      int idx = 0;

      if (ln->buf[idx] == '>') {
	/* skip to first non-white-space character */
	while ( isspace( (int) ln->buf[++idx] ) );
	name = &(ln->buf[idx++]); 
	/* skip to end of name */
	while ( ln->buf[idx] != '\0' && !isspace( (int) ln->buf[++idx] ) );
	ln->buf[idx] = '\0';

	/* locate the sequence of interest */
	interested = ((g_seq_idx = dict_lookup( glist->seq_id_dict, name )) >= 0) ? TRUE : FALSE;
	
	if (interested) {
	  g_seq = glist->seq_list[g_seq_idx];

	  if (g_seq->dna_seq != NULL)
	    /* we've seen this dna before - error / ignore */
	    interested = FALSE;
	  else {

	    g_seq->dna_seq = (char *) malloc_util (ALLOC_STEP * sizeof( char ) );
	    
	    if (g_seq->seq_region.s == 0)
	      g_seq->seq_region.s = g_seq->offset_dna;

	  }
	  
	  num_bases = 0;
	}
      }
      else {
	/* this is a DNA line */
	if (interested) {
	  for(i=0; i < line_len; i++) {
	    if (! isspace( (int) ln->buf[i] )) {
	      c = tolower( (int) ln->buf[i]);
	      
	      if ( g_seq->offset_dna < g_seq->seq_region.s ) {

		/* do nothing */
	      }
	      else if (g_seq->seq_region.e && (g_seq->offset_dna > g_seq->seq_region.e)) {
		/* we've gone past the end of the region of interest for this sequence */
		interested = FALSE;
	      }
	      else {
		/* store the base - increase memory if necessary */
		if ( (num_bases % ALLOC_STEP) == 0 ) 
		  if (num_bases != 0)
		    g_seq->dna_seq = (char *) realloc_util( g_seq->dna_seq, 
							    (num_bases + ALLOC_STEP) * sizeof(char) );

		g_seq->dna_seq[num_bases++] = c;
	      }

	      g_seq->offset_dna++;
	    }
	  }
	}
      }
    }
  }

  /* finally, clean up and check */
  
  for(i=0; i < glist->num_seqs; i++) {
    g_seq = glist->seq_list[i];
    
    if (g_seq->dna_seq != NULL) {
      if (g_seq->seq_region.e == 0)
	g_seq->seq_region.e = g_seq->offset_dna - 1;
      
      num_bases = g_seq->seq_region.e - g_seq->seq_region.s + 1;
      
      g_seq->dna_seq = (char *) realloc_util( g_seq->dna_seq, num_bases + 1 );
      g_seq->dna_seq[num_bases] = '\0';

    }
    else {
      /* we have a sequence for which ther was no DNA in any of the files */
      fprintf(stderr, "Warning: no DNA found for %s\n", g_seq->seq_name );
    } 

    /* make the begin and end features consistent with the seq_region */
    g_seq->beg_ft->real_pos.s = g_seq->seq_region.s;
    g_seq->beg_ft->real_pos.e = g_seq->seq_region.s;
    g_seq->end_ft->real_pos.s = g_seq->seq_region.e;
    g_seq->end_ft->real_pos.e = g_seq->seq_region.e;
  }

  free_Line( ln );

}






/*********************************************************************
 FUNCTION: get_features_from_gff
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void get_features_from_gff( Gaze_Sequence_list *glist,
			    Array *file_list,
			    Array *gff2fts,
			    boolean listen_select) {

  int i, j, k, f, g_seq_idx;
  boolean selected, antiselected;
  Gaze_Sequence *g_seq;

  GFF_line *gff_line = new_GFF_line();
  
  for(f=0; f < file_list->len; f++) {
    FILE *file = index_Array( file_list, FILE *, f );

    while( read_GFF_line( file, gff_line ) != 0 ) {
      /* First check that we reading annotation for the same sequence that we've
	 doing so all along */

      if ( (g_seq_idx = dict_lookup( glist->seq_id_dict, gff_line->seqname )) >= 0 )
	g_seq = glist->seq_list[g_seq_idx];
      else
	continue;
      
      selected = antiselected = FALSE;
      
      if (listen_select) {
	if (gff_line->group != NULL) {
	  /* attribute field present */
	  do {
	    for (i=strlen(gff_line->group); gff_line->group[i] != ';' && i > 0; i--);
	    for (j=i; isspace( (int) gff_line->group[j] ) || gff_line->group[j] == ';' ; j++); 
	    for (k=j; gff_line->group[k] != '\0' && ! isspace( (int) gff_line->group[k] ); k++); 
	    gff_line->group[k] = '\0';
	    if (strcmp( &(gff_line->group[j]), "Selected") == 0) {
	      selected = TRUE;
	    }
	    else if (strcmp( &(gff_line->group[j]), "Antiselected") == 0) {
	      antiselected = TRUE;
	    }
	    if (i > 0)
	      gff_line->group[i] = '\0';
	  } while (i > 0);
	}
      }
      
      if ((gff_line->end >= g_seq->seq_region.s) && (gff_line->start <= g_seq->seq_region.e)) {
	/* we have a partial overlap */
	for(i=0; i < gff2fts->len; i++) {
	  GFF_to_features *con = index_Array( gff2fts, GFF_to_features *, i);
	  
	  if (((con->gff_source == NULL) || strcmp(con->gff_source, gff_line->source) == 0) &&
	      ((con->gff_feature == NULL) || strcmp(con->gff_feature, gff_line->type) == 0) &&
	      ((con->gff_strand == NULL) || strcmp(con->gff_strand, gff_line->strand) == 0) &&
	      ((con->gff_frame == NULL) || strcmp(con->gff_frame, gff_line->frame) == 0)) {

	    if ((gff_line->start >= g_seq->seq_region.s) && (gff_line->end <= g_seq->seq_region.e)) {
	      /* only features completly within the given region are considered */
	      for(j=0; j < con->features->len; j++) {
		
		Feature *ft = clone_Feature( index_Array( con->features,
							  Feature *,
							  j ));
		ft->real_pos.s = gff_line->start;
		ft->real_pos.e = gff_line->end;
		ft->score = gff_line->score;
		ft->is_selected = selected;
		ft->is_antiselected = antiselected;

		if (ft->score < index_Array( g_seq->min_scores, double, ft->feat_idx ))
		  index_Array( g_seq->min_scores, double, ft->feat_idx ) = ft->score;

		append_val_Array( g_seq->features, ft );
	      }
	    }
	    
	    /* ...whereas all overlapping segments are added */
	    for(j=0; j < con->segments->len; j++) {
	      Segment *sg1, *sg2;
	      Segment_lists *correct_list;

	      sg1 = clone_Segment( index_Array( con->segments, Segment *, j ));
	      sg1->pos.s = gff_line->start;
	      sg1->pos.e = gff_line->end;
	      sg1->score = gff_line->score;
	      sg2 = clone_Segment( sg1 );
	      
	      correct_list = index_Array( g_seq->segment_lists, Segment_lists *, sg1->seg_idx );
	      append_val_Array( index_Array( correct_list->orig, Array *, sg1->pos.s % 3 ), sg1 );
	      append_val_Array( index_Array( correct_list->orig, Array *, 3 ), sg2 );
	    }
	  }
	}
      }
    }
  }

  free_GFF_line( gff_line );
} 




/*********************************************************************
 FUNCTION: read_in_paths
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
   This function converts the given GFF file into a list
   of features that GAZE knows about (i.e. those defined in
   the given feature dictionary). Returns true iff everything went okay
 *********************************************************************/
boolean read_in_paths( Gaze_Sequence_list *glist,
		       Array *file_list,
		       Dict *feat_dict ) {
  
  int i, j, seq_idx, feat_idx;
  GFF_line *gff_line;
  boolean match;
  boolean no_problem = TRUE;
  Feature *ft;

  /* It is important that the overall feature list and the local path features list 
     are ordered by the same scheme for feature lookup/verification. We therefore sort
     the given feature lists, but they will be sorted again prior to D.P. anyway */

  gff_line = new_GFF_line();

  for (i=0; i < file_list->len; i++) {
    FILE *file = index_Array( file_list, FILE *, i);

    while( read_GFF_line( file, gff_line ) != 0 ){
      
      if ( (seq_idx = dict_lookup( glist->seq_id_dict, gff_line->seqname )) >= 0) {
	Gaze_Sequence *g_seq = glist->seq_list[seq_idx];

	if (g_seq->path == NULL) {
	  g_seq->path = new_Array( sizeof( Feature * ), TRUE );
	  ft = clone_Feature( g_seq->beg_ft );
	  append_val_Array( g_seq->path, ft );
	  ft = clone_Feature( g_seq->end_ft );
	  append_val_Array( g_seq->path, ft );
	}

	if ((feat_idx = dict_lookup( feat_dict, gff_line->type )) >= 0) {
	  if (gff_line->start >= g_seq->seq_region.s && gff_line->end <= g_seq->seq_region.e) {

	    Feature *feat = new_Feature();
	    
	    feat->feat_idx = feat_idx;
	    feat->real_pos.s = gff_line->start;
	    feat->real_pos.e = gff_line->end;
	    /* other fields are irrelevant for now */
	    
	    append_val_Array( g_seq->path, feat );
	  }
	}
	else {
	  fprintf(stderr, "In path for %s, feature '%s' not recognised\n", 
		  g_seq->seq_name, 
		  gff_line->type );
	  no_problem = FALSE;
	}
      }
    }
  }

  /* now locate "correct path" features in the complete list */

  if (no_problem ) {

    for(seq_idx=0; seq_idx < glist->num_seqs; seq_idx++) { 
      Gaze_Sequence *g_seq = glist->seq_list[seq_idx];

      if (g_seq->path != NULL) {

	qsort( g_seq->features->data, 
	       g_seq->features->len, 
	       sizeof(Feature *), 
	       &order_features_standard );
	
	qsort( g_seq->path->data, 
	       g_seq->path->len, 
	       sizeof(Feature *), 
	       &order_features_standard );

	for (i=j=0; i < g_seq->path->len; i++) {
	  Feature *f1 = index_Array( g_seq->path, Feature *, i);
	  match = FALSE;
	  for (; j < g_seq->features->len && ! match; j++) {
	    Feature *f2 = index_Array( g_seq->features, Feature *, j );
	    if ((f1->real_pos.s == f2->real_pos.s) &&  
		(f1->real_pos.e == f2->real_pos.e) && 
		(f1->feat_idx == f2->feat_idx)) {
	      
	      for (j++; j < g_seq->features->len; j++) {
		Feature *f3 = index_Array(g_seq->features, Feature *, j );
		
		if ((f3->feat_idx == f2->feat_idx) && 
		    (f3->real_pos.s == f2->real_pos.s) &&
		    (f3->real_pos.e == f2->real_pos.e)) { 
		  
		  if (f3->score > f2->score)
		    f2 = f3;
		}
		else {
		  j--;
		  break;
		}
	      }
	      f2->is_correct = TRUE;
	      
	      free_Feature( index_Array( g_seq->path, Feature *, i ) );
	      index_Array( g_seq->path, Feature *, i ) = f2;
	      
	      match = TRUE;
	    }
	  }
	  
	  if (! match) {
	    fprintf (stderr, 
		     "In sequence %s, feature in 'correct' path not found in complete list: %s %d %d\n",
		     g_seq->seq_name,
		     index_Array( feat_dict, char *, f1->feat_idx ),
		     f1->real_pos.s, f1->real_pos.e);
	    no_problem = FALSE;
	    j = 0;
	  }
	}
      }
    }
  }

  free_GFF_line( gff_line );

  return no_problem;

}

