/*  Last edited: Aug  3 15:53 2002 (klh) */
/**********************************************************************
 ** File: sequence.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#include "sequence.h"

#define ALLOC_STEP 100

/********************************************************************/
/**************** static functions **********************************/
/********************************************************************/


/*********************************************************************
 FUNCTION: convert_gff_line_to_Gaze_entities
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES: Helper to:
   - convert_gff_Gaze_Sequence
   - convert_gff_Gaze_Sequence_list
 *********************************************************************/
static void convert_gff_line_to_Gaze_entities(Gaze_Sequence *g_seq,
					      GFF_line *gff_line,
					      Array *gff2fts) {
  int i, j;

  if ((gff_line->start >= g_seq->seq_region.s) && (gff_line->end <= g_seq->seq_region.e)) {
    /* we have a partial overlap */
    for(i=0; i < gff2fts->len; i++) {
      GFF_to_Gaze_entities *con = index_Array( gff2fts, GFF_to_Gaze_entities *, i);

      if (((con->gff_source == NULL) || strcmp(con->gff_source, gff_line->source) == 0) &&
	  ((con->gff_feature == NULL) || strcmp(con->gff_feature, gff_line->type) == 0) &&
	  ((con->gff_strand == NULL) || strcmp(con->gff_strand, gff_line->strand) == 0) &&
	  ((con->gff_frame == NULL) || strcmp(con->gff_frame, gff_line->frame) == 0)) {
	
	if ((gff_line->start >= g_seq->seq_region.s) && (gff_line->end <= g_seq->seq_region.e)) {
	  /* only features completly within the given region are considered */
	  for(j=0; j < con->features->len; j++) {
	    Gaze_entity *ge = index_Array( con->features, Gaze_entity *, j );
	    Feature *ft = new_Feature();

	    ft->feat_idx = ge->entity_idx;
	    ft->real_pos.s = gff_line->start + ge->offsets.s;
	    ft->real_pos.e = gff_line->end - ge->offsets.e;
	    ft->score = gff_line->score;
	    if (ge->has_score)
	      ft->score = ge->score;
	    
	    if (ft->score < index_Array( g_seq->min_scores, double, ft->feat_idx ))
	      index_Array( g_seq->min_scores, double, ft->feat_idx ) = ft->score;
	    
	    append_val_Array( g_seq->features, ft );
	  }
	}
	
	/* ...whereas all overlapping segments are added */
	for(j=0; j < con->segments->len; j++) {
	  Gaze_entity *ge = index_Array( con->segments, Gaze_entity *, j );
	  Segment *seg = new_Segment();

	  seg->seg_idx = ge->entity_idx;
	  seg->pos.s = gff_line->start + ge->offsets.s;
	  seg->pos.e = gff_line->end - ge->offsets.e;
	  seg->score = gff_line->score;
	  if (ge->has_score)
	    seg->score = ge->score;

	  if (seg->pos.s < g_seq->seq_region.s) {
	    int trimmed = g_seq->seq_region.s - seg->pos.s;
	    double trimmed_score = trimmed * (seg->score / (seg->pos.e - seg->pos.s + 1));
	    seg->pos.s = g_seq->seq_region.s;
	    seg->score -= trimmed_score;
	  }
	  if (seg->pos.e > g_seq->seq_region.e) {
	    int trimmed = seg->pos.e - g_seq->seq_region.s;
	    double trimmed_score = trimmed * (seg->score / (seg->pos.e - seg->pos.s + 1));
	    seg->pos.e = g_seq->seq_region.e;
	    seg->score -= trimmed_score;
	  }
	  if (seg->pos.e <= g_seq->seq_region.e &&
	      seg->pos.s >= g_seq->seq_region.s &&
	      seg->pos.e >= seg->pos.s)
	    append_to_Segment_list( index_Array( g_seq->segment_lists, Segment_list *, seg->seg_idx ),
				    seg );
	  
	  /********************************************************************************/
	  /*  the following was some attempt to get per_base scoring working. It 
	      worked, bit it was not general enough to fit in with the GAZE way 
	      of doing things. In particular, it makes heavy use of the real meaning
	      frame and strand fields of the gff file, which is not how GAZE
	      has done things in the past. Things to point out:
	      
	      1. By convention, I assume that for non unity-width per-base scores, 
	      the last base of the region is the "base" (first base on rev strand), 
	      and the previous n bases are dependency bases. This allows hexamers, 
	      pentamers etc to be given. However, the given "frame" arguement is the 
	      codon position of the *first* base of the region, so we need to convert 
	      this into a codon position of the base of interest. 
	      
	      2. A trick is employed to convert a reverse strand codon-position, which 
	      run 2->1->0 along the length of the sequence when viewing from the forward 
	      strand, to a "virtual forward strand codon-postion", which run 0->1->2 along 
	      the length of the sequence. This is a semantinc convenience that makes the 
	      lookup of the correct list consistent with that of normal segments
	      
	      3. By using the modulo position of the start (or end for reverse strand) of 
	      the segments to determine where in the list they go, they can be retrieved 
	      in a simple fashion in the segment score calculation; this way, the position 
	      (mod 3) gives the index of the zeroth codon position, and the other two codon 
	      positions are obtained by subtraction
	  */
	  
	  /*
	    int cod_pos = atoi( gff_line->frame );
	    int pos = gff_line->end;
	    
	    if ((gff_line->end > g_seq->seq_region.e) || (gff_line->start < g_seq->seq_region.s))
	    continue;
	    
	    correct_list = index_Array( g_seq->segment_lists, 
					Segment_list *, 
					index_Array( con->segments, Segment *, j)->seg_idx );
					  
	    cod_pos = (cod_pos + gff_line->end - gff_line->start) % 3;
					  
	    if (! strcmp(gff_line->strand, "-")) {
		 cod_pos = 3 - cod_pos - 1;
		 pos = gff_line->start;
	    }
	    correct_list->per_base[(pos - cod_pos + 3)%3][pos - g_seq->seq_region.s] = gff_line->score;
	  */
	  /*********************************************************************************/
	}
      }
    }
  }
}



/*********************************************************************
 FUNCTION: get_correct_feature_from_gff_line
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES: Helper to:
   - get_correct_features_Gaze_Sequence
   - get_correct_features_Gaze_Sequence_list
 *********************************************************************/
static boolean get_correct_feature_from_gff_line( Gaze_Sequence *g_seq,
						  GFF_line *gff_line,
						  Array *given_list,
						  Dict *feat_dict,
						  boolean define_paths ) {
  int feat_idx, l, j, k;
  
  boolean no_problem = TRUE;

  if ((feat_idx = dict_lookup( feat_dict, gff_line->type )) >= 0) {	  
    if (gff_line->start >= g_seq->seq_region.s && gff_line->end <= g_seq->seq_region.e) {
      
      Feature *feat = new_Feature();
      
      feat->feat_idx = feat_idx;
      feat->real_pos.s = gff_line->start;
      feat->real_pos.e = gff_line->end;
      
      if (! define_paths) {
	/* 9th field is used to determine what sort of "correct" feature this is :
	   Selected   => this features should be used all paths during the DP
	   Antiselected => this feature should be used in NO paths during the DP
	   Path => this feature should be trated as belonging to the single correct
	   path (along with all other "Path" features given for the sequence)
	   If nothing appears in the 9th field, "Path" is assumed 
	*/
	boolean is_selected = TRUE;
	boolean is_antiselected = FALSE;
	
	if (gff_line->group != NULL) {
	  do {
	    for (l=strlen(gff_line->group); gff_line->group[l] != ';' && l > 0; l--);
	    for (j=l; isspace( (int) gff_line->group[j] ) || gff_line->group[j] == ';' ; j++); 
	    for (k=j; gff_line->group[k] != '\0' && ! isspace( (int) gff_line->group[k] ); k++); 
	    gff_line->group[k] = '\0';
	    if (strcmp( &(gff_line->group[j]), "Selected") == 0)
	      ; /* default - do nothing */
	    else if (strcmp( &(gff_line->group[j]), "Antiselected") == 0) {
	      is_selected = FALSE;
	      is_antiselected = TRUE;
	  }
	    if (l > 0)
	      gff_line->group[l] = '\0';
	  } while (l > 0);
	}
	
	feat->is_selected = is_selected;
	feat->is_antiselected = is_antiselected;
      }
      
      append_val_Array( given_list, feat );
    }
  }	
  else {
    fprintf(stderr, "In path for %s, feature '%s' not recognised\n", 
	    g_seq->seq_name, 
	    gff_line->type );
    no_problem = FALSE;
  }

  return no_problem;
}



/*********************************************************************
 FUNCTION: transfer_correct_features_Gaze_Sequence
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES: Helper to:
   - get_correct_features_Gaze_Sequence
   - get_correct_features_Gaze_Sequence_list
 *********************************************************************/

static boolean transfer_correct_features_Gaze_Sequence( Gaze_Sequence *g_seq,
							Array *this_list,
							Dict *feat_dict,
							boolean define_paths) {
  int i, j;
  boolean no_problem = TRUE;

  if ( this_list->len != 0) {
    
    qsort( this_list->data,
	   this_list->len, 
	   sizeof(Feature *), 
	   &order_features );
    
    for (i=j=0; i < this_list->len; i++) {
      Feature *f1 = index_Array( this_list, Feature *, i);
      boolean match = FALSE;
      
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
	  
	  if (define_paths) {
	    f2->is_correct = TRUE;
	    
	    if (g_seq->path == NULL)
	      g_seq->path = new_Array( sizeof( Feature *), TRUE );
	    append_val_Array( g_seq->path, f2 );
	  }
	  else {
	    f2->is_selected = f1->is_selected;
	    f2->is_antiselected = f1->is_antiselected;
	    
	    /* all selected features must belong the correct path, and furthermore
	       no selected feature can. Otherwise, forward score is meaningless, and
	       correct path could be invalidated (for example) by a forced in-frame
	       stop codon! */
	    
	    if ( g_seq->path != NULL && 
		 ( (f2->is_selected && ! f2->is_correct) ||
		   (f2->is_antiselected && f2->is_correct) )) {
	      
	      fprintf (stderr, 
		       "In sequence %s, inconsistency between selected/antiselected feat and path\n",
		       g_seq->seq_name );
	      no_problem = FALSE;
	    }
	  }
	  
	  free_Feature( f1 );
	  match = TRUE;
	}
      }
      
      if (! match) {
	fprintf (stderr, 
		 "In sequence %s, selected/antselected/correct feature not in list: %s %d %d\n",
		 g_seq->seq_name,
		 index_Array( feat_dict, char *, f1->feat_idx ),
		 f1->real_pos.s, f1->real_pos.e);
	no_problem = FALSE;
	j = 0;
      }
    }
  }
  
  if (define_paths) {
    if (no_problem) {
      /* begin and end features to the path */
      Feature *ft;
      
      if (g_seq->path == NULL)
	g_seq->path = new_Array( sizeof( Feature *), TRUE );
      
      ft = g_seq->end_ft;
      ft->is_correct = TRUE;
      append_val_Array( g_seq->path, ft );
      ft = g_seq->beg_ft;
      ft->is_correct = TRUE;
      prepend_val_Array( g_seq->path, ft );
    }
    else {
      if (g_seq->path != NULL) {
	free_Array( g_seq->path, TRUE );
	g_seq->path = NULL;
      }
    }
  }
  
  return no_problem;
}


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
void free_Gaze_Sequence( Gaze_Sequence *g_seq, boolean free_all ) {
  int i;
  
  if (g_seq != NULL) {
    if (g_seq->seq_name != NULL && free_all)
      free_util( g_seq->seq_name );

    if (g_seq->features != NULL) {
      for(i=0; i < g_seq->features->len; i++)
	free_Feature( index_Array( g_seq->features, Feature *, i));
      free_Array( g_seq->features, TRUE);
      g_seq->features = NULL;
    }
    if (g_seq->segment_lists != NULL) {
      for(i=0; i < g_seq->segment_lists->len; i++)
	free_Segment_list( index_Array( g_seq->segment_lists, Segment_list *, i ));
      free_Array( g_seq->segment_lists, TRUE);
      g_seq->segment_lists = NULL;
    }

    if (g_seq->min_scores != NULL) {
      free_Array( g_seq->min_scores, TRUE );
      g_seq->min_scores = NULL;
    }

    if (g_seq->path != NULL) {
      /* the features themselved belong to the main feature
	 list, so we don't free them */
      free_Array( g_seq->path, TRUE );
      g_seq->path = NULL;
    }

    if (g_seq->gff_file_names != NULL && free_all)
      free_Array( g_seq->gff_file_names, TRUE );

    if (g_seq->selected_file_names != NULL && free_all)
      free_Array( g_seq->selected_file_names, TRUE );

    /* freed elsewhere: dna_seq */ 
    /* freed elsewhere: beg_ft */    
    /* freed elsewhere: end_ft */    

    if (free_all)
      free_util( g_seq );
  }
}



/*********************************************************************
 FUNCTION: initialise_Gaze_Sequence
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void initialise_Gaze_Sequence( Gaze_Sequence *g_seq,
			       Gaze_Structure *gs ) {
  int i;

  g_seq->features = new_Array( sizeof(Feature *), TRUE);

  g_seq->beg_ft = new_Feature();
  g_seq->beg_ft->feat_idx = dict_lookup( gs->feat_dict, "BEGIN" );
  g_seq->beg_ft->real_pos.s = g_seq->seq_region.s;  
  g_seq->beg_ft->real_pos.e = g_seq->seq_region.s;  
  g_seq->beg_ft->is_selected = TRUE;
  append_val_Array( g_seq->features, g_seq->beg_ft );
    
  g_seq->end_ft = new_Feature();
  g_seq->end_ft->feat_idx = dict_lookup( gs->feat_dict, "END" );
  g_seq->end_ft->real_pos.s = g_seq->seq_region.e;  
  g_seq->end_ft->real_pos.e = g_seq->seq_region.e; 
  g_seq->end_ft->is_selected = TRUE;
  append_val_Array( g_seq->features, g_seq->end_ft );

  g_seq->segment_lists = new_Array( sizeof(Segment_list *), TRUE);
  set_size_Array( g_seq->segment_lists, gs->seg_dict->len );

  for(i=0; i < g_seq->segment_lists->len; i++) 
    index_Array( g_seq->segment_lists, Segment_list *, i) = 
      new_Segment_list( g_seq->seq_region.s, g_seq->seq_region.e ); 

  g_seq->min_scores = new_Array( sizeof( double ), TRUE );
  set_size_Array( g_seq->min_scores, gs->feat_dict->len );
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
				  int end ) {

  Gaze_Sequence *g_seq = (Gaze_Sequence *) malloc_util( sizeof( Gaze_Sequence ) );

  g_seq->seq_name = strdup_util( seq_name );

  g_seq->offset_dna = 1;
  g_seq->seq_region.s = sta;
  g_seq->seq_region.e = end;
  g_seq->dna_seq = NULL;
  g_seq->path = NULL;
  g_seq->features = NULL;
  g_seq->segment_lists = NULL;
  g_seq->min_scores = NULL;
  g_seq->beg_ft = NULL;
  g_seq->end_ft = NULL;
  g_seq->dna_file_name = NULL;
  g_seq->gene_file_name = NULL;
  g_seq->gff_file_names = NULL;
  g_seq->selected_file_names = NULL;

  return g_seq;
}


/*********************************************************************
 FUNCTION: convert_dna_Gaze_Sequence
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void convert_dna_Gaze_Sequence ( Gaze_Sequence *g_seq,
				 Array *dna2fts, 
				 Array *offsets, 
				 Dict *motif_dict ) {

  int i,j,k,offset;
  StartEnd *off;

  /* dna_str[i] = residue (dna_off + i) */
  
  if (g_seq->dna_seq == NULL)
    return;

  /* first get the features from the DNA... */

  for (i=0; i < dna2fts->len; i++) {
    DNA_to_Gaze_entities *con = index_Array(dna2fts, DNA_to_Gaze_entities *, i);
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
	Gaze_entity *ge = index_Array( con->features, Gaze_entity *, j );
	Feature *ft = new_Feature();

	ft->feat_idx = ge->entity_idx;
	ft->real_pos.s = start_match + ge->offsets.s;
	ft->real_pos.e = end_match - ge->offsets.e;

	ft->score = ge->has_score ? ge->score : index_Array( g_seq->min_scores, double, ft->feat_idx );

	/* only add the feature if its adjusted position lies within the sequence */
	if (ft->real_pos.s >= g_seq->seq_region.s && ft->real_pos.e <= g_seq->seq_region.e)
	  append_val_Array( g_seq->features, ft );
      }

      for(j=0; j < con->segments->len; j++) {
	Gaze_entity *ge = index_Array( con->segments, Gaze_entity *, j );
	Segment *seg = new_Segment(); 

	seg->seg_idx = ge->entity_idx;
	seg->pos.s = start_match + ge->offsets.s; 
	seg->pos.e = end_match - ge->offsets.e;
	seg->score = ge->has_score ? ge->score : 0.0;

	/* May need to trim back the segment so that it fits inside the sequence */
	if (seg->pos.s < g_seq->seq_region.s) {
	  int trimmed = g_seq->seq_region.s - seg->pos.s;
	  double trimmed_score = trimmed * (seg->score / (seg->pos.e - seg->pos.s + 1));
	  seg->pos.s = g_seq->seq_region.s;
	  seg->score -= trimmed_score;
	}
	if (seg->pos.e > g_seq->seq_region.e) {
	  int trimmed = seg->pos.e - g_seq->seq_region.s;
	  double trimmed_score = trimmed * (seg->score / (seg->pos.e - seg->pos.s + 1));
	  seg->pos.e = g_seq->seq_region.e;
	  seg->score -= trimmed_score;
	}
	if (seg->pos.e <= g_seq->seq_region.e &&
	    seg->pos.s >= g_seq->seq_region.s &&
	    seg->pos.e >= seg->pos.s)
	  append_to_Segment_list( index_Array( g_seq->segment_lists, Segment_list *, seg->seg_idx ),
				  seg );
      }

      *match = save;
      ptr = match + 1; offset++;
    }
  }

  /* and then get the DNA for the features... */

  if (offsets != NULL) {
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
}


/*********************************************************************
 FUNCTION: convert_gff_Gaze_Sequence
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void convert_gff_Gaze_Sequence( Gaze_Sequence *g_seq,
				Array *full_file_list,
				Array *gff2fts ) {
  
  int f_idx;

  GFF_line *gff_line = new_GFF_line();
  Array *useful_files = NULL;

  Array *file_list = g_seq->gff_file_names == NULL ? full_file_list : g_seq->gff_file_names;

  for(f_idx=0; f_idx < file_list->len; f_idx++) {
    char *this_file_name = index_Array( file_list, char *, f_idx );
    FILE *file = fopen (this_file_name, "r" );
    boolean useful_data = FALSE;

    while( read_GFF_line( file, gff_line ) != 0 ) {
      /* First check that we reading annotation for the same sequence that we've
	 doing so all along */
      if (strcmp( g_seq->seq_name, gff_line->seqname))
	continue;
      else {
	useful_data = TRUE;
	convert_gff_line_to_Gaze_entities( g_seq,
					   gff_line,
					   gff2fts );
	
      }
    }
    fclose( file );

    if (useful_data && g_seq->gff_file_names == NULL) {
      /* register this file as useful */
      if (useful_files == NULL)
	useful_files = new_Array( sizeof( char *), TRUE );
      append_val_Array( useful_files, this_file_name );
    }
  }

  if (g_seq->gff_file_names == NULL)
    g_seq->gff_file_names = useful_files;

  free_GFF_line( gff_line );
} 






/*********************************************************************
 FUNCTION: get_correct_feats_Gaze_Sequence
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
   This function converts the given GFF file into a list
   of features that GAZE knows about (i.e. those defined in
   the given feature dictionary). Returns true iff everything went okay
 *********************************************************************/
boolean get_correct_feats_Gaze_Sequence( Gaze_Sequence *g_seq,
					 Array *total_file_list,
					 Dict *feat_dict,
					 boolean define_paths) {

  int f_idx = 0;
  GFF_line *gff_line;
  boolean no_problem = TRUE;
  boolean no_more_files = FALSE;
  
  Array *feat_list = new_Array( sizeof( Feature *), TRUE );
  Array *useful_files = NULL;

  gff_line = new_GFF_line();

  while (! no_more_files) {
    FILE *gff_file;
    char *this_file_name;
    boolean useful_data = FALSE;

    if ( define_paths) {
      this_file_name = g_seq->gene_file_name;

      if (this_file_name == NULL) {
	if (f_idx > total_file_list->len - 1) {
	  no_more_files = TRUE;
	  continue;
	}
	else
	  this_file_name = index_Array( total_file_list, char *, f_idx++ );
      }
      else
	no_more_files = TRUE;
    }
    else {
      if (g_seq->selected_file_names != NULL) {
	if (f_idx > g_seq->selected_file_names->len - 1) {
	  no_more_files = TRUE;
	  continue;
	}
	else
	  this_file_name = index_Array( g_seq->selected_file_names, char *, f_idx++ );
      }
      else {
	if (f_idx > total_file_list->len - 1) {
	  no_more_files = TRUE;
	  continue;
	}
	else
	  this_file_name = index_Array( total_file_list, char *, f_idx++ );
      }
    }

    gff_file = fopen( this_file_name, "r");

    while( read_GFF_line( gff_file, gff_line ) != 0 ){

      if (strcmp( g_seq->seq_name, gff_line->seqname))
	continue;
      else {
	useful_data = TRUE;
	no_problem &= get_correct_feature_from_gff_line( g_seq,
							 gff_line,
							 feat_list,
							 feat_dict,
							 define_paths );
      }      
    }

    fclose( gff_file );

    /* need to register this file as having the path/correct features for this */
    if (useful_data) {
      if (define_paths) {
	if (g_seq->gene_file_name == NULL) {
	  g_seq->gene_file_name = this_file_name;
	  no_more_files = TRUE;
	}
      }
      else {	  
	if (g_seq->selected_file_names == NULL) {
	  if (useful_files == NULL)
	    useful_files = new_Array( sizeof( char *), TRUE );
      
	  append_val_Array( useful_files, this_file_name );
	}
      }
    }
  }

  /* now locate "correct path" features in the complete list */

  if (no_problem )
    no_problem = transfer_correct_features_Gaze_Sequence( g_seq,
							   feat_list,
							   feat_dict,
							   define_paths );

  if (! define_paths && g_seq->selected_file_names == NULL) 
    g_seq->selected_file_names = useful_files;

  /* features themselves will have been freed during the transfer */
  free_Array(feat_list, TRUE );
  free_GFF_line( gff_line );

  return no_problem;

}



/*********************************************************************
 FUNCTION: read_dna_Gaze_Sequence
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
   This routine changes the value of offset_dna for each Gaze_Sequence,
   but since thisis the only routine that makes use of it, this is
   not harmful
 *********************************************************************/
void read_dna_Gaze_Sequence( Gaze_Sequence *g_seq,
			     Array *total_file_list ) {
  
  char *name, c;
  Line *ln = new_Line(); 
  int num_bases = 0;
  int line_len, i, f_idx = 0;
  boolean no_more_files = FALSE;
  int dna_offset = g_seq->offset_dna;

  if (g_seq->dna_seq != NULL) {
    free_util( g_seq->dna_seq );
    g_seq->dna_seq = NULL;
  }

  while (! no_more_files) {
    FILE *dna_file;
    char *this_file_name;

    this_file_name = g_seq->dna_file_name;
    if (this_file_name == NULL) {
      if (f_idx > total_file_list->len - 1) {
	no_more_files = TRUE;
	continue;
      }
      else
	this_file_name = index_Array( total_file_list, char *, f_idx++ );
    }
    else
      no_more_files = TRUE;

    dna_file = fopen( this_file_name, "r");

    while( (line_len = read_Line( dna_file, ln )) != 0) {
      int idx = 0;

      if (ln->buf[idx] == '>') {
	/* skip to first non-white-space character */
	while ( isspace( (int) ln->buf[++idx] ) );
	name = &(ln->buf[idx++]); 
	/* skip to end of name */
	while ( ln->buf[idx] != '\0' && !isspace( (int) ln->buf[++idx] ) );
	ln->buf[idx] = '\0';

	if (strcmp( name, g_seq->seq_name ) == 0) {

	  g_seq->dna_seq = (char *) malloc_util (ALLOC_STEP * sizeof( char ) );
	    
	  if (g_seq->seq_region.s == 0)
	    g_seq->seq_region.s = dna_offset;
	  	  
	  num_bases = 0;
	}
	else {
	  /* we've come across another sequence, so if we've already read our sequence, we're done */
	  if (g_seq->dna_seq != NULL)
	    break;
	}
      }
      else {
	/* this is a DNA line */
	if (g_seq->dna_seq != NULL) {
	  for(i=0; i < line_len; i++) {
	    if (! isspace( (int) ln->buf[i] )) {
	      c = tolower( (int) ln->buf[i]);
	      
	      if ( dna_offset < g_seq->seq_region.s )
		; /* do nothing */
	      else if (g_seq->seq_region.e && (dna_offset > g_seq->seq_region.e)) {
		/* we've gone past the end of the region of interest for this sequence */
		dna_offset++;
		break;
	      }
	      else {
		/* store the base - increase memory if necessary */
		if ( (num_bases % ALLOC_STEP) == 0 ) 
		  if (num_bases != 0)
		    g_seq->dna_seq = (char *) realloc_util( g_seq->dna_seq, 
							    (num_bases + ALLOC_STEP) * sizeof(char) );

		g_seq->dna_seq[num_bases++] = c;
	      }

	      dna_offset++;
	    }
	  }
	}
      }
    }
    fclose( dna_file );

    if (g_seq->dna_seq != NULL && g_seq->dna_file_name == NULL) {
      /* register this file as having the DNA for this seq */
      g_seq->dna_file_name = this_file_name;
      no_more_files = TRUE;
    }
  }

  /* finally, clean up and check */
    
  if (g_seq->dna_seq != NULL) {
    if (g_seq->seq_region.e == 0)
      g_seq->seq_region.e = dna_offset - 1;
      
    num_bases = g_seq->seq_region.e - g_seq->seq_region.s + 1;
    
    g_seq->dna_seq = (char *) realloc_util( g_seq->dna_seq, num_bases + 1 );
    g_seq->dna_seq[num_bases] = '\0';
    
  }
  else {
    /* we have a sequence for which ther was no DNA in any of the files */
    fprintf(stderr, "Warning: no DNA found for %s\n", g_seq->seq_name );
  } 

  free_Line( ln );

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

  /* now condense the array */

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
      free_Gaze_Sequence( list->seq_list[i], TRUE );
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
Gaze_Sequence_list *new_Gaze_Sequence_list ( Array *names ) {
  int i;

  Gaze_Sequence_list *list = (Gaze_Sequence_list *) malloc_util( sizeof(Gaze_Sequence_list));

  list->num_seqs = names->len;
  list->seq_list = (Gaze_Sequence **) malloc_util( list->num_seqs * sizeof( Gaze_Sequence *) );
  
  list->seq_id_dict = new_Array( sizeof( char *), TRUE );
  set_size_Array( list->seq_id_dict, list->num_seqs );

  for(i=0; i < names->len; i++)
    index_Array( list->seq_id_dict, char *, i) = index_Array( names, char *, i);

  return list;
}



/*********************************************************************
 FUNCTION: read_dna_Gaze_Sequence_list
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
   DEPRECATED: does not handle several regions from the same sequence
 *********************************************************************/
void read_dna_Gaze_Sequence_list( Gaze_Sequence_list *glist,
				  Array *file_list ) {
  
  char *name, c;
  Line *ln = new_Line(); 
  boolean interested = FALSE;
  Gaze_Sequence *g_seq = NULL;
  int num_bases = 0;
  int line_len, g_seq_idx, i, f;

  for (f=0; f < file_list->len; f++) {
    FILE *dna_file = fopen( index_Array( file_list, char *, f ), "r");

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
    fclose( dna_file );
  }

  /* finally, clean up and check */
  
  for(i=0; i < glist->num_seqs; i++) {
    g_seq = glist->seq_list[i];
    
    if (g_seq->dna_seq != NULL) {
      if (g_seq->seq_region.e == 0)
	g_seq->seq_region.e = g_seq->offset_dna - 1;
      
      num_bases = g_seq->seq_region.e - g_seq->seq_region.s + 1;
      /* The following is not technically necessary but is done for consistency */
      g_seq->offset_dna = g_seq->seq_region.s;

      g_seq->dna_seq = (char *) realloc_util( g_seq->dna_seq, num_bases + 1 );
      g_seq->dna_seq[num_bases] = '\0';

    }
    else {
      /* we have a sequence for which ther was no DNA in any of the files */
      fprintf(stderr, "Warning: no DNA found for %s\n", g_seq->seq_name );
    } 
  }

  free_Line( ln );

}




/*********************************************************************
 FUNCTION: convert_gff_Gaze_Sequence_list
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void convert_gff_Gaze_Sequence_list( Gaze_Sequence_list *glist,
				     Array *file_list,
				     Array *gff2fts ) {
  
  int f, s;
  Gaze_Sequence *g_seq;

  GFF_line *gff_line = new_GFF_line();
  
  for(f=0; f < file_list->len; f++) {
    FILE *file = fopen (index_Array( file_list, char *, f ), "r" );

    while( read_GFF_line( file, gff_line ) != 0 ) {
      /* First check that we reading annotation for the same sequence that we've
	 doing so all along */

      for (s=0; s < glist->num_seqs; s++) {
	g_seq = glist->seq_list[s];

	if ( ! strcmp (g_seq->seq_name,  gff_line->seqname ))
	  convert_gff_line_to_Gaze_entities( g_seq,
					     gff_line,
					     gff2fts );
	
      }
    }
    fclose( file );
  }

  free_GFF_line( gff_line );
} 




/*********************************************************************
 FUNCTION: get_correct_feats_Gaze_Sequence_list
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
   This function converts the given GFF file into a list
   of features that GAZE knows about (i.e. those defined in
   the given feature dictionary). Returns true iff everything went okay
 *********************************************************************/
boolean get_correct_feats_Gaze_Sequence_list( Gaze_Sequence_list *glist,
					      Array *file_list,
					      Dict *feat_dict,
					      boolean define_paths) {
  
  int i, seq_idx;
  GFF_line *gff_line;
  boolean no_problem = TRUE;

  Array *all_lists = new_Array( sizeof( Array *), TRUE );
  set_size_Array( all_lists, glist->num_seqs );
  for(seq_idx=0; seq_idx < glist->num_seqs; seq_idx++)
    index_Array( all_lists, Array *, seq_idx ) = new_Array( sizeof( Feature * ), TRUE );

  gff_line = new_GFF_line();

  for (i=0; i < file_list->len; i++) {
    FILE *file = fopen( index_Array( file_list, char *, i), "r" );

    while( read_GFF_line( file, gff_line ) != 0 ){
      for(seq_idx=0; seq_idx < glist->num_seqs; seq_idx++) {
	if (! strcmp( glist->seq_list[seq_idx]->seq_name, gff_line->seqname))
	  no_problem &= get_correct_feature_from_gff_line( glist->seq_list[seq_idx], 
							   gff_line,
							   index_Array( all_lists, Array *, seq_idx ),
							   feat_dict,
							   define_paths );
      }	
    }

    fclose( file );
  }

  /* now locate "correct path" features in the complete list */

  if (no_problem ) {

    for(seq_idx=0; seq_idx < glist->num_seqs; seq_idx++) 
      no_problem &= transfer_correct_features_Gaze_Sequence( glist->seq_list[seq_idx],
							     index_Array( all_lists, 
									  Array *, 
									  seq_idx ),
							     feat_dict,
							     define_paths );
  }

  for(i=0; i < all_lists->len; i++)
    free_Array( index_Array( all_lists, Array *, i), TRUE );
  free_Array( all_lists, TRUE );
  free_GFF_line( gff_line );

  return no_problem;

}


/********************************************************************/
/**************** Segment_list **************************************/
/********************************************************************/

/*********************************************************************
 FUNCTION: append_to_Segment_list
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void append_to_Segment_list( Segment_list *sl, Segment *seg ) {
  Segment *sg2;
  
  /* make the segment score per-base */
  seg->score /= (seg->pos.e - seg->pos.s + 1);
  sg2 = clone_Segment( seg );
  
  append_val_Array( index_Array( sl->orig, Array *, seg->pos.s % 3 ), seg );
  append_val_Array( index_Array( sl->orig, Array *, 3 ), sg2 );
 
}

/*********************************************************************
 FUNCTION: free_Segment_list
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Segment_list( Segment_list *sl ) {
  int i, j;

  if (sl != NULL) {
    if (sl->orig != NULL) {
      for(i=0; i < sl->orig->len; i++) {
	Array *segs = index_Array( sl->orig, Array *, i);
	if (segs != NULL) {
	  for ( j=0; j < segs->len; j++)
	    free_Segment( index_Array( segs, Segment *, j));
	  free_Array( segs, TRUE );
	}
      }
      free_Array( sl->orig, TRUE );
    }

    if (sl->proj != NULL) {
      for(i=0; i < sl->proj->len; i++) {
	Array *segs = index_Array( sl->proj, Array *, i);
	if (segs != NULL) {
	  for ( j=0; j < segs->len; j++)
	    free_Segment( index_Array( segs, Segment *, j));
	  free_Array( segs, TRUE );
	}
      }
      free_Array( sl->proj, TRUE );
    }

    for (i=0; i < 3; i++)
      if (sl->per_base[i] != NULL)
	free_util( sl->per_base[i] );

    free_util( sl );
  }
}

/*********************************************************************
 FUNCTION: new_Segment_list
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Segment_list *new_Segment_list( int start_reg, int end_reg ) {
  int i; 

  Segment_list *sl = (Segment_list *) malloc_util( sizeof( Segment_list ) );

  sl->orig = new_Array( sizeof( Array * ), TRUE );
  /* one list for each frame, plus a final list for the complete seg list (all frames) */
  set_size_Array( sl->orig, 4 );
  for (i=0; i < 4; i++)
    index_Array( sl->orig, Array *, i) = new_Array( sizeof( Segment * ), TRUE );

  /* in future, the following will only be performed for segments types that
     require it. For now though, do it by default */

  sl->reg_len = end_reg - start_reg + 1;


  sl->proj = NULL;
  for(i=0; i < 3; i++) {
    /*
    sl->per_base[i] = (double *) malloc_util( sl->reg_len * sizeof( double ) );
    for (j=0; j < sl->reg_len; j++)
      sl->per_base[i][j] = 0.0;
    */
    sl->per_base[i] = NULL;
  }


  return sl;
}

/*********************************************************************
 FUNCTION: index_Segment_list
 DESCRIPTION:
   This function pre-processes the semgnet_list structure to make it 
   more amenable to segment lookup. 
   Firstly, it makes the per_base element (if it exists) cumulative,
   allowing calculation by simple subtraction.

   Secondly, it stores, for each segment, the maximum right-hand position
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
void index_Segment_list(Segment_list *sl) {
  int tp_idx, i, j;

  /* first, cumulativeise the per_base element */
  
  /*
  for (i=0; i < 3; i++) {
    double score_so_far = 0.0;
    
    for(j=0; j < sl->reg_len; j++) {
      sl->per_base[i][j] += score_so_far;
      score_so_far = sl->per_base[i][j];
    }
  }
  */

  /* and now index the normal segments */

  for (tp_idx = 0; tp_idx < 2; tp_idx++) {
    Array *tp = (tp_idx % 2 == 0) ? sl->orig : sl->proj;

    if (tp != NULL) {
      for (i=0; i < tp->len; i++) {
	Array *this_list = index_Array( tp, Array *, i );

	if (this_list != NULL) {
	  int max_end_upstream = -1;
	  int idx_of_max_end_upstream = -1;
	  
	  for (j=0; j < this_list->len; j++) {
	    Segment *this_seg = index_Array( this_list, Segment *, j);
	    
	    if (this_seg->pos.e > max_end_upstream) {
	      max_end_upstream = this_seg->pos.e;
	      idx_of_max_end_upstream = i;
	    }
	    
	    this_seg->max_end_up = max_end_upstream;
	    this_seg->max_end_up_idx = idx_of_max_end_upstream;
	  }    
	}
      }
    }
  }
}


/*********************************************************************
 FUNCTION: project_Segment_list
 DESCRIPTION:
   This function takes a list of possible overlapping scored segments,
   and returns a non-overlapping list of segment that results from
   the projection of the original list onto the sequence.
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void project_Segment_list( Segment_list *sl) {
  boolean found_match;
  Segment *last_seg;
  int i,j, k;

  if (sl->orig == NULL || sl->proj != NULL)
    return;
  
  /* set up the projected segment list */
  sl->proj = new_Array( sizeof( Array * ), TRUE );
  set_size_Array( sl->proj, sl->orig->len );
  for (i=0; i < sl->proj->len; i++)
    index_Array( sl->proj, Array *, i) = NULL;
  
  /* and now populate it */
  for (k=0; k < sl->orig->len; k++) {
    Array *segs = index_Array( sl->orig, Array *, k);
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
    
    index_Array( sl->proj, Array *, k) = new_Array( sizeof( Segment *), TRUE );

    last_seg = NULL;
    for (i=0; i < proj->len; i++) {
      Segment *this_seg = index_Array( proj, Segment *, i );

      if (last_seg == NULL || 
	  this_seg->pos.s - last_seg->pos.e != 1 || 
	  ABS(this_seg->score - last_seg->score) > 1.0e-10) {
	
	/* add the seg */
	append_val_Array(index_Array( sl->proj, Array *, k), this_seg );
	last_seg = this_seg;
      }
      else {
	/* adjust the end of last seg and don't push this one */
	last_seg->pos.e = this_seg->pos.e;
	free_Segment( this_seg );
      }
    }

    free_Array( proj, TRUE );
  }

}



/*********************************************************************
 FUNCTION: scale_Segment_list
 DESCRIPTION:
   This function takes a segment list and a sclaling factor, and
   scales the scores of the various elements of the segment list
   by the scaling factor
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void scale_Segment_list( Segment_list *sl, double scale ) {
  int i, j;

  /* first the per-base element, if it exists */

  /*
  for(i=0; i < 3; i++) 
    if (sl->per_base[i] != NULL)
      for (j=0; j < sl->reg_len; j++)
	sl->per_base[i][j] *= scale;
  */

  /* only need to scale the orig list, and this scaling
     will be propagated when the segments are projected */

  if (sl->orig != NULL) {
    for( i=0; i < sl->orig->len; i++) {      
      Array *list = index_Array( sl->orig, Array *, i );
      
      if (list != NULL) {
	for(j=0; j < list->len; j++) {
	  Segment *seg = index_Array( list, Segment *, j );
	  seg->score *= scale;
	}
      }
    }
  }
}



/*********************************************************************
 FUNCTION: sort_Segment_list
 DESCRIPTION:
   Sorts a segment list by the standard method for sorting segments
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void sort_Segment_list( Segment_list *sl ) {
  int i;

  /* per_base element does not need sorting; sorted by construction */

  if (sl->orig != NULL) {
    for (i=0; i < sl->orig->len; i++) {
      Array *segs = index_Array( sl->orig, Array *, i );

      if (segs != NULL)
	qsort( segs->data, segs->len, sizeof( Segment *), &order_segments );
    }
  }

  if (sl->proj != NULL) {
    for (i=0; i < sl->proj->len; i++) {
      Array *segs = index_Array( sl->proj, Array *, i );

      if (segs != NULL)
	qsort( segs->data, segs->len, sizeof( Segment *), &order_segments );
    }
  }

}
