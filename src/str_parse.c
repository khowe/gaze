/*  Last edited: Jul 24 14:49 2002 (klh) */
/**********************************************************************
 ** File: str_parse.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/


#include "str_parse.h"


struct Parse_context {
  Array *tag_stack;
  boolean finished_parsing;
  int error;
  XML_Parser the_parser;
  int current_idx;
  int current_idx_2;  /* sometimes require two levels of context */
  Gaze_Structure *gs;
};


/* expat call-back routines */

static void start_tag_structure(void *, const char *, const char **);
static void end_tag_structure(void *, const char *);

/* parsing routines */

static void parse_Gaze_Structure_declarations( struct Parse_context *, const char **);
static void parse_Gaze_Structure_dna2gaze( struct Parse_context *, const char **);
static void parse_Gaze_Structure_dnafeat( struct Parse_context *, const char **);
static void parse_Gaze_Structure_feat( struct Parse_context *, const char **);
static void parse_Gaze_Structure_feature( struct Parse_context *, const char **);
static void parse_Gaze_Structure_model( struct Parse_context *, const char **);
static void parse_Gaze_Structure_gaze( struct Parse_context *, const char **);
static void parse_Gaze_Structure_gff2gaze( struct Parse_context *, const char **);
static void parse_Gaze_Structure_gffline( struct Parse_context *, const char **);
static void parse_Gaze_Structure_killdna( struct Parse_context *, const char **);
static void parse_Gaze_Structure_killfeat( struct Parse_context *, const char **);
static void parse_Gaze_Structure_lengthfunc( struct Parse_context *, const char **);
static void parse_Gaze_Structure_lengthfunction( struct Parse_context *,const char **);
static void parse_Gaze_Structure_lengthfunctions( struct Parse_context *, const char **);
static void parse_Gaze_Structure_output( struct Parse_context *, const char **);
static void parse_Gaze_Structure_point( struct Parse_context *, const char **);
static void parse_Gaze_Structure_seg( struct Parse_context *, const char **);
static void parse_Gaze_Structure_segment( struct Parse_context *, const char **);
static void parse_Gaze_Structure_source( struct Parse_context *, const char **);
static void parse_Gaze_Structure_takedna( struct Parse_context *, const char **);
static void parse_Gaze_Structure_target( struct Parse_context *, const char **);
static void parse_Gaze_Structure_useseg( struct Parse_context *, const char **);


/* Indices into the array of XML tags. Nesting partially
   reflects structure */

enum Tags { 
  GAZE,
    DECLARATIONS,
      FEATURE,
      SEGMENT,
      LENGTHFUNCTION,
    GFF2GAZE,
      GFFLINE,
        FEAT,
        SEG,
    DNA2GAZE,
      DNAFEAT,
      TAKEDNA,
    MODEL,
      TARGET,
        USESEG,
        KILLDNA,
        KILLFEAT,
        OUTPUT,
        SOURCE,
    LENGTHFUNCTIONS,
      LENGTHFUNC,
        POINT, 
  MAXTAG
};

		 
/* XML tags used in Gaze stucture files. Order should correspond
   with that above */

static struct {
  char *tag;
  void (*func)( struct Parse_context *, const char **);
} tag_map[] = {

  { "gaze", &parse_Gaze_Structure_gaze },
  { "declarations", &parse_Gaze_Structure_declarations},
  { "feature", &parse_Gaze_Structure_feature },
  { "segment", &parse_Gaze_Structure_segment },
  { "lengthfunction", &parse_Gaze_Structure_lengthfunction },
  { "gff2gaze", &parse_Gaze_Structure_gff2gaze },
  { "gffline", &parse_Gaze_Structure_gffline },
  { "feat", &parse_Gaze_Structure_feat },
  { "seg", &parse_Gaze_Structure_seg },
  { "dna2gaze", &parse_Gaze_Structure_dna2gaze },
  { "dnafeat", &parse_Gaze_Structure_dnafeat },
  { "takedna", &parse_Gaze_Structure_takedna },
  { "model", &parse_Gaze_Structure_model },
  { "target", &parse_Gaze_Structure_target },
  { "useseg", &parse_Gaze_Structure_useseg },
  { "killdna", &parse_Gaze_Structure_killdna },
  { "killfeat", &parse_Gaze_Structure_killfeat },
  { "output", &parse_Gaze_Structure_output },
  { "source", &parse_Gaze_Structure_source },
  { "lengthfunctions", &parse_Gaze_Structure_lengthfunctions },
  { "lengthfunc", &parse_Gaze_Structure_lengthfunc },
  { "point", &parse_Gaze_Structure_point }
};


#define PARSE_BUFFER_SIZE 8000





/********************************************************************/
/**************** expat callback  routines **************************/
/********************************************************************/


/*********************************************************************
 FUNCTION: start_tag_structure
 DESCRIPTION:
   This function is called whenver an opening XML tag is parsed
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void start_tag_structure(void *data, const char *el, const char **attr) {
  int i;
  char *t;
  boolean match = FALSE;

  struct Parse_context *state = (struct Parse_context *) data;

  if (! state->error) {
    if (state->finished_parsing) {
      /* extra data we are not expecting. Deal with it later */
    }
    else {
      for (i=0; i < MAXTAG && ! match; i++) {
	if (! strcmp( el, tag_map[i].tag )) {
	  match = TRUE;
	  break;
	}
      }
      
      if (! match) {
	fprintf(stderr, "Unrecognised tag: %s\n", el);
	state->error =  XML_GetCurrentLineNumber( state->the_parser );
	/* raise an error some how and quit */      
      }
      else {
	(*tag_map[i].func)(state, attr);
	t =  strdup_util( el );
	append_val_Array(state->tag_stack, t);
      }
    }
  } 
}



/*********************************************************************
 FUNCTION: end_tag_structure
 DESCRIPTION:
   This function is called whenever an XML closing tag is parsed
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void end_tag_structure(void *data, const char *el) {
  struct Parse_context *state = (struct Parse_context *) data;

  if (! state->error) {
    if (state->finished_parsing) {
      /* extra data we are not expecting. Deal with it later */
    }
    else {
      if (state->tag_stack->len) {
	int index =  state->tag_stack->len - 1;
	char *temp = index_Array( state->tag_stack, char *, index );
	remove_index_Array( state->tag_stack, index );
	free_util( temp );
	
	if (! state->tag_stack->len) 
	  state->finished_parsing = TRUE;
      }
    }
  }    
}  /* End of end handler */


/********************************************************************/
/********************** parsing  routines ***************************/
/********************************************************************/




/*********************************************************************
 FUNCTION: parse_Gaze_Structure_declarations
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
*********************************************************************/
static void parse_Gaze_Structure_declarations( struct Parse_context *state, 
					       const char **attr ) {


  if (! state->tag_stack->len ||
      strcmp(index_Array(state->tag_stack, 
			   char *, 
			   state->tag_stack->len - 1), 
	     tag_map[GAZE].tag)) {
    /* flag error and return */
    fprintf(stderr, "Tag 'declarations' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}



/*********************************************************************
 FUNCTION: parse_Gaze_Structure_dna2gaze
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
*********************************************************************/
static void parse_Gaze_Structure_dna2gaze( struct Parse_context *state, 
					   const char **attr ) {

  if (! state->tag_stack->len ||
      strcmp(index_Array(state->tag_stack, 
			   char *, 
			   state->tag_stack->len - 1), 
	     tag_map[GAZE].tag)) {
    /* top of tag stack should be 'gaze' */
    fprintf(stderr, "Tag 'dna2gaze' not expected in this context\n");
    state->error =  XML_GetCurrentLineNumber( state->the_parser );
  }
}




/*********************************************************************
 FUNCTION: parse_Gaze_Structure_dnafeat
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_dnafeat( struct Parse_context *state, 
					  const char **attr ) {

  char *pattern = NULL;
  boolean has_score = FALSE;
  double score = 0.0;
  int i;
  DNA_to_features *new_dna2ft;
  
  if (state->tag_stack->len && 
      ! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[DNA2GAZE].tag)) {
    
    for( i=0; attr[i]; i += 2 ) {
      if (! strcmp( attr[i], "pattern" )) {
	pattern = strdup_util( attr[i+1] );
      }
      else if (! strcmp( attr[i], "score" )) {
	has_score = TRUE;
	score = atof( attr[i+1] );
      }
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'dnafeat', attr '%s' not recognised\n", attr[i]); 
	state->error =  XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    if (! state->error) {
      if (pattern != NULL) {
      	new_dna2ft = new_DNA_to_features();
	new_dna2ft->dna_motif = pattern;
	new_dna2ft->has_score = has_score;
	new_dna2ft->score = score;
	append_val_Array( state->gs->dna_to_feats, new_dna2ft );
      }
      else {
	/* no pattern attribute - bad news */
	fprintf(stderr, "In tag 'dnafeat' attr 'pattern' is not defined\n");
	state->error =  XML_GetCurrentLineNumber( state->the_parser );
      }
    }
  }
  else {
    fprintf(stderr, "Tag 'dnafeat' not expected in this context\n");
    state->error =  XML_GetCurrentLineNumber( state->the_parser );
  }
}




/*********************************************************************
 FUNCTION: parse_Gaze_Structure_feat
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_feat( struct Parse_context *state, 
				       const char **attr ) {
  int i, index = 0;
  Feature *feat;
  double score = 0.0;

  for( i=0; attr[i]; i += 2 ) {
    if (! strcmp( attr[i], "id" )) {
      if ((index = dict_lookup( state->gs->feat_dict, attr[i+1]))< 0) {
	/* undeclared feature */
	fprintf(stderr, "In tag 'feat' attr 'id' has illegal value\n");
	state->error =  XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    else if (! strcmp( attr[i], "score")) {
      score = atof( attr[i+1] );
    }
    else {
      /* unrecognised attribute */
      fprintf(stderr, "In tag 'feat', attr '%s' not recognised\n", attr[i]); 
      state->error =  XML_GetCurrentLineNumber( state->the_parser );
    }
  }  
  
  if (! state->error) {
    if (state->tag_stack->len) {
      if (! strcmp(index_Array(state->tag_stack, 
				 char *, 
				 state->tag_stack->len - 1), 
		   tag_map[GFFLINE].tag)) {
	
	GFF_to_features *gff2fts;
	gff2fts = index_Array( state->gs->gff_to_feats, 
				 GFF_to_features *,
				 state->gs->gff_to_feats->len - 1);
	feat = new_Feature();
	feat->feat_idx = index;
	append_val_Array( gff2fts->features, feat );
      }
      else if (! strcmp(index_Array(state->tag_stack, 
				      char *, 
				      state->tag_stack->len - 1), 
			tag_map[DNAFEAT].tag)) {
	
	DNA_to_features *dna2fts;
	dna2fts = index_Array( state->gs->dna_to_feats, 
				 DNA_to_features *,
				 state->gs->dna_to_feats->len - 1);
	feat = new_Feature();
	feat->feat_idx = index;
	feat->score = score;
	append_val_Array( dna2fts->features, feat );	
      }
      else {
	/* Tag out of context */
	state->error =  XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    else {
      /* Tag out of context */
      fprintf(stderr, "Tag 'feat' not expected in this context\n");
      state->error =  XML_GetCurrentLineNumber( state->the_parser );
    }
  }
}




/*********************************************************************
 FUNCTION: parse_Gaze_Structure_feature
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
*********************************************************************/
static void parse_Gaze_Structure_feature( struct Parse_context *state, 
					  const char **attr ) {
  char *id = NULL;
  int start_offset = 0;
  int end_offset = 0;
  double mul = 1.0;
  int i;
  Feature_Info *new_feat;
  
  if (state->tag_stack->len && 
      ! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[DECLARATIONS].tag)) {
    /* top of tag stack should be 'declarations' */
    
    for( i=0; attr[i]; i += 2 ) {
      if (! strcmp( attr[i], "id" )) {
	/* check this is not a re-defintion */
	if ( dict_lookup(state->gs->feat_dict, attr[i+1]) >= 0) {
	  fprintf(stderr, "In tag 'feature' attr 'id' has already been seen\n");
	  state->error =  XML_GetCurrentLineNumber( state->the_parser );
	}
	else {
	  id = strdup_util( attr[i+1] );
	}
      }
      else if (! strcmp( attr[i], "st_off" )) {
	start_offset = atoi( attr[i+1] ); 
      }
      else if (! strcmp( attr[i], "en_off" )) {
	end_offset = atoi( attr[i+1] );
      }
      else if (! strcmp( attr[i], "mul" )) {
	mul = atof( attr[i+1] ); 
      }
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'feature', attr '%s' not recognised\n", attr[i]); 
	state->error =  XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    if (! state->error) {
      if (id != NULL) { 
      	new_feat = new_Feature_Info( start_offset, end_offset, mul );
	/* we have to place this into the dictionary at the penultimate
	   position, just before the END feature. This is so that the
	   features can be sorted correctly */

	insert_val_Array( state->gs->feat_dict, state->gs->feat_dict->len - 1, id );
	insert_val_Array( state->gs->feat_info, state->gs->feat_info->len - 1, new_feat );
      }
      else {
	/* no id attribute - bad news */
	fprintf(stderr, "In tag 'feature' attr 'id' is not defined\n");
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
  }
  else {
    fprintf(stderr, "Tag 'feature' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}






/*********************************************************************
 FUNCTION: parse_Gaze_Structure_model
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
*********************************************************************/
static void parse_Gaze_Structure_model( struct Parse_context *state, 
					const char **attr ) {

  if (! state->tag_stack->len ||
      strcmp(index_Array(state->tag_stack, 
			   char *, 
			   state->tag_stack->len - 1), 
	     tag_map[GAZE].tag)) {
    fprintf(stderr, "Tag 'model' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}




/*********************************************************************
 FUNCTION: parse_Gaze_Structure_gaze
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
*********************************************************************/
static void parse_Gaze_Structure_gaze( struct Parse_context *state, 
				       const char **attr ) {

  if (state->tag_stack->len ) {
    /* tag stack should be empty */
    fprintf(stderr, "Tag 'gaze' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}



/*********************************************************************
 FUNCTION: parse_Gaze_Structure_gff2gaze
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
*********************************************************************/
static void parse_Gaze_Structure_gff2gaze( struct Parse_context *state, 
					   const char **attr ) {

  if (! state->tag_stack->len ||
      strcmp(index_Array(state->tag_stack, 
			   char *, 
			   state->tag_stack->len - 1), 
	     tag_map[GAZE].tag)) {
    fprintf(stderr, "Tag 'gff2gaze' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}



/*********************************************************************
 FUNCTION: parse_Gaze_Structure_gffline
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_gffline( struct Parse_context *state, 
					  const char **attr ) {

  char *feature = NULL;
  char *source = NULL;
  char *strand = NULL;
  char *frame = NULL;
  int i;


  if (state->tag_stack->len && 
      ((! strcmp(index_Array(state->tag_stack, 
			       char *, 
			       state->tag_stack->len - 1), 
		 tag_map[GFF2GAZE].tag)))) {
    
    for( i=0; attr[i]; i += 2 ) {
      if (! strcmp( attr[i], "feature" )) {
	feature = strdup_util( attr[i+1] );
      }
      else if (! strcmp( attr[i], "source" )) {
	source = strdup_util( attr[i+1] );
      }
      else if (! strcmp( attr[i], "strand" )) {
	strand = strdup_util( attr[i+1] );
      }
      else if (! strcmp( attr[i], "frame" )) {
	frame = strdup_util( attr[i+1] );
      }
      else {
	/* unrecognised attribute */
	state->error = XML_GetCurrentLineNumber( state->the_parser );
	fprintf(stderr, "In tag 'gffline', attr '%s' not recognised\n", attr[i]); 
      }
    }
    if (! state->error) {
      if (feature != NULL || source != NULL || strand != NULL) { 	
	GFF_to_features *new_gff2ft = new_GFF_to_features();
	new_gff2ft->gff_feature = feature;
	new_gff2ft->gff_source = source;
	new_gff2ft->gff_strand = strand;
	new_gff2ft->gff_frame = frame;

	append_val_Array( state->gs->gff_to_feats, new_gff2ft );
      }
      else {
	  fprintf(stderr, "In tag 'gffline' you must define some attributes!\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
  }
  else {
    fprintf(stderr, "Tag 'gffline' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}



/*********************************************************************
 FUNCTION: parse_Gaze_Structure_killdna
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_killdna( struct Parse_context *state, 
					  const char **attr ) {

  int i;
  int j;
  int src_dna_idx = -1;
  int tgt_dna_idx = -1;
  char *new_motif = NULL;

  if (state->tag_stack->len && 
      ! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[SOURCE].tag)) {

    for( i=0; attr[i] && ! state->error; i += 2 ) {
      if (! strcmp( attr[i], "source_dna" )) {
	if (src_dna_idx >= 0) {
	  /* more than one src - error */
	  fprintf(stderr, "In tag 'killdna', attribute 'source_dna' cannot be defined more than once\n"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
	else {
	  new_motif = strdup_util( attr[i+1] );
	  for(j=0; new_motif[j] != '\0'; j++)
	    new_motif[j] = tolower( (int) new_motif[j]);
	  if ((src_dna_idx = dict_lookup(state->gs->motif_dict, attr[i+1])) < 0) {
	    append_val_Array( state->gs->motif_dict, new_motif );
	    src_dna_idx = state->gs->motif_dict->len - 1;
	  }
	  else 
	    free_util( new_motif );
	}
      }
      else if (! strcmp( attr[i], "target_dna" )) {
	if (tgt_dna_idx >= 0) {
	  /* error */
	  fprintf(stderr, "In tag 'killdna', attribute 'target_dna' cannot be defined more than once\n"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
	else {
	  new_motif = strdup_util( attr[i+1] );
	  for(j=0; new_motif[j] != '\0'; j++)
	    new_motif[j] = tolower( (int) new_motif[j]);
	  if ((tgt_dna_idx = dict_lookup(state->gs->motif_dict, attr[i+1])) < 0) {
	    append_val_Array( state->gs->motif_dict, new_motif );
	    tgt_dna_idx = state->gs->motif_dict->len - 1;
	  }
	  else 
	    free_util( new_motif );
	}
      }
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'killdna', attr '%s' not recognised\n", attr[i]); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }


    if (! state->error) {
      Feature_Info *tgt_info = index_Array( state->gs->feat_info, Feature_Info *, state->current_idx);
      Feature_Relation *ft_rel = index_Array( tgt_info->sources, Feature_Relation *, state->current_idx_2);
      Killer_DNA_Qualifier *kdq = new_Killer_DNA_Qualifier();

      kdq->src_dna = src_dna_idx;
      kdq->tgt_dna = tgt_dna_idx;

      if (ft_rel->kill_dna_quals == NULL) 
	ft_rel->kill_dna_quals = new_Array( sizeof( Killer_DNA_Qualifier *), TRUE );
      
      append_val_Array( ft_rel->kill_dna_quals, kdq );
    }
  }
  else {
    /* flag error and return */
    fprintf(stderr, "Tag 'killdna' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}



/*********************************************************************
 FUNCTION: parse_Gaze_Structure_killfeat
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_killfeat( struct Parse_context *state, 
					 const char **attr ) {

  int i;
  boolean in_source_tag = FALSE;
  
  if (state->tag_stack->len && 
      ((! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[TARGET].tag)) ||
      (! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[SOURCE].tag)))) {
    /* top of tag stack should be 'target' or 'source' */

    Killer_Feature_Qualifier *kq = new_Killer_Feature_Qualifier();

    in_source_tag =  strcmp(index_Array(state->tag_stack, 
					  char *, 
					  state->tag_stack->len - 1), 
			    tag_map[TARGET].tag);

    for( i=0; attr[i] && ! state->error; i += 2 ) {
      if (! strcmp( attr[i], "id" )) {
	if ((kq->feat_idx = dict_lookup(state->gs->feat_dict, attr[i+1])) < 0) {
	  fprintf(stderr, "In tag 'killfeat' attr 'id' has an illegal value\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "source_phase" )) {
	if (kq->has_tgt_phase) {
	  fprintf(stderr, "In tag 'killfeat' attrs 'target_phase' and 'source_phase' cannot both be defined\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
	else {
	  kq->has_src_phase = TRUE;
	  kq->phase = atoi( attr[i+1] );
	  if (kq->phase < 0 || kq->phase > 2) {
	    fprintf(stderr, "In tag 'killfeat' attr 'source_phase' has an illegal value\n");
	    state->error = XML_GetCurrentLineNumber( state->the_parser );
	  }
	}
      }
      else if (! strcmp( attr[i], "target_phase" )) {
	if (kq->has_src_phase) {
	  fprintf(stderr, "In tag 'killfeat' attrs 'target_phase' and 'source_phase' cannot both be defined\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
	else {
	  kq->has_tgt_phase = TRUE;
	  kq->phase = atoi( attr[i+1] );
	  if (kq->phase < 0 || kq->phase > 2) {
	    fprintf(stderr, "In tag 'killfeat' attr 'target_phase' has an illegal value\n");
	    state->error = XML_GetCurrentLineNumber( state->the_parser );
	  }
	}
      }
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'killfeat', attr '%s' not recognised\n", attr[i]); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }

    if (! state->error) {
      Feature_Info *target = index_Array( state->gs->feat_info, Feature_Info *, state->current_idx);

      /* if we are in the 'target' tag, we need to update the general kill_off_up
	 array (the one that applies to all sources for the given target)
	 If not, we update the Feature_Relation for the source-target directly */

      if (!in_source_tag) {
	/* must be in target tag - the killer is global to all sources for this target */
	if (target->kill_feat_quals == NULL) {
	  target->kill_feat_quals = new_Array( sizeof( Killer_Feature_Qualifier * ), TRUE);
	  /* set_size_Array(target->kill_feat_quals, state->gs->feat_dict->len); */
	}
	/* index_Array( target->kill_feat_quals, Killer_Feature_Qualifier *, feat_idx ) = kq; */
	append_val_Array( target->kill_feat_quals, kq );
      }
      else {
	/* must be in source tag - the killer is local to this particular source and target */
	
	Feature_Relation *ft_rel = index_Array(target->sources, 
						 Feature_Relation *, 
						 state->current_idx_2);
	if (ft_rel->kill_feat_quals == NULL) {
	  ft_rel->kill_feat_quals = new_Array(sizeof( Killer_Feature_Qualifier *), TRUE);
	  /* set_size_Array( ft_rel->kill_feat_quals, state->gs->feat_dict->len ); */
	}
	/*
	if ( index_Array( ft_rel->kill_feat_quals, Killer_Feature_Qualifier *, feat_idx ) == NULL )
	  index_Array( ft_rel->kill_feat_quals, Killer_Feature_Qualifier *, feat_idx ) = kq;
	else {
	  fprintf(stderr, "In tag 'killfeat', redefintion of killer qualifier\n"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
	*/
	append_val_Array( ft_rel->kill_feat_quals, kq );
      }

      /* need to mark the feat_idx feature as a killer feature */
      index_Array( state->gs->feat_info, Feature_Info *, kq->feat_idx)->is_killer_feat = TRUE;
    }
  }
  else {
    /* flag error and return */
    fprintf(stderr, "Tag 'killfeat' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}




/*********************************************************************
 FUNCTION: parse_Gaze_Structure_lengthfunc
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_lengthfunc( struct Parse_context *state, 
				      const char **attr ) {

  int i;
  char *file = NULL;
  
  if (state->tag_stack->len && 
      ! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[LENGTHFUNCTIONS].tag)) {
    /* top of tag stack should be 'lengthfunctions' */
    
    state->current_idx = -1;
    for( i=0; attr[i] && ! state->error; i += 2 ) {
      if (! strcmp( attr[i], "id" )) {
	if ((state->current_idx = 
	     dict_lookup( state->gs->len_fun_dict, attr[i+1])) < 0) {
	  fprintf(stderr, "In tag 'lengthfunc' attr 'id' has illegal value\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "file" ))
	file = strdup_util(attr[i+1]);
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'lengthfunc', attr '%s' not recognised\n", attr[i]); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    if (state->current_idx < 0) {
      state->error = XML_GetCurrentLineNumber( state->the_parser );
    }
    else if (file != NULL) {
	/* open file attr[i+1] */
	/* we need to append vals to the length function, but which one? */
	/* read in x-y pairs. Flag error if bad format */
      FILE *fp;

      if ( (fp= fopen( file, "r")) == NULL) {
	fprintf(stderr, "lengthfunction file '%s' could not be opened\n", file); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
      else {
	int  x;
	double y;
	Length_Function *lf = index_Array( state->gs->length_funcs, 
					     Length_Function *, 
					     state->current_idx );
	
	while( fscanf( fp, "%d %lf", &x, &y ) != EOF ) {
	  append_val_Array( lf->raw_x_vals, x );
	  append_val_Array( lf->raw_y_vals, y );
	}
	calc_Length_Function( lf );
      }

      free_util( file );
    }
  }
  else {
    fprintf(stderr, "Tag 'lengthfunc' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}



/*********************************************************************
 FUNCTION: parse_Gaze_Structure_lengthfunction
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_lengthfunction( struct Parse_context *state, 
						 const char **attr ) {

  char *id = NULL;
  double mul = 1.0;
  int i;
  
  if (state->tag_stack->len && 
      ! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[DECLARATIONS].tag)) {
    
    for( i=0; attr[i]; i += 2 ) {
      if (! strcmp( attr[i], "id" )) {
	/* check this is not a re-defintion */
	if ( dict_lookup(state->gs->len_fun_dict, attr[i+1]) >= 0) {
	fprintf(stderr, "In tag 'lengthfunction' attr 'id' has illeagl value\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
	else {
	  id = strdup_util( attr[i+1] );
	}
      }
      else if (! strcmp( attr[i], "mul" )) {
	mul = atof( attr[i+1] ); 
      }
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'lengthfunction', attr '%s' not recognised\n", attr[i]); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    if (! state->error) {
      if (id != NULL) { 
	Length_Function *lf = new_Length_Function( mul );
	append_val_Array( state->gs->len_fun_dict, id);
	append_val_Array( state->gs->length_funcs, lf);
      }
      else {
	/* no id attribute - bad news */
	fprintf(stderr, "In tag 'lengthfunction' attr 'id' is not defined\n");
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
  }
  else {
    fprintf(stderr, "Tag 'lengthfunction' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}



/*********************************************************************
 FUNCTION: parse_Gaze_Structure_lengthfunctions
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
*********************************************************************/
static void parse_Gaze_Structure_lengthfunctions( struct Parse_context *state, 
						  const char **attr ) {

  if (! state->tag_stack->len ||
      strcmp(index_Array(state->tag_stack, 
			   char *, 
			   state->tag_stack->len - 1), 
	     tag_map[GAZE].tag)) {
    fprintf(stderr, "Tag 'lengthfunctions' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}



/*********************************************************************
 FUNCTION: parse_Gaze_Structure_output
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_output( struct Parse_context *state, 
					 const char **attr ) {
  int i;

  Output_Qualifier *oq = new_Output_Qualifier();
  boolean in_source_tag;
  

  if (state->tag_stack->len && 
      ((! strcmp(index_Array(state->tag_stack, 
			       char *, 
			       state->tag_stack->len - 1), 
		 tag_map[SOURCE].tag)) ||
       (! strcmp(index_Array(state->tag_stack, 
			       char *, 
			       state->tag_stack->len - 1), 
		 tag_map[TARGET].tag)))) {
    
    in_source_tag =  ! strcmp(index_Array(state->tag_stack, 
					  char *, 
					  state->tag_stack->len - 1), 
			    tag_map[SOURCE].tag);

    for( i=0; attr[i]; i += 2 ) {
      if (! strcmp( attr[i], "feature" )) {
	oq->feature = strdup_util( attr[i+1] );
      }
      else if (! strcmp( attr[i], "strand" )) {
	oq->strand = strdup_util( attr[i+1] );
      }
      else if (! strcmp( attr[i], "frame" )) {
	oq->frame = strdup_util( attr[i+1] );
      }
      else if ( ! strcmp( attr[i], "all_regions" )) {
	if (! strcmp( attr[i+1], "TRUE" ))
	  oq->need_to_print = TRUE;
      }
      else {
	/* unrecognised attribute */
	state->error = XML_GetCurrentLineNumber( state->the_parser );
	fprintf(stderr, "In tag 'output', attr '%s' not recognised\n", attr[i]); 
      }
    }
    if (! state->error) {
      if (oq->feature != NULL || oq->strand != NULL || oq->frame != NULL) { 	
	Feature_Info *target = index_Array( state->gs->feat_info, Feature_Info *, state->current_idx);
	if ( in_source_tag )
	  index_Array(target->sources, Feature_Relation *, state->current_idx_2)->out_qual = oq;
	else 
	  target->out_qual = oq;
      }
      else {
	fprintf(stderr, "In tag 'output' you must define some output attributes!\n");
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
  }
  else {
    fprintf(stderr, "Tag 'output' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}






/*********************************************************************
 FUNCTION: parse_Gaze_Structure_point
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_point( struct Parse_context *state, 
					const char **attr ) {
  int i, x;
  double y;
  boolean parsed_x = FALSE;
  boolean parsed_y = FALSE;
  
  if (state->tag_stack->len && 
      ! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[LENGTHFUNC].tag)) {
    /* top of tag stack should be 'lengthfunctions' */

    Length_Function *len_fun = 
      index_Array( state->gs->length_funcs, Length_Function *, state->current_idx );

    /* check if this length function has already been determined (by reading a file) */
    
    if (len_fun->value_map == NULL) {

      for( i=0; attr[i] && ! state->error; i += 2 ) {
	if (! strcmp( attr[i], "x" )) {
	  if (! parsed_x) {
	  x = atoi( attr[i+1] );
	  parsed_x = TRUE;
	  }
	  else {
	    fprintf(stderr, "In tag 'point' attr 'x' must only be defined once\n");
	    state->error = XML_GetCurrentLineNumber( state->the_parser );
	  }
	}
	else if (! strcmp( attr[i], "y" )) {
	  if (! parsed_y) {
	    y = atof( attr[i+1] );
	    parsed_y = TRUE;
	  }
	  else {
	    fprintf(stderr, "In tag 'point' attr 'y' must only be defined once\n");
	    state->error = XML_GetCurrentLineNumber( state->the_parser );
	  }
	}
	else {
	  /* unrecognised attribute */
	  fprintf(stderr, "In tag 'point', attr '%s' not recognised\n", attr[i]); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      
      if (parsed_x && parsed_y) {
	append_val_Array( len_fun->raw_x_vals, x );
	append_val_Array( len_fun->raw_y_vals, y );
      }
      else {
	fprintf(stderr, "In tag 'point' attrs 'x' and 'y' must be defined\n");
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    else {
      fprintf(stderr, "Attempt to define both a file and 'point' for length function\n");
      state->error = XML_GetCurrentLineNumber( state->the_parser );
    }
  }
  else {
    fprintf(stderr, "Tag 'point' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}



/*********************************************************************
 FUNCTION: parse_Gaze_Structure_seg
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_seg( struct Parse_context *state, 
				      const char **attr ) {

  int i, index = 0;
  Segment *seg;

  for( i=0; attr[i]; i += 2 ) {
    if (! strcmp( attr[i], "id" )) {
      if ((index = dict_lookup( state->gs->seg_dict, attr[i+1])) < 0) {
	/* undeclared feature */
	fprintf(stderr, "In ag 'seg' attr 'id' has illegal valuet\n");
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    else {
      /* unrecognised attribute */
      fprintf(stderr, "In tag 'seg', attr '%s' not recognised\n", attr[i]); 
      state->error = XML_GetCurrentLineNumber( state->the_parser );
    }
  }  
  
  if (! state->error) {
    if (state->tag_stack->len) {
      if (! strcmp(index_Array(state->tag_stack, 
				 char *, 
				 state->tag_stack->len - 1), 
		   tag_map[GFFLINE].tag)) {
	
	GFF_to_features *gff2fts;
	gff2fts = index_Array( state->gs->gff_to_feats, 
				 GFF_to_features *,
				 state->gs->gff_to_feats->len - 1);
	seg = new_Segment();
	seg->seg_idx = index;
	append_val_Array( gff2fts->segments, seg );
      }
      else if (! strcmp(index_Array(state->tag_stack, 
				      char *, 
				      state->tag_stack->len - 1), 
			tag_map[DNAFEAT].tag)) {
	
	DNA_to_features *dna2fts;
	dna2fts = index_Array( state->gs->dna_to_feats, 
				 DNA_to_features *,
				 state->gs->dna_to_feats->len - 1);
	seg = new_Segment();
	seg->seg_idx = index;
	append_val_Array( dna2fts->segments, seg );	
      }
      else {
	fprintf(stderr, "Tag 'seg' not expected in this context\n");
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    else {
      /* Tag out of context */
      fprintf(stderr, "Tag 'seg' not expected in this context\n");
      state->error = XML_GetCurrentLineNumber( state->the_parser );
    }
  }
}




/*********************************************************************
 FUNCTION: parse_Gaze_Structure_segment
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_segment( struct Parse_context *state, 
					  const char **attr ) {
  char *id = NULL;
  double mul = 1.0;
  int i;
  Segment_Info *new_seg;
  boolean use_projected = TRUE;  /* default */
  boolean score_sum = TRUE;      /* default */
  boolean partial = TRUE;        /* default */

  if (state->tag_stack->len && 
      ! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[DECLARATIONS].tag)) {
    
    for( i=0; attr[i]; i += 2 ) {
      if (! strcmp( attr[i], "id" )) {
	/* check this is not a re-defintion */
	if ( dict_lookup(state->gs->seg_dict, attr[i+1]) >= 0) {
	  fprintf(stderr, "In tag 'segment',  this 'id' has been defined before\n"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
	else {
	  id = strdup_util( attr[i+1] );
	}
      }
      else if (! strcmp( attr[i], "mul" )) {
	mul = atof( attr[i+1] ); 
      }
      else if (! strcmp( attr[i], "partial")) {
	if (! strcmp( attr[i+1], "TRUE" )) {
	  /* this one is the default */
	  partial = TRUE;
	}
	else if (! strcmp( attr[i+1], "FALSE" )) {
	  partial = FALSE;
	}
	else {
	  fprintf(stderr, "In tag 'segment', attr 'partial' has an illegal value\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "scoring" )) {

	if (! strcmp( attr[i+1], "project_sum" )) {
	  /* this one is the default */
	  use_projected = TRUE;
	  score_sum = TRUE;
	}
	else if (! strcmp( attr[i+1], "project_max" )) {
	  use_projected = TRUE;
	  score_sum = FALSE;
	}
	else if (! strcmp( attr[i+1], "standard_sum" )) {
	  use_projected = FALSE;
	  score_sum = TRUE;
	}
	else if (! strcmp( attr[i+1], "standard_max" )) {
	  use_projected = FALSE;
	  score_sum = FALSE;
	}
	else {
	  fprintf(stderr, "In tag 'segment', attr 'scoring' has an illegal value\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'segment', attr '%s' not recognised\n", attr[i]); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    if (! state->error) {
      if (id != NULL) { 
      	new_seg = new_Segment_Info( mul );
	new_seg->use_projected = use_projected;
	new_seg->score_sum = score_sum;
	new_seg->partial = partial;

	append_val_Array( state->gs->seg_dict, id );
	append_val_Array( state->gs->seg_info, new_seg );
      }
      else {
	/* no id attribute - bad news */
	fprintf(stderr, "In tag 'segment' no 'id' attr defined\n"); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
  }
  else {
    /* flag error and return */
    fprintf(stderr, "Tag 'segment' not expected in this context\n"); 
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }

}



/*********************************************************************
 FUNCTION: parse_Gaze_Structure_source
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_source( struct Parse_context *state, 
					 const char **attr ) {

  int  i, source_idx, len_fun_idx;
  int *phase = NULL;
  int *mindist = NULL;
  int *maxdist = NULL;
  int *len_fun = NULL;

  if (state->tag_stack->len && 
      ! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[TARGET].tag)) {
    /* top of tag stack should be 'target' */

    source_idx = -1;
    for( i=0; attr[i] && ! state->error; i += 2 ) {
      if (! strcmp( attr[i], "id" )) {
	if ((source_idx = dict_lookup(state->gs->feat_dict, attr[i+1])) < 0) {
	  fprintf(stderr, "In tag 'source' attr 'id' has illegal value"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "phase" )) {
	phase = (int *) malloc_util( sizeof(int) );
        *phase = atoi( attr[i+1] ); 
	if (*phase < 0 || *phase > 2) {
	  fprintf(stderr, "In tag 'source' attr 'phase' has illegal value"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "mindis" )) {
	mindist = (int *) malloc_util( sizeof(int) );
        *mindist = atoi( attr[i+1] ); 
      }
      else if (! strcmp( attr[i], "maxdis" )) {
	maxdist = (int *) malloc_util( sizeof(int) );
        *maxdist = atoi( attr[i+1] ); 
      }
      else if (! strcmp( attr[i], "len_fun" )) {
	if ((len_fun_idx = dict_lookup( state->gs->len_fun_dict, attr[i+1])) >= 0) {
	  len_fun = (int *) malloc_util( sizeof(int) );
	  *len_fun = len_fun_idx;
	}
	else {
	  fprintf(stderr, "In tag 'source' attr 'len_fun' has illegal value"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'source', attr '%s' not recognised\n", attr[i]); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }

    if (! state->error) {
      Feature_Info *this_one;
      Feature_Relation *ft_src;

      this_one =
	index_Array( state->gs->feat_info, Feature_Info *, state->current_idx);
      if (source_idx >= 0) {
	ft_src = new_Feature_Relation();
	ft_src->target = state->current_idx;
	ft_src->source = source_idx;
	ft_src->min_dist = mindist;
	ft_src->max_dist = maxdist;
	ft_src->phase = phase;
	ft_src->len_fun = len_fun;

	index_Array( this_one->sources, Feature_Relation *, source_idx ) = ft_src;
	state->current_idx_2 = source_idx;
      }
      else {
	fprintf(stderr, "In tag 'source' no 'id' attr defined"); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
  }
  else {
    /* flag error and return */
    fprintf(stderr, "Tag 'source' not expected in this context"); 
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}




/*********************************************************************
 FUNCTION: parse_Gaze_Structure_takedna
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_takedna( struct Parse_context *state, 
					  const char **attr ) {
  
  int i, feat_idx = 0;
  int st_off = 0;
  int en_off = 0;
  boolean has_st_off = FALSE;
  boolean has_en_off= FALSE;

  if (state->tag_stack->len && 
      ! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[DNA2GAZE].tag)) {
    
    for( i=0; attr[i]; i += 2 ) {
      if (! strcmp( attr[i], "id" )) {
	if ((feat_idx = dict_lookup(state->gs->feat_dict, attr[i+1])) < 0) {
	  fprintf(stderr, "In tag 'takedna' attr 'id' has illegal value"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "st_off" )) {
	if (! has_st_off) {
	  st_off = atoi( attr[i+1] ); 
	  has_st_off = TRUE;
	}
	else {
	  fprintf(stderr, "In tag 'takedna', attr 'st_off' should only be defind once\n"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "en_off" )) {
	if (! has_en_off) {
	  en_off = atoi( attr[i+1] ); 
	  has_en_off = TRUE;
	}
	else {
	  fprintf(stderr, "In tag 'takedna', attr 'en_off' should only be defind once\n"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'source', attr '%s' not recognised\n", attr[i]); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }

    if (! state->error) {
      if (state->gs->take_dna == NULL) {
	state->gs->take_dna = new_Array( sizeof( StartEnd * ), TRUE);
	set_size_Array(state->gs->take_dna, state->gs->feat_dict->len);
      } 
      if (index_Array( state->gs->take_dna, StartEnd *, feat_idx ) == NULL)
	index_Array( state->gs->take_dna, StartEnd *, feat_idx ) = new_StartEnd( st_off, en_off );
      else {
	/* unrecognised attribute */
	fprintf(stderr, "Redef of feat In tag 'takedna'; canonly take one piece of dna for each feature\n");
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
  }
  else {
    fprintf(stderr, "Tag 'dnafeat' not expected in this context\n");
    state->error =  XML_GetCurrentLineNumber( state->the_parser );
  }
}




/*********************************************************************
 FUNCTION: parse_Gaze_Structure_target
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_target( struct Parse_context *state, 
					 const char **attr ) {
  int i;
  
  if (state->tag_stack->len && 
      ! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[MODEL].tag)) {
    /* top of tag stack should be 'model' */
    
    state->current_idx = -1;
    for( i=0; attr[i] && ! state->error; i += 2 ) {
      if (! strcmp( attr[i], "id" )) {
	if ((state->current_idx = 
	     dict_lookup( state->gs->feat_dict, attr[i+1])) < 0) {
	  fprintf(stderr, "In tag 'target' no 'id' attr has illegal value"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
	else {
	  /* set up the memory for the source information for
	     this feature. Target information will be filled in 
	     when the structure is symmetricalzed. */
	  Feature_Info *fi = index_Array( state->gs->feat_info, 
					    Feature_Info *, 
					    state->current_idx);
	  if (fi->sources == NULL) {
	    fi->sources = new_Array( sizeof( Feature_Relation *), TRUE);
	    set_size_Array(fi->sources, state->gs->feat_dict->len); 
	  }
	  else {
	    fprintf(stderr, "Error: Redefinition of target '%s'\n", 
		    index_Array(state->gs->feat_dict, char *, state->current_idx) ); 
	    state->error = XML_GetCurrentLineNumber( state->the_parser );
	  }
	}
      }
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'target', attr '%s' not recognised\n", attr[i]); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    if (state->current_idx < 0) {
      fprintf(stderr, "In tag 'target' no 'id' attr defined"); 
      state->error = XML_GetCurrentLineNumber( state->the_parser );
    }
  }
  else {
    fprintf(stderr, "Tag 'target' not expected in this context"); 
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}




/*********************************************************************
 FUNCTION: parse_Gaze_Structure_useseg
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_useseg( struct Parse_context *state, 
					 const char **attr ) {
  int  i;
  boolean in_source_tag = FALSE;
  boolean user_specified_scoring = FALSE;
  boolean user_specified_partial = FALSE;

  if (state->tag_stack->len && 
      ((! strcmp(index_Array(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[TARGET].tag)) ||
      (! strcmp(index_Array(state->tag_stack, 
			      char *, 
			      state->tag_stack->len - 1), 
		tag_map[SOURCE].tag)))) {
    /* top of tag stack should be 'target' or 'source' */
    Segment_Qualifier *seg_qual = new_Segment_Qualifier();

    in_source_tag =  strcmp(index_Array(state->tag_stack, 
					  char *, 
					  state->tag_stack->len - 1), 
			    tag_map[TARGET].tag);

    for( i=0; attr[i] && ! state->error; i += 2 ) {
      if (! strcmp( attr[i], "id" )) {
	/* check this is not a re-defintion */
	if ((seg_qual->seg_idx = dict_lookup(state->gs->seg_dict, attr[i+1])) < 0) {
	  fprintf(stderr, "In tag 'useseg', attr 'id' has an illegal value\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "target_phase" )) {
	if (seg_qual->has_src_phase) {
	  fprintf(stderr, 
		  "In tag 'useseg', attrs 'source_phase' and 'target_phase' cannot both be defined\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
	else {
	  seg_qual->has_tgt_phase = TRUE;
	  seg_qual->phase = atoi( attr[i+1] );
	  if (seg_qual->phase < 0 || seg_qual->phase > 2) {
	    fprintf(stderr, "In tag 'useseg', attr 'tgt_phase' has an illegal value\n");
	    state->error = XML_GetCurrentLineNumber( state->the_parser );
	  }
	}
      }
      else if (! strcmp( attr[i], "source_phase" )) {
	if (seg_qual->has_tgt_phase) {
	  fprintf(stderr, 
		  "In tag 'useseg', attrs 'source_phase' and 'target_phase' cannot both be defined\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
	else {
	  seg_qual->has_src_phase = TRUE;
	  seg_qual->phase = atoi( attr[i+1] );
	  if (seg_qual->phase < 0 || seg_qual->phase > 2) {
	    fprintf(stderr, "In tag 'useseg', attr 'src_phase' has an illegal value\n");
	    state->error = XML_GetCurrentLineNumber( state->the_parser );
	  }
	}
      }
      else if (! strcmp( attr[i], "exact" )) {
	if (! strcmp( attr[i+1], "source" )) {
	  seg_qual->is_exact_src = TRUE;
	}
	else if (! strcmp( attr[i+1], "target" )) {
	  seg_qual->is_exact_tgt = TRUE;
	}
	else if (! strcmp( attr[i+1], "both" )) {
	  seg_qual->is_exact_src = TRUE;
	  seg_qual->is_exact_tgt = TRUE;
	}
	else {
	  fprintf(stderr, "In tag 'useseg', attr 'exact' has an illegal value\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "scoring" )) {
	user_specified_scoring = TRUE;

	if (! strcmp( attr[i+1], "project_sum" )) {
	  /* this one is the default */
	  seg_qual->use_projected = TRUE;
	  seg_qual->score_sum = TRUE;
	}
	else if (! strcmp( attr[i+1], "project_max" )) {
	  seg_qual->use_projected = TRUE;
	  seg_qual->score_sum = FALSE;
	}
	else if (! strcmp( attr[i+1], "standard_sum" )) {
	  seg_qual->use_projected = FALSE;
	  seg_qual->score_sum = TRUE;
	}
	else if (! strcmp( attr[i+1], "standard_max" )) {
	  seg_qual->use_projected = FALSE;
	  seg_qual->score_sum = FALSE;
	}
	else {
	  fprintf(stderr, "In tag 'useseg', attr 'scoring' has an illegal value\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "partial" )) {
	user_specified_partial = TRUE;

	if (! strcmp( attr[i+1], "FALSE" ))
	  seg_qual->partial = FALSE;
	else if (! strcmp( attr[i+1], "TRUE" ))
	  seg_qual->partial = TRUE;
	else {
	  fprintf(stderr, "In tag 'useseg', attr 'partial' has an illegal value\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'useseg', attr '%s' not recognised\n", attr[i]); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    if (! state->error) {
      Feature_Info *target = index_Array( state->gs->feat_info, Feature_Info *, state->current_idx);

      /* firstly, if user did not specify a scoring scheme, then we go back
	 to the segment info to look for the scoring scheme */

      if (! user_specified_partial)
	seg_qual->partial = index_Array( state->gs->seg_info, Segment_Info *, seg_qual->seg_idx )->partial;

      if (! user_specified_scoring) {
	seg_qual->use_projected = index_Array( state->gs->seg_info, Segment_Info *, seg_qual->seg_idx )->use_projected;
	seg_qual->score_sum = index_Array( state->gs->seg_info, Segment_Info *, seg_qual->seg_idx )->score_sum;
      }

      /* if we came across an "exact" qualifier, then we cannot use projected segments,
	 so, die and tell the user */

      if  (seg_qual->is_exact_src || seg_qual->is_exact_tgt) {
	if (user_specified_scoring && (seg_qual->use_projected || seg_qual->score_sum)) {
	  /* flag error and return */
	  fprintf(stderr, "In tag 'useseg': 'exact' qualifier cannot be used with this seg scoring\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
	else {
	  seg_qual->use_projected = FALSE;
	  seg_qual->score_sum = FALSE;
	}
      }
      
      if (! state->error) {
	/* if we are in the 'target' tag, we need to update the general segment 
	   qualifiers array (the one that is copied to each source as it is encountered).
	   If not, we update the Feature_Relation for the source-target directly */
	
	if (! in_source_tag) { 
	  if (target->seg_quals == NULL) {
	    target->seg_quals = new_Array( sizeof( Segment_Qualifier *), TRUE);
	    /* set_size_Array( target->seg_quals, state->gs->seg_dict->len ); */
	  }
	  /* index_Array( target->seg_quals, Segment_Qualifier *, seg_idx ) = seg_qual; */
	  append_val_Array( target->seg_quals, seg_qual );
	}
	else {
	  Feature_Relation *ft_rel = index_Array(target->sources, Feature_Relation *, state->current_idx_2);
	  if (ft_rel->seg_quals == NULL) {
	    ft_rel->seg_quals = new_Array(sizeof( Segment_Qualifier *), TRUE);
	  }
	  append_val_Array( ft_rel->seg_quals, seg_qual );
	}
      }
    }
  }
  else {
    /* flag error and return */
    fprintf(stderr, "Not expecting 'useseg' in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}



/********************************************************************/
/********************* main parsing routine *************************/
/********************************************************************/


/*********************************************************************
 FUNCTION: parse_Gaze_Structure
 DESCRIPTION:
   This function is called whenever an XML closing tag is parsed
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/

Gaze_Structure *parse_Gaze_Structure( char *structure_file_nm ) {
  int i;
  boolean end_of_file = FALSE;
  struct Parse_context *state;
  Gaze_Structure *gs;
  FILE *structure_file;

  XML_Parser p = XML_ParserCreate(NULL);
  if (! p) {
    fprintf(stderr, "Couldn't allocate memory for parser\n");
    exit(-1);
  }

  gs = new_Gaze_Structure();
  state = (struct Parse_context *) malloc_util ( sizeof( struct Parse_context ) );
  state->finished_parsing = FALSE;
  state->error = 0;
  state->the_parser = p;
  state->current_idx = -1;
  state->tag_stack = new_Array( sizeof( char *), TRUE);
  state->gs = gs;

  XML_SetElementHandler(p, &start_tag_structure, &end_tag_structure);
  XML_SetUserData( p,  state );

  structure_file = fopen( structure_file_nm, "r" );

  while (! end_of_file && ! state->finished_parsing && ! state->error ) {
    int len;
    void *buffer = XML_GetBuffer(p, PARSE_BUFFER_SIZE );
    if (buffer == NULL) {
      fprintf(stderr, "Could not allocate memory for buffer\n");
      exit(1);
    }

    len = fread(buffer, 1, PARSE_BUFFER_SIZE, structure_file );
    if (ferror(stdin)) {
      fprintf(stderr, "Read error\n");
      exit(-1);
    }
    end_of_file = feof(structure_file);

    if (! XML_Parse(p, buffer, len, end_of_file)) {
      fprintf(stderr, "Parse error at line %d:\n%s\n",
              XML_GetCurrentLineNumber(p),
              XML_ErrorString(XML_GetErrorCode(p)));
      exit(-1);
    }
  }
  fclose( structure_file );

  if (state->error) {
    fprintf(stderr, "Error occurred at line %d\n", state->error);
    free_Gaze_Structure( gs );
    gs = NULL;
  }
  else {
    /* calculate derived information for length functions, and structure
       in general */
    
    for(i=0; i < gs->length_funcs->len; i++) {
      Length_Function *lf = index_Array( gs->length_funcs, Length_Function *, i );
      /* Some length functions will have already been read in from files, so no need
	 to calculate them */
      if (lf->value_map == NULL) 
	calc_Length_Function( lf );
    }
    fill_in_Gaze_Structure( gs );
  }    

  for(i=0; i < state->tag_stack->len; i++) {
    /* stack should be empty, but just in case ... */
    free_util( index_Array( state->tag_stack, char *, i));
  }
  free_Array( state->tag_stack, TRUE);
  free_util( state );

  /* free the parse state object */
    
  XML_ParserFree( p );

  return gs;
}


