/*  Last edited: Nov 22 18:26 2001 (klh) */
/**********************************************************************
 ** File: str_parse.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/


#include "str_parse.h"


struct Parse_context {
  GArray *tag_stack;
  gboolean finished_parsing;
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
static void parse_Gaze_Structure_gfffeat( struct Parse_context *, const char **);
static void parse_Gaze_Structure_killdna( struct Parse_context *, const char **);
static void parse_Gaze_Structure_killfeat( struct Parse_context *, const char **);
static void parse_Gaze_Structure_lengthfunc( struct Parse_context *, const char **);
static void parse_Gaze_Structure_lengthfunction( struct Parse_context *,const char **);
static void parse_Gaze_Structure_lengthfunctions( struct Parse_context *, const char **);
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
      GFFFEAT,
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
  { "gfffeat", &parse_Gaze_Structure_gfffeat },
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
  gboolean match = FALSE;

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
	t =  g_strdup( el );
	g_array_append_val(state->tag_stack, t);
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
	char *temp = g_array_index( state->tag_stack, char *, index );
	g_array_remove_index( state->tag_stack, index );
	g_free( temp );
	
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
      strcmp(g_array_index(state->tag_stack, 
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
      strcmp(g_array_index(state->tag_stack, 
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
  gboolean has_score = FALSE;
  double score = 0.0;
  int i;
  DNA_to_features *new_dna2ft;
  
  if (state->tag_stack->len && 
      ! strcmp(g_array_index(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[DNA2GAZE].tag)) {
    
    for( i=0; attr[i]; i += 2 ) {
      if (! strcmp( attr[i], "pattern" )) {
	pattern = g_strdup( attr[i+1] );
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
	g_array_append_val( state->gs->dna_to_feats, new_dna2ft );
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
  int i, index;
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
      if (! strcmp(g_array_index(state->tag_stack, 
				 char *, 
				 state->tag_stack->len - 1), 
		   tag_map[GFFFEAT].tag)) {
	
	GFF_to_features *gff2fts;
	gff2fts = g_array_index( state->gs->gff_to_feats, 
				 GFF_to_features *,
				 state->gs->gff_to_feats->len - 1);
	feat = new_Feature();
	feat->feat_idx = index;
	g_array_append_val( gff2fts->features, feat );
      }
      else if (! strcmp(g_array_index(state->tag_stack, 
				      char *, 
				      state->tag_stack->len - 1), 
			tag_map[DNAFEAT].tag)) {
	
	DNA_to_features *dna2fts;
	dna2fts = g_array_index( state->gs->dna_to_feats, 
				 DNA_to_features *,
				 state->gs->dna_to_feats->len - 1);
	feat = new_Feature();
	feat->feat_idx = index;
	feat->score = score;
	g_array_append_val( dna2fts->features, feat );	
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
      ! strcmp(g_array_index(state->tag_stack, 
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
	  id = g_strdup( attr[i+1] );
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

	g_array_insert_val( state->gs->feat_dict, state->gs->feat_dict->len - 1, id );
	g_array_insert_val( state->gs->feat_info, state->gs->feat_info->len - 1, new_feat );
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
      strcmp(g_array_index(state->tag_stack, 
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
      strcmp(g_array_index(state->tag_stack, 
			   char *, 
			   state->tag_stack->len - 1), 
	     tag_map[GAZE].tag)) {
    fprintf(stderr, "Tag 'gff2gaze' not expected in this context\n");
    state->error = XML_GetCurrentLineNumber( state->the_parser );
  }
}



/*********************************************************************
 FUNCTION: parse_Gaze_Structure_gfffeat
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static void parse_Gaze_Structure_gfffeat( struct Parse_context *state, 
					  const char **attr ) {

  char *type = NULL;
  char *source = NULL;
  char *strand = NULL;
  int i;
  GFF_to_features *new_gff2ft;
  
  if (state->tag_stack->len && 
      ! strcmp(g_array_index(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[GFF2GAZE].tag)) {

    for( i=0; attr[i]; i += 2 ) {
      if (! strcmp( attr[i], "feature" )) {
	type = g_strdup( attr[i+1] );
      }
      else if (! strcmp( attr[i], "source" )) {
	source = g_strdup( attr[i+1] );
      }
      else if (! strcmp( attr[i], "strand" )) {
	strand = g_strdup( attr[i+1] );
      }
      else {
	/* unrecognised attribute */
	state->error = XML_GetCurrentLineNumber( state->the_parser );
	fprintf(stderr, "In tag 'gfffeat', attr '%s' not recognised\n", attr[i]); 
      }
    }
    if (! state->error) {
      if (type != NULL || source != NULL || strand != NULL) { 
      	new_gff2ft = new_GFF_to_features();
	new_gff2ft->gff_feature = type;
	new_gff2ft->gff_source = source;
	new_gff2ft->gff_strand = strand;
	g_array_append_val( state->gs->gff_to_feats, new_gff2ft );

      }
      else {
	/* no id attribute - bad news */
	fprintf(stderr, "In tag 'gfffeat' attr 'id' is not defined\n");
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
  }
  else {
    fprintf(stderr, "Tag 'gfffeat' not expected in this context\n");
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

  int i,j;
  int src_dna_idx = -1;
  int tgt_dna_idx = -1;
  char *new_motif = NULL;

  if (state->tag_stack->len && 
      ! strcmp(g_array_index(state->tag_stack, 
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
	  new_motif = g_strdup( attr[i+1] );
	  for(j=0; new_motif[j] != '\0'; j++)
	    new_motif[j] = tolower(new_motif[j]);
	  if ((src_dna_idx = dict_lookup(state->gs->motif_dict, attr[i+1])) < 0) {
	    g_array_append_val( state->gs->motif_dict, new_motif );
	    src_dna_idx = state->gs->motif_dict->len - 1;
	  }
	  else 
	    g_free( new_motif );
	}
      }
      else if (! strcmp( attr[i], "target_dna" )) {
	if (tgt_dna_idx >= 0) {
	  /* error */
	  fprintf(stderr, "In tag 'killdna', attribute 'target_dna' cannot be defined more than once\n"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
	else {
	  new_motif = g_strdup( attr[i+1] );
	  for(j=0; new_motif[j] != '\0'; j++)
	    new_motif[j] = tolower(new_motif[j]);
	  if ((tgt_dna_idx = dict_lookup(state->gs->motif_dict, attr[i+1])) < 0) {
	    g_array_append_val( state->gs->motif_dict, new_motif );
	    tgt_dna_idx = state->gs->motif_dict->len - 1;
	  }
	  else 
	    g_free( new_motif );
	}
      }
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'killdna', attr '%s' not recognised\n", attr[i]); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }


    if (! state->error) {
      Feature_Info *tgt_info = g_array_index( state->gs->feat_info, Feature_Info *, state->current_idx);
      Feature_Relation *ft_rel = g_array_index( tgt_info->sources, Feature_Relation *, state->current_idx_2);
      Killer_DNA_Qualifier *kdq = new_Killer_DNA_Qualifier();

      kdq->src_dna = src_dna_idx;
      kdq->tgt_dna = tgt_dna_idx;

      if (ft_rel->kill_dna_quals == NULL) 
	ft_rel->kill_dna_quals = g_array_new(FALSE, TRUE, sizeof( Killer_DNA_Qualifier *));
      
      g_array_append_val( ft_rel->kill_dna_quals, kdq );
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

  int  i, feat_idx;
  gboolean in_source_tag = FALSE;
  
  if (state->tag_stack->len && 
      ((! strcmp(g_array_index(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[TARGET].tag)) ||
      (! strcmp(g_array_index(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[SOURCE].tag)))) {
    /* top of tag stack should be 'target' or 'source' */

    Killer_Feature_Qualifier *kq = new_Killer_Feature_Qualifier();

    in_source_tag =  strcmp(g_array_index(state->tag_stack, 
					  char *, 
					  state->tag_stack->len - 1), 
			    tag_map[TARGET].tag);

    for( i=0; attr[i] && ! state->error; i += 2 ) {
      if (! strcmp( attr[i], "id" )) {
	if ((feat_idx = dict_lookup(state->gs->feat_dict, attr[i+1])) < 0) {
	  fprintf(stderr, "In tag 'killfeat' attr 'id' has an illegal value\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "phase" )) {
	kq->has_phase = TRUE;
	kq->phase = atoi( attr[i+1] );
	if (kq->phase < 0 || kq->phase > 2) {
	  fprintf(stderr, "In tag 'killfeat' attr 'tgt_phase' has an illegal value\n");
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'killfeat', attr '%s' not recognised\n", attr[i]); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }

    if (! state->error) {
      Feature_Info *target = g_array_index( state->gs->feat_info, Feature_Info *, state->current_idx);

      /* if we are in the 'target' tag, we need to update the general kill_off_up
	 array (the one that applies to all sources for the given target)
	 If not, we update the Feature_Relation for the source-target directly */

      if (!in_source_tag) {
	/* must be in target tag - the killer is global to all sources for this target */
	if (target->kill_feat_quals_up == NULL) {
	  target->kill_feat_quals_up = g_array_new( FALSE, TRUE, sizeof( Killer_Feature_Qualifier * ));
	  g_array_set_size(target->kill_feat_quals_up, state->gs->feat_dict->len); 
	}
	g_array_index( target->kill_feat_quals_up, Killer_Feature_Qualifier *, feat_idx ) = kq;
      }
      else {
	/* must be in source tag - the killer is local to this particular source and target */
	
	Feature_Relation *ft_rel = g_array_index(target->sources, 
						 Feature_Relation *, 
						 state->current_idx_2);
	if (ft_rel->kill_feat_quals == NULL) {
	  ft_rel->kill_feat_quals = g_array_new(FALSE, TRUE, sizeof( Killer_Feature_Qualifier *));
	  g_array_set_size( ft_rel->kill_feat_quals, state->gs->feat_dict->len );
	}
	if ( g_array_index( ft_rel->kill_feat_quals, Killer_Feature_Qualifier *, feat_idx ) == NULL )
	  g_array_index( ft_rel->kill_feat_quals, Killer_Feature_Qualifier *, feat_idx ) = kq;
	else {
	  fprintf(stderr, "In tag 'killfeat', redefintion of killer qualifier\n"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }

      /* need to mark the feat_idx feature as a killer feature */
      g_array_index( state->gs->feat_info, Feature_Info *, feat_idx)->is_killer_feat = TRUE;
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
      ! strcmp(g_array_index(state->tag_stack, 
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
	file = g_strdup(attr[i+1]);
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
	Length_Function *lf = g_array_index( state->gs->length_funcs, 
					     Length_Function *, 
					     state->current_idx );
	
	while( fscanf( fp, "%d %lf", &x, &y ) != EOF ) {
	  g_array_append_val( lf->raw_x_vals, x );
	  g_array_append_val( lf->raw_y_vals, y );
	}
	calc_Length_Function( lf );
      }

      g_free( file );
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
      ! strcmp(g_array_index(state->tag_stack, 
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
	  id = g_strdup( attr[i+1] );
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
	g_array_append_val( state->gs->len_fun_dict, id);
	g_array_append_val( state->gs->length_funcs, lf);
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
      strcmp(g_array_index(state->tag_stack, 
			   char *, 
			   state->tag_stack->len - 1), 
	     tag_map[GAZE].tag)) {
    fprintf(stderr, "Tag 'lengthfunctions' not expected in this context\n");
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
  gboolean parsed_x = FALSE;
  gboolean parsed_y = FALSE;
  
  if (state->tag_stack->len && 
      ! strcmp(g_array_index(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[LENGTHFUNC].tag)) {
    /* top of tag stack should be 'lengthfunctions' */

    Length_Function *len_fun = 
      g_array_index( state->gs->length_funcs, Length_Function *, state->current_idx );

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
	g_array_append_val( len_fun->raw_x_vals, x );
	g_array_append_val( len_fun->raw_y_vals, y );
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

  int i, index;
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
      if (! strcmp(g_array_index(state->tag_stack, 
				 char *, 
				 state->tag_stack->len - 1), 
		   tag_map[GFFFEAT].tag)) {
	
	GFF_to_features *gff2fts;
	gff2fts = g_array_index( state->gs->gff_to_feats, 
				 GFF_to_features *,
				 state->gs->gff_to_feats->len - 1);
	seg = new_Segment();
	seg->seg_idx = index;
	g_array_append_val( gff2fts->segments, seg );
      }
      else if (! strcmp(g_array_index(state->tag_stack, 
				      char *, 
				      state->tag_stack->len - 1), 
			tag_map[DNAFEAT].tag)) {
	
	DNA_to_features *dna2fts;
	dna2fts = g_array_index( state->gs->dna_to_feats, 
				 DNA_to_features *,
				 state->gs->dna_to_feats->len - 1);
	seg = new_Segment();
	seg->seg_idx = index;
	g_array_append_val( dna2fts->segments, seg );	
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
  gboolean use_projected = TRUE;  /* default */
  gboolean score_sum = TRUE;      /* default */

  if (state->tag_stack->len && 
      ! strcmp(g_array_index(state->tag_stack, 
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
	  id = g_strdup( attr[i+1] );
	}
      }
      else if (! strcmp( attr[i], "mul" )) {
	mul = atof( attr[i+1] ); 
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

	g_array_append_val( state->gs->seg_dict, id );
	g_array_append_val( state->gs->seg_info, new_seg );
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
  char *out_feature = NULL;
  char *out_strand = NULL;
  char *out_frame = NULL;
  int *len_fun = NULL;

  if (state->tag_stack->len && 
      ! strcmp(g_array_index(state->tag_stack, 
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
	phase = (int *) g_malloc( sizeof(int) );
        *phase = atoi( attr[i+1] ); 
	if (*phase < 0 || *phase > 2) {
	  fprintf(stderr, "In tag 'source' attr 'phase' has illegal value"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "mindis" )) {
	mindist = (int *) g_malloc( sizeof(int) );
        *mindist = atoi( attr[i+1] ); 
      }
      else if (! strcmp( attr[i], "maxdis" )) {
	maxdist = (int *) g_malloc( sizeof(int) );
        *maxdist = atoi( attr[i+1] ); 
      }
      else if (! strcmp( attr[i], "len_fun" )) {
	if ((len_fun_idx = dict_lookup( state->gs->len_fun_dict, attr[i+1])) >= 0) {
	  len_fun = (int *) g_malloc( sizeof(int) );
	  *len_fun = len_fun_idx;
	}
	else {
	  fprintf(stderr, "In tag 'source' attr 'len_fun' has illegal value"); 
	  state->error = XML_GetCurrentLineNumber( state->the_parser );
	}
      }
      else if (! strcmp( attr[i], "out_feat" )) {
	out_feature = g_strdup( attr[i+1] );
      }
      else if (! strcmp( attr[i], "out_str" )) {
	out_strand = g_strdup( attr[i+1] );
      }
      else if (! strcmp( attr[i], "out_frm" )) {
	out_frame = g_strdup( attr[i+1] );
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
	g_array_index( state->gs->feat_info, Feature_Info *, state->current_idx);
      if (source_idx >= 0) {
	ft_src = new_Feature_Relation();
	ft_src->target = state->current_idx;
	ft_src->source = source_idx;
	ft_src->min_dist = mindist;
	ft_src->max_dist = maxdist;
	ft_src->phase = phase;
	ft_src->out_feature = out_feature;
	ft_src->out_strand = out_strand;
	ft_src->out_frame = out_frame;
	ft_src->len_fun = len_fun;
	/*
	if (this_one->seg_quals != NULL) {
	  ft_src->seg_quals = g_array_new( FALSE, TRUE, sizeof(Segment_Qualifier *));
	  g_array_set_size( ft_src->seg_quals, this_one->seg_quals->len );
	  for (i=0; i < this_one->seg_quals->len; i++) {
	    g_array_index( ft_src->seg_quals, Segment_Qualifier *, i) =
	      clone_Segment_Qualifier( g_array_index( this_one->seg_quals, Segment_Qualifier *, i));
	  }
	}
	*/
	g_array_index( this_one->sources, Feature_Relation *, source_idx ) = ft_src;
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
  
  int i, feat_idx;
  int st_off = 0;
  int en_off = 0;
  gboolean has_st_off = FALSE;
  gboolean has_en_off= FALSE;

  if (state->tag_stack->len && 
      ! strcmp(g_array_index(state->tag_stack, 
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
	state->gs->take_dna = g_array_new( FALSE, TRUE, sizeof( StartEnd * ));
	g_array_set_size(state->gs->take_dna, state->gs->feat_dict->len);
      } 
      if (g_array_index( state->gs->take_dna, StartEnd *, feat_idx ) == NULL) 
	g_array_index( state->gs->take_dna, StartEnd *, feat_idx ) = new_StartEnd( st_off, en_off );
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
      ! strcmp(g_array_index(state->tag_stack, 
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
	  Feature_Info *fi = g_array_index( state->gs->feat_info, 
					    Feature_Info *, 
					    state->current_idx);
	  fi->sources = g_array_new( FALSE, TRUE, sizeof( Feature_Relation *));
	  g_array_set_size(fi->sources, state->gs->feat_dict->len); 
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
  int  i, seg_idx;
  gboolean in_source_tag = FALSE;
  gboolean user_specified_scoring = FALSE;

  if (state->tag_stack->len && 
      ((! strcmp(g_array_index(state->tag_stack, 
			     char *, 
			     state->tag_stack->len - 1), 
	       tag_map[TARGET].tag)) ||
      (! strcmp(g_array_index(state->tag_stack, 
			      char *, 
			      state->tag_stack->len - 1), 
		tag_map[SOURCE].tag)))) {
    /* top of tag stack should be 'target' or 'source' */
    Segment_Qualifier *seg_qual = new_Segment_Qualifier();

    in_source_tag =  strcmp(g_array_index(state->tag_stack, 
					  char *, 
					  state->tag_stack->len - 1), 
			    tag_map[TARGET].tag);

    for( i=0; attr[i] && ! state->error; i += 2 ) {
      if (! strcmp( attr[i], "id" )) {
	/* check this is not a re-defintion */
	if ((seg_idx = dict_lookup(state->gs->seg_dict, attr[i+1])) < 0) {
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
      else {
	/* unrecognised attribute */
	fprintf(stderr, "In tag 'useseg', attr '%s' not recognised\n", attr[i]); 
	state->error = XML_GetCurrentLineNumber( state->the_parser );
      }
    }
    if (! state->error) {
      Feature_Info *target = g_array_index( state->gs->feat_info, Feature_Info *, state->current_idx);

      /* firstly, if user did not specify a scoring scheme, then we go back
	 to the segment info to look for the scoring scheme */

      if (! user_specified_scoring) {
	Segment_Info *si = g_array_index( state->gs->seg_info, Segment_Info *, seg_idx );
	seg_qual->use_projected = si->use_projected;
	seg_qual->score_sum = si->score_sum;
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
	    target->seg_quals = g_array_new(FALSE, TRUE, sizeof( Segment_Qualifier *));
	    g_array_set_size( target->seg_quals, state->gs->seg_dict->len );
	  }
	  g_array_index( target->seg_quals, Segment_Qualifier *, seg_idx ) = seg_qual;
	}
	else {
	  Feature_Relation *ft_rel = g_array_index(target->sources, Feature_Relation *, state->current_idx_2);
	  if (ft_rel->seg_quals == NULL) {
	    ft_rel->seg_quals = g_array_new(FALSE, TRUE, sizeof( Segment_Qualifier *));
	    g_array_set_size( ft_rel->seg_quals, state->gs->seg_dict->len );
	  }
	  if ( g_array_index( ft_rel->seg_quals, Segment_Qualifier *, seg_idx ) == NULL )
	    g_array_index( ft_rel->seg_quals, Segment_Qualifier *, seg_idx ) = seg_qual;
	  else {
	    fprintf(stderr, "In tag 'useseg', redefintion of segment qualifier\n"); 
	    state->error = XML_GetCurrentLineNumber( state->the_parser );
	  }
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

Gaze_Structure *parse_Gaze_Structure( FILE *structure_file ) {
  int i;
  gboolean end_of_file = FALSE;
  struct Parse_context *state;
  Gaze_Structure *gs;

  XML_Parser p = XML_ParserCreate(NULL);
  if (! p) {
    fprintf(stderr, "Couldn't allocate memory for parser\n");
    exit(-1);
  }

  gs = new_Gaze_Structure();
  state = (struct Parse_context *) g_malloc ( sizeof( struct Parse_context ) );
  state->finished_parsing = FALSE;
  state->error = 0;
  state->the_parser = p;
  state->current_idx = -1;
  state->tag_stack = g_array_new( FALSE, TRUE, sizeof( char *));
  state->gs = gs;

  XML_SetElementHandler(p, &start_tag_structure, &end_tag_structure);
  XML_SetUserData( p,  state );

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

  if (state->error) {
    fprintf(stderr, "Error occurred at line %d\n", state->error);
    free_Gaze_Structure( gs );
    gs = NULL;
  }
  else {
    /* calculate derived information for length functions, and structure
       in general */
    
    for(i=0; i < gs->length_funcs->len; i++) {
      Length_Function *lf = g_array_index( gs->length_funcs, Length_Function *, i );
      /* Some length functions will have already been read in from files, so no need
	 to calculate them */
      if (lf->value_map == NULL) 
	calc_Length_Function( lf );
    }
    fill_in_Gaze_Structure( gs );
  }    

  for(i=0; i < state->tag_stack->len; i++) {
    /* stack should be empty, but just in case ... */
    g_free( g_array_index( state->tag_stack, char *, i));
  }
  g_array_free( state->tag_stack, TRUE);
  g_free( state );

  /* free the parse state object */
    
  XML_ParserFree( p );

  return gs;
}


