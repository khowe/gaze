/*  Last edited: Jul 22 11:03 2002 (klh) */
/**********************************************************************
 ** File: structure.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description :
     This file takes care of the representation and parsing of the 
     XML encoded gaze structure file. 

     The Gaze structure file encodes three classes of information:

     1. Info about Features, and how they relate to each other
     2. Info about making features from DNA sequence
     3. Info about making features from a GFF file

 **********************************************************************/


#include "structure.h"

/*********************************************************************
 FUNCTION: free_Gaze_Structure
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Gaze_Structure( Gaze_Structure *gs ) {
  int i;

  if (gs != NULL) {
    if (gs->feat_dict != NULL) {
      for(i=0; i < gs->feat_dict->len; i++) 
	free_util( index_Array( gs->feat_dict, char *, i));
      free_Array( gs->feat_dict, TRUE );
    }
    if (gs->seg_dict != NULL) {
      for(i=0; i < gs->seg_dict->len; i++) 
	free_util( index_Array( gs->seg_dict, char *, i));
      free_Array( gs->seg_dict, TRUE );
    }
    if (gs->len_fun_dict != NULL) {
      for(i=0; i < gs->len_fun_dict->len; i++) 
	free_util( index_Array( gs->len_fun_dict, char *, i));
      free_Array( gs->len_fun_dict, TRUE );
    }
    if (gs->motif_dict != NULL) {
      for(i=0; i < gs->motif_dict->len; i++) 
	free_util( index_Array( gs->motif_dict, char *, i));
      free_Array( gs->motif_dict, TRUE );
    }
    if (gs->feat_info != NULL) {
      for(i=0; i < gs->feat_info->len; i++) 
	free_Feature_Info( index_Array( gs->feat_info, Feature_Info *, i));
      free_Array( gs->feat_info, TRUE );
    }
    if (gs->seg_info != NULL) {
      for(i=0; i < gs->seg_info->len; i++) 
	free_Segment_Info( index_Array(gs->seg_info, Segment_Info *, i));
      free_Array( gs->seg_info, TRUE );
    }
    if (gs->length_funcs != NULL) {
      for(i=0; i < gs->length_funcs->len; i++) 
	free_Length_Function( index_Array(gs->length_funcs, Length_Function *, i));
      free_Array( gs->length_funcs, TRUE );
    }
    if (gs->take_dna != NULL) {
      for(i=0; i < gs->take_dna->len; i++) 
	free_StartEnd( index_Array( gs->take_dna, StartEnd *, i ));
      free_Array( gs->take_dna, TRUE );
    }
    if (gs->dna_to_feats != NULL) {
      for(i=0; i < gs->dna_to_feats->len; i++) 
	free_DNA_to_features( index_Array(gs->dna_to_feats, DNA_to_features *, i));
      free_Array( gs->dna_to_feats, TRUE );
    }
    if (gs->gff_to_feats != NULL) {
      for(i=0; i < gs->gff_to_feats->len; i++) 
	free_GFF_to_features( index_Array(gs->gff_to_feats, GFF_to_features *, i));
      free_Array( gs->gff_to_feats, TRUE );
    }

    free_util( gs );
  }
}


/*********************************************************************
 FUNCTION: new_Gaze_Structure
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Gaze_Structure *new_Gaze_Structure( void ) {
  Gaze_Structure *g_str;
  char *begin = strdup_util("BEGIN");
  char *end = strdup_util("END");
  Feature_Info *begin_info = new_Feature_Info( 0, 0, 0.0 );
  Feature_Info *end_info = new_Feature_Info( 0, 0, 0.0 );

  g_str = (Gaze_Structure *) malloc_util( sizeof( Gaze_Structure ) );
  g_str->feat_dict = new_Array( sizeof( char *), TRUE );
  g_str->seg_dict = new_Array( sizeof( char *), TRUE );
  g_str->len_fun_dict = new_Array( sizeof( char *), TRUE );
  g_str->motif_dict = new_Array( sizeof( char *), TRUE );
  g_str->feat_info = new_Array( sizeof( Feature_Info *), TRUE);
  g_str->seg_info = new_Array( sizeof( Segment_Info *), TRUE);
  g_str->length_funcs = new_Array( sizeof( Length_Function *), TRUE);
  g_str->dna_to_feats = new_Array( sizeof( DNA_to_features *), TRUE);
  g_str->gff_to_feats = new_Array( sizeof( GFF_to_features *), TRUE);
  g_str->take_dna = NULL;

  /* need to add BEGIN and END features to the feature dictionary, 
     and create dummy Feature_Info objects for them */

  append_val_Array( g_str->feat_dict, begin );
  append_val_Array( g_str->feat_dict, end );
  append_val_Array( g_str->feat_info, begin_info );
  append_val_Array( g_str->feat_info, end_info );

  return g_str;
}


/*********************************************************************
 FUNCTION: write_Gaze_Structure
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void write_Gaze_Structure( Gaze_Structure *gs, FILE *out ) {
  int i,j,k;

  fprintf(out,"Gaze Internal Structure\n");
  fprintf( out,  "-----------------------\n\n");
  
  for( i=0; i < gs->feat_info->len; i++) {    
    Feature_Info *fi = index_Array( gs->feat_info, Feature_Info *, i);
    fprintf(out, "Feature definition: %s\n", index_Array(gs->feat_dict, char*, i));
    fprintf(out, "  Start offset:%d\n", fi->start_offset);
    fprintf(out, "  End   offset:%d\n", fi->end_offset);
    if (gs->take_dna != NULL) {
      if (index_Array(gs->take_dna, StartEnd *, i ) != NULL)
	fprintf(out, "  Take DNA: st_off:%d en_off:%d\n", 
		index_Array(gs->take_dna, StartEnd *, i )->s,
		index_Array(gs->take_dna, StartEnd *, i )->e);
    }
    if (fi->kill_feat_quals != NULL) {
      fprintf(out, "  Killers for all upstream sources:\n");
      for( j=0; j < fi->kill_feat_quals->len; j++) {
	Killer_Feature_Qualifier *kq = index_Array(fi->kill_feat_quals, Killer_Feature_Qualifier *, j);
	if ( kq != NULL) {
	  fprintf(out, "    id:%s", index_Array( gs->feat_dict, char *, j));
	  if (kq->has_src_phase)
	    fprintf(out, ", src_phase:%d", kq->phase );
	  else if (kq->has_tgt_phase)
	    fprintf(out, ", tgt_phase:%d", kq->phase );
	  fprintf( out, "\n" );
	}
      }
    }

    if (fi->sources != NULL) {
      fprintf(out, "  Sources:\n");
      for( j=0; j < fi->sources->len; j++) {
	if (index_Array(fi->sources, Feature_Relation *, j) != NULL) {
	  Feature_Relation *fs = index_Array(fi->sources, Feature_Relation *, j);
	  int *phase = fs->phase;
	  int *mindis = fs->min_dist;
	  int *maxdis = fs->max_dist;
	  int *len_fun = fs->len_fun;
	  fprintf(out, "    src:%s", index_Array(gs->feat_dict, char *, j));
	  if (mindis != NULL) 
	    fprintf(out, " min:%d", *mindis);
	  if (maxdis != NULL) 
	    fprintf(out, " max:%d", *maxdis);
	  if (phase != NULL) 
	    fprintf(out, " phase:%d", *phase);
	  if (len_fun != NULL) {
	    char *name = index_Array( gs->len_fun_dict, char *, *len_fun );
	    fprintf(out, " len_fun:%s", name);
	  }
	  if (fs->out_qual != NULL && fs->out_qual->feature != NULL)
	    fprintf(out, " label:%s", fs->out_qual->feature);
	  fprintf(out, "\n");
	  
	  if (fs->seg_quals != NULL ) {
	    fprintf( out, "      Use segments:\n");
	    for (k=0; k < fs->seg_quals->len; k++) {
	      Segment_Qualifier *sq = index_Array( fs->seg_quals, Segment_Qualifier *, k);
	      if (sq != NULL) {
		fprintf(out, "       id:%s ", index_Array( gs->seg_dict, char *, k));
		if (sq->has_tgt_phase)
		  fprintf(out, "tgt_phase:%d ", sq->phase);
		if (sq->has_src_phase)
		  fprintf(out, "src_phase:%d ", sq->phase);
		if (sq->is_exact_src)
		  fprintf(out, "src_exact ");
		if (sq->is_exact_tgt)
		  fprintf(out, "tgt_exact");
		if (sq->use_projected)
		  fprintf(out, "projected ");
		if (sq->score_sum)
		  fprintf(out, "score_sum");
		fprintf(out, "\n");
	      }
	    }
	  }
	  if (fs->kill_feat_quals != NULL ) {
	    fprintf( out, "      Killer features specific this source-target:\n");
	    for (k=0; k < fs->kill_feat_quals->len; k++) {
	      Killer_Feature_Qualifier *kq = index_Array( fs->kill_feat_quals, 
							    Killer_Feature_Qualifier *, 
							    k);
	      if (kq != NULL) {
		fprintf(out, "       id:%s ", index_Array( gs->feat_dict, char *, k));
		if (kq->has_src_phase)
		  fprintf(out, "src_phase:%d ", kq->phase);
		else if (kq->has_tgt_phase)
		  fprintf(out, "tgt_phase:%d ", kq->phase);
		fprintf(out, "\n");
	      }
	    }
	  }
	  if (fs->kill_dna_quals != NULL ) {
	    fprintf( out, "      Killer DNA motifs for this source-target:\n");
	    for (k=0; k < fs->kill_dna_quals->len; k++) {
	      Killer_DNA_Qualifier *kdq = index_Array( fs->kill_dna_quals, 
							 Killer_DNA_Qualifier *, 
							 k);
	      fprintf(out, "       src dna:%s tgt dna:%s\n", 
		      index_Array( gs->motif_dict, char *, kdq->src_dna ),
		      index_Array( gs->motif_dict, char *, kdq->tgt_dna ));

	    }
	  }
	}
      }
    }

    if (fi->is_killer_feat) {
      fprintf(out, "  This is a Feature Killer\n");
    }
    fprintf( out,  "\n");

  }
  fprintf( out,  "\n");  

  for( i=0; i < gs->seg_info->len; i++) {    
    Segment_Info *si = index_Array( gs->seg_info, Segment_Info *, i);
    fprintf(out, "Segment definition: %s\n", index_Array(gs->seg_dict, char*, i));
    fprintf(out, "  Multiplier:%f\n\n", si->multiplier);
  }
  fprintf( out,  "\n");  

  fprintf( out,  "Conversion information\n");
  fprintf( out,  "----------------------\n\n");

  fprintf(out, "From GFF:\n\n");
  for(i=0; i < gs->gff_to_feats->len; i++ ) {
    GFF_to_features *gff2fts = 
      index_Array( gs->gff_to_feats, GFF_to_features *, i);

    fprintf(out, "GFF: type:%s, source:%s, strand:%s\n", 
	    (gff2fts->gff_feature)?gff2fts->gff_feature:"n/a",
	    (gff2fts->gff_source)?gff2fts->gff_source:"n/a",
	    (gff2fts->gff_strand)?gff2fts->gff_strand:"n/a");

    for(j=0; j < gff2fts->features->len; j++) {
      Feature *ft = index_Array( gff2fts->features, Feature *, j);
      int idx = ft->feat_idx;
      
      fprintf(out, "  Feature: %s", index_Array( gs->feat_dict, char *, idx));
    }
    for(j=0; j < gff2fts->segments->len; j++) {
      Segment *ft = index_Array( gff2fts->segments, Segment *, j);
      int idx = ft->seg_idx;
      
      fprintf(out, "  Segment: %s\n", index_Array( gs->seg_dict, char *, idx));
    }
    fprintf(out,"\n");
  }
  fprintf(out, "\n");

  fprintf(out, "From DNA:\n\n");
  for(i=0; i < gs->dna_to_feats->len; i++ ) {
    DNA_to_features *dna2fts = 
      index_Array( gs->dna_to_feats, DNA_to_features *, i);

    fprintf(out, "DNA: pattern:%s\n", dna2fts->dna_motif);
    for(j=0; j < dna2fts->features->len; j++) {
      Feature *ft = index_Array( dna2fts->features, Feature *, j);
      int idx = ft->feat_idx;
      
      fprintf(out, "  Feature: %s", index_Array( gs->feat_dict, char *, idx));
    }
    for(j=0; j < dna2fts->segments->len; j++) {
      Segment *ft = index_Array( dna2fts->segments, Segment *, j);
      int idx = ft->seg_idx;
      
      fprintf(out, "  Segment: %s\n", index_Array( gs->seg_dict, char *, idx));
    }
    fprintf(out,"\n");
  }
  fprintf(out, "\n");

}


/*********************************************************************
 FUNCTION: fill_in_Gaze_Structure
 DESCRIPTION:
   This function, called after parsing, fills in all the structure
   information that is not know until parsing is finished (e.g. 
   propagating target usesegs to their sources) and also attempts to make 
   the structure in a sense symmetrical, by deriving upstream killers
   for sources from downstream ones for targets (not that segments are
   not made symmetrical in this way - their directionality is controlled
   by the user via src_phase and tgt_phase
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void fill_in_Gaze_Structure( Gaze_Structure *gs) {
  int tgt_idx, src_idx, i;

  for (tgt_idx=0; tgt_idx < gs->feat_info->len; tgt_idx++) {
    Feature_Info *tgt_inf = index_Array( gs->feat_info, Feature_Info *, tgt_idx);
    
    if (tgt_inf->sources != NULL) {
      for(src_idx=0; src_idx < tgt_inf->sources->len; src_idx++) {
	Feature_Relation *src_tgt = index_Array( tgt_inf->sources, Feature_Relation *, src_idx);

	if (src_tgt != NULL) {
	  /* Firstly, propagate global information locally to local sources */

	  if (tgt_inf->seg_quals != NULL) {
	    if (src_tgt->seg_quals == NULL)
	      src_tgt->seg_quals = new_Array( sizeof(Segment_Qualifier *), TRUE);
	    for (i=0; i < tgt_inf->seg_quals->len; i++) {
	      Segment_Qualifier *new_sq = 
		clone_Segment_Qualifier( index_Array( tgt_inf->seg_quals, Segment_Qualifier *, i));
	      append_val_Array( src_tgt->seg_quals, new_sq );
	    }
	  }

	  if (tgt_inf->kill_feat_quals != NULL) {
	    if (src_tgt->kill_feat_quals == NULL) 
	      src_tgt->kill_feat_quals = new_Array( sizeof(Killer_Feature_Qualifier *), TRUE);
	    for (i=0; i < tgt_inf->kill_feat_quals->len; i++) {
	      Killer_Feature_Qualifier *new_kq = 
		clone_Killer_Feature_Qualifier( index_Array( tgt_inf->kill_feat_quals, Killer_Feature_Qualifier *, i));
	      append_val_Array( src_tgt->kill_feat_quals, new_kq );
	    }
	  }

	  if (src_tgt->out_qual == NULL) {
	    if (tgt_inf->out_qual != NULL)
	      src_tgt->out_qual = clone_Output_Qualifier( tgt_inf->out_qual );
	    else
	      src_tgt->out_qual = new_Output_Qualifier();
	  }
	  else {
	    /* copy across global information where there is no conflict */
	    if (tgt_inf->out_qual != NULL) {
	      if (src_tgt->out_qual->feature == NULL && tgt_inf->out_qual->feature == NULL)
		src_tgt->out_qual->feature = strdup_util (tgt_inf->out_qual->feature );
	      if (src_tgt->out_qual->strand == NULL && tgt_inf->out_qual->strand == NULL)
		src_tgt->out_qual->strand = strdup_util (tgt_inf->out_qual->strand );
	      if (src_tgt->out_qual->frame == NULL && tgt_inf->out_qual->frame == NULL)
		src_tgt->out_qual->frame = strdup_util (tgt_inf->out_qual->frame );
	      if (tgt_inf->out_qual->need_to_print == TRUE) 
		src_tgt->out_qual->need_to_print = TRUE;
	    }
	  }	  
	}
      }
    }

    /* Finally, destroy the global qualifiers, beacause they are never used */
    if (tgt_inf->seg_quals != NULL) {
      for (i=0; i < tgt_inf->seg_quals->len; i++) {
	/* This ensures that we don't overwrite more specific useseg specified in a source */
	if (index_Array( tgt_inf->seg_quals, Segment_Qualifier *, i) != NULL) {
	  free_Segment_Qualifier( index_Array( tgt_inf->seg_quals, Segment_Qualifier *, i));
	  index_Array(  tgt_inf->seg_quals, Segment_Qualifier *, i) = NULL;
	}
      }
      free_Array( tgt_inf->seg_quals, TRUE );
      tgt_inf->seg_quals = NULL;
    }
    if (tgt_inf->kill_feat_quals != NULL) {
      for (i=0; i < tgt_inf->kill_feat_quals->len; i++) {
	/* This ensures that we don't overwrite more specific useseg specified in a source */
	if (index_Array( tgt_inf->kill_feat_quals, Killer_Feature_Qualifier *, i) != NULL) {
	  free_Killer_Feature_Qualifier( index_Array( tgt_inf->kill_feat_quals, 
							Killer_Feature_Qualifier *, 
							i));
	  index_Array(  tgt_inf->kill_feat_quals, Killer_Feature_Qualifier *, i) = NULL;
	}
      }
      free_Array( tgt_inf->kill_feat_quals, TRUE );
      tgt_inf->kill_feat_quals = NULL;
    }
    if (tgt_inf->out_qual != NULL) {
      free_Output_Qualifier( tgt_inf->out_qual );
      tgt_inf->out_qual = NULL;
    }
  }
}
