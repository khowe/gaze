/*  Last edited: Apr 23 15:45 2002 (klh) */
/**********************************************************************
 ** File: output.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description :
 **********************************************************************/


#include "output.h"


/*********************************************************************
 FUNCTION: Gaze_Output
 DESCRIPTION:
   Thisd 
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Gaze_Output *new_Gaze_Output( void ) {
  Gaze_Output *out = (Gaze_Output *) g_malloc( sizeof( Gaze_Output ) );

  out->fh = NULL;
  out->seq_name = NULL;    
  out->posterior = FALSE;
  out->use_threshold = FALSE;
  out->threshold = 0.0;

  return out;
}


/*********************************************************************
 FUNCTION: free_Gaze_Output
 DESCRIPTION:
   An Output_line is a shopping bag of things needed to produce
   a line of GAZE output
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Gaze_Output( Gaze_Output *out ) {
  if (out != NULL) {
    if (out->seq_name != NULL)
      g_free( out->seq_name );

    g_free( out );
  }
}



/*********************************************************************
 FUNCTION: print_gff_path
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void print_GFF_path( Gaze_Output *out,
		     GArray *fts,
		     Gaze_Structure *gs) {

  Feature *f1, *f2 = NULL;
  Feature_Info *f1_info, *f2_info;
  Feature_Relation *src;
  double feature_score, region_score;
  int i;


  fprintf( out->fh, "##  Score of path : %.6f\n", g_array_index(fts,Feature *,fts->len-1)->path_score);
  fprintf( out->fh, "##  Probability of path : %.6f\n", 
	   exp ( g_array_index(fts, Feature *,fts->len-1)->path_score -
		 g_array_index(fts, Feature *, fts->len-1)->forward_score) );

  for( i=0; i < fts->len - 1; i++) {
    f1 = g_array_index( fts, Feature *, i);
    f2 = g_array_index( fts, Feature *, i+1);

    f1_info = g_array_index( gs->feat_info, Feature_Info *, f1->feat_idx);
    f2_info = g_array_index( gs->feat_info, Feature_Info *, f2->feat_idx);

    src = g_array_index( f2_info->sources, Feature_Relation *, f1->feat_idx );

    if (src == NULL) {
      fprintf( stderr, "Error: Illegal feature list\n");
      return;
    }

    /* first print feature itself, then print the region */
    feature_score = f1->score;
    if (out->posterior)
      feature_score = exp( f1->forward_score + 
			   f1->backward_score - 
			   g_array_index( fts, Feature *, 0)->backward_score);
    

    fprintf(out->fh, "%s\tGAZE\t%s\t%d\t%d\t%.3f\t.\t.\n", 
	    out->seq_name, 
	    g_array_index( gs->feat_dict, char *, f1->feat_idx ),
	    f1->real_pos.s, 
	    f1->real_pos.e,
	    feature_score);
    
    region_score = f2->path_score - f1->path_score - f2->score; 
    
    if (out->posterior) 
      region_score = exp( f1->forward_score +
			  f2->backward_score +
			  region_score + 
			  f2->score -
			  g_array_index( fts, Feature *, 0)->backward_score );

    fprintf(out->fh, "%s\tGAZE\t%s\t%d\t%d\t%.3f\t%s\t%s\n", 
	    out->seq_name, 
	    src->out_qual->feature != NULL ? src->out_qual->feature : "Not_Given",
	    f1->adj_pos.s, 
	    f2->adj_pos.e,
	    region_score,
	    src->out_qual->strand != NULL ? src->out_qual->strand : ".",
	    src->out_qual->frame != NULL ? src->out_qual->frame : "." );
  }

  /* finally, the last feature */

  feature_score = f2->score;
  if (out->posterior)
    feature_score = exp( f2->forward_score + 
			 f2->backward_score - 
			 g_array_index( fts, Feature *, 0)->backward_score);
  if (f2 != NULL) {
    fprintf(out->fh, "%s\tGAZE\t%s\t%d\t%d\t%.3f\t.\t.\n", 
	    out->seq_name, 
	    g_array_index( gs->feat_dict, char *, f2->feat_idx ),
	    f2->real_pos.s, 
	    f2->real_pos.e, 
	    feature_score);
  }
}



/*********************************************************************
 FUNCTION: print_GFF_Gaze_Features
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void print_GFF_Gaze_Features( Gaze_Output *out,
			      GArray *fts,
			      Gaze_Structure *gs) {


  int i, discarded = 0;

  if (out->posterior) {
    fprintf( out->fh, "## GAZE feature scored by posterior probability\n");
    fprintf( out->fh, "##     Fend = %.10f,     Bbegin = %.10f\n", 
	     g_array_index( fts, Feature *, fts->len - 1)->forward_score,
	     g_array_index( fts, Feature *, 0)->backward_score );
  }

  if (out->use_threshold)
    fprintf( out->fh, "##   (only features scoring above %.2f are shown)\n", out->threshold );

  for(i=0; i < fts->len; i++) {
    Feature *f = g_array_index( fts, Feature *, i);
    double score = f->score;

    if (out->posterior) {
      score = exp( f->forward_score + 
		   f->backward_score - 
		   g_array_index( fts, Feature *, 0)->backward_score );
    }

    if (! out->use_threshold || score > out->threshold) 
      fprintf(out->fh, "%s\tGAZE\t%s\t%d\t%d\t%.8f\t.\t.\t%s\n", 
	      out->seq_name, 
	      g_array_index( gs->feat_dict,
			     char *,
			     f->feat_idx ),
	      f->real_pos.s, f->real_pos.e, score, f->is_correct?"TRUE":"");
    else
      discarded++;
  }

  if (out->use_threshold)
    fprintf( out->fh, "## Discarded %d features with post. probs. %.2f or below\n", discarded, out->threshold );
}



/*********************************************************************
 FUNCTION: write_GFF_header
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void write_GFF_header( Gaze_Output *out, int seq_begin, int seq_end ) {

  fprintf( out->fh, "##gff-version 2\n");
  fprintf( out->fh, "##sequence-region %s %d %d\n", out->seq_name, seq_begin, seq_end );
  fprintf( out->fh, "##source-version GAZE 1.0\n");
}

