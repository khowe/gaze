/*  Last edited: Jan 23 16:55 2002 (klh) */
/**********************************************************************
 ** File: output.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description :
 **********************************************************************/


#include "output.h"




/*********************************************************************
 FUNCTION: print_gff_Gaze_Features
 DESCRIPTION:
   Prints out the given path as a list of Gaze features (rather than
   as a list of derived feature names and the regions between them,
   as in print_GFF_path)
 RETURNS:
 ARGS: 
 NOTES:
   1. This function does not check whether the given feature list
   represents a true path.
   2. Score, Strand and frame are given uninformative values. However,
   the output generated will mainly be used for reading PATHs back
   into GAZE, which does not require score, strand, or frame 
   information
 *********************************************************************/
void print_GFF_Gaze_Features( FILE *fh,
			      GArray *fts, 
			      Gaze_Structure *gs,
			      char *seq_name) {
  Feature *f1;
  int i;

  fprintf( fh, "##gff-version 2\n");
  fprintf( fh, "## List of GAZE features implying a gene structure of %s\n", seq_name);
  fprintf( fh, "##  Score of path : %.6f\n", g_array_index(fts,Feature *,fts->len-1)->path_score);
  fprintf( fh, "##  Forward score : %.6f\n", g_array_index(fts,Feature *,fts->len-1)->forward_score);
  fprintf( fh, "##  Probability of path : %.6f\n", exp(g_array_index(fts,Feature *,fts->len-1)->path_score -
						       g_array_index(fts,Feature *,fts->len-1)->forward_score));

  /* leave off BEGIN and END */
  
  for( i=1; i < fts->len - 2; i++) {
    f1 = g_array_index( fts, Feature *, i);

    fprintf(fh, "%s\tGAZE\t%s\t%d\t%d\t%.3f\t.\t.\n", 
	    seq_name, g_array_index( gs->feat_dict,
				     char *,
				     f1->feat_idx ),
	    f1->real_pos.s, f1->real_pos.e,
	    f1->score);
  }
}




/*********************************************************************
 FUNCTION: print_gff_path
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void print_GFF_path( FILE *fh,
		     GArray *fts, 
		     Gaze_Structure *gs,
		     char *seq_name) {
  Feature *f1, *f2 = NULL;
  Feature_Info *f1_info, *f2_info;
  Feature_Relation *src;
  char empty = '.';
  char strand, frame;
  int i;

  fprintf( fh, "##gff-version 2\n");
  fprintf( fh, "## GAZE gene structure of %s\n", seq_name);
  fprintf( fh, "##  Score of path : %.6f\n", g_array_index(fts,Feature *,fts->len-1)->path_score);
  /*
  fprintf( fh, "##  Forward score : %.6f\n", g_array_index(fts,Feature *,fts->len-1)->forward_score);
  fprintf( fh, "##  Probability of path : %.6f\n", exp(g_array_index(fts,Feature *,fts->len-1)->path_score -
						       g_array_index(fts,Feature *,fts->len-1)->forward_score));
  */

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

    fprintf(fh, "%s\tGAZE\t%s\t%d\t%d\t%.3f\t.\t.\n", 
	    seq_name, g_array_index( gs->feat_dict,
				     char *,
				     f1->feat_idx ),
	    f1->real_pos.s, f1->real_pos.e,
	    f1->score);

    strand = src->out_strand ? src->out_strand : empty;
    frame = src->out_frame ? src->out_frame : empty;
    
    fprintf(fh, "%s\tGAZE\t%s\t%d\t%d\t%.3f\t%c\t%c\n", 
	    seq_name, src->out_feature,
	    /*
	    f1->real_pos.s + f1_info->start_offset,
	    f2->real_pos.e - f2_info->end_offset,*/
	    f1->adj_pos.s, f2->adj_pos.e,
	    f2->path_score - f1->path_score - f2->score,
	    strand, frame);
  }

  /* finally, the last feature */

  if (f2 != NULL) {
    fprintf(fh, "%s\tGAZE\t%s\t%d\t%d\t%.3f\t.\t.\n", 
	    seq_name, g_array_index( gs->feat_dict,
				     char *,
				     f2->feat_idx ),
	    f2->real_pos.s, f2->real_pos.e, f2->score);
  }
}




/*********************************************************************
 FUNCTION: print_post_probs
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void print_post_probs( FILE *fh,
		       GArray *fts, 
		       double thresh,
		       Gaze_Structure *gs, 
		       char *seq_name) {

  int i;

  fprintf( fh, "##gff-version 2\n");
  fprintf( fh, "## GAZE posterior probabilities of features from %s\n", seq_name);
  fprintf( fh, "##     Fend = %.10f,     Bbegin = %.10f\n", 
	   g_array_index( fts, Feature *, fts->len - 1)->forward_score,
	   g_array_index( fts, Feature *, 0)->backward_score );


  for(i=0; i < fts->len; i++) {
    Feature *f = g_array_index( fts, Feature *, i);
    
    double fw = f->forward_score;
    double bk = f->backward_score;
    double sum = g_array_index(fts, Feature *, 0)->backward_score;
    double prob = exp( fw + bk - sum );
    
    fprintf(fh, "%s\tGAZE\t%s\t%d\t%d\t%.8f\t.\t.\t%s\n", 
	    seq_name, g_array_index( gs->feat_dict,
				     char *,
				     f->feat_idx ),
	    f->real_pos.s, f->real_pos.e, prob, f->is_correct?"TRUE":"");
  }
}
