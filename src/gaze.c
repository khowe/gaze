/*  Last edited: May  7 09:44 2002 (klh) */
/**********************************************************************
 ** File: gaze.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#include <glib.h>
#include <stdio.h>
#include <stdlib.h>

#include "options.h"
#include "info.h"
#include "str_parse.h"
#include "structure.h"
#include "features.h"
#include "g_engine.h"
#include "output.h"



static char gaze_usage_string[] = "\
Usage: gaze <options>\n\
Options are:\n\
\n\
Input files:\n\
\
 -structure_file <s>    XML file containing the gaze structure\n\
 -gff_file <s>          name of the GFF file containing the features\n\
 -dna_file <s>          file containing the DNA sequence\n\
 -gene_file <s>         name of a GFF file containing a user-specified gene structure\n\
 -defaults <s>          name of the file of default options (def: './gaze.defaults')\n\
\n\
Output files:\n\
 -out_file <s>       print gene structure to given file (def: stdout)\n\
\n\
Where to look for genes:\n\
 -begin_dna <n>         residue number to start looking for genes (def: 1)\n\
 -end_dna <n>           residue number to stop looking for genes (def: sequence length)\n\
 -offset_dna <n>        residue number of the first residue in the DNA file (def: 1)\n\
\n\
Output format:\n\
 -output <c>            The output type. Should be one of the following:\n\
                          B  (B)est (highest scoring) gene structure (default)\n\
                          S  (S)ampled gene structure based on forward score\n\
                          F  (F)eature candidates (most sensibly used with -posterior)\n\
                          R  (R)egion candidates (for regions indicated in structure file\n\
\n\
Other options:\n\
\
 -posterior             Show element scores as posteriror probs. rather than raw scores\n\
 -selected              look out for Selected features in input\n\
 -full_calc             perform full dynamic programming (as opposed to faster heurstic method)\n\
 -verbose               write basic progess information to stderr\n\
 -help                  show this message\n" ;

static Option options[] = {
  { "-begin_dna", INT_ARG },
  { "-end_dna", INT_ARG },
  { "-offset_dna", INT_ARG },
  { "-dna_file", STRING_ARG },
  { "-structure_file", STRING_ARG },
  { "-gff_file", STRING_ARG },
  { "-out_file", STRING_ARG },
  { "-gene_file", STRING_ARG },
  { "-defaults_file", STRING_ARG },
  { "-output", CHAR_ARG },
  { "-help", NO_ARGS },
  { "-verbose", NO_ARGS },
  { "-posterior", NO_ARGS },
  { "-full_calc", NO_ARGS },
  { "-sigma", FLOAT_ARG }
};


static struct {
  int begin_dna;
  int end_dna;
  int offset_dna;
  double sigma;
  char *structure_file_name;
  FILE *structure_file;
  GArray *gff_file_names; /* of string */
  GArray *gff_files;      /* of FILE */
  char *dna_file_name;
  FILE *dna_file;
  char *out_file_name;
  FILE *out_file;
  char *gene_file_name;
  FILE *gene_file;
  enum {
    BEST_PATH,
    SAMPLE_PATH,
    ALL_FEATURES,
    ALL_REGIONS 
  } output;    
  gboolean full_calc;
  gboolean use_selected;
  gboolean verbose;

  gboolean posterior;
  gboolean use_threshold;
  double threshold;

} gaze_options;



/*********************************************************************
 FUNCTION: process_Gaze_Options
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static gboolean process_Gaze_Options(char *optname,
				     char *optarg ) {
  int i;
  gboolean options_error = FALSE;

  if (strcmp(optname, "-begin_dna") == 0) gaze_options.begin_dna = atoi( optarg );
  else if (strcmp(optname, "-end_dna") == 0) gaze_options.end_dna = atoi( optarg );	     
  else if (strcmp(optname, "-offset_dna") == 0) gaze_options.offset_dna = atoi( optarg );
  else if (strcmp(optname, "-sigma") == 0) gaze_options.sigma = atof( optarg );
  else if (strcmp(optname, "-selected") == 0) gaze_options.use_selected = TRUE;	     
  else if (strcmp(optname, "-verbose") == 0) gaze_options.verbose = TRUE;
  else if (strcmp(optname, "-posterior") == 0) gaze_options.posterior = TRUE;  
  else if (strcmp(optname, "-full_calc") == 0) gaze_options.full_calc = TRUE;
  else if (strcmp(optname, "-output") == 0) {
    char output = (char) toupper( (int) optarg[0] );
    if (output == 'B')
      gaze_options.output = BEST_PATH;
    else if (output == 'S')
      gaze_options.output = SAMPLE_PATH;
    else if (output == 'F')
      gaze_options.output = ALL_FEATURES;
    else if (output == 'R')
      gaze_options.output = ALL_REGIONS;
    else {
      fprintf( stderr, "Unrecognised qualifier for -output (%s)\n", optarg );
      options_error = TRUE;
    }
  }
  else if (strcmp(optname, "-threshold") == 0) {
    gaze_options.use_threshold = TRUE;
    gaze_options.threshold = atof( optarg );
  }
  else if (strcmp(optname, "-out_file") == 0) {
    if ((gaze_options.out_file = fopen( optarg, "w")) == NULL) {
      fprintf( stderr, "Could not open output file %s for writing\n", optarg );
      options_error = TRUE;
    }
    else {
      if (gaze_options.out_file_name != NULL)
	g_free( gaze_options.out_file_name );
      gaze_options.out_file_name = g_strdup( optarg );
    }
  }
  else if (strcmp(optname, "-dna_file") == 0) {
    if ((gaze_options.dna_file = fopen( optarg, "r")) == NULL) {
      fprintf( stderr, "Could not open dna file %s for reading\n", optarg );
      options_error = TRUE;
    }
    else {
      if (gaze_options.dna_file_name != NULL)
	g_free( gaze_options.dna_file_name );
      gaze_options.dna_file_name = g_strdup( optarg );
    }
  }
  else if (strcmp(optname, "-structure_file") == 0) {
    if ((gaze_options.structure_file = fopen( optarg, "r")) == NULL) {
      fprintf( stderr, "Could not open structure file %s for reading\n", optarg );
      options_error = TRUE;
    }
    else {
      if (gaze_options.structure_file_name != NULL)
	g_free( gaze_options.structure_file_name );
      gaze_options.structure_file_name = g_strdup( optarg );
    }
  }
  else if (strcmp(optname, "-gene_file") == 0) {
    if ((gaze_options.gene_file = fopen( optarg, "r")) == NULL) {
      fprintf( stderr, "Could not open gene file %s for reading\n", optarg );
      options_error = TRUE;
    }
    else {
      if (gaze_options.gene_file_name != NULL)
	g_free( gaze_options.gene_file_name );
      gaze_options.gene_file_name = g_strdup( optarg );
    }
  }
  else if (strcmp(optname, "-gff_file") == 0) {
    /* I allow multiple gff files to be specified. So, fo each one, 
       we have to check that it hasn't already been opened, and
       if it hasn't, open it and store it in the file name list */
    gboolean match = FALSE;
    for (i=0; i < gaze_options.gff_file_names->len; i++) {
      if (! strcmp( optarg, g_array_index( gaze_options.gff_file_names, char *, i)))
	match = TRUE;
    }
    
    if (match)
      /* need to warn about duplicate feature file */
	fprintf( stderr, "Warning: feature file %s was given more than once\n", optarg );
    else {
      /* Try to open the file, and if success, store both the file name
	 and the file handle */
      FILE *tmp_f = fopen( optarg, "r");
      if (tmp_f == NULL) {
	fprintf( stderr, "Could not open feature file %s for reading\n", optarg );
	options_error = TRUE;
      }
      else {
	char *tmp_f_name = g_strdup( optarg );
	g_array_append_val( gaze_options.gff_file_names, tmp_f_name );
	g_array_append_val( gaze_options.gff_files, tmp_f );
      }
    } 
  }

  /* else do nothing */

  return options_error;
}



/*********************************************************************
 FUNCTION: parse_command_line
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
static int parse_command_line( int argc, char *argv[] ) {
  gboolean options_error = FALSE;
  gboolean help_wanted = FALSE;
  int optindex;
  char *optname, *optarg;
  FILE *defaults_fh = NULL;

  while (get_option( argc, argv, options,
		     sizeof(options) / sizeof( Option ),
		     &optindex, &optname, &optarg, &options_error )){
    if (strcmp(optname, "-help") == 0) {
      help_wanted = TRUE;
    }
  }

  if (optindex != argc)
    options_error = TRUE;

  if (help_wanted) { 
    fprintf( stderr, gaze_usage_string );
    return 0;
  }
  else if (options_error) {
    fprintf( stderr, "(gaze -h gives usage).\n" );
    return 0;
  }

  /* at this point, we can be sure that all options on the
     command line are valid, so we need not check 
     options_error any more */

  while (get_option( argc, argv, options,
		     sizeof(options) / sizeof( Option ),
		     &optindex, &optname, &optarg, &options_error )){
    if (strcmp(optname, "-defaults_file") == 0) {
      /* open the defaults file and process it */
      if ((defaults_fh = fopen( optarg, "r" )) == NULL) {
	fprintf(stderr,  "Could not open defaults file %s for reading\n", optarg );
	return 0;
      }
    }
  }

  if (defaults_fh == NULL) {
    defaults_fh = fopen( "defaults.gaze", "r" );
  } 

  gaze_options.begin_dna = 1;
  gaze_options.end_dna = 300000000;   /* need to do something about this */
  gaze_options.offset_dna = 1;
  gaze_options.sigma = 1.0;
  gaze_options.dna_file_name = NULL;
  gaze_options.dna_file = NULL;
  gaze_options.gff_file_names = g_array_new( FALSE, TRUE, sizeof( char *) );
  gaze_options.gff_files = g_array_new( FALSE, TRUE, sizeof( FILE *) );
  gaze_options.structure_file_name = NULL;
  gaze_options.structure_file = NULL;
  gaze_options.out_file_name = g_strdup( "stdout ");
  gaze_options.out_file = stdout;
  gaze_options.gene_file_name = NULL;
  gaze_options.gene_file = NULL;;
  gaze_options.use_selected = FALSE;
  gaze_options.verbose = FALSE;
  gaze_options.full_calc = FALSE;
  gaze_options.posterior = FALSE;
  gaze_options.use_threshold = FALSE;
  gaze_options.threshold = 0.0;
  gaze_options.output = BEST_PATH;

  if (process_default_Options( defaults_fh, &process_Gaze_Options  )) {
    return 0;
  }

  /* anything left on the command line has priority, so will
     overwrite settings made so far */

  
  while (! options_error && get_option( argc, argv, options,
					sizeof(options) / sizeof( Option ),
					&optindex, &optname, &optarg, 
					&options_error )) {
    
    options_error = process_Gaze_Options( optname, optarg );
  }
  /* check that compulsory args were actually given */

  if (! options_error ) {
    if (gaze_options.structure_file == NULL) {
      fprintf( stderr, "You have not specified a structure file\n");
      options_error = TRUE;
    }
    else if (gaze_options.dna_file == NULL) {
      fprintf( stderr, "Warning: You have not specified a DNA file\n");
    }
    else if (gaze_options.gff_files->len == 0) {
      fprintf( stderr, "You have not specified a GFF feature file\n");
      options_error = TRUE;
    }
    else if (gaze_options.begin_dna > gaze_options.end_dna) {
      fprintf( stderr, "You have given an illegal DNA start/end range\n");
      options_error = TRUE;
    }
  }

  return (!options_error);
}


/*********************************************************************
                        MAIN
 *********************************************************************/

int main (int argc, char *argv[]) {

  Gaze_Structure *gs;
  Gaze_Output *g_out;
  GArray *features, *segments, *min_scores, *feature_path = NULL;
  Feature *beg_ft, *end_ft;
  char *dna_seq;
  int i,j,k, num_segs = 0;
  enum DP_Calc_Mode calc_mode;

  if (! parse_command_line(argc, argv) )
    exit(1);
  
  if(gaze_options.verbose)
    fprintf(stderr, "Parsing structure file\n");
  
  if ((gs = parse_Gaze_Structure( gaze_options.structure_file )) == NULL)
    exit(1);

  /** set up the output object */
  g_out = new_Gaze_Output();
  g_out->fh = gaze_options.out_file;
  g_out->posterior = gaze_options.posterior;
  g_out->use_threshold = gaze_options.use_threshold;
  g_out->threshold = gaze_options.threshold;
  /***/

  features = g_array_new( FALSE, TRUE, sizeof(Feature *));

  segments = g_array_new( FALSE, TRUE, sizeof(Segment_lists *));
  g_array_set_size( segments, gs->seg_dict->len );

  for(i=0; i < segments->len; i++) 
    g_array_index( segments, Segment_lists *, i) = new_Segment_lists(); 

  /* need to record the minimum score seen for each feature type, so 
     that we can obtain a score for features extracted from the DNA */

  min_scores = g_array_new( FALSE, TRUE, sizeof(double) );
  g_array_set_size( min_scores, gs->feat_dict->len );
  for(i=0; i < min_scores->len; i++)
    g_array_index( min_scores, double, i ) = 0.0;

  /******************************************************************/
  /* get all the features *******************************************/
  /******************************************************************/

  /* start by adding BEGIN and END by hand... */ 
  
  beg_ft = new_Feature();
  beg_ft->feat_idx = dict_lookup( gs->feat_dict, "BEGIN" );
  beg_ft->real_pos.s = gaze_options.begin_dna;
  beg_ft->real_pos.e = gaze_options.begin_dna;
  g_array_append_val( features, beg_ft );
  
  end_ft = new_Feature();
  end_ft->feat_idx = dict_lookup( gs->feat_dict, "END" );
  end_ft->real_pos.s = gaze_options.end_dna;
  end_ft->real_pos.e = gaze_options.end_dna;
  g_array_append_val( features, end_ft );

  /* Get the features from the GFF files... */
  
  if (gaze_options.verbose)
    fprintf(stderr, "Reading the gff files...\n");
  g_out->seq_name = get_features_from_gff( gaze_options.gff_files, features, segments, 
					   gs->gff_to_feats, min_scores, 
					   gaze_options.begin_dna, gaze_options.end_dna, 
					   gaze_options.use_selected ); 
  
  if (g_out->seq_name == NULL)
    exit(1);

  /* ...and from the DNA files */

  if (gaze_options.dna_file != NULL) {
    if (gaze_options.verbose)
      fprintf(stderr, "Reading the dna file...\n");
    dna_seq = read_dna_seq( gaze_options.dna_file, 
			    gaze_options.begin_dna, gaze_options.end_dna, gaze_options.offset_dna );
    get_features_from_dna( dna_seq, features, segments, 
			   gs->dna_to_feats, min_scores,
			   gaze_options.begin_dna); 
    if (gs->take_dna != NULL) {
      if (gaze_options.verbose)
	fprintf(stderr, "Getting dna for features...\n");
      get_dna_for_features( dna_seq, features, gs->take_dna, gs->motif_dict,
			    gaze_options.begin_dna, gaze_options.end_dna );
    }
    g_free( dna_seq );
  }
  
  g_array_free( min_scores, TRUE );


  /***********************************************/
  /* Scale, sort, and remove duplicates features */
  /***********************************************/

  if (gaze_options.verbose)
    fprintf(stderr, "Features: sorting, scaling %d feats and removing duplicates...", features->len );

  for( i=0; i < features->len; i++ ) {
    Feature *ft = g_array_index( features, Feature *, i );
    ft->score *= g_array_index( gs->feat_info, Feature_Info *, ft->feat_idx )->multiplier;
    ft->score *= gaze_options.sigma;

    /* for sorting, we need to calculate the effective "position" of each feature, 
       using the start_offset and end_offset */
    ft->adj_pos.s = ft->real_pos.s 
      + g_array_index( gs->feat_info, Feature_Info *, ft->feat_idx )->start_offset;

    ft->adj_pos.e = ft->real_pos.e 
      - g_array_index( gs->feat_info, Feature_Info *, ft->feat_idx )->end_offset;
  }

  qsort( features->data, features->len, sizeof(Feature *), &order_features_forwards); 
  features = remove_duplicate_features( features );

  if (gaze_options.verbose)
    fprintf(stderr, "%d features left\n", features->len);

  /***************************************/
  /* scale, sort and index segments ******/
  /***************************************/

  if (gaze_options.verbose)
    fprintf(stderr, "Segments: sorting, scaling, indexing and projecting...");

  for( i=0; i < segments->len; i++ ) {
    double multiplier = g_array_index( gs->seg_info, Segment_Info *, i )->multiplier;
    Segment_lists *seg_lists = g_array_index( segments, Segment_lists *, i);

    /* fourth element holds segments of this type in all frames */
    num_segs += g_array_index( seg_lists->orig, GArray *, 3 )->len;
    
    for (j=0; j < seg_lists->orig->len; j++) {
      GArray *o = g_array_index( seg_lists->orig, GArray *, j);
      GArray *p;
      
      for (k=0; k < o->len; k++) {
	Segment *seg = g_array_index( o, Segment *, k );
	seg->score *= multiplier;
	seg->score *= gaze_options.sigma;
	/* the following makes it a per-residue score */
	seg->score /= (seg->pos.e - seg->pos.s + 1);
      }
      
      qsort( o->data, o->len, sizeof(Segment *), &order_segments); 
      index_Segments( o );
      p = project_Segments( o );
      index_Segments( p );
      
      g_array_index( seg_lists->proj, GArray *, j) = p;
    }      
  }

  if (gaze_options.verbose)
    fprintf(stderr, " %d non-projected segs\n", num_segs);

  /******************************/
  /* Scale the length penalties */
  /******************************/
  for(i=0; i < gs->length_funcs->len; i++) {
    Length_Function *lf = g_array_index( gs->length_funcs, Length_Function *, i );
    for( j=0; j < lf->value_map->len; j++) {
      g_array_index( lf->value_map, double, j ) = 
	g_array_index( lf->value_map, double, j ) * lf->multiplier;
      g_array_index( lf->value_map, double, j ) = 
	g_array_index( lf->value_map, double, j ) * gaze_options.sigma;
    }
  }
  

  
  /************************************************************************
         Dynamic programming
  ************************************************************************/

  /* need to write the head first because forwards_calc produces 
     the output of all candidate regions, for space-saving reasons */
  write_GFF_header( g_out, gaze_options.begin_dna, gaze_options.end_dna );
  
  calc_mode = gaze_options.full_calc ? STANDARD_SUM : PRUNED_SUM;
  
  if (g_out->posterior) {
    if (gaze_options.verbose)
      fprintf(stderr, "Doing backward calculation...\n"); 
    backwards_calc( features, 
		    segments, 
		    gs, 
		    calc_mode);
  }
  
  if (gaze_options.verbose)
    fprintf(stderr, "Doing forward calculation...\n");
  
  forwards_calc( features, 
		 segments, 
		 gs, 
		 calc_mode, 
		 gaze_options.output == ALL_REGIONS ? g_out : NULL );
  
  if (gaze_options.gene_file != NULL) {

    /* read in and score the path */
    if (gaze_options.verbose)
      fprintf(stderr, "Reading the gff correct path file...\n");
    
    if ( (feature_path = read_in_path(gaze_options.gene_file, 
				      gs->feat_dict, 
				      features, 
				      stderr)) == NULL || 
	 ! is_legal_path( feature_path, gs, stderr ))
      exit(1);
    
    /* the following is called for its side effect of filling in path
       score up-to-and-including each feature in the path */ 
    calculate_path_score( feature_path, segments, gs );    
    print_GFF_path( g_out, feature_path, gs );      
  }
  else if (gaze_options.output == BEST_PATH || gaze_options.output == SAMPLE_PATH) {
      /* obtain a path by traceback */
      
      if (gaze_options.verbose)
	fprintf( stderr, "Tracing back...\n");
      feature_path = trace_back_general(features, 
					segments, 
					gs,
					MAX_TRACEBACK ); 
      
      print_GFF_path( g_out, feature_path, gs );        
  }
  else if (gaze_options.output == ALL_FEATURES) {
    /* output must be simple list. In one case, regions will have already been
       printed during the forwards calculation. In the other case, we need
       to print a feature list */
    
    /* before printing the posterior probabilities, re-sort the features in the standard
       way. The method of sorting used for the D.P. will not list the complete set of 
       features in an order that is intuitive */ 
    qsort( features->data, features->len, sizeof(Feature *), &order_features_standard); 
    print_GFF_Gaze_Features( g_out, features, gs );
  }


  free_Gaze_Output( g_out );
  for(i=0; i < features->len; i++)
    free_Feature( g_array_index( features, Feature *, i));
  g_array_free( features, TRUE);
  for(i=0; i < segments->len; i++)
    free_Segment_lists( g_array_index( segments, Segment_lists *, i ));
  g_array_free( segments, TRUE);
  if (feature_path != NULL)
    g_array_free (feature_path, TRUE );
  free_Gaze_Structure( gs );
  
  return 0;
}

