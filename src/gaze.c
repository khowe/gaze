/*  Last edited: Aug  3 15:55 2002 (klh) */
/**********************************************************************
 ** File: gaze.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "options.h"
#include "info.h"
#include "str_parse.h"
#include "g_engine.h"
#include "output.h"
#include "sequence.h"


static Gaze_Structure *gazeStructure;
static Gaze_Output *gazeOutput;
static Gaze_Sequence_list *allGazeSequences;


/*********************************************************************/
/***** Processing the command-line ***********************************/
/*********************************************************************/

static char gaze_usage_string[] = "\
Usage: gaze <options> seq_name1/start-end seq_name2/start-end ... \n\
Options are:\n\
\n\
Input files:\n\
\n\
 -structure_file <s>    XML file containing the gaze structure\n\
 -gff_file <s>          name of a GFF file containing the features (can give many)\n\
 -dna_file <s>          name of a DNA file in fasta format (can give many)\n\
 -gene_file <s>         name of a GFF file containing a user-specified gene structure\n\
 -selected_file <s>     name of file containing features that should be Selected/Not selected\n\
 -id_file <s>           name of a file containing a list of name/start-ends to process\n\
 -defaults <s>          name of the file of default options (def: './gaze.defaults')\n\
\n\
Output files:\n\
 -out_file <s>          print gene structure to given file (def: stdout)\n\
\n\
Output options:\n\
 -sample_gene           use posterior sampling to obtain gene structure (default: max)\n\
 -regions               show region candidates (for regions indicated in structure file)\n\
 -feature               show feature candidates (most sensibly used with -probability)\n\
    (If non of the above is suuplied, the highest scoring gene structure is obtained)\n\
\n\
Other options:\n\
\
 -probability           show output scores as posterior probabilities\n\
 -full_calc             perform full dynamic programming (as opposed to faster heurstic method)\n\
 -verbose               write basic progess information to stderr\n\
 -help                  show this message\n" ;

static Option options[] = {
  { "-dna_file", STRING_ARG },
  { "-structure_file", STRING_ARG },
  { "-gff_file", STRING_ARG },
  { "-out_file", STRING_ARG },
  { "-id_file", STRING_ARG },
  { "-gene_file", STRING_ARG },
  { "-selected_file", STRING_ARG },
  { "-defaults_file", STRING_ARG },
  { "-sample_gene", NO_ARGS },
  { "-regions", NO_ARGS },
  { "-features", NO_ARGS },
  { "-help", NO_ARGS },
  { "-verbose", NO_ARGS },
  { "-probability", NO_ARGS },
  { "-full_calc", NO_ARGS },
  { "-sigma", FLOAT_ARG }
};


static struct {
  char *structure_file_name;
  char *id_file_name;
  char *out_file_name;
  FILE *out_file;

  Array *sequence_names;   /* of char */
  Array *sequence_starts;  /* of int  */ 
  Array *sequence_ends;    /* of int */

  Array *gff_file_names;       /* of string */
  Array *dna_file_names;       /* of string */
  Array *gene_file_names;      /* of string */
  Array *selected_file_names;  /* of string */

  boolean sample_gene;
  boolean output_regions;
  boolean output_features;

  boolean full_calc;
  boolean use_selected;
  boolean verbose;
  boolean probability;
  boolean use_threshold;

  double threshold;
  double sigma;

} gaze_options;


/*********************************************************************
 FUNCTION: process_Gaze_Options
 *********************************************************************/
static boolean process_Gaze_Options(char *optname,
				     char *optarg ) {
  boolean options_error = FALSE;

  if (strcmp(optname, "-sigma") == 0) gaze_options.sigma = atof( optarg );
  else if (strcmp(optname, "-selected") == 0) gaze_options.use_selected = TRUE;	     
  else if (strcmp(optname, "-verbose") == 0) gaze_options.verbose = TRUE;
  else if (strcmp(optname, "-probability") == 0) gaze_options.probability = TRUE;  
  else if (strcmp(optname, "-full_calc") == 0) gaze_options.full_calc = TRUE;
  else if (strcmp(optname, "-sample_gene") == 0) gaze_options.sample_gene = TRUE;
  else if (strcmp(optname, "-regions") == 0) gaze_options.output_regions = TRUE;
  else if (strcmp(optname, "-features") == 0) gaze_options.output_features = TRUE;
  else if (strcmp(optname, "-threshold") == 0) {
    gaze_options.use_threshold = TRUE;
    gaze_options.threshold = atof( optarg );
  }
  else if (strcmp(optname, "-structure_file") == 0) {
    FILE *test = fopen( optarg, "r");
    if (test == NULL) {
      fprintf( stderr, "Could not open structure file %s for reading\n", optarg );
      options_error = TRUE;
    }
    else {
      fclose(test);
      if (gaze_options.structure_file_name != NULL)
	free_util( gaze_options.structure_file_name );
      gaze_options.structure_file_name = strdup_util( optarg );
    }
  }
  else if (strcmp(optname, "-id_file") == 0) {
    FILE *test = fopen( optarg, "r");
    if (test == NULL) {
      fprintf( stderr, "Could not open ID file %s for reading\n", optarg );
      options_error = TRUE;
    }
    else {
      fclose(test);
      if (gaze_options.id_file_name != NULL)
	free_util( gaze_options.id_file_name );
      gaze_options.id_file_name = strdup_util( optarg );
    }
  }
  else if (strcmp(optname, "-out_file") == 0) {
    if ((gaze_options.out_file = fopen( optarg, "w")) == NULL) {
      fprintf( stderr, "Could not open output file %s for writing\n", optarg );
      options_error = TRUE;
    }
    else {
      if (gaze_options.out_file_name != NULL)
	free_util( gaze_options.out_file_name );
      gaze_options.out_file_name = strdup_util( optarg );
    }
  }

  else if (strcmp(optname, "-gff_file") == 0 || 
	   strcmp(optname, "-dna_file") == 0 ||
	   strcmp(optname, "-selected_file") == 0 ||
	   strcmp(optname, "-gene_file") == 0 ) {
    Array *file_names = NULL;

    if (strcmp(optname, "-gff_file") == 0)
      file_names = gaze_options.gff_file_names;
    else if (strcmp(optname, "-dna_file") == 0)
      file_names = gaze_options.dna_file_names;
    else if (strcmp(optname, "-gene_file") == 0)
      file_names = gaze_options.gene_file_names;
    else if (strcmp(optname, "-selected_file") == 0)
      file_names = gaze_options.selected_file_names;

    if (dict_lookup( file_names, optarg ) >= 0)
      /* need to warn about duplicate feature file */
	fprintf( stderr, "Warning: file %s was given more than once\n", optarg );
    else {
      /* Try to open the file, and if success, store both the file name
	 and the file handle */
      FILE *test_f = fopen( optarg, "r");
      if (test_f == NULL) {
	fprintf( stderr, "Could not open file %s for reading\n", optarg );
	options_error = TRUE;
      }
      else {
	char *tmp_f_name = strdup_util( optarg );
	append_val_Array( file_names, tmp_f_name );
	fclose( test_f );
      }
    }
  }
  else if (strcmp(optname, "-help") == 0 ||
	   strcmp(optname, "-defaults_file") == 0 ) {
    /* ignore these options because they have already been dealt with */

  }
  else {
    fprintf(stderr, "Unrecognised option: %s\n", optname );
    options_error = TRUE;
  }

  return options_error;
}



/*********************************************************************
 FUNCTION: parse_command_line
 *********************************************************************/
static int parse_command_line( int argc, char *argv[] ) {
  boolean options_error = FALSE;
  boolean help_wanted = FALSE;
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

  gaze_options.sigma = 1.0;

  gaze_options.sequence_names = new_Array (sizeof( char * ), TRUE);
  gaze_options.sequence_starts = new_Array (sizeof( int ), TRUE);
  gaze_options.sequence_ends = new_Array (sizeof( int ), TRUE);

  gaze_options.structure_file_name = NULL;
  gaze_options.out_file_name = strdup_util( "stdout ");
  gaze_options.out_file = stdout;

  gaze_options.dna_file_names = new_Array( sizeof( char *), TRUE );
  gaze_options.gff_file_names = new_Array( sizeof( char *), TRUE );
  gaze_options.gene_file_names = new_Array( sizeof( char *), TRUE );
  gaze_options.selected_file_names = new_Array( sizeof( char *), TRUE );

  gaze_options.use_selected = FALSE;
  gaze_options.verbose = FALSE;
  gaze_options.full_calc = FALSE;
  gaze_options.probability = FALSE;
  gaze_options.use_threshold = FALSE;
  gaze_options.threshold = 0.0;

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
    if (gaze_options.structure_file_name == NULL) {
      fprintf( stderr, "Error: You have not specified a structure file\n");
      options_error = TRUE;
    }
    if (gaze_options.gff_file_names->len == 0) {
      fprintf( stderr, "Warning: You have not specified any GFF files\n");
      options_error = TRUE;
    }
    if (gaze_options.dna_file_names->len == 0) {
      fprintf( stderr, "Warning: You have given not any DNA files\n");
      options_error = TRUE;
    }
  }

  /* check for legal combinations of options */

  if (gaze_options.sample_gene || gaze_options.output_regions || gaze_options.output_features) {
    if ( (gaze_options.sample_gene && gaze_options.output_regions) ||
	 (gaze_options.output_regions && gaze_options.output_features) ||
	 (gaze_options.sample_gene && gaze_options.output_features) ) {
      
      fprintf( stderr, "Error: You can only give one of -sample_gene -regions and -features\n");
      options_error = TRUE;
    }
  }

  /* the list of file names to process is the combination of thise names
     in the -id_file file, and those left on the command-line */

  if (! options_error) {
    /* first, form list of names */
    Array *nmstends = new_Array( sizeof( char *), TRUE );
    boolean dup_found;
    int nmidx;

    /* first from the file */
    if (gaze_options.id_file_name != NULL) {
      FILE *ids = fopen( gaze_options.id_file_name, "r" );
      Line *ln = new_Line();
      while( read_Line( ids, ln ) > 0 ) {
	if (ln->buf[0] != '#') { 
	  char *temp = strdup_util( ln->buf );
	  append_val_Array( nmstends, temp );
	}
      }
      fclose(ids);
      free_Line(ln);
    }
    /* the from the command line */
    for (; optindex < argc; optindex++) {
      char *temp = strdup_util( argv[optindex] );
      append_val_Array( nmstends, temp );
    }

    for (nmidx = 0; nmidx < nmstends->len; nmidx++) {
      char *st, *en, *ptr;
      int start = 0;
      int i, end = 0;

      char *seq_id = index_Array( nmstends, char *, nmidx );
      
      if ((ptr = strchr( seq_id, '/')) != NULL) {
	*ptr = '\0';
	st = ptr+1;
	if ((ptr = strchr( st, '-')) != NULL) {
	  boolean all_digits = TRUE;

	  *ptr = '\0';
	  en = ptr+1;

	  /* ideally, allow the user to specify their own name/start-end format */
	  for (i=0; all_digits && i < strlen( st ); i++)
	    if (! isdigit((int)st[i])) 
	      all_digits = FALSE;
	  for (i=0; all_digits && i < strlen( en ); i++)
	    if (! isdigit((int)en[i])) 
	      all_digits = FALSE;

	  if (all_digits) {
	    start = atoi( st );
	    end = atoi ( en );

	    if (start > end) {
	      fprintf(stderr, "For %s, start (%d) is greater than end (%d)\n", seq_id, start, end);
	      options_error = TRUE;
	    }
	  }
	  else {
	    fprintf(stderr, "For %s, start-end of \"%s-%s\" is illegal\n", seq_id, st, en);
	    options_error = TRUE;
	  }
	}
      }
      
      dup_found = FALSE;
      for (i=0; i < gaze_options.sequence_names->len; i++)
	if (! strcmp( seq_id, index_Array( gaze_options.sequence_names, char *, i))) {
	  fprintf( stderr, "Warning: sequence %s given more than one - ignoring\n", seq_id );
	  dup_found = TRUE;
	  break;
	}
      
      if (! dup_found) {	    
	append_val_Array( gaze_options.sequence_names, seq_id );
	append_val_Array( gaze_options.sequence_starts, start );
	append_val_Array( gaze_options.sequence_ends, end );
      }
    }
    
    free_Array( nmstends, TRUE );
  }
  
  if ( gaze_options.sequence_names->len == 0 ) {
    fprintf( stderr, "Error: you have not give any sequences for GAZE to work on\n");
    options_error = TRUE;
  } 

  return (!options_error);
}




/*********************************************************************
 FUNCTION: prepare_Gaze_Sequence_for_work
    This function basically fills in the Sequence by reading the
    GFF and dna files, and then sorting and scaling the features etc.

 *********************************************************************/
static void prepare_Gaze_Sequence_for_work( Gaze_Sequence *g_seq ) {
  int i;

  /*******************************************************************/
  /* get the dna sequences *******************************************/
  /*******************************************************************/  
    
  if (gaze_options.verbose)
    fprintf(stderr, "Getting for DNA for %s...\n", g_seq->seq_name);
  
  read_dna_Gaze_Sequence( g_seq,
			  gaze_options.dna_file_names );
  
  /* sequences are intialised after reading the DNA, just in case
     the user did not supply start-end information in which case
     we have to serive it from the DNA */
  initialise_Gaze_Sequence( g_seq, gazeStructure );
  
  /******************************************************************/
  /* First, obtain and set up all the Features and Segments *********/
  /******************************************************************/
  
  if (gaze_options.verbose)
    fprintf(stderr, "Reading the gff files...\n");
  convert_gff_Gaze_Sequence( g_seq,
			       gaze_options.gff_file_names,
			     gazeStructure->gff_to_feats ); 
  
  if (gaze_options.verbose)
    fprintf(stderr, "Getting features from dna...\n");
  
  if (g_seq->dna_seq != NULL) {
    convert_dna_Gaze_Sequence( g_seq,
			       gazeStructure->dna_to_feats,
			       gazeStructure->take_dna, 
			       gazeStructure->motif_dict );
    
    /* we never need the sequence itself again */
    free_util( g_seq->dna_seq );
    g_seq->dna_seq = NULL;
  } 
  
  /******************************************************************/
  /*** Sorting and scaling of feature and segments ******************/
  /******************************************************************/
  
  if (gaze_options.verbose)
    fprintf(stderr, "Sorting, and scaling of features and segments...\n");
  
  /* first the features */
  qsort( g_seq->features->data, g_seq->features->len, sizeof(Feature *), &order_features); 
  remove_duplicate_features( g_seq );     
  
  for( i=0; i < g_seq->features->len; i++ ) {
    Feature *ft = index_Array( g_seq->features, Feature *, i );
    ft->score *= index_Array( gazeStructure->feat_info, Feature_Info *, ft->feat_idx )->multiplier;
      ft->score *= gaze_options.sigma;
      
      ft->adj_pos.s = ft->real_pos.s 
	+ index_Array( gazeStructure->feat_info, Feature_Info *, ft->feat_idx )->start_offset;
      
      ft->adj_pos.e = ft->real_pos.e 
	- index_Array( gazeStructure->feat_info, Feature_Info *, ft->feat_idx )->end_offset;
  }
    
  /* now the segments..*/
  for( i=0; i < g_seq->segment_lists->len; i++ ) {
    Segment_list *seg_list = index_Array( g_seq->segment_lists, Segment_list *, i);
    double multiplier = index_Array( gazeStructure->seg_info, Segment_Info *, i )->multiplier;
    
    scale_Segment_list( seg_list, multiplier * gaze_options.sigma );
    sort_Segment_list ( seg_list );
    project_Segment_list( seg_list );
    index_Segment_list( seg_list );
  }
  
  /******************************************************************/
  /** Obtain the given paths, if there are any **********************/
  /******************************************************************/
  
  if ( gaze_options.gene_file_names->len > 0) {
    if (gaze_options.verbose)
      fprintf(stderr, "Reading the gff correct path files...\n");
    
    if (! get_correct_feats_Gaze_Sequence( g_seq,
					     gaze_options.gene_file_names, 
					   gazeStructure->feat_dict, 
					   TRUE))
      fatal_util( "There was a problem reading in the correct paths\n" );
    
    /* check that any paths that were given are actually legal paths */
    
    if (g_seq->path != NULL && ! is_legal_path( g_seq->path, gazeStructure ))
      fatal_util( "For sequence %s, the given \"correct\" path was illegal according to the model",
		  g_seq->seq_name);
  }

  
  if ( gaze_options.selected_file_names->len > 0) {
    if (gaze_options.verbose)
      fprintf(stderr, "Reading selected feature files...\n");
    
    if (! get_correct_feats_Gaze_Sequence( g_seq,
					   gaze_options.selected_file_names, 
					   gazeStructure->feat_dict, 
					   FALSE))
      fatal_util( "There was a problem reading in the selected features\n" );
  }
}


/*********************************************************************
 FUNCTION: cleanup_Gaze_Sequence_after_work
    This function frees the parts of the Gaze_Sequence that
    are not needed any more

 *********************************************************************/
static void cleanup_Gaze_Sequence_after_work ( Gaze_Sequence *g_seq ) {

  free_Gaze_Sequence( g_seq, FALSE );
}



/*********************************************************************
 *********************************************************************
                        MAIN
 *********************************************************************
 *********************************************************************/
int main (int argc, char *argv[]) {

  int i = 0;
  Gaze_Sequence *g_seq;

  if (! parse_command_line(argc, argv) )
    fatal_util( "use \"gaze -h\" to find out about usage");
  
  if(gaze_options.verbose)
    fprintf(stderr, "Parsing structure file\n");
  
  if ((gazeStructure = parse_Gaze_Structure( gaze_options.structure_file_name )) == NULL)
    exit(1);
	    
  /******************************/
  /* Scale the length penalties */
  /******************************/
  for(i=0; i < gazeStructure->length_funcs->len; i++) {
    Length_Function *lf = index_Array( gazeStructure->length_funcs, Length_Function *, i );
    scale_Length_Function( lf, lf->multiplier * gaze_options.sigma );
  }

  gazeOutput = new_Gaze_Output(gaze_options.out_file,
			       gaze_options.probability,
			       gaze_options.sample_gene,
			       gaze_options.output_features,
			       gaze_options.output_regions,
			       gaze_options.use_threshold,
			       gaze_options.threshold);

  allGazeSequences = new_Gaze_Sequence_list( gaze_options.sequence_names );
  for (i=0; i < allGazeSequences->num_seqs; i++) 
    allGazeSequences->seq_list[i] = new_Gaze_Sequence( index_Array( gaze_options.sequence_names, char *, i),
						       index_Array( gaze_options.sequence_starts, int, i),
						       index_Array( gaze_options.sequence_ends, int, i) );

  for (i=0; i < allGazeSequences->num_seqs; i++) {
    g_seq = allGazeSequences->seq_list[i];

    prepare_Gaze_Sequence_for_work ( g_seq );
      
    if(gaze_options.verbose)
      fprintf(stderr, "Running GAZE for sequence %s (%d-%d), %d feats\n", 
	      g_seq->seq_name, 
	      g_seq->seq_region.s, 
	      g_seq->seq_region.e,
	      g_seq->features->len);

    if (gazeOutput->probability) {
      if (gaze_options.verbose)
	fprintf(stderr, "Doing backward calculation...\n"); 
      backwards_calc( g_seq,
		      gazeStructure, 
		      ! gaze_options.full_calc );
    }
    
    if (gaze_options.verbose)
      fprintf(stderr, "Doing forward calculation...\n");

    /* need to write the head first because forwards_calc produces 
       the output of all candidate regions, for space-saving reasons */
    write_Gaze_header( gazeOutput, g_seq );
    
    forwards_calc( g_seq,
		   gazeStructure, 
		   ! gaze_options.full_calc,
		   gazeOutput );
  
    if (gaze_options.output_features)
       write_Gaze_Features( gazeOutput, g_seq, gazeStructure );
    else if (!gaze_options.output_regions) {
      if (g_seq->path == NULL) {
	if (gaze_options.verbose)
	  fprintf( stderr, "Tracing back...\n");
	trace_back_general(g_seq );
      }
      
      calculate_path_score( g_seq, gazeStructure );
      write_Gaze_path( gazeOutput, g_seq, gazeStructure );
    }

    cleanup_Gaze_Sequence_after_work( g_seq );
  }

  free_Gaze_Output( gazeOutput );
  free_Gaze_Structure( gazeStructure );
  free_Gaze_Sequence_list( allGazeSequences );
  
  return 0;
}

