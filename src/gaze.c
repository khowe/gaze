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
\
 -structure_file <s>    XML file containing the gaze structure\n\
 -gff_file <s>          name of a GFF file containing the features (can give many)\n\
 -dna_file <s>          name of a DNA file in fasta format (can give many)\n\
 -gene_file <s>         name of a GFF file containing a user-specified gene structure\n\
 -id_file <s>           name of a file containing a list of name/start-ends to process\n\
 -defaults <s>          name of the file of default options (def: './gaze.defaults')\n\
\n\
Output files:\n\
 -out_file <s>       print gene structure to given file (def: stdout)\n\
\n\
Output format:\n\
 -output <c>            The output type. Should be one of the following:\n\
                          B  (B)est (highest scoring) gene structure (default)\n\
                          S  (S)ampled gene structure based on forward score\n\
                          F  (F)eature candidates (most sensibly used with -posterior)\n\
                          R  (R)egion candidates (for regions indicated in structure file)\n\
\n\
Other options:\n\
\
 -posterior             Show element scores as posteriror probs. rather than raw scores\n\
 -selected              look out for Selected features in input\n\
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
  { "-defaults_file", STRING_ARG },
  { "-output", CHAR_ARG },
  { "-help", NO_ARGS },
  { "-verbose", NO_ARGS },
  { "-posterior", NO_ARGS },
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

  Array *gff_file_names;  /* of string */
  Array *dna_file_names;  /* of string */
  Array *gene_file_names; /* of string */

  enum {
    BEST_PATH,
    SAMPLE_PATH,
    ALL_FEATURES,
    ALL_REGIONS 
  } output;    

  boolean full_calc;
  boolean use_selected;
  boolean verbose;
  boolean posterior;
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
	   strcmp(optname, "-gene_file") == 0 ) {
    Array *file_names, *files;

    if (strcmp(optname, "-gff_file") == 0)
      file_names = gaze_options.gff_file_names;
    else if (strcmp(optname, "-dna_file") == 0)
      file_names = gaze_options.dna_file_names;
    else if (strcmp(optname, "-gene_file") == 0)
      file_names = gaze_options.gene_file_names;
    
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
	char *temp = strdup( ln->buf );
	append_val_Array( nmstends, temp );
      }
      fclose(ids);
      free_Line(ln);
    }
    /* the from the command line */
    for (; optindex < argc; optindex++) {
      char *temp = strdup( argv[optindex] );
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
	  *ptr = '\0';
	  en = ptr+1;
	  
	    /* ideally, need a check here that start ends are digits, and to allow
	       the user to specify their own name/start-end format */
	    start = atoi( st );
	    end = atoi ( en );
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
 *********************************************************************
                        MAIN
 *********************************************************************
 *********************************************************************/
int main (int argc, char *argv[]) {

  int i,s = 0;

  if (! parse_command_line(argc, argv) )
    fatal_util( "use \"gaze -h\" to find out about usage");
  
  if(gaze_options.verbose)
    fprintf(stderr, "Parsing structure file\n");
  
  if ((gazeStructure = parse_Gaze_Structure( gaze_options.structure_file_name )) == NULL)
    exit(1);
	    

  allGazeSequences = new_Gaze_Sequence_list( gaze_options.sequence_names );

  /*******************************************************************/
  /* get the dna sequences *******************************************/
  /*******************************************************************/  

  for (s=0; s < allGazeSequences->num_seqs; s++)
    allGazeSequences->seq_list[s] = new_Gaze_Sequence( index_Array( gaze_options.sequence_names, char *, s),
						       index_Array( gaze_options.sequence_starts, int, s),
						       index_Array( gaze_options.sequence_ends, int, s) );
  if (gaze_options.verbose)
    fprintf(stderr, "Reading the dna files...\n");
  read_dna_seqs(allGazeSequences,
		gaze_options.dna_file_names );
    
  for (s=0; s < allGazeSequences->num_seqs; s++)
    initialise_Gaze_Sequence( allGazeSequences->seq_list[s], gazeStructure );

  /******************************************************************/
  /* First, obtain and set up all the Features and Segments *********/
  /******************************************************************/
    
  if (gaze_options.verbose)
    fprintf(stderr, "Reading the gff files...\n");
  get_features_from_gff( allGazeSequences,
			 gaze_options.gff_file_names,
			 gazeStructure->gff_to_feats, 
			 gaze_options.use_selected ); 

      
  if (gaze_options.verbose)
    fprintf(stderr, "Getting features from dna and dna for features...\n");
  
  for (i=0; i < allGazeSequences->num_seqs; i++) {
    if (allGazeSequences->seq_list[i]->dna_seq != NULL) {

      get_features_from_dna( allGazeSequences->seq_list[i], 
			     gazeStructure->dna_to_feats );

      if (gazeStructure->take_dna != NULL) {
	get_dna_for_features( allGazeSequences->seq_list[i], 
			      gazeStructure->take_dna, 
			      gazeStructure->motif_dict );

      }
      
      /* we never need the sequence itself again */
      free_util( allGazeSequences->seq_list[i]->dna_seq );
    }  
  }

  /******************************************************************/
  /** Obtain the given paths, if there are any **********************/
  /******************************************************************/

  if ( gaze_options.gene_file_names->len > 0) {
    if (gaze_options.verbose)
      fprintf(stderr, "Reading the gff correct path files...\n");
  
    if (! read_in_paths( allGazeSequences,
			 gaze_options.gene_file_names, 
			 gazeStructure->feat_dict ))
	fatal_util( "There was a problem reading in the correct paths\n" );

    /* check that any paths that were given are actually legal paths */
    for (i=0; i < allGazeSequences->num_seqs; i++) {
      Array *path = allGazeSequences->seq_list[i]->path;
      if (path != NULL && ! is_legal_path( path, gazeStructure ))
	fatal_util( "For sequence %s, the given \"correct\" path was illegal according to the model",
		     allGazeSequences->seq_list[i]->seq_name);
    }
  }

  /******************************/
  /* Scale the length penalties */
  /******************************/

  for(i=0; i < gazeStructure->length_funcs->len; i++) {
    Length_Function *lf = index_Array( gazeStructure->length_funcs, Length_Function *, i );
    scale_Length_Function( lf, lf->multiplier * gaze_options.sigma );
  }
  
  /***********************************************/
  /* Scale, sort, and remove duplicates features */
  /***********************************************/

  if (gaze_options.verbose)
    fprintf(stderr, "Sorting, scaling etc of features and segments...\n");
  
  for (s=0; s < allGazeSequences->num_seqs; s++) {
    Gaze_Sequence *g_seq = allGazeSequences->seq_list[s];

    for( i=0; i < g_seq->features->len; i++ ) {
      Feature *ft = index_Array( g_seq->features, Feature *, i );
      ft->score *= index_Array( gazeStructure->feat_info, Feature_Info *, ft->feat_idx )->multiplier;
      ft->score *= gaze_options.sigma;
      
      ft->adj_pos.s = ft->real_pos.s 
	+ index_Array( gazeStructure->feat_info, Feature_Info *, ft->feat_idx )->start_offset;
      
      ft->adj_pos.e = ft->real_pos.e 
	- index_Array( gazeStructure->feat_info, Feature_Info *, ft->feat_idx )->end_offset;
    }
    
    qsort( g_seq->features->data, g_seq->features->len, sizeof(Feature *), &order_features_for_dp); 
    remove_duplicate_features( g_seq );
    
    /***************************************/
    /* scale, sort and index segments ******/
    /***************************************/
    
    for( i=0; i < g_seq->segment_lists->len; i++ ) {
      Segment_list *seg_list = index_Array( g_seq->segment_lists, Segment_list *, i);
      double multiplier = index_Array( gazeStructure->seg_info, Segment_Info *, i )->multiplier;

      scale_Segment_list( seg_list, multiplier * gaze_options.sigma );
      sort_Segment_list ( seg_list );
      project_Segment_list( seg_list );
      index_Segment_list( seg_list );
    }
  }
    
  /************************************************************************/
  /* Finally, do the work                                                 */
  /************************************************************************/

  gazeOutput = new_Gaze_Output(gaze_options.out_file,
			       gaze_options.posterior,
			       gaze_options.use_threshold,
			       gaze_options.threshold);

  for (s=0; s < allGazeSequences->num_seqs; s++) {
    Gaze_Sequence *g_seq = allGazeSequences->seq_list[s];

    if(gaze_options.verbose)
      fprintf(stderr, "Running GAZE for sequence %s (%d-%d), %d feats\n", 
	      g_seq->seq_name, 
	      g_seq->seq_region.s, 
	      g_seq->seq_region.e,
	      g_seq->features->len);

    if (gazeOutput->posterior) {
      if (gaze_options.verbose)
	fprintf(stderr, "Doing backward calculation...\n"); 
      backwards_calc( g_seq,
		      gazeStructure, 
		      gaze_options.full_calc ? STANDARD_SUM : PRUNED_SUM );
    }
    
    if (gaze_options.verbose)
      fprintf(stderr, "Doing forward calculation...\n");


    /* need to write the head first because forwards_calc produces 
       the output of all candidate regions, for space-saving reasons */
    write_Gaze_header( gazeOutput, g_seq );
    
    forwards_calc( g_seq,
		   gazeStructure, 
		   gaze_options.full_calc ? STANDARD_SUM : PRUNED_SUM,
		   gaze_options.output == SAMPLE_PATH ? SAMPLE_TRACEBACK : MAX_TRACEBACK,
		   gaze_options.output == ALL_REGIONS ? gazeOutput : NULL );
  
    if (g_seq->path != NULL) {
    
      /* the following is called for its side effect of filling in path
	 score up-to-and-including each feature in the path */ 
      calculate_path_score( g_seq, gazeStructure );    
      write_Gaze_path( gazeOutput, g_seq, gazeStructure );      
    }
    else if (gaze_options.output == BEST_PATH || gaze_options.output == SAMPLE_PATH) {
      /* obtain a path by traceback */
      if (gaze_options.verbose)
	fprintf( stderr, "Tracing back...\n");
      trace_back_general(g_seq, 
			 gazeStructure ); 

      /* for non-strandard traceback, we need to recalculate the path score */
      calculate_path_score( g_seq, gazeStructure );
      write_Gaze_path( gazeOutput, g_seq, gazeStructure );        
    }
    else if (gaze_options.output == ALL_FEATURES) {
      /* before printing the posterior probabilities, re-sort the features in the standard
	 way. The method of sorting used for the D.P. will not list the complete set of 
	 features in an order that is intuitive */ 
      qsort( g_seq->features->data, g_seq->features->len, sizeof(Feature *), &order_features_standard); 
      write_Gaze_Features( gazeOutput, g_seq, gazeStructure );
    }
  }

  free_Gaze_Output( gazeOutput );
  free_Gaze_Structure( gazeStructure );
  free_Gaze_Sequence_list( allGazeSequences );
  
  return 0;
}

