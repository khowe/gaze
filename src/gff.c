/*  Last edited: Jul 26 15:59 2002 (klh) */
/**********************************************************************
 ** File: gff.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#include "gff.h"


/*********************************************************************
 FUNCTION: free_gff_line
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_GFF_line( GFF_line *gff_line ) {
  if (gff_line != NULL) {
    if (gff_line->ln != NULL)
      free_Line( gff_line->ln );
    free_util(gff_line);
  }
}



/*********************************************************************
 FUNCTION: new_gff_line
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
GFF_line *new_GFF_line( void ) {
  GFF_line *line = (GFF_line *) malloc_util( sizeof( GFF_line ));

  line->ln = new_Line();
  line->seqname = NULL;
  line->source = NULL;
  line->type = NULL;
  line->start = line->end = 0;
  line->score = 0.0;
  line->strand = NULL;
  line->frame = NULL;
  line->group = NULL;

  return line;
}


/*********************************************************************
 FUNCTION: read_GFF_line
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
int read_GFF_line(FILE *file,
		  GFF_line *line) {

  char *fields[10];
  int line_len, i, j;
  boolean got_line = FALSE;

  while( ! got_line && (line_len = read_Line( file, line->ln )) != 0)  {

    if (line->ln->buf[0] != '#') {
      fields[0] = line->ln->buf;
      /* nullify 9th field, because it might not be present - first 8 will be */
      fields[8] = NULL;

      for (i=0, j=1; i < line_len && j < 9; i++) {
	if (line->ln->buf[i] == '\t') {
	  line->ln->buf[i] = '\0';
	  fields[j++] = &(line->ln->buf[i+1]);
	}
      }

      line->seqname = fields[0];
      line->source = fields[1];
      line->type = fields[2];
      line->start = atoi(fields[3]);
      line->end = atoi(fields[4]);
      line->score = atof(fields[5]);
      line->strand = fields[6];
      line->frame = fields[7];
      line->group = fields[8];

      got_line = TRUE;
    }
  }

  return line_len;
}



/*********************************************************************
 FUNCTION: write_GFF_header
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void write_GFF_comment( FILE *fh,
			char *fmt,
			...) {
  va_list args;

  va_start( args, fmt );
  fprintf( fh, "##" );
  vfprintf( fh, fmt, args);
  fprintf( fh,"\n");
  va_end( args );
}



/*********************************************************************
 FUNCTION: write_GFF_line
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void write_GFF_line( FILE *fh,
		     char *seq,
		     char *source,
		     char *feature,
		     int start,
		     int end,
		     double score,
		     char *strand,
		     char *frame,
		     char *group) {

  fprintf( fh, "%s\t%s\t%s\t%d\t%d\t%.4f\t%s\t%s", 
	   seq != NULL ? seq : "Not_given",
	   source != NULL ? source : "Not_given",
	   feature != NULL ? feature : "Not_given",
	   start,
	   end,
	   score,
	   strand != NULL ? strand : ".",
	   frame != NULL ? frame : "." );
  if ( group != NULL ? group : "" )
    fprintf( fh, "\t%s\n", group );
  fprintf( fh, "\n");

}


/*********************************************************************
 FUNCTION: write_GFF_header
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void write_GFF_header( FILE *fh,
		       char *name, 
		       int start,
		       int end ) {

  fprintf( fh, "##gff-version 2\n");
  fprintf( fh, "##sequence-region %s %d %d\n", name, start, end );
}
