/*  Last edited: Jul 22 11:42 2002 (klh) */
/**********************************************************************
 ** FILE: util.c
 ** NOTES:
 **   This file contains general utiliy functions used throughout the 
 **   application, such as those for memory management and error
 **   messaging. I have used it as a place to put other general stuff
 **   until I have somewhere better to put it.
 **********************************************************************/

#include "util.h"


/*********************************************************************/
/********************** Debug functions ******************************/
/*********************************************************************/

static long int total_bytes;

long int how_many_bytes (void) {
  return total_bytes;
}



/*********************************************************************/
/********************** Error reporting ******************************/
/*********************************************************************/


/********************************************************************* 
 FUNCTION: fatal_util
 DESCRIPTION: 
   Prints the given formatted error message and exits
 RETURNS:
 ARGS:
   A format string + args, c.f. printf
 NOTES:
 *********************************************************************/

void fatal_util( char *fmt, ... ) {
  va_list args;
  
  va_start( args, fmt );
  fprintf( stderr, "\nA Fatal Error occurred: ");
  vfprintf( stderr, fmt, args);
  fprintf( stderr,"\n");
  va_end( args );
  exit(1);
}



/********************************************************************* 
 FUNCTION: warning_util
 DESCRIPTION: 
   Prints the given formatted warning to stderr
 RETURNS:
 ARGS:
   A format string + args, c.f. printf
 NOTES:
 *********************************************************************/

void warning_util( char *fmt, ...) {
  va_list args;
	
  va_start( args, fmt );
  fprintf( stderr, "\nWARNING: " );
  vfprintf( stderr, fmt, args );
  fprintf( stderr, "\n" );
  va_end( args );
}




/*********************************************************************/
/********************** String functions *****************************/
/*********************************************************************/

#define LINE_ALLOC_STEP 50

/********************************************************************* 
 FUNCTION: new_Line
 DESCRIPTION: 
    Duplicates the given string and returns it
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
Line *new_Line( void ) {
  Line *l = (Line *) malloc_util( sizeof( Line ) );

  l->buf = (char *) malloc_util( LINE_ALLOC_STEP * sizeof( char ) );
  l->buf_size = LINE_ALLOC_STEP;

  return l;

}

/********************************************************************* 
 FUNCTION: free_Line
 DESCRIPTION: 
    Duplicates the given string and returns it
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
void free_Line( Line *l ) {
  if (l != NULL) {
    if (l->buf != NULL) 
      free_util( l->buf );

    free_util(l);
  }
}


/********************************************************************* 
 FUNCTION: read_Line
 DESCRIPTION: 
    Duplicates the given string and returns it
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
int read_Line( FILE *file, Line *ln ) {

  boolean got_line = FALSE;
  int c, line_len = 0;

  while(!got_line && (c = fgetc(file)) != EOF) {

    if (c == '\n')
      /* blank line, no good */
      continue;

    do {
      if (line_len >= ln->buf_size) {
	ln->buf = (char *) realloc_util( ln->buf, (line_len + LINE_ALLOC_STEP) * sizeof(char));
	ln->buf_size += LINE_ALLOC_STEP;
      }
      ln->buf[line_len++] = c;

    } while((c = fgetc( file)) != '\n');

    if (line_len >= ln->buf_size) {
      ln->buf = (char *) realloc_util( ln->buf, (line_len + 1) * sizeof(char) );
      ln->buf_size++;
    }
    ln->buf[line_len] = '\0';
    got_line = TRUE;
  }

  return line_len;

}


/********************************************************************* 
 FUNCTION: strdup_util
 DESCRIPTION: 
    Duplicates the given string and returns it
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
char *strdup_util(const char *str) {
  char *new_str;

  if (str) {
    new_str = (char *) malloc_util ( strlen (str) + 1 );
    strcpy (new_str, str);
  }
  else
    new_str = NULL;
  
  return new_str;
}




/**********************************************************************/
/*************** Core memory allocation wrappers **********************/
/**********************************************************************/


void *calloc_util( size_t numobjs, size_t size) {
  void *ret;
	
  if ((ret = calloc( numobjs, size )) == NULL)
    fatal_util("calloc_util: Out of memory");

  return ret;	
}



void *malloc0_util(size_t numbytes) {
  void *ret;

  if ((ret = malloc( numbytes )) == NULL)
    fatal_util("malloc_util: out of memory when requesting %d bytes", numbytes);
  memset( ret, 0, numbytes );

  return ret;	
}

void *malloc_util(size_t numbytes) {
  void *ret;

  if ((ret = malloc( numbytes )) == NULL)
    fatal_util("malloc_util: out of memory when requesting %d bytes", numbytes);

  return ret;	
}



void *realloc_util(void *ptr, size_t bytes) {
  void *ret = NULL;

  if (ptr == NULL)
    fatal_util("Call to realloc_util with a null pointer");
  else { 
    if ((ret = realloc(ptr, bytes)) == NULL)
      fatal_util("realloc_util: out of memory when requesting %d bytes", bytes);
  }
  return ret;
}  


void free_util( void *ptr ) {
  if (ptr == NULL)
    fatal_util("Call to free_util with null pointer");
  else {
    free(ptr);
    ptr = NULL;
  }
}



/**********************************************************************/
/*************** dynamically growable arrays **************************/
/**********************************************************************/


#define MIN_ARRAY_SIZE  16

struct _RealArray
{
  char *data;
  int   len;
  int   alloc;
  int   elt_size;
  int   clear;
};

typedef struct _RealArray  RealArray;

static int _next_power_of_two (int num);
static void _array_expand_if_necessary (RealArray *array,
				  int len);

/********************************************************************* 
 FUNCTION: new_Array
 DESCRIPTION: 
   Creates a new dynamically growable array, with no elements
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
Array *new_Array ( int elt_size,
		   boolean clear ) {
  RealArray *array;

  array = (RealArray *) malloc_util ( sizeof (RealArray) );
  array->data            = NULL;
  array->len             = 0;
  array->alloc           = 0;
  array->clear           = (clear ? 1 : 0);
  array->elt_size        = elt_size;

  return (Array*) array;
}



/********************************************************************* 
 FUNCTION: free_Array
 DESCRIPTION: 
   Frees the dynamically growable array, including the data
   segment if the second argument is true.
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
void free_Array (Array  *array, boolean free_segment) {
  if (free_segment && array->data != NULL)
    free_util (array->data);

  free_util( array );
}



/********************************************************************* 
 FUNCTION: append_vals_Array
 DESCRIPTION: 
    Copys "len" elements from the given data block, to the end 
    of the given array, resizing if necessary
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
Array* append_vals_Array (Array *farray,
			  const void *data,
			  int len) {

  RealArray *array = (RealArray*) farray;

  _array_expand_if_necessary (array, len);
  memcpy (array->data + array->elt_size * array->len, data, array->elt_size * len);
  array->len += len;

  return farray;
}


/********************************************************************* 
 FUNCTION: prepend_vals_Array
 DESCRIPTION: 
    Copys "len" elements from the given data block, to the start
    of the given array, resizing if necessary
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
Array* prepend_vals_Array (Array *farray,
			   const void *data,
			   int len) {

  RealArray *array = (RealArray*) farray;

  _array_expand_if_necessary (array, len);
  memmove (array->data + array->elt_size * len, array->data, array->elt_size * array->len);
  memcpy (array->data, data, len * array->elt_size);
  array->len += len;

  return farray;
}


/********************************************************************* 
 FUNCTION: insert_vals_Array
 DESCRIPTION: 
    Copys "len" elements from the given data block, to "index"
    of the given array, resizing if necessary
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
Array* insert_vals_Array (Array *farray,
			  int index,
			  const void *data,
			  int len) {

  RealArray *array = (RealArray*) farray;

  _array_expand_if_necessary (array, len);
  memmove (array->data + array->elt_size * (len + index), 
	     array->data + array->elt_size * index, 
	     array->elt_size * (array->len - index));
  memcpy (array->data + array->elt_size * index, data, len * array->elt_size);
  array->len += len;

  return farray;
}


/********************************************************************* 
 FUNCTION: set_size_Array
 DESCRIPTION: 
   Called after new_Array if the array size is known, to set the
   size of the array
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
Array *set_size_Array (Array *farray,
		       int length) {

  RealArray *array = (RealArray*) farray;

  if (array->len != 0)
    warning_util("Attempt to set the size of an array that is already populated\n");
  else {
    array->alloc = MAX ( MIN_ARRAY_SIZE, length * array->elt_size );

    array->data = malloc_util( array->alloc );

    if (array->clear)
      memset (array->data, 0, array->alloc );
  }

  array->len = length;

  return farray;
}



/********************************************************************* 
 FUNCTION: remove_index_Array
 DESCRIPTION: 
    Removes the element from the array at the given index, shiftng
    the others along
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
Array *remove_index_Array (Array* farray,
			   int index) {

  RealArray* array = (RealArray*) farray;

  if (index != array->len - 1)
      memmove (array->data + array->elt_size * index, 
	       array->data + array->elt_size * (index + 1), 
	       array->elt_size * (array->len - index - 1));
  
  array->len -= 1;

  return farray;
}



/********************************************************************* 
 FUNCTION: _next_power_of_two
 DESCRIPTION: 
    Returns the nearest power of two to the given number
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
static int _next_power_of_two (int num)
{
  int n = 1;

  while (n < num)
    n <<= 1;

  return n;
}


/********************************************************************* 
 FUNCTION: _array_expand_if_necessary
 DESCRIPTION: 
   Reallocates the given array to the given size, if the given 
   size is bigger than the current size.
 RETURNS:
 ARGS:
 NOTES:
 *********************************************************************/
static void _array_expand_if_necessary (RealArray *array,
					int len)
{
  int want_alloc = (array->len + len) * array->elt_size;

  if (want_alloc > array->alloc) {
    int old_alloc = array->alloc;

    array->alloc = _next_power_of_two(want_alloc);
    array->alloc = MAX (array->alloc, MIN_ARRAY_SIZE);

    if (old_alloc == 0) {
      array->data = malloc_util (array->alloc);
    }
    else {
      array->data = realloc_util (array->data, array->alloc);
    }
    
    if (array->clear)
      memset (array->data + old_alloc, 0, array->alloc - old_alloc);
  }
}



/**********************************************************************/
/*************** Dictionaries *****************************************/
/**********************************************************************/


/*********************************************************************
 FUNCTION: dict_lookup
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
signed char dict_lookup( Dict *dict, const char *name ) {
  boolean match = FALSE;
  signed char j;

  for( j=0; j < ((Array *)dict)->len; j++ ) { 
    if (! strcmp( index_Array( (Array *)dict, char *, j), name )) {
      match = TRUE;
      break;
    }
  }

  return (match)?(signed char)j:-1;
}


