/*  Last edited: Nov 22 16:20 2001 (klh) */
/**********************************************************************
 ** File: info.c
 ** Author : Kevin Howe
 ** E-mail : klh@sanger.ac.uk
 ** Description : 
 **********************************************************************/

#include "info.h"


/***************** Feature Info ************************/


/*********************************************************************
 FUNCTION: empty_Feature_Info
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Feature_Info *empty_Feature_Info(void) {
  Feature_Info *temp;

  temp = (Feature_Info *) g_malloc( sizeof(Feature_Info) );
  temp->start_offset = 0;
  temp->end_offset = 0;
  temp->is_killer_feat = FALSE;
  temp->sources = NULL;
  temp->targets = NULL;
  temp->kill_feat_quals_up = NULL;
  temp->kill_feat_quals_down = NULL;
  temp->seg_quals = NULL;

  return temp;
}


/*********************************************************************
 FUNCTION: free_Feature_Info
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Feature_Info(Feature_Info *ft_info) {
  int i;

  if (ft_info != NULL) {
    if (ft_info->sources != NULL) {
      for(i=0; i < ft_info->sources->len; i++) {
	free_Feature_Relation( g_array_index(ft_info->sources, Feature_Relation *, i) );
      }
      g_array_free( ft_info->sources, TRUE);
    }
    if (ft_info->targets != NULL) {
      for(i=0; i < ft_info->targets->len; i++) {
	free_Feature_Relation( g_array_index(ft_info->targets, Feature_Relation *, i) );
      }
      g_array_free( ft_info->targets, TRUE);
    }
    if (ft_info->kill_feat_quals_up != NULL) {
      for(i=0; i < ft_info->kill_feat_quals_up->len; i++) {
	if (g_array_index(ft_info->kill_feat_quals_up, Killer_Feature_Qualifier *, i)) {
	 free_Killer_Feature_Qualifier( g_array_index( ft_info->kill_feat_quals_up, Killer_Feature_Qualifier *, i));
	}
      }
      g_array_free( ft_info->kill_feat_quals_up, TRUE);
    }
    if (ft_info->kill_feat_quals_down != NULL) {
      for(i=0; i < ft_info->kill_feat_quals_down->len; i++) {
	if (g_array_index(ft_info->kill_feat_quals_down, Killer_Feature_Qualifier *, i)) {
	  free_Killer_Feature_Qualifier( g_array_index( ft_info->kill_feat_quals_down, Killer_Feature_Qualifier *, i));
	}
      }
      g_array_free( ft_info->kill_feat_quals_down, TRUE);
    }
    if (ft_info->seg_quals != NULL) {
      for(i=0; i < ft_info->seg_quals->len; i++) {
	if (g_array_index(ft_info->seg_quals, Segment_Qualifier *, i)) 
	  free_Segment_Qualifier( g_array_index( ft_info->seg_quals, Segment_Qualifier *, i));
      }
      g_array_free( ft_info->seg_quals, TRUE);
    }
  
    g_free( ft_info );
  }
}


/*********************************************************************
 FUNCTION: new_Feature_Info
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Feature_Info *new_Feature_Info(int start_off, int end_off, double multiplier) {
  Feature_Info *temp;

  temp = empty_Feature_Info();
  temp->multiplier = multiplier;
  temp->start_offset = start_off;
  temp->end_offset = end_off;

  return temp;
}


/****************** Feature_Relation ***********************/


/*********************************************************************
 FUNCTION: clone_Feature_Relation
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Feature_Relation *clone_Feature_Relation(Feature_Relation *src) {
  int i;
  Feature_Relation *dest = NULL;

  if (src != NULL) {
    dest = (Feature_Relation *) g_malloc( sizeof(Feature_Relation) );
    dest->target = src->target;
    dest->source = src->source;

    if (src->min_dist != NULL) {
      dest->min_dist = (int *) g_malloc( sizeof( int ) );
      *(dest->min_dist) = *(src->min_dist);
    }
    else
      dest->min_dist = NULL;

    if (src->max_dist != NULL) {
      dest->max_dist = (int *) g_malloc( sizeof( int ) );
      *(dest->max_dist) = *(src->max_dist);
    }
    else
      dest->max_dist = NULL;
    
    if (src->phase != NULL) {
      dest->phase = (int *) g_malloc( sizeof( int ) );
      *(dest->phase) = *(src->phase);
    }
    else
      dest->phase =NULL;
 
    if (src->len_fun != NULL) {
      dest->len_fun = (int *) g_malloc( sizeof( int ) );
      *(dest->len_fun) = *(src->len_fun);
    }
    else
      dest->len_fun =NULL;
 
   
    dest->out_feature = (src->out_feature != NULL)?g_strdup( src->out_feature ):NULL;
    dest->out_strand = (src->out_strand != NULL)?g_strdup( src->out_strand ):NULL;
    dest->out_frame = (src->out_frame != NULL)?g_strdup( src->out_frame ):NULL;
    
    if (src->seg_quals != NULL) {
      dest->seg_quals = g_array_new( FALSE, TRUE, sizeof(Segment_Qualifier *));
      g_array_set_size( dest->seg_quals, src->seg_quals->len );
      for (i=0; i < src->seg_quals->len; i++)
	g_array_index( dest->seg_quals, Segment_Qualifier *, i) =
	  clone_Segment_Qualifier( g_array_index(src->seg_quals,  Segment_Qualifier *, i) );
    }
    else
     dest->seg_quals = NULL;

    if (src->kill_feat_quals != NULL) {
      dest->kill_feat_quals = g_array_new( FALSE, TRUE, sizeof(Killer_Feature_Qualifier *));
      g_array_set_size( dest->kill_feat_quals, src->kill_feat_quals->len );
      for (i=0; i < src->kill_feat_quals->len; i++)
	g_array_index( dest->kill_feat_quals, Killer_Feature_Qualifier *, i) =
	  clone_Killer_Feature_Qualifier( g_array_index(src->kill_feat_quals, Killer_Feature_Qualifier *, i) );
    }
    else
     dest->kill_feat_quals = NULL;

    if (src->kill_dna_quals != NULL) {
      dest->kill_dna_quals = g_array_new( FALSE, TRUE, sizeof(Killer_DNA_Qualifier *));
      g_array_set_size( dest->kill_dna_quals, src->kill_dna_quals->len );
      for (i=0; i < src->kill_dna_quals->len; i++) 
	g_array_index( dest->kill_dna_quals, Killer_DNA_Qualifier *, i) = 
	  clone_Killer_DNA_Qualifier( g_array_index( src->kill_dna_quals, Killer_DNA_Qualifier *, i ) ); 
    }
    else
      dest->kill_dna_quals = NULL;

  }    

  return dest;
}




/*********************************************************************
 FUNCTION: free_Feature_Relation
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Feature_Relation(Feature_Relation *ft_src) {
  int i;

  if (ft_src != NULL) {
    if (ft_src->seg_quals != NULL) {
      for (i=0; i < ft_src->seg_quals->len; i++)
	  free_Segment_Qualifier( g_array_index( ft_src->seg_quals, Segment_Qualifier  *, i));
      g_array_free( ft_src->seg_quals, TRUE );
    }

    if (ft_src->kill_feat_quals != NULL) {
      for (i=0; i < ft_src->kill_feat_quals->len; i++)
	  free_Killer_Feature_Qualifier( g_array_index( ft_src->kill_feat_quals, Killer_Feature_Qualifier  *, i));
      g_array_free( ft_src->kill_feat_quals, TRUE );
    }

    if (ft_src->kill_dna_quals != NULL) {
      for (i=0; i < ft_src->kill_dna_quals->len; i++) 
	free_Killer_DNA_Qualifier( g_array_index( ft_src->kill_dna_quals, Killer_DNA_Qualifier *, i ) );
      g_array_free( ft_src->kill_dna_quals, TRUE );
    }

    if (ft_src->min_dist != NULL)
      g_free(ft_src->min_dist);
    if (ft_src->max_dist != NULL)
      g_free(ft_src->max_dist);
    if (ft_src->phase != NULL) 
      g_free(ft_src->phase);
    if (ft_src->len_fun != NULL)
      g_free( ft_src->len_fun );
    if (ft_src->out_feature != NULL)
      g_free( ft_src->out_feature );
    if (ft_src->out_strand != NULL)
      g_free( ft_src->out_strand );
    if (ft_src->out_frame != NULL)
      g_free( ft_src->out_frame );

    g_free( ft_src );
  }
}


/*********************************************************************
 FUNCTION: new_Feature_Relation
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Feature_Relation *new_Feature_Relation(void) {
  Feature_Relation *temp;

  temp = (Feature_Relation *) g_malloc( sizeof(Feature_Relation) );
  temp->target = 0;
  temp->source = 0;
  temp->min_dist = NULL;
  temp->max_dist = NULL;
  temp->phase = NULL;
  temp->seg_quals = NULL;
  temp->kill_feat_quals = NULL;
  temp->kill_dna_quals = NULL;
  temp->len_fun = NULL;
  temp->out_feature = NULL;
  temp->out_strand = NULL;
  temp->out_frame = NULL;
  
  return temp;
}




/************** Length_Function **********************/

/*********************************************************************
 FUNCTION: calc_Length_Function
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void calc_Length_Function(Length_Function *len_fun) {
  int last_x;
  int point_ctr, fctr;
  int x1, x2;
  double y1, y2, dy_by_dx;

  if (len_fun->raw_x_vals->len &&
      (len_fun->raw_x_vals->len == len_fun->raw_y_vals->len)) {

    last_x = g_array_index( len_fun->raw_x_vals,
			    int,
			    len_fun->raw_x_vals->len - 1);
    len_fun->value_map = g_array_new( FALSE, TRUE, sizeof( double ) );
    g_array_set_size( len_fun->value_map, last_x + 1 );
    

    point_ctr = 1;
    x1 = g_array_index( len_fun->raw_x_vals, int, point_ctr - 1);
    y1 = g_array_index( len_fun->raw_y_vals, double, point_ctr - 1);
	
    /* first, initialize all points up to x1 */

    for (fctr=0; fctr < x1; fctr++) {
      g_array_index( len_fun->value_map, double, fctr ) = 0.0;
    }
	
    while (fctr < last_x) {
      x2 = g_array_index( len_fun->raw_x_vals, int, point_ctr );
      y2 = g_array_index( len_fun->raw_y_vals, double, point_ctr );

      dy_by_dx = (y2 - y1) / ((double) x2 - (double) x1);
      while (fctr < x2) {
	g_array_index( len_fun->value_map, double, fctr ) = 
	  y1 + ((fctr - x1) * dy_by_dx);
	fctr++;
      }
	    
      x1 = x2;
      y1 = y2;
      point_ctr++;
    }

    /* finally, initialise the remaining point */

    g_array_index( len_fun->value_map, double, fctr ) = y1; 
      
  }
}



/*********************************************************************
 FUNCTION: free_Length_Function
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Length_Function(Length_Function *len_fun) {
  if (len_fun != NULL) {
    if (len_fun->value_map != NULL) {
      g_array_free( len_fun->value_map, TRUE );
    }
    if (len_fun->raw_x_vals != NULL) {
      g_array_free( len_fun->raw_x_vals, TRUE );
    }
    if (len_fun->raw_y_vals != NULL) {
      g_array_free( len_fun->raw_y_vals, TRUE );
    }
    
    g_free( len_fun );
  }
}



/*********************************************************************
 FUNCTION: new_Length_Function
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Length_Function *new_Length_Function( double multiplier) {
  Length_Function *temp;

  temp = (Length_Function *) g_malloc( sizeof(Length_Function));
  temp->value_map = NULL;
  temp->raw_x_vals = g_array_new( FALSE, TRUE, sizeof(int) );
  temp->raw_y_vals = g_array_new( FALSE, TRUE, sizeof(double) );
  temp->multiplier = multiplier;

  return temp;
}


/*********** Killer_DNA_Qualifier ********************/


/*********************************************************************
 FUNCTION: clone_Killer_DNA_Qualifier
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Killer_DNA_Qualifier *clone_Killer_DNA_Qualifier(Killer_DNA_Qualifier *src) {
  Killer_DNA_Qualifier *dest = NULL;

  if (src != NULL) {
    dest = new_Killer_DNA_Qualifier();
    dest->src_dna = src->src_dna;
    dest->tgt_dna = src->tgt_dna;
  }    

  return dest;
}


/*********************************************************************
 FUNCTION: free_Killer_DNA_Qualifier
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Killer_DNA_Qualifier(Killer_DNA_Qualifier *kill_qual) {
  if (kill_qual != NULL)
    g_free( kill_qual );
}



/*********************************************************************
 FUNCTION: new_Killer_DNA_Qualifier
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Killer_DNA_Qualifier *new_Killer_DNA_Qualifier(void) {
  Killer_DNA_Qualifier *temp;
  
  temp = (Killer_DNA_Qualifier *) g_malloc( sizeof(Killer_DNA_Qualifier) );
  temp->src_dna = -1;
  temp->tgt_dna = -1;
  
  return temp;
}





/*********** Killer_Feature_Qualifier ********************/

/*********************************************************************
 FUNCTION: clone_Killer_Feature_Qualifier
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Killer_Feature_Qualifier *clone_Killer_Feature_Qualifier(Killer_Feature_Qualifier *src) {
  Killer_Feature_Qualifier *dest = NULL;

  if (src != NULL) {
    dest = new_Killer_Feature_Qualifier();
    dest->has_phase = src->has_phase;
    dest->phase = src->phase;
  }    

  return dest;
}


/*********************************************************************
 FUNCTION: free_Killer_Feature_Qualifier
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Killer_Feature_Qualifier(Killer_Feature_Qualifier *kill_qual) {
  if (kill_qual != NULL)
    g_free( kill_qual );
}



/*********************************************************************
 FUNCTION: new_Killer_Feature_Qualifier
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Killer_Feature_Qualifier *new_Killer_Feature_Qualifier(void) {
  Killer_Feature_Qualifier *temp;
  
  temp = (Killer_Feature_Qualifier *) g_malloc( sizeof(Killer_Feature_Qualifier) );
  temp->has_phase = FALSE;
  temp->phase = 0;
  
  return temp;
}



/*********** Segment_Info ********************/


/*********************************************************************
 FUNCTION: empty_Segment_Info
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Segment_Info *empty_Segment_Info(void) {
  Segment_Info *temp;
  
  temp = (Segment_Info *) g_malloc( sizeof(Segment_Info) );

  return temp;
}




/*********************************************************************
 FUNCTION: free_Segment_Info
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Segment_Info(Segment_Info *seg_info) {
  if (seg_info != NULL)
    g_free( seg_info );
}



/*********************************************************************
 FUNCTION: new_Segment_Info
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Segment_Info *new_Segment_Info(double mul) {
  Segment_Info *temp;
  
  temp = empty_Segment_Info();
  temp->multiplier = mul;
  temp->use_projected = FALSE;
  temp->score_sum = FALSE;

  return temp;
}



/*********************************************************************
 FUNCTION: clone_Segment_Qualifier
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Segment_Qualifier *clone_Segment_Qualifier(Segment_Qualifier *src) {
  Segment_Qualifier *dest = NULL;

  if (src != NULL) {
    dest = new_Segment_Qualifier();
    dest->use_projected = src->use_projected;
    dest->score_sum = src->score_sum;
    dest->is_exact_src = src->is_exact_src;
    dest->is_exact_tgt = src->is_exact_tgt;
    dest->has_tgt_phase = src->has_tgt_phase;
    dest->has_src_phase = src->has_src_phase;
    dest->phase = src->phase;
  }    

  return dest;
}


/*********************************************************************
 FUNCTION: free_Segment_Qualifier
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Segment_Qualifier(Segment_Qualifier *seg_qual) {
  if (seg_qual != NULL)
    g_free( seg_qual );
}



/*********************************************************************
 FUNCTION: new_Segment_Qualifier
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Segment_Qualifier *new_Segment_Qualifier(void) {
  Segment_Qualifier *temp;
  
  temp = (Segment_Qualifier *) g_malloc( sizeof(Segment_Qualifier) );

  temp->use_projected = FALSE;
  temp->score_sum = FALSE;
  temp->is_exact_src = FALSE;
  temp->is_exact_tgt = FALSE;
  temp->has_tgt_phase = FALSE;
  temp->has_src_phase = FALSE;
  temp->phase = 0;
  
  return temp;
}
