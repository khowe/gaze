/*  Last edited: Jul 24 14:06 2002 (klh) */
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

  temp = (Feature_Info *) malloc_util( sizeof(Feature_Info) );
  temp->start_offset = 0;
  temp->end_offset = 0;
  temp->is_killer_feat = FALSE;
  temp->sources = NULL;
  temp->targets = NULL;
  temp->seg_quals = NULL;
  temp->kill_feat_quals = NULL;
  temp->out_qual = NULL;

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
	free_Feature_Relation( index_Array(ft_info->sources, Feature_Relation *, i) );
      }
      free_Array( ft_info->sources, TRUE);
    }
    if (ft_info->targets != NULL) {
      for(i=0; i < ft_info->targets->len; i++) {
	free_Feature_Relation( index_Array(ft_info->targets, Feature_Relation *, i) );
      }
      free_Array( ft_info->targets, TRUE);
    }
    if (ft_info->kill_feat_quals != NULL) {
      for(i=0; i < ft_info->kill_feat_quals->len; i++) {
	if (index_Array(ft_info->kill_feat_quals, Killer_Feature_Qualifier *, i) != NULL) {
	 free_Killer_Feature_Qualifier( index_Array( ft_info->kill_feat_quals, Killer_Feature_Qualifier *, i));
	}
      }
      free_Array( ft_info->kill_feat_quals, TRUE);
    }
    if (ft_info->seg_quals != NULL) {
      for(i=0; i < ft_info->seg_quals->len; i++) {
	if (index_Array(ft_info->seg_quals, Segment_Qualifier *, i) != NULL) 
	  free_Segment_Qualifier( index_Array( ft_info->seg_quals, Segment_Qualifier *, i));
      }
      free_Array( ft_info->seg_quals, TRUE);
    }
  
    free_util( ft_info );
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
    dest = (Feature_Relation *) malloc_util( sizeof(Feature_Relation) );
    dest->target = src->target;
    dest->source = src->source;

    if (src->min_dist != NULL) {
      dest->min_dist = (int *) malloc_util( sizeof( int ) );
      *(dest->min_dist) = *(src->min_dist);
    }
    else
      dest->min_dist = NULL;

    if (src->max_dist != NULL) {
      dest->max_dist = (int *) malloc_util( sizeof( int ) );
      *(dest->max_dist) = *(src->max_dist);
    }
    else
      dest->max_dist = NULL;
    
    if (src->phase != NULL) {
      dest->phase = (int *) malloc_util( sizeof( int ) );
      *(dest->phase) = *(src->phase);
    }
    else
      dest->phase =NULL;
 
    if (src->len_fun != NULL) {
      dest->len_fun = (int *) malloc_util( sizeof( int ) );
      *(dest->len_fun) = *(src->len_fun);
    }
    else
      dest->len_fun =NULL;
 
    if (src->out_qual != NULL)
      dest->out_qual = clone_Output_Qualifier( src->out_qual );

    if (src->seg_quals != NULL) {
      dest->seg_quals = new_Array( sizeof(Segment_Qualifier *), TRUE);
      set_size_Array( dest->seg_quals, src->seg_quals->len );
      for (i=0; i < src->seg_quals->len; i++)
	index_Array( dest->seg_quals, Segment_Qualifier *, i) =
	  clone_Segment_Qualifier( index_Array(src->seg_quals,  Segment_Qualifier *, i) );
    }
    else
     dest->seg_quals = NULL;

    if (src->kill_feat_quals != NULL) {
      dest->kill_feat_quals = new_Array( sizeof(Killer_Feature_Qualifier *), TRUE);
      set_size_Array( dest->kill_feat_quals, src->kill_feat_quals->len );
      for (i=0; i < src->kill_feat_quals->len; i++)
	index_Array( dest->kill_feat_quals, Killer_Feature_Qualifier *, i) =
	  clone_Killer_Feature_Qualifier( index_Array(src->kill_feat_quals, Killer_Feature_Qualifier *, i) );
    }
    else
     dest->kill_feat_quals = NULL;

    if (src->kill_dna_quals != NULL) {
      dest->kill_dna_quals = new_Array( sizeof(Killer_DNA_Qualifier *), TRUE);
      set_size_Array( dest->kill_dna_quals, src->kill_dna_quals->len );
      for (i=0; i < src->kill_dna_quals->len; i++) 
	index_Array( dest->kill_dna_quals, Killer_DNA_Qualifier *, i) = 
	  clone_Killer_DNA_Qualifier( index_Array( src->kill_dna_quals, Killer_DNA_Qualifier *, i ) ); 
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
    if (ft_src->min_dist != NULL)
      free_util(ft_src->min_dist);
    if (ft_src->max_dist != NULL)
      free_util(ft_src->max_dist);
    if (ft_src->phase != NULL) 
      free_util(ft_src->phase);
    if (ft_src->len_fun != NULL)
      free_util( ft_src->len_fun );

    if (ft_src->seg_quals != NULL) {
      for (i=0; i < ft_src->seg_quals->len; i++)
	  free_Segment_Qualifier( index_Array( ft_src->seg_quals, Segment_Qualifier  *, i));
      free_Array( ft_src->seg_quals, TRUE );
    }

    if (ft_src->kill_feat_quals != NULL) {
      for (i=0; i < ft_src->kill_feat_quals->len; i++)
	  free_Killer_Feature_Qualifier( index_Array( ft_src->kill_feat_quals, Killer_Feature_Qualifier  *, i));
      free_Array( ft_src->kill_feat_quals, TRUE );
    }

    if (ft_src->kill_dna_quals != NULL) {
      for (i=0; i < ft_src->kill_dna_quals->len; i++) 
	free_Killer_DNA_Qualifier( index_Array( ft_src->kill_dna_quals, Killer_DNA_Qualifier *, i ) );
      free_Array( ft_src->kill_dna_quals, TRUE );
    }

    if (ft_src->out_qual != NULL)
      free_Output_Qualifier( ft_src->out_qual );

    free_util( ft_src );
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

  temp = (Feature_Relation *) malloc_util( sizeof(Feature_Relation) );
  temp->target = 0;
  temp->source = 0;
  temp->min_dist = NULL;
  temp->max_dist = NULL;
  temp->phase = NULL;
  temp->len_fun = NULL;
  temp->seg_quals = NULL;
  temp->kill_feat_quals = NULL;
  temp->kill_dna_quals = NULL;
  temp->out_qual = NULL;
  
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

    last_x = index_Array( len_fun->raw_x_vals,
			    int,
			    len_fun->raw_x_vals->len - 1);
    len_fun->value_map = new_Array( sizeof( double ), TRUE );
    set_size_Array( len_fun->value_map, last_x + 1 );
    

    point_ctr = 1;
    x1 = index_Array( len_fun->raw_x_vals, int, point_ctr - 1);
    y1 = index_Array( len_fun->raw_y_vals, double, point_ctr - 1);
	
    /* first, initialize all points up to x1 */

    for (fctr=0; fctr < x1; fctr++) {
      index_Array( len_fun->value_map, double, fctr ) = 0.0;
    }
	
    while (fctr < last_x) {
      x2 = index_Array( len_fun->raw_x_vals, int, point_ctr );
      y2 = index_Array( len_fun->raw_y_vals, double, point_ctr );

      dy_by_dx = (y2 - y1) / ((double) x2 - (double) x1);
      while (fctr < x2) {
	index_Array( len_fun->value_map, double, fctr ) = 
	  y1 + ((fctr - x1) * dy_by_dx);
	fctr++;
      }
	
      if (y2 >= y1) {
	if (! len_fun->becomes_monotonic ) {
	  len_fun->becomes_monotonic = TRUE;
	  len_fun->monotonic_point = x1;
	}
      }
      else 
	len_fun->becomes_monotonic = FALSE;
      
      x1 = x2;
      y1 = y2;
      point_ctr++;
    }

    /* finally, initialise the remaining point */

    index_Array( len_fun->value_map, double, fctr ) = y1; 
      
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
      free_Array( len_fun->value_map, TRUE );
    }
    if (len_fun->raw_x_vals != NULL) {
      free_Array( len_fun->raw_x_vals, TRUE );
    }
    if (len_fun->raw_y_vals != NULL) {
      free_Array( len_fun->raw_y_vals, TRUE );
    }
    
    free_util( len_fun );
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

  temp = (Length_Function *) malloc_util( sizeof(Length_Function));
  temp->value_map = NULL;
  temp->raw_x_vals = new_Array( sizeof(int), TRUE );
  temp->raw_y_vals = new_Array( sizeof(double), TRUE );
  temp->multiplier = multiplier;
  temp->becomes_monotonic = TRUE;
  temp->monotonic_point = 0;

  return temp;
}



/*********************************************************************
 FUNCTION: scale_Length_Function
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void scale_Length_Function(Length_Function *lf, double scaling) {
  int i;

  for( i=0; i < lf->value_map->len; i++) {
    index_Array( lf->value_map, double, i ) = 
      index_Array( lf->value_map, double, i ) * scaling;
  }

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
    free_util( kill_qual );
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
  
  temp = (Killer_DNA_Qualifier *) malloc_util( sizeof(Killer_DNA_Qualifier) );
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
    dest->feat_idx = src->feat_idx;
    dest->has_tgt_phase = src->has_tgt_phase;
    dest->has_src_phase = src->has_src_phase;
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
    free_util( kill_qual );
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
  
  temp = (Killer_Feature_Qualifier *) malloc_util( sizeof(Killer_Feature_Qualifier) );
  temp->feat_idx = 0;
  temp->has_tgt_phase = FALSE;
  temp->has_src_phase = FALSE;
  temp->phase = 0;
  
  return temp;
}


/*********** Output_Qualifier **********************************/

/*********************************************************************
 FUNCTION: clone_Output_Qualifier
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
Output_Qualifier *clone_Output_Qualifier( Output_Qualifier *src ) {
  Output_Qualifier *tgt = NULL;

  if (src != NULL) {
    if (src->feature != NULL)
      tgt->feature = strdup_util( src->feature );
    if (src->strand != NULL)
      tgt->strand = strdup_util( src->strand );
    if (src->frame != NULL) 
      tgt->frame = strdup_util( src->frame );
      
    tgt->need_to_print = src->need_to_print;
  }
  
  return tgt;
}


/*********************************************************************
 FUNCTION: free_Output_Qualifier
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/
void free_Output_Qualifier( Output_Qualifier *oq ) {
  if (oq != NULL) {
    if (oq->feature != NULL)
      free_util( oq->feature );
    if (oq->strand != NULL) 
      free_util( oq->strand );
    if (oq-> frame != NULL) 
      free_util( oq->frame );

    free_util( oq );
  }
}


/*********************************************************************
 FUNCTION: new_Output_Qualifier
 DESCRIPTION:
 RETURNS:
 ARGS: 
 NOTES:
 *********************************************************************/

Output_Qualifier *new_Output_Qualifier( void ) {
  Output_Qualifier *oq = (Output_Qualifier *) malloc_util( sizeof( Output_Qualifier ) );

  oq->feature = NULL;
  oq->strand = NULL;
  oq->frame = NULL;
  oq->need_to_print = FALSE;

  return oq;
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
  
  temp = (Segment_Info *) malloc_util( sizeof(Segment_Info) );

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
    free_util( seg_info );
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
  temp->partial = FALSE;

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
    dest->seg_idx = src->seg_idx;
    dest->use_projected = src->use_projected;
    dest->score_sum = src->score_sum;
    dest->is_exact_src = src->is_exact_src;
    dest->is_exact_tgt = src->is_exact_tgt;
    dest->has_tgt_phase = src->has_tgt_phase;
    dest->has_src_phase = src->has_src_phase;
    dest->phase = src->phase;
    dest->partial = src->partial;
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
    free_util( seg_qual );
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
  
  temp = (Segment_Qualifier *) malloc_util( sizeof(Segment_Qualifier) );

  temp->seg_idx = 0;
  temp->use_projected = FALSE;
  temp->score_sum = FALSE;
  temp->is_exact_src = FALSE;
  temp->is_exact_tgt = FALSE;
  temp->has_tgt_phase = FALSE;
  temp->has_src_phase = FALSE;
  temp->phase = 0;
  temp->partial = FALSE;
  
  return temp;
}
