#include "boltzmann_structs.h"
#include "boltzmann_flatten_alloc.h"
#include "boltzmann_flatten_external.h"
#include "boltzmann_flatten_aux_name.h"
#include "boltzmann_flatten_scalars.h"
#include "boltzmann_flatten_alloc1.h"
#include "boltzmann_flatten_oneway_data.h"
#include "boltzmann_flatten_twoway_data.h"
#include "boltzmann_flatten_state.h"
int boltzmann_flatten_state (struct state_struct **state_p,
			     void **fstate_p,
			     int direction, FILE *lfp) {
  /*
    Flatten a state struct into an allocated 
    flattened contiguous piece of memory if direction = 0 (save).

    allocate and load state_struct from a contiguous piece of memory
    if direction = 1 (load).

    While we could do this in separate routines, for the 
    sake of consistency and maintainability it is useful to do this in
    one routine, with the flatten and unflatten lines adjacent in the code.
    
    In order that the flattened state be somewhat self describing we use
    the notion of having sections each section starts with two eight byte
    words that are int64_t's.
    Word 0 is word count (words are 8 bytes) for the section including the
    2 header bytes.
    Word_1 is >0 if the section has subsections (each also with its own 
    2 word  header) or if word_1 is < 0,
       the section is composed 8 * word_0 bytes of type indicated by
n       word_1:  -1 doubles   (1 double/word)
                -2 int64_t's (8 bytes/int64_t  = 1 int64_t/word)
		-3 int's     (4 bytes/int  =  2 int's/word) 
		-4 char's    (8 chars/word)
		-5 compressed struct of 8 byte words of varying types -
		   basically just a memory block.

    Called by: boltzmann_save_state, boltzmann_load_state
    Calls:     boltzmann_flatten_alloc,
               boltzmann_flatten_external,
	       boltzmann_flatten_aux_name,
	       boltzmann_flatten_scalars,
	       boltzmann_flatten_alloc1,
	       boltzmann_flatten_oneway_data,
	       boltzmann_flatten_twoway_data,
	       fprintf,fflush

		
  */       
  struct  state_struct *state;
  void    *fstate;
  int64_t *lfstate;
  int64_t num_top_elements;
  int64_t oneway_pos;
  int64_t twoway_pos;

  int success;
  int total_size;
  int word_pos;



  /*
    Allocate flattened_state (fstate)  or state depending on direction
    being 0 or 1.
  */
  success = boltzmann_flatten_alloc(state_p,
				    fstate_p,
				    direction, &word_pos, lfp);
  if (success) {
    state = (struct state_struct *)*state_p;
    fstate = (void *)*fstate_p;
    lfstate = (int64_t *)fstate;
    /* 
      Transfer external constants.
      if (direction == 0) fstate <= state, elxe state<-fstate
      
    */
    success = boltzmann_flatten_external(state,fstate,direction,
					 &word_pos,lfp);
  }
  if (success) {
    /*
      Store or retrieve Name of auxiliary data file.
    */
    success = boltzmann_flatten_aux_name(state,fstate,direction,&word_pos,
					 lfp);
  }
  /*
    At this point word_pos should be 31.
    Flatten or load the scalar fields of state.
  */
  if (success) {
    success = boltzmann_flatten_scalars(state,fstate,
					direction,&word_pos,lfp);
  } 
  /*
    At this point word_pos should be 127.
  */
  /*
    Next if direction == 1 we need to allocate all the nonscalar
    data fields of state, and the work fields as well.
  */
  if (direction > 0) {
    success = boltzmann_flatten_alloc1(state);
  }
  /*
    Flatten or load the one way fields of state.
  */
  if (success) {
    oneway_pos = word_pos;
    success = boltzmann_flatten_oneway_data(state,fstate,direction,
					    &word_pos,lfp);
  }
  if (success) {
    twoway_pos = word_pos;
    success = boltzmann_flatten_twoway_data(state,fstate,direction,
					    &word_pos,lfp);
  }
  if (success) {
    total_size = word_pos + 1;
    if (direction == 0) {
      lfstate[0] = total_size;
    } else {
      if (lfstate[0] != total_size) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"boltzmann_flatten_state: state sizes do not match\n");
	  fprintf(lfp,"        lfstate[0] = %ld, total_size = %d\n",
		  lfstate[0],total_size);
	  fflush(lfp);
	}
      }
    }
  }
  return(success);
}
