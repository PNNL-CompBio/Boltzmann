#include "boltzmann_structs.h"

#include "boot_alloc0.h"
#include "boot_alloc1.h"
#include "boot_io_init.h"
#include "size_rxns_list.h"
#include "boot_alloc2.h"
#include "boot_init.h"

int boot_init(struct state_struct *local_state, struct boot_state_struct **boot_state_p) {
  /*
    Allocate a boot_state_structure and its fields, 
    Allocate the fixed size of the superstate struct,
    open the output files,
    size the reactions_list file,
    and allocate the per_reaction file pointers.

    Called by: boltzmann_boot
    Calls:     boot_alloc0,
               boot_alloc1,
	       boot_io_ini,
	       size_rxns_list,
	       boot_alloc2
    
  */
  struct boot_state_struct *boot_state;
  int64_t page_size;
  int64_t page_mask;
  int64_t one_l;
  int64_t zero_l;
  int64_t filename_len;
  int64_t rxn_list_buffer_len;
  int success;
  int log2_page_size;

  zero_l  = (int64_t)0;
  one_l   = (int64_t)1;
  success = boot_alloc0(boot_state_p);
  if (success) {
    boot_state     = *boot_state_p;
    log2_page_size = 12;
    page_size      = one_l << log2_page_size;
    page_mask      = page_size - one_l;
    filename_len   = (int64_t)1024;
    rxn_list_buffer_len             = (int64_t)1024;
    boot_state->page_size 	    = page_size;
    boot_state->page_mask 	    = page_mask;
    boot_state->filename_len        = filename_len;
    boot_state->rxn_list_buffer_len = rxn_list_buffer_len;
    boot_state->align_len           = local_state->align_len;
    boot_state->align_mask          = local_state->align_mask;
    success = boot_alloc1(boot_state);
  }
  if (success) {
    success = boot_io_init(local_state,boot_state);
  }
  if (success) {
    success = size_rxns_list(boot_state);
  }
  if (success) {
    success = boot_alloc2(boot_state);
  }
  if (success) {
    /*
      Initialize the global counters.
    */
    boot_state->work_offset             = zero_l;
    boot_state->molecule_text_length    = zero_l;
    boot_state->compartment_text_length = zero_l;
    boot_state->unique_molecules_max    = zero_l;
    boot_state->unique_compartments_max = zero_l;
    boot_state->io_buff_size_in_pages   = (int64_t)16;
  }
  return(success);
}
