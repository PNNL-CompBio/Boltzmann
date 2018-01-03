#include "boltzmann_structs.h"
#include "catenate_compartments_and_molecules.h"
int catenate_compartments_and_molecules(struct boot_state_struct *boot_state) {
  /*
    Read in and catenate the molecules and compartments texts and meta data
    from the boot_work file.
    Called by: boltzmann_boot
    Calls      pread, fprinf, fflush, strerror
  */
  struct state_struct local_state_instance;
  struct state_struct *local_state;
  struct compartment_struct compartment_instance;
  struct compartment_struct *compartment_sort_ws;
  struct compartment_struct *compartments;
  struct molecule_struct molecule_instance;
  struct molecule_struct *molecule_sort_ws;
  struct molecule_struct *molecules;
  int64_t *offset_size_pairs;
  int64_t *compartment_list_starts;
  int64_t *molecule_map_starts;
  char    *compartment_text_ws;
  char    *molecule_text_ws;
  int64_t state_size;
  int64_t molecule_text_pos;
  int64_t compartment_text_pos;
  int64_t i;
  int64_t j;
  int64_t nr;
  int64_t offset;
  int64_t num_reaction_files;
  int64_t lseek_pos;
  int64_t compartment_meta_length;
  int64_t molecule_meta_length;
  int64_t compartment_text_length;
  int64_t molecule_text_length;
  int64_t nunique_molecules;
  int64_t nunique_compartments;
  int64_t mi_base;
  int64_t ci_base;
  int success;
  int boot_work_fd;
  int save_errno;
  int padi;
  FILE *lfp;
  FILE *efp;
  success            	  = 1;
  lfp                	  = boot_state->lfp;
  offset_size_pairs  	  = boot_state->state_offsets_sizes;
  num_reaction_files 	  = boot_state->num_reaction_files;
  boot_work_fd       	  = boot_state->boot_work_fd;
  compartment_sort_ws     = boot_state->compartment_sort_ws;
  molecule_sort_ws        = boot_state->molecule_sort_ws;
  compartment_list_starts = boot_state->compartment_list_starts;
  molecule_map_starts     = boot_state->molecule_map_starts;
  local_state             = &local_state_instance;
  state_size              = (int64_t)sizeof(local_state_instance);
  molecule_text_ws        = boot_state->molecule_text_ws;
  compartment_text_ws     = boot_state->compartment_text_ws;
  molecule_text_pos       = (int64_t)0;
  compartment_text_pos    = (int64_t)0;
  for (i=0;((i<num_reaction_files) && success);i++) {
    /*
      Seek to the position of the metadata for reaction file i.
    */
    offset = offset_size_pairs[i+i];
    /*
      lseek(tmp_state_fd,offset,whence);
      get the local state pointers.
    */
    nr = pread(boot_work_fd,local_state,state_size,offset);
    save_errno=errno;
    if (nr != state_size) {
      success =0;
      if (lfp) {
	fprintf(lfp,"catenate_compartments_and_molecules: Error i = %lld, could not read "
		"%lld bytes for a local_state, errno= %d, %s\n",i,state_size,
		save_errno,strerror(save_errno));
	fflush(lfp);
	/*print some error message.*/
      }
    }
    if (success) {
      /*
	read in the compartment text 
      */
      compartment_text_length = local_state->compartment_text_length;
      lseek_pos = offset + local_state->compartments_text_offset_in_bytes;
      /*
      lseek(boot_work_fd,lseek_pos,whence);
      */
      nr = pread(boot_work_fd,(char *)&compartment_text_ws[compartment_text_pos],
		 compartment_text_length,lseek_pos);
      save_errno = errno;
      if (nr != compartment_text_length) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"catenate_compartments_and_molecules: "
		  "Error i = %lld, could not read "
		  "%lld bytes for compartment text, errno = %d,%s\n",
		  i,compartment_text_length,save_errno,strerror(save_errno));
	  fflush(lfp);
	}
      }
    }
    if (success) {
      /*
	read in molecules text.
      */
      molecule_text_length    = local_state->molecule_text_length;
      lseek_pos = offset + local_state->molecules_text_offset_in_bytes;
      /*
      lseek(boot_work_fd,lseek_pos,whence);
      */
      nr = pread(boot_work_fd,(char *)&molecule_text_ws[molecule_text_pos],
		 molecule_text_length,lseek_pos);
      save_errno = errno;
      if (nr != molecule_text_length) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"catenate_compartments_and_molecules: Error i = %lld,"
		  " could not read %lld bytes for molecules text, "
		  "errno = %d,%s\n",i,molecule_text_length,save_errno,
		  strerror(save_errno));
	  fflush(lfp);
	}
      }
    }
    if (success) {
      /*
	read in the meta data  for the compartments.
      */
      nunique_compartments = local_state->nunique_compartments;
      compartment_meta_length = nunique_compartments * sizeof(compartment_instance);
      lseek_pos = offset + local_state->sorted_compartments_offset_in_bytes;
      /*
	lseek(boot_work_fd,lseek_pos,whence);
      */
      ci_base = compartment_list_starts[i];
      compartments = (struct compartment_struct *)&compartment_sort_ws[ci_base];
      nr = pread(boot_work_fd,(void *)compartments,compartment_meta_length,lseek_pos);
      save_errno = errno;
      if (nr != compartment_text_length) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"catenate_compartments_and_molecules: Error rxn, i = %lld,"
		  " could not read %lld bytes for sorted compartments meta "
		  "data, errno = %d,%s\n",i,compartment_text_length,
		  save_errno,strerror(save_errno));
	  fflush(lfp);
	} 
      } else {
	/*
	  Set up global_compartment string pointers, from the
	  local_compartment string field.
	*/
	for (j=0;j<nunique_compartments;j++) {
	  compartments->c_index = ci_base + j;
	  compartments->g_index = i;
	  compartments->string += compartment_text_pos;
	  compartments++; /* Caution address arithmetic. */
	}
	compartment_text_pos += compartment_text_length;
      }
    }
    if (success) {
      /*
	read in the metadata for the molecules;
      */
      nunique_molecules = local_state->nunique_molecules;
      molecule_meta_length = nunique_molecules * sizeof(molecule_instance);
      lseek_pos = offset + local_state->sorted_molecules_offset_in_bytes;
      /*
	lseek(boot_work_fd,lseek_pos,whence);
      */
      mi_base = molecule_map_starts[i];
      molecules = (struct molecule_struct *)&molecule_sort_ws[mi_base];
      nr = pread(boot_work_fd,
		 (void *)molecules,molecule_meta_length,lseek_pos);
      save_errno = errno;
      if (nr != molecule_text_length) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"catenate_compartments_and_molecules: Error i = %lld,"
		  " could not read %lld bytes for sorted molecules meta data, "
		  " errno = %d, %s\n",i,molecule_text_length,
		  save_errno,strerror(save_errno));
	  fflush(lfp);
	}
      } else {
	/*
	  Set up global_molecule string pointers, from the
	  local_molecule string field.
	*/
	for (j=0;j<nunique_molecules;j++) {
	  molecules->m_index = mi_base + j;
	  molecules->string  += molecule_text_pos;
	  molecules->g_index = i;
	  molecules++; /* Caution address arithmetic. */
	}
	molecule_text_pos += molecule_text_length;
      }
    }
    offset_size_pairs += 2;   /* Caution address arithmetic. */
  } /* end for(i...) */
  return(success);
}
