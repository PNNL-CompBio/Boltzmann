#include "boltzmann_structs.h"
#include "copy_local_states.h"
int copy_local_states(struct boot_state_struct *boot_state) {
  /*
    Copy the local states into the global(super) state.
    Called by: write_super_state;
    Calls:     fseek, read, write, fprintf, fflush, fclose
  */
  int64_t nw;
  int64_t ntr;
  int64_t nr;
  int64_t read_pos;
  int64_t sum_local_state_lengths;
  int64_t page_size;
  int64_t page_buffer_length;
  int64_t io_buff_size_in_pages;
  int64_t lseek_pos;
  char    *page_buffer;

  int boot_work_fd;
  int global_state_fd;
  int success;
  int whence;

  FILE *lfp;
  FILE *boot_work_fp;

  success                 = 1;
  sum_local_state_lengths = boot_state->work_offset;
  io_buff_size_in_pages   = boot_state->io_buff_size_in_pages;
  page_size               = boot_state->page_size;
  page_buffer             = boot_state->io_buff;
  boot_work_fd            = boot_state->boot_work_fd;
  global_state_fd         = boot_state->global_state_fd; 
  lfp                     = boot_state->lfp;
  boot_work_fp            = boot_state->boot_work_fp;
  page_buffer_length      = io_buff_size_in_pages * page_size;
  whence     = SEEK_SET;
  lseek_pos = (int64_t)0;
  lseek(boot_work_fd,lseek_pos,whence);
  /*
    Now copy the local states into the global state file.
    Do this in largest chunks as possible.
  */
  for (read_pos=(int64_t)0;((read_pos<sum_local_state_lengths) && success);
       read_pos+=page_buffer_length) {
    ntr = page_buffer_length;
    if (read_pos + ntr >= sum_local_state_lengths) {
      ntr = sum_local_state_lengths - read_pos;
    }
    nr = read(boot_work_fd,page_buffer,ntr);
    if (nr != ntr) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"copy_local_state_to_super_state: Error reading temporary local state file\n");
	fflush(lfp);
      }
    }
    if (success) {
      nw = write(global_state_fd,page_buffer,ntr);
      if (nw != ntr) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"write_super_state: Error writing global state file\n");
	  fflush(lfp);
	}
      }
    }
  } /* end for (read_pos...) */
  if (boot_work_fp) {
    fclose(boot_work_fp);
  }
  return(success);
}
