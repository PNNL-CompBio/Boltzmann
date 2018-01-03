#include "boltzmann_structs.h"

#include "boot_alloc0.h"
int boot_alloc0(struct boot_state_struct **boot_state_p) {
  /*
    Allocate space for the boot_state_struct.
    Called by: boot_init
    Calls:     calloc,fprintf,fflush
    Returns 1 on success, 0 on failure.
  */
  struct boot_state_struct *boot_state;
  struct boot_state_struct bss;
  int64_t ask_for;
  int64_t one_l;
  int success;
  int padi;
  success = 1;
  one_l   = (int64_t)1;
  ask_for = (int64_t)sizeof(bss);
  boot_state = (struct boot_state_struct *)calloc(one_l,ask_for);
  if (boot_state == NULL) {
    success = 0;
    fprintf(stderr,"boot_alloc0: Error unable to allocate %lld bytes for boot_state\n",ask_for);
    fflush(stderr);
  }
  *boot_state_p = boot_state;
  return(success);
}
