#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "compute_flattened_state_size.h"
int64_t compute_flattened_state_size(struct state_struct *state) {
  /*
    Compute the size of the flattened state.
    Called by: boltzmann_flatten_alloc, boltzmann_boot
  */
  struct cvodes_params_struct cv_p_s;
  int64_t fsize;
  int nr;
  int nu;
  int nc;
  int nz;
  int cvps;
  cvps = (sizeof (cv_p_s) + 7) >> 3;
  /*
    Make cvps be a multiple of 8.
  */
  if ((cvps & 7) != 0) {
    cvps += (8-(cvps&7));
  }
  nr = state->number_reactions;
  nu = state->nunique_molecules;
  nc = state->nunique_compartments;
  nz = state->number_molecules;
  /*
                       words    ps
    Goblal meta data:  2        2

    External meta data  2       4 
    External data      10      14 (1 spare)

    Aux filename m.d.   2      16
    Aux filename data  16      32

    Scalars meta        2      34
      Scalar I8  m.d.   2      36
      Scalar I8  data  60      96 (4 spare)
      Scalar R8  md     2      98
      Scalar R8  data  38     136 (7 spare)

    Oneway Data meta    2     138
      rxns m.d          2     140
      rxns data        18*nr  140 + 18*nr

      molc m.d          2     142 + 18*nr
      molc data         4*nu  142 + 18*nr + 4*nu

      comp m.d.         2     144 + 18*nr + 4*nu
      comp data         8*nc  144 + 18*nr + 4*nu + 8*nc

      rmtrx m.d.        2     146 + 18*nr + 4*nu + 8*nc
      rmtrx data        1+3*nr + 4*nz
                              147 + 21*nr + 4*nu + 8*nc + 4*nz

      mmtrx m.d.	2     149 + 21*nr + 4*nu + 8*nc + 4*nz
      mmtrx data        1+nu+3*nz
                              150 + 21*nr + 5*nu + 8*nz + 7*nz

      1way vecs m.f.    2     152 + 21*nr + 5*nu + 8*nz + 7*nz
      1way_vecs data    13*nr + 7*nu      
                              152 + 34*nr + 12*nu + 8*nc + 7*nz

    Tweway m.d.         2     154 + 34*nr + 12*nu + 8*nc + 7*nz
      vgrng_state*2    32     186 + 34*nr + 12*nu + 8*nc + 7*nz
      cvodes_params    cvps (=49 I think)
      vectors          3*nu + nr
                              186 + 35*nr + 15*nu + 8*nc + 7*nz + cvps;

  */
  fsize = (int64_t)186 + (35*nr) + (15*nu) + (8*nc) + (7*nz) + cvps;
  return(fsize);
}  
