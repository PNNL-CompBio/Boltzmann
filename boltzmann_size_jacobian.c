#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "boltzmann_size_jacobian.h"
void boltzmann_size_jacobian(struct state_struct *state) {
  /*
    Set the nnz, nnzm, nnzl, and nnzu fields of the cvodes_params struct
    from the reactions matrix.
    nnz is an upper bound for the number of nonzero entries in the jacobian.
    nnzm is the upper bound for the number of entries in the netwon 
    interation matrix, M = I - gamma J, ans is  <= nnz + ny,
    nnzl is the bound on the number of entries in the approximate lower
    triangular factor of M and nnzu is the similar bound on the approximate
    upper triangular factor of M
    These numbers are used in allocating space for working with the 
    jacobian in boltzmann_covdes.
    
    Called by: boltzmann_cvodes
  */
  struct cvodes_params_struct *cvodes_params;
  struct reactions_matrix_struct *reactions_matrix;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;

  int jacobian_choice;
  int ny;

  int nrxns;
  int i;

  int j;
  int rowlen;

  int nnz;
  int nnzm;

  int nnzl;
  int nnzu;

  int fill;
  int ny2;
  

  ny                = state->nunique_molecules;
  nrxns             = state->number_reactions;
  fill              = state->cvodes_prec_fill;
  cvodes_params     = state->cvodes_params;
  reactions_matrix  = state->reactions_matrix;
  rxn_ptrs          = reactions_matrix->rxn_ptrs;
  molecules_indices = reactions_matrix->molecules_indices;
  ny2               = ny * ny;
  /*
    loop over reactions.
  */
  nnz = 0;
  for (i=0;i<nrxns;i++) {
    rowlen = rxn_ptrs[i+1] - rxn_ptrs[i];
    nnz += (rowlen * rowlen);
  }
  if (nnz > ny2) {
    nnz = ny2; /* Dense Jacobian, might want to  set a flag here. */
  }
  /*
    I don't think we need to add ny to nnz for nnzm because every row
    will have a diagonal element.
  */
  nnzm = nnz; 
  /*
    Now it could be that all elements in M are below or above the diagonal
  */
  nnzl = nnzm + (fill * ny);
  nnzu = nnzm + (fill * ny);
  cvodes_params->nnz  = nnz;
  cvodes_params->nnzm = nnzm;
  cvodes_params->nnzl = nnzl;
  cvodes_params->nnzu = nnzu;
}
