#include "boltzmann_structs.h"

/*
#include "compute_flux_scaling.h"
*/
#include "update_rxn_likelihoods.h"

#include "lr6_approximate_delta_concs.h"
int lr6_approximate_delta_concs(struct state_struct *state, 
				double *concs,
				double *flux, 
				int choice) {
  /*
    Compute approximations to concentration changes wrt time, 
    Based on likelihood computational alchemy
    for computing ratios of forward reaction rates.

    so kf_rel[i] = 


       Kf[i]/Kf[base] = exp ((sum(nu[j,base]*mu0[j]) - sum(nu[k,i]*mu0[k]))/ rt)
                             for species j              for species k
			     a reactant in              a reactant in 
			     reaction base              reaction i

   nu[m,n] = stoichiometric coefficient for species m in reaction n
   mu0[m]  = delta g0 of solvation for species m

   We will want vectors of length num_reactions to store kf_rel and
   rxn_quotient vectors.


    reverse reaction likelihoods are the same as the base reaction.
    Get reference from Bill Cannon
    Called by: approximate_delta_concs
    Calls:     

                                TMF
    state                       *SI   Boltzmant state structure.
                                      uses number_reactions,
				           unique_moleules,
                                           molecules_matrix,
					   and lfp,
				      

    concs			D1I   molecule concentrations vector of length 
                                      nunique_moleucles

    flux                        D1O   vector of length  unique_molecules
                                      of concentration change per unit time.
				      Set by this routine.

    choice                      IOI   Not used by this routine.

  */
  struct  molecule_struct *molecules;
  struct  molecule_struct *molecule;
  struct  molecules_matrix_struct *molecules_matrix;
  struct  rxn_matrix_struct *rxn_matrix;
  /*
  double  *forward_rxn_likelihoods;
  double  *reverse_rxn_likelihoods; 
  double  *counts;
  double  *reactant_term;
  
  double  *conc_to_count;
  */
  double  *c;
  double  *log_kf_rel;
  double  *log_kr_rel;
  double  *ke;
  double  *rke;
  double  b;
  double  base;
  double  log_c;
  double  log_m_c;
  double  dg0_scale_factor;
  /*
  double  flux_scaling;
  double  frb;
  double  lrb;
  double  recip_frb;
  double  recip_base_rt;
  */
  double  pt;
  double  rt;
  double  max_log_g0_sum;
  double  sflux;
  int64_t *molecules_ptrs;
  int64_t *rxn_indices;
  int64_t *coefficients;
  int64_t *rxn_ptrs;
  int64_t *molecule_indices;
  int64_t *rcoefficients;
  int num_species;
  int num_rxns;
  int rxn;
  int padi;
  /*
  int base_rxn;
  */
  int i;
  int j;

  int mi;
  int success;

  FILE *lfp;
  FILE *efp;
  /*
#define DBG 1
  */
  /*
    Check that base_rxn is in range.
  */
  success = 1;
  num_rxns = state->number_reactions;
  num_species = state->nunique_molecules;
  molecules   = state->sorted_molecules;
  molecules_matrix = state->molecules_matrix;
  molecules_ptrs   = molecules_matrix->molecules_ptrs;
  rxn_indices      = molecules_matrix->rxn_indices;
  coefficients     = molecules_matrix->coefficients;
  rxn_matrix       = state->reactions_matrix;
  rxn_ptrs         = rxn_matrix->rxn_ptrs;
  molecule_indices = rxn_matrix->molecules_indices;
  rcoefficients    = rxn_matrix->coefficients;
  c                = state->product_term;
  /*
  dg0_scale_factor = .001;
  */
  dg0_scale_factor = state->dg0_scale_factor;
  /*
  product_term     = state->product_term;
  reactant_term    = state->reactant_term;
  base_rxn         = (int)state->base_reaction;
  counts           = state->ode_counts;
  conc_to_count    = state->conc_to_count;
  */
  ke               = state->ke;
  rke              = state->rke;
  log_kf_rel       = state->log_kf_rel;
  log_kr_rel       = state->log_kr_rel;
  max_log_g0_sum   = state->max_log_g0_sum;
  /*
  forward_rxn_likelihoods = state->ode_forward_lklhds;
  reverse_rxn_likelihoods = state->ode_reverse_lklhds;
  */
  if (success) {
    /*
    for (i=0;i<num_species;i++) {
      counts[i] = concs[i] * conc_to_count[i];
    }
    success = update_rxn_likelihoods(state,counts,forward_rxn_likelihoods,
				     reverse_rxn_likelihoods);
    */
  }
  /*
  state->flux_scaling = 0.0;
  */
  /*
    kf(base) * concs[base_reaction_reactants]^nu)
  */
  /*
  flux_scaling     = compute_flux_scaling(state,concs);
  */
  lfp      = state->lfp;
  /*
    The flux for a species is given by the sum of contributions from
    the reactions it is present in.

    The contribution, c, from a reaction with forward rate constant kf;
    reverse rate constant kr; reciprocal equilibrium constant, rke (kr/kf);
    product of reactant concentrations, R; 
    and product of product concentrations, P; for a product species is
    
    c = kf * R - kr * P
    
    but kf and kr may not be computable even though c might be quite reasonable
    ad kf and kr are nearly the same size.
    But we have rke = kr/kf, and we have log(kf) and log(kr) available.

    then we write c = kf * (R - rke * P)

    and then let b = R - rke * P

    then if b > 0 then 
    
     we may write 

           log(c) = log(kf) + log(b);

	   and c = exp(log(c));

     if b = 0, c = 0

     if b < 0 then
         log(-c) = log(kf) + log(-b)
	  c     = -exp(log(-c))

     so c is a product_per_reaction_flux_contribution,
     and -c is a reactant_per_reaction_flux_contribution.

     As we won't need to use the product_term vector used in other routines
     here we use it to store the vector c.
 
     For a reactant species the sign is merely flipped.    

     NB one might want to factor out kr instead of kb above then we 
     have 
              c = kr * (ke * R - P)

	      let b = ke * R - P then

	      log_c = log(kr) + log(b) when b> 0

	      log_m_c = log(kr) + log(-b) when b < 0
    
              so we might want to choose the kr or kf with the smaller 
	      absolute value so that the b term is larger (and as its a
	      subtraction would have less cancellation error).
         
  */
  if (success) {
    for (i=0;i<num_rxns;i++) {
      pt = 1.0;
      rt = 1.0;
      for (j=rxn_ptrs[i];j<rxn_ptrs[i+1];j++) {
	mi = molecule_indices[j];
	if (rcoefficients[j] < 0) {
	  /*
	  rt = rt * counts[mi];
	  */
	  rt = rt * concs[mi];
	} else {
	  if (rcoefficients[j] > 0) {
	    /*
	    pt = pt * counts[mi];
	    */
	    pt = pt * concs[mi];
	  }
	}
      } /* end for j */
      if (abs(log_kf_rel[i]) < abs(log_kr_rel[i])) {
	b = rt - (rke[i] * pt);
	base = log_kf_rel[i] * dg0_scale_factor;
      } else {
	b = (ke[i] * rt) - pt;
	base = log_kr_rel[i] * dg0_scale_factor;
      }
      if (b == 0.0) {
	/*
	  reaction is at equilbrium
	*/
	c[i] = 0.0;
      } else {
	if (b > 0) {
	  log_c = base + log(b);
	  /* 
	    here we should check on suitability of  log_c for 
	    an arguemt to exp (should be in [-704:704])
	  */
	  if (log_c > max_log_g0_sum) {
	    if (lfp) {
	      fprintf(lfp,
		      "lr6_approximate_delta_concs: log_c > %le, set to %le\n",
		      max_log_g0_sum,max_log_g0_sum);
	      fflush(lfp);
	    }
	    /*
	    success = 0;
	    */
	    log_c = max_log_g0_sum;
	  }
	  c[i] = exp(log_c);
	} else {
	  log_m_c = base + log(-b);
	  /* 
	    here we should check on suitability of  log_c for 
	    an arguemt to exp (should be in [-704:704])
	  */
	  if (log_m_c > max_log_g0_sum) {
	    if (lfp) {
	      fprintf(lfp,
	      "lr6_approximate_delta_concs: log_m_c > %le, set to %le\n",
		      max_log_g0_sum,max_log_g0_sum);
	      fflush(lfp);
	    }
	    /*
	    success = 0;
	    */
	    log_m_c = max_log_g0_sum;
	  }
	  c[i] = -exp(log_m_c);
	}
      }
    } /* end for (i...) */
  } /* end if (success) */
  //------------------------------------
  if (success) {
    molecule = molecules;
    for (i=0;i<num_species;i++) {
      sflux = 0.0;
      if (molecule->variable == 1) {
	for (j=molecules_ptrs[i];j<molecules_ptrs[i+1];j++) {
	  rxn = rxn_indices[j];
	  if (coefficients[j] < 0) {
	    sflux -= c[rxn];
	  } else {
	    if (coefficients[j] > 0) {
	      sflux += c[rxn];
	    }
	  }
	} /* end for(j...) */
	/*
	flux[i] = flux_scaling * ((forward * recip_frb) + (backward * lrb));
	*/
	flux[i] = sflux;
      } else {
	flux[i] = 0.0;
      }
      molecule += 1; /* Caution address arithmetic heare. */
    } /* end for (i...) */
  } /* end if success */
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"Mol_index\t   conc   \t    flux\n");
    for (i=0;i<num_species;i++) {
      fprintf(lfp,"%d\t%le\t%le\n",
	      i,concs[i],flux[i]);
    }
    fflush(lfp);
  }
#endif
  return (success);
}
