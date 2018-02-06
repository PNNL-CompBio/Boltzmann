/* echo_params.c
*******************************************************************************
boltzmann

Pacific Northwest National Laboratory, Richland, WA 99352.

Copyright (c) 2010 Battelle Memorial Institute.

Publications based on work performed using the software should include 
the following citation as a reference:


Licensed under the Educational Community License, Version 2.0 (the "License"); 
you may not use this file except in compliance with the License. 
The terms and conditions of the License may be found in 
ECL-2.0_LICENSE_TERMS.TXT in the directory containing this file.
        
Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
CONDITIONS OF ANY KIND, either express or implied. See the License for the 
specific language governing permissions and limitations under the License.
******************************************************************************/
#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "echo_params.h"
int echo_params (FILE *lfp, struct state_struct *state) {
  /*
    Echo the state paramaters for the boltzmann code to determine equilbrium
    concentrations of a set of reactions via Monte Carlo methods.
    
    Called by: echo_inputs
    Calls:     fprintf (intrinsic)
  */
  struct cvodes_params_struct *cvodes_params;
  struct ode23tb_params_struct *ode23tb_params;
  int success;
  int pad1;
  cvodes_params = state->cvodes_params;
  ode23tb_params = state->ode23tb_params;
  success = 1;
  if (lfp) {
    fprintf(lfp,"state->params_file    	       = %s\n",state->params_file);
    fprintf(lfp,"state->reaction_file  	       = %s\n",state->reaction_file);
    fprintf(lfp,"state->init_conc_file 	       = %s\n",state->init_conc_file);
    fprintf(lfp,"state->input_dir      	       = %s\n",state->input_dir);
    fprintf(lfp,"state->output_file    	       = %s\n",state->output_file);
    fprintf(lfp,"state->log_file      	       = %s\n",state->log_file);
    fprintf(lfp,"state->counts_out_file        = %s\n",state->counts_out_file);
    fprintf(lfp,"state->concs_out_file         = %s\n",state->concs_out_file);
    fprintf(lfp,"state->rxn_lklhd_file         = %s\n",state->rxn_lklhd_file);
    fprintf(lfp,"state->free_energy_file       = %s\n",state->free_energy_file);
    fprintf(lfp,"state->restart_file           = %s\n",state->restart_file);
    fprintf(lfp,"state->rxn_view_file          = %s\n",state->rxn_view_file);
    fprintf(lfp,"state->bndry_flux_file        = %s\n",state->bndry_flux_file);
    fprintf(lfp,"state->pseudoisomer_file      = %s\n",state->pseudoisomer_file);
    fprintf(lfp,"state->compartment_file       = %s\n",state->compartment_file);
    fprintf(lfp,"state->sbml_file              = %s\n",state->sbml_file);
    fprintf(lfp,"state->ms2js_file             = %s\n",state->ms2js_file);
    fprintf(lfp,"state->kg2js_file             = %s\n",state->kg2js_file);
    fprintf(lfp,"state->rxn_echo_file          = %s\n",state->rxn_echo_file);
    fprintf(lfp,"state->rxn_mat_file           = %s\n",state->rxn_mat_file);
    fprintf(lfp,"state->arxn_mat_file          = %s\n",state->arxn_mat_file);
    fprintf(lfp,"state->dg0ke_file             = %s\n",state->dg0ke_file);
    fprintf(lfp,"state->dictionary_file        = %s\n",state->dictionary_file);
    fprintf(lfp,"state->ode_concs_file         = %s\n",state->ode_concs_file);
    fprintf(lfp,"state->ode_lklhd_file         = %s\n",state->ode_lklhd_file);
    fprintf(lfp,"state->ode_dconcs_file        = %s\n",state->ode_dconcs_file);
    fprintf(lfp,"state->ode_bflux_file         = %s\n",state->ode_bflux_file);
    fprintf(lfp,"state->net_lklhd_file         = %s\n",state->net_lklhd_file);
    fprintf(lfp,"state->ode_sens_file          = %s\n",state->ode_sens_file);
    fprintf(lfp,"state->ode_dsens_file         = %s\n",state->ode_dsens_file);
    fprintf(lfp,"state->dfdke_file             = %s\n",state->dfdke_file);
    fprintf(lfp,"state->dfdmu0_file            = %s\n",state->dfdmu0_file);
    fprintf(lfp,"state->nl_bndry_flx_file      = %s\n",state->nl_bndry_flx_file);
    
    fprintf(lfp,"state->solvent_string         = %s\n",
	    state->solvent_string);
    
    fprintf(lfp,"state->output_dir     	       = %s\n",state->output_dir);
    fprintf(lfp,"state->align_len              = %ld\n",state->align_len);
    fprintf(lfp,"state->max_filename_len       = %ld\n",state->max_filename_len);
    fprintf(lfp,"state->max_param_line_len     = %ld\n",state->max_param_line_len);
    fprintf(lfp,"state->warmup_steps           = %ld\n",state->warmup_steps);
    fprintf(lfp,"state->record_steps           = %ld\n",state->record_steps);
    fprintf(lfp,"state->free_energy_format     = %ld\n",state->free_energy_format);
    fprintf(lfp,"state->rxn_view_freq          = %ld\n",state->rxn_view_freq);
    fprintf(lfp,"state->count_view_freq        = %ld\n",state->count_view_freq);
    fprintf(lfp,"state->lklhd_view_freq        = %ld\n",state->lklhd_view_freq);
    fprintf(lfp,"state->fe_view_freq           = %ld\n",state->fe_view_freq);
    fprintf(lfp,"state->ode_rxn_view_freq      = %ld\n",state->ode_rxn_view_freq);
    fprintf(lfp,"state->print_ode_concs        = %ld\n",state->print_ode_concs);
    fprintf(lfp,"state->adjust_steady_state    = %ld\n",state->adjust_steady_state);
    fprintf(lfp,"state->print_output           = %ld\n",state->print_output);
    fprintf(lfp,"state->print_concs_or_counts  = %ld\n",state->print_concs_or_counts);
    fprintf(lfp,"state->use_bulk_water         = %ld\n",state->use_bulk_water);
    fprintf(lfp,"state->use_activities         = %ld\n",state->use_activities);
    fprintf(lfp,"state->use_deq                = %ld\n",state->use_deq);
    fprintf(lfp,"state->use_pseudoisomers      = %ld\n",state->use_pseudoisomers);
    fprintf(lfp,"state->use_metropolis         = %ld\n",state->use_metropolis);
    fprintf(lfp,"state->use_regulation         = %ld\n",state->use_regulation);
    fprintf(lfp,"state->max_regs_per_rxn       = %ld\n",state->max_regs_per_rxn);
    fprintf(lfp,"state->compute_sensitivities  = %ld\n",state->compute_sensitivities);
    fprintf(lfp,"state->base_reaction          = %ld\n",state->base_reaction);
    fprintf(lfp,"state->ode_solver_choice      = %ld\n",state->ode_solver_choice);
    fprintf(lfp,"state->delta_concs_choice     = %ld\n",state->delta_concs_choice);
    fprintf(lfp,"state->ode_jacobian_choice    = %ld\n",state->ode_jacobian_choice);
    fprintf(lfp,"state->ode_stop_norm          = %ld\n",state->ode_stop_norm);
    fprintf(lfp,"state->ode_stop_rel           = %ld\n",state->ode_stop_rel);
    fprintf(lfp,"state->ode_stop_style         = %ld\n",state->ode_stop_style);
    fprintf(lfp,"state->cvodes_rhs_choice      = %ld\n",state->cvodes_rhs_choice);
    fprintf(lfp,"state->cvodes_jtimes_choice   = %ld\n",state->cvodes_jtimes_choice);
    fprintf(lfp,"state->cvodes_prec_choice     = %ld\n",state->cvodes_prec_choice);
    fprintf(lfp,"state->cvodes_prec_fill       = %ld\n",state->cvodes_prec_fill);
    fprintf(lfp,"cvodes_params->linear_multistep_method = %d\n",
	    cvodes_params->linear_multistep_method);
    fprintf(lfp,"cvodes_params->linear_solver_method    = %d\n",
	    cvodes_params->linear_solver_method);
    fprintf(lfp,"cvodes_params->iterative_method        = %d\n",
	    cvodes_params->iterative_method);
    fprintf(lfp,"cvodes_params->adams_q_max             = %d\n",
	    cvodes_params->adams_q_max);
    fprintf(lfp,"cvodes_params->bdf_q_max               = %d\n",
	    cvodes_params->bdf_q_max);
    fprintf(lfp,"cvodes_params->max_ord                 = %d\n",
	    cvodes_params->max_ord);
    fprintf(lfp,"cvodes_params->mxsteps                 = %d\n",
	    cvodes_params->mxsteps);
    fprintf(lfp,"cvodes_params->mxhnil                  = %d\n",
	    cvodes_params->mxhnil);
    fprintf(lfp,"cvodes_params->use_stab_lim_det        = %d\n",
	    cvodes_params->use_stab_lim_det);
    fprintf(lfp,"cvodes_params->maxnef                  = %d\n",
	    cvodes_params->maxnef);
    fprintf(lfp,"cvodes_params->maxcor                  = %d\n",
	    cvodes_params->maxcor);
    fprintf(lfp,"cvodes_params->maxncf                  = %d\n",
	    cvodes_params->maxncf);
    fprintf(lfp,"cvodes_params->maxl                    = %d\n",
	    cvodes_params->maxl);
    fprintf(lfp,"cvodes_params->pretype                 = %d\n",
	    cvodes_params->pretype);
    fprintf(lfp,"cvodes_params->gstype                  = %d\n",
	    cvodes_params->gstype);
    fprintf(lfp,"cvodes_params->num_cvode_steps         = %d\n",
	    cvodes_params->num_cvode_steps);
    fprintf(lfp,"cvodes_params->reltol                  = %le\n",
	    cvodes_params->reltol);
    fprintf(lfp,"cvodes_params->abstol                  = %le\n",
	    cvodes_params->abstol);
    fprintf(lfp,"cvodes_params->hin                     = %le\n",
	    cvodes_params->hin);
    fprintf(lfp,"cvodes_params->hmin                    = %le\n",
	    cvodes_params->hmin);
    fprintf(lfp,"cvodes_params->hmax                    = %le\n",
	    cvodes_params->hmax);
    fprintf(lfp,"cvodes_params->nlscoef                 = %le\n",
	    cvodes_params->nlscoef);
    fprintf(lfp,"cvodes_params->eplifac                 = %le\n",
	    cvodes_params->eplifac);

    fprintf(lfp,"state->ideal_gas_r            = %le\n",state->ideal_gas_r);
    fprintf(lfp,"state->temp_kelvin            = %le\n",state->temp_kelvin);
    fprintf(lfp,"state->ph                     = %le\n",state->ph);
    fprintf(lfp,"state->ionic_strength         = %le\n",state->ionic_strength);
    fprintf(lfp,"state->default_volume         = %le\n",state->default_volume);
    fprintf(lfp,"state->flux_scaling           = %le\n",state->flux_scaling);
    fprintf(lfp,"state->kf_base_reaction       = %le\n",state->kf_base_reaction);
    fprintf(lfp,"state->ode_t_final            = %le\n",state->ode_t_final);
    fprintf(lfp,"state->min_conc               = %le\n",state->min_conc);
    
    
  }
  return (success);
}
