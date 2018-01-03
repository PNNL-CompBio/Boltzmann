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
#include "echo_params.h"
int echo_params (FILE *lfp, struct state_struct *state) {
  /*
    Echo the state paramaters for the boltzmann code to determine equilbrium
    concentrations of a set of reactions via Monte Carlo methods.
    
    Called by: echo_inputs
    Calls:     fprintf (intrinsic)
 */
  int success;
  int pad1;
  success = 1;
  if (lfp) {
    fprintf(lfp,"state->params_file    	       = %s\n",state->params_file);
    fprintf(lfp,"state->reaction_file  	       = %s\n",state->reaction_file);
    fprintf(lfp,"state->init_conc_file 	       = %s\n",state->init_conc_file);
    fprintf(lfp,"state->input_dir      	       = %s\n",state->input_dir);
    fprintf(lfp,"state->output_file    	       = %s\n",state->output_file);
    fprintf(lfp,"state->log_file      	       = %s\n",state->log_file);
    fprintf(lfp,"state->counts_out_file        = %s\n",state->counts_out_file);
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
    fprintf(lfp,"state->dg0ke_file             = %s\n",state->dg0ke_file);
    fprintf(lfp,"state->dictionary_file        = %s\n",state->dictionary_file);
    fprintf(lfp,"state->ode_concs_file         = %s\n",state->ode_concs_file);
    fprintf(lfp,"state->ode_lklhd_file         = %s\n",state->ode_lklhd_file);
    fprintf(lfp,"state->ode_flux_file          = %s\n",state->ode_flux_file);
    fprintf(lfp,"state->ode_bflux_file         = %s\n",state->ode_bflux_file);
    fprintf(lfp,"state->net_lklhd_file         = %s\n",state->net_lklhd_file);
    fprintf(lfp,"state->nl_bndry_flx_file      = %s\n",state->nl_bndry_flx_file);
    
    fprintf(lfp,"state->solvent_string         = %s\n",
	    state->solvent_string);
    
    fprintf(lfp,"state->output_dir     	       = %s\n",state->output_dir);
    fprintf(lfp,"state->align_len              = %lld\n",state->align_len);
    fprintf(lfp,"state->max_filename_len       = %lld\n",state->max_filename_len);
    fprintf(lfp,"state->max_param_line_len     = %lld\n",state->max_param_line_len);
    fprintf(lfp,"state->warmup_steps           = %lld\n",state->warmup_steps);
    fprintf(lfp,"state->record_steps           = %lld\n",state->record_steps);
    fprintf(lfp,"state->free_energy_format     = %lld\n",state->free_energy_format);
    fprintf(lfp,"state->rxn_view_freq          = %lld\n",state->rxn_view_freq);
    fprintf(lfp,"state->count_view_freq        = %lld\n",state->count_view_freq);
    fprintf(lfp,"state->lklhd_view_freq        = %lld\n",state->lklhd_view_freq);
    fprintf(lfp,"state->fe_view_freq           = %lld\n",state->fe_view_freq);
    fprintf(lfp,"state->ode_rxn_view_freq      = %lld\n",state->ode_rxn_view_freq);
    fprintf(lfp,"state->adjust_steady_state    = %lld\n",state->adjust_steady_state);
    fprintf(lfp,"state->print_output           = %lld\n",state->print_output);
    fprintf(lfp,"state->use_activities         = %lld\n",state->use_activities);
    fprintf(lfp,"state->use_deq                = %lld\n",state->use_deq);
    fprintf(lfp,"state->use_pseudoisomers      = %lld\n",state->use_pseudoisomers);
    fprintf(lfp,"state->use_metropolis         = %lld\n",state->use_metropolis);
    fprintf(lfp,"state->use_regulation         = %lld\n",state->use_regulation);
    fprintf(lfp,"state->max_regs_per_rxn       = %lld\n",state->max_regs_per_rxn);
    fprintf(lfp,"state->base_reaction          = %lld\n",state->base_reaction);
    fprintf(lfp,"state->ode_solver_choice      = %lld\n",state->ode_solver_choice);
    fprintf(lfp,"state->delta_concs_choice      = %lld\n",state->delta_concs_choice);
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
