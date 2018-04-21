#ifndef _BOLTZMANN_WATCH_H_
#define _BOLTZMANN_WATCH_H_ 1
extern void boltzmann_watch(struct state_struct *state,
			    int64_t *choice_view_step_p,
			    int64_t *count_view_step_p,
			    int64_t *lklhd_view_step_p,
			    int64_t *rxn_view_step_p,
			    int64_t *fe_view_step_p,
			    int64_t *rxn_view_pos_p,
			    double  r_sum_likelihood,
			    double  dg_forward,
			    double  entropy,
			    int64_t iter,
			    int64_t rxn_choice);
#endif
