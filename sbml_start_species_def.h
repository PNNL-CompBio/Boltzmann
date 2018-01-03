/* sbml_start_species_def.h
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
extern void  sbml_start_species_def(char *species,
				    char *species_name,
				    char *kegg_id,
				    char *comp,
				    double default_recip_volume,
				    double *multiplier,
				    double *recip_volume,
				    double *init_conc,
				    double *init_amt,
				    int *variable,
				    int *in_species_tag,
				    int *in_species_group);
