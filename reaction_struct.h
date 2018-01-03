/* reaction_struct.h
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
#ifndef _REACTION_STRUCT_H_
#define _REACTION_STRUCT_H_ 1
struct reaction_struct {
  int64_t title;
  int64_t pathway;
  int64_t lcompartment;
  int64_t rcompartment;

  double  delta_g0;
  /*
    unit_v is the delta_g0 units, 1.0 for calories 4.184 for Joules.
  */
  double  unit_v;
  double  k_epsilon;
  /*
    activity is to be used as a likelihood multiplier in determining
    which reaction fires. It is in [0,1], 0 means the reaction does not 
    fire, one means its fully activated, things in between represent 
    partial activations, this is for later use.
  */
  double  activity;
  double  enzyme_level;
  /*
    Lines added by DGT on April 17, 2013
    ph is the pH
    temp_kelvin is the temperature in degrees Kelvin
    ionic_strength is the ionic string in molar (M) units.
  */
  double ph;
  double temp_kelvin;
  double ionic_strength;

  int  num_reactants;
  int  num_products;

  int  self_id;
  /*
    unit_i: 0 for Kcalories, 1 for KJoules.
  */
  int  unit_i;
  /*
    Left compartment, and right compartment indicator, -1 = not set.
  */
  int  left_compartment;
  int  right_compartment;

  /*
    Added by DGT on 4/22/2013
    deltag0_computed = 1 if deltag0 of the reaction was successfully computed.
    Otherwise it is set to 0.
  */
  int deltag0_computed;
  int num_regulators;
  
}
;
#endif
