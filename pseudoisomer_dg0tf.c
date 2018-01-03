/* pseudoisomer_dg0tf.c
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
/*
 * compute_molecule_dg0tf.c
 *
 *  Created on: Apr 18, 2013
 *      Author:  Dennis G. Thomas
        Modified by Doug Baxter May 31, 2013
 */
#include "boltzmann_structs.h"
#include "pseudoisomer_dg0tf.h"

double pseudoisomer_dg0tf(double ph, 
			  double mrt, 
			  double ionic_str, 
			  double nh, 
			  double z, 
			  double deltag0){
  /*
    pseudoisomer_dg0tf calcualtes the standard transformed Gibbs energy of 
    formation of a pseudospecies of a compound at a
    given pH and ionic strength, based on Eq. 4.4-10 in Alberty's book.
    Called by: compute_molecule_dg0tf
    Calls:     log,sqrt
   
  */
  double deltag_tf;
  double is_sqrt;
  /*
  deltag_tf = deltag0 + nh*IDEAL_GAS_CONST_KJ_PER_KELVIN_MOL*temp*log(10)*pH-(2.91482*(z*z-nH)*sqrt(ionic_str))/
			(1+1.6*sqrt(ionic_str));
  */
  is_sqrt   = sqrt(ionic_str);
  deltag_tf = deltag0 - (((nh*mrt)*log(10))*ph)-((2.91482*((z*z)-nh)*is_sqrt)/
					   (1.0+1.6*is_sqrt));

  /*
    cout << "[CalcFreeEnergy(speciesDeltaG_tf)] species deltaG0 = " << deltaG0 << " kJ/mol" << endl;
    
    cout << "[CalcFreeEnergy(speciesDeltaG_tf)] species deltaG_tf = " << deltaG_tf << " kJ/mol" << endl;
  */
  return deltag_tf;
}
