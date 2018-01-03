/* compute_molecule_dg0tf.c
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

#include "compute_molecule_dg0tf.h"
int compute_molecule_dg0tf(double ph, 
			   double m_rt, 
			   double m_r_rt,
			   double ionic_str,
			   char *cpd_name, 
			   struct pseudoisomer_struct *pseudoisomers,
			   char *pseudoisomer_strings,
			   int num_cpds,
			   double *deltag_tf) {
  /*
    compute the standard transformed Gibbs energy of formation of the
    pseudoisomer group of the molecule. 
    Equation 4.5-1
    Called by: compute_molecule_dg0tfs.
    Calls      pseudoisomer_dg0tf, log, exp

    Dougs comment: This process takes time proportional to the product
    of the number of pseudoisomers and the number of unique molecules.
    By adding a sort of the pseudoisomers we can make the whole process
    take time proportional to the sum of the number of pseudoisomers and
    the number of unique molecules. Will work on that later.
  */
  struct pseudoisomer_struct *pseudoisomer;
  char *json_cpd_name;

  int success;
  int cmp1;

  int i;
  int index;


  int priority_level;
  int count1,count2;
  
  int check;
  int padi;

  double deltag0;
  double nh, z;
  double min;
  double sum;
  
  double dgtc;
  double dgt1;
  double dgtmin;
  
  /*
  double deltag_temp[50];
  */

  success = 1;

  priority_level = 0;
  check = 0;
  *deltag_tf = 0.0;
  count1 = 0;
  count2 = 0;
  pseudoisomer = pseudoisomers;
  sum = 0;
  for(i=0;i<num_cpds;i++){

    json_cpd_name = (char*)&pseudoisomer_strings[pseudoisomer->json_cpd_name];
    cmp1 = strcmp(cpd_name,json_cpd_name); // compare reaction cpd names and json cpd names

    if(cmp1==0){
      nh = pseudoisomer->nh;
      z = pseudoisomer->z;
      deltag0 = pseudoisomer->dg0_f;
      if((pseudoisomer->priority_level==1) && (priority_level !=2)){
	priority_level = 1;
	check = 1;
	/*
	deltaG_temp[count1] = pseudoisomer_dg0tf(ph,m_rt,
						 ionic_str,nh,z,deltag0);
	*/
	dgtc = pseudoisomer_dg0tf(ph,m_rt,ionic_str,nh,z,deltag0);
	if (count1 == 0) {
	  dgt1 = dgtc;
	  dgtmin = dgtc;
	  sum    = 1.0;
	} else {
	  /*
	  sum = sum + exp((deltag_temp[i]-deltag_temp[index])*m_r_rt);
	  */
	  sum += exp((dgtc - dgt1)*m_r_rt);
	  if (dgtc < dgtmin) {
	    dgtmin = dgtc;
	  }
	}
	count1 = count1+1;
      }	else if((pseudoisomer->priority_level==2) && (priority_level!=1)){
	check = 1;
	/*
	deltag_temp[count2] = pseudoisomer_dg0tf(ph,m_rt,
						 ionic_str,nh,z,deltag0);
	*/
	dgtc = pseudoisomer_dg0tf(ph,m_rt,ionic_str,nh,z,deltag0);
	if (count2 == 0) {
	  dgt1 = dgtc;
	  dgtmin = dgtc;
	  sum = 1.0;
	} else {
	  sum += exp((dgtc - dgt1)*m_r_rt);
	  if (dgtc < dgtmin) {
	    dgtmin = dgtc;
	  }
	}
	count2 = count2+1;
	priority_level = 2;
      }
    }
    pseudoisomer += 1; /* Caution address arithmetic */
  }
  if (dgt1 != dgtmin) {
    sum *= exp((dgt1 - dgtmin)*m_r_rt);
  }
  if (check) {
    /*
    *deltaG_tf = deltaG_temp[index]+m_rt*log(sum);
    */
    *deltag_tf = dgtmin + (m_rt * log(sum));
  } else {
    *deltag_tf = 0.0;
  }
  /*
  min = 0;
  if(count1!=0){
    num = count1;
  }
  else if(count2!=0){
    num = count2;
  }
  */
  /* 
    find the index of the deltaG_temp array element that has the minimum value
  */
  /*
  for(i=0;i<num;i++){
    if(deltag_temp[i]<=min){
      min = deltaG_temp[i];
      index = i;
    }
  }
  sum = 0;
  for(i=0;i<num;i++){
  */
    /*
    sum = sum + exp((deltag_temp[index]-deltag_temp[i])/(IDEAL_GAS_CONST_KJ_PER_KELVIN_MOL*temp));
    */
  /*
    sum = sum + exp((deltag_temp[i]-deltag_temp[index])*m_r_rt);

  }
  */
  /*
  *deltaG_tf = deltaG_temp[index]-IDEAL_GAS_CONST_KJ_PER_KELVIN_MOL*temp*log(sum);
  *deltaG_tf = deltaG_temp[index]+m_rt*log(sum);

  if (*check==0) {
    *deltaG_tf = 0.0;
  }
  */
  return (success);
}
