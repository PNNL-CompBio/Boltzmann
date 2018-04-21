/* merge_compartments.c
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

#include "merge_compartments.h"
int merge_compartments(struct compartment_struct *list1,
		       struct compartment_struct *list2,
		       struct compartment_struct *mlist,
		       char *compartment_text,
		       int l1,
		       int l2) {
  /*
    Merge two sorted arrays of compartments of length l1 and l2 respectively.
    Called by: sort_compartments, sort_global_compartments
    Calls:     strcmp (intrinsic)
  */
  struct compartment_struct ces;
  struct compartment_struct *p1;
  struct compartment_struct *p2;
  struct compartment_struct *p3;
  size_t move_size;
  size_t e_size;
  char *string1;
  char *string2;
  int success;
  int n;

  int c;
  int j1;

  int j2;
  int j3;

  success = 1;
  j1  = 0;
  j2  = 0;
  j3  = 0;
  n   = l1 + l2;
  p1  = list1;
  p2  = list2;
  p3  = mlist;
  string1 = NULL;
  string2 = NULL;
  e_size = (size_t)sizeof(ces);
  if (p1->string >= 0) {
    string1 = (char*)&compartment_text[p1->string];
  }
  if (p2->string >= 0) {
    string2 = (char*)&compartment_text[p2->string];
  }
  for (j3 = 0;j3 < n; j3++) {
    /*
      Compare first element in list1 to first element in list2.
    */
    if (string1) {
      if (string2) {
	c = strcmp(string1,string2);
      } else {
	c = 1;
      }
    } else {
      c = -1;
    }
    if (c == 0) {
      c = -1;
      if (p1->c_index > p2->c_index) {
	c = 1;
      }
    }
    if (c < 0) {
      /*
	Smaller value was in p1.
      */
      memcpy((void *)p3,(void*)p1,e_size);
      /*
      p3->volume   = p1->volume;
      p3->string = p1->string;
      p3->c_index  = p1->c_index;
      p3->g_index  = p1->g_index;
      */
      j1++;
      p1 += 1; /* Caution Address arithmetic here. */
      p3 += 1; /* Caution Address arithmetic here. */
      /* 
	If we have seen all of list 1, catenate the rest of list 2.
      */
      if (j1 == l1) {
	move_size = (size_t)(l2 - j2) * e_size;
	if (move_size > 0) {
	  memcpy((void *)p3,(void *)p2,move_size);
	}
	/*
	for (j = j2;j<l2;j++) {
	  p3->volume   = p2->volume;
	  p3->string   = p2->string;
	  p3->c_index  = p2->c_index;
	  p3->g_index  = p2->g_index;
	  p2 += 1; // Caution Address arithmetic here. 
	  p3 += 1; // Caution Address arithmetic here. 
	}
        */
	break;
      } else {
	if (p1->string >= 0) {
	  string1 = (char*)&compartment_text[p1->string];
	} else {
	  string1 = NULL;
	}
      }
    } else {
      /*
	Smaller value was in p2.
      */
      memcpy((void *)p3,(void*)p2,e_size);
      /*
      p3->volume   = p2->volume;
      p3->string   = p2->string;
      p3->c_index  = p2->c_index;
      p3->g_index  = p2->g_index;
      */
      j2++;
      p2 += 1; /* Caution Address arithmetic here. */
      p3 += 1; /* Caution Address arithmetic here. */
      /* 
	If we have seen all of list 2, catenate the rest of list 1.
      */
      if (j2 == l2) {
	move_size = (size_t)(l1 - j1) * e_size;
	if (move_size > 0) {
	  memcpy((void *)p3,(void *)p1,move_size);
	}
	/*
	for (j = j1;j<l1;j++) {
	  p3->volume   = p1->volume;
	  p3->string   = p1->string;
	  p3->c_index  = p1->c_index;
	  p3->g_index  = p1->g_index;
	  p1 += 1; // Caution Address arithmetic here. 
	  p3 += 1; // Caution Address arithmetic here. 
	}
	*/
	break;
      } else {
	if (p2->string >= 0) {
	  string2 = (char*)&compartment_text[p2->string];
	} else {
	  string2 = NULL;
	}
      }
    }
  }
  return(success);
}
