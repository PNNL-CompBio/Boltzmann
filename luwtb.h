/* luwtb.h
*******************************************************************************
MSPolygraph

Pacific Northwest National Laboratory, Richland, WA 99352.

MSPolygraph is a mass spectrogram analysis tool that identifies likely 
matching candidate peptides from a reference database.

Copyright (c) 2010 Battelle Memorial Institute.

Publications based on work performed using the software should include 
the following citation as a reference:

    Cannon WR, KH Jarman, BM Webb-Robertson, DJ Baxter, CS Oehmen, 
    KD Jarman, A Heredia-Langner, GA Anderson, and KJ Auberry.
    A Comparison of Probability and Likelihood Models for Peptide 
    Identification from Tandem Mass Spectrometry Data, 
    J. Proteome Res., 2005 4(5):1687-1698

Licensed under the Educational Community License, Version 2.0 (the "License"); 
you may not use this file except in compliance with the License. 
The terms and conditions of the License may be found in 
ECL-2.0_LICENSE_TERMS.TXT in the directory containing this file.
        
Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
CONDITIONS OF ANY KIND, either express or implied. See the License for the 
specific language governing permissions and limitations under the License.
******************************************************************************/
void show_stack(void)
{

    int err;
    char foo[1024];
    unw_context_t ucp;

    unw_cursor_t cursor;
    unw_word_t ip, sp;
    unw_word_t offp;
    unw_proc_info_t pip;

    err = unw_getcontext(&ucp);

    if(err) {
        fprintf(stderr, "unw_getcontext(&ucp) failed\n");
    }
    
    unw_init_local(&cursor, &ucp);

    /* FIXME: Function call error handling */
    do {
        /* Get procedure name */
        err = unw_get_proc_name(&cursor, foo, sizeof(foo), &offp);
        /* Procedure info */
        err = unw_get_proc_info(&cursor, &pip);
        /* Instruction Ptr */
        err = unw_get_reg(&cursor, UNW_REG_IP, &ip);
        /* Stact Ptr */
        err = unw_get_reg(&cursor, UNW_REG_SP, &sp);

        fprintf (stderr,"Func = %s\toffset = 0x%lx\tip = 0x%lx\tsp = 0x%lx\n",
               foo, (long) offp, (long) ip, (long) sp);
    } while (unw_step(&cursor) > 0);
    
}

void handler(int signal, siginfo_t *info, void *unused)
{
    fprintf(stderr, "Handler Received signal %d, errno %d code = %d\n",
        signal, info->si_errno, info->si_code);

    if(signal == SIGSEGV && info->si_code == SEGV_MAPERR)
        fprintf(stderr, "SEGV: address not mapped to object\n");
    if(signal == SIGSEGV && info->si_code == SEGV_ACCERR)
        fprintf(stderr, "SEGV: invalid permissions for mapped object\n");

    show_stack();
    
    exit(0);
}

