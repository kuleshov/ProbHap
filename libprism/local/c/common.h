/****** CONVERSION FUNCTIONS ******/

/* Returns 1 if i is in the s-th subset of the powerset of {0,...,N-1}
   and 0 otherwise. */
int i_in_s(int i, unsigned int s){
    return ((s >> i) % 2 == 1);
}

unsigned long ys_to_m(int y, unsigned int s, unsigned int N_subsets)
{
    unsigned long m = N_subsets * y + s;
    return m;
}

/****** SET MANIPULATION FUNCTIONS ******/

/* Takes a set s of active indices at the current position. Puts into
*active_set the subset of active indices at the current position that also
appear at the previous position, in the coordinates of the previous position.
 Does the same with *inactive_set and inactive indices. */

// This function is currently not used anywhere.
void set_active_and_inactive(unsigned int *active_set,
                             unsigned int *inactive_set,
                             unsigned int s, int N, int* current_to_prev_pos)
{
    int i0 = 0;
    *active_set = 0;
    *inactive_set = 0;

    for (int i=0; i<N; i++)
    {
        i0 = current_to_prev_pos[i];
        if (i0 == -1) continue;

        // if i is in s:
        if ((s >> i) % 2 == 1) *active_set |= 1 << i0; // add to active indices
        else *inactive_set |= 1 << i0;              // add to inactive indices
    }
}

unsigned int set_mask(unsigned int s, int N, int* current_to_prev_pos)
{
    int i0 = 0;
    unsigned int mask = 0;

    for (int i=0; i<N; i++)
    {
        i0 = current_to_prev_pos[i];
        if (i0 == -1) continue;

        // if i is in s, set it in the mask
        if ((s >> i) % 2 == 1) mask |= 1 << i0;
    }

    return mask;
}

/* Returns 1 if set s0 is consistent with set s (i.e. all active indices of
   s are active in s0, and all inactive indices are also inactive).
   Indices of s must be already converted to the numbering of s0. */
int is_consistent(unsigned int s0, unsigned int s_active,
                  unsigned int s_inactive)
{
    return (((s0 & s_active) == s_active) && ((s0 & s_inactive) == 0));
}

unsigned int full_prev_subset(unsigned int s_free, int N_free,
                              int* free_indices, unsigned int mask)
{
    unsigned int s_prev = 0;
    int i = 0;

    // set clamped indices
    s_prev = (s_prev | mask);

    for (int i_free = 0; i_free<N_free; i_free++)
    {
        i = free_indices[i_free];

        // if index is in s_free, set it in s_prev
        if ((s_free >> i_free) % 2 == 1)
            s_prev |= 1 << i;

    }

    return s_prev;
}

/****** FUNCTIONS FOR CALCULATING PROBABILITIES ******/

/* Returns the emission likelihood at a heterozygous site.
        y: allele at strand zero
        Nclouds: number of clouds
        clouds_at_position: array of cloud position
        prob0: for each cloud, probability of emitting 0
        prob1: for each cloud, probability of emitting 1
*/

double logaddexp(double log_p, double log_q)
{
	if (log_p == -INFINITY)
		return log_q;
	else if (log_q == -INFINITY)
		return log_p;
	else if (abs(log_p - log_q) > 30)
		return fmax(log_p, log_q);
	else
		return log_p + log(1 + exp(log_q - log_p));
}

double emission_lik(int y, unsigned int s, int Nclouds,
                    int* clouds_at_position, double* prob0, double* prob1)
{
    double log_prob = 0.0;
    int strand0 = (0 + y) % 2;
    int strand1 = (1 + y) % 2;
    int i = 0;

    for (int c=0; c<Nclouds; c++)
    {
        i = clouds_at_position[c];

        // if i is in the set indexed by s
        if ((s >> i) % 2)
        {
            // segment associates with strand 1
            if (strand1 == 0) log_prob += prob0[i];
            else if (strand1 == 1) log_prob += prob1[i];
        }
        else
        {
            // segment associates with strand 0
            if (strand0 == 0) log_prob += prob0[i];
            else if (strand0 == 1) log_prob += prob1[i];
        }

        if (isnan(log_prob)){
            printf("~ %d %d %d %f %f\n", y, s, i, prob0[i], prob1[i]);
            exit(1);
        }
    }

    return log_prob;
}