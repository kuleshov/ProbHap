/* In this script, we mention two types of coordinates: global and local.

    When we refer to a cloud by its position among all clouds,
    we call that global coordinate indexing (e.g., it is the 78th cloud among
    the 1636 clouds we are considering).

    When we refer to a cloud by its position among clouds that cover
    heterozygous site j, we call that local coordinate indexing (e.g.,
    there are 8 clouds that cover site j, and the 78th cloud is the 3rd cloud
    among these 8.)


    We receive from python the following objects (all are either variables or
    numpy arrays):

	bwd: current bwd values at position j
	bwd_next: bwd values at the next j

    np_cloud_values: value of each cloud at site j

    np_prob0: P(0|cloud[j]) for all clouds at position j,
              in local coordinates
    np_prob1: P(1|cloud[j] for all clouds at position j,
              in local coordinates

    np_current_to_next_pos: for each cloud at site j, local index of that cloud
                            at site j-1 (undefined if cloud doesn't appear at
                            j-1)
    np_clouds_at_position: array of cloud indices that actually cover
                           position j

    N: number of segments at current position
    N_bwd_next: number of segments at next position
    Nclouds: number of clouds at that current positions (i.e. clouds that
             actually have reads)

    end_position: 1 if this is the last position, 0 otherwise


    NOTE: This code represents sets of active segments by unsigned ints.
    Therefore it only works when coverage is less than 32.
*/

unsigned int Nsubsets = pow(2, N);
unsigned int Nsubsets_next = pow(2, N_next);
unsigned int Nsubsets_free_next = pow(2, N_free);
double t = 0;

// for speed, convert numpy arrays to C arrays
double* prob0 = (double *) malloc(sizeof(double) * N_next);
double* prob1 = (double *) malloc(sizeof(double) * N_next);
int* current_to_next_pos = (int *) malloc(sizeof(int) * N);
int* clouds_at_position = (int *) malloc(sizeof(int) * Nclouds);
int* free_indices = (int *) malloc(sizeof(int) * N_free);
double* bwd_next_arr = (double *) malloc(sizeof(double) * 2 * Nsubsets_next);

// test by printing the results
for (int i=0; i<N; i++)
    current_to_next_pos[i] = (int) np_current_to_next_pos(i);

for (int i=0; i<N_next; i++){
    prob0[i] = np_prob0(i);
    prob1[i] = np_prob1(i);
}

for (int i=0; i<2*Nsubsets_next; i++)
    bwd_next_arr[i] = bwd_next(i);

for (int i=0; i<Nclouds; i++){
    clouds_at_position[i] = np_clouds_at_position(i);
}

for (int i=0; i<N_free; i++)
    free_indices[i] = np_free_indices(i);

double Nsubsets_new = pow(2.0, N_free);

// start forwards:

for (int y=0; y<2; y++)
{
    for (unsigned int s=0; s<Nsubsets; s++)
    {
        unsigned long m = ys_to_m(y, s, Nsubsets);

        if (end_position)
        {
            bwd(m) = 0.0;
        }
        else
        {
            double curr_M = -INFINITY;
            unsigned int mask = set_mask(s, N, current_to_next_pos);

            for (int y_next=0; y_next<2; y_next++)
            {
                for (unsigned int s_free=0; s_free<Nsubsets_free_next;
                     s_free++)
                {
                    unsigned int s_next = full_prev_subset(s_free, N_free,
                                                           free_indices, mask);
                    unsigned long m_next = ys_to_m(y_next, s_next,
                                                   Nsubsets_next);

                    double E = emission_lik(y_next, s_next, Nclouds,
                                            clouds_at_position, prob0, prob1);

                    t = log(0.5) + // transition probability for y variable
                        log(1.0/Nsubsets_new) + // tr. prob. for z variable
                        E +
                        bwd_next_arr[m_next];

                    curr_M = logaddexp(curr_M, t);
                    if (isnan(curr_M)){
                        printf("> %d %d %d %d %d\n", j, y, s, y_next, s_next);
                        exit(1);

                    }
                }
            }

            bwd(m) = curr_M;
        }
    }
}

free(prob0);
free(prob1);
free(current_to_next_pos);
free(clouds_at_position);
free(free_indices);
free(bwd_next_arr);
