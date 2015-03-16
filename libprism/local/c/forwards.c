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

	fwd: current fwd values at position j
	fwd_prev: fwd values at the previous j

    np_cloud_values: value of each cloud at site j

    np_prob0: P(0|cloud[j]) for all clouds at position j,
              in local coordinates
    np_prob1: P(1|cloud[j] for all clouds at position j,
              in local coordinates

    np_current_to_prev_pos: for each cloud at site j, local index of that cloud
                            at site j-1 (undefined if cloud doesn't appear at
                            j-1)
    np_clouds_at_position: array of cloud indices that actually cover
                           position j
    np_free_indices: indices at previous position that are not bound by the
                         clouds at the current position

    N: number of segments at current position
    N_prev: number of segments at previous position
    Nclouds: number of clouds at that current positions (i.e. clouds that
             actually have reads)

    start_position: 1 if this is the first position, 0 otherwise


    NOTE: This code represents sets of active segments by unsigned ints.
    Therefore it only works when coverage is less than 32.
*/

unsigned int Nsubsets = pow(2, N);
unsigned int Nsubsets_prev = pow(2, N_prev);
unsigned int Nsubsets_free_prev = pow(2, N_free);
double t = 0;

// for speed, convert numpy arrays to C arrays
double* prob0 = (double *) malloc(sizeof(double) * N);
double* prob1 = (double *) malloc(sizeof(double) * N);
int* current_to_prev_pos = (int *) malloc(sizeof(int) * N);
int* clouds_at_position = (int *) malloc(sizeof(int) * Nclouds);
int* free_indices = (int *) malloc(sizeof(int) * N_free);
double* fwd_prev_arr = (double *) malloc(sizeof(double) * 2 * Nsubsets_prev);

// test by printing the results
for (int i=0; i<N; i++)
{
    prob0[i] = np_prob0(i);
    prob1[i] = np_prob1(i);
    current_to_prev_pos[i] = (int) np_current_to_prev_pos(i);
}

for (int i=0; i<2*Nsubsets_prev; i++)
    fwd_prev_arr[i] = fwd_prev(i);

for (int i=0; i<Nclouds; i++){
    clouds_at_position[i] = np_clouds_at_position(i);
}

for (int i=0; i<N_free; i++)
    free_indices[i] = np_free_indices(i);

// count number of new indices at current position:
int N_new = 0;
for (int i=0; i<N; i++)
{
    if (current_to_prev_pos[i] == -1) N_new++;
}
double Nsubsets_new = pow(2, N_new);

// start forwards:

for (int y=0; y<2; y++)
{
    for (unsigned int s=0; s<Nsubsets; s++)
    {
        unsigned long m = ys_to_m(y, s, Nsubsets);

        double E = emission_lik(y, s, Nclouds, clouds_at_position, prob0,
                                prob1);

        if (start_position)
        {
            if (y==0)
                fwd(m) = log(1.0 / ((float) Nsubsets)) + E;
            else
                fwd(m) = -INFINITY;
        }
        else
        {
            double curr_M = -INFINITY;

            unsigned int mask = set_mask(s, N, current_to_prev_pos);

            for (int y_prev=0; y_prev<2; y_prev++)
            {
                for (unsigned int s_free=0; s_free<Nsubsets_free_prev;
                     s_free++)
                {
                    unsigned int s_prev = full_prev_subset(s_free, N_free,
                                                           free_indices, mask);
                    unsigned long m_prev = ys_to_m(y_prev, s_prev,
                                                   Nsubsets_prev);


                    t = log(0.5) + // transition probability for y variable
                        log(1.0/Nsubsets_new) +
                        E +
                        fwd_prev_arr[m_prev];

                    curr_M = logaddexp(curr_M, t);
                }
            }

            fwd(m) = curr_M;
        }
    }
}



free(prob0);
free(prob1);
free(current_to_prev_pos);
free(clouds_at_position);
free(free_indices);
free(fwd_prev_arr);
