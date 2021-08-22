import numpy as np

n = 2
c_0 = 1
D = 2**(c_0 * n)
K = 1000 #Make sure this value works, described in Lemma 6

def each_vec(B):
    gcds = []

    for i in range(len(B)):
        values = B[i].tolist()
        spot = 0
        for num in values:
            values[spot] = int(num)
            spot += 1

        gcds.append(np.gcd.reduce(values))

    vec_points = []
    for vecs in range(gcds[0] + 1):
        new_vec = B[0] * (vecs / gcds[0])
        vec_points.append(new_vec)
    for basis in range(1, len(B)):
        new_vecs = []
        for vecs in range(1, gcds[basis] + 1):
            new_vec = B[basis] * (vecs/gcds[basis])
            new_vecs.append(new_vec)
        new_comb_vecs = []
        for v in range(len(vec_points)):
            for n in range(len(new_vecs)):
                new_comb_vecs.append(vec_points[v] + new_vecs[n])
        for vec in new_comb_vecs:
            vec_points.append(vec)

    return vec_points

def add_perturbation(S, K, n):
    perturbed = []
    v_perturbations = []
    for biggoIndex in range(len(S)):
        perturbation = np.random.normal(0, 1/((K*n)**(1/2)))
        new_v = np.array([S[biggoIndex][0] + perturbation])
        v_perturbation = np.array([perturbation])
        for index in range(1, len(S[biggoIndex])):
            more_perturbs = np.random.normal(0, 1/((K*n)**(1/2)))
            new_boi = S[biggoIndex][index] + more_perturbs
            new_v = np.insert(new_v, len(new_v), new_boi)
            v_perturbation = np.insert(v_perturbation, len(v_perturbation), more_perturbs)
        perturbed.append(new_v)
        v_perturbations.append(v_perturbation)
    return perturbed, v_perturbations



def sample_vectors(L, n, D, K):
    orthonormal_basis = np.identity(n)
    fi_basis = orthonormal_basis * D
    fi_in_lattice_basis = np.matmul(fi_basis, np.linalg.inv(L))
    fi_int_L_basis = np.rint(fi_in_lattice_basis)
    each_sample_point = each_vec(fi_int_L_basis)
    perturbed_samples, perturbations = add_perturbation(each_sample_point, K, n)







L = np.array([[1, 0], [1, 1]])
sample_vectors(L, n, D, K)
