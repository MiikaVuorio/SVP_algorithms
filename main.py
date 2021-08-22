import numpy as np
import random
import math

n = 2
c_0 = 12
c_3 = 2
c_5 = 4 * c_3
D = 2 ** (c_0 * n)
print(D)
K = 1000  # Make sure this value works, described in Lemma 6
start_vector_count = -1  #check what this number should be from the paper
min_num_vectors = 2 * 8 ** n


def unitary_random_vectors(B):
    while True:
        rando_vec = np.zeros(2)
        for d in B:
            new_vec = np.array([])
            for v in d:
                np.insert(new_vec, random.randint(0, v))
            rando_vec += new_vec

def each_vec(B, dimensions):
    gcds = []
    print(B)
    for i in range(len(B)):
        values = B[i].tolist()
        spot = 0
        for num in values:
            values[spot] = int(abs(num))
            spot += 1
        print(values)
        gcds.append(np.gcd.reduce(values))
    print(gcds)
    vec_points = []
    for vecs in range(gcds[0] + 1):
        new_vec = B[0] * (vecs / gcds[0])
        vec_points.append(new_vec)
    for basis in range(1, len(B)):
        new_vecs = []
        for vecs in range(1, gcds[basis] + 1):
            new_vec = B[basis] * (vecs / gcds[basis])
            new_vecs.append(new_vec)
        new_comb_vecs = []
        for v in range(len(vec_points)):
            for n in range(len(new_vecs)):
                new_comb_vecs.append(vec_points[v] + new_vecs[n])
        for vec in new_comb_vecs:
            vec_points.append(vec)

    min_num = 2 * 8 ** dimensions
    print(len(vec_points))
    if len(vec_points) < min_num:
        raise Exception("Not enuf vectars")

    return vec_points


def add_perturbation(S, K, n):
    perturbed = []
    v_perturbations = []
    for biggoIndex in range(len(S)):
        perturbation = np.random.normal(0, 1 / ((K * n) ** (1 / 2)))
        new_v = np.array([S[biggoIndex][0] + perturbation])
        v_perturbation = np.array([perturbation])
        for index in range(1, len(S[biggoIndex])):
            more_perturbs = np.random.normal(0, 1 / ((K * n) ** (1 / 2)))
            new_boi = S[biggoIndex][index] + more_perturbs
            new_v = np.insert(new_v, len(new_v), new_boi)
            v_perturbation = np.insert(v_perturbation, len(v_perturbation), more_perturbs)
        perturbed.append(new_v)
        v_perturbations.append(v_perturbation)
    return perturbed, v_perturbations


def packing_balls(paralellepiped_basis):
    circ_centres = []
    for divisor in range(8):
        new_cent = paralellepiped_basis[0] * (divisor / 8) + paralellepiped_basis[0] / 16
        circ_centres.append(new_cent)
    for basis in range(1, len(paralellepiped_basis)):
        new_vecs = []
        for vecs in range(1, 8):
            new_vec = paralellepiped_basis[basis] * (vecs / 8) + paralellepiped_basis[basis] / 16
            new_vecs.append(new_vec)
        new_comb_vecs = []
        for v in range(len(circ_centres)):
            for n in range(len(new_vecs)):
                new_comb_vecs.append(circ_centres[v] + new_vecs[n])
        for vec in new_comb_vecs:
            circ_centres.append(vec)

    return circ_centres


# def euclidean_length(M):
#     sum = 0
#     for vector in M:
#

def assign_vectors(circ_centres, vectors, R):
    r = R / 4
    assigned_centres = []
    vec_index = 0
    centre_dic = {}
    for i in range(len(circ_centres)):
        centre_dic[i] = []
    for v in vectors:
        centre_index = 0
        for centre in circ_centres:
            if np.linalg.norm(v - centre) <= r:
                assigned_centres.append(centre_index)
                centre_dic[centre_index].append(vec_index)
                vec_index += 1
                break
            centre_index += 1

    return assigned_centres, centre_dic


def change_from_basis(B, S):
    basis_changed = []
    for v in S:
        basis_changed.append(np.matmul(v, B))
    return basis_changed


def choose_reps(centre_dic, min_num_vectors):
    reps = []
    for i in range(len(centre_dic.keys())):
        if len(centre_dic[i]) == 0:
            del centre_dic[i]
            continue
        index = len(centre_dic[i]) // 2
        reps.append(centre_dic[i][index])
        centre_dic[i].pop(index)
    if len([item for subl in centre_dic.values() for item in subl]) < min_num_vectors:
        raise Exception("not enuf vectors")
    return centre_dic, reps


def to_dic(vecs):
    xi_dic = {}
    for i in range(len(vecs)):
        xi_dic[i] = vecs[i]
    return xi_dic


def sieve(current_xis, current_ais, zis, no_reps, reps):
    new_ais = {}
    for i in range(len(reps)):
        for index in range(len(no_reps[i])):
            new_ais[no_reps[i][index]] = []
            new_ais[no_reps[i][index]].append(current_ais[no_reps[i][index]] + zis[reps[i]] - current_ais[reps[i]])

        del current_xis[reps[i]]
        del zis[reps[i]]
    return current_xis, new_ais, zis


def sample_vectors(L, n, D, K, c_5, min_num_vectors):
    orthonormal_basis = np.identity(n)
    fi_basis = orthonormal_basis * D
    fi_in_lattice_basis = np.matmul(fi_basis, np.linalg.inv(L))
    fi_int_L_basis = np.rint(fi_in_lattice_basis)
    each_sample_point = each_vec(fi_int_L_basis, n)
    zis_in_standard_basis = change_from_basis(L, each_sample_point)
    zis = to_dic(zis_in_standard_basis)
    perturbed_samples, perturbations = add_perturbation(each_sample_point, K, n)
    yis_in_standard_basis = change_from_basis(L, perturbations)
    perturbed_in_standard_basis = change_from_basis(L, perturbed_samples)
    current_xis = to_dic(perturbed_in_standard_basis)
    parallelepiped_in_standard_basis = np.matmul(fi_int_L_basis, L)
    circ_centres = packing_balls(parallelepiped_in_standard_basis)
    R = np.linalg.norm(parallelepiped_in_standard_basis)
    centre_indexes, centre_dic = assign_vectors(circ_centres, perturbed_in_standard_basis, R)
    no_reps, reps = choose_reps(centre_dic, min_num_vectors)

    current_ais = {}
    for i in current_xis.keys():
        current_ais[i] = np.zeros(2)

    current_xis, current_ais, current_zis = sieve(current_xis, current_ais, zis, no_reps, reps)

    while True:
        survivors = current_zis.keys()
        smalls = 0
        for key in survivors:
            if np.linalg.norm(current_zis[key] - current_ais[key]) <= c_5:
                smalls += 1
        if smalls == len(survivors):
            break
        else:
            no_reps, reps = choose_reps(no_reps, min_num_vectors)
            current_xis, current_ais, current_zis = sieve(current_xis, current_ais, current_zis, no_reps, reps)

    possible_vectors = []
    for key in current_zis.keys():
        possible_vectors.append(current_zis[key]-current_ais[key])

    sh_vec = np.array([999999])
    sh_len = c_5 + 1
    it_index = 0
    for vector in possible_vectors:
        vec_len = np.linalg.norm(vector)
        if vec_len < sh_len:
            sh_len = vec_len
            sh_vec = vector
        it_index += 1
    print(sh_vec)







L = np.array([[4, 1], [5, 2]])
sample_vectors(L, n, D, K, c_5, min_num_vectors)
