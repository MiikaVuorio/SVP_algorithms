import numpy as np
import random
import math

test = np.array([[1/math.sqrt(2), 1/math.sqrt(2)], [1,-1]])
L = np.array([[4, 1, 1, 1], [5, 2, 2, 1], [3, 1, 3, 1], [2,5,6,1]])
n = len(L[0])
c_0 = 2
#c_1 must be large enough that (N-O(n*2^(2*n)))*2^-O(log(n)) by lemmas 10, 11. Maybe add a check for this, three should work for 4d for 5d, 2 might work
c_1 = 2 + math.log(n, 2)/n + 0.000001
c_3 = 2
c_5 = 4 * c_3
D = 2 ** (c_0 * n)
K = 100  # Make sure this value works, described in Lemma 6
min_num_vectors = 2 * 8 ** n #Do not know why I at some point thought this'd be the minimum required amount of vectors, it isn't, the exception arising from this has been removed
N = math.ceil(2**(n*c_1))
print(N)

def unitary_random_vectors(B, L, N, n):
    standard_basis_parallelepiped = np.matmul(B, L)
    uni_rand_vectors = []
    while True:
        rando_vec = np.zeros(n)
        for d in B:
            values = []
            for v in d:
                if v < 0:
                    min = v
                    max = 0
                else:
                    min = 0
                    max = v
                values.append(random.randint(min, max))
            new_vec = np.array(values)
            rando_vec += new_vec

        rand_vec_standard_basis = np.matmul(rando_vec, L)

        # checking whether vector is in parralelepiped
        in_parallelepiped = True
        for v in standard_basis_parallelepiped:
            dot = np.dot(v, rand_vec_standard_basis)
            v_norm = np.linalg.norm(v)
            if dot < 0 or dot/v_norm > v_norm or np.linalg.norm(rand_vec_standard_basis) == 0:
                in_parallelepiped = False
                break

        # ADD removal of duplicates
        if in_parallelepiped:
            uni_rand_vectors.append(rando_vec)
            # print("ACCEPTED")
            # print(rando_vec)
        else:
            # print("REJECTED")
            # print(rando_vec)
            pass
        #print(len(uni_rand_vectors))
        if len(uni_rand_vectors) >= N:
            break

    return uni_rand_vectors

# Unused function, thought of a cool way to uniformly choose random vectors from inside, turns out it don't work
def each_vec(B, dimensions):
    gcds = []
    for i in range(len(B)):
        values = B[i].tolist()
        spot = 0
        for num in values:
            values[spot] = int(abs(num))
            spot += 1
        # print(values)
        gcds.append(np.gcd.reduce(values))
    # print(gcds)
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
    # print(len(vec_points))
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
    for numerator in range(8):
        new_cent = paralellepiped_basis[0] * (numerator / 8) + paralellepiped_basis[0] / 16
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


def assign_vectors(circ_centres, vectors, R, n):
    r = R / 4

    # print(len(circ_centres))

    assigned_centres = []
    vec_index = 0
    centre_dic = {}
    for i in range(len(circ_centres)):
        centre_dic[i] = []
    #print(centre_dic)
    #print(len(circ_centres))
    for v in vectors:
        centre_index = 0
        for centre in circ_centres:
            if np.linalg.norm(v - centre) <= r:
                assigned_centres.append(centre_index)
                centre_dic[centre_index].append(vec_index)
                vec_index += 1
                break
            centre_index += 1
            if centre_index == 8**n:
                print(v)

                raise Exception("This ain't supposed to happen")
    #print(len(centre_dic))
    return assigned_centres, centre_dic

def assign_dict_vectors(circ_centres, vectors, R, n):
    r = R / 4

    # print(len(circ_centres))

    assigned_centres = []
    centre_dic = {}
    for i in range(len(circ_centres)):
        centre_dic[i] = []
    #print(centre_dic)
    #print(len(circ_centres))
    for v in vectors.keys():
        centre_index = 0
        for centre in circ_centres:
            if np.linalg.norm(vectors[v] - centre) <= r:
                assigned_centres.append(centre_index)
                centre_dic[centre_index].append(v)
                break
            centre_index += 1
            if centre_index == 8**n:
                print(8**n)
                print(len(circ_centres))
                print(R)
                print(vectors[v])

                raise Exception("This ain't supposed to happen")
    #print(len(centre_dic))
    return assigned_centres, centre_dic


def change_from_basis(B, S):
    basis_changed = []
    for v in S:
        basis_changed.append(np.matmul(v, B))
    return basis_changed


def choose_reps(centre_dic, min_num_vectors):
    reps = []
    saved_from_execution = []
    marked_for_execution = []
    # print(len(centre_dic.keys()))
    for i in centre_dic.keys():
        if len(centre_dic[i]) != 0:
            index = len(centre_dic[i]) // 2
            reps.append(centre_dic[i][index])
            centre_dic[i].pop(index)
            if len(centre_dic[i]) == 0:
                saved_from_execution.append(i)
    for i in centre_dic.keys():
        if len(centre_dic[i]) == 0 and i not in saved_from_execution:
            marked_for_execution.append(i)
    for identity in marked_for_execution:
        del centre_dic[identity]



    # if len([item for subl in centre_dic.values() for item in subl]) < min_num_vectors:
    #     raise Exception("not enuf vectors")

    #print(len(centre_dic.keys()))
    #print(len(reps))
    return centre_dic, reps


def to_dic(vecs):
    xi_dic = {}
    for i in range(len(vecs)):
        xi_dic[i] = vecs[i]
    return xi_dic


def sieve(current_xis, current_ais, zis, no_reps, reps):
    new_ais = {}
    temp_fix = 0 #better fix is to make reps into a dictionary form as well
    # print(len([item for subl in no_reps.values() for item in subl]))
    # print(len(no_reps.keys()))
    # print(len(reps))
    for i in no_reps.keys():
        for index in range(len(no_reps[i])):
            new_ais[no_reps[i][index]] = current_ais[no_reps[i][index]] + zis[reps[temp_fix]] - current_ais[reps[temp_fix]]

        del current_xis[reps[temp_fix]]
        del zis[reps[temp_fix]]
        temp_fix += 1


    return current_xis, new_ais, zis

def post_sieve_ball_fun(R, n):

    ball_centres_per_dimension = []
    ball_centres = []
    for dimension in range(n):
        ball_centres_per_dimension.append([])
        for multiplier in range(8):
            ball_place = np.zeros(n)
            ball_place[dimension] = (-R*7/16 + R * multiplier/8)*2
            ball_centres_per_dimension[dimension].append(ball_place)

    final_ball_centres = []
    for bb in ball_centres_per_dimension[0]:
        ball_centres.append(bb)
    for bingobing in range(n-1):
        final_ball_centres = []
        for b in ball_centres_per_dimension[bingobing + 1]:
            #print(len(ball_centres))
            for ball in ball_centres:
                king = ball+b
                final_ball_centres.append(king)
                # print(king)
                # print(len(final_ball_centres))
        ball_centres = []
        for biko in final_ball_centres:
            ball_centres.append(biko)
        # print(ball_centres)
        # print(len(ball_centres))
        # print(bingobing)
    #print(len(final_ball_centres))
    #print(final_ball_centres)
    return final_ball_centres


def sample_vectors(L, n, D, K, c_5, min_num_vectors):
    orthonormal_basis = np.identity(n)
    fi_basis = orthonormal_basis * D
    fi_in_lattice_basis = np.matmul(fi_basis, np.linalg.inv(L))
    fi_rint_L_basis = np.rint(fi_in_lattice_basis)
    # each_sample_point = each_vec(fi_rint_L_basis, n)
    each_sample_point = unitary_random_vectors(fi_rint_L_basis, L, N, n)
    zis_in_standard_basis = change_from_basis(L, each_sample_point)
    zis = to_dic(zis_in_standard_basis)
    perturbed_samples, perturbations = add_perturbation(each_sample_point, K, n)
    yis_in_standard_basis = change_from_basis(L, perturbations)
    perturbed_in_standard_basis = change_from_basis(L, perturbed_samples)
    current_xis = to_dic(perturbed_in_standard_basis)
    parallelepiped_in_standard_basis = np.matmul(fi_rint_L_basis, L)
    circ_centres = packing_balls(parallelepiped_in_standard_basis)
    R = np.linalg.norm(parallelepiped_in_standard_basis)
    # print(parallelepiped_in_standard_basis)
    # print(R)
    # long_line = parallelepiped_in_standard_basis[0] + parallelepiped_in_standard_basis[1]
    # print(long_line)
    # print(np.linalg.norm(long_line))
    centre_indexes, centre_dic = assign_vectors(circ_centres, perturbed_in_standard_basis, R, n)
    no_reps, reps = choose_reps(centre_dic, min_num_vectors)

    current_ais = {}
    for i in current_xis.keys():
        current_ais[i] = np.zeros(n)
    # print(current_ais.keys())
    # print(zis.keys())
    current_xis, current_ais, current_zis = sieve(current_xis, current_ais, zis, no_reps, reps)
    # print(current_ais.keys())
    while True:
        survivors = current_zis.keys()
        smalls = 0
        # print(current_zis)
        # print(current_ais)
        # print(current_xis)
        tot_sum = 0
        for key in survivors:
            # print(current_zis[key])
            # print("zis ABOVE")
            # print(current_ais[key])
            diff = np.linalg.norm(current_zis[key] - current_ais[key])
            if diff <= c_5:
                smalls += 1
            tot_sum += diff

        print(tot_sum/len(survivors))


        if smalls == len(survivors):
            break
        else:
            R = R/2
            current_xis_minus_ais = {}
            for key in current_xis.keys():
                current_xis_minus_ais[key] = current_xis[key] - current_ais[key]

            # pp = []
            # for b in range(len(parallelepiped_in_standard_basis)):
            #     basis_vec = np.zeros(n)
            #     basis_vec[b] = R
            #     pp = pp.append(basis_vec)

            circ_centres = post_sieve_ball_fun(R, n)
            centr_int, centre_dictionary = assign_dict_vectors(circ_centres, current_xis_minus_ais, R, n)
            no_reps, reps = choose_reps(centre_dictionary, min_num_vectors)
            current_xis, current_ais, current_zis = sieve(current_xis, current_ais, current_zis, no_reps, reps)

    possible_vectors = []
    #original_vector = []
    for key in current_zis.keys():
        possible_vectors.append(current_zis[key] - current_ais[key])
        #original_vector.append(key)

    sh_vec = np.array([999999])
    sh_len = c_5 + 1
    #it_index = 0
    #og_key = -1
    for vector in possible_vectors:
        vec_len = np.linalg.norm(vector)
        if vec_len < sh_len and vec_len != 0:
            sh_len = vec_len
            sh_vec = vector
            #og_key = original_vector[it_index]
        #it_index += 1
    print(sh_vec)
    print(sh_len)
    #print(each_sample_point[og_key])


sample_vectors(L, n, D, K, c_5, min_num_vectors)
