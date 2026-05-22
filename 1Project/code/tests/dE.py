import numpy as np

N = 10
J = 1.0
Hz = 0.5
D = 1.0

def random_spin():
    v = np.random.normal(size=3)
    return v / np.linalg.norm(v)

def dot(a,b):
    return np.dot(a,b)

def cross_x(a,b):
    return a[1]*b[2] - a[2]*b[1]

def cross_y(a,b):
    return a[2]*b[0] - a[0]*b[2]

# full lattice
def init():
    return np.array([[random_spin() for _ in range(N)] for _ in range(N)])

# FULL ENERGY (truth)
def full_energy(S):
    E = 0.0
    for i in range(N):
        for j in range(N):
            ip = (i+1)%N
            jp = (j+1)%N

            E += -Hz * S[i,j,2]
            E += -J * dot(S[i,j], S[ip,j])
            E += -J * dot(S[i,j], S[i,jp])
            E += -D * cross_x(S[i,j], S[ip,j])
            E += -D * cross_y(S[i,j], S[i,jp])
    return E

# PATCH energy (your structure)
def patch_energy(S, n1, n2):
    E = 0.0
    NC = 1

    for di in range(-NC, NC+1):
        for dj in range(-NC, NC+1):
            i = (n1 + di) % N
            j = (n2 + dj) % N

            ip = (i+1)%N
            im = (i-1)%N
            jp = (j+1)%N
            jm = (j-1)%N

            E += -Hz * S[i,j,2]

            E += -0.5*J*(dot(S[i,j],S[ip,j]) + dot(S[i,j],S[im,j]))
            E += -0.5*J*(dot(S[i,j],S[i,jp]) + dot(S[i,j],S[i,jm]))

            E += -0.5*D*(cross_x(S[i,j],S[ip,j]) + cross_x(S[i,j],S[im,j]))
            E += -0.5*D*(cross_y(S[i,j],S[i,jp]) + cross_y(S[i,j],S[i,jm]))

    return E

# correct local ΔE
def local_delta(S, n1, n2, new_spin):
    S_old = S.copy()
    S_new = S.copy()

    S_new[n1,n2] = new_spin

    sites = [
        (n1,n2),
        ((n1+1)%N,n2), ((n1-1)%N,n2),
        (n1,(n2+1)%N), (n1,(n2-1)%N)
    ]

    def site_energy(T,i,j):
        ip = (i+1)%N
        im = (i-1)%N
        jp = (j+1)%N
        jm = (j-1)%N

        E = -Hz * T[i,j,2]
        E += -J*(dot(T[i,j],T[ip,j]) + dot(T[i,j],T[im,j]) +
                 dot(T[i,j],T[i,jp]) + dot(T[i,j],T[i,jm]))
        E += -D*(cross_x(T[i,j],T[ip,j]) + cross_x(T[i,j],T[im,j]) +
                 cross_y(T[i,j],T[i,jp]) + cross_y(T[i,j],T[i,jm]))
        return E

    dE = 0.0
    for (i,j) in sites:
        dE += site_energy(S_new,i,j) - site_energy(S_old,i,j)

    return dE


# run experiment
S = init()

err_patch = []
err_local = []

for _ in range(10000):

    n1 = np.random.randint(N)
    n2 = np.random.randint(N)

    old_spin = S[n1,n2].copy()
    new_spin = random_spin()

    # GLOBAL ΔE (truth)
    E_before = full_energy(S)
    S[n1,n2] = new_spin
    E_after = full_energy(S)
    dE_global = E_after - E_before

    # PATCH ΔE
    S[n1,n2] = old_spin
    E_before_p = patch_energy(S, n1, n2)
    S[n1,n2] = new_spin
    E_after_p = patch_energy(S, n1, n2)
    dE_patch = E_after_p - E_before_p

    # LOCAL ΔE
    S[n1,n2] = old_spin
    dE_loc = local_delta(S, n1, n2, new_spin)

    err_patch.append(dE_patch - dE_global)
    err_local.append(dE_loc - dE_global)

print("PATCH mean error:", np.mean(err_patch))
print("PATCH std error:", np.std(err_patch))

print("LOCAL mean error:", np.mean(err_local))
print("LOCAL std error:", np.std(err_local))

#%%

import numpy as np

N = 10
J = 1.0
Hz = 0.5
D = 1.0

def random_spin():
    v = np.random.normal(size=3)
    return v / np.linalg.norm(v)

def dot(a,b):
    return np.dot(a,b)

def cross_x(a,b):
    return a[1]*b[2] - a[2]*b[1]

def cross_y(a,b):
    return a[2]*b[0] - a[0]*b[2]

def init():
    return np.array([[random_spin() for _ in range(N)] for _ in range(N)])

# FULL energy (truth)
def full_energy(S):
    E = 0.0
    for i in range(N):
        for j in range(N):
            ip = (i+1)%N
            jp = (j+1)%N

            E += -Hz * S[i,j,2]
            E += -J * dot(S[i,j], S[ip,j])
            E += -J * dot(S[i,j], S[i,jp])
            E += -D * cross_x(S[i,j], S[ip,j])
            E += -D * cross_y(S[i,j], S[i,jp])
    return E

# PATCH energy with variable NC
def patch_energy(S, n1, n2, NC):
    E = 0.0

    for di in range(-NC, NC+1):
        for dj in range(-NC, NC+1):

            i = (n1 + di) % N
            j = (n2 + dj) % N

            ip = (i+1)%N
            im = (i-1)%N
            jp = (j+1)%N
            jm = (j-1)%N

            E += -Hz * S[i,j,2]

            E += -0.5*J*(dot(S[i,j],S[ip,j]) + dot(S[i,j],S[im,j]))
            E += -0.5*J*(dot(S[i,j],S[i,jp]) + dot(S[i,j],S[i,jm]))

            E += -0.5*D*(cross_x(S[i,j],S[ip,j]) + cross_x(S[i,j],S[im,j]))
            E += -0.5*D*(cross_y(S[i,j],S[i,jp]) + cross_y(S[i,j],S[i,jm]))

    return E

# experiment
def run_test(NC_list=[0,1,2,3]):
    errors = {NC: [] for NC in NC_list}

    for _ in range(2000):

        S = init()

        n1 = np.random.randint(N)
        n2 = np.random.randint(N)

        old_spin = S[n1,n2].copy()
        new_spin = random_spin()

        # GLOBAL ΔE (truth)
        E_before = full_energy(S)
        S[n1,n2] = new_spin
        E_after = full_energy(S)
        dE_global = E_after - E_before

        # reset
        S[n1,n2] = old_spin

        for NC in NC_list:
            E_before_p = patch_energy(S, n1, n2, NC)
            S[n1,n2] = new_spin
            E_after_p = patch_energy(S, n1, n2, NC)

            dE_patch = E_after_p - E_before_p

            S[n1,n2] = old_spin

            errors[NC].append(dE_patch - dE_global)

    return errors


errors = run_test()

for NC, vals in errors.items():
    print("NC =", NC,
          "mean =", np.mean(vals),
          "std =", np.std(vals))