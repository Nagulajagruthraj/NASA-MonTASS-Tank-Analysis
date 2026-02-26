import numpy as np
from numpy.polynomial.legendre import legder, legval
from scipy.linalg import eigh


# -------------------------------
# Legendre Polynomial Generator
# -------------------------------
def legendre_poly(n, x):
    coeffs = np.zeros(n + 1)
    coeffs[n] = 1
    return legval(x, coeffs)


def legendre_derivative(n, x, order=1):
    coeffs = np.zeros(n + 1)
    coeffs[n] = 1
    for _ in range(order):
        coeffs = legder(coeffs)
    return legval(x, coeffs)


# -------------------------------
# Plate Properties
# -------------------------------
E = 70e9         # Young's modulus (Pa)
nu = 0.33        # Poisson's ratio
rho = 2700       # Density (kg/m3)
t = 0.005        # Thickness (m)
a = 2.0          # Span (m)
b = 0.5          # Chord (m)

D = (E * t ** 3) / (12 * (1 - nu ** 2))

# -------------------------------
# Basis Order
# -------------------------------
N = 3
ndof = (N + 1) ** 2

# -------------------------------
# Gaussian Quadrature
# -------------------------------
gp = np.polynomial.legendre.leggauss(5)

K = np.zeros((ndof, ndof))
M = np.zeros((ndof, ndof))

# -------------------------------
# Assembly
# -------------------------------
index = lambda i, j: i * (N + 1) + j

for xi, wx in zip(*gp):
    for eta, wy in zip(*gp):
        weight = wx * wy

        Nvec = np.zeros(ndof)
        Bvec = np.zeros(ndof)

        for i in range(N + 1):
            for j in range(N + 1):
                idx = index(i, j)

                Pi = legendre_poly(i, xi)
                Pj = legendre_poly(j, eta)

                d2Pi = legendre_derivative(i, xi, 2)
                d2Pj = legendre_derivative(j, eta, 2)

                Nvec[idx] = Pi * Pj

                kx = (2 / a) ** 2 * d2Pi * Pj
                ky = (2 / b) ** 2 * Pi * d2Pj

                Bvec[idx] = kx + ky

        K += D * np.outer(Bvec, Bvec) * weight
        M += rho * t * np.outer(Nvec, Nvec) * weight

# -------------------------------
# Solve Eigenvalue Problem
# -------------------------------
eigvals, eigvecs = eigh(K, M)

frequencies = np.sqrt(np.abs(eigvals)) / (2 * np.pi)

print("Natural Frequencies (Hz):")
print(frequencies[:5])
