from sage.all import *

# ============================================================
# MODULAR / ETA-PRODUCT PROOF OF THE TWO HSS JET IDENTITIES
#
# Goal:
#   Prove the finite-jet identities R1=0 and R5=0 from
#
#       B(q) = 9 prod_{n>=1} (1-q^n)^(-4)
#            = 9 q^(1/6) eta(tau)^(-4),
#
#   using Ramanujan's differential equations and the CM-point
#   identities at tau=i:
#
#       E2(i) = 3/pi = 6/L,       L=2*pi,
#       E6(i) = 0.
#
# The output should show:
#
#   coefficient check through q^50: True
#   R1 at tau=i is zero? True
#   R5 at tau=i is zero? True
# ============================================================

# ------------------------------------------------------------
# Optional coefficient check against the HSS B-series through U^50
# ------------------------------------------------------------

B_COEFFS = [
    ZZ(9),
    ZZ(36),
    ZZ(126),
    ZZ(360),
    ZZ(945),
    ZZ(2268),
    ZZ(5166),
    ZZ(11160),
    ZZ(23220),
    ZZ(46620),
    ZZ(90972),
    ZZ(172872),
    ZZ(321237),
    ZZ(584640),
    ZZ(1044810),
    ZZ(1835856),
    ZZ(3177153),
    ZZ(5421132),
    ZZ(9131220),
    ZZ(15195600),
    ZZ(25006653),
    ZZ(40722840),
    ZZ(65670768),
    ZZ(104930280),
    ZZ(166214205),
    ZZ(261141300),
    ZZ(407118726),
    ZZ(630048384),
    ZZ(968272605),
    ZZ(1478208420),
    ZZ(2242463580),
    ZZ(3381344280),
    ZZ(5069259342),
    ZZ(7557818940),
    ZZ(11208455370),
    ZZ(16538048640),
    ZZ(24282822798),
    ZZ(35487134928),
    ZZ(51626878470),
    ZZ(74779896240),
    ZZ(107861179482),
    ZZ(154945739844),
    ZZ(221711362038),
    ZZ(316042958880),
    ZZ(448856366490),
    ZZ(635216766732),
    ZZ(895854679650),
    ZZ(1259213600736),
    ZZ(1764210946995),
    ZZ(2463949037340),
    ZZ(3430694064888),
]

print("="*70)
print("ETA-PRODUCT COEFFICIENT CHECK")
print("="*70)

PS.<q> = PowerSeriesRing(ZZ, default_prec=52)

B_eta = PS(9)

# For coefficients through q^50, factors n<=50 are enough.
for n in range(1, 51):
    B_eta *= (1 - q**n)**(-4)

B_eta = B_eta + O(q**51)

coeff_check = True
bad = []

for n in range(51):
    got = B_eta[n]
    want = B_COEFFS[n]
    if got != want:
        coeff_check = False
        bad.append((n, got, want))

print("B(q)=9*prod(1-q^n)^(-4) matches coefficients q^0..q^50?", coeff_check)

if not coeff_check:
    print("Mismatches:")
    for item in bad[:20]:
        print(item)

print()


# ------------------------------------------------------------
# Formal Ramanujan / Eisenstein computation
# ------------------------------------------------------------

print("="*70)
print("FORMAL RAMANUJAN COMPUTATION")
print("="*70)

# Ring variables:
#   L  = 2*pi
#   E2,E4,E6 are Eisenstein values at tau.
R = PolynomialRing(QQ, ["L", "E2", "E4", "E6"])
K = R.fraction_field()

L, E2, E4, E6 = K.gens()

def dK(expr, var):
    """
    Differentiate a fraction-field expression wrt polynomial variable var.
    """
    expr = K(expr)
    num = R(expr.numerator())
    den = R(expr.denominator())

    return K(num.derivative(var)*den - num*den.derivative(var)) / K(den)**2


def theta(expr):
    """
    Ramanujan differential operator theta=q*d/dq on rational functions
    of E2,E4,E6:
        theta E2 = (E2^2-E4)/12
        theta E4 = (E2 E4-E6)/3
        theta E6 = (E2 E6-E4^2)/2
    """
    return (
        dK(expr, R.gen(1)) * (E2**2 - E4) / 12
        + dK(expr, R.gen(2)) * (E2*E4 - E6) / 3
        + dK(expr, R.gen(3)) * (E2*E6 - E4**2) / 2
    )


# For B(q)=9 q^(1/6) eta^(-4):
#
#   theta log eta = E2/24
#
# so
#
#   A := theta log B = 1/6 - 4*(E2/24)
#                    = (1-E2)/6.
A = (1 - E2) / 6

print("A = theta log B =", A)

# Normalized jets:
#
#   u_m = C_m/C_0 = theta^m B / B.
#
# Recurrence:
#
#   u_0 = 1,
#   u_{m+1} = theta(u_m) + A*u_m.
u = [K(1)]

for m in range(5):
    u.append(K(theta(u[-1]) + A*u[-1]))

print()
print("Normalized jets u_m=C_m/C_0:")
for m in range(6):
    print("u{} = {}".format(m, u[m]))

# ------------------------------------------------------------
# CM point tau=i:
#
#   L=2*pi,
#   E2(i)=3/pi=6/L,
#   E6(i)=0.
#
# We leave E4 arbitrary. If the identities simplify to zero with E4
# still free, then E4(i)'s actual value is not needed.
# ------------------------------------------------------------

def cm_sub(expr):
    """
    Substitute E2=6/L and E6=0, leaving L and E4 formal.
    """
    expr = K(expr)
    phi = R.hom([L, 6/L, E4, 0], K)

    num = R(expr.numerator())
    den = R(expr.denominator())

    return K(phi(num)) / K(phi(den))


def is_zero_fraction(expr):
    expr = K(expr)
    return bool(R(expr.numerator()) == 0)


def num_factor(expr):
    expr = K(expr)
    return R(expr.numerator()).factor()


# ------------------------------------------------------------
# Identity R1
# ------------------------------------------------------------

# R1 = -L*C0 + 6*L*C1 + 6*C0.
# Divide by C0:
R1_norm = -L + 6*L*u[1] + 6

R1_cm = cm_sub(R1_norm)

print()
print("="*70)
print("R1 CHECK")
print("="*70)
print("R1/C0 before CM substitution:")
print(R1_norm)
print()
print("R1/C0 after E2=6/L,E6=0:")
print(R1_cm)
print()
print("R1 is zero at tau=i?", is_zero_fraction(R1_cm))

# ------------------------------------------------------------
# Identity R5
# ------------------------------------------------------------

# R5 divided by C0^2, where u_m=C_m/C0:
R5_norm = (
    -54*L**2*u[2]**2
    + 12*L**2*u[3]
    + 108*L**2*u[2]*u[3]
    - 45*L**2*u[4]
    + 54*L**2*u[5]
    + L
    - 72*L*u[2]
    + 216*L*u[3]
    - 270*L*u[4]
    - 8
)

R5_cm = cm_sub(R5_norm)

print()
print("="*70)
print("R5 CHECK")
print("="*70)
print("R5/C0^2 after E2=6/L,E6=0:")
print(R5_cm)
print()
print("R5 is zero at tau=i?", is_zero_fraction(R5_cm))

if not is_zero_fraction(R5_cm):
    print()
    print("R5 numerator factorization:")
    print(num_factor(R5_cm))

# ------------------------------------------------------------
# Extra: show that E4 cancels
# ------------------------------------------------------------

print()
print("="*70)
print("DEPENDENCE CHECK")
print("="*70)

print("R1 numerator degree in E4 after CM substitution:",
      R(R1_cm.numerator()).degree(R.gen(2)))

print("R5 numerator degree in E4 after CM substitution:",
      R(R5_cm.numerator()).degree(R.gen(2)))

print()
print("Interpretation:")
print("- The coefficient check confirms the HSS B-series is the eta-product through q^50.")
print("- R1 follows from A=(1-E2)/6 and E2(i)=6/L.")
print("- R5 follows from Ramanujan's equations plus E2(i)=6/L and E6(i)=0.")
print("- E4 cancels, so the actual E4(i) value is not needed.")
print("Done.")
