"""Independent exact check for the square-lattice closure.

Requires SymPy.  All calculations are exact over Q(L,b,X), where
L=2*pi formally and b=E4(i).  The script verifies:
  * the normalized square-lattice jets u_0,...,u_5;
  * R1=R3=R5=0;
  * the exact alpha and beta;
  * the generic-X closure residual is identically zero.
"""

import sympy as sp

L, b, X = sp.symbols("L b X", nonzero=True)

u0 = sp.Integer(1)
u1 = (L - 6) / (6 * L)
u2 = (L**2*b + 2*L**2 - 24*L + 36) / (72*L**2)
u3 = (3*L**2*b + 2*L**2 - 36*L + 108) / (432*L**2)
u4 = (3*L**2*b**2 + 3*L**2*b + L**2 - 24*L + 108) / (1296*L**2)
u5 = (
    15*L**2*b**2 + 5*L**2*b + L**2 + 90*L*b**2 - 30*L + 180
) / (7776*L**2)

R1 = sp.factor(-L + 6*L*u1 + 6)
R3 = sp.factor((L - 9) - 54*L*u2 + 108*L*u3)
R5 = sp.factor(
    -54*L**2*u2**2
    + 12*L**2*u3
    + 108*L**2*u2*u3
    - 45*L**2*u4
    + 54*L**2*u5
    + L
    - 72*L*u2
    + 216*L*u3
    - 270*L*u4
    - 8
)

assert R1 == 0
assert R3 == 0
assert R5 == 0

sigma = (
    216*L**3*X**3 + 216*L**3*X**2 - 9*L**3*X*b
    + 72*L**3*X - 3*L**3*b + 8*L**3
    - 216*L**2*X**2 - 144*L**2*X - 6*L**2*b - 24*L**2
    - 756*L*X - 252*L - 648
) / (5832*L**3*(3*X + 1)**2)

Dchan = (
    (L**2*b + 36)*(3*L*X + L + 3)
) / (1944*L**3*(3*X + 1)**2)

Scurv = (
    864*L**5*X**3*b + 864*L**5*X**2*b - 99*L**5*X*b**2
    + 288*L**5*X*b - 33*L**5*b**2 + 32*L**5*b
    - 864*L**4*X**2*b - 576*L**4*X*b - 87*L**4*b**2
    - 96*L**4*b + 10368*L**3*X**3 + 10368*L**3*X**2
    + 1512*L**3*X*b + 3456*L**3*X + 504*L**3*b + 384*L**3
    - 10368*L**2*X**2 - 6912*L**2*X + 2088*L**2*b
    - 1152*L**2 + 84240*L*X + 28080*L + 89424
) / (1679616*L**3*(3*X + 1)**2)

alpha = (L**2*b + 12) / 72
beta = -(7*L**4*b**2 - 552*L**2*b - 13392) / (288*(L**2*b + 36))

residual = sp.factor(sp.cancel(Scurv - alpha*sigma - beta*Dchan))
assert residual == 0

print("R1 = R3 = R5 = 0: verified")
print("generic-X closure residual: 0")
print("alpha =", sp.factor(alpha))
print("beta  =", sp.factor(beta))
