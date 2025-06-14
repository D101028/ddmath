from ddmath.matrix import Matrix
from fractions import Fraction

A = Matrix([
    [0,-2,0,4],
    [-2,0,1,1],
    [0,1,1,2],
    [4,1,2,0],
], field=Fraction)

D, E = A.bilinear_diag()

if D is None or E is None:
    exit(0)

print("A =")
print(A)
print()
print("RESULT: D = E^T * A * E, where")
print()
print("D =")
print(D)
print()
print("E =")
print(E)
