from typing import Self
from ddmath import Polyn, PolynQuotientRing, IntField, generate_Zn

from fractions import Fraction

Z7 = generate_Zn(7)

a = Z7(5)
b = Z7(4)
c = a / b

print(c ** 114514)



