from __future__ import annotations
from typing import Self, TypeVar, Generic, Union, List, Any, Type, overload

from math import gcd

from .quotient import QuotientRing
from .ring import Ring

class Field(Ring):
    """Every field inherits this class"""
    def __init__(self, ele) -> None:
        super().__init__(ele)
    
    def __add__(self, other: Self) -> Self:
        ...

    def __sub__(self, other: Self) -> Self:
        ...

    def __mul__(self, other: Self) -> Self:
        ...

    def __neg__(self) -> Self:
        ...

    def __eq__(self, other: Any) -> bool:
        ...

    def __mod__(self, other: Self) -> Self:
        ...

    def __pow__(self, exponent: int) -> Self:
        ...

    def __truediv__(self, other: Self) -> Self:
        ...
    
    def __floordiv__(self, other: Self) -> Self:
        ...

class IntField(QuotientRing, Field):
    def __init__(self, ele: int, mod_ele: int) -> None:
        if gcd(ele, mod_ele) != 1:
            raise ValueError("the element and the quotient must be relatively prime")
        self.ele: int
        self.mod_ele: int
        super().__init__(ele, mod_ele, int)
    
    def _construct(self, ele: int) -> Self:
        return type(self)(ele, self.mod_ele)

    def mul_inv(self) -> Self:
        a: int = self.ele
        m: int = self.mod_ele
        if a == 0:
            raise ZeroDivisionError("0 has no multiplicative inverse")
        # Extended Euclidean Algorithm
        t: int = 0
        new_t: int = 1
        r: int = m
        new_r: int = a
        while new_r != 0:
            quotient: int = r // new_r
            t, new_t = new_t, t - quotient * new_t
            r, new_r = new_r, r - quotient * new_r
        if r > 1:
            raise ValueError(f"{a} has no inverse modulo {m}")
        if t < 0:
            t += m
        return self._construct(t)

    def __pow__(self, exponent: int) -> Self:
        if exponent >= 0:
            return self._construct(pow(self.ele, exponent, self.mod_ele))
        else:
            return self._construct(pow(self.ele, exponent, self.mod_ele)).mul_inv()

    def __truediv__(self, other: Self) -> Self:
        if self.mod_ele != other.mod_ele:
            raise ValueError("Cannot divide elements from different fields")
        return self * other.mul_inv()
    
    def __floordiv__(self, other: Self) -> Self:
        if self.mod_ele != other.mod_ele:
            raise ValueError("Cannot divide elements from different fields")
        return (self - (self % other)) / other

def generate_Zn(n: int):
    if n <= 0:
        raise ValueError("n should be a positive integer")
    class Zn(IntField):
        def __init__(self, ele: int) -> None:
            super().__init__(ele, n)
        
        def _construct(self, ele: int) -> Self:
            return type(self)(ele)
    
    return Zn
