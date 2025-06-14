from __future__ import annotations
from typing import Any

def gcd(a: int, b: int):
    assert isinstance(a, int) and isinstance(b, int) and a >= 0 and b >= 0
    while True:
        if a >= b:
            if b == 0:
                return a
            a %= b
        else:
            if a == 0:
                return b
            b %= a

def reduced(fenzu: int, fenmu: int):
    assert isinstance(fenzu, int) and isinstance(fenmu, int)
    if fenmu == 0:
        raise ValueError("int can not be devided by zero")
    g = gcd(abs(fenzu), abs(fenmu))
    if g == 1:
        return (fenzu, fenmu)
    fenzu //= g
    fenmu //= g
    return (fenzu, fenmu)
    
class Frac:
    def __init__(self, p1: int | Frac | list[int] | tuple[int], p2: int | None = None):
        """p1/p2"""
        self.tup: tuple[int, int]
        if not (isinstance(p2, int) or p2 is None):
            raise TypeError(f"p2 should be int or None, not {type(p2).__name__}")
        if p2 == 0:
            raise ValueError("int can not be devided by zero")
        if isinstance(p1, int):
            if p2 is None:
                self.tup = (p1, 1)
            else:
                tup = reduced(p1, p2)
                if tup[1] < 0:
                    self.tup = (-tup[0], -tup[1])
                else:
                    self.tup = tup
        elif isinstance(p1, Frac):
            if p2 is None:
                self.tup = p1.tup
            else:
                self.tup = (p1 / p2).tup
        elif isinstance(p1, list) or isinstance(p1, tuple):
            if len(p1) != 2:
                raise ValueError(f"doesn't have correct length: {len(p1)}")
            if p2 is None:
                tup = reduced(p1[0], p1[1])
                if tup[1] < 0:
                    self.tup = (-tup[0], -tup[1])
                else:
                    self.tup = tup
            else:
                self.tup = (Frac(p1[0], p1[1]) / Frac(p2, 1)).tup
        else:
            raise TypeError(f"incorrect type: {type(p1).__name__}")
    
    def __float__(self):
        return self.tup[0] / self.tup[1]
    def __int__(self):
        return self.tup[0] // self.tup[1]
    def __add__(self, other: Frac | int):
        if isinstance(other, int):
            tup2 = (other, 1)
        elif isinstance(other, Frac):
            tup2 = other.tup
        else:
            raise TypeError(f"incorrect type: {type(other).__name__}")
        return Frac(self.tup[0] * tup2[1] + self.tup[1] * tup2[0], self.tup[1] * tup2[1])
    def __sub__(self, other: Frac | int):
        if isinstance(other, int):
            tup2 = (other, 1)
        elif isinstance(other, Frac):
            tup2 = other.tup
        else:
            raise TypeError(f"incorrect type: {type(other).__name__}")
        return Frac(self.tup[0] * tup2[1] - self.tup[1] * tup2[0], self.tup[1] * tup2[1])
    def __mul__(self, other: Frac | int):
        if isinstance(other, int):
            tup2 = (other, 1)
        elif isinstance(other, Frac):
            tup2 = other.tup
        else:
            raise TypeError(f"incorrect type: {type(other).__name__}")
        return Frac(self.tup[0] * tup2[0], self.tup[1] * tup2[1])
    def __truediv__(self, other: Frac | int):
        if isinstance(other, int):
            tup2 = (other, 1)
        elif isinstance(other, Frac):
            tup2 = other.tup
        else:
            raise TypeError(f"incorrect type: {type(other).__name__}")
        if tup2[0] == 0:
            raise ValueError("Frac should not be divided by zero")
        return Frac(self.tup[0] * tup2[1], self.tup[1] * tup2[0])
    def __radd__(self, other: Frac | int):
        if isinstance(other, int):
            tup2 = (other, 1)
        elif isinstance(other, Frac):
            tup2 = other.tup
        else:
            raise TypeError(f"incorrect type: {type(other).__name__}")
        return Frac(self.tup[0] * tup2[1] + self.tup[1] * tup2[0], self.tup[1] * tup2[1])
    def __rsub__(self, other: Frac | int):
        if isinstance(other, int):
            tup2 = (other, 1)
        elif isinstance(other, Frac):
            tup2 = other.tup
        else:
            raise TypeError(f"incorrect type: {type(other).__name__}")
        return Frac(self.tup[1] * tup2[0] - self.tup[0] * tup2[1], self.tup[1] * tup2[1])
    def __rmul__(self, other: Frac | int):
        if isinstance(other, int):
            tup2 = (other, 1)
        elif isinstance(other, Frac):
            tup2 = other.tup
        else:
            raise TypeError(f"incorrect type: {type(other).__name__}")
        return Frac(self.tup[0] * tup2[0], self.tup[1] * tup2[1])
    def __rtruediv__(self, other: Frac | int):
        if isinstance(other, int):
            tup2 = (other, 1)
        elif isinstance(other, Frac):
            tup2 = other.tup
        else:
            raise TypeError(f"incorrect type: {type(other).__name__}")
        if self.tup[0] == 0:
            raise ValueError("Frac should not be divided by zero")
        return Frac(tup2[0] * self.tup[1], tup2[1] * self.tup[0])
    def __pow__(self, other: int):
        if not isinstance(other, int):
            raise TypeError(f"incorrect type: {type(other).__name__}")
        if other >= 0:
            return Frac(self.tup[0] ** other, self.tup[1] ** other)
        else:
            if self.tup[0] == 0:
                return ValueError("Frac should not be devided by zero")
            other = -other
            return Frac(self.tup[1] ** other, self.tup[0] ** other)
    def __eq__(self, other: Frac | int):
        if isinstance(other, int):
            tup2 = (other, 1)
        elif isinstance(other, Frac):
            tup2 = other.tup
        else:
            raise TypeError(f"incorrect type: {type(other).__name__}")
        return self.tup == tup2
    def __gt__(self, other: Frac | Any):
        if not isinstance(other, Frac):
            return self.__float__() > other
        return self.tup[0] * other.tup[1] > self.tup[1] * other.tup[0]
    def __ge__(self, other: Frac | Any):
        if not isinstance(other, Frac):
            return self.__float__() >= other
        return self.tup[0] * other.tup[1] >= self.tup[1] * other.tup[0]
    def __lt__(self, other: Frac | Any):
        if not isinstance(other, Frac):
            return self.__float__() < other
        return self.tup[0] * other.tup[1] < self.tup[1] * other.tup[0]
    def __le__(self, other: Frac | Any):
        if not isinstance(other, Frac):
            return self.__float__() <= other
        return self.tup[0] * other.tup[1] <= self.tup[1] * other.tup[0]
    def __str__(self):
        if self.tup[1] == 1:
            return str(self.tup[0])
        return f"{self.tup[0]}/{self.tup[1]}"

