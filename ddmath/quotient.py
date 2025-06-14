from __future__ import annotations
from typing import TypeVar, Generic, Union, List, Any, Type, overload, Self

from .typesetting import Ring, Field

T = TypeVar('T', bound=Ring)
U = TypeVar('U', bound=Field)

class QuotientRing(Ring, Generic[T]):
    def __init__(self, ele: T, mod_ele: T, ring: Type[T]) -> None:
        if mod_ele == ring.add_idn():
            raise ValueError("you can't quotient a zero element")
        # check whether % is available
        if not hasattr(ring, '__mod__') \
            or not callable(getattr(ring, '__mod__')) \
            or ring.__mod__ is Ring.__mod__:
            raise TypeError("ring does not provide a % method")

        self.mod_ele = mod_ele
        self.ele = ele % mod_ele
        self.ring = ring

    def __str__(self) -> str:
        return f"({self.ele} mod {self.mod_ele})"
    def __repr__(self) -> str:
        return f"QuotientRing({self.ele} mod {self.mod_ele})"
    
    def _construct(self, ele: T) -> Self:
        """預設建構方式（讓子類覆寫）"""
        return type(self)(ele, self.mod_ele, self.ring)

    def __eq__(self, other: Self) -> bool:
        if not isinstance(other, QuotientRing):
            return False
        if self.ring != other.ring:
            return False
        if self.mod_ele != other.mod_ele:
            return False
        return self.ele == other.ele
    
    def __ne__(self, other: Self) -> bool:
        return not self.__eq__(other)

    def __add__(self, other: Self) -> Self:
        if not isinstance(other, QuotientRing):
            raise ValueError("It is not a QuotientRing")
        if self.ring != other.ring:
            raise ValueError("They do not under the same ring")
        if self.mod_ele != other.mod_ele:
            raise ValueError("They do not quotient the same element")
        return self._construct((self.ele + other.ele) % self.mod_ele)
        
    def __radd__(self, other: Self) -> Self:
        if not isinstance(other, QuotientRing):
            raise ValueError("It is not a QuotientRing")
        if self.ring != other.ring:
            raise ValueError("They do not under the same ring")
        if self.mod_ele != other.mod_ele:
            raise ValueError("They do not quotient the same element")
        return self._construct((other.ele + self.ele) % self.mod_ele)

    def __sub__(self, other: Self) -> Self:
        if not isinstance(other, QuotientRing):
            raise ValueError("It is not a QuotientRing")
        if self.ring != other.ring:
            raise ValueError("They do not under the same ring")
        if self.mod_ele != other.mod_ele:
            raise ValueError("They do not quotient the same element")
        return self._construct((self.ele - other.ele) % self.mod_ele)

    def __rsub__(self, other: Self) -> Self:
        if not isinstance(other, QuotientRing):
            raise ValueError("It is not a QuotientRing")
        if self.ring != other.ring:
            raise ValueError("They do not under the same ring")
        if self.mod_ele != other.mod_ele:
            raise ValueError("They do not quotient the same element")
        return self._construct((other.ele - self.ele) % self.mod_ele)

    def __mul__(self, other: Self) -> Self:
        if not isinstance(other, QuotientRing):
            raise ValueError("It is not a QuotientRing")
        if self.ring != other.ring:
            raise ValueError("They do not under the same ring")
        if self.mod_ele != other.mod_ele:
            raise ValueError("They do not quotient the same element")
        return self._construct((self.ele * other.ele) % self.mod_ele)

    def __rmul__(self, other: Self) -> Self:
        if not isinstance(other, QuotientRing):
            raise ValueError("It is not a QuotientRing")
        if self.ring != other.ring:
            raise ValueError("They do not under the same ring")
        if self.mod_ele != other.mod_ele:
            raise ValueError("They do not quotient the same element")
        return self._construct((other.ele * self.ele) % self.mod_ele)

    def __neg__(self) -> Self:
        return self._construct((-self.ele) % self.mod_ele)

    def __pow__(self, exponent: int) -> Self:
        if exponent < 0:
            raise ValueError("Exponent must be non-negative")
        result = self.ring.mul_idn()
        base = self.ele
        exp = exponent
        mod = self.mod_ele
        while exp > 0:
            if exp % 2 == 1:
                result = (result * base) % mod
            base = (base * base) % mod
            exp //= 2
        return self._construct(result)

class QuotientField(QuotientRing, Field, Generic[U]):
    def __init__(self, ele: U, mod_ele: U, field: Type[U]) -> None:
        if mod_ele == field.add_idn():
            raise ValueError("you can't quotient a zero element")
        # check whether % is available
        if not hasattr(field, '__mod__') \
            or not callable(getattr(field, '__mod__')) \
            or field.__mod__ is Field.__mod__:
            raise TypeError("field does not provide a % method")
        super().__init__(ele, mod_ele, field)
        self.field = field

