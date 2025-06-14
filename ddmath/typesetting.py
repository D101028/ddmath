from __future__ import annotations
from typing import Self, Any

NonNegativeInt = int  # type alias for non-negative integers

class Ring:
    """Every ring inherits this class."""
    def __init__(self, ele: Any) -> None:
        raise NotImplementedError

    @classmethod
    def add_idn(cls) -> Self:
        return cls(0)
    
    @classmethod
    def mul_idn(cls) -> Self:
        return cls(1)
    
    def __add__(self, other: Self) -> Self:
        raise NotImplementedError

    def __sub__(self, other: Self) -> Self:
        raise NotImplementedError

    def __mul__(self, other: Self) -> Self:
        raise NotImplementedError

    def __neg__(self) -> Self:
        raise NotImplementedError

    def __eq__(self, other: Any) -> bool:
        raise NotImplementedError
    
    def __neq__(self, other: Any) -> bool:
        raise NotImplementedError
        
    def __pow__(self, exponent: NonNegativeInt) -> Self:
        raise NotImplementedError
    
    def __mod__(self, other: Self) -> Self:
        raise NotImplementedError(f"{type(self).__name__} does not support '%' operation.")

class Field(Ring):
    """Every field inherits this class"""
    def __pow__(self, exponent: int) -> Self:
        raise NotImplementedError

    def __truediv__(self, other: Self) -> Self:
        raise NotImplementedError
    
    def __floordiv__(self, other: Self) -> Self:
        raise NotImplementedError

    def __mod__(self, other: Self) -> Self:
        raise NotImplementedError(f"{type(self).__name__} does not support '%' operation.")
