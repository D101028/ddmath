from __future__ import annotations
from typing import Self, Any

NonNegativeInt = int  # type alias for non-negative integers

class Ring:
    """Every ring inherits this class."""
    def __init__(self, ele: Any) -> None:
        ...
    
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
    
    def __neq__(self, other: Any) -> bool:
        ...

    def __mod__(self, other: Self) -> Self:
        ...
        
    def __pow__(self, exponent: NonNegativeInt) -> Self:
        ...

class Field(Ring):
    """Every field inherits this class"""
    def __init__(self, ele: Any) -> None:
        ...
    
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
