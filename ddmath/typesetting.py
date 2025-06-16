from __future__ import annotations
from typing import Self, Any

NonNegativeInt = int  # type alias for non-negative integers

class Group:
    """Every group inherits this class."""
    def __init__(self, ele: Any) -> None:
        """You should overload this method."""
        raise NotImplementedError
    
    def __add__(self, other: Self) -> Self:
        """You should overload this method."""
        raise NotImplementedError
    
    def __neg__(self) -> Self:
        """You should overload this method."""
        raise NotImplementedError

    def __eq__(self, other: Any) -> bool:
        """You should overload this method."""
        raise NotImplementedError
    
    @classmethod
    def add_idn(cls) -> Self:
        """
        Creates and returns the additive identity (usually zero) for the class.

        Returns:
            Self: An instance of the class initialized to its additive identity (e.g., zero).

        Example:
            >>> MyNumber.add_idn()
            MyNumber(0)

        Note:
            This method assumes that the class constructor accepts a single numeric argument representing the additive identity.
        """
        return cls(0)
    
    def __sub__(self, other: Self) -> Self:
        """
        Implements the subtraction operator for the class.

        This method returns the result of subtracting another instance (`other`) from the current instance (`self`).
        It achieves this by negating `other` and adding it to `self`, assuming that the class defines both `__add__` and unary negation (`__neg__`).

        Args:
            other (Self): The instance to subtract from `self`.

        Returns:
            Self: The result of the subtraction (`self - other`).
        """
        """"""
        return self + (-other)
    
    def __ne__(self, other: Any) -> bool:
        return not self.__eq__(other)
    
    def __lmul__(self, other: int) -> Self: 
        """
        Implements left multiplication for the object, allowing expressions like `3 * g` where `g` is an instance of this class.
        This method should return a new instance representing the sum of `g` added to itself `3` times (i.e., repeated addition).
        Args:
            other (int): The integer multiplier on the left side of the multiplication.
        Returns:
            Self: A new instance representing the result of the repeated addition.
        Raises:
            TypeError: If `other` is not an integer.
        Example:
            If `g` is an instance, `3 * g` will call `g.__lmul__(3)` and return `g + g + g`.
        """
        """3 * g means g + g + g"""
        if not isinstance(other, int):
            raise TypeError(f"Left multiplier must be int, got {type(other).__name__}")
        if other < 0:
            result = -self.add_idn()
        else:
            result = self.add_idn()
        base = self
        exp = other
        while exp > 0:
            if exp & 1:
                result = result + base
            base = base + base
            exp >>= 1
        return result

    def __mul__(self, other: int) -> Self:
        """
        Implements right multiplication for the object, allowing expressions like `g * 3` where `g` is an instance of this class.
        This method should return a new instance representing the sum of `g` added to itself `3` times (i.e., repeated addition).
        Args:
            other (int): The integer multiplier on the right side of the multiplication.
        Returns:
            Self: A new instance representing the result of the repeated addition.
        Raises:
            TypeError: If `other` is not an integer.
        Example:
            If `g` is an instance, `g * 3` will call `g.__mul__(3)` and return `g + g + g`.
        """
        if not isinstance(other, int):
            raise TypeError(f"Right multiplier must be int, got {type(other).__name__}")
        return self.__lmul__(other)

class Ring(Group):
    """Every ring inherits this class."""
    def __init__(self, ele: Any) -> None:
        """You should overload this method."""
        raise NotImplementedError

    def __mul__(self, other: Self | int) -> Self:
        """You should overload some parts of this method."""
        if isinstance(other, int):
            return super().__mul__(other)
        elif isinstance(other, type(self)):
            # You should overload this method in subclasses to define multiplication between two elements
            raise NotImplementedError(f"{type(self).__name__} does not implement multiplication with its own type.")
        else:
            raise TypeError(f"Unsupported operand type(s) for *: '{type(self).__name__}' and '{type(other).__name__}'")
    
    def __lmul__(self, other: Self | int) -> Self:
        """You should overload some parts of this method."""
        if isinstance(other, int):
            return super().__lmul__(other)
        elif isinstance(other, type(self)):
            # You should overload this method in subclasses to define multiplication between two elements
            raise NotImplementedError(f"{type(self).__name__} does not implement multiplication with its own type.")
        else:
            raise TypeError(f"Unsupported operand type(s) for *: '{type(self).__name__}' and '{type(other).__name__}'")
    
    def __mod__(self, other: Self) -> Self:
        """You may overload this method."""
        raise NotImplementedError(f"{type(self).__name__} does not support '%' operation.")

    def __truediv__(self, other: Self) -> Self:
        """You may overload this method."""
        raise NotImplementedError(f"{type(self).__name__} does not support '//' operation.")

    def __floordiv__(self, other: Self) -> Self:
        """You may overload this method."""
        raise NotImplementedError(f"{type(self).__name__} does not support '//' operation.")

    def mul_inv(self) -> Self:
        """You may overload this method."""
        raise NotImplementedError(f"{type(self).__name__} does not support 'mul_inv'.")

    @classmethod
    def mul_idn(cls) -> Self:
        """
        Creates and returns the multiplicative identity (usually one) for the class.

        Returns:
            Self: An instance of the class initialized to its multiplicative identity (e.g., one).

        Example:
            >>> MyNumber.mul_idn()
            MyNumber(1)

        Note:
            This method assumes that the class constructor accepts a single numeric argument representing the multiplicative identity.
        """
        return cls(1)

    def _mul_int(self, other: int) -> Self:
        if not isinstance(other, int):
            raise TypeError(f"Left multiplier must be int, got {type(other).__name__}")
        if other < 0:
            result = -self.add_idn()
        else:
            result = self.add_idn()
        base = self
        exp = other
        while exp > 0:
            if exp & 1:
                result = result + base
            base = base + base
            exp >>= 1
        return result
    
    def __pow__(self, exponent: int) -> Self:
        if not isinstance(exponent, int):
            raise TypeError(f"Exponent must be int, got {type(exponent).__name__}")
        if exponent < 0:
            return self.mul_inv().__pow__(-exponent)
        result = self.mul_idn()
        base = self
        exp = exponent
        while exp > 0:
            if exp & 1:
                result = result * base
            base = base * base
            exp >>= 1
        return result
    
class Field(Ring):
    """Every field inherits this class"""
    def mul_inv(self) -> Self:
        """You should overload this method."""
        raise NotImplementedError

    def __truediv__(self, other: Self) -> Self:
        """You should overload this method."""
        raise NotImplementedError
    
    def __floordiv__(self, other: Self) -> Self:
        """You may overload this method."""
        raise NotImplementedError
    