from __future__ import annotations
from fractions import Fraction
from typing import TypeVar, Generic, Union, List, Any, Type, overload

from math import gcd as _gcd
from collections.abc import Iterable

from .quotient import QuotientRing
from .typesetting import Ring, Field

# 定義係數類型的 TypeVar
T = TypeVar('T', bound=Field)

def _modulus(p_coeffs: list[T], mod_coeffs: list[T], field: Type[T]) -> list[T]:
        """多項式取模運算 (p_coeffs % mod_coeffs)，返回餘式"""
        if len(mod_coeffs) == 0 or all(c == field.add_idn() for c in mod_coeffs):
            raise ZeroDivisionError("modulus by zero polynomial")
        result_coeffs = p_coeffs[:]
        mod_lead = mod_coeffs[-1]
        while len(result_coeffs) >= len(mod_coeffs):
            deg_diff = len(result_coeffs) - len(mod_coeffs)
            lead_coeff = result_coeffs[-1] / mod_lead
            for i in range(len(mod_coeffs)):
                result_coeffs[deg_diff + i] -= lead_coeff * mod_coeffs[i]
            while len(result_coeffs) > 0 and result_coeffs[-1] == field.add_idn():
                result_coeffs.pop()
            if not result_coeffs:
                result_coeffs = [field.add_idn()]
        return result_coeffs

def gcd(p1: Polyn[T] | int, p2: Polyn[T] | int) -> Polyn[T] | int:
    if isinstance(p1, int) and isinstance(p2, int):
        return _gcd(p1, p2)
    if isinstance(p1, int) and isinstance(p2, Polyn):
        p1 = Polyn([p1], p2.field)
    if isinstance(p2, int) and isinstance(p1, Polyn):
        p2 = Polyn([p2], p1.field)
    if not (isinstance(p1, Polyn) and isinstance(p2, Polyn)):
        raise TypeError("Both arguments must be either integers or one integer and one Polyn instance")
    if p1.field != p2.field:
        raise TypeError("Polynomials must have the same coefficient field")
    a, b = p1, p2
    zero = a.field.add_idn()
    while not (len(b.coeffs) == 1 and b.coeffs[0] == zero):
        a, b = b, a % b
    # Normalize so leading coefficient is 1 if possible
    lead = a.coeffs[-1]
    if hasattr(lead, '__truediv__') and lead != zero:
        a = a * (a.field.mul_idn() / lead)
    return a

class Polyn(Ring, Generic[T]):
    """
    Polyn is a generic polynomial class supporting coefficients of any numeric type.
    This class provides arithmetic operations, evaluation, differentiation, integration, and construction from roots for polynomials with coefficients of a specified numeric type (such as int, float, Fraction, etc.).
        coefficients (List[T]): List of coefficients, from constant term to highest degree term [a₀, a₁, a₂, ...], representing the polynomial a₀ + a₁x + a₂x² + ...
        field (Type[T]): The numeric type used for coefficients (e.g., int, float, Fraction).
    Attributes:
        coeffs (List[T]): The normalized list of coefficients.
        field (Type[T]): The numeric type used for coefficients.
    Methods:
        degree: Returns the degree of the polynomial.
        __str__: Returns a human-readable string representation of the polynomial.
        __repr__: Returns a string representation suitable for debugging.
        __eq__, __ne__: Equality and inequality comparison with another polynomial or scalar.
        __add__, __radd__: Addition with another polynomial or scalar.
        __sub__, __rsub__: Subtraction with another polynomial or scalar.
        __mul__, __rmul__: Multiplication with another polynomial or scalar.
        __pow__: Exponentiation (non-negative integer powers).
        __neg__: Negation of the polynomial.
        __call__: Evaluates the polynomial at a given value using Horner's method.
        __mod__: Polynomial modulus operation, returning the remainder after division by another polynomial.
        derivative: Returns the derivative of the polynomial.
        integrate: Returns the indefinite integral (antiderivative) of the polynomial, with an optional constant of integration.
        from_roots: Class method to construct a polynomial given its roots.
    Raises:
        TypeError: If coefficients is not a list or if modulus is not a Polyn instance.
        ZeroDivisionError: If modulus by the zero polynomial is attempted.
        ValueError: If negative exponent is used in __pow__.
    Examples:
        >>> p = Polyn([1, 2, 3], int)  # Represents 1 + 2x + 3x^2
        >>> str(p)
        '(1 + 2x + 3x^2)'
        >>> p(2)
        17
        >>> p.derivative()
        Polyn(int, [2, 6])
    """
    def __init__(self, coefficients: List[Any] | Any, field: Type[T]) -> None:
        self.field: Type[T] = field

        if not isinstance(coefficients, list):
            if not isinstance(coefficients, Iterable) or isinstance(coefficients, (str, bytes)):
                coefficients = [coefficients]
            else:
                coefficients = list(coefficients)
        
        # 將所有係數轉換為指定的 field 類型
        self.coeffs: List[T] = [field(c) if not isinstance(c, field) else c for c in coefficients]
        self._normalize()
    
    def _normalize(self) -> None:
        """移除最高次項的零係數"""
        zero = self.field.add_idn()
        while len(self.coeffs) > 1 and self.coeffs[-1] == zero:
            self.coeffs.pop()
    
    @property
    def degree(self) -> int:
        """回傳多項式的次數"""
        return len(self.coeffs) - 1
    
    def __str__(self) -> str:
        """字串表示法"""
        zero = self.field.add_idn()
        one = self.field.mul_idn()
        neg_one = self.field(-1)
        
        if all(c == zero for c in self.coeffs):
            return "0"
        
        terms: List[str] = []
        for i, coeff in enumerate(self.coeffs):
            if coeff == zero:
                continue
            
            # 係數部分
            if i == 0:  # 常數項
                terms.append(str(coeff))
            elif coeff == one:
                if i == 1:
                    terms.append("x")
                else:
                    terms.append(f"x^{i}")
            elif coeff == neg_one:
                if i == 1:
                    terms.append("-x")
                else:
                    terms.append(f"-x^{i}")
            else:
                if i == 1:
                    terms.append(f"{coeff}x")
                else:
                    terms.append(f"{coeff}x^{i}")
        
        if not terms:
            return "0"
        
        # 組合項目，處理正負號
        result = terms[0]
        for term in terms[1:]:
            if term.startswith('-'):
                result += " - " + term[1:]
            else:
                result += " + " + term
        
        return f"({result})"
    
    def __repr__(self) -> str:
        return f"Polyn({self.field.__name__}, {self.coeffs})"
    
    def __eq__(self, other: Any) -> bool:
        """相等比較"""
        if isinstance(other, Polyn):
            return self.coeffs == other.coeffs
        elif isinstance(other, (int, float, Fraction)):
            zero = self.field.add_idn()
            converted_other = self.field(other)
            return len(self.coeffs) == 1 and self.coeffs[0] == converted_other
        return False
    
    def __ne__(self, other: Any) -> bool:
        """不等比較"""
        return not self.__eq__(other)
    
    @overload
    def __add__(self, other: 'Polyn[T]') -> 'Polyn[T]': ...
    
    @overload
    def __add__(self, other: Union[T, int, float]) -> 'Polyn[T]': ...
    
    def __add__(self, other: Union['Polyn[T]', T, int, float]) -> 'Polyn[T]':
        """加法運算"""
        if isinstance(other, Polyn):
            # 確保兩個多項式有相同的長度
            max_len = max(len(self.coeffs), len(other.coeffs))
            new_coeffs: List[T] = []
            zero = self.field.add_idn()
            for i in range(max_len):
                a = self.coeffs[i] if i < len(self.coeffs) else zero
                b = other.coeffs[i] if i < len(other.coeffs) else zero
                new_coeffs.append(a + b)
            return Polyn(new_coeffs, self.field)
        else:
            # 與純量相加
            new_coeffs = self.coeffs.copy()
            scalar_value = self.field(other)
            new_coeffs[0] = new_coeffs[0] + scalar_value
            return Polyn(new_coeffs, self.field)
    
    def __radd__(self, other: Union[T, int, float]) -> 'Polyn[T]':
        """右加法運算"""
        return self.__add__(other)
    
    @overload
    def __sub__(self, other: 'Polyn[T]') -> 'Polyn[T]': ...
    
    @overload
    def __sub__(self, other: Union[T, int, float]) -> 'Polyn[T]': ...
    
    def __sub__(self, other: Union['Polyn[T]', T, int, float]) -> 'Polyn[T]':
        """減法運算"""
        if isinstance(other, Polyn):
            max_len = max(len(self.coeffs), len(other.coeffs))
            new_coeffs: List[T] = []
            zero = self.field.add_idn()
            for i in range(max_len):
                a = self.coeffs[i] if i < len(self.coeffs) else zero
                b = other.coeffs[i] if i < len(other.coeffs) else zero
                new_coeffs.append(a - b)
            return Polyn(new_coeffs, self.field)
        else:
            new_coeffs = self.coeffs.copy()
            scalar_value = self.field(other)
            new_coeffs[0] = new_coeffs[0] - scalar_value
            return Polyn(new_coeffs, self.field)
    
    def __rsub__(self, other: Union[T, int, float]) -> 'Polyn[T]':
        """右減法運算"""
        new_coeffs: List[T] = [-c for c in self.coeffs]
        scalar_value = self.field(other)
        new_coeffs[0] = new_coeffs[0] + scalar_value
        return Polyn(new_coeffs, self.field)
    
    @overload
    def __mul__(self, other: 'Polyn[T]') -> 'Polyn[T]': ...
    
    @overload
    def __mul__(self, other: Union[T, int, float]) -> 'Polyn[T]': ...
    
    def __mul__(self, other: Union['Polyn[T]', T, int, float]) -> 'Polyn[T]':
        """乘法運算"""
        if isinstance(other, Polyn):
            # 多項式乘法
            zero = self.field.add_idn()
            result_coeffs: List[T] = [zero] * (len(self.coeffs) + len(other.coeffs) - 1)
            for i, a in enumerate(self.coeffs):
                for j, b in enumerate(other.coeffs):
                    product = a * b
                    result_coeffs[i + j] = result_coeffs[i + j] + product
            return Polyn(result_coeffs, self.field)
        else:
            # 與純量相乘
            scalar_value = self.field(other)
            new_coeffs: List[T] = [c * scalar_value for c in self.coeffs]
            return Polyn(new_coeffs, self.field)
    
    def __rmul__(self, other: Union[T, int, float]) -> 'Polyn[T]':
        """右乘法運算"""
        return self.__mul__(other)
    
    def __pow__(self, n: int) -> 'Polyn[T]':
        """冪運算"""
        if n < 0:
            raise ValueError("不支援負指數")
        if n == 0:
            one = self.field.mul_idn()
            return Polyn([one], self.field)
        
        one = self.field.mul_idn()
        result: 'Polyn[T]' = Polyn([one], self.field)
        base: 'Polyn[T]' = self
        
        # 快速冪演算法
        while n > 0:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        
        return result
    
    def __neg__(self) -> 'Polyn[T]':
        """負號運算"""
        new_coeffs: List[T] = [-c for c in self.coeffs]
        return Polyn(new_coeffs, self.field)
    
    def __call__(self, x: Union[T, int, float]) -> T:
        """求值運算，使用 Horner 方法"""
        if not self.coeffs:
            return self.field.add_idn()
        
        x_value = self.field(x)
        result: T = self.coeffs[-1]
        for i in range(len(self.coeffs) - 2, -1, -1):
            result = result * x_value + self.coeffs[i]
        return result
    
    def __mod__(self, mod_poly: 'Polyn[T]') -> 'Polyn[T]':
        """多項式取模運算 (self % mod_poly)，返回餘式"""
        if not isinstance(mod_poly, Polyn):
            raise TypeError("modulus must be a Polyn instance")
        if len(mod_poly.coeffs) == 0 or all(c == self.field.add_idn() for c in mod_poly.coeffs):
            raise ZeroDivisionError("modulus by zero polynomial")
        result_coeffs = _modulus(self.coeffs, mod_poly.coeffs, self.field)
        return Polyn(result_coeffs, self.field)

    def derivative(self) -> 'Polyn[T]':
        """求導數"""
        if len(self.coeffs) <= 1:
            zero = self.field.add_idn()
            return Polyn([zero], self.field)
        
        new_coeffs: List[T] = []
        for i in range(1, len(self.coeffs)):
            multiplier = self.field(i)
            new_coeffs.append(self.coeffs[i] * multiplier)
        
        return Polyn(new_coeffs, self.field)
    
    def integrate(self, constant: Union[T, int, float] = 0) -> 'Polyn[T]':
        """求不定積分"""
        constant_value = self.field(constant)
        new_coeffs: List[T] = [constant_value]
        for i, coeff in enumerate(self.coeffs):
            divisor = self.field(i + 1)
            new_coeffs.append(coeff / divisor)
        
        return Polyn(new_coeffs, self.field)
    
    @classmethod
    def from_roots(
        cls, 
        field: Type[T], 
        roots: List[Union[T, int, float]]
    ) -> 'Polyn[T]':
        """從根建立多項式"""
        one = field.mul_idn()
        result: 'Polyn[T]' = cls([one], field)  # 從 1 開始
        for root in roots:
            # 乘以 (x - root)
            neg_root = field((-field.mul_idn()) * field(root))
            one_coeff = field.mul_idn()
            factor: 'Polyn[T]' = cls([neg_root, one_coeff], field)
            result = result * factor
        return result

class PolynQuotientRing(QuotientRing, Generic[T]):
    def __init__(self, coefficients: List[T], mod_polyn: 'Polyn[T]') -> None:
        polyn = Polyn(coefficients, mod_polyn.field)
        class TempPolynRing(Polyn):
            def __init__(self, coefficients: List[Any] | Any) -> None:
                super().__init__(coefficients, mod_polyn.field)
        super().__init__(polyn, mod_polyn, TempPolynRing)


