from __future__ import annotations
from fractions import Fraction
from typing import TypeVar, Generic, Union, List, Any, Type, overload

from .field import Field
from .quotient import QuotientRing
from .ring import Ring

# 定義係數類型的 TypeVar
T = TypeVar('T', bound=Any)

def _modulus(p_coeffs: list[T], mod_coeffs: list[T], field: Type[T]) -> list[T]:
        """多項式取模運算 (p_coeffs % mod_coeffs)，返回餘式"""
        if len(mod_coeffs) == 0 or all(c == field(0) for c in mod_coeffs):
            raise ZeroDivisionError("modulus by zero polynomial")
        result_coeffs = p_coeffs[:]
        mod_lead = mod_coeffs[-1]
        while len(result_coeffs) >= len(mod_coeffs):
            deg_diff = len(result_coeffs) - len(mod_coeffs)
            lead_coeff = result_coeffs[-1] / mod_lead  # type: ignore
            for i in range(len(mod_coeffs)):
                result_coeffs[deg_diff + i] -= lead_coeff * mod_coeffs[i]  # type: ignore
            while len(result_coeffs) > 0 and result_coeffs[-1] == field(0):
                result_coeffs.pop()
            if not result_coeffs:
                result_coeffs = [field(0)]
        return result_coeffs

class Polyn(Ring, Generic[T]):
    """
    多項式類別，支援泛型係數類型
    
    Args:
        coefficients: 係數列表，從常數項到最高次項 [a₀, a₁, a₂, ...]
                     表示多項式 a₀ + a₁x + a₂x² + ...
        field: 係數使用的數值類型 (如 int, float, Fraction 等)
    """
    
    def __init__(self, coefficients: List[T], field: Type[T]) -> None:
        self.field: Type[T] = field

        if not isinstance(coefficients, list):
            raise TypeError(f"coefficients should be a list, not {type(coefficients).__name__}")
        
        # 將所有係數轉換為指定的 field 類型
        self.coeffs: List[T] = [field(c) if not isinstance(c, field) else c for c in coefficients]  # type: ignore
        self._normalize()
    
    def _normalize(self) -> None:
        """移除最高次項的零係數"""
        zero = self.field(0)  # type: ignore
        while len(self.coeffs) > 1 and self.coeffs[-1] == zero:
            self.coeffs.pop()
    
    @property
    def degree(self) -> int:
        """回傳多項式的次數"""
        return len(self.coeffs) - 1
    
    def __str__(self) -> str:
        """字串表示法"""
        zero = self.field(0)  # type: ignore
        one = self.field(1)  # type: ignore
        neg_one = self.field(-1)  # type: ignore
        
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
            zero = self.field(0)  # type: ignore
            converted_other = self.field(other)  # type: ignore
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
            zero = self.field(0)  # type: ignore
            for i in range(max_len):
                a = self.coeffs[i] if i < len(self.coeffs) else zero
                b = other.coeffs[i] if i < len(other.coeffs) else zero
                new_coeffs.append(a + b)  # type: ignore
            return Polyn(new_coeffs, self.field)
        else:
            # 與純量相加
            new_coeffs = self.coeffs.copy()
            scalar_value = self.field(other)  # type: ignore
            new_coeffs[0] = new_coeffs[0] + scalar_value  # type: ignore
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
            zero = self.field(0)  # type: ignore
            for i in range(max_len):
                a = self.coeffs[i] if i < len(self.coeffs) else zero
                b = other.coeffs[i] if i < len(other.coeffs) else zero
                new_coeffs.append(a - b)  # type: ignore
            return Polyn(new_coeffs, self.field)
        else:
            new_coeffs = self.coeffs.copy()
            scalar_value = self.field(other)  # type: ignore
            new_coeffs[0] = new_coeffs[0] - scalar_value  # type: ignore
            return Polyn(new_coeffs, self.field)
    
    def __rsub__(self, other: Union[T, int, float]) -> 'Polyn[T]':
        """右減法運算"""
        new_coeffs: List[T] = [-c for c in self.coeffs]  # type: ignore
        scalar_value = self.field(other)  # type: ignore
        new_coeffs[0] = new_coeffs[0] + scalar_value  # type: ignore
        return Polyn(new_coeffs, self.field)
    
    @overload
    def __mul__(self, other: 'Polyn[T]') -> 'Polyn[T]': ...
    
    @overload
    def __mul__(self, other: Union[T, int, float]) -> 'Polyn[T]': ...
    
    def __mul__(self, other: Union['Polyn[T]', T, int, float]) -> 'Polyn[T]':
        """乘法運算"""
        if isinstance(other, Polyn):
            # 多項式乘法
            zero = self.field(0)  # type: ignore
            result_coeffs: List[T] = [zero] * (len(self.coeffs) + len(other.coeffs) - 1)
            for i, a in enumerate(self.coeffs):
                for j, b in enumerate(other.coeffs):
                    product = a * b  # type: ignore
                    result_coeffs[i + j] = result_coeffs[i + j] + product  # type: ignore
            return Polyn(result_coeffs, self.field)
        else:
            # 與純量相乘
            scalar_value = self.field(other)  # type: ignore
            new_coeffs: List[T] = [c * scalar_value for c in self.coeffs]  # type: ignore
            return Polyn(new_coeffs, self.field)
    
    def __rmul__(self, other: Union[T, int, float]) -> 'Polyn[T]':
        """右乘法運算"""
        return self.__mul__(other)
    
    def __pow__(self, n: int) -> 'Polyn[T]':
        """冪運算"""
        if n < 0:
            raise ValueError("不支援負指數")
        if n == 0:
            one = self.field(1)  # type: ignore
            return Polyn([one], self.field)
        
        one = self.field(1)  # type: ignore
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
        new_coeffs: List[T] = [-c for c in self.coeffs]  # type: ignore
        return Polyn(new_coeffs, self.field)
    
    def __call__(self, x: Union[T, int, float]) -> T:
        """求值運算，使用 Horner 方法"""
        if not self.coeffs:
            return self.field(0)  # type: ignore
        
        x_value = self.field(x)  # type: ignore
        result: T = self.coeffs[-1]
        for i in range(len(self.coeffs) - 2, -1, -1):
            result = result * x_value + self.coeffs[i]  # type: ignore
        return result
    
    def __mod__(self, mod_poly: 'Polyn[T]') -> 'Polyn[T]':
        """多項式取模運算 (self % mod_poly)，返回餘式"""
        if not isinstance(mod_poly, Polyn):
            raise TypeError("modulus must be a Polyn instance")
        if len(mod_poly.coeffs) == 0 or all(c == self.field(0) for c in mod_poly.coeffs):
            raise ZeroDivisionError("modulus by zero polynomial")
        result_coeffs = _modulus(self.coeffs, mod_poly.coeffs, self.field)
        return Polyn(result_coeffs, self.field)

    def derivative(self) -> 'Polyn[T]':
        """求導數"""
        if len(self.coeffs) <= 1:
            zero = self.field(0)  # type: ignore
            return Polyn([zero], self.field)
        
        new_coeffs: List[T] = []
        for i in range(1, len(self.coeffs)):
            multiplier = self.field(i)  # type: ignore
            new_coeffs.append(self.coeffs[i] * multiplier)  # type: ignore
        
        return Polyn(new_coeffs, self.field)
    
    def integrate(self, constant: Union[T, int, float] = 0) -> 'Polyn[T]':
        """求不定積分"""
        constant_value = self.field(constant)  # type: ignore
        new_coeffs: List[T] = [constant_value]
        for i, coeff in enumerate(self.coeffs):
            divisor = self.field(i + 1)  # type: ignore
            new_coeffs.append(coeff / divisor)  # type: ignore
        
        return Polyn(new_coeffs, self.field)
    
    @classmethod
    def from_roots(
        cls, 
        field: Type[T], 
        roots: List[Union[T, int, float]]
    ) -> 'Polyn[T]':
        """從根建立多項式"""
        one = field(1)  # type: ignore
        result: 'Polyn[T]' = cls([one], field)  # 從 1 開始
        for root in roots:
            # 乘以 (x - root)
            neg_root = field((-field(1)) * field(root))  # type: ignore
            one_coeff = field(1)  # type: ignore
            factor: 'Polyn[T]' = cls([neg_root, one_coeff], field)
            result = result * factor
        return result

class PolynQuotientRing(QuotientRing, Generic[T]):
    def __init__(self, coefficients: List[T], mod_polyn: 'Polyn[T]') -> None:
        polyn = Polyn(coefficients, mod_polyn.field)
        super().__init__(polyn, mod_polyn, mod_polyn.field)
    
    