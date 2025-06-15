from __future__ import annotations
from typing import Any, Callable, Type, TypeVar, Generic, Iterable, Sized

from .typesetting import Field

# Define the TypeVar for matrix elements
T = TypeVar('T', bound=Field)

def myLen(obj: Iterable) -> int:
    if isinstance(obj, Sized):
        return len(obj)
    else:
        count = 0
        for _ in obj:
            count += 1
        return count

def check_arr_form(arr: Iterable[Iterable[Any]]) -> bool:
    length1 = myLen(arr)
    if length1 == 0:
        return True
    it = iter(arr)
    a0 = next(it)
    length = myLen(a0)
    for _ in range(length1-1):
        if myLen(next(it)) != length:
            return False
    return True

def matrix_output(arg: list[list]) -> str:
    pos = 0
    length = len(arg)
    lengths = [len(a) for a in arg]
    temp: list[list[str]] = [[] for _ in range(length)]
    while True:
        width = 0
        is_end = True
        for i in range(length):
            if lengths[i] <= pos:
                continue
            is_end = False
            ch = str(arg[i][pos])
            if len(ch) > width:
                new_width = len(ch)
                for j in range(i):
                    if lengths[j] <= pos:
                        continue
                    temp[j][pos] = " "*(new_width - width) + temp[j][pos]
                width = new_width
            else:
                ch = " "*(width - len(ch)) + ch
            temp[i].append(ch)
        if is_end:
            break
        else:
            pos += 1
    return "\n".join(" ".join(a) for a in temp)

class Matrix(Generic[T]):
    def __init__(self, 
                 arr: Iterable[Iterable[Any]], 
                 field: Type[T]) -> None:
        if not check_arr_form(arr):
            raise ValueError("incorrect matrix elements input")

        self.field = field

        self.arr: list[list[T]] = [[self.field(j) for j in a] for a in arr]
        self.add_idn: T = self.field(0) if not hasattr(self.field, 'add_idn') else self.field.add_idn()
        self.mul_idn: T = self.field(1) if not hasattr(self.field, 'mul_idn') else self.field.mul_idn()

    def __str__(self) -> str:
        return matrix_output(self.arr)
    def __getitem__(self, index: int) -> list[Any]:
        if not isinstance(index, int):
            raise TypeError(f"incorrect type: {type(index).__name__}")
        return self.arr[index]
    
    def __add__(self, other: Matrix[T]) -> Matrix[T]:
        if not isinstance(other, Matrix):
            raise TypeError(f"incorrect type: {type(other).__name__}")
        if self.get_shape() != other.get_shape():
            raise ValueError("error shapes for adding")
        if self.add_idn != other.add_idn or self.mul_idn != other.mul_idn:
            raise ValueError("different fields")
        rows, cols = self.get_shape()
        arr: list[list[T | None]] = [[None]*cols for _ in range(rows)]
        for r in range(rows):
            for c in range(cols):
                arr[r][c] = self.arr[r][c] + other.arr[r][c]
        return Matrix(arr, self.field)
    def __sub__(self, other: Matrix[T]) -> Matrix[T]:
        if not isinstance(other, Matrix):
            raise TypeError(f"incorrect type: {type(other).__name__}")
        if self.get_shape() != other.get_shape():
            raise ValueError("error shapes for adding")
        if self.add_idn != other.add_idn or self.mul_idn != other.mul_idn:
            raise ValueError("different fields")
        rows, cols = self.get_shape()
        arr: list[list[T | None]] = [[None]*cols for _ in range(rows)]
        for r in range(rows):
            for c in range(cols):
                arr[r][c] = self.arr[r][c] - other.arr[r][c]
        return Matrix(arr, self.field)
    def __mul__(self, other: Matrix[T] | Any) -> Matrix[T]:
        """scalar multiplication or matrix multiplication (self * other)"""
        if not isinstance(other, Matrix):
            # scalar multiplication
            m = self.copy()
            for r, _ in enumerate(self.arr):
                m.row_operation_2(r, other, True)
            return m
        if self.add_idn != other.add_idn or self.mul_idn != other.mul_idn:
            raise ValueError("different fields")
        rows1, cols1 = self.get_shape()
        rows2, cols2 = other.get_shape()
        if cols1 != rows2:
            raise ValueError("error shapes for multiplying")
        arr: list[list[Any]] = [[None] * cols2 for _ in range(rows1)]
        for r in range(rows1):
            for c in range(cols2):
                a = self.add_idn
                for i in range(cols1):
                    a += self.arr[r][i] * other.arr[i][c]
                arr[r][c] = a
        return Matrix(arr, self.field)
    def __rmul__(self, other: Any) -> Matrix[T]:
        if isinstance(other, Matrix):
            return NotImplemented
        return self * other
    def __pow__(self, other: int) -> Matrix[T]:
        if not isinstance(other, int):
            raise TypeError(f"incorrect type: {type(other).__name__}")
        rows, cols = self.get_shape()
        if rows != cols:
            raise ValueError("It should be a square matrix")
        result = self.identity()
        # result = Matrix([[0]*i + [1] + [0]*(cols - i - 1) for i in range(rows)])
        if other >= 0:
            current_base = self.copy()
            exponent = other
            while exponent > 0:
                if exponent % 2 == 1:
                    result *= current_base
                current_base *= current_base
                exponent //= 2
            return result
        else:
            current_base = self.inverse()
            exponent = -other
            while exponent > 0:
                if exponent % 2 == 1:
                    result *= current_base
                current_base *= current_base
                exponent //= 2
            return result
    def __eq__(self, other: Matrix[T]) -> bool:
        if not isinstance(other, Matrix):
            return False
        if self.get_shape() != other.get_shape():
            return False
        if self.add_idn != other.add_idn or self.mul_idn != other.mul_idn:
            return False
        rows, cols = self.get_shape()
        for r in range(rows):
            for c in range(cols):
                if self.arr[r][c] != other.arr[r][c]:
                    return False
        return True
    def __ne__(self, other: Matrix[T]) -> bool:
        return not self == other
    
    def identity(self, n: int | None = None) -> Matrix[T]:
        if n is None:
            rows, cols = self.get_shape()
            if rows != cols:
                raise ValueError("It should be a square matrix")
            n = rows
        if not isinstance(n, int):
            raise TypeError(f"incorrect type: {type(n).__name__}")
        if n <= 0:
            raise ValueError("n should be a positive integer")
        return Matrix([[self.add_idn]*i + [self.mul_idn] + [self.add_idn]*(n - i - 1) for i in range(n)], self.field)

    def append_row(self, arr: list[Any] | tuple[Any]) -> None:
        if not isinstance(arr, list) and not isinstance(arr, tuple):
            raise TypeError("incorrect type for row")
        rows, cols = self.get_shape()
        if len(arr) != cols:
            raise ValueError("incorrect columns number for the matrix")
        new_row = []
        for i in arr:
            if not isinstance(i, type(self.add_idn)):
                raise ValueError("incorrect type in arr for matrix")
            new_row.append(i)
        self.arr.append(new_row)
    
    def append_col(self, arr: list[Any] | tuple[Any]) -> None:
        if not isinstance(arr, list) and not isinstance(arr, tuple):
            raise TypeError("incorrect type for row")
        rows, cols = self.get_shape()
        if len(arr) != rows:
            raise ValueError("incorrect rows number for the matrix")
        for i in range(rows):
            if not isinstance(arr[i], type(self.add_idn)):
                raise ValueError("incorrect type in arr for matrix")
            self.arr.append(arr[i])

    def pop_row(self, index: int = -1) -> list[Any]:
        if not isinstance(index, int):
            raise TypeError(f"incorrect type for index: {type(index).__name__}")
        rows, cols = self.get_shape()
        if index >= rows or index < -rows:
            raise IndexError(f"index out of range: {index}")
        return self.arr.pop(index)
    
    def pop_col(self, index: int = -1) -> list[Any]:
        if not isinstance(index, int):
            raise TypeError(f"incorrect type for index: {type(index).__name__}")
        rows, cols = self.get_shape()
        if index >= cols or index < -cols:
            raise IndexError(f"index out of range: {index}")
        result = []
        for row in self.arr:
            result.append(row.pop(index))
        return result

    def get_shape(self) -> tuple[int, int]:
        if len(self.arr) == 0:
            return (0, 0)
        return (len(self.arr), len(self.arr[0]))

    def minor(self, i: int, j: int) -> Matrix[T]:
        if not isinstance(i, int) or not isinstance(j, int):
            raise TypeError("incorrect type for minor parameters")
        rows, cols = self.get_shape()
        if not all((0 <= i < rows, 0 <= j < cols)):
            raise IndexError("out of range")
        if rows == 1 and cols == 1:
            raise ValueError("cannot get minor from a 1x1 matrix")
        arr = self.arr[:i] + self.arr[i+1:]
        for k, a in enumerate(arr):
            arr[k] = a[:j] + a[j+1:]
        return Matrix(arr, self.field)

    def sub_matrix(self, r1: int, c1: int, r2: int | None = None, c2: int | None = None) -> Matrix[T]:
        rows, cols = self.get_shape()
        if r2 is None:
            r2 = rows
        if c2 is None:
            c2 = cols
        if not(isinstance(r1, int) or isinstance(r2, int) or isinstance(c1, int) or isinstance(c2, int)):
            raise TypeError("incorrect type for sub_matrix parameters")
        rows, cols = self.get_shape()
        if any((r1 > rows, r2 > rows, c1 > cols, c2 > cols, 
               r1 < -rows-1, r2 < -rows-1, c1 < -cols-1, c2 < -cols-1)):
            raise IndexError("out of range")
        m: list[list[Any]] = [[None] * abs((c2 - c1)) for _ in range(abs(r2 - r1))]
        for i, r in enumerate(range(r1, r2)):
            for j, c in enumerate(range(c1, c2)):
                m[i][j] = self.field(self.arr[r][c])
        return Matrix(m, self.field)

    def substitude(self, r: int, c: int, sub_matrix: Matrix[T]) -> Matrix[T]:
        rows, cols = self.get_shape()
        sub_rows, sub_cols = sub_matrix.get_shape()
        for idx_i, i in enumerate(range(r, min(rows, r + sub_rows))):
            for idx_j, j in enumerate(range(c, min(cols, c + sub_cols))):
                self.arr[i][j] = sub_matrix[idx_i][idx_j]
        return self

    def change_field(self, f: Callable) -> None:
        for i, x in enumerate(self.arr):
            for j, y in enumerate(x):
                self.arr[i][j] = f(y)
        self.add_idn = f(self.add_idn)
        self.mul_idn = f(self.mul_idn)

    def row_operation_1(self, r1: int, r2: int) -> None:
        """exchange two rows"""
        if not isinstance(r1, int):
            raise TypeError(f"incorrect type: {type(r1).__name__}")
        if not isinstance(r2, int):
            raise TypeError(f"incorrect type: {type(r2).__name__}")
        if r1 == r2:
            return 
        rows, cols = self.get_shape()
        for c in range(cols):
            temp = self.arr[r1][c]
            self.arr[r1][c] = self.arr[r2][c]
            self.arr[r2][c] = temp

    def row_operation_2(self, r: int, mul: Any, is_ignore_add_idn_error: bool = False) -> None:
        """multiply a row by a nonzero scalar"""
        if not isinstance(r, int):
            raise TypeError(f"incorrect type: {type(r).__name__}")
        if mul == self.add_idn and not is_ignore_add_idn_error:
            raise ValueError(f"mul should not be zero")
        rows, cols = self.get_shape()
        for c in range(cols):
            self.arr[r][c] *= mul

    def row_operation_3(self, mul: Any, r1: int, r2: int) -> None:
        """multiply r1 by `mul` and add it to r2 """
        if not isinstance(r1, int):
            raise TypeError(f"incorrect type: {type(r1).__name__}")
        if not isinstance(r2, int):
            raise TypeError(f"incorrect type: {type(r2).__name__}")
        if mul == self.add_idn:
            raise ValueError(f"mul should not be zero")
        rows, cols = self.get_shape()
        for c in range(cols):
            self.arr[r2][c] += self.arr[r1][c] * mul

    def col_operation_1(self, c1: int, c2: int) -> None:
        """exchange two columns"""
        if not isinstance(c1, int):
            raise TypeError(f"incorrect type: {type(c1).__name__}")
        if not isinstance(c1, int):
            raise TypeError(f"incorrect type: {type(c2).__name__}")
        if c1 == c2:
            return 
        rows, cols = self.get_shape()
        for r in range(rows):
            temp = self.arr[r][c1]
            self.arr[r][c1] = self.arr[r][c2]
            self.arr[r][c2] = temp

    def col_operation_2(self, c: int, mul: Any, is_ignore_add_idn_error: bool = False) -> None:
        """multiply a column by a nonzero scalar"""
        if not isinstance(c, int):
            raise TypeError(f"incorrect type: {type(c).__name__}")
        if mul == self.add_idn and not is_ignore_add_idn_error:
            raise ValueError(f"mul should not be zero")
        rows, cols = self.get_shape()
        for r in range(rows):
            self.arr[r][c] *= mul

    def col_operation_3(self, mul: Any, c1: int, c2: int) -> None:
        """multiply c1 by `mul` and add it to c2 """
        if not isinstance(c1, int):
            raise TypeError(f"incorrect type: {type(c1).__name__}")
        if not isinstance(c2, int):
            raise TypeError(f"incorrect type: {type(c2).__name__}")
        if mul == self.add_idn:
            raise ValueError(f"mul should not be zero")
        rows, cols = self.get_shape()
        for r in range(rows):
            self.arr[r][c2] += self.arr[r][c1] * mul

    def col_extend(self, matrix2: Matrix[T]) -> None:
        if not isinstance(matrix2, Matrix):
            raise TypeError(f"incorrect type: {type(matrix2).__name__}")
        if self.get_shape()[0] != matrix2.get_shape()[0]:
            raise ValueError("the number of rows of matrix2 is not corresponding")
        for i, row in enumerate(self.arr):
            row += matrix2[i]
    
    def row_extend(self, matrix2: Matrix[T]) -> None:
        if not isinstance(matrix2, Matrix):
            raise TypeError(f"incorrect type: {type(matrix2).__name__}")
        if self.get_shape()[1] != matrix2.get_shape()[1]:
            raise ValueError("the number of columns of matrix2 is not corresponding")
        self.arr += matrix2.arr

    def is_square(self) -> bool:
        """empty matrix is not square"""
        shape = self.get_shape()
        return shape[0] == shape[1] != 0

    def copy(self) -> Matrix[T]:
        return Matrix(self.arr, self.field)

    def rref(self) -> Matrix[T]:
        """reduced row echolen form"""
        A = self.copy()
        rows, cols = A.get_shape()
        leading = 0
        for c in range(cols):
            for r in range(leading, rows):
                if A.arr[r][c] != self.add_idn:
                    mul = self.mul_idn / A.arr[r][c]
                    # make that row lead by 1
                    A.row_operation_2(r, mul)
                    # exchange it to leading row
                    A.row_operation_1(r, leading)
                    # make that col be 0 except itself
                    for i in range(rows):
                        if i == leading:
                            continue
                        if A.arr[i][c] != self.add_idn:
                            mul = (self.add_idn - self.mul_idn) * A.arr[i][c]
                            A.row_operation_3(mul, leading, i)
                    leading += 1
                    break
        return A

    def nullity(self) -> int:
        temp = 0
        for row in reversed(self.rref().arr):
            if any(c != self.add_idn for c in row):
                break
            else:
                temp += 1
        return temp

    def rank(self) -> int:
        temp = 0
        for row in self.rref():
            if self.mul_idn in row:
                temp += 1
            else:
                break
        return temp

    @property
    def T(self) -> Matrix[T]:
        rows, cols = self.get_shape()
        arr: list[list[T | None]] = [[None] * rows for _ in range(cols)]
        for i, x in enumerate(self.arr):
            for j, y in enumerate(x):
                arr[j][i] = y
        return Matrix(arr, self.field)

    def trace(self) -> Any:
        result = self.add_idn
        for i in range(max(self.get_shape())):
            result += self.arr[i][i]
        return result

    def is_invertible(self) -> bool:
        if not self.is_square():
            raise ValueError("It must be a square matrix.")
        if self.mul_idn in self.rref()[-1]:
            return True
        else:
            return False

    def det(self) -> Any:
        if not self.is_square():
            raise ValueError("It must be a square matrix.")
        
        result = self.mul_idn
        negative_mul_idn = self.add_idn - self.mul_idn
        A = self.copy()
        rows, cols = A.get_shape()
        leading = 0
        for c in range(cols):
            for r in range(leading, rows):
                if A.arr[r][c] != self.add_idn:
                    # change it to the leading row
                    if r != leading:
                        A.row_operation_1(r, leading)
                        result *= negative_mul_idn
                    leading_scalar = A.arr[leading][c]
                    result *= leading_scalar
                    # make the col entries be 0 except itself
                    for r in range(leading, rows):
                        if r == leading:
                            continue
                        if A.arr[r][c] != self.add_idn:
                            mul = negative_mul_idn * A.arr[r][c] / leading_scalar
                            A.row_operation_3(mul, leading, r)
                    leading += 1
                    break
            else:
                return self.add_idn
        return result

    def inverse(self) -> Matrix[T]:
        if not self.is_square():
            raise ValueError("It must be a square matrix.")
        A = self.copy()
        rows, cols = A.get_shape()
        I = Matrix([[0]*i + [1] + [0]*(cols - i - 1) for i in range(rows)], self.field)
        leading = 0
        for c in range(cols):
            for r in range(leading, rows):
                if A.arr[r][c] != self.add_idn:
                    mul = self.mul_idn / A.arr[r][c]
                    # make that row lead by 1
                    A.row_operation_2(r, mul)
                    I.row_operation_2(r, mul)
                    # exchange it to leading row
                    A.row_operation_1(r, leading)
                    I.row_operation_1(r, leading)
                    # make that col be 0 except itself
                    for i in range(rows):
                        if i == leading:
                            continue
                        if A.arr[i][c] != self.add_idn:
                            mul = (self.add_idn - self.mul_idn) * A.arr[i][c]
                            A.row_operation_3(mul, leading, i)
                            I.row_operation_3(mul, leading, i)
                    leading += 1
                    break
            else:
                raise ValueError("This matrix is not invertible")
        return I

    def adjoint(self) -> Matrix[T]:
        rows, cols = self.get_shape()
        if rows != cols:
            raise ValueError("This is not a square matrix")
        if rows == 1 and cols == 1:
            return Matrix([[self.mul_idn]], self.field)
        arr: list[list[Any]] = [[None] * rows for _ in range(cols)]
        negative_mul_idn = self.add_idn - self.mul_idn
        for i in range(rows):
            for j in range(cols):
                arr[j][i] = negative_mul_idn ** (i+j) * self.minor(i, j).det()
        return Matrix(arr, self.field)

    def is_symmetric(self) -> bool:
        rows, cols = self.get_shape()
        if rows != cols:
            return False
        for i in range(rows):
            for j in range(i + 1, cols):
                if self[i][j] != self[j][i]:
                    return False
        return True

    def bilinear_diag(self) -> tuple[Matrix[T], Matrix[T]] | tuple[None, None]:
        """
        Diagonalizes a symmetric bilinear form (matrix) using congruence transformations.
        Returns a tuple (D, E) where D is the diagonalized matrix and E is the transformation matrix,
        or (None, None) if diagonalization is not possible.
        """
        if not self.is_symmetric():
            raise ValueError("It must be a symmetric matrix.")
        rows, cols = self.get_shape()
        A = self.copy()
        E0 = A.identity()
        def reduce(M: Matrix) -> Matrix | None:
            _, n = M.get_shape()
            layer = cols - n
            if n == 1:
                return M
            if M[0][0] == self.add_idn:
                for j in range(1, n):
                    if M[0][j] == 0:
                        continue
                    if M[j][j] == 0:
                        E0.col_operation_3(self.mul_idn, layer + j, layer)
                        M.col_operation_3(self.mul_idn, j, 0)
                        M.row_operation_3(self.mul_idn, j, 0)
                    else:
                        E0.col_operation_1(layer, layer + j)
                        M.col_operation_1(0, j)
                        M.row_operation_1(0, j)
                    break
                else:
                    sub_matrix = reduce(M.sub_matrix(1, 1))
                    if sub_matrix is None:
                        return None
                    return M.substitude(1, 1, sub_matrix)

            for j in range(1, n):
                if M[0][j] == 0:
                    continue
                mul = -(M[0][j] / M[0][0])
                E0.col_operation_3(mul, layer, layer + j)
                M.col_operation_3(mul, 0, j)
                M.row_operation_3(mul, 0, j)
            sub_matrix = reduce(M.sub_matrix(1, 1))
            if sub_matrix is None:
                return None
            return M.substitude(1, 1, sub_matrix)
            
        result_A = reduce(A)
        if result_A is None:
            return None, None
        return result_A, E0
