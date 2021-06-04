<?php

declare(strict_types=1);

namespace Np;

use Np\core\nd;
use Np\core\blas;
use Np\core\lapack;
use Np\reductions\ref;
use Np\reductions\rref;
use Np\decompositions\lu;
use Np\decompositions\svd;
use Np\decompositions\eigen;
use Np\decompositions\cholesky;

/**
 * Matrix
 * A fast lite memory efficient Scientific Computing for php
 * 
 * @package   NumPhp
 * @category  Scientific Computing
 * @author    ghost (Shubham Chaudhary)
 * @email     ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */
class matrix extends nd{

    public $row, $col;
    /**
     * create empty 2d matrix for given data type
     * @param int $row  num of rows 
     * @param int $col  num of cols
     * @param int $dtype matrix data type float|double
     * @return \Np\matrix
     */
    public static function factory(int $row, int $col, int $dtype = self::FLOAT): matrix {
        return new self($row, $col, $dtype);
    }

    /**
     * create 2d matrix using php array
     * @param array $data
     * @param int $dtype matrix data type float|double
     * @return \Np\matrix
     */
    public static function ar(array $data, int $dtype = self::FLOAT): matrix {
        if (is_array($data) && is_array($data[0])) {
            $ar = self::factory(count($data), count($data[0]), $dtype);
            $ar->setData($data);
            unset($data);
            return $ar;
        } else {
            self::_err('given array is not rank-2 or given is not an array');
        }
    }

    /**
     * create one like 2d matrix
     * @param int $row
     * @param int $col
     * @return \Np\matrix
     */
    public static function ones(int $row, int $col, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        for ($i = 0; $i < $ar->ndim; ++$i) {
            $ar->data[$i] = 1;
        }
        return $ar;
    }

    /**
     * Create Matrix with random values
     * @param int $row
     * @param int $col
     * @param int $dtype  Float|Double
     * @return \Np\matrix
     */
    public static function randn(int $row, int $col, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        $max = getrandmax();
        for ($i = 0; $i < $ar->ndim; ++$i) {
            $ar->data[$i] = rand() / $max;
        }
        return $ar;
    }

    /**
     * Return 2d matrix with uniform values
     * @param int $row
     * @param int $col
     * @param int $dtype
     * @return \Np\matrix
     */
    public static function uniform(int $row, int $col, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        $max = getrandmax();
        for ($i = 0; $i < $ar->ndim; ++$i) {
            $ar->data[$i] = rand(-$max, $max) / $max;
        }
        return $ar;
    }

    /**
     * Return a zero matrix with the given dimensions.
     * @param int $row
     * @param int $col
     * @param int $dtype
     * @return \Np\matrix
     */
    public static function zeros(int $row, int $col, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        for ($i = 0; $i < $ar->ndim; ++$i) {
            $ar->data[$i] = 0.0;
        }
        return $ar;
    }

    /**
     * create a null like 2d matrix
     * @param int $row
     * @param int $col
     * @return \Np\matrix
     */
    public static function null(int $row, int $col, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        for ($i = 0; $i < $ar->ndim; ++$i) {
            $ar->data[$i] = null;
        }
        return $ar;
    }

    /**
     * create a 2d matrix with given scalar value
     * @param int $row
     * @param int $col
     * @param int|float|double $val
     * @return \Np\matrix
     */
    public static function full(int $row, int $col, $val, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        for ($i = 0; $i < $ar->ndim; ++$i) {
            $ar->data[$i] = $val;
        }
        return $ar;
    }

    /**
     * create a diagonal 2d matrix with given 1d array;
     * @param array $elements
     * @return \Np\matrix
     */
    public static function diagonal(array $elements, int $dtype = self::FLOAT): matrix {
        $n = count($elements);
        $ar = self::factory($n, $n, $dtype);
        for ($i = 0; $i < $n; ++$i) {
            $ar->data[$i * $n + $i] = $elements[$i]; #for ($j = 0; $j < $n; ++$j) {$i === $j ? $elements[$i] : 0;#} 
        }
        return $ar;
    }

    /**
     * Generate a m x n matrix with elements from a Poisson distribution.
     *
     * @param int $row
     * @param int $col
     * @param float $lambda
     * @param int $dtype 
     * @return \Np\matrix
     */
    public static function poisson(int $row, int $col, float $lambda = 1.0, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        $max = getrandmax();
        $l = exp(-$lambda);
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $k = 0;
                $p = 1.0;
                while ($p > $l) {
                    ++$k;
                    $p = $p * rand() / $max;
                }
                $ar->data[$i * $col + $j] = $k - 1;
            }
        }
        return $ar;
    }

    /**
     * Return a standard normally distributed random matrix i.e values
     * between -1 and 1.
     * @param int $row
     * @param int $col
     * @param int $dtype Description
     * @return \Np\matrix
     */
    public static function gaussian(int $row, int $col, int $dtype = self::FLOAT): matrix {
        $max = getrandmax();
        $a = $extras = [];

        while (count($a) < $row) {
            $rowA = [];

            if (!empty($extras)) {
                $rowA[] = array_pop($extras);
            }

            while (count($rowA) < $col) {
                $r = sqrt(-2.0 * log(rand() / $max));

                $phi = rand() / $max * self::TWO_PI;

                $rowA[] = $r * sin($phi);
                $rowA[] = $r * cos($phi);
            }

            if (count($rowA) > $col) {
                $extras[] = array_pop($rowA);
            }

            $a[] = $rowA;
        }

        return self::ar($a, $dtype);
    }

    /**
     * create an identity matrix with the given dimensions.
     * @param int $n
     * @param int $dtype
     * @return matrix
     * @throws \InvalidArgumentException
     */
    public static function identity(int $n, int $dtype = self::FLOAT): matrix {
        if ($n < 1) {
            self::_invalidArgument('Dimensionality must be greater than 0 on all axes.');
        }

        $ar = self::factory($n, $n, $dtype);
        for ($i = 0; $i < $n; ++$i) {
            for ($j = 0; $j < $n; ++$j) {
                $ar->data[$i * $n + $j] = $i === $j ? 1 : 0;
            }
        }
        return $ar;
    }

    /**
     * Return the element-wise minimum of two matrices.
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    public function minimum(matrix $m): matrix {
        if ($this->checkShape($m)) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->ndim; ++$i) {
                $ar->data[$i] = min($this->data[$i], $m->data[$i]);
            }
            return $ar;
        }
    }

    /**
     * Return the element-wise maximum of two matrices.
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    public function maximum(matrix $m): matrix {
        if ($this->checkShape($m)) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->ndim; ++$i) {
                $ar->data[$i] = max($this->data[$i], $m->data[$i]);
            }
            return $ar;
        }
    }

    /**
     * 2D convolution between a matrix ma and kernel kb, with a given stride.
     * @param \Np\matrix $m
     * @param int $stride
     * @return matrix
     */
    public function convolve(matrix $m, int $stride = 1): matrix {
        return convolve::conv2D($this, $m, $stride);
    }

    /**
     * Calculate the determinant of the matrix.
     * @return float
     */
    public function det(): float {
        if (!$this->isSquare()) {
            self::_err('determinant is undefined for a non square matrix');
        }
        $lu = $this->lu();
        $nSwaps = $lu->p()->diagonalAsVector()->subtract($lu->p()->diagonalAsVector()->sum())->col - 1;
        $detP = (-1) ** $nSwaps;
        $detL = $lu->l()->diagonalAsVector()->product();
        $detU = $lu->u()->diagonalAsVector()->product();
        unset($lu);
        return ($detP * $detL * $detU);
    }

    /**
     * Return the trace of the matrix i.e the sum of all diagonal elements of a square matrix.
     * @return float
     */
    public function trace(): float {
        if (!$this->isSquare()) {
            self::_err('Error::matrix is not a squared can not Trace!');
        }
        $trace = 0.0;
        for ($i = 0; $i < $this->row; ++$i) {
            for ($j = 0; $j < $this->col; ++$j) {
                if ($i == $j) {
                    $trace += $this->data[$i * $this->col + $i];
                }
            }
        }
        return $trace;
    }

    /**
     * dignoalInterChange
     */
    public function dignoalInterChange() {
        for ($i = 0; $i < $this->row; ++$i) {
            for ($j = 0; $j < $this->col; ++$j) {
                $tmp = $this->data[$i * $this->col - $j];
                $this->data[$i * $this->col - $j] = $tmp;
            }
        }
    }

    //----------------Linear Algebra Opreations-------------------------------

    /**
     *  
     * get dot product of m.m or m.v
     * 
     * @param \Np\matrix|\Np\vector $d
     * @return matrix|vector
     */
    public function dot(matrix|vector $d): matrix|vector {
        if ($d instanceof self) {
            return $this->dotMatrix($d);
        } else {
            return $this->dotVector($d);
        }
    }

    /**
     * get matrix & matrix dot product
     * @param \Np\matrix $matrix
     * @return \Np\matrix
     */
    protected function dotMatrix(matrix $matrix): matrix {
        if ($this->checkDtype($matrix) && $this->checkDimensions($matrix)) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            blas::gemm($this, $matrix, $ar);
            return $ar;
        }
    }

    /**
     * get dot product of matrix & a vector
     * @param \Np\vector $vector
     * @return \Np\vector
     */
    protected function dotVector(vector $vector): vector {
        if ($this->dtype != $vector->dtype) {
            self::_dtypeErr('Mismatch Dtype of given vector');
        }
        if ($this->row != $vector->col) {
            self::_dimensionaMisMatchErr('Mismatch row and col of matrix and vector');
        }
        $mvr = vector::factory($this->col, $this->dtype);
        blas::gemv($this, $vector, $mvr);
        return $mvr;
    }

    //---------------Arthmetic Opreations-----------------------------------

    /**
     * multiply this matrix with another matrix|scalar element-wise
     * Matrix Scalar\Matrix multiplication
     * @param int|float|matrix|vector $m
     * @return matrix|vector
     */
    public function multiply(int|float|matrix|vector $m): matrix|vector {
        if ($m instanceof self) {
            return $this->multiplyMatrix($m);
        } else if ($m instanceof vector) {
            return $this->multiplyVector($m);
        } else {
            return $this->scale($m);
        }
    }

    /**
     * 
     * @param \Np\vector $v
     * @return matrix
     */
    protected function multiplyVector(vector $v): matrix {
        if ($this->row == $v->col && $this->dtype == $v->dtype) {
            $ar = matrix::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] = $v->data[$j] * $this->data[$i * $this->col + $j];
                }
            }
            return $ar;
        }
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function multiplyMatrix(matrix $m): matrix {
        if ($this->checkDtype($m) && $this->checkShape($m)) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] = $this->data[$i * $this->col + $j] * $m->data[$i * $this->col + $j];
                }
            }
            return $ar;
        }
    }

    /**
     * 
     * @param int|float $scalar
     * @return matrix
     */
    public function scale(int|float $scalar): matrix {
        if ($scalar == 0) {
            return self::zeros($this->row, $this->col, $this->dtype);
        }

        $ar = $this->copyMatrix();
        for ($i = 0; $i < $this->ndim; ++$i) {
            $ar->data[$i] *= $scalar;
        }

        return $ar;
    }

    /**
     * Sum of Rows of matrix
     * @return vector
     */
    public function sumRows(): vector {
        $vr = vector::factory($this->row, $this->dtype);
        for ($i = 0; $i < $this->row; ++$i) {
            $sum = 0.0;
            for ($j = 0; $j < $this->col; ++$j) {
                $sum += $this->data[$i * $this->col + $j];
            }
            $vr->data[$i] = $sum;
        }
        return $vr;
    }

    /**
     * Sum of two matrix, vector or a scalar to current matrix
     * 
     * @param int|float|matrix|vector $m
     * @return matrix
     */
    public function sum(int|float|matrix|vector $m): matrix {
        if ($m instanceof self) {
            return $this->sumMatrix($m);
        } elseif ($m instanceof vector) {
            return $this->sumVector($m);
        } else {
            return $this->sumScalar($m);
        }
    }

    protected function sumScalar(int|float $s): matrix {
        $ar = self::factory($this->row, $this->col, $this->dtype);
        for ($i = 0; $i < $this->ndim; ++$i) {
            $ar->data[$i] = $this->data[$i] + $s;
        }
        return $ar;
    }

    protected function sumMatrix(matrix $m): matrix {
        if ($this->checkShape($m) && $this->checkDtype($m)) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->ndim; ++$i) {
                $ar->data[$i] = $this->data[$i] + $m->data[$i];
            }
            return $ar;
        }
    }

    protected function sumVector(vector $v): matrix {
        if ($this->row == $v->col && $this->dtype == $v->dtype) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] = $v->data[$j] + $this->data[$i * $this->col + $j];
                }
            }
            return $ar;
        }
    }

    /**
     * subtract another matrix, vector or a scalar to this matrix
     * @param int|float|matrix|vector $d matrix|$scalar to subtract this matrix
     * @return \Np\matrix
     */
    public function subtract(int|float|matrix|vector $d): matrix {
        if ($d instanceof self) {
            return $this->subtractMatrix($d);
        } elseif ($d instanceof vector) {
            return $this->subtractVector($d);
        } else {
            return $this->subtractScalar($d);
        }
    }

    protected function subtractScalar(int|float $s): matrix {
        $ar = self::factory($this->row, $this->col, $this->dtype);
        for ($i = 0; $i < $this->ndim; ++$i) {
            $ar->data[$i] = $this->data[$i] - $s;
        }
        return $ar;
    }

    /**
     * 
     * @param matrix $m
     * @return matrix
     */
    protected function subtractMatrix(matrix $m): matrix {
        if ($this->checkShape($m) && $this->checkDtype($m)) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->ndim; ++$i) {
                $ar->data[$i] = $this->data[$i] - $m->data[$i];
            }
            return $ar;
        }
    }

    /**
     * 
     * @param vector $v
     * @return matrix
     */
    protected function subtractVector(vector $v): matrix {
        if ($this->row == $v->col && $this->dtype == $v->dtype) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] = $this->data[$i * $this->col + $j] - $v->data[$j];
                }
            }
            return $ar;
        }
    }

    /**
     * 
     * @param vector $v
     * @return matrix
     */
    public function subtractColumnVector(vector $v): matrix {
        if ($this->row == $v->col && $this->dtype == $v->dtype) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($j = 0; $j < $this->col; ++$j) {
                for ($i = 0; $i < $this->row; ++$i) {
                    $ar->data[$i * $this->col + $j] = $this->data[$i * $this->col + $j] - $v->data[$i];
                }
            }
            return $ar;
        }
    }

    /**
     * Return the division of two elements, element-wise.
     * @param int|float|matrix $d
     * @return matrix
     */
    public function divide(int|float|matrix|vector $d): matrix {
        if ($d instanceof self) {
            return $this->divideMatrix($d);
        } elseif ($d instanceof vector) {
            return $this->divideVector($d);
        } else {
            return $this->divideScalar($d);
        }
    }

    protected function divideMatrix(matrix $m): matrix {
        if ($this->checkShape($m) && $this->checkDtype($m)) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->ndim; ++$i) {
                $ar->data[$i] = $this->data[$i] / $m->data[$i];
            }
            return $ar;
        }
    }

    protected function divideVector(vector $v): matrix {
        if ($this->row == $v->col && $this->dtype == $v->dtype) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] = $this->data[$i * $this->col + $j] / $v->data[$j];
                }
            }
            return $ar;
        }
    }

    protected function divideScalar(int|float $s): matrix {
        $ar = self::factory($this->row, $this->col, $this->dtype);
        for ($i = 0; $i < $this->ndim; ++$i) {
            $ar->data[$i] = $this->data[$i] / $s;
        }
        return $ar;
    }

    /**
     * 
     * Raise this matrix to the power of the element-wise entry in another matrix.
     * 
     * @param int|float|matrix $m
     * @return matrix
     */
    public function pow(int|float|matrix|vector $d): matrix {
        if ($d instanceof self) {
            return $this->powMatrix($d);
        } else if ($d instanceof vector) {
            return $this->powVector($d);
        } else {
            return $this->powScalar($d);
        }
    }

    protected function powMatrix(matrix $m): matrix {
        if ($this->checkShape($m) && $this->checkDtype($m)) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->ndim; ++$i) {
                $ar->data[$i] = $this->data[$i] ** $m->data[$i];
            }
            return $ar;
        }
    }

    protected function powVector(vector $v): matrix {
        if ($this->row == $v->col && $this->dtype == $v->dtype) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] = $this->data[$i * $this->col + $j] ** $v->data[$j];
                }
            }
            return $ar;
        }
    }

    protected function powScalar(int|float $s): matrix {
        $ar = $this->copyMatrix();
        for ($i = 0; $i < $this->ndim; ++$i) {
            $ar->data[$i] **= $s;
        }
        return $ar;
    }

    /**
     * Calculate the modulus i.e remainder of division between this matrix and another matrix.
     * @param int|float|matrix|vector $d
     * @return matrix
     */
    public function mod(int|float|matrix|vector $d): matrix {
        if ($d instanceof self) {
            $this->modMatrix($d);
        } else if ($d instanceof vector) {
            $this->modVector($d);
        } else {
            $this->modScalar($d);
        }
    }

    protected function modMatrix(matrix $m): matrix {
        if ($this->checkShape($m) && $this->checkDtype($m)) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->ndim; ++$i) {
                $ar->data[$i] = $this->data[$i] % $m->data[$i];
            }
            return $ar;
        } 
    }

    protected function modVector(vector $v): matrix {
        if ($this->row == $v->col && $this->dtype == $v->dtype) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] = $this->data[$i * $this->col + $j] % $v->data[$j];
                }
            }
            return $ar;
        }
    }

    protected function modScalar(int|float $s): matrix {
        $ar = $this->copyMatrix();
        for ($i = 0; $i < $this->ndim; ++$i) {
            $ar->data[$i] %= $s;
        }
        return $ar;
    }

    /**
     * Return the element-wise reciprocal of the matrix.
     * 
     * @return matrix
     */
    public function reciprocal(): matrix {
        return self::ones($this->row, $this->col, $this->dtype)->divideMatrix($this);
    }

    /**
     * 
     * @param int|float $d
     * @return bool
     */
    public static function is_zero($d): bool {
        if (abs($d) < self::EPSILON) {
            return true;
        }
        return false;
    }

    /**
     *  is row zero
     * @param int $row
     * @return bool
     */
    public function is_rowZero(int $row): bool {
        for ($i = 0; $i < $this->col; ++$i) {
            if ($this->data[$row * $this->col + $i] != 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * 
     * @return bool
     */
    public function has_ZeroRow(): bool {
        for ($i = 0; $i < $this->row; ++$i) {
            if ($this->is_rowZero($i)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Transpose the matrix i.e row become cols and cols become rows.
     * @return \Np\matrix
     */
    public function transpose(): matrix {
        $ar = self::factory($this->col, $this->row, $this->dtype);
        for ($i = 0; $i < $ar->row; ++$i) {
            for ($j = 0; $j < $ar->col; ++$j) {
                $ar->data[$i * $ar->col + $j] = $this->data[$j * $ar->col + $i];
            }
        }
        return $ar;
    }

    /**
     * swap specific values in matrix
     * @param int $i1
     * @param int $i2
     */
    public function swapValue(int $i1, int $i2) {
        $tmp = $this->data[$i1];
        $this->data[$i1] = $this->data[$i2];
        $this->data[$i2] = $tmp;
    }

    /**
     * swap specific rows in matrix
     * @param int $r1
     * @param int $r2
     */
    public function swapRows(int $r1, int $r2) {
        for ($i = 0; $i < $this->col; ++$i) {
            $tmp = $this->data[$r1 * $this->col + $i];
            $this->data[$r1 * $this->col + $i] = $this->data[$r2 * $this->col + $i];
            $this->data[$r2 * $this->col + $i] = $tmp;
        }
    }

    /**
     * swap specific cols in matrix
     * @param int $c1
     * @param int $c2
     */
    public function swapCols(int $c1, int $c2) {
        for ($i = 0; $i < $this->row; ++$i) {
            $tmp = $this->data[$i * $this->row + $c1];
            $this->data[$i * $this->row + $c1] = $this->data[$i * $this->row + $c2];
            $this->data[$i * $this->row + $c2] = $tmp;
        }
    }

    /**
     * scale all the elements of a row 
     * @param int $row
     * @param float $c
     */
    public function scaleRow(int $row, float $c) {
        for ($i = 0; $i < $this->col; ++$i) {
            $this->data[$row * $this->col + $i] *= $c;
        }
    }

    /**
     * 
     * @param int $r1
     * @param int $r2
     * @param float $c
     */
    public function addScaleRow(int $r1, int $r2, float $c) {
        for ($i = 0; $i < $this->col; ++$i) {
            $this->data[$r2 * $this->col + $i] += $this->data[$r1 * $this->col + $i] * $c;
        }
    }

    /**
     * Attach given matrix to the left of this matrix.
     * 
     * @param \Np\matrix $m
     * @return \Np\matrix
     */
    public function joinLeft(matrix $m): matrix {
        if ($this->row != $m->row && $this->checkDtype($m)) {
            self::_err('Error::Invalid size! or DataType!');
        }
        $col = $this->col + $m->col;
        $ar = self::factory($this->row, $col, $this->dtype);
        for ($i = 0; $i < $this->row; ++$i) {
            for ($j = 0; $j < $this->col; ++$j) {
                $ar->data[$i * $col + $j] = $this->data[$i * $this->col + $j];
            }
            for ($j = 0; $j < $m->col; ++$j) {
                $ar->data[$i * $col + ($this->col + $j)] = $m->data[$i * $m->col + $j];
            }
        }
        return $ar;
    }

    /**
     * Join matrix m to the Right of this matrix.
     * @param \Np\matrix $m
     * @return matrix
     */
    public function joinRight(matrix $m): matrix {
        if ($this->row != $m->row && $this->checkDtype($m)) {
            self::_err('Error::Invalid size! or DataType!');
        }
        $col = $this->col + $m->col;
        $ar = self::factory($this->row, $col, $this->dtype);
        for ($i = 0; $i < $m->row; ++$i) {
            for ($j = 0; $j < $m->col; ++$j) {
                $ar->data[$i * $col + $j] = $m->data[$i * $m->col + $j];
            }
            for ($j = 0; $j < $this->col; ++$j) {
                $ar->data[$i * $col + ($this->col + $j)] = $this->data[$i * $this->col + $j];
            }
        }
        return $ar;
    }

    /**
     * Join matrix m Above this matrix.
     * @param \Np\matrix $m
     * @return matrix
     */
    public function joinAbove(matrix $m): matrix {
        if ($this->col !== $m->col && $this->checkDtype($m)) {
            self::_err('Error::Invalid size! or DataType!');
        }
        $row = $this->row + $m->row;
        $ar = self::factory($row, $this->col, $this->dtype);
        for ($i = 0; $i < $m->row; ++$i) {
            for ($j = 0; $j < $m->col; ++$j) {
                $ar->data[$i * $m->col + $j] = $m->data[$i * $m->col + $j];
            }
            for ($j = 0; $j < $this->col; ++$j) {
                $ar->data[($i + $this->row) * $this->col + $j] = $this->data[$i * $this->col + $j];
            }
        }
        return $ar;
    }

    /**
     * Join matrix m below this matrix.
     * @param \Np\matrix $m
     * @return matrix
     */
    public function joinBelow(matrix $m): matrix {
        if ($this->col !== $m->col && $this->checkDtype($m)) {
            self::_err('Error::Invalid size! or DataType!');
        }
        $row = $this->row + $m->row;
        $ar = self::factory($row, $this->col, $this->dtype);
        for ($i = 0; $i < $this->row; ++$i) {
            for ($j = 0; $j < $this->col; ++$j) {
                $ar->data[$i * $this->col + $j] = $this->data[$i * $this->col + $j];
            }
            for ($j = 0; $j < $m->col; ++$j) {
                $ar->data[($i + $m->row) * $m->col + $j] = $m->data[$i * $m->col + $j];
            }
        }
        return $ar;
    }

    /**
     * Calculate the row echelon form of the matrix. 
     * Return the reduced matrix.
     *
     * @return matrix|null
     */
    public function ref(): matrix|null {
        return ref::factory($this);
    }

    /**
     * Return the lower triangular matrix of the Cholesky decomposition.
     *
     * @return matrix|null
     */
    public function cholesky(): matrix|null {
        return cholesky::factory($this);
    }

    /**
     * FIXME--------------
     * RREF
     * The reduced row echelon form (RREF) of a matrix.
     * @return \Np\matrix
     */
    public function rref(): matrix {
        return rref::factory($this);
    }

    /**
     * make copy of the matrix
     * @return \Np\matrix
     */
    public function copyMatrix(): matrix {
        return clone $this;
    }

    /**
     * 
     * @param int $cols
     * @return \Np\matrix
     */
    public function diminish_left(int $cols): matrix {
        $ar = self::factory($this->row, $cols, $this->dtype);
        for ($i = 0; $i < $ar->row; ++$i) {
            for ($j = 0; $j < $ar->col; ++$j) {
                $ar->data[$i * $ar->col + $j] = $this->data[$i * $this->col + $j];
            }
        }
        return $ar;
    }

    /**
     * 
     * @param int $cols
     * @return \Np\matrix
     */
    public function diminish_right(int $cols): matrix {
        $ar = self::factory($this->row, $cols, $this->dtype);
        for ($i = 0; $i < $ar->row; ++$i) {
            for ($j = 0; $j < $ar->col; ++$j) {
                $ar->data[$i * $ar->col + $j] = $this->data[$i * $this->col - $cols + $j];
            }
        }
        return $ar;
    }

    /**
     * Return the index of the maximum element in every row of the matrix.
     * @return \Np\vector int
     */
    public function argMax(): vector {
        $v = vector::factory($this->row, vector::INT);
        for ($i = 0; $i < $this->row; ++$i) {
            $v->data[$i] = blas::max($this->rowAsVector($i));
        }
        return $v;
    }

    /**
     * Return the index of the minimum element in every row of the matrix.
     * @return \Np\vector int
     */
    public function argMin(): vector {
        $v = vector::factory($this->row, vector::INT);
        for ($i = 0; $i < $this->row; ++$i) {
            $v->data[$i] = blas::min($this->rowAsVector($i));
        }

        return $v;
    }

    /**
     * Set given data in matrix
     * @param int|float|array $data
     * @param bool $dignoal
     * @return void
     */
    public function setData(int|float|array $data, bool $dignoal = false): void {
        if ($dignoal == false) {
            if (is_array($data) && is_array($data[0])) {
                $f = $this->flattenArray($data);
                foreach ($f as $k => $v) {
                    $this->data[$k] = $v;
                }
            } elseif (is_numeric($data) && $dignoal == false) {
                for ($i = 0; $i < $this->ndim; ++$i) {
                    $this->data[$i] = $data;
                }
            } elseif (is_numeric($data) && $dignoal == true) {
                for ($i = 0; $i < $this->row; ++$i) {
                    $this->data[$i * $this->col * $i] = $data;
                }
            }
        }
    }

    /**
     * get the shape of matrix
     * @return object
     */
    public function getShape(): object {
        return (object) ['m' => $this->row, 'n' => $this->col];
    }

    /**
     * get the number of elements in the matrix.
     * @return int
     */
    public function getSize(): int {
        return $this->row * $this->col;
    }

    /**
     * 
     * @return bool
     */
    public function isSquare(): bool {
        if ($this->row === $this->col) {
            return true;
        }
        return false;
    }

    public function getDtype() {
        return $this->dtype;
    }

    /**
     * Return a row as vector from the matrix.
     * @param int $index
     * @return \Np\vector
     */
    public function rowAsVector(int $index): vector {
        $vr = vector::factory($this->col, $this->dtype);
        for ($j = 0; $j < $this->col; ++$j) {
            $vr->data[$j] = $this->data[$index * $this->col + $j];
        }
        return $vr;
    }

    /**
     * Return a col as vector from the matrix.
     * @param int $index
     * @return \Np\vector
     */
    public function colAsVector(int $index): vector {
        $vr = vector::factory($this->row, $this->dtype);
        for ($i = 0; $i < $this->row; ++$i) {
            $vr->data[$i] = $this->data[$i * $this->row + $index];
        }
        return $vr;
    }

    /**
     * Return the diagonal elements of a square matrix as a vector.
     * @return \Np\vector
     */
    public function diagonalAsVector(): vector {
        if (!$this->isSquare()) {
            self::_err('Can not trace of a none square matrix');
        }
        $vr = vector::factory($this->row, $this->dtype);
        for ($i = 0; $i < $this->row; ++$i) {
            $vr->data[$i] = $this->getDiagonalVal($i);
        }
        return $vr;
    }

    /**
     * Flatten i.e unravel the matrix into a vector.
     *
     * @return \Np\vector
     */
    public function asVector(): vector {
        $vr = vector::factory($this->ndim, $this->dtype);
        for ($i = 0; $i < $this->ndim; ++$i) {
            $vr->data[$i] = $this->data[$i];
        }
        return $vr;
    }

    /**
     * Return the elements of the matrix in a 2-d array.
     * @return array
     */
    public function asArray(): array {
        $ar = array_fill(0, $this->row, array_fill(0, $this->col, null));
        for ($i = 0; $i < $this->row; ++$i) {
            for ($j = 0; $j < $this->col; ++$j) {
                $ar[$i][$j] = $this->data[$i * $this->col + $j];
            }
        }
        return $ar;
    }

    /**
     * get a diagonal value from matrix
     * @param int $i
     * @return float
     */
    public function getDiagonalVal(int $i) {
        if ($this->isSquare()) {
            return $this->data[$i * $this->row + $i];
        }
    }

    /**
     *
     * Compute the multiplicative inverse of the matrix.
     * @return matrix
     */
    public function inverse(): matrix {
        if (!$this->isSquare()) {
            self::_err('Error::invalid Size of matrix!');
        }
        $imat = $this->copyMatrix();
        $ipiv = vector::factory($this->row, vector::INT);
        $lp = lapack::getrf($imat, $ipiv);
        if ($lp != 0) {
            return null;
        }
        $lp = lapack::getri($imat, $ipiv);
        if ($lp != 0) {
            return null;
        }
        unset($ipiv);
        unset($lp);
        return $imat;
    }
    
    /**
     * Compute the (Moore-Penrose) pseudo inverse of the general matrix.
     * @return matrix|null
     */
    public function pseudoInverse(): matrix|null {
        $k = min($this->row, $this->col);
        $s = vector::factory($k, $this->dtype);
        $u = self::factory($this->row, $this->row, $this->dtype);
        $vt = self::factory($this->col, $this->col, $this->dtype);
        $imat = $this->copyMatrix();
        $lp = lapack::gesdd($imat, $s, $u, $vt);
        if ($lp != 0) {
            return null;
        }
        for ($i = 0; $i < $k; ++$i) {
            blas::scale(1.0 / $s->data[$i], $vt->rowAsVector($i));
        }
        unset($imat);
        unset($k);
        unset($lp);
        unset($s);
        $mr = self::factory($this->col, $this->row, $this->dtype);
        blas::gemm($vt, $u, $mr);
        unset($u);
        unset($vt);
        return $mr;
    }

    /**
     * Compute the singular value decomposition of a matrix and 
     * return an object of the singular values and unitary matrices
     *
     * @return object (u,s,v)
     */
    public function svd(): svd {
        return svd::factory($this);
    }

    /**
     * Compute the eigen decomposition of a general matrix.
     * return the eigenvalues and eigenvectors as object
     * 
     * @param bool $symmetric
     * @return eigen
     */
    public function eign(bool $symmetric = false): eigen {
        return eigen::factory($this, $symmetric);
    }

    /**
     *  
     * Compute the LU factorization of matrix.
     * return lower, upper, and permutation matrices as object.
     * 
     * @return lu
     */
    public function lu(): lu {
        return lu::factory($this);
    }

    /**
     * Return the L1 norm of the matrix.
     * @return float
     */
    public function normL1(): float {
        return lapack::lange('l', $this);
    }

    /**
     * Return the L2 norm of the matrix.
     * @return float
     */
    public function normL2(): float {
        return lapack::lange('f', $this);
    }

    /**
     * Return the L1 norm of the matrix.
     * @return float
     */
    public function normINF(): float {
        return lapack::lange('i', $this);
    }

    /**
     * Return the Frobenius norm of the matrix.
     * @return float
     */
    public function normFrob(): float {
        return $this->normL2();
    }

    /**
     * Run a function over all of the elements in the matrix. 
     * @param callable $func
     * @return \Np\matrix
     */
    public function map(callable $func): matrix {
        $ar = self::factory($this->row, $this->col, $this->dtype);
        for ($i = 0; $i < $this->ndim; ++$i) {
            $ar->data[$i] = $func($this->data[$i]);
        }
        return $ar;
    }

    public function abs(): matrix {
        return $this->map('abs');
    }

    public function sqrt(): matrix {
        return $this->map('sqrt');
    }

    public function exp(): matrix {
        return $this->map('exp');
    }

    public function exp1(): matrix {
        return $this->map('exp1');
    }

    public function log(float $b = M_E): matrix {
        $ar = $this->copyMatrix();
        for ($i = 0; $i < $ar->ndim; ++$i) {
            log($ar->data[$i], $b);
        }
        return $ar;
    }

    public function log1p(): matrix {
        return $this->map('log1p');
    }

    public function sin(): matrix {
        return $this->map('sin');
    }

    public function asin(): matrix {
        return $this->map('asin');
    }

    public function cos(): matrix {
        return $this->map('cos');
    }

    public function acos(): matrix {
        return $this->map('acos');
    }

    public function tan(): matrix {
        return $this->map('tan');
    }

    public function atan(): matrix {
        return $this->map('atan');
    }

    public function radToDeg(): matrix {
        return $this->map('rad2deg');
    }

    public function degToRad(): matrix {
        return $this->map('deg2rad');
    }

    public function floor(): matrix {
        return $this->map('floor');
    }

    public function ceil(): matrix {
        return $this->map('ceil');
    }

    /**
     * Compute the means of each row and return them in a vector.
     *
     * @return vector
     */
    public function mean(): vector {
        return $this->sumRows()->divide($this->col);
    }

    /**
     * Compute the row variance of the matrix.
     * 
     * @param vector|null $mean
     * @return vector
     */
    public function variance(vector|null $mean = null): vector {
        if (isset($mean)) {
            if (!$mean instanceof vector) {
                self::_invalidArgument('mean must be a vector!');
            }
            if ($this->row !== $mean->col) {
                self::_err('Err:: given mean vector dimensionality mismatched!');
            }
        } else {
            $mean = $this->mean();
        }
        return $this->subtractColumnVector($mean)->square()
                        ->sumRows()->divide($this->row);
    }

    /**
     *  Return the median vector of this matrix.
     * @return vector
     */
    public function median(): vector {
        $mid = intdiv($this->col, 2);
        $odd = $this->col % 2 === 1;
        $vr = vector::factory($this->row, $this->dtype);
        for ($i = 0; $i < $this->row; ++$i) {
            $a = $this->rowAsVector($i)->sort();
            if ($odd) {
                $median = $a->data[$mid];
            } else {
                $median = ($a->data[$mid - 1] + $a->data[$mid]) / 2.0;
            }
            $vr->data[$i] = $median;
        }
        unset($a);
        return $vr;
    }

    /**
     * Compute the covariance matrix.
     * 
     * @param vector|null $mean
     * @return matrix
     */
    public function covariance(vector|null $mean = null): matrix {
        if (isset($mean)) {
            if ($mean->col !== $this->row) {
                self::_err('Err:: given mean vector dimensionality mismatched!');
            }
        } else {
            $mean = $this->mean();
        }

        $b = $this->subtractColumnVector($mean);

        return $b->dot($b->transpose())
                        ->divideScalar($this->row);
    }

    /**
     * Clip the elements in the matrix to be between given minimum and maximum
     * and return a new matrix.
     * 
     * @param float $min
     * @param float $max
     * @return matrix
     */
    public function clip(float $min, float $max): matrix {
        $ar = self::factory($this->row, $this->col, $this->dtype);
        for ($i = 0; $i < $this->ndim; ++$i) {
            if ($this->data[$i] > $max) {
                $ar->data[$i] = $max;
                continue;
            }
            if ($this->data[$i] < $min) {
                $ar->data[$i] = $min;
                continue;
            }
            $ar->data[$i] = $this->data[$i];
        }
        return $ar;
    }

    /**
     * Clip the matrix to be lower bounded by a given minimum.
     * @param float $min
     * @return matrix
     */
    public function clipLower(float $min): matrix {
        $ar = self::factory($this->row, $this->col, $this->dtype);
        for ($i = 0; $i < $this->ndim; ++$i) {
            if ($this->data[$i] < $min) {
                $ar->data[$i] = $min;
                continue;
            }
            $ar->data[$i] = $this->data[$i];
        }
        return $ar;
    }

    /**
     * Clip the matrix to be upper bounded by a given maximum.
     *
     * @param float $max
     * @return matrix
     */
    public function clipUpper(float $max): matrix {
        $ar = self::factory($this->row, $this->col, $this->dtype);
        for ($i = 0; $i < $this->ndim; ++$i) {
            if ($this->data[$i] > $max) {
                $ar->data[$i] = $max;
                continue;
            }
            $ar->data[$i] = $this->data[$i];
        }
        return $ar;
    }

    /**
     * Square of matrix
     * @return matrix
     */
    public function square(): matrix {
        return $this->multiplyMatrix($this);
    }

    /**
     * 
     * @param int|float|matrix|vector $d
     * @return matrix
     */
    public function equal(int|float|matrix|vector $d): matrix {
        if ($d instanceof self) {
            return $this->equalMatrix($d);
        } elseif ($d instanceof vector) {
            return $this->equalVector($d);
        } else {
            return $this->equalScalar($d);
        }
    }

    protected function equalMatrix(matrix $m): matrix {
        if ($this->checkShape($m) && $this->checkDtype($m)) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->ndim; ++$i) {
                $ar->data[$i] = $this->data[$i] == $m->data[$i] ? 1 : 0;
            }
            return $ar;
        }
    }

    protected function equalVector(vector $v): matrix {
        if ($this->row == $v->col && $this->dtype == $v->dtype) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] = $this->data[$i * $this->col + $j] == $v->data[$j] ? 1 : 0;
                }
            }
            return $ar;
        }
    }

    protected function equalScalar(int|float $s): matrix {
        $ar = self::factory($this->row, $this->col, $this->dtype);
        for ($i = 0; $i < $this->ndim; ++$i) {
            $ar->data[$i] = $this->data[$i] == $s ? 1 : 0;
        }
        return $ar;
    }

    /**
     * 
     * @param int|float|matrix|vector $d
     * @return matrix
     */
    public function greater(int|float|matrix|vector $d): matrix {
        if ($d instanceof self) {
            return $this->greaterMatrix($d);
        } elseif ($d instanceof vector) {
            return $this->greaterVector($d);
        } else {
            return $this->greaterScalar($d);
        }
    }

    protected function greaterMatrix(matrix $m): matrix {
        if ($this->checkShape($m) && $this->checkDtype($m)) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->ndim; ++$i) {
                $ar->data[$i] = $this->data[$i] > $m->data[$i] ? 1 : 0;
            }
            return $ar;
        }
    }

    protected function greaterVector(vector $v): matrix {
        if ($this->row == $v->col && $this->dtype == $v->dtype) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] = $this->data[$i * $this->col + $j] > $v->data[$j] ? 1 : 0;
                }
            }
            return $ar;
        }
    }

    protected function greaterScalar(int|float $s): matrix {
        $ar = self::factory($this->row, $this->col, $this->dtype);
        for ($i = 0; $i < $this->ndim; ++$i) {
            $ar->data[$i] = $this->data[$i] > $s ? 1 : 0;
        }
        return $ar;
    }

    /**
     * 
     * @param int|float|matrix $m
     * @return matrix
     */
    public function less(int|float|matrix $m): matrix {
        $ar = self::factory($this->row, $this->col, $this->dtype);
        if ($m instanceof self) {
            if ($this->checkShape($m)) {
                for ($i = 0; $i < $this->ndim; ++$i) {
                    $ar->data[$i] = $this->data[$i] < $m->data[$i] ? 1 : 0;
                }
                return $ar;
            }
        } else {
            for ($i = 0; $i < $this->ndim; ++$i) {
                $ar->data[$i] = $this->data[$i] < $m ? 1 : 0;
            }
            return $ar;
        }
    }

    /**
     * Is the matrix symmetric i.e. is it equal to its own transpose?
     *
     * @return bool
     */
    public function isSymmetric(): bool {
        if (!$this->isSquare()) {
            return false;
        }
        $ar = $this->transpose();
        for ($i = 0; $i < $ar->ndim; ++$i) {
            if ($ar->data[$i] != $this->data[$i]) {
                unset($ar);
                return false;
            }
        }
        unset($ar);
        return true;
    }

    /**
     * print the matrix in consol
     */
    public function printMatrix() {
        echo __CLASS__ . PHP_EOL;
        for ($i = 0; $i < $this->row; ++$i) {
            for ($j = 0; $j < $this->col; ++$j) {
                printf('%lf  ', $this->data[$i * $this->col + $j]);
            }
            echo PHP_EOL;
        }
    }

    public function __toString() {
        return (string) $this->printMatrix();
    }

    protected function flattenArray(array $ar) {
        if (is_array($ar) && is_array($ar[0])) {
            $a = [];
            foreach ($ar as $y => $value) {
                foreach ($value as $k => $v) {
                    $a[] = $v;
                }
            }
            return $a;
        }
    }

    protected function checkShape(matrix $matrix) {
        if ($this->row != $matrix->row || $this->col != $matrix->col) {
            self::_dimensionaMisMatchErr('mismatch Dimensions of given matrix');
        }
        return true;
    }

    protected function checkDimensions(matrix $matrix) {
        if ($this->col != $matrix->row) {
            self::_dimensionaMisMatchErr('Mismatch Dimensions of given matrix! Matrix-A col & Matrix-B row amount need to be the same');
        }
        return true;
    }

    protected function checkDtype(matrix $matrix) {
        if ($this->dtype != $matrix->dtype) {
            self::_dtypeErr('mismatch dtype of given matrix');
        }
        return true;
    }
    
    protected function __construct(int $row, int $col, int $dtype = self::Float) {
        if ($row < 1 || $col < 1) {
            self::_invalidArgument('* To create Numphp/Matrix row & col must be greater than 0!, Op Failed! * ');
        }
        parent::__construct($row*$col, $dtype);
        $this->row = $row;
        $this->col = $col;
        return $this;
    }
}
