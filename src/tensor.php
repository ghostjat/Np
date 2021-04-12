<?php

declare(strict_types=1);

namespace numphp;

use FFI,
    InvalidArgumentException;

if (!extension_loaded('FFI')) {
    die('FFI extension required');
}

/** A fast lite memory efficient Scientific Computing in php
 * 
 * @category Scientific Computing
 * @package numphp/tensor
 * @author ghostjat
 * @version 0.0.1.alpha
 */
class tensor {

    const TWO_PI = 2. * M_PI, EPSILON = 1e-8;
    const CblasLeft = 141, CblasRight = 142;
    const CblasUpper = 121, CblasLower = 122;
    const CblasNonUnit = 131, CblasUnit = 132;
    const CblasRowMajor = 101, CblasColMajor = 102;
    const CblasNoTrans = 111, CblasTrans = 112, CblasConjTrans = 113;
    const Float = 1, Double = 2;

    public $tp, $row, $col, $view, $type;
    public static $ffi_blas = null, $_time = null, $_mem = null;

    /**
     * create empty vector/2d matrix
     * @param int $row
     * @param int $col
     * @return \numphp\tensor
     */
    public static function factory(int $row, int $col = null): tensor {
        return new self($row, $col);
    }

    /**
     * create vector/2d matrix using php array
     * @param array $data
     * @return \numphp\tensor
     */
    public static function ar(array $data): tensor {
        if (is_array($data) && is_array($data[0])) {
            $ar = self::factory(count($data), count($data[0]));
            $ar->setData($data);
        } else {
            $row = count($data);
            $ar = self::factory($row, null);
            for ($i = 0; $i < $row; ++$i) {
                $ar->tp[$i] = $data[$i];
            }
        }
        return $ar;
    }

    /**
     * Return Tensor with random values
     * @param int $row
     * @param int $col
     * @param const $type
     * @return \blas\tensor
     */
    public static function randn(int $row, int $col): tensor {
        $ar = self::factory($row, $col);
        $max = getrandmax();
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->tp[$i * $col + $j] = rand() / $max;
            }
        }
        return $ar;
    }

    /**
     * Return vector with random values
     * @param int $row
     * @return \numphp\tensor
     */
    public static function randnVector(int $row): tensor {
        $ar = self::factory($row);
        $max = getrandmax();
        for ($i = 0; $i < $ar->row; ++$i) {
            $ar->tp[$i] = rand() / $max;
        }
        return $ar;
    }

    /**
     * create a 2d uniform matrix
     * @param int $row
     * @param int $col
     * @return \blas\tensor
     */
    public static function uniform(int $row, int $col): tensor {
        $ar = self::factory($row, $col);
        $max = getrandmax();
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->tp[$i * $col + $j] = rand(-$max, $max) / $max;
            }
        }
        return $ar;
    }

    /**
     * cretae zero like 2d matrix
     * @param int $row
     * @param int $col
     * @return \blas\tensor
     */
    public static function zero(int $row, int $col): tensor {
        $ar = self::factory($row, $col);
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->tp[$i * $col + $j] = 0;
            }
        }
        return $ar;
    }

    /**
     * create one like 2d matrix
     * @param int $row
     * @param int $col
     * @return \blas\tensor
     */
    public static function ones(int $row, int $col): tensor {
        $ar = self::factory($row, $col, self::Double);
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->tp[$i * $col + $j] = 1;
            }
        }
        return $ar;
    }

    /**
     * create a null like 2d matrix
     * @param int $row
     * @param int $col
     * @return \blas\tensor
     */
    public static function null(int $row, int $col): tensor {
        $ar = self::factory($row, $col, self::Double);
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->tp[$i * $col + $j] = null;
            }
        }
        return $ar;
    }

    /**
     * create a 2d matrix with given scalar value
     * @param int $row
     * @param int $col
     * @param int|float|double $val
     * @return \blas\tensor
     */
    public static function full(int $row, int $col, $val): tensor {
        $ar = self::factory($row, $col);
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->tp[$i * $col + $j] = $val;
            }
        }
        return $ar;
    }

    /**
     * 
     * @param array $elements
     * @return \numphp\tensor
     */
    public static function diagonal(array $elements): tensor {
        $n = count($elements);
        $ar = self::factory($n, $n);
        for ($i = 0; $i < $n; ++$i) {
            for ($j = 0; $j < $n; ++$j) {
                $ar->tp[$i * $n + $j] = $i === $j ? $elements[$i] : 0;
            }
        }

        return $ar;
    }

    /**
     * 
     * @param int $row
     * @param int $col
     * @param float $lambda
     * @return \numphp\tensor
     */
    public static function poisson(int $row, int $col, float $lambda = 1.0): tensor {
        $max = getrandmax();
        $l = exp(-$lambda);
        $a = [];
        while (count($a) < $row) {
            $rowA = [];

            while (count($rowA) < $col) {
                $k = 0;
                $p = 1.0;
                while ($p > $l) {
                    ++$k;
                    $p *= rand() / $max;
                }
                $rowA[] = $k - 1;
            }
            $a[] = $rowA;
        }

        return self::ar($a);
    }

    /**
     * 
     * @param int $row
     * @param int $col
     * @return \blas\tensor
     */
    public static function gaussian(int $row, int $col): tensor {
        $max = getrandmax();
        $a = $extras = [];

        while (count($a) < $row) {
            $rowA = [];

            if ($extras) {
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

        return self::ar($a);
    }

    /**
     * 
     * @param int $n
     * @return \numphp\tensor
     * @throws InvalidArgumentException
     */
    public static function identity(int $n): tensor {
        if ($n < 1) {
            throw new InvalidArgumentException('Dimensionality must be greater than 0 on all axes.');
        }

        $ar = self::factory($n, $n);
        for ($i = 0; $i < $n; ++$i) {
            for ($j = 0; $j < $n; ++$j) {
                $ar->tp[$i * $n + $j] = $i === $j ? 1 : 0;
            }
        }
        return $ar;
    }

    public function setData($data) {
        if (!is_null($data)) {
            if ($this->type == 0) {
                if (is_array($data) && is_array($data[0]) && $this->type == 0) {
                    for ($i = 0; $i < $this->row; ++$i) {
                        for ($j = 0; $j < $this->col; ++$j) {
                            $this->tp[$i * $this->col + $j] = $data[$i][$j];
                        }
                    }
                } elseif (is_numeric($data)) {
                    for ($i = 0; $i < $this->row; ++$i) {
                        for ($j = 0; $j < $this->col; ++$j) {
                            $this->tp[$i * $this->col + $j] = $data[$i][$j];
                        }
                    }
                }
            } elseif ($this->type == 1) {
                if (is_array($data)) {
                    for ($i = 0; $i < $this->col; ++$i) {
                        $this->tp[$i] = $data[$i];
                    }
                } elseif (is_numeric($data)) {
                    for ($i = 0; $i < $this->col; ++$i) {
                        $this->tp[$i] = $data;
                    }
                }
            }
        }
    }

    /**
     * matrix dot product
     * @param \numphp\tensor $matrix
     * @return \numphp\tensor
     */
    public function dotMatrix(\numphp\tensor $matrix): tensor {
        $ar = self::factory($this->row, $this->col);
        $this->_init();
        self::$ffi_blas->cblas_dgemm(self::CblasRowMajor, self::CblasNoTrans,
                self::CblasNoTrans, $this->row, $this->col,
                $this->col, 1.0, $this->tp,
                $this->col, $matrix->tp, $this->row,
                0.0, $ar->tp, $this->row);
        return $ar;
    }

    /**
     * Matrix vector multiplication
     * @param \numphp\tensor $vecotr
     * @return \numphp\tensor
     */
    public function mul_MatrixVector(\numphp\tensor $vecotr): tensor {
        $this->_init();
        $ar = self::factory($this->row, null);
        self::$ffi_blas->cblas_dgemv(self::CblasRowMajor,
                self::CblasNoTrans, $this->row, $this->col,
                0.5, $this->tp, $this->row,
                $vecotr->tp, 1, 1.0,
                $ar->tp, 1);
        return $ar;
    }

    /**
     * Matrix Scalar multiplication
     * @param int|float $scalar
     * @return \numphp\tensor
     */
    public function mul_MatrixScalar($scalar): tensor {
        if (!is_int($scalar) and!is_float($scalar)) {
            $this->_err('Scalar must be an integer or float, ' . gettype($scalar) . ' found.');
        }

        if ($scalar == 0) {
            return self::zeros($this->row, $this->col);
        }

        $ar = self::factory($this->row, $this->col);
        for ($i = 0; $i < $this->row; ++$i) {
            for ($j = 0; $j < $this->col; ++$j) {
                $ar->tp[$i * $this->col + $j] = $this->tp[$i * $this->col + $j] * $scalar;
            }
        }

        return $ar;
    }

    /**
     * sum given matrix 
     * @param \numphp\tensor $matrix
     * @return \numphp\tensor
     */
    public function sum_Matrix(\numphp\tensor $matrix): tensor {
        if ($this->row != $matrix->row || $this->col != $matrix->col) {
            $this->_err('Inavlid matrix size');
        }
        $ar = self::factory($this->row, $this->col);
        for ($i = 0; $i < $this->row; ++$i) {
            for ($j = 0; $j < $this->col; ++$j) {
                $ar->tp[$i * $this->col + $j] = $this->tp[$i * $this->col + $j] + $matrix->tp[$i * $this->col + $j];
            }
        }
        return $ar;
    }

    /**
     * subtract given matrix
     * @param \numphp\tensor $matrix
     * @return \numphp\tensor
     */
    public function subtract_Matrix(\numphp\tensor $matrix): tensor {
        if ($this->row != $matrix->row || $this->col != $matrix->col) {
            $this->_err('Inavlid matrix size');
        }
        $ar = self::factory($this->row, $this->col);
        for ($i = 0; $i < $this->row; ++$i) {
            for ($j = 0; $j < $this->col; ++$j) {
                $ar->tp[$i * $this->col + $j] = $this->tp[$i * $this->col + $j] - $matrix->tp[$i * $this->col + $j];
            }
        }
        return $ar;
    }

    /**
     * vector-vector dot product
     * @param \numphp\tensor $vector
     * @param int $incX
     * @param int $incY
     * @return type
     */
    public function dot_vector(\numphp\tensor $vector, int $incX = 1, int $incY = 1) {
        $this->_init();
        return self::$ffi_blas->cblas_ddot($this->col, $this->tp, $incX,
                        $vector->tp, $incY);
    }

    /**
     * 
     * @param int $cols
     * @return \numphp\tensor
     */
    public function diminish_left(int $cols): tensor {
        $ar = self::factory($this->row, $cols);
        for ($i = 0; $i < $ar->row; ++$i) {
            for ($j = 0; $j < $ar->col; ++$j) {
                $ar->tp[$i * $ar->col + $j] = $this->tp[$i * $ar->col + $j];
            }
        }
        return $ar;
    }

    /**
     * 
     * @param int $cols
     * @return \numphp\tensor
     */
    public function diminish_right(int $cols): tensor {
        $ar = self::factory($this->row, $cols);
        for ($i = 0; $i < $ar->row; ++$i) {
            for ($j = 0; $j < $ar->col; ++$j) {
                $ar->tp[$i * $ar->col + $j] = $this->tp[$i * $this->cols - $cols + $j];
            }
        }
        return $ar;
    }

    /**
     *  is row zero
     * @param int $row
     * @return bool
     */
    public function is_rowZero(int $row): bool {
        for ($i = 0; $i < $this->col; ++$i) {
            if ($this->tp[$row * $this->col + $i] != 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * transpose the matrix
     * @return \numphp\tensor
     */
    public function transpose(): tensor {
        $ar = self::factory($this->col, $this->row);
        for ($i = 0; $i < $ar->row; ++$i) {
            for ($j = 0; $j < $ar->col; ++$j) {
                $ar->tp[$i * $ar->col + $j] = $this->tp[$j * $ar->col + $i];
            }
        }
        return $ar;
    }
    
    /**
     * swap specific values in tensor
     * @param int $index1
     * @param int $index2
     */
    public function swapValue(int $index1, int $index2){
        $tmp = $this->tp[$index1];
        $this->tp[$index1] = $this->tp[$index2];
        $this->tp[$index2] = $tmp;
    }
    
    /**
     * 
     * @param float $c
     * @return \numphp\tensor
     */
    public function scale(float $c) :tensor {
        $ar = $this->copyTensor();
        for($i = 0; $i < $ar->row; ++$i) {
            for($j = 0; $j < $ar->col; ++$j) {
                $ar->tp[$i * $ar->col + $j] *= $c;
            }
        }
        return $ar;
    }

    /**
     * make copy the tensor
     * @return \numphp\tensor
     */
    public function copyTensor(): tensor {
        return clone $this;
    }

    /**
     * get the shape of tensor
     * @return object
     */
    public function getShape(): object {
        return (object) ['m' => $this->row, 'n' => $this->col];
    }

    /**
     * get the type of tensor
     * @return int
     */
    public function getType(): int {
        return $this->type;
    }

    /**
     * set Timer, get total time 
     */
    public static function time() {
        if (is_null(self::$_time)) {
            self::$_time = microtime(true);
        } else {
            echo 'Time-Consumed:- ' . (microtime(true) - self::$_time) . PHP_EOL;
        }
    }

    /**
     * set memory dog, get total memory
     */
    public static function getMemory() {
        if (is_null(self::$_mem)) {
            self::$_mem = memory_get_usage();
        } else {
            $memory = memory_get_usage() - self::$_mem;
            $unit = ['b', 'kb', 'mb', 'gb', 'tb', 'pb'];
            echo round($memory / pow(1024, ($i = floor(log($memory, 1024)))), 2) . $unit[$i] . PHP_EOL;
        }
    }

    /**
     * print the tensor in consol
     */
    public function printTensor() {
        echo __CLASS__ . PHP_EOL;
        if ($this->type == 0) {
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    printf('%lf  ', $this->tp[$i * $this->col + $j]);
                }
                echo PHP_EOL;
            }
        } elseif ($this->type == 1) {
            for ($j = 0; $j < $this->col; ++$j) {
                printf('%lf  ', $this->tp[$j]);
            }
            echo PHP_EOL;
        }
    }

    public function __toString() {
        return (string) $this->printTensor();
    }

    /**
     * c double* matrix
     * @param int $row
     * @param int $col
     * @return \ffi\cdata
     */
    protected static function c_Matrix(int $row, int $col) {
        return \FFI::cast('double *', FFI::new("double[$row][$col]"));
    }

    /**
     * c double* vector
     * @param int $col
     * @return \ffi\cdata
     */
    protected static function c_Vector(int $col) {
        return \FFI::cast('double *', FFI::new("double[$col]"));
    }

    protected static function castDouble($cdata) {
        return \FFI::cast('double *', $cdata);
    }

    protected static function castInt($cdata) {
        return FFI::cast('int*', $cdata);
    }

    protected static function castFloat($cdata) {
        return FFI::cast('float*', $cdata);
    }

    protected function free(FFI\CData $ptr) {
        return FFI::free($ptr);
    }

    protected function __construct(int $row, int $col = null) {
        if ($row < 1) {
            throw new InvalidArgumentException(' * To create Numphp/tensor row must be greater than 0!, Op Failed! * ');
        }
        if (is_null($col)) {
            $this->type = 1;
            $this->row = $row;
            $this->tp = self::c_Vector($this->row);
        } else {
            $this->type = 0;
            $this->row = $row;
            $this->col = $col;
            $this->tp = self::c_Matrix($this->row, $this->col);
        }
        return $this;
    }

    protected function _init() {
        if (is_null(self::$ffi_blas)) {
            self::$ffi_blas = \FFI::load(__DIR__ . '/blas.h');
        }
        return self::$ffi_blas;
    }

    private static function _err($msg) {
        throw new \Exception($msg);
    }

}
