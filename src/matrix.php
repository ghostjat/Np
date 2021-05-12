<?php

declare(strict_types=1);

namespace numphp;

/**
 * A fast lite memory efficient Scientific Computing for php
 * Matrix
 * @package NumPhp\Matrix
 * @category Scientific Computing
 * @author ghost (Shubham Chaudhary)
 * @email ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */
class matrix {

    const TWO_PI = 2. * M_PI, EPSILON = 1e-8;
    const FLOAT = 1, DOUBLE = 2, INT = 3;
    protected $l,$u,$p;
    public $data, $row, $col, $dtype;
    public static $_time = null, $_mem = null;

    /**
     * create empty 2d matrix for given data type
     * @param int $row  num of rows 
     * @param int $col  num of cols
     * @param const $dtype matrix data type float|double
     * @return \numphp\matrix
     */
    public static function factory(int $row, int $col, int $dtype = self::FLOAT): matrix {
        return new self($row, $col, $dtype);
    }

    /**
     * create 2d matrix using php array
     * @param array $data
     * @param const $dtype matrix data type float|double
     * @return \numphp\matrix
     */
    public static function ar(array $data, int $dtype = self::FLOAT): matrix {
        if (is_array($data) && is_array($data[0])) {
            $ar = self::factory(count($data), count($data[0]), $dtype);
            $ar->setData($data);
        } else {
            self::_err('data must be of same dimensions');
        }
        return $ar;
    }

    /**
     * Create Matrix with random values
     * @param int $row
     * @param int $col
     * @param int $dtype  Float|Double
     * @return \numphp\matrix
     */
    public static function randn(int $row, int $col, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        $max = getrandmax();
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->data[$i * $col + $j] = rand() / $max;
            }
        }
        return $ar;
    }

    /**
     * Return 2d matrix with uniform values
     * @param int $row
     * @param int $col
     * @param int $dtype
     * @return \numphp\matrix
     */
    public static function uniform(int $row, int $col, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        $max = getrandmax();
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->data[$i * $col + $j] = rand(-$max, $max) / $max;
            }
        }
        return $ar;
    }

    /**
     * 
     * @param int $row
     * @param int $col
     * @param int $dtype
     * @return \numphp\matrix
     */
    public static function zeros(int $row, int $col, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->data[$i * $col + $j] = 0;
            }
        }
        return $ar;
    }

    /**
     * create one like 2d matrix
     * @param int $row
     * @param int $col
     * @return \numphp\matrix
     */
    public static function ones(int $row, int $col, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->data[$i * $col + $j] = 1;
            }
        }
        return $ar;
    }

    /**
     * create a null like 2d matrix
     * @param int $row
     * @param int $col
     * @return \numphp\matrix
     */
    public static function null(int $row, int $col, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->data[$i * $col + $j] = null;
            }
        }
        return $ar;
    }

    /**
     * create a 2d matrix with given scalar value
     * @param int $row
     * @param int $col
     * @param int|float|double $val
     * @return \numphp\matrix
     */
    public static function full(int $row, int $col, $val, int $dtype = self::FLOAT): matrix {
        $ar = self::factory($row, $col, $dtype);
        for ($i = 0; $i < $row; ++$i) {
            for ($j = 0; $j < $col; ++$j) {
                $ar->data[$i * $col + $j] = $val;
            }
        }
        return $ar;
    }

    /**
     * create a diagonal 2d matrix with given 1d array;
     * @param array $elements
     * @return \numphp\matrix
     */
    public static function diagonal(array $elements, int $dtype = self::FLOAT): matrix {
        $n = count($elements);
        $ar = self::factory($n, $n,$dtype);
        for ($i = 0; $i < $n; ++$i) {
            $ar->data[$i * $n + $i] = $elements[$i]; #for ($j = 0; $j < $n; ++$j) {$i === $j ? $elements[$i] : 0;#} 
        }
        return $ar;
    }

    /**
     * 
     * @param int $row
     * @param int $col
     * @param float $lambda
     * @return \numphp\matrix
     */
    public static function poisson(int $row, int $col, float $lambda = 1.0): matrix {
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
     * @return \numphp\matrix
     */
    public static function gaussian(int $row, int $col): matrix {
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
     * create an identity matrix with the given dimensions.
     * @param int $n
     * @return \numphp\matrix
     * @throws InvalidArgumentException
     */
    public static function identity(int $n): matrix {
        if ($n < 1) {
            throw new InvalidArgumentException('Dimensionality must be greater than 0 on all axes.');
        }

        $ar = self::factory($n, $n);
        for ($i = 0; $i < $n; ++$i) {
            for ($j = 0; $j < $n; ++$j) {
                $ar->data[$i * $n + $j] = $i === $j ? 1 : 0;
            }
        }
        return $ar;
    }

    /**
     * trace
     * @return float
     */
    public function trace() {
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

    /**
     * get matrix dot product
     * @param \numphp\matrix $matrix
     * @return \numphp\matrix
     */
    public function dotMatrix(\numphp\matrix $matrix): matrix {
        if ($this->checkDtype($matrix) && $this->checkDimensions($matrix)) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            if ($this->dtype == self::FLOAT) {
                core\blas::sgemm($this, $matrix, $ar);
            } else {
                core\blas::dgemm($this, $matrix, $ar);
            }

            return $ar;
        }
    }

    /**
     * get multiplication of matrix vector
     * @param \numphp\vector $vector
     * @return \numphp\matrix
     */
    public function mulMatrixVector(\numphp\vector $vector): matrix {
        if ($this->dtype != $vector->dtype) {
            self::_err('Mismatch Dtype of given vector');
        }
        $mvr = self::factory($this->row, $this->col, $this->dtype);
        if ($this->dtype == self::FLOAT) {
            core\blas::sgemv($this, $vector, $mvr);
        } else {
            core\blas::dgemv($this, $vector, $mvr);
        }

        return $mvr;
    }

    /**
     * multiply this matrix with another matrix|scalar element-wise
     * Matrix Scalar\Matrix multiplication
     * @param int|float|matrix|vector $value
     * @return mixed (matrix,vector)
     */
    public function multiply($value) {
        if ($value instanceof self) {
            if ($this->checkDtype($value) && $this->checkDimensions($value)) {
                $ar = $this->copyMatrix();
                for ($i = 0; $i < $this->row; ++$i) {
                    for ($j = 0; $j < $this->col; ++$j) {
                        $ar->data[$i * $this->col + $j] *= $value->data[$i * $this->col + $j];
                    }
                }
                return $ar;
            }
        } else if ($value instanceof vector) {
            return $value->mulVectorMatrix($this);
        } else {
            if (!is_int($value) || !is_float($value)) {
                self::_err('Scalar must be an integer or float, ' . gettype($value) . ' found.');
            }

            if ($value == 0) {
                return self::zeros($this->row, $this->col, $this->dtype);
            }

            $ar = $this->copyMatrix();
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] *= $value;
                }
            }

            return $ar;
        }
    }

    /**
     * Add another matrix or a scalar to current matrix
     * 
     * @param mixed(scalar,\numphp\matrix) $value matrix|$scalar to add this matrix
     * @return \numphp\matrix
     */
    public function add($value): matrix {
        if ($value instanceof self) {
            if ($this->row != $value->row || $this->col != $value->col) {
                self::_err('Inavlid matrix size');
            }
            $ar = $this->copyMatrix();
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] +=  $value->data[$i * $this->col + $j];
                }
            }
            return $ar;
        }
        else {
            $ar = $this->copyMatrix();
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] += $value;
                }
            }
            return $ar;
        }
    }

    /**
     * subtract another matrix or a scalar to this matrix
     * @param mixed(scalar,\numphp\matrix) $value matrix|$scalar to add this matrix
     * @return \numphp\matrix
     */
    public function subtract($value): matrix {
        if ($value instanceof self) {
            if ($this->row != $value->row || $this->col != $value->col) {
                self::_err('Inavlid matrix size');
            }
            $ar = $this->copyMatrix();
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] -= $value->data[$i * $this->col + $j];
                }
            }
            return $ar;
        }
        else{
            $ar = $this->copyMatrix();
            for ($i = 0; $i < $this->row; ++$i) {
                for ($j = 0; $j < $this->col; ++$j) {
                    $ar->data[$i * $this->col + $j] -= $value;
                }
            }
            return $ar;
        }
    }
    
    /**
     * 
     * @param int|float $d
     * @return bool
     */
    public static function is_zero($d): bool {
        if(abs($d) < self::EPSILON){
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
     * transpose the matrix
     * @return \numphp\matrix
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
     * swap specific values in tensor
     * @param int $index1
     * @param int $index2
     */
    public function swapValue(int $index1, int $index2) {
        $tmp = $this->data[$index1];
        $this->data[$index1] = $this->data[$index2];
        $this->data[$index2] = $tmp;
    }

    /**
     * swap specific rows in tensor
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
     * swap specific cols in tensor
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
     * 
     * @param float $c
     * @return \numphp\matrix
     */
    public function scale(float $c): matrix {
        $ar = $this->copyMatrix();
        for ($i = 0; $i < $ar->row; ++$i) {
            for ($j = 0; $j < $ar->col; ++$j) {
                $ar->data[$i * $ar->col + $j] *= $c;
            }
        }
        return $ar;
    }

    /**
     * 
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
     * 
     * @param \numphp\matrix $matrix
     * @return \numphp\matrix
     */
    public function augment(\numphp\matrix $matrix): matrix {
        if ($this->row != $matrix->row) {
            self::_err('Error::Invalid size!');
        }
        $col = $this->col + $matrix->col;
        $ar = self::factory($this->row, $col, $this->dtype);
        for ($i = 0; $i < $this->row; ++$i) {
            for ($j = 0; $j < $this->col; ++$j) {
                $ar->data[$i * $col + $j] = $this->data[$i * $this->col + $j];
            }
            for ($j = 0; $j < $matrix->col; ++$j) {
                $ar->data[$i * $col + ($this->col + $j)] = $matrix->data[$i * $matrix->col + $j];
            }
        }
        return $ar;
    }

    /**
     * RREF
     * The reduced row echelon form (RREF) of a matrix.
     * @return \numphp\matrix
     */
    public function rref(): matrix {
        $lead = 0;
        $ar = $this->copyMatrix();
        for ($r = 0; $r < $ar->row; ++$r) {
            if ($lead >= $ar->col)
                break;
            {
                $i = $r;
                while ($ar->data[$i * $ar->col + $lead] == 0) {
                    $i++;
                    if ($i == $ar->row) {
                        $i = $r;
                        $lead++;
                        if ($lead == $ar->col) {
                            return $ar;
                        }
                    }
                }
                $ar->swapRows($r, $i);
            }
            {
                $lv = $ar->data[$r * $ar->col +$lead];
                for ($j = 0; $j < $ar->col; ++$j) {
                    $ar->data[$r * $ar->col + $j] = $ar->data[$r * $ar->col + $j] / $lv;
                }
            }
            for ($i = 0; $i < $ar->row; ++$i) {
                if ($i != $r) {
                    $lv = $ar->data[$i * $ar->col + $lead];
                    for ($j = 0; $j < $ar->col; ++$j) {
                        $ar->data[$i * $ar->col + $j] -= $lv * $ar->data[$r * $ar->col + $j];
                    }
                }
            }
            $lead++;
        }
        return $ar;
    }

    public function determinant() {
        if(!$this->isSquare()) {
            self::_err('Error:: Determinant of non-square matrix!');
        }
        elseif ($this->row == 1) {
            return $this->data[0];
        }
        $smaller = array_fill(0, ($this->row -1) * ($this->col -1), null);
        $b = self::factory($this->row - 1, $this->col - 1);
        $det = 0.0; $sign = 1;
        for($col = 0; $col < $this->col ; ++$col)  {
            $n =0;
            for($i = 0; $i < $this->row; ++$i) {
                for($j = 0; $j < $this->col; ++$j) {
                    if($j == $i) {
                        continue;
                    }
                    $smaller[$n] = $this->data[$i * $this->col +$j];
                    ++$n;
                }
            }
            $b->setData($smaller);
            $det += $sign * $this->data[0 * $this->col + $col];
        }
    }

    /**
     * make copy of the matrix
     * @return \numphp\matrix
     */
    public function copyMatrix(): matrix {
        return clone $this;
    }
    

    /**
     * 
     * @param int $cols
     * @return \numphp\matrix
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
     * @return \numphp\matrix
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
     * @return \numphp\vector int
     */
    public function argMax(): vector {
        $v = vector::factory($this->row, vector::INT);
        if($this->dtype === self::DOUBLE ){
            for($i = 0; $i < $this->row; ++$i) {
                $v->data[$i] = core\blas::dmax($this->rowAsVector($i));
            }
        }
        elseif($this->dtype === self::FLOAT) {
            for($i = 0; $i < $this->row; ++$i) {
                $v->data[$i] = core\blas::smax($this->rowAsVector($i));
            }
        }
        return $v;
    }

    /**
     * Return the index of the minimum element in every row of the matrix.
     * @return \numphp\vector int
     */
    public function argMin(): vector {
        $v = vector::factory($this->row, vector::INT);
        if($this->dtype === self::DOUBLE ){
            for($i = 0; $i < $this->row; ++$i) {
                $v->data[$i] = core\blas::dmin($this->rowAsVector($i));
            }
        }
        elseif($this->dtype === self::FLOAT) {
            for($i = 0; $i < $this->row; ++$i) {
                $v->data[$i] = core\blas::smin($this->rowAsVector($i));
            }
        }
        return $v;
    }

    /**
     * set data to matrix
     * @param type $data
     */
    public function setData($data,$dignoal=false) {
        if (!is_null($data) && $dignoal == false) {
            if (is_array($data) && is_array($data[0])) {
                for ($i = 0; $i < $this->row; ++$i) {
                    for ($j = 0; $j < $this->col; ++$j) {
                        $this->data[$i * $this->col + $j] = $data[$i][$j];
                    }
                }
            } elseif (is_numeric($data) && $dignoal == false) {
                for ($i = 0; $i < $this->row; ++$i) {
                    for ($j = 0; $j < $this->col; ++$j) {
                        $this->data[$i * $this->col + $j] = $data;
                    }
                }
            }
            elseif(is_numeric($data) && $dignoal == true) {
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
     * @return \numphp\vector
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
     * @return \numphp\vector
     */
    public function colAsVector(int $index) : vector {
        $vr = vector::factory($this->row, $this->dtype);
        for($i = 0; $i< $this->row; ++$i) {
            $vr->data[$i] = $this->data[$i * $this->row + $index];
        }
        return $vr;
    }

    /**
     * Return the diagonal elements of a square matrix as a vector.
     * @return \numphp\vector
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
     * get a diagonal value from matrix
     * @param int $i
     * @return type
     */
    public function getDiagonalVal(int $i) {
        if($this->isSquare()){
            return $this->data[$i * $this->row + $i];
        }
    }

    /**
     * Compute the multiplicative inverse of the matrix.
     * @return \numphp\matrix
     */
    public function inverse():matrix|null {
        if(!$this->isSquare()){
            self::_err('Error::invalid Size of matrix!');
        }
        $imat = $this->copyMatrix();
        $ipiv = vector::factory($this->row, vector::INT);
        $lp = core\lapack::sgetrf($imat, $ipiv);
        if($lp != 0) {
            return null;
        }
        $lp = core\lapack::sgetri($imat, $ipiv);
        if($lp != 0) {
            return null;
        }
        unset($ipiv);
        unset($lp);
        return $imat;
    }
    
    
    public function determinate() {
        if(!$this->isSquare()) {
            self::_err('determinant is undefined for a non square matrix');
        }
    }
    
    
    /**
     *Run a function over all of the elements in the matrix. 
     * @param callable $func
     * @return \numphp\matrix
     */
    public function map(callable $func) : matrix {
        if($func instanceof \Closure) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
            for($i=0;$i<$this->row;++$i) {
                for($j=0;$j<$this->col;++$j) {
                    $ar->data[$i * $this->col + $j] = $func($this->data[$i * $this->col + $j]); 
                }
            }
            return $ar;
        }
    }
    
     /**
     * Is the matrix symmetric i.e. is it equal to its own transpose?
     *
     * @return bool
     */
    public function isSymmetric() : bool {
        if(!$this->isSquare()){
            return false;
        }
        $ar = $this->transpose();
        for($i = 0; $i < $this->row; ++$i) {
            for($j = 0 ; $j < $this->col; ++ $j) {
                if($ar->data[$i * $this->col + $j] != $this->data[$i * $this->col +$j]) {
                    unset($ar);
                    return false;
                }
            }
        }
        unset($ar);
        return true;
    }


    /**
     * Flatten i.e unravel the matrix into a vector.
     *
     * @return \numphp\vector
     */
    public function flatten():vector {
        $vr = vector::factory($this->row * $this->col, $this->dtype);
        for($i = 0; $i< $this->row; ++$i) {
            for($j = 0; $j < $this->col; ++$j) {
                $vr->data[$i* $this->col + $j] = $this->data[$i * $this->col +$j];
            }
        }
        return $vr;
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

    protected function checkDimensions(\numphp\matrix $matrix) {
        if ($this->row != $matrix->row || $this->col != $matrix->col) {
            self::_err('Mismatch Dimensions of given matrix');
        }
        return true;
    }

    protected function checkDtype(\numphp\matrix $matrix) {
        if ($this->dtype != $matrix->dtype) {
            self::_err('Mismatch Dtype of given matrix');
        }
        return true;
    }

    protected static function c_FloatMatrix(int $row, int $col) {
        return \FFI::cast('float *', \FFI::new("float[$row][$col]"));
    }

    protected static function c_DoubleMatrix(int $row, int $col) {
        return \FFI::cast('double *', \FFI::new("double[$row][$col]"));
    }
    
    protected static function c_IntMatrix(int $row, int $col) {
        return \FFI::cast('int *', \FFI::new("int[$row][$col]"));
    }
    
    protected function __construct(int $row, int $col, int $dtype = self::Float) {
        if ($row < 1 || $col < 1) {
            $this->_invalidArgument('* To create Numphp/Matrix row & col must be greater than 0!, Op Failed! * ');
        }
        $this->row = $row;
        $this->col = $col;
        $this->dtype = $dtype;
        switch ($dtype) {
            case self::FLOAT:
                $this->data = self::c_FloatMatrix($this->row, $this->col);
                break;
            case self::DOUBLE:
                $this->data = self::c_DoubleMatrix($this->row, $this->col);
                break;
            case self::INT:
                $this->data = self::c_IntMatrix($this->row, $this->col);
                break;
            default :
                $this->_invalidArgument('given dtype is not supported by numphp');
                break;
        }
        return $this;
    }
    

    private static function _err($msg): \Exception {
        throw new \Exception($msg);
    }

    private function _invalidArgument($argument): \InvalidArgumentException {
        throw new InvalidArgumentException($argument);
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

}
