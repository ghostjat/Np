<?php

declare(strict_types=1);

namespace Np;

use Np\core\{nd,blas,lapack};
use Np\exceptions\invalidArgumentException;

/** 
 * Vector (rank-1)
 * 
 * @package Np
 * @category Scientific Library for PHP
 * @author ghost (Shubham Chaudhary)
 * @email ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 * 
 */
 class vector extends nd {
     use ops,linAlgb\linAlg;

    /**
     * Factory method to build a new vector.
     * 
     * @param int $col
     * @param int $dtype
     * @return vector
     */
    public static function factory(int $col, int $dtype = self::DOUBLE): vector {
        return new self($col, $dtype);
    }

    /**
     * Build a new vector from a php array.
     * 
     * @param array $data
     * @return vector
     */
    public static function ar(array $data): vector {
        if (is_array($data) && !is_array($data[0])) {
            $ar = self::factory(count($data));
            $ar->setData($data);
            return $ar;
        } else {
            self::_err('data must be of same dimensions');
        }
    }

    /**
     * Return vector with random values
     * @param int $col
     * @return vector
     */
    public static function randn(int $col): vector {
        $ar = self::factory($col);
        $max = getrandmax();
        for ($i = 0; $i < $ar->col; ++$i) {
            $ar->data[$i] = rand() / $max;
        }
        return $ar;
    }

    /**
     * Return vector with uniform values
     * @param int $col
     * @return vector
     */
    public static function uniform(int $col): vector {
        $ar = self::factory($col);
        $max = getrandmax();
        for ($i = 0; $i < $col; ++$i) {
            $ar->data[$i] = rand(-$max, $max) / $max;
        }
        return $ar;
    }

    /**
     * Build a vector of zeros with n elements.
     * 
     * @param int $col
     * @return vector
     */
    public static function zeros(int $col): vector {
        $ar = self::factory($col);
        for ($i = 0; $i < $col; ++$i) {
            $ar->data[$i] = 0;
        }
        return $ar;
    }

    /**
     * create one like vector
     * 
     * @param int $col
     * @return vector
     */
    public static function ones(int $col): vector {
        $ar = self::factory($col);
        for ($i = 0; $i < $col; ++$i) {
            $ar->data[$i] = 1;
        }
        return $ar;
    }

    /**
     * create a null like vector
     * @param int $col
     * @return vector
     */
    public static function null(int $col): vector {
        $ar = self::factory($col);
        for ($i = 0; $i < $col; ++$i) {
            $ar->data[$i] = null;
        }
        return $ar;
    }

    /**
     * create a vector with given scalar value
     * @param int $col
     * @param int|float|double $val
     * @return vector
     */
    public static function full(int $col, int|float $val): vector {
        $ar = self::factory($col);
        for ($i = 0; $i < $col; ++$i) {
            $ar->data[$i] = $val;
        }
        return $ar;
    }

    /**
     * Return evenly spaced values within a given interval.
     *
     * @param int|float $start
     * @param int|float $end
     * @param int|float $interval
     * @return vector
     */
    public static function range(int|float $start, int|float $end, int|float $interval = 1): vector {
        return self::ar(range($start, $end, $interval));
    }

    /**
     * Return a Gaussian random vector with mean 0
     * and unit variance.
     *
     * @param int $n
     * @return self
     */
    public static function gaussian(int $n): vector {
        $max = getrandmax();
        $a = [];
        while (count($a) < $n) {
            $r = sqrt(-2.0 * log(rand() / $max));
            $phi = rand() / $max * (2. * M_PI);
            $a[] = $r * sin($phi);
            $a[] = $r * cos($phi);
        }
        if (count($a) > $n) {
            $a = array_slice($a, 0, $n);
        }
        return self::ar($a);
    }

    /**
     * Generate a vector with n elements from a Poisson distribution.
     *
     * @param int $n
     * @param float $lambda
     * @return vector
     */
    public static function poisson(int $n, float $lambda = 1.0): vector {
        $max = getrandmax();
        $l = exp(-$lambda);
        $a = new self($n);
        for ($i = 0; $i < $n; ++$i) {
            $k = 0;
            $p = 1.0;
            while ($p > $l) {
                ++$k;
                $p *= rand() / $max;
            }
            $a->data[$i] = $k - 1;
        }
        return $a;
    }

    /**
     * Return a vector of n evenly spaced numbers between minimum and maximum.
     *
     * @param float $min
     * @param float $max
     * @param int $n
     * @throws invalidArgumentException
     * @return vector
     */
    public static function linspace(float $min, float $max, int $n): vector {
        if ($min > $max) {
            throw new invalidArgumentException('Minimum must be less than maximum.');
        }
        if ($n < 2) {
            throw new invalidArgumentException('Number of elements must be greater than 1.');
        }
        $k = $n - 1;
        $interval = abs($max - $min) / $k;
        $a = [$min];
        while (count($a) < $k) {
            $a[] = end($a) + $interval;
        }
        $a[] = $max;
        return self::ar($a);
    }
    
    /**
     * Return the index of the minimum element in the vector.
     * 
     * @return int
     */
    public function argMin(): int {
        return blas::min($this);
    }

    /**
     * Return the index of the maximum element in the vector.
     * 
     * @return int
     */
    public function argMax(): int {
        return blas::max($this);
    }

    /**
     * The sum of the vector.
     * @return float
     */
    public function sum(): float {
        return blas::asum($this);
    }

    /**
     * Return the product of the vector.
     * @return int|float
     */
    public function product(): float {
        $r = 1.0;
        for ($i = 0; $i < $this->col; ++$i) {
            $r *= $this->data[$i];
        }
        return $r;
    }

    /**
     * Compute the vector-matrix dot product of this vector and matrix .
     * @param \Np\matrix $m
     * @return vector
     */
    public function dotMatrix(\Np\matrix $m): vector {
        if ($this->checkDtype($this, $m)) {
            $mvr = self::factory($this->col);
            core\blas::gemv($m, $this, $mvr);
            return $mvr;
        }
    }

    /**
     * 
     * @param int|float|matrix|vector $d
     * @return matrix|vector
     */
    public function divide(int|float|matrix|vector $d): matrix|vector {
        if ($d instanceof matrix) {
            return $this->divideMatrix($d);
        }
        if ($d instanceof self) {
            return $this->divideVector($d);
        }
        return $this->divideScalar($d);
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function divideMatrix(\Np\matrix $m): matrix {
        if ($this->checkShape($this, $m)) {
            $vr = matrix::factory($m->row, $m->col);
            for ($i = 0; $i < $m->row; ++$i) {
                for ($j = 0; $j < $m->col; ++$j) {
                    $vr->data[$i * $m->col + $j] = $this->data[$j] / $m->data[$i * $m->col + $j];
                }
            }
            return $vr;
        }
    }

    /**
     * 
     * @param vector $v
     * @return vector
     */
    protected function divideVector(vector $v): vector {
        if ($this->checkShape($this, $v)) {
            $vr = self::factory($this->col);
            for ($i = 0; $i < $this->col; ++$i) {
                $vr->data[$i] = $this->data[$i] / $v->data[$i];
            }
            return $vr;
        }
    }

    /**
     * 
     * @param int|float $s
     * @return vector
     */
    protected function divideScalar(int|float $s): vector {
        $vr = self::factory($this->col);
        for ($i = 0; $i < $this->col; ++$i) {
            $vr->data[$i] = $this->data[$i] / $s;
        }
        return $vr;
    }

    /**
     * 
     * @param int|float|matrix|vector $d
     * @return matrix|vector
     */
    public function multiply(int|float|matrix|vector $d): matrix|vector {
        if ($d instanceof matrix) {
            return $this->multiplyMatrix($d);
        }
        if ($d instanceof self) {
            return $this->multiplyVector($d);
        }
        return $this->multiplyScalar($d);
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function multiplyMatrix(\Np\matrix $m): matrix {
        if ($this->checkShape($this, $m)) {
            $vr = matrix::factory($m->row, $m->col);
            for ($i = 0; $i < $m->row; ++$i) {
                for ($j = 0; $j < $m->col; ++$j) {
                    $vr->data[$i * $m->col + $j] = $this->data[$j] * $m->data[$i * $m->col + $j];
                }
            }
            return $vr;
        }
    }

    /**
     * 
     * @param \Np\vector $vector
     * @return vector
     */
    protected function multiplyVector(\Np\vector $vector): vector {
        if ($this->checkShape($this, $vector)) {
            $vr = self::factory($this->col);
            for ($i = 0; $i < $this->col; ++$i) {
                $vr->data[$i] = $this->data[$i] * $vector->data[$i];
            }
            return $vr;
        }
    }

    /**
     * 
     * @param int|float $s
     * @return vector
     */
    protected function multiplyScalar(int|float $s): vector {
        $vr = $this->copy();
        blas::scale($s, $vr);
        return $vr;
    }

    /**
     * 
     * @param int|float|matrix|vector $d
     * @return matrix|vector
     */
    public function add(int|float|matrix|vector $d): matrix|vector {
        if ($d instanceof matrix) {
            return $this->addMatrix($d);
        }
        if ($d instanceof self) {
            return $this->addVector($d);
        }
        return $this->addScalar($d);
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function addMatrix(\Np\matrix $m): matrix {
        if ($this->checkShape($this, $m)) {
            $vr = matrix::factory($m->row, $m->col);
            for ($i = 0; $i < $m->row; ++$i) {
                for ($j = 0; $j < $m->col; ++$j) {
                    $vr->data[$i * $m->col + $j] = $this->data[$j] + $m->data[$i * $m->col + $j];
                }
            }
            return $vr;
        }
    }

    /**
     * 
     * @param \Np\vector $vector
     * @return vector
     */
    protected function addVector(\Np\vector $vector): vector {
        if ($this->checkShape($this, $vector)) {
            $vr = self::factory($this->col);
            for ($i = 0; $i < $this->col; ++$i) {
                $vr->data[$i] = $this->data[$i] + $vector->data[$i];
            }
            return $vr;
        }
    }

    /**
     * 
     * @param int|float $s
     * @return vector
     */
    protected function addScalar(int|float $s): vector {
        $vr = $this->copy();
        for ($i = 0; $i < $this->col; ++$i) {
            $vr->data[$i] += $s;
        }
        return $vr;
    }

    /**
     * 
     * @param int|float|\Np\matrix|\Np\vector $d
     * @return matrix|vector
     */
    public function pow(int|float|\Np\matrix|\Np\vector $d): matrix|vector {
        if ($d instanceof matrix) {
            return $this->powMatrix($d);
        }
        if ($d instanceof vector) {
            return $this->powVector($d);
        }
        return $this->powScalar($d);
        
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function powMatrix(\Np\matrix $m): matrix {
        if ($this->checkDimensions($this, $m)) {
            $ar = matrix::factory($m->row, $m->col);
            for ($i = 0; $i < $m->row; ++$i) {
                for ($j = 0; $j < $m->col; ++$j) {
                    $ar->data[$i * $m->col + $j] = $m->data[$i * $m->col + $j] ** $this->data[$j];
                }
            }
            return $ar;
        }
    }

    /**
     * 
     * @param \Np\vector $vector
     * @return vector
     */
    protected function powVector(\Np\vector $vector): vector {
        if ($this->checkShape($this, $vector)) {
            $vr = self::factory($this->col);
            for ($i = 0; $i < $this->col; ++$i) {
                $vr->data[$i] = $this->data[$i] ** $vector->data[$i];
            }
            return $vr;
        }
    }

    /**
     * 
     * @param int|float $s
     * @return vector
     */
    protected function powScalar(int|float $s): vector {
        $v = $this->copy();
        for ($i = 0; $i < $this->col; ++$i) {
            $v->data[$i] = $v->data[$i] ** $s;
        }
        return $v;
    }

    /**
     * 
     * @param int|float|\Np\matrix|\Np\vector $d
     * @return matrix|vector
     */
    public function mod(int|float|\Np\matrix|\Np\vector $d): matrix|vector {
        if ($d instanceof matrix) {
            return $this->powMatrix($d);
        }
        if ($d instanceof vector) {
            return $this->powVector($d);
        }
        return $this->powScalar($d);    
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function modMatrix(\Np\matrix $m): matrix {
        if ($this->checkDimensions($this, $m)) {
            $ar = matrix::factory($m->row, $m->col);
            for ($i = 0; $i < $m->row; ++$i) {
                for ($j = 0; $j < $m->col; ++$j) {
                    $ar->data[$i * $m->col + $j] = $m->data[$i * $m->col + $j] % $this->data[$j];
                }
            }
            return $ar;
        }
    }

    /**
     * 
     * @param \Np\vector $vector
     * @return vector
     */
    protected function modVector(\Np\vector $vector): vector {
        if ($this->checkShape($this, $vector)) {
            $vr = self::factory($this->col);
            for ($i = 0; $i < $this->col; ++$i) {
                $vr->data[$i] = $this->data[$i] % $vector->data[$i];
            }
            return $vr;
        }
    }

    /**
     * 
     * @param int|float $s
     */
    protected function modScalar(int|float $s) {
        $v = $this->copy();
        for ($i = 0; $i < $this->col; ++$i) {
            $v->data[$i] = $v->data[$i] % $s;
        }
    }

    /**
     * 
     * @param int|float|matrix|vector $d
     * @return matrix|vector
     */
    public function subtract(int|float|matrix|vector $d): matrix|vector {
        if ($d instanceof matrix) {
            return $this->subtractMatrix($d);
        }
        if ($d instanceof self) {
            return $this->subtractVector($d);
        }
        return $this->substractScalar($d);
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function subtractMatrix(\Np\matrix $m): matrix {
        if ($this->checkShape($this, $m)) {
            $vr = matrix::factory($m->row, $m->col);
            for ($i = 0; $i < $m->row; ++$i) {
                for ($j = 0; $j < $m->col; ++$j) {
                    $vr->data[$i * $m->col + $j] = $this->data[$j] - $m->data[$i * $m->col + $j];
                }
            }
            return $vr;
        }
    }

    /**
     * 
     * @param \Np\vector $vector
     * @return vector
     */
    protected function subtractVector(\Np\vector $vector): vector {
        if ($this->checkShape($this, $vector) ) {
            $vr = self::factory($this->col);
            for ($i = 0; $i < $this->col; ++$i) {
                $vr->data[$i] = $this->data[$i] - $vector->data[$i];
            }
            return $vr;
        }
    }

    /**
     * 
     * @param \Np\vector $scalar
     * @return \Np\vector
     */
    protected function substractScalar(int|float $scalar): vector {
        $vr = self::factory($this->col);
        for ($i = 0; $i < $this->col; ++$i) {
            $vr->data[$i] = $this->data[$i] - $scalar;
        }
        return $vr;
    }

    /**
     * 
     * @param \Np\vector $v
     * @param int $stride
     * @return vector
     */
    public function convolve(\Np\vector $v, int $stride = 1): vector {
        return convolve::conv1D($this, $v, $stride);
    }

    public function max() {
        return $this->data[blas::max($this)];
    }

    public function min() {
        $this->data[blas::min($this)];
    }
    
    /**
     * Return the inner product of two vectors.
     *
     * @param \Np\vector $vector
     * 
     */
    public function inner(\Np\vector $vector) {
        return $this->dotVector($vector);
    }

    /**
     * Calculate the L1 norm of the vector.
     * @return float
     */
    public function normL1(): float {
        return $this->abs()->sum();
    }

    public function normL2() {
        return sqrt($this->square()->sum());
    }

    public function normMax() {
        return $this->abs()->max();
    }

    public function normP(float $p = 2.5) {
        if ($p <= 0.0) {
            self::_invalidArgument('P must be greater than 0.0 !');
        }
        return $this->abs()->powScalar($p)->sum() ** (1.0 / $p);
    }

    /**
     * Return the reciprocal of the vector element-wise.
     *
     * @return self
     */
    public function reciprocal(): vector {
        return self::ones($this->col)->divideVector($this);
    }
    
    /**
     * 
     * @return int|float
     */
    public function mean():int|float {
        return $this->sum()/ $this->col;
    }
    
    /**
     * 
     * @return int|float
     */
    public function median():int|float {
        $mid = intdiv($this->col, 2);

        $a = $this->copy()->sort();
        if ($this->col % 2 === 1) {
            $median = $a->data[$mid];
        } else {
            $median = ($a->data[$mid - 1] + $a->data[$mid]) / 2.;
        }
        return $median;
    }
    
    public function variance($mean = null)
    {
        if (is_null($mean)) {
            $mean = $this->mean();
        }

        $sd = $this->substractScalar($mean)->square()->sum();

        return $sd / $this->col;
    }

    /**
     * 
     * @return vector
     */
    public function square(): vector {
        return $this->multiplyVector($this);
    }

    /**
     * sort the vector 
     * @param string $type i or d
     * 
     */
    public function sort($type = 'i') {
        lapack::sort($this, $type);
        return $this;
    }

    /**
     * set data to vector
     * @param int|float|array $data
     */
    public function setData(int|float|array $data) {
        if (is_array($data) && !is_array($data[0])) {
            for ($i = 0; $i < $this->col; ++$i) {
                $this->data[$i] = $data[$i];
            }
        }
        if (is_numeric($data)) {
            for ($i = 0; $i < $this->col; ++$i) {
                $this->data[$i] = $data;
            }
        }
    }

    /**
     * get the size of vector
     * @return int
     */
    public function getSize(): int {
        return $this->col;
    }

    public function getDtype() {
        return $this->dtype;
    }

    public function asArray() {
        $ar = array_fill(0, $this->col, null);
        for ($i = 0; $i < $this->col; ++$i) {
            $ar[$i] = $this->data[$i];
        }
        return $ar;
    }

    public function printVector() {
        echo __CLASS__ . PHP_EOL;
        for ($j = 0; $j < $this->col; ++$j) {
            printf('%lf  ', $this->data[$j]);
        }
        echo PHP_EOL;
    }

    public function __toString() {
        return (string) $this->printVector();
    }

    protected function __construct(public int $col, int $dtype = self::DOUBLE) {
        if ($this->col < 1) {
            throw new invalidArgumentException('* To create Numphp/Vector col must be greater than 0!, Op Failed! * ');
        }
        parent::__construct($this->col, $dtype);
        return $this;
    }

}
