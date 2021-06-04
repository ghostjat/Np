<?php

declare(strict_types=1);

namespace Np;

use Np\core\{
    nd,
    blas,
    lapack
};
use Np\exceptions\{
    invalidArgumentException,
};

/** A fast lite memory efficient Scientific Computing in php
 * Vector (rank-1)
 * 
 * @package NumPhp
 * @version V0.0.alpha
 * @category Php Scientific Library
 * @author ghost (Shubham Chaudhary)
 * @email ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 * 
 */
class vector extends nd {

    /**
     * Factory method to build a new vector.
     * 
     * @param int $col
     * @param int $dtype
     * @return vector
     */
    public static function factory(int $col, int $dtype = self::FLOAT): vector {
        return new self($col, $dtype);
    }

    /**
     * Build a new vector from a php array.
     * 
     * @param array $data
     * @param int $dtype
     * @return vector
     */
    public static function ar(array $data, int $dtype = self::FLOAT): vector {
        if (is_array($data) && !is_array($data[0])) {
            $ar = self::factory(count($data), $dtype);
            $ar->setData($data);
            return $ar;
        } else {
            self::_err('data must be of same dimensions');
        }
    }

    /**
     * Return vector with random values
     * @param int $col
     * @param int $dtype
     * @return vector
     */
    public static function randn(int $col, int $dtype = self::FLOAT): vector {
        $ar = self::factory($col, $dtype);
        $max = getrandmax();
        for ($i = 0; $i < $ar->col; ++$i) {
            $ar->data[$i] = rand() / $max;
        }
        return $ar;
    }

    /**
     * Return vector with uniform values
     * @param int $col
     * @param int $dtype
     * @return vector
     */
    public static function uniform(int $col, int $dtype = self::FLOAT): vector {
        $ar = self::factory($col, $dtype);
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
     * @param int $dtype
     * @return vector
     */
    public static function zeros(int $col, int $dtype = self::FLOAT): vector {
        $ar = self::factory($col, $dtype);
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
    public static function ones(int $col, int $dtype = self::FLOAT): vector {
        $ar = self::factory($col, $dtype);
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
    public static function null(int $col, int $dtype = self::FLOAT): vector {
        $ar = self::factory($col, $dtype);
        for ($i = 0; $i < $col; ++$i) {
            $ar->data[$i] = null;
        }
        return $ar;
    }

    /**
     * create a vector with given scalar value
     * @param int $col
     * @param int|float|double $val
     * @param int $dtype
     * @return vector
     */
    public static function full(int $col, int|float $val, int $dtype = self::FLOAT): vector {
        $ar = self::factory($col, $dtype);
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
     * @param int $dtype 
     * @return vector
     */
    public static function range(int|float $start, int|float $end, int|float $interval = 1, int $dtype = self::FLOAT): vector {
        return self::ar(range($start, $end, $interval), $dtype);
    }

    /**
     * Return a Gaussian random vector with mean 0
     * and unit variance.
     *
     * @param int $n
     * @param int $dtype
     * @return self
     */
    public static function gaussian(int $n, int $dtype = self::FLOAT): vector {
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
        return self::ar($a, $dtype);
    }

    /**
     * Generate a vector with n elements from a Poisson distribution.
     *
     * @param int $n
     * @param float $lambda
     * @param int $dtype 
     * @return vector
     */
    public static function poisson(int $n, float $lambda = 1.0, int $dtype = self::FLOAT): vector {
        $max = getrandmax();
        $l = exp(-$lambda);
        $a = new self($n, $dtype);
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
     * @param int $dtype
     * @throws invalidArgumentException
     * @return vector
     */
    public static function linspace(float $min, float $max, int $n, int $dtype = self::FLOAT): vector {
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
        return self::ar($a, $dtype);
    }

    /**
     * make a copy of vector
     * @return vector
     */
    public function copyVector(): vector {
        return clone $this;
    }

    /**
     * Return the element-wise maximum of given vector with current vector
     * 
     * @param \Np\vector $vector
     * @return vector
     */
    public function maximum(\Np\vector $vector): vector {
        if ($this->checkShape($this, $vector) && $this->checkDtype($this, $vector)) {
            $v = new self($this->ndim, $this->dtype);
            for ($i = 0; $i < $v->ndim; ++$i) {
                $v->data[$i] = max($this->data[$i], $vector->data[$i]);
            }
            return $v;
        }
    }

    /**
     * Return the element-wise minium of given vector with current vector
     * 
     * @param \Np\vector $vector
     * @return vector
     */
    public function minium(\Np\vector $vector): vector {
        if ($this->checkShape($this, $vector) && $this->checkDtype($this, $vector)) {
            $v = new self($this->ndim, $this->dtype);
            for ($i = 0; $i < $v->ndim; ++$i) {
                $v->data[$i] = min($this->data[$i], $vector->data[$i]);
            }
            return $v;
        }
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
    public function argMx(): int {
        return blas::max($this);
    }

    /**
     * vector-vector dot product
     * @param \Np\vector $vector
     * @param int $incX
     * @param int $incY
     * @return vector
     */
    public function dotVector(\Np\vector $v) {
        if ($this->checkDtype($this, $v)) {
            return blas::dot($this, $v);
        }
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
            $mvr = self::factory($this->col, $this->dtype);
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
        } elseif ($d instanceof self) {
            return $this->divideVector($d);
        } else {
            return $this->divideScalar($d);
        }
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function divideMatrix(\Np\matrix $m): matrix {
        if ($this->checkShape($this, $m) && $this->checkDtype($this, $m)) {
            $vr = matrix::factory($m->row, $m->col, $m->dtype);
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
        if ($this->checkShape($this, $v) && $this->checkDtype($this, $v)) {
            $vr = self::factory($this->col, $this->dtype);
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
        $vr = self::factory($this->col, $this->dtype);
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
        } elseif ($d instanceof self) {
            return $this->multiplyVector($d);
        } else {
            return $this->multiplyScalar($d);
        }
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function multiplyMatrix(\Np\matrix $m): matrix {
        if ($this->checkShape($this, $m) && $this->checkDtype($this, $m)) {
            $vr = matrix::factory($m->row, $m->col, $m->dtype);
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
        if ($this->checkShape($this, $vector) && $this->checkDtype($this, $vector)) {
            $vr = self::factory($this->col, $this->dtype);
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
        $vr = $this->copyVector();
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
        } elseif ($d instanceof self) {
            return $this->addVector($d);
        } else {
            return $this->addScalar($d);
        }
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function addMatrix(\Np\matrix $m): matrix {
        if ($this->checkShape($this, $m) && $this->checkDtype($this, $m)) {
            $vr = matrix::factory($m->row, $m->col, $m->dtype);
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
        if ($this->checkShape($this, $vector) && $this->checkDtype($this, $vector)) {
            $vr = self::factory($this->col, $this->dtype);
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
        $vr = $this->copyVector();
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
        } elseif ($d instanceof vector) {
            return $this->powVector($d);
        } else {
            return $this->powScalar($d);
        }
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function powMatrix(\Np\matrix $m): matrix {
        if ($this->checkDimensions($this, $m) && $this->checkDtype($this, $m)) {
            $ar = matrix::factory($m->row, $m->col, $this->dtype);
            for ($i = 0; $i < $m->row; ++$i) {
                for ($j = 0; $j < $m - col; ++$j) {
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
        if ($this->checkShape($this, $vector) && $this->checkDtype($this, $vector)) {
            $vr = self::factory($this->col, $this->dtype);
            for ($i = 0; $i < $this->col; ++$i) {
                $vr->data[$i] = $this->data[$i] ** $vector->data[$i];
            }
            return $vr;
        }
    }

    /**
     * 
     * @param int|float $s
     */
    protected function powScalar(int|float $s) {
        $v = $this->copyVector();
        for ($i = 0; $i < $this->col; ++$i) {
            $v->data[$i] = $v->data[$i] ** $s;
        }
    }

    /**
     * 
     * @param int|float|\Np\matrix|\Np\vector $d
     * @return matrix|vector
     */
    public function mod(int|float|\Np\matrix|\Np\vector $d): matrix|vector {
        if ($d instanceof matrix) {
            return $this->powMatrix($d);
        } elseif ($d instanceof vector) {
            return $this->powVector($d);
        } else {
            return $this->powScalar($d);
        }
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function modMatrix(\Np\matrix $m): matrix {
        if ($this->checkDimensions($this, $m) && $this->checkDtype($this, $m)) {
            $ar = matrix::factory($m->row, $m->col, $this->dtype);
            for ($i = 0; $i < $m->row; ++$i) {
                for ($j = 0; $j < $m - col; ++$j) {
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
        if ($this->checkShape($this, $vector) && $this->checkDtype($this, $vector)) {
            $vr = self::factory($this->col, $this->dtype);
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
        $v = $this->copyVector();
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
        } elseif ($d instanceof self) {
            return $this->subtractVector($d);
        } else {
            return $this->substractScalar($d);
        }
    }

    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    protected function subtractMatrix(\Np\matrix $m): matrix {
        if ($this->checkShape($this, $m) && $this->checkDtype($this, $m)) {
            $vr = matrix::factory($m->row, $m->col, $m->dtype);
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
        if ($this->checkShape($this, $vector) && $this->checkDtype($this, $vector)) {
            $vr = self::factory($this->col, $this->dtype);
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
        $vr = self::factory($this->col, $this->dtype);
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

    /**
     * Run a function over all of the elements in the vector. 
     * 
     * @param callable $func
     * @return vector
     */
    public function map(callable $func): vector {
        $vr = self::factory($this->col, $this->dtype);
        for ($i = 0; $i < $this->col; ++$i) {
            $vr->data[$i] = $func($this->data[$i]);
        }
        return $vr;
    }

    public function log(float $b = M_E): vector {
        $vr = $this->copyVector();
        for ($i = 0; $i < $vr->col; ++$i) {
            log($vr->data[$i], $b);
        }
        return $vr;
    }

    public function max() {
        return $this->data[blas::max($this)];
    }

    public function min() {
        $this->data[blas::min($this)];
    }

    public function abs(): vector {
        return $this->map('abs');
    }

    public function sqrt(): vector {
        return $this->map('sqrt');
    }

    public function exp(): vector {
        return $this->map('exp');
    }

    public function exp1(): vector {
        return $this->map('exp1');
    }

    public function log1p(): vector {
        return $this->map('log1p');
    }

    public function sin(): vector {
        return $this->map('sin');
    }

    public function asin(): vector {
        return $this->map('asin');
    }

    public function cos(): vector {
        return $this->map('cos');
    }

    public function acos(): vector {
        return $this->map('acos');
    }

    public function tan(): vector {
        return $this->map('tan');
    }

    public function atan(): vector {
        return $this->map('atan');
    }

    public function radToDeg(): vector {
        return $this->map('rad2deg');
    }

    public function degToRad(): vector {
        return $this->map('deg2rad');
    }

    public function floor(): vector {
        return $this->map('floor');
    }

    public function ceil(): vector {
        return $this->map('ceil');
    }
    
    /**
     * 
     * @param float $min
     * @param float $max
     * @return vector
     */
    public function clip(float $min, float $max) : vector {
        if ($min > $max) {
            self::_invalidArgument('Minimum cannot be greater than maximum.');
        }

        $vr = self::factory($this->col, $this->dtype);
        
        for($i = 0; $i < $this->col; ++$i) {
            if ($this->data[$i] > $max) {
                $vr->data[$i] = $max;
                continue;
            }
            if ($this->data[$i] < $min) {
                $vr->data[$i] = $min;
                continue;
            }
        }

        return $vr;
    }
    
    public function clipUpper(float $min) : vector {
        $vr = self::factory($this->col, $this->dtype);
        for($i = 0; $i < $this->col; ++$i) {
            if ($this->data[$i] > $min) {
                $vr->data[$i] = $min;
                continue;
            }
        }
        return $vr;
    }
    
    public function clipLower(float $min) : vector {
        $vr = self::factory($this->col, $this->dtype);
        for($i = 0; $i < $this->col; ++$i) {
            if ($this->data[$i] < $min) {
                $vr->data[$i] = $min;
                continue;
            }
        }
        return $vr;
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
        return self::ones($this->col, $this->dtype)
                        ->divideVector($this);
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

        $a = $this->copyVector()->sort();
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

        $sd = $this->substractScalar($mean)
            ->square()
            ->sum();

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
        } elseif (is_numeric($data)) {
            for ($i = 0; $i < $this->col; ++$i) {
                $this->data[$i] = $data;
            }
        }
    }

    public function asMatrix(): matrix {
        $size = (int) sqrt($this->col);
        $ar = matrix::factory($size, $size, $this->dtype);
        for ($i = 0; $i < $ar->ndim; ++$i) {
            $ar->data[$i] = $this->data[$i];
        }
        return $ar;
    }

    /**
     * get the shape of matrix
     * @return int
     */
    public function getShape(): int {
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
        for ($j = 0; $j < $this->col; ++$j) {
            printf('%lf  ', $this->data[$j]);
        }
        echo PHP_EOL;
    }

    public function __toString() {
        return (string) $this->printVector();
    }

    protected function __construct(public int $col, int $dtype = self::FLOAT) {
        if ($this->col < 1) {
            throw new invalidArgumentException('* To create Numphp/Vector col must be greater than 0!, Op Failed! * ');
        }
        parent::__construct($this->col, $dtype);
        return $this;
    }

}
