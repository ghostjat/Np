<?php
declare(strict_types=1);

namespace numphp;

/** A fast lite memory efficient Scientific Computing in php
 * Vector (rank-1)
 * 
 * @package NumPhp
 * @version 0.0.1.alpha
 * @category Php Scientific Library
 * @author ghost (Shubham Chaudhary)
 * @email ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 * 
 */
class vector {
    const INT=0, FLOAT = 1, DOUBLE = 2;
    public $data,$col,$dtype;

    /**
     * 
     * @param int $col
     * @param int $dtype
     * @return vector
     */
    public static function factory(int $col, int $dtype = self::FLOAT) : vector {
        return new self($col, $dtype);
    }
    
    /**
     * create vector using php array
     * @param array $data
     * @param int $dtype
     * @return vector
     */
    public static function ar(array $data,int $dtype= self::FLOAT): vector {
        if (is_array($data) && !is_array($data[0])) {
            $ar = self::factory(count($data),$dtype);
            $ar->setData($data);
        } else {
            self::_err('data must be of same dimensions');
        }
        return $ar;
    }
    
    
    /**
     * Return vector with random values
     * @param int $col
     * @param int $dtype
     * @return vector
     */
    public static function randn(int $col, int $dtype= self::FLOAT): vector {
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
     * @return vector
     */
    public static function full(int $col, $val, int $dtype = self::FLOAT): vector {
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
     * @return vector
     */
    public static function range($start, $end, $interval = 1) : vector {
        return self::ar(range($start, $end, $interval));
    }
    
    /**
     * Return a Gaussian random vector with mean 0
     * and unit variance.
     *
     * @param int $n
     * @return self
     */
    public static function gaussian(int $n) : vector {

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
    public static function poisson(int $n, float $lambda = 1.) : vector
    {
        $max = getrandmax();

        $l = exp(-$lambda);

        $a = [];

        while (count($a) < $n) {
            $k = 0;
            $p = 1.0;

            while ($p > $l) {
                ++$k;

                $p *= rand() / $max;
            }

            $a[] = $k - 1;
        }

        return self::ar($a);
    }
    
    /**
     * Return a vector of n evenly spaced numbers between minimum and maximum.
     *
     * @param float $min
     * @param float $max
     * @param int $n
     * @throws \InvalidArgumentException
     * @return vector
     */
    public static function linspace(float $min, float $max, int $n) : vector {
        if ($min > $max) {
            self::_invalidArgument('Minimum must be less than maximum.');
        }
        if ($n < 2) {
            self::_invalidArgument('Number of elements must be greater than 1.');
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
     * make a copy of vector
     * @return vector
     */
    public function copyVector() : vector {
        return clone $this;
    }


    /**
     * vector-vector dot product
     * @param \numphp\vector $vector
     * @param int $incX
     * @param int $incY
     * @return vector
     */
    public function dotVector(\numphp\vector $vector ) {
        if($this->checkDtype($vector)) {
            if($this->dtype == self::FLOAT) {
                return core\blas::sdot($this, $vector);
            }
            if($this->dtype == self::DOUBLE ) {
                return core\blas::ddot($this, $vector);
            }
        }
    }
    
    public function sum():float {
        $r = 0.0;
        for($i = 0; $i < $this->col; ++$i) {
            $r += $this->data[$i];
        }
        return $r;
    }
    
    /**
     * Return the product of the vector.
     * @return int|float
     */
    public function product():float {
        $r = 1.0;
        for($i = 0; $i < $this->col; ++$i) {
            $r *= $this->data[$i];
        }
        return $r;
    }
    
    /**
     * Compute the vector-matrix dot product of this vector and matrix .
     * @param \numphp\matrix $matrix
     * @return vector
     */
    public function dotMatrix(\numphp\matrix $matrix): vector {
        if($this->dtype != $matrix->dtype) {
            self::_err('Mismatch Dtype of given matrix');
        }
        $mvr = self::factory($this->col, $this->dtype);
        if($this->dtype == self::FLOAT) {
            core\blas::sgemv($matrix, $this, $mvr);
        }
        else {
            core\blas::dgemv($matrix, $this, $mvr);
        }
        return $mvr;
    }
    
    /**
     * 
     * @param \numphp\matrix $matrix
     * @return matrix
     */
    public function addMatrix(\numphp\matrix $matrix):matrix {
        if($this->col == $matrix->col && $this->dtype == $matrix->dtype) {
            $vr = matrix::factory($matrix->row,$matrix->col, $matrix->dtype);
            for($i = 0; $i < $matrix->row; ++$i) {
                for($j = 0; $j < $matrix->col; ++$j) {
                    $vr->data[$i * $matrix->col +$j] = $this->data[$j] + $matrix->data[$i * $matrix->col + $j];
                }
            }
            return $vr;
        }
        self::_invalidArgument('');
    }
    
    /**
     * 
     * @param \numphp\matrix $matrix
     * @return matrix
     */
    public function multiplyMatrix(\numphp\matrix $matrix):matrix {
        if($this->col == $matrix->col && $this->dtype == $matrix->dtype) {
            $vr = matrix::factory($matrix->row,$matrix->col, $matrix->dtype);;
            for($i = 0; $i < $matrix->row; ++$i) {
                for($j = 0; $j < $matrix->col; ++$j) {
                    $vr->data[$i * $matrix->col +$j] = $this->data[$j] * $matrix->data[$i * $matrix->col + $j];
                }
            }
            return $vr;
        }
        self::_invalidArgument('');
    }
    
    /**
     * 
     * @param int|float|matrix|vector $d
     * @return matrix|vector
     */
    public function divide(int|float|matrix|vector $d):matrix|vector {
        if($d instanceof matrix){
            return $this->divideMatrix($d);
        } elseif ($d instanceof self) {
            return $this->divideVector($d);
        } else {
            return $this->divideScalar($d);
        }
    }
    
    /**
     * 
     * @param \numphp\matrix $matrix
     * @return matrix
     */
    protected function divideMatrix(\numphp\matrix $matrix):matrix {
        if($this->col == $matrix->col && $this->dtype == $matrix->dtype) {
            $vr = matrix::factory($matrix->row,$matrix->col, $matrix->dtype);
            for($i = 0; $i < $matrix->row; ++$i) {
                for($j = 0; $j < $matrix->col; ++$j) {
                    $vr->data[$i * $matrix->col +$j] = $this->data[$j] / $matrix->data[$i * $matrix->col + $j];
                }
            }
            return $vr;
        }
        self::_invalidArgument('');
    }
    
    /**
     * 
     * @param vector $v
     * @return vector
     */
    protected function divideVector(vector $v) :vector {
        if($this->checkDimensions($v) && $this->checkDtype($v)) {
            $vr = self::factory($this->col, $this->dtype);
            for($i = 0; $i < $this->col; ++$i) {
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
     * @param \numphp\vector $vector
     * @return vector
     */
    public function multiplyVector(\numphp\vector $vector):vector {
        if($this->checkDimensions($vector) && $this->checkDtype($vector)) {
            $vr = self::factory($this->col, $this->dtype);
            for($i = 0; $i < $this->col; ++$i) {
                $vr->data[$i] = $this->data[$i] * $vector->data[$i];
            }
            return $vr;
        }
    }
    
    /**
     * 
     * @param \numphp\vector $vector
     * @return vector
     */
    public function sumVector(\numphp\vector $vector):vector {
        if($this->checkDimensions($vector) && $this->checkDtype($vector)) {
            $vr = self::factory($this->col, $this->dtype);
            for($i = 0; $i < $this->col; ++$i) {
                $vr->data[$i] = $this->data[$i] + $vector->data[$i];
            }
            return $vr;
        }
    }
    
    
    
    /**
     * 
     * @param \numphp\vector $vector
     * @return vector
     */
    public function powVector(\numphp\vector $vector):vector {
        if($this->checkDimensions($vector) && $this->checkDtype($vector)) {
            $vr = self::factory($this->col, $this->dtype);
            for($i = 0; $i < $this->col; ++$i) {
                $vr->data[$i] = $this->data[$i] ** $vector->data[$i];
            }
            return $vr;
        }
    }
    
    /**
     * 
     * @param \numphp\vector $vector
     * @return vector
     */
    public function modVector(\numphp\vector $vector):vector {
        if($this->checkDimensions($vector) && $this->checkDtype($vector)) {
            $vr = self::factory($this->col, $this->dtype);
            for($i = 0; $i < $this->col; ++$i) {
                $vr->data[$i] = $this->data[$i] % $vector->data[$i];
            }
            return $vr;
        }
    }
    
    /**
     * 
     * @param \numphp\vector $vector
     * @return vector
     */
    public function subtractVector(\numphp\vector $vector):vector {
        if($this->checkDimensions($vector) && $this->checkDtype($vector)) {
            $vr = self::factory($this->col, $this->dtype);
            for($i = 0; $i < $this->col; ++$i) {
                $vr->data[$i] = $this->data[$i] - $vector->data[$i];
            }
            return $vr;
        }
    }
    
    /**
     * 
     * @param \numphp\vector $vector
     * @return \numphp\vector
     */
    public function substractScalar(int|float $scalar): vector {
        $vr = self::factory($this->col, $this->dtype);
        for ($i = 0; $i < $this->col; ++$i) {
            $vr->data[$i] = $this->data[$i] - $scalar;
        }
        return $vr;
    }

    public function convolve(\numphp\vector $v, int $stride = 1): vector {
        return convolve::conv1D($this, $v, $stride);
    }

    /**
     * Return the inner product of two vectors.
     *
     * @param \numphp\vector $vector
     * @return float
     */
    public function inner(\numphp\vector $vector) {
        return $this->dotVector($vector);
    }

    public function l1_norm() {
        
    }
    
    public function l2_norm() {
        
    }

    /**
     *  Return the element-wise maximum of given vector with current vector
     * @param \numphp\vector $vector
     * @return vector
     */
    public function maximum(\numphp\vector $vector) :vector {
        if($this->checkDimensions($vector)) {
            
        }
    }
    
    /**
     *  Return the element-wise minium of given vector with current vector
     * @param \numphp\vector $vector
     * @return vector
     */
    public function minium(\numphp\vector $vector) :vector {
        if($this->checkDimensions($vector)) {
            
        }
    }


    /**
     * set data to vector
     * @param type $data
     */
    public function setData($data) {
        if (!is_null($data)) {
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
    }
    
    public function asMatrix():matrix {
        $ar = matrix::factory(1, $this->col, $this->dtype);
        for($i = 0; $i < 1; ++$i) {
            for($j = 0; $j < $this->col; ++$j) {
                $ar->data[$i * $this->col + $j] = $this->data[$j];
            }
        }
        return $ar;
    }
    
    /**
     * get the shape of matrix
     * @return int
     */
    public function getShape() :int {
        return  $this->col;
    }
    
    public function getDtype() {
        return $this->dtype;
    }
    
    public function asArray() {
        $ar = array_fill(0, $this->col, null);
        for($i = 0; $i < $this->col; ++$i) {
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
    
    protected function checkDimensions(vector $vector) {
        if($this->col != $vector->col) {
            self::_err('Mismatch Dimensions of given vector');
        }
        return true;
    }
    
     protected function checkDtype(vector $vector) {
        if($this->dtype != $vector->dtype) {
            self::_err('Mismatch Dtype of given vector');
        }
        return true;
    }
    
    protected function __construct(int $col, int $dtype = self::FLOAT) {
        if($col < 1 ) {
            throw new InvalidArgumentException('* To create Numphp/Vector col must be greater than 0!, Op Failed! * ');
        }
        $this->col = $col;
        $this->dtype = $dtype;
        switch ($dtype) {
            case self::FLOAT:
                $this->data = self::c_FloatVector($this->col);
                break;
            case self::DOUBLE:
                $this->data = self::c_DoubleVector($this->col);
                break;
            case self::INT:
                $this->data = self::c_IntVector($this->col);
                break;
            default :
                self::_invalidArgument('given dtype is not supported by numphp');
                break;
        }
        return $this;
    }
    
    private static function c_FloatVector(int $col): \FFI\CData {
        return \FFI::cast('float *', \FFI::new("float[$col]"));
    }
    
    private static function c_IntVector(int $col): \FFI\CData {
        return \FFI::cast('int *', \FFI::new("int[$col]"));
    }
    
    private static function c_DoubleVector(int $col): \FFI\CData {
        return \FFI::cast('double *', \FFI::new("double[$col]"));
    }
    
    private static function _err($msg): \Exception {
        throw new \Exception($msg);
    }
    
    private static function _invalidArgument($argument) : \InvalidArgumentException{
        throw new \InvalidArgumentException($argument);
    }
    
}
