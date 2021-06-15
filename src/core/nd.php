<?php

declare(strict_types=1);

namespace Np\core;

use Np\exceptions\{
    dtypeException,
    invalidArgumentException,
    runtimeException,
    dimensionalityMismatch
};

/**
 * ND
 * 
 * A fast lite memory efficient Scientific Computing for php
 * 
 * @package   Np
 * @version   V0.0.1
 * @category  Scientific Computing
 * @author    ghost (Shubham Chaudhary)
 * @email     ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 * 
 */
class nd {

    const TWO_PI = 2. * M_PI, EPSILON = 1e-8;
    const FLOAT = 1, DOUBLE = 2, INT = 3;
    public $data;
    protected $_time = null, $_mem = null;
    
    public function checkDimensions(\Np\matrix|\Np\vector $Obj1, \Np\matrix $Obj2) {
        if($Obj1->col == $Obj2->row){
            return true;
        }
        self::_dimensionaMisMatchErr('Mismatch Dimensions of given Objects! Obj-A col & Obj-B row amount need to be the same!');
    }
    
    public function checkDtype(\Np\matrix|\Np\vector $Obj1, \Np\matrix|\Np\vector $Obj2){
        if($Obj1->dtype == $Obj2->dtype) {
            return true;
        }
        self::_dtypeErr('mismatch data type of given Np\Objects!');
    }
    
    public function checkShape(\Np\matrix|\Np\vector $Obj1, \Np\matrix|\Np\vector $Obj2) {
        if ($Obj1 instanceof \Np\vector && $Obj2 instanceof \Np\vector) {
            if ($Obj1->col == $Obj2->col) {
                return true;
            }
            self::_dimensionaMisMatchErr('mismatch Dimensions of given vectors!');
        } elseif ($Obj1 instanceof \Np\vector && $Obj2 instanceof \Np\matrix) {
            if ($Obj1->col == $Obj2->col) {
                return true;
            }
            self::_dimensionaMisMatchErr('mismatch Dimensions of given vectors & matrix!');
        } else {
            if ($Obj1->row == $Obj2->row || $Obj1->col == $Obj2->col) {
                return true;
            }
            self::_dimensionaMisMatchErr('mismatch Dimensions of given matrix!');
        }
    }
    
    public function asType(int $dtype){
        switch ($dtype){
            case self::FLOAT:
                \FFI::cast('float *', $this->data);
                break;
            case self::DOUBLE:
                \FFI::cast('double *', $this->data);
                break;
            case self::INT:
                \FFI::cast('int *', $this->data);
                break;
        }
    }

    protected function __construct(public int $ndim, public int $dtype = self::DOUBLE) {
        $this->getMemory();
        $this->time();
        $this->_nd();
    }
    
    protected function _nd() {
        switch ($this->dtype) {
            case self::FLOAT:
                $this->data = self::_ndFloat($this->ndim);
                break;
            case self::DOUBLE:
                $this->data = self::_ndDouble($this->ndim);
                break;
            case self::INT:
                $this->data = self::_ndInt($this->ndim);
                break;
            default :
                throw new dtypeException('given dtype is not supported by Np');
        }
    }

    protected static function _ndFloat(int $size) {
        return \FFI::cast('float *', \FFI::new("float[$size]"));
    }

    protected static function _ndDouble(int $size) {
        return \FFI::cast('double *', \FFI::new("double[$size]"));
    }

    protected static function _ndInt(int $size) {
        return \FFI::cast('int *', \FFI::new("int[$size]"));
    }

    protected static function _err($msg): runtimeException {
        throw new runtimeException($msg);
    }

    protected static function _invalidArgument($argument): invalidArgumentException {
        throw new invalidArgumentException($argument);
    }
    
    protected static function _dtypeErr($msg) : dtypeException {
        throw new dtypeException($msg);
    }
    
    protected static function _dimensionaMisMatchErr($msg) :dimensionalityMismatch {
        throw new dimensionalityMismatch($msg);
    }

    /**
     * set Timer, get total time 
     */
    public function time() {
        if (is_null($this->_time)) {
            $this->_time = microtime(true);
        } else {
            echo 'Time-Consumed:- ' . (microtime(true) - $this->_time) . PHP_EOL;
        }
    }

    /**
     * set memory dog, get total memory
     */
    public function getMemory() {
        if (is_null($this->_mem)) {
            $this->_mem = memory_get_usage();
        } else {
            $memory = memory_get_usage() - $this->_mem;
            $unit = ['b', 'kb', 'mb', 'gb', 'tb', 'pb'];
            echo round($memory / pow(1024, ($i = floor(log($memory, 1024)))), 2) . $unit[$i] . PHP_EOL;
        }
    }

}
