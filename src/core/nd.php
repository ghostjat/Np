<?php

declare(strict_types=1);

namespace Np\core;

/**
 * ND
 * A fast lite memory efficient Scientific Computing for php
 * 
 * @package   NumPhp
 * @category  Scientific Computing
 * @author    ghost (Shubham Chaudhary)
 * @email     ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 * 
 */
class nd {

    const TWO_PI = 2. * M_PI, EPSILON = 1e-8;
    const FLOAT = 1, DOUBLE = 2, INT = 3;

    public $data, $ndim, $dtype;
    public static $_time = null, $_mem = null;

    protected function __construct(int $size, int $dtype = self::FLOAT) {
        $this->ndim = $size;
        $this->dtype = $dtype;
        switch ($dtype) {
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
                self::_invalidArgument('given dtype is not supported by Np');
                break;
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

    protected static function _err($msg): \Exception {
        throw new \Exception($msg);
    }

    protected static function _invalidArgument($argument): \InvalidArgumentException {
        throw new \InvalidArgumentException($argument);
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
