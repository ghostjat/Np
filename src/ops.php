<?php

namespace Np;

/**
 * Operations
 * A fast lite & memory efficient Scientific Computing for php
 * 
 * @package   NP
 * @category  Scientific Computing
 * @author    ghost (Shubham Chaudhary)
 * @email     ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */
trait ops {

    /**
     * Return the element-wise maximum of two matrices| two vectors.
     * @param matrix|vector $d
     * @return matrix|vector
     */
    public function max(matrix|vector $d): matrix|vector {
        if ($this instanceof matrix && $d instanceof matrix && $this->checkShape($this, $d) && $this->checkDtype($this, $d)) {
            $r = self::factory($this->row, $this->col, $this->dtype);
        } elseif ($this instanceof vector && $d instanceof vector && $this->checkShape($this, $d) && $this->checkDtype($this, $d)) {
            $r = self::factory($this->col, $this->dtype);
        }
        for ($i = 0; $i < $this->ndim; ++$i) {
            $r->data[$i] = max($this->data[$i], $d->data[$i]);
        }
        return $r;
    }

    /**
     * Return the element-wise minimum of two matrices|two vectors.
     * @param matrix|vector $d
     * @return matrix|vector
     */
    public function min(matrix|vector $d): matrix|vector {
        if ($this instanceof matrix && $d instanceof matrix && $this->checkShape($this, $d) && $this->checkDtype($this, $d)) {
            $r = self::factory($this->row, $this->col, $this->dtype);
        } elseif ($this instanceof vector && $d instanceof vector && $this->checkShape($this, $d) && $this->checkDtype($this, $d)) {
            $r = self::factory($this->col, $this->dtype);
        }
        for ($i = 0; $i < $this->ndim; ++$i) {
            $r->data[$i] = max($this->data[$i], $d->data[$i]);
        }
        return $r;
    }

    /**
     * Run a function over all of the elements in the matrix|vector. 
     * @param callable $func
     * @return matrix|vector
     */
    public function map(callable $func): matrix|vector {
        if ($this instanceof matrix) {
            $r = self::factory($this->row, $this->col, $this->dtype);
        } else {
            $r = self::factory($this->col, $this->dtype);
        }
        for ($i = 0; $i < $this->ndim; ++$i) {
            $r->data[$i] = $func($this->data[$i]);
        }
        return $r;
    }

    public function abs(): matrix|vector {
        return $this->map('abs');
    }

    public function sqrt(): matrix|vector {
        return $this->map('sqrt');
    }

    public function exp(): matrix|vector {
        return $this->map('exp');
    }

    public function exp1(): matrix|vector {
        return $this->map('exp1');
    }

    public function log(float $b = M_E): matrix|vector {
        $ar = $this->copy();
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
     * make a copy of matrix|vector;
     * @return matrix|vector
     */
    public function copy(): matrix|vector {
        return clone $this;
    }

    /**
     * Clip the elements in the matrix to be between given minimum and maximum
     * and return a new matrix.
     * 
     * @param float $min
     * @param float $max
     * @return matrix
     */
    public function clip(float $min, float $max): matrix|vector {
        if ($this instanceof matrix) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
        } else {
            $ar = self::factory($this->col, $this->dtype);
        }
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
     * Clip the matrix|vector to be lower bounded by a given minimum.
     * @param float $min
     * @return matrix
     */
    public function clipLower(float $min): matrix {
        if ($this instanceof matrix) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
        } else {
            $ar = self::factory($this->col, $this->dtype);
        }
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
     * Clip the matrix|vector to be upper bounded by a given maximum.
     *
     * @param float $max
     * @return matrix
     */
    public function clipUpper(float $max): matrix|vector {
        if ($this instanceof matrix) {
            $ar = self::factory($this->row, $this->col, $this->dtype);
        } else {
            $ar = self::factory($this->col, $this->dtype);
        }
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
     * return a reshaped data buffer as matrix
     * @param int $row
     * @param int $col
     * @return matrix
     */
    public function reshape(int $row, int $col): matrix {
        if ($this->ndim != $row * $col) {
            self::_dimensionaMisMatchErr('given dimenssion is not valid for current bufferData');
        }
        if ($this instanceof vector) {
            $ar = matrix::factory($row, $col, $this->dtype);
            $ar->data = $this->data;
            return $ar;
        }
        $this->row = $row;
        $this->col = $col;
        return $this;
    }

}
