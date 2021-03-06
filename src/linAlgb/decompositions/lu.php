<?php

declare (strict_types=1);

namespace Np\linAlgb\decompositions;

use Np\matrix;
use Np\vector;
use Np\core\lapack;
use Np\exceptions\invalidArgumentException;

/**
 * LU
 *
 * The LU decomposition is a factorization of a Matrix as the product of a
 * lower and upper triangular matrix as well as a permutation matrix.
 *
 * @package Np
 * @category Scientific Library
 * @author ghost (Shubham Chaudhary)
 * @email ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */
class lu {

    /**
     * 
     * @param matrix $m
     * @return self
     * @throws InvalidArgumentException
     */
    public static function factory(matrix $m): self {
        if (!$m->isSquare()) {
            throw new invalidArgumentException('Matrix must be given.');
        }
        $ipiv = vector::factory($m->col, vector::INT);
        $ar = $m->copy();
        $lp = lapack::getrf($ar, $ipiv);
        if ($lp != 0) {
            return null;
        }
        $l = matrix::factory($m->col, $m->col);
        $u = matrix::factory($m->col, $m->col);
        $p = matrix::factory($m->col, $m->col);
        for ($i = 0; $i < $m->col; ++$i) {
            for ($j = 0; $j < $i; ++$j) {
                $l->data[$i * $m->col + $j] = $ar->data[$i * $m->col + $j];
            }
            $l->data[$i * $m->col + $i] = 1.0;
            for ($j = $i + 1; $j < $m->col; ++$j) {
                $l->data[$i * $m->col + $j] = 0.0;
            }
        }
        for ($i = 0; $i < $m->col; ++$i) {
            for ($j = 0; $j < $i; ++$j) {
                $u->data[$i * $m->col + $j] = 0.0;
            }
            for ($j = $i; $j < $m->col; ++$j) {
                $u->data[$i * $m->col + $j] = $ar->data[$i * $m->col + $j];
            }
        }
        for ($i = 0; $i < $m->col; ++$i) {
            for ($j = 0; $j < $m->col; ++$j) {
                if ($j == $ipiv->data[$i] - 1) {
                    $p->data[$i * $m->col + $j] = 1;
                } else {
                    $p->data[$i * $m->col + $j] = 0;
                }
            }
        }
        unset($ar);
        unset($ipiv);
        return new self($l, $u, $p);
    }

    /**
     * 
     * @param matrix $l
     * @param matrix $u
     * @param matrix $p
     */
    protected function __construct(protected matrix $l, protected matrix $u, protected matrix $p) {
        
    }

    /**
     * Return the lower triangular matrix.
     *
     * @return matrix
     */
    public function l(): matrix {
        return $this->l;
    }

    /**
     * Return the upper triangular matrix.
     *
     * @return matrix
     */
    public function u(): matrix {
        return $this->u;
    }

    /**
     * Return the permutation matrix.
     * 
     * @return matrix
     */
    public function p(): matrix {
        return $this->p;
    }

}
