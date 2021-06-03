<?php

declare (strict_types=1);

namespace Np\decompositions;

use Np\matrix;
use Np\vector;
use Np\core\lapack;

/**
 * SVD
 * Compute the singular value decomposition of a matrix and 
 * return an object of the singular values and unitary matrices
 * 
 * @package NumPhp
 * @category Scientific Library
 * @author ghost (Shubham Chaudhary)
 * @email ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */
class svd {

    protected $u, $s, $v;

    /**
     * 
     * @param matrix $m
     * @return self
     */
    public static function factory(\Np\matrix $m): self {
        $k = min($m->row, $m->col);
        $ar = $m->copyMatrix();
        $s = vector::factory($k, $m->dtype);
        $u = matrix::factory($m->row, $m->row, $m->dtype);
        $v = matrix::factory($m->col, $m->col, $m->dtype);
        $lp = lapack::gesdd($ar, $s, $u, $v);
        if ($lp != 0) {
            return null;
        }
        unset($ar);
        unset($k);
        unset($lp);
        return new self($u, $v, $s);
    }

    /**
     * 
     * @param \Np\matrix $u
     * @param \Np\matrix $v
     * @param \Np\vector $s
     */
    protected function __construct(\Np\matrix $u, \Np\matrix $v, \Np\vector $s) {
        $this->u = $u;
        $this->s = $s;
        $this->v = $v;
    }

    /**
     * Return the U matrix.
     * 
     * @return matrix
     */
    public function u(): matrix {
        return $this->u;
    }

    /**
     * Return the singular values of matrix.
     * 
     * @return vector
     */
    public function s(): vector {
        return $this->s;
    }

    /**
     * Return the V matrix.
     * 
     * @return matrix
     */
    public function v(): matrix {
        return $this->v;
    }

    /**
     * Return the V transposed matrix.
     *
     * @return matrix
     */
    public function vt(): matrix {
        return $this->v->transpose();
    }

}
