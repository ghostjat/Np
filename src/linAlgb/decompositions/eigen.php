<?php

declare (strict_types=1);

namespace Np\linAlgb\decompositions;

use Np\matrix;
use Np\vector;
use Np\core\lapack;
use Np\exceptions\invalidArgumentException;

/**
 * Eigen
 *
 * The Eigen decompositon or (Spectral decomposition) is a matrix factorization resulting in a matrix of eigenvectors and a
 * corresponding array of eigenvalues.
 * 
 * @package Np
 * @category Scientific Library
 * @author ghost (Shubham Chaudhary)
 * @email ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */
class eigen {

    protected $eignVal;
    protected $eignVec;

    public static function factory(\Np\matrix $m, bool $symmetric = false): self {
        if (!$m->isSquare()) {
            throw new invalidArgumentException('A Non Square Matrix is given!');
        }
        $wr = vector::factory($m->col);
        $ar = $m->copy();
        if ($symmetric) {
            $lp = lapack::syev($ar, $wr);
            if ($lp != 0) {
                return null;
            }

            return new self($wr, $ar);
        } else {
            $wi = vector::factory($m->col);
            $vr = matrix::factory($m->col, $m->col);

            $lp = lapack::geev($ar, $wr, $wi, $vr);
            if ($lp != 0) {
                return null;
            }

            return new self($wr, $vr);
        }
    }

    /**
     * 
     * @param vector $eignVal
     * @param matrix $eignVec
     */
    protected function __construct(vector $eignVal, matrix $eignVec) {
        $this->eignVal = $eignVal;
        $this->eignVec = $eignVec;
    }

    /**
     * Return the eigenvalues of the eigen decomposition.
     * 
     * @return vector
     */
    public function eigenVal(): vector {
        return $this->eignVal;
    }

    /**
     * Return the eigen vectors of the eigen decomposition.
     * 
     * @return matrix
     */
    public function eigenVec(): matrix {
        return $this->eignVec->transpose();
    }

}
