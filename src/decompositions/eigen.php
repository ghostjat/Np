<?php

declare (strict_types=1);

namespace Np\decompositions;

use Np\matrix;
use Np\vector;
use Np\core\lapack;
use InvalidArgumentException;

/**
 * Eigen
 *
 * The Eigen decompositon or (Spectral decomposition) is a matrix factorization resulting in a matrix of eigenvectors and a
 * corresponding array of eigenvalues.
 * 
 * @package NumPhp
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
            throw new InvalidArgumentException('A Non Square Matrix is given!');
        }
        $wr = vector::factory($m->col, $m->dtype);
        $ar = $m->copyMatrix();
        if ($symmetric) {
            $lp = lapack::syev($ar, $wr);
            if ($lp != 0) {
                return null;
            }

            return new self($wr, $ar);
        } else {
            $wi = vector::factory($m->col, $m->dtype);
            $vr = matrix::factory($m->col, $m->col, $m->dtype);

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
