<?php
declare(strict_types=1);
namespace Np\linAlgb;

use Np\core\{blas,lapack};
use Np\matrix, Np\vector;
/**
 * Linear Algebra
 * 
 * 
 * 
 * @package   Np
 * @category  Scientific Computing
 * @author    ghost (Shubham Chaudhary)
 * @email     ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */
trait linAlg {
    
    /**
     *  
     * get dot product of m.m | m.v | v.v
     * 
     * @param \Np\matrix|\Np\vector $d
     * @return matrix|vector
     */
    public function dot(matrix|vector $d): matrix|vector {
        if ($this instanceof matrix) {
            if ($d instanceof matrix) {
                return $this->dotMatrix($d);
            }
            return $this->dotVector($d);
        }
        return blas::dot($this, $d);
    }

    /**
     * get matrix & matrix dot product
     * @param \Np\matrix $matrix
     * @return \Np\matrix
     */
    protected function dotMatrix(matrix $matrix): matrix {
        if ($this->checkDimensions($this,$matrix)) {
            $ar = self::factory($this->row, $matrix->col);
            blas::gemm($this, $matrix, $ar);
            return $ar;
        }
    }

    /**
     * get dot product of matrix & a vector
     * @param \Np\vector $vector
     * @return \Np\vector
     */
    protected function dotVector(vector $vector): vector {
        if ($this->checkDimensions($vector, $this)) {
            $mvr = vector::factory($this->col);
            blas::gemv($this, $vector, $mvr);
            return $mvr;
        }
    }
    
    /**
     * 
     * Compute the multiplicative inverse of the matrix.
     * @return matrix|null
     */
    public function inverse(): matrix|null {
        if ($this->isSquare()) {
            $imat = $this->copy();
            $ipiv = vector::factory($this->row, vector::INT);
            $lp = lapack::getrf($imat, $ipiv);
            if ($lp != 0) {
                return null;
            }
            $lp = lapack::getri($imat, $ipiv);
            if ($lp != 0) {
                return null;
            }
            unset($ipiv);
            unset($lp);
            return $imat;
        }
        self::_err('Error::invalid Size of matrix!');
    }

    /**
     * FIXEME:-Bug noticed on 10/06/21
     * Compute the (Moore-Penrose) pseudo inverse of the general matrix.
     * @return matrix|null
     */
    public function pseudoInverse(): matrix|null {
        $k = min($this->row, $this->col);
        $s = vector::factory($k);
        $u = self::factory($this->row, $this->row);
        $vt = self::factory($this->col, $this->col);
        $imat = $this->copy();
        $lp = lapack::gesdd($imat, $s, $u, $vt);
        if ($lp != 0) {
            return null;
        }
        for ($i = 0; $i < $k; ++$i) {
            blas::scale(1.0 / $s->data[$i], $vt->rowAsVector($i));
        }
        unset($imat);
        unset($k);
        unset($lp);
        unset($s);
        $mr = self::factory($this->col, $this->row, $this->dtype);
        blas::gemm($vt, $u, $mr);
        unset($u);
        unset($vt);
        return $mr;
    }
    
}
