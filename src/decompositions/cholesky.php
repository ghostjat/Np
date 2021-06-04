<?php

declare (strict_types=1);

namespace Np\decompositions;

use Np\matrix;
use Np\core\lapack;

/**
 * Cholesky
 *
 * An efficient decomposition of a square positive definite matrix into a
 * lower triangular matrix and its conjugate transpose.
 * 
 * @package NumPhp
 * @category Scientific Library
 * @author ghost (Shubham Chaudhary)
 * @email ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */
class cholesky {

    /**
     * 
     * @param matrix $m
     * @return matrix|null
     * 
     */
    public static function factory(matrix $m): matrix|null {
        if ($m->isSquare()) {
            $ar = $m->copyMatrix();
            $lp = lapack::potrf($ar);
            if ($lp != 0) {
                return null;
            }
            for ($i = 0; $i < $m->col; ++$i) {
                for ($j = $i + 1; $j < $m->col; ++$j) {
                    $ar->data[$i * $m->col + $j] = 0.0;
                }
            }
            unset($lp);
            return $ar;
        }
    }

}
