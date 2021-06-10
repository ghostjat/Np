<?php
declare(strict_types=1);
namespace Np\reductions;

use Np\matrix;
use Np\vector;
use Np\core\lapack;

/**
 * REF
 *
 * The row echelon form (REF) of a matrix.
 * 
 * @package Np
 * @category Scientific Computing
 * @author ghost (Shubham Chaudhary)
 * @email ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */

class ref { 
    
     /**
      * 
      * @param \Np\matrix $m
      * @return matrix|null
      */
    public static function factory(\Np\matrix $m): matrix|null {
        $ipiv = vector::factory(min($m->row, $m->col), vector::INT);
        $ar = $m->copy();
        $lp = lapack::getrf($ar, $ipiv);
        if ($lp != 0) {
            return null;
        }
        
        return $ar;
    }
}
