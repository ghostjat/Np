<?php
declare(strict_types=1);
namespace Np\reductions;

use Np\matrix;

/**
 * RREF
 *
 * The reduced row echelon form (RREF) of a matrix.
 *
 * @category Scientific Computing
 * @author ghost (Shubham Chaudhary)
 * @email ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */

class rref { 
     
    /**
     * 
     * @param \Np\matrix $m
     * @return matrix
     */
    public static function factory(\Np\matrix $m): matrix {
        $lead = 0;
        $ar = $m->copyMatrix();
        for ($r = 0; $r < $ar->row; ++$r) {
            if ($lead >= $ar->col)
                break; {
                $i = $r;
                while ($ar->data[$i * $ar->col + $lead] == 0) {
                    $i++;
                    if ($i == $ar->row) {
                        $i = $r;
                        $lead++;
                        if ($lead == $ar->col) {
                            return $ar;
                        }
                    }
                }
                $ar->swapRows($r, $i);
            } {
                $lv = $ar->data[$r * $ar->col + $lead];
                for ($j = 0; $j < $ar->col; ++$j) {
                    $ar->data[$r * $ar->col + $j] = $ar->data[$r * $ar->col + $j] / $lv;
                }
            }
            for ($i = 0; $i < $ar->row; ++$i) {
                if ($i != $r) {
                    $lv = $ar->data[$i * $ar->col + $lead];
                    for ($j = 0; $j < $ar->col; ++$j) {
                        $ar->data[$i * $ar->col + $j] -= $lv * $ar->data[$r * $ar->col + $j];
                    }
                }
            }
            $lead++;
        }
        return $ar;
    }

}
