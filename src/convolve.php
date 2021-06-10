<?php
declare(strict_types=1);

namespace Np;
use Np\exceptions\dtypeException;

/**
 * Convolve
 * 
 * 1D & 2D SignalProcessing in pure php
 * 
 * @package Np
 * @category  Scientific Computing
 * @author    ghost (Shubham Chaudhary)
 * @email     ghost.jat@gmail.com
 * @copyright (c) 2020-2021, Shubham Chaudhary
 */
class convolve {
    
    /**
     * 1D convolution between a vector v and kernel k, with a given stride.
     * 
     * @param \Np\vector $v
     * @param \Np\vector $k
     * @param int $stride
     * @return vector
     * @throws \Exception
     */
    public static function conv1D(\Np\vector $v, \Np\vector $k, int $stride = 1): vector {
        if ($v->dtype == $k->dtype) {
            $nc = $v->col + $k->col - 1;
            $r = vector::factory($nc / $stride);
            for ($i = 0; $i < $nc; $i += $stride) {
                $jmin = $i >= $k->col - 1 ? $i - ($k->col - 1) : 0;
                $jmax = $i <= $v->col ? $i : $v->col - 1;
                $sigma = 0.0;
                for ($j = $jmin; $j <= $jmax; ++$j) {
                    $sigma += $v->data[$j] * $k->data[$i - $j];
                }
                $r->data[$i] = $sigma;
            }
            return $r;
        }
        else {
            throw new dtypeException('Err::given vectors has diffrent data type!');
        }
    }
    
    /**
     * 2D convolution between a matrix ma and kernel kb, with a given stride.
     * 
     * @param \Np\matrix $ma
     * @param \Np\matrix $kb
     * @param int $stride
     * @return matrix
     * @throws \Exception
     */
    public static function conv2D(\Np\matrix $ma, \Np\matrix $kb, int $stride = 1): matrix{
        if ($ma->dtype == $kb->dtype) {
            $p = $kb->row / 2;
            $q = $kb->col / 2;
            $rc = matrix::factory($ma->row / $stride, $ma->col / $stride, $ma->dtype);
            for ($i = 0; $i < $ma->row; $i += $stride) {

                for ($j = 0; $j < $ma->col; $j += $stride) {
                    $sgima = 0.0;
                    for ($k = 0; $k < $kb->row; ++$k) {
                        $x = $i + $p - $k;
                        if ($x < 0 || $x >= $ma->row) {
                            continue;
                        }
                        for ($l = 0; $l < $kb->col; ++$l) {
                            $y = $j + $q - $l;
                            if ($y >= 0 && $y < $ma->col) {
                                $sgima += $ma->data[$x * $ma->col + $y] * $kb->data[$k * $kb->col + $l];
                            }
                        }
                    }
                    $rc->data[$i * $ma->col + $j] = $sgima;
                }
            }
            return $rc;
        }
        else {
            throw new \Exception('Err::given matrixes has different data type!');
        }
    }
}