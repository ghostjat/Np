<?php
declare(strict_types=1);

namespace numphp;

class convolve {
    
    /**
     * 1D convolution between a vector va and kernel kb, with a given stride.
     * @param \numphp\vector $va
     * @param \numphp\vector $kb
     * @param int $stride
     */
    public static function conv1D(\numphp\vector $va, \numphp\vector $kb, int $stride = 1): vector {
        if ($va->dtype == $kb->dtype) {
            $nc = $va->col + $kb->col - 1;
            $rc = vector::factory($nc / $stride);
            for ($i = 0; $i < $nc; $i += $stride) {
                $jmin = $i >= $kb->col - 1 ? $i - ($kb->col - 1) : 0;
                $jmax = $i <= $va->col ? $i : $va->col - 1;
                $sigma = 0.0;
                for ($j = $jmin; $j <= $jmax; ++$j) {
                    $sigma += $va->data[$j] * $kb->data[$i - $j];
                }
                $rc->data[$i] = $sigma;
            }
            return $rc;
        }
        else {
            throw new \Exception('Err::given vectors has diffrent data type!');
        }
    }
    
    /**
     * 2D convolution between a matrix ma and kernel kb, with a given stride.
     * @param \numphp\matrix $ma
     * @param \numphp\matrix $kb
     * @param int $stride
     * @return matrix
     */
    public static function conv2D(\numphp\matrix $ma, \numphp\matrix $kb, int $stride = 1): matrix{
        if ($ma->dtype == $kb->dtype) {
            $p = $kb->row / 2;
            $q = $kb->col / 2;
            $rc = matrix::factory($ma->row / $stride, $ma->col / $stride, $ma->dtype);
            #core\blas::sgemm($ma, $kb, $rc);

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