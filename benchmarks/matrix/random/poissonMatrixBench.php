<?php

namespace numphp\benchmarks\matrix\random;

use numphp\matrix;

/**
 * @Groups({"Random"})
 */
class poissonMatrixBench
{
    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function poisson() : void
    {
        matrix::poisson(500, 500);
    }
}
