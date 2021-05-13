<?php

namespace numphp\benchmarks\matrix\random;

use numphp\matrix;

/**
 * @Groups({"Random"})
 */
class gaussianMatrixBench
{
    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function gaussian() : void
    {
        matrix::gaussian(500, 500);
    }
}
