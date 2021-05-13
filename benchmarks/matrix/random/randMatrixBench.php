<?php

namespace numphp\benchmarks\matrix\random;

use numphp\matrix;

/**
 * @Groups({"Random"})
 */
class randMatrixBench
{
    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function randn() : void
    {
        matrix::randn(500, 500);
    }
}
