<?php

namespace numphp\benchmarks\matrix\random;

use numphp\matrix;

/**
 * @Groups({"Random"})
 */
class uniformMatrixBench
{
    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function uniform() : void
    {
        matrix::uniform(500, 500);
    }
}
