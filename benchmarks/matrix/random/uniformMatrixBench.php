<?php

namespace Np\benchmarks\matrix\random;

use Np\matrix;

/**
 * @Groups({"Random"})
 */
class uniformMatrixBench
{
    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function uniform() : void
    {
        matrix::uniform(500, 500);
    }
}
