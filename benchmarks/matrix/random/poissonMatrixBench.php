<?php

namespace Np\benchmarks\matrix\random;

use Np\matrix;

/**
 * @Groups({"Random"})
 */
class poissonMatrixBench {

    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function poisson(): void {
        matrix::poisson(500, 500);
    }

}
