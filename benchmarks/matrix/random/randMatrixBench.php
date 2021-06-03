<?php

namespace Np\benchmarks\matrix\random;

use Np\matrix;

/**
 * @Groups({"Random"})
 */
class randMatrixBench {

    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function randn(): void {
        matrix::randn(500, 500);
    }

}
