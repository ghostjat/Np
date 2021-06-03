<?php

namespace Np\benchmarks\matrix\random;

use Np\matrix;

/**
 * @Groups({"Random"})
 */
class gaussianMatrixBench {

    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function gaussian(): void {
        matrix::gaussian(500, 500);
    }

}
