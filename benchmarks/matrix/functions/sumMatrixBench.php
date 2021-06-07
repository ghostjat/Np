<?php

namespace Np\benchmarks\matrix\functions;

use Np\matrix;

/**
 * @Groups({"Functions"})
 * @BeforeMethods({"setUp"})
 */
class sumMatrixBench {

    /**
     * @var \Np\matrix
     */
    protected $a;

    public function setUp(): void {
        $this->a = matrix::uniform(500, 500);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("milliseconds", precision=3)
     */
    public function sum(): void {
        $this->a->sumRows();
    }

}
