<?php

namespace Np\benchmarks\matrix\functions;

use Np\matrix;

/**
 * @Groups({"Functions"})
 * @BeforeMethods({"setUp"})
 */
class expMatrixBench {

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
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function exp(): void {
        $this->a->exp();
    }

}
