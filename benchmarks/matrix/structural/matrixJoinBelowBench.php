<?php

namespace Np\benchmarks\matrix\structural;

use Np\matrix;

/**
 * @Groups({"Structural"})
 * @BeforeMethods({"setUp"})
 */
class matrixJoinBelowBench {

    /**
     * @var \Np\matrix
     */
    protected $a;

    /**
     * @var \Np\matrix
     */
    protected $b;

    public function setUp(): void {
        $this->a = matrix::uniform(500, 500);

        $this->b = matrix::uniform(500, 500);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function joinBelow(): void {
        $this->a->joinBelow($this->b);
    }

}
