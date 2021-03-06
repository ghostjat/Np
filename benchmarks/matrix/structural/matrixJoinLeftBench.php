<?php

namespace Np\benchmarks\matrix\structural;

use Np\matrix;

/**
 * @Groups({"Structural"})
 * @BeforeMethods({"setUp"})
 */
class matrixJoinLeftBench {

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
    public function joinLeft(): void {
        $this->a->joinLeft($this->b);
    }

}
