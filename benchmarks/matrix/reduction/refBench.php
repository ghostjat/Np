<?php

namespace Np\benchmarks\matrix\reduction;

use Np\matrix;

/**
 * @Groups({"Reductions"})
 * @BeforeMethods({"setUp"})
 */
class refBench {

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
    public function ref(): void {
        $this->a->ref();
    }

}
