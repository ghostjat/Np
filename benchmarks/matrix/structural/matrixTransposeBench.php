<?php

namespace Np\benchmarks\matrix\structural;

use Np\matrix;

/**
 * @Groups({"Structural"})
 * @BeforeMethods({"setUp"})
 */
class matrixTransposeBench
{
    /**
     * @var \Np\matrix
     */
    protected $a;

    public function setUp() : void
    {
        $this->a = matrix::uniform(1000, 1000);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function transpose() : void
    {
        $this->a->transpose();
    }
}
