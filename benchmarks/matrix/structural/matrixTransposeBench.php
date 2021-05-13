<?php

namespace numphp\benchmarks\matrix\structural;

use numphp\matrix;

/**
 * @Groups({"Structural"})
 * @BeforeMethods({"setUp"})
 */
class matrixTransposeBench
{
    /**
     * @var \Tensor\Matrix
     */
    protected $a;

    public function setUp() : void
    {
        $this->a = matrix::uniform(1000, 1000);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function transpose() : void
    {
        $this->a->transpose();
    }
}
