<?php

namespace numphp\benchmarks\matrix\structural;

use numphp\matrix;

/**
 * @Groups({"Structural"})
 * @BeforeMethods({"setUp"})
 */
class matrixJoinBelowBench
{
    /**
     * @var \numphp\matrix
     */
    protected $a;

    /**
     * @var \numphp\matrix
     */
    protected $b;

    public function setUp() : void
    {
        $this->a = matrix::uniform(500, 500);

        $this->b = matrix::uniform(500, 500);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function joinBelow() : void
    {
        $this->a->joinBelow($this->b);
    }
}
