<?php

namespace numphp\benchmarks\matrix\reduction;

use numphp\matrix;

/**
 * @Groups({"Reductions"})
 * @BeforeMethods({"setUp"})
 */
class rrefBench
{
    /**
     * @var \numphp\matrix
     */
    protected $a;

    public function setUp() : void
    {
        $this->a = matrix::uniform(500, 500);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function rref() : void
    {
        $this->a->rref();
    }
}
