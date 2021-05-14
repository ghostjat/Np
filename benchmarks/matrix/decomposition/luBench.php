<?php

namespace numphp\benchmarks\matrix\decomposition;

use numphp\matrix;

/**
 * @Groups({"Decompositions"})
 * @BeforeMethods({"setUp"})
 */
class luBench
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
     * @Revs(10)
     * @Iterations(10)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function lu() : void
    {
        $this->a->lu();
    }
}
