<?php

namespace Np\benchmarks\matrix\decomposition;

use Np\matrix;

/**
 * @Groups({"Decompositions"})
 * @BeforeMethods({"setUp"})
 */
class luBench
{
    /**
     * @var \Np\matrix
     */
    protected $a;

    public function setUp() : void
    {
        $this->a = matrix::uniform(500, 500);
    }

    /**
     * @Subject
     * @Revs(5)
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function lu() : void
    {
        $this->a->lu();
    }
}
