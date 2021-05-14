<?php

namespace numphp\benchmarks\matrix\decomposition;

use numphp\matrix;

/**
 * @Groups({"Decompositions"})
 * @BeforeMethods({"setUp"})
 */
class choleskyBench
{
    /**
     * @var \numphp\matrix
     */
    protected $a;

    public function setUp() : void
    {
        $this->a = matrix::randn(500, 500);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function cholesky() : void
    {
        $this->a->cholesky();
    }
}
