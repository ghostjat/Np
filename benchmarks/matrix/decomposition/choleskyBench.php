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
     * @var \Tensor\Matrix
     */
    protected $a;

    public function setUp() : void
    {
        $this->a = matrix::rand(500, 500);
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
