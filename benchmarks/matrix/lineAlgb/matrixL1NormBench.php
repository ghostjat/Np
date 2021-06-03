<?php

namespace Np\benchmarks\matrix\lineAlgb;

use Np\matrix;

/**
 * @Groups({"LinearAlgebra"})
 * @BeforeMethods({"setUp"})
 */
class matrixL1NormBench
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
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function normL1() : void
    {
        $this->a->normINF();
    }
}
