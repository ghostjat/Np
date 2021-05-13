<?php

namespace numphp\benchmarks\matrix\lineAlgb;

use numphp\matrix;

/**
 * @Groups({"LinearAlgebra"})
 * @BeforeMethods({"setUp"})
 */
class matrixPseudoInverseBench
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
    public function inverse() : void
    {
        $this->a->pseudoInverse();
    }
}
