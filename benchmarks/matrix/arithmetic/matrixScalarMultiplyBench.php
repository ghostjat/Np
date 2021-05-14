<?php

namespace numphp\benchmarks\matrix\arithmetic;

use numphp\matrix;

/**
 * @Groups({"Arithmetic"})
 * @BeforeMethods({"setUp"})
 */
class matrixScalarMultiplyBench
{
    /**
     * @var \numphp\matrix
     */
    protected $a;

    /**
     * @var float
     */
    protected $b = 3.14;

    public function setUp() : void
    {
        $this->a = matrix::uniform(1000, 1000);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function multiply() : void
    {
        $this->a->multiply($this->b);
    }
}
