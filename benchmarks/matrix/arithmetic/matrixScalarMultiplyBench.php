<?php

namespace Np\benchmarks\matrix\arithmetic;

use Np\matrix;

/**
 * @Groups({"Arithmetic"})
 * @BeforeMethods({"setUp"})
 */
class matrixScalarMultiplyBench
{
    /**
     * @var \Np\matrix
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
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function multiply() : void
    {
        $this->a->multiply($this->b);
    }
}
