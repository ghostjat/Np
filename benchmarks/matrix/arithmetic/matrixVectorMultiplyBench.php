<?php

namespace Np\benchmarks\matrix\arithmetic;

use Np\{matrix,vector};

/**
 * @Groups({"Arithmetic"})
 * @BeforeMethods({"setUp"})
 */
class matrixVectorMultiplyBench
{
    /**
     * @var \Np\matrix
     */
    protected $a;

    /**
     * @var \Np\vector
     */
    protected $b;

    public function setUp() : void
    {
        $this->a = matrix::uniform(1000, 1000);

        $this->b = vector::uniform(1000);
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
