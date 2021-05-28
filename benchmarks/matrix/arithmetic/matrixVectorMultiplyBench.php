<?php

namespace numphp\benchmarks\matrix\arithmetic;

use numphp\{matrix,vector};

/**
 * @Groups({"Arithmetic"})
 * @BeforeMethods({"setUp"})
 */
class matrixVectorMultiplyBench
{
    /**
     * @var \numphp\matrix
     */
    protected $a;

    /**
     * @var \numphp\vector
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
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function multiply() : void
    {
        $this->a->multiply($this->b);
        #$this->b->multiplyMatrix($this->a);
    }
}
