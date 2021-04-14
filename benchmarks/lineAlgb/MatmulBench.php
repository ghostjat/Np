<?php

namespace numphp\benchmarks\lineAlgb;

use numphp\matrix;

/**
 * @Groups({"LinearAlgebra"})
 * @BeforeMethods({"setUp"})
 */
class dotMatrixBench
{
    /**
     * @var \numphp\matrix
     */
    protected $a;

    /**
     * @var \numphp\matrix
     */
    protected $b;

    public function setUp() : void
    {
        $this->a = matrix::uniform(1000, 1000);

        $this->b = matrix::uniform(1000, 1000);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @Revs(1)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function dotMatrix() : void
    {
        $this->a->dotMatrix($this->b);
    }
}
