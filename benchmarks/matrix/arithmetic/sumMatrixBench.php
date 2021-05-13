<?php

namespace numphp\benchmarks\matrix\arithmetic;

use numphp\matrix;

/**
 * @Groups({"Arithmetic"})
 * @BeforeMethods({"setUp"})
 */
class sumMatrixBench
{
    /**
     * @var \numphp\matrix
     */
    protected $a,$b;

    public function setUp() : void {
        $this->a = matrix::uniform(500, 500);
        $this->b = matrix::uniform(500, 500);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function sum() : void
    {
        $this->a->sum($this->b);
    }
}
