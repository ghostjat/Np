<?php

namespace numphp\benchmarks\matrix\trigonometric;

use numphp\matrix;

/**
 * @Groups({"Trigonometric"})
 * @BeforeMethods({"setUp"})
 */
class matrixTanBench
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
    public function tangent() : void
    {
        $this->a->tan();
    }
}
