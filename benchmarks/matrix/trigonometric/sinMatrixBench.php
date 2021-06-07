<?php

namespace Np\benchmarks\matrix\trigonometric;

use Np\matrix;

/**
 * @Groups({"Trigonometric"})
 * @BeforeMethods({"setUp"})
 */
class sinMatrixBench
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
    public function sine() : void
    {
        $this->a->sin();
    }
}
