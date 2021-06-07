<?php

namespace Np\benchmarks\matrix\trigonometric;

use Np\matrix;

/**
 * @Groups({"Trigonometric"})
 * @BeforeMethods({"setUp"})
 */
class cosMatrixBench
{
    /**
     * @var \Np\matrix
     */
    protected $a;

    public function setUp() : void {
        $this->a = matrix::uniform(500, 500);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function cos() : void {
        $this->a->cos();
    }
}
