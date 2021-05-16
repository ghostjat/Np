<?php

namespace numphp\benchmarks\vector\signalProcessing;

use numphp\matrix;

/**
 * @Groups({"Signal Processing"})
 * @BeforeMethods({"setUp"})
 */
class convolveBench
{
    /**
     * @var \numphp\matrix
     */
    protected $a;

    /**
     * @var \numphp\matrix
     */
    protected $b;

    public function setUp() : void {
        $this->a = matrix::uniform(500,500);

        $this->b = matrix::uniform(50,50);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function convolve() : void
    {
        $this->a->convolve($this->b);
    }
}
