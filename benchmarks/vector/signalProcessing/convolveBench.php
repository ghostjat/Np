<?php

namespace Np\benchmarks\vector\signalProcessing;

use Np\matrix;

/**
 * @Groups({"Signal Processing"})
 * @BeforeMethods({"setUp"})
 */
class convolveBench
{
    /**
     * @var \Np\matrix
     */
    protected $a;

    /**
     * @var \Np\matrix
     */
    protected $b;

    public function setUp() : void {
        $this->a = matrix::uniform(500,500);

        $this->b = matrix::uniform(50,50);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function convolve() : void
    {
        $this->a->convolve($this->b);
    }
}
