<?php

namespace Np\benchmarks\vector\signalProcessing;

use Np\vector;


/**
 * @Groups({"Signal Processing"})
 * @BeforeMethods({"setUp"})
 */
class vectorConvolveBench
{
    /**
     * @var \Np\vector
     */
    protected $a;

    /**
     * @var \Np\vector
     */
    protected $kernel;

    public function setUp() : void
    {
        $this->a = vector::uniform(250000);

        $this->kernel = vector::uniform(100);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function convolve() : void {
        $this->a->convolve($this->kernel);
    }
}
