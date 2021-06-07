<?php

namespace Np\benchmarks\matrix\signalProcessing;

use Np\matrix;

/**
 * @Groups({"Signal Processing"})
 * @BeforeMethods({"setUp"})
 */
class matrixConvolveBench {

    /**
     * @var \Np\matrix
     */
    protected $a;

    /**
     * @var \Np\Matrix
     */
    protected $kernel;

    public function setUp(): void {
        $this->a = matrix::uniform(500, 500);

        $this->kernel = matrix::uniform(10, 10);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function convolve(): void {
        $this->a->convolve($this->kernel);
    }

}
