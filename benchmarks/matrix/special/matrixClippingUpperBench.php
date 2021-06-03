<?php
namespace Np\benchmarks\matrix\special;

use Np\matrix;

/**
 * @Groups({"Special"})
 * @BeforeMethods({"setUp"})
 */
class matrixClippingUpperBench {
    
    /**
     * @var \Np\matrix
     */
    protected $a;

    /**
     * @var \Np\matrix
     */
    protected $kernel;

    public function setUp(): void {
        $this->a = matrix::uniform(1000, 1000);
    }
    
    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function clipUpper(): void {
        $this->a->clipUpper(0);
    }
}

