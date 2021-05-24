<?php

namespace numphp\benchmarks\matrix\special;

use numphp\matrix;

/**
 * @Groups({"Special"})
 * @BeforeMethods({"setUp"})
 */
class matrixClippingBench
{
    /**
     * @var \numphp\matrix
     */
    protected $a;

    /**
     * @var \numphp\matrix
     */
    protected $kernel;

    public function setUp() : void
    {
        $this->a = matrix::uniform(1000, 1000);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function clip() : void
    {
        $this->a->clip(0, 1);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function clipUpper() : void
    {
        $this->a->clipUpper(0);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function clipLower() : void
    {
        $this->a->clipLower(0);
    }
}
