<?php

namespace numphp\benchmarks\matrix\decomposition;

use numphp\matrix;

/**
 * @Groups({"Decompositions"})
 */
class eignBench
{
    /**
     * @var \Tensor\Matrix
     */
    protected $a;

    public function setUp() : void
    {
        $this->a = matrix::uniform(500, 500);

        $this->a = $this->a->dotMatrix($this->a);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @BeforeMethods({"setUp"})
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function eign() : void
    {
        $this->a->eign();
    }
}
