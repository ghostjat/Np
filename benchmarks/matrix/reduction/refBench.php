<?php

namespace numphp\benchmarks\matrix\reduction;

use numphp\matrix;

/**
 * @Groups({"Reductions"})
 * @BeforeMethods({"setUp"})
 */
class refBench
{
    /**
     * @var \Tensor\Matrix
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
    public function ref() : void
    {
        $this->a->ref();
    }
}
