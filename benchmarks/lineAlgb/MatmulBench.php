<?php

namespace numphp\benchmarks\lineAlgb;

use numphp\tensor;

/**
 * @Groups({"LinearAlgebra"})
 * @BeforeMethods({"setUp"})
 */
class MatmulBench
{
    /**
     * @var \numphp\tensor
     */
    protected $a;

    /**
     * @var \numphp\tensor
     */
    protected $b;

    public function setUp() : void
    {
        $this->a = tensor::uniform(1000, 1000);

        $this->b = tensor::uniform(1000, 1000);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @Revs(1)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function matmul() : void
    {
        $this->a->dotMatrix($this->b);
    }
}
