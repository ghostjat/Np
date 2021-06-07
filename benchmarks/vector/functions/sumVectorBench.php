<?php

namespace Np\benchmarks\vector\functions;

use Np\vector;

/**
 * @Groups({"Functions"})
 * @BeforeMethods({"setUp"})
 */
class sumVectorBench {
    /**
     * @var \Np\vector
     */
    protected $a;

    public function setUp() : void {
        $this->a = Vector::uniform(100000);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("milliseconds", precision=3)
     */
    public function sum() : void
    {
        $this->a->sum();
    }
}
