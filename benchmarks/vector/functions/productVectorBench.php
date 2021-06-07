<?php

namespace Np\benchmarks\vector\functions;

use Np\vector;

/**
 * @Groups({"Functions"})
 * @BeforeMethods({"setUp"})
 */
class productVectorBench
{
    /**
     * @var \Np\vector
     */
    protected $a;

    public function setUp() : void {
        $this->a = vector::uniform(100000);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("milliseconds", precision=3)
     */
    public function product() : void {
        $this->a->product();
    }
}
