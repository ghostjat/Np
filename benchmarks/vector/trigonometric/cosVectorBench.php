<?php

namespace Np\benchmarks\vector\trigonometric;

use Np\vector;

/**
 * @Groups({"Trigonometric"})
 * @BeforeMethods({"setUp"})
 */
class cosVectorBench {

    /**
     * @var \Np\vector
     */
    protected $a;

    public function setUp(): void {
        $this->a = vector::uniform(100000);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @OutputTimeUnit("milliseconds", precision=3)
     */
    public function cos(): void {
        $this->a->cos();
    }

}
