<?php

namespace numphp\benchmarks\vector\trigonometric;

use numphp\vector;

/**
 * @Groups({"Trigonometric"})
 * @BeforeMethods({"setUp"})
 */
class vectorCosBench
{
    /**
     * @var \numphp\vector
     */
    protected $a;

    public function setUp() : void
    {
        $this->a = vector::uniform(100000);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @OutputTimeUnit("milliseconds", precision=3)
     */
    public function cosine() : void
    {
        $this->a->cos();
    }
}
