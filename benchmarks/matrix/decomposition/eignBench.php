<?php

namespace Np\benchmarks\matrix\decomposition;

use Np\matrix;

/**
 * @Groups({"Decompositions"})
 */
class eignBench
{
    /**
     * @var \Np\matrix
     */
    protected $a;

    public function setUp() : void
    {
        $this->a = matrix::uniform(500, 500);

        $this->a = $this->a->dot($this->a);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @revs(5)
     * @BeforeMethods({"setUp"})
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function eign() : void
    {
        $this->a->eign();
    }
}
