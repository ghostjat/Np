<?php

namespace Np\benchmarks\matrix\lineAlgb;

use Np\matrix;

/**
 * @Groups({"LinearAlgebra"})
 * @BeforeMethods({"setUp"})
 */
class dotMatrixBench {

    /**
     * @var \Np\matrix
     */
    protected $a;

    /**
     * @var \Np\matrix
     */
    protected $b;

    public function setUp(): void {
        $this->a = matrix::uniform(500, 500);

        $this->b = matrix::uniform(500, 500);
    }

    /**
     * @Subject
     * @Iterations(5)
     * @Revs(5)
     * @OutputTimeUnit("seconds", precision=3)
     */
    public function dot(): void {
        $this->a->dot($this->b);
    }

}
