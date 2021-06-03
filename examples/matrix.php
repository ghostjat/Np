<?php

require __DIR__ . '/../vendor/autoload.php';

use Np\matrix;

matrix::time();
matrix::getMemory();
$a = matrix::randn(1000, 1000);
$b = matrix::randn(1000, 1000);
$a->dot($b);
matrix::getMemory();           // get memory use
matrix::time();               // get time
/**
 * Memory-Consumed 7.7mb
 * Time-Consumed:- 0.18390893936157
 */