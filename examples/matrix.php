<?php

require __DIR__ . '/../vendor/autoload.php';

use numphp\matrix;

matrix::time();
matrix::getMemory();
$a = matrix::ar([[23,56,45],[89,98,55],[56,75,65]]);
$b = \numphp\vector::ar([56,75,65]);

echo $a->clip(0,1);
echo PHP_EOL;
matrix::getMemory();           // get memory use
matrix::time();               // get time
/**
 * Memory-Consumed 7.7mb
 * Time-Consumed:- 0.19370794296265
 */