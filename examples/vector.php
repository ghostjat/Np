<?php

require __DIR__ . '/../vendor/autoload.php';

use numphp\{vector,matrix};

matrix::time();
matrix::getMemory();
$v = vector::ar([1,2,3,4]);        // to genrate random vector

echo $v->product() . PHP_EOL;
matrix::getMemory();           // get memory use
matrix::time();               // get time
/**
 * Memory-Consumed 7.7mb
 * Time-Consumed:- 0.19370794296265
 */