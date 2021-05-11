<?php

require __DIR__ . '/../vendor/autoload.php';

use numphp\matrix;

matrix::time();
matrix::getMemory();
$ta = matrix::ar([[21,12,32],[432,322,112,],[22,342,21]]); // to genrate random 2d tensor
$c = $ta->multiply(\numphp\vector::ar([2,2,2]));          // do a dot opreation on given tensor
echo $c;
matrix::getMemory();           // get memory use
matrix::time();               // get time
/**
 * Memory-Consumed 7.7mb
 * Time-Consumed:- 0.19370794296265
 */