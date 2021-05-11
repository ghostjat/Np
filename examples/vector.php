<?php

require __DIR__ . '/../vendor/autoload.php';

use numphp\{vector,matrix};

matrix::time();
matrix::getMemory();
$ta = vector::randn(3);        // to genrate random vector
$tb = matrix::randn(3,3);
$c = $ta->mulVectorMatrix($tb);            // do a dot opreation on given tensor
echo $c;
unset($c);
$c = $ta->multiplyMatrix($tb);
echo $c;
matrix::getMemory();           // get memory use
matrix::time();               // get time
/**
 * Memory-Consumed 7.7mb
 * Time-Consumed:- 0.19370794296265
 */