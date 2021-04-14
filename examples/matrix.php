<?php

require __DIR__ . '/../vendor/autoload.php';

use numphp\matrix;

matrix::time();
matrix::getMemory();
$ta = matrix::randn(1000,1000); // to genrate random 2d tensor
$tb = matrix::randn(1000,1000);
$ta->dotMatrix($tb);            // do a dot opreation on given tensor

matrix::getMemory();           // get memory use
matrix::time();               // get time
/**
 * Memory-Consumed 7.7mb
 * Time-Consumed:- 0.19370794296265
 */