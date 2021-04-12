<?php

require __DIR__ . '/../vendor/autoload.php';

use numphp\tensor;

tensor::time();
tensor::getMemory();
$ta = tensor::randn(1000, 1000); // to genrate random 2d matrix
$tb = tensor::randn(1000, 1000);
$ta->dotMatrix($tb);            // do a dot opreation on given matrix
tensor::getMemory();           // get memory use
tensor::time();               // get time

/**
 * Memory-Consumed 15.3mb
 * Time-Consumed:- 0.26119709014893
 */