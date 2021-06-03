<?php

require __DIR__ . '/../vendor/autoload.php';

use Np\vector;
vector::time();
vector::getMemory();
$v = vector::ar(range(random_int(1,2), random_int(99999,9999999))); 

echo $v->sum() . PHP_EOL;
vector::getMemory();
vector::time();
