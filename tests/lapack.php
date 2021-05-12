<?php
require __DIR__ . '/../vendor/autoload.php';

use numphp\matrix;
matrix::time();
matrix::getMemory();
$mat = matrix::ar([[21,12,32],[432,322,112,],[22,342,21]]);
$i = $mat->inverse();
echo $i;
matrix::time();
matrix::getMemory();


