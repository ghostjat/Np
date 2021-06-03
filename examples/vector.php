<?php

require __DIR__ . '/../vendor/autoload.php';

use Np\vector;
$unit = ['b', 'kb', 'mb', 'gb', 'tb', 'pb'];
$time = microtime(1);
$mem = memory_get_usage();
$v = vector::ar(range(-100000, 100000));        // to genrate random vector

echo PHP_EOL;
$memory = memory_get_usage() - $mem . PHP_EOL;
echo round($memory / pow(1024, ($i = floor(log($memory, 1024)))), 2) . $unit[$i] . PHP_EOL;

echo microtime(1) - $time . PHP_EOL;
