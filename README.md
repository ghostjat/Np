[![Scrutinizer Code Quality](https://scrutinizer-ci.com/g/ghostjat/numphp/badges/quality-score.png?b=main)](https://scrutinizer-ci.com/g/ghostjat/numphp/?branch=main)
[![Minimum PHP version: 7.4.0](https://img.shields.io/badge/php-7.4.0%2B-blue.svg)](https://php.net)
[![Build Status](https://scrutinizer-ci.com/g/ghostjat/numphp/badges/build.png?b=main)](https://scrutinizer-ci.com/g/ghostjat/numphp/build-status/main)
[![Code Intelligence Status](https://scrutinizer-ci.com/g/ghostjat/numphp/badges/code-intelligence.svg?b=main)](https://scrutinizer-ci.com/code-intelligence)

<p align="center">
  <img src="https://github.com/ghostjat/numphp/blob/main/numphp.png">
</p>
## Description
   -----------
a lite &amp; memory efficient  php library for scientific computing

numphp is a library that provides objects for computing large sets of numbers in [PHP](https://php.net).

## Installation
Install Numphp into your project with [Composer](https://getcomposer.org/):
```sh
$ composer require ghostjat/numphp
```
##Sample Code
```php
require __DIR__ . '/../vendor/autoload.php';

use numphp\matrix;

matrix::time();
matrix::getMemory();
$ta = matrix::randn(1000, 1000);    
$tb = matrix::randn(1000, 1000); // to generate random 2d matrix
$ta->dotMatrix($tb);            // do a dot operation on given matrix
matrix::getMemory();           // get memory use
matrix::time();               // get time

/**
 * Memory-Consumed 7.7mb
 * Time-Consumed:- 0.19370794296265
 */
```
Synopsis
--------
WARNING:  
This module is in its early stages and should be considered a Work in Progress.
The interface is not final and may change in the future. 

Requirements
------------
- [PHP](https://php.net) 7.4+ with ffi & blas

Make sure you have all the necessary tools installed such as libblas,ffi.

Performance
-----------
PhpBench @git_tag@. Running benchmarks.
Using configuration file: /home/ghost/projects/git/numphp/phpbench.json

\numphp\benchmarks\lineAlgb\dotMatrixBench (#0 dotMatrix)

#0  0.046 0.046 0.046 0.049 0.046 (s) [μ Mo]/r: 0.046 0.046 μRSD/r: 2.67%
-----------------------------------------------------------------------
1 subjects, 5 iterations, 1 revs, 0 rejects, 0 failures, 0 warnings

(best [mean mode] worst) = 0.046 [0.046 0.046] 0.049 (s)
--------------------------------------------------------
⅀
T: 0.232s μSD/r 0.001s μRSD/r: 2.672%

suite: 13462ee50fa758e45d873b58a32184e5a5f9a052, date: 2021-04-14, stime: 12:00:47

| benchmark      | subject   | set | revs | its | mem_peak    | best   | mean   | mode   | worst  | stdev  | rstdev | diff  |
|----------------|-----------|-----|------|-----|-------------|--------|--------|--------|--------|--------|--------|-------|
| dotMatrixBench | dotMatrix | 0   | 1    | 5   | 79,603,120b | 0.046s | 0.046s | 0.046s | 0.049s | 0.001s | 2.67%  | 1.00x |


License
-------
The code is licensed [MIT](LICENSE) and the documentation is licensed [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).

Author
------
Shubham Chaudhary <ghost.jat@gmail.com>