[![Scrutinizer Code Quality](https://scrutinizer-ci.com/g/ghostjat/numphp/badges/quality-score.png?b=main)](https://scrutinizer-ci.com/g/ghostjat/numphp/?branch=main)
![Packagist PHP Version Support](https://img.shields.io/packagist/php-v/ghostjat/numphp)
[![Build Status](https://scrutinizer-ci.com/g/ghostjat/numphp/badges/build.png?b=main)](https://scrutinizer-ci.com/g/ghostjat/numphp/build-status/main)
[![Code Intelligence Status](https://scrutinizer-ci.com/g/ghostjat/numphp/badges/code-intelligence.svg?b=main)](https://scrutinizer-ci.com/code-intelligence)
![GitHub contributors](https://img.shields.io/github/contributors/ghostjat/numphp)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/ghostjat/numphp)
![GitHub last commit](https://img.shields.io/github/last-commit/ghostjat/numphp)
![Packagist Version](https://img.shields.io/packagist/v/ghostjat/numphp)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/ghostjat/numphp)
![Lines of code](https://img.shields.io/tokei/lines/github/ghostjat/numphp)
![GitHub top language](https://img.shields.io/github/languages/top/ghostjat/numphp)
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
- [PHP](https://php.net) 7.4+ with ffi & #libblas, #liblapacke   

Make sure you have all the necessary tools installed such as FFI, libblas, liblapacke.

Performance
-----------
System-Conf:- 
CPU:- Intel(R) Core(TM) i3-2370M CPU @ 2.40GHz 64bit
MEM:- 8GB
PhpBench @git_tag@. Running benchmarks.
Using configuration file: /home/ghost/projects/git/numphp/phpbench.json

| benchmark                 | subject   | set | revs | its | mem_peak | mode    | rstdev   |
|---------------------------|-----------|-----|------|-----|----------|---------|----------|
| eignBench                 | eign      | 0   | 1    | 5   | 2.699mb  | 0.309s  | ±4.51%   |
| svdBench                  | svd       | 0   | 1    | 5   | 3.604mb  | 0.148s  | ±3.60%   |
| poissonMatrixBench        | poisson   | 0   | 1    | 5   | 11.738mb | 0.105s  | ±7.07%   |
| gaussianMatrixBench       | gaussian  | 0   | 1    | 5   | 11.738mb | 0.112s  | ±17.12%  |
| randMatrixBench           | randn     | 0   | 1    | 5   | 1.429mb  | 0.048s  | ±2.37%   |
| uniformMatrixBench        | uniform   | 0   | 1    | 5   | 1.429mb  | 0.063s  | ±8.16%   |
| matrixTransposeBench      | transpose | 0   | 1    | 5   | 8.431mb  | 0.120s  | ±1.32%   |
| rrefBench                 | rref      | 0   | 1    | 5   | 1.501mb  | 28.513s | ±1.90%   |
| refBench                  | ref       | 0   | 1    | 5   | 1.731mb  | 0.023s  | ±7.24%   |
| sumMatrixBench            | sum       | 0   | 1    | 5   | 2.434mb  | 0.051s  | ±3.59%   |
| matrixPseudoInverseBench  | inverse   | 0   | 1    | 5   | 4.775mb  | 0.222s  | ±13.76%  |
| matrixInverseBench        | inverse   | 0   | 1    | 5   | 1.731mb  | 0.032s  | ±127.50% |
| dotMatrixBench            | dotMatrix | 0   | 1    | 5   | 3.656mb  | 0.013s  | ±27.94%  |


License
-------
The code is licensed [MIT](LICENSE) and the documentation is licensed [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).

Author
------
Shubham Chaudhary <ghost.jat@gmail.com>