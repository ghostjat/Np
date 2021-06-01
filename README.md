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
Lite, Fast &amp; Memory Efficient php-FFI library for scientific computing

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
$ta->dot($tb);                  // do a dot operation on given matrix
matrix::getMemory();           // get memory use
matrix::time();               // get time

/**
 * 7.7mb
 * Time-Consumed:- 0.18390893936157
 */
```
Synopsis
--------
WARNING:  
This module is in its early stages and should be considered a Work in Progress.
The interface is not final and may change in the future. 

Requirements
------------
- [PHP](https://php.net) 8+ with ffi & #libblas, #liblapacke   

Make sure you have all the necessary tools installed such as FFI, libblas, liblapacke.

Performance
-----------
System-Conf:- 
CPU:- Intel(R) Core(TM) i3-2370M CPU @ 2.40GHz 64bit
MEM:- 8GB
PhpBench @git_tag@. Running benchmarks.
Using configuration file: /home/ghost/projects/git/numphp/phpbench.json


<p align="center">
  <img src="https://github.com/ghostjat/numphp/blob/main/Numphp Benchmarks.png">
</p>

Current Benchmark

| benchmark                 | subject  | set | revs | its | mem_peak | best   | mode   | mean   | worst  | stdev  | rstdev |
|---------------------------|----------|-----|------|-----|----------|--------|--------|--------|--------|--------|--------|
|sumMatrixBench             | sum      | 0   | 5    | 5   | 3.606mb  | 0.014s | 0.014s | 0.015s | 0.015s | 0.000s | ±1.57% |
|matrixVectorMultiplyBench  | multiply | 0   | 5    | 5   | 8.589mb  | 0.070s | 0.071s | 0.071s | 0.071s | 0.000s | ±0.23% |
|luBench                    | lu       | 0   | 5    | 5   | 4.648mb  | 0.064s | 0.065s | 0.065s | 0.068s | 0.001s | ±1.73% |
|eignBench                  | eign     | 0   | 5    | 5   | 2.801mb  | 0.085s | 0.086s | 0.086s | 0.088s | 0.001s | ±1.20% |
|choleskyBench              | cholesky | 0   | 5    | 5   | 1.621mb  | 0.001s | 0.001s | 0.001s | 0.001s | 0.000s | ±0.91% |
|svdBench                   | svd      | 0   | 5    | 5   | 3.706mb  | 0.126s | 0.126s | 0.127s | 0.133s | 0.002s | ±1.72% |
|matrixL2NormBench          | normL2   | 0   | 5    | 5   | 1.621mb  | 0.003s | 0.003s | 0.003s | 0.003s | 0.000s | ±0.18% |
|matrixPseudoInverseBench   | inverse  | 0   | 5    | 5   | 4.903mb  | 0.156s | 0.156s | 0.158s | 0.163s | 0.003s | ±1.87% |
|matrixInverseBench         | inverse  | 0   | 5    | 5   | 1.819mb  | 0.016s | 0.016s | 0.016s | 0.017s | 0.000s | ±1.75% |
|matrixL1NormBench          | normL1   | 0   | 5    | 5   | 1.621mb  | 0.001s | 0.001s | 0.001s | 0.001s | 0.000s | ±0.35% |
|dotMatrixBench             | dotMatrix| 0   | 5    | 5   | 3.769mb  | 0.006s | 0.006s | 0.006s | 0.006s | 0.000s | ±1.40% |
|matrixDeterminantBench     | det      | 0   | 5    | 5   | 4.662mb  | 0.066s | 0.066s | 0.067s | 0.067s | 0.000s | ±0.56% |
|rrefBench                  | rref     | 0   | 5    | 5   | 1.529mb  | 9.227s | 9.271s | 9.309s | 9.427s | 0.072s | ±0.77% |
|refBench                   | ref      | 0   | 5    | 5   | 1.818mb  | 0.007s | 0.008s | 0.008s | 0.008s | 0.000s | ±1.79% |
|matrixClippingBench        | clip     | 0   | 5    | 5   | 8.516mb  | 0.073s | 0.076s | 0.075s | 0.077s | 0.002s | ±2.42% |
|matrixClippingBench        | clipUpper| 0   | 5    | 5   | 8.516mb  | 0.055s | 0.056s | 0.057s | 0.059s | 0.002s | ±3.10% |
|matrixClippingBench        | clipLower| 0   | 5    | 5   | 8.516mb  | 0.055s | 0.058s | 0.057s | 0.059s | 0.002s | ±2.82% |
|matrixJoinBelowBench       | joinBelow| 0   | 5    | 5   | 4.517mb  | 0.027s | 0.027s | 0.027s | 0.028s | 0.000s | ±0.99% |
|matrixTransposeBench       | transpose| 0   | 5    | 5   | 8.504mb  | 0.057s | 0.057s | 0.058s | 0.059s | 0.001s | ±0.92% |
|matrixJoinLeftBench        | joinLeft | 0   | 5    | 5   | 4.511mb  | 0.025s | 0.025s | 0.026s | 0.027s | 0.001s | ±2.18% |
|poissonMatrixBench         | poisson  | 0   | 5    | 5   | 1.590mb  | 0.029s | 0.029s | 0.029s | 0.030s | 0.000s | ±1.63% |
|gaussianMatrixBench        | gaussian | 0   | 5    | 5   | 20.203mb | 0.056s | 0.056s | 0.056s | 0.056s | 0.000s | ±0.65% |
|randMatrixBench            | randn    | 0   | 5    | 5   | 1.528mb  | 0.017s | 0.017s | 0.017s | 0.017s | 0.000s | ±0.48% |
|uniformMatrixBench         | uniform  | 0   | 5    | 5   | 1.528mb  | 0.021s | 0.021s | 0.021s | 0.022s | 0.000s | ±1.16% |
|matrixScalarMultiplyBench  | multiply | 0   | 5    | 5   | 4.507mb  | 0.042s | 0.042s | 0.043s | 0.045s | 0.001s | ±2.01% |

Previous BenchMark

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
| matrixL1NormBench         | normL1    | 0   | 1    | 10  | 1.525mb  | 0.001s  | ±0.80%   |
| matrixL2NormBench         | normL2    | 0   | 1    | 10  | 1.525mb  | 0.003s  | ±1.63%   |

License
-------
The code is licensed [MIT](LICENSE) and the documentation is licensed [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).

Author
------
Shubham Chaudhary <ghost.jat@gmail.com>