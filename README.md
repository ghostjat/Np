[![Scrutinizer Code Quality](https://scrutinizer-ci.com/g/ghostjat/Np/badges/quality-score.png?b=main)](https://scrutinizer-ci.com/g/ghostjat/Np/?branch=main)
![Packagist PHP Version Support](https://img.shields.io/packagist/php-v/ghostjat/Np)
[![Build Status](https://scrutinizer-ci.com/g/ghostjat/Np/badges/build.png?b=main)](https://scrutinizer-ci.com/g/ghostjat/Np/build-status/main)
[![Code Intelligence Status](https://scrutinizer-ci.com/g/ghostjat/Np/badges/code-intelligence.svg?b=main)](https://scrutinizer-ci.com/code-intelligence)
![GitHub contributors](https://img.shields.io/github/contributors/ghostjat/Np)
![GitHub commit activity](https://img.shields.io/github/commit-activity/m/ghostjat/Np)
![GitHub last commit](https://img.shields.io/github/last-commit/ghostjat/Np)
![Packagist Version](https://img.shields.io/packagist/v/ghostjat/Np)
![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/ghostjat/Np)
![Lines of code](https://img.shields.io/tokei/lines/github/ghostjat/Np)
![GitHub top language](https://img.shields.io/github/languages/top/ghostjat/Np)
![Np](https://github.com/ghostjat/numphp/blob/main/np.png)

## Description
   -----------
Lite, Fast &amp; Memory Efficient *Mathematical PHP library for scientific computing*

Np(numphp) is a library that provides objects for computing large sets of numbers in [PHP](https://php.net).

## Installation
Install [Np](https://packagist.org/packages/ghostjat/np) into your project with [Composer](https://getcomposer.org/):

```sh
$ composer require ghostjat/np
```
##Sample Code
```php
require __DIR__ . '/../vendor/autoload.php';
use Np\matrix;

$ta = matrix::randn(1000, 1000);    
$tb = matrix::randn(1000, 1000); // to generate random 2d matrix
$ta->dot($tb);                  // do a dot operation on given matrix
$ta->getMemory();              // get memory use
$ta->time();                  // get time
/**
 * 7.7mb
 * Time-Consumed:- 0.18390893936157
 */
```
*Synopsis*
--------
WARNING:  
This module is in its early stages and should be considered a Work in Progress.The interface is not final and may change in the future. 

*Requirements*
------------
- [PHP](https://php.net) 8+ 64bit with ffi & #libblas, #liblapacke   

Make sure you have all the necessary tools installed such as FFI, libblas, liblapacke.

*Performance*
-----------

System Conf:- Intel(R) Core(TM) i3-2370M CPU @ 2.40GHz 64bit 
Memory:- 8GB
php:- 8.0.5 64bit

*Current Benchmarks of this library*
-----------------------------------
![Benckmark](https://github.com/ghostjat/numphp/blob/main/npbm.png)

Data Size :- [500x500] Revolutions:- 5 Iterations:- 5

| subject  | mem_peak | best   | mode   | mean   | worst  | stdev |  
|----------|----------|--------|--------|--------|--------|-------|
| sum      | 3.606mb  | 0.014s | 0.014s | 0.015s | 0.015s | 0.000s| 
| multiply | 8.589mb  | 0.070s | 0.071s | 0.071s | 0.071s | 0.000s|
| lu       | 4.648mb  | 0.064s | 0.065s | 0.065s | 0.068s | 0.001s|
| eign     | 2.801mb  | 0.085s | 0.086s | 0.086s | 0.088s | 0.001s|
| cholesky | 1.621mb  | 0.001s | 0.001s | 0.001s | 0.001s | 0.000s|
| svd      | 3.706mb  | 0.126s | 0.126s | 0.127s | 0.133s | 0.002s|
| normL2   | 1.621mb  | 0.003s | 0.003s | 0.003s | 0.003s | 0.000s|
| Pinverse | 4.903mb  | 0.156s | 0.156s | 0.158s | 0.163s | 0.003s|
| inverse  | 1.819mb  | 0.016s | 0.016s | 0.016s | 0.017s | 0.000s|
| normL1   | 1.621mb  | 0.001s | 0.001s | 0.001s | 0.001s | 0.000s|
| dotMatrix| 3.769mb  | 0.006s | 0.006s | 0.006s | 0.006s | 0.000s|
| det      | 4.662mb  | 0.066s | 0.066s | 0.067s | 0.067s | 0.000s|
| rref     | 1.529mb  | 9.227s | 9.271s | 9.309s | 9.427s | 0.072s|
| ref      | 1.818mb  | 0.007s | 0.008s | 0.008s | 0.008s | 0.000s|
| clip     | 8.516mb  | 0.073s | 0.076s | 0.075s | 0.077s | 0.002s|
| clipUpper| 8.516mb  | 0.055s | 0.056s | 0.057s | 0.059s | 0.002s|
| clipLower| 8.516mb  | 0.055s | 0.058s | 0.057s | 0.059s | 0.002s|
| joinBelow| 4.517mb  | 0.027s | 0.027s | 0.027s | 0.028s | 0.000s|
| transpose| 8.504mb  | 0.057s | 0.057s | 0.058s | 0.059s | 0.001s|
| joinLeft | 4.511mb  | 0.025s | 0.025s | 0.026s | 0.027s | 0.001s|
| poisson  | 1.590mb  | 0.029s | 0.029s | 0.029s | 0.030s | 0.000s|
| gaussian | 20.203mb | 0.056s | 0.056s | 0.056s | 0.056s | 0.000s|
| randn    | 1.528mb  | 0.017s | 0.017s | 0.017s | 0.017s | 0.000s|
| uniform  | 1.528mb  | 0.021s | 0.021s | 0.021s | 0.022s | 0.000s|
| multiply | 4.507mb  | 0.042s | 0.042s | 0.043s | 0.045s | 0.001s|

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