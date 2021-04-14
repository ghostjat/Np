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

## Requirements
- [PHP](https://php.net) 7.4+ with ffi & blas

Make sure you have all the necessary tools installed such as libblas,ffi.

## Performance

## License
The code is licensed [MIT](LICENSE) and the documentation is licensed [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).

Author
------
Shubham Chaudhary <ghost.jat@gmail.com>