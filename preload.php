<?php

FFI::load(__DIR__ . '/src/core/blas.h');
FFI::load(__DIR__ . '/src/core/lapack.h');
if (!opcache_is_script_cached(__DIR__ . '/src/core/nd.php')) {
    opcache_compile_file(__DIR__ . '/src/core/nd.php');
}
if (!opcache_is_script_cached(__DIR__ . '/src/core/blas.php')) {
    opcache_compile_file(__DIR__ . '/src/core/blas.php');
}
if (!opcache_is_script_cached(__DIR__ . '/src/core/lapack.php')) {
    opcache_compile_file(__DIR__ . '/src/core/lapack.php');
}
if (!opcache_is_script_cached(__DIR__ . '/src/ops.php')) {
    opcache_compile_file(__DIR__ . '/src/ops.php');
}
if (!opcache_is_script_cached(__DIR__ . '/src/linAlg.php')) {
    opcache_compile_file(__DIR__ . '/src/linAlg.php');
}
if (!opcache_is_script_cached(__DIR__ . '/src/matrix.php')) {
    opcache_compile_file(__DIR__ . '/src/matrix.php');
}
if (!opcache_is_script_cached(__DIR__ . '/src/vector.php')) {
    opcache_compile_file(__DIR__ . '/src/vector.php');
}
if (!opcache_is_script_cached(__DIR__ . '/src/reductions/rref.php')) {
    opcache_compile_file(__DIR__ . '/src/reductions/rref.php');
}
if (!opcache_is_script_cached(__DIR__ . '/src/convolve.php')) {
    opcache_compile_file(__DIR__ . '/src/convolve.php');
}
