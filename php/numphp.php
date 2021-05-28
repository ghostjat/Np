<?php 
$unit = ['b', 'kb', 'mb', 'gb', 'tb', 'pb'];
$time = microtime(1);
$mem = memory_get_usage(); 
$ffi = FFI::cdef('
    typedef struct _matrix {
        int rows;
        int cols;
        double *data;
    } matrix;
    typedef struct _lup lup;
    matrix *mat(int rows, int cols);
    matrix *matIdentity(int rows, int cols);
    matrix *matRandom(int rows, int cols);
    matrix *matUniform(int rows, int cols);
    matrix *matPoisson(int rows, int cols, double lambda);
    matrix *matMinimum(matrix *m1, matrix *m2);
    matrix *matMultiply(matrix* m1, matrix* m2);
    matrix *matDgemm(matrix* m1, matrix* m2);
    matrix *matJoinLeft(matrix *m1, matrix *m2);
    matrix *matJoinRight(matrix *m1, matrix *m2);
    matrix *matJoinAbove(matrix *m1, matrix *m2);
    matrix *matJoinBelow(matrix *m1, matrix *m2);
    matrix *matRef(matrix *m);
    matrix *matInverse(matrix *m);
    double matTrace(matrix *m);
    void dignoalInterChange(matrix *m);
    int isSquare(matrix* m);
    int checkDimensions(matrix *m1, matrix *m2);
    int checkShape(matrix *m1, matrix *m2);
    void matPrint(matrix* mat);
        ','/home/ghost/projects/ffi/c/numphp/libnumphp.so');


$row = $col = 3;
$np1 = $ffi->matRandom($row,$col);
#$ffi->matPrint($np1);
#echo PHP_EOL;
$a =$ffi->matInverse($np1);
$ffi->matPrint($a);
#echo PHP_EOL;
echo ' time:-' . (microtime(1) - $time) .PHP_EOL;
$memory = memory_get_usage() - $mem;
echo round($memory / pow(1024, ($i = floor(log($memory, 1024)))), 2) . $unit[$i] . PHP_EOL;
