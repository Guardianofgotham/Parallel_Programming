#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include "hclib_cpp.h"
#include <inttypes.h>

using namespace hclib;

 struct t  {
    int a ;
    int b ;
}__attribute__( ( aligned ( 64 ) ) ) typedef t;


int main(){
    char const *deps[] = { "system" };
    int count = 0;
    t d;
    printf("%lu\n",sizeof(d ));
    //printf("size of int is: %lu bytes\n",sizeof(int));
    return 0;
}