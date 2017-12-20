#include <stdio.h>
#include <stdlib.h>

int32_t deletefile ( char *dfile)
{
   FILE *fp;
   int ret;

   fp = fopen(dfile, "r");

   if( fp == NULL )
     { return -1; }

   fclose(fp);

   ret = remove(dfile);

   if(ret != 0) {
      return -1;
   } 
   return 0;
}
