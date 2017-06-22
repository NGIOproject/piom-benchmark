#include <stdio.h>
#include <stdlib.h>
#include "arralloc.h"

#define NCHAR 10
#define CSIZE  8

#define M 8
#define N 6

#define PAD 2
#define SCALE 1000

#define CBLKSIZE (SCALE*CSIZE)

#define MBLKSIZE (SCALE*2*(CSIZE+PAD))
#define NBLKSIZE (SCALE*(CSIZE+2*PAD)) 

#define MSIZE (M*MBLKSIZE)
#define NSIZE (N*NBLKSIZE)

char pixmap[NCHAR][CSIZE][CSIZE] =
{
{
{'.','.','X','X','X','.','.','.'},
{'.','X','.','.','.','X','.','.'},
{'X','.','.','.','.','.','X','.'},
{'X','.','.','.','.','.','X','.'},
{'X','.','.','.','.','.','X','.'},
{'.','X','.','.','.','X','.','.'},
{'.','.','X','X','X','.','.','.'},
{'.','.','.','.','.','.','.','.'},
},
{
{'.','.','.','X','.','.','.','.'},
{'.','.','X','X','.','.','.','.'},
{'.','X','.','X','.','.','.','.'},
{'.','.','.','X','.','.','.','.'},
{'.','.','.','X','.','.','.','.'},
{'.','.','.','X','.','.','.','.'},
{'.','X','X','X','X','X','.','.'},
{'.','.','.','.','.','.','.','.'},
},
{
{'.','X','X','X','X','X','.','.'},
{'X','.','.','.','.','.','X','.'},
{'.','.','.','.','.','.','X','.'},
{'.','X','X','X','X','X','.','.'},
{'X','.','.','.','.','.','.','.'},
{'X','.','.','.','.','.','.','.'},
{'X','X','X','X','X','X','X','.'},
{'.','.','.','.','.','.','.','.'},
},
{
{'.','X','X','X','X','X','.','.'},
{'X','.','.','.','.','.','X','.'},
{'.','.','.','.','.','.','X','.'},
{'.','X','X','X','X','X','.','.'},
{'.','.','.','.','.','.','X','.'},
{'X','.','.','.','.','.','X','.'},
{'.','X','X','X','X','X','.','.'},
{'.','.','.','.','.','.','.','.'},
},
{
{'X','.','.','.','.','.','.','.'},
{'X','.','.','.','.','X','.','.'},
{'X','.','.','.','.','X','.','.'},
{'X','.','.','.','.','X','.','.'},
{'X','X','X','X','X','X','X','.'},
{'.','.','.','.','.','X','.','.'},
{'.','.','.','.','.','X','.','.'},
{'.','.','.','.','.','.','.','.'},
},
{
{'X','X','X','X','X','X','X','.'},
{'X','.','.','.','.','.','.','.'},
{'X','.','.','.','.','.','.','.'},
{'X','X','X','X','X','X','.','.'},
{'.','.','.','.','.','.','X','.'},
{'X','.','.','.','.','.','X','.'},
{'.','X','X','X','X','X','.','.'},
{'.','.','.','.','.','.','.','.'},
},
{
{'.','X','X','X','X','X','.','.'},
{'X','.','.','.','.','.','X','.'},
{'X','.','.','.','.','.','.','.'},
{'X','X','X','X','X','X','.','.'},
{'X','.','.','.','.','.','X','.'},
{'X','.','.','.','.','.','X','.'},
{'.','X','X','X','X','X','.','.'},
{'.','.','.','.','.','.','.','.'},
},
{
{'X','X','X','X','X','X','X','.'},
{'X','.','.','.','.','X','.','.'},
{'.','.','.','.','X','.','.','.'},
{'.','.','.','X','.','.','.','.'},
{'.','.','X','.','.','.','.','.'},
{'.','.','X','.','.','.','.','.'},
{'.','.','X','.','.','.','.','.'},
{'.','.','.','.','.','.','.','.'},
},
{
{'.','X','X','X','X','X','.','.'},
{'X','.','.','.','.','.','X','.'},
{'X','.','.','.','.','.','X','.'},
{'.','X','X','X','X','X','.','.'},
{'X','.','.','.','.','.','X','.'},
{'X','.','.','.','.','.','X','.'},
{'.','X','X','X','X','X','.','.'},
{'.','.','.','.','.','.','.','.'},
},
{
{'.','X','X','X','X','X','.','.'},
{'X','.','.','.','.','.','X','.'},
{'X','.','.','.','.','.','X','.'},
{'.','X','X','X','X','X','X','.'},
{'.','.','.','.','.','.','X','.'},
{'X','.','.','.','.','.','X','.'},
{'.','X','X','X','X','X','.','.'},
{'.','.','.','.','.','.','.','.'},
}
};


#define MAXCHAR 128

int main(void){
  int ichar, i, j, im, in, i0, j0, is, js, ip, jp;
  char c0, c1;

  int zero=0, one=1;

  char filename[MAXCHAR];

  float **array;

  FILE *fp;

  array = arralloc(sizeof(float), 2, MSIZE, NSIZE);

  printf("CSIZE = %d, M = %d, N = %d\n", CSIZE, M, N);
  printf("PAD = %d, SCALE = %d\n", PAD, SCALE);
  printf("CLKSIZE = %d\n", CBLKSIZE);
  printf("MBLKSIZE = %d, NBLKSIZE = %d\n", MBLKSIZE, NBLKSIZE);
  printf("MSIZE = %d, NSIZE = %d\n", MSIZE, NSIZE);


  for (ichar=0; ichar < NCHAR; ichar++)
  {
    for (i=0; i < CSIZE; i++)
    {
      for (j=0; j < CSIZE; j++)
      {
        if (pixmap[ichar][i][j] != 'X')
        {
          pixmap[ichar][i][j] = '.';
        }
      }
    }
  }

  for (i=0; i < MSIZE; i++)
  {
    for (j=0; j < NSIZE; j++)
    {
      array[i][j] = 0.0;
    }
  }

  for (im = 0; im < M; im++)
  {
    for (in = 0; in < N; in++)
    {
      for (i=0; i < CSIZE; i++)
      {
        for (j=0; j < CSIZE; j++)
        {
          c0 = pixmap[im][CSIZE-j-1][i];
          c1 = pixmap[in][CSIZE-j-1][i];

          i0 = SCALE*PAD+im*MBLKSIZE+i*SCALE;
          j0 = SCALE*PAD+in*NBLKSIZE+j*SCALE;

          for (is=0; is < SCALE; is++)
          {
            for (js=0; js < SCALE; js++)
            {
              ip = i0+is;
              jp = j0+js;

              if (ip < 0 || ip >= MSIZE || jp < 0 || jp >= NSIZE)
              {
                printf("Invalid ip = %d, jp = %d\n", ip, jp);
                exit(1);
              }

              if (c0 == 'X')
              {
                array[ip][jp] = 1.0;
              }

              ip = ip+CBLKSIZE;

              if (ip < 0 || ip >= MSIZE || jp < 0 || jp >= NSIZE)
              {
                printf("Invalid ip = %d, jp = %d\n", ip, jp);
                exit(1);
              }

              if (c1 == 'X')
              {
                array[ip][jp] = 1.0;
              } 
            }
          }
        }
      }
    }
  }          

/* This commented out code will produce a test PGM file for debugging

  fp = fopen("test.pgm", "w");
  
  fprintf(fp, "P2\n");
  fprintf(fp, "# Written by mkinfile\n");
  fprintf(fp, "%d %d\n", MSIZE, NSIZE);
  fprintf(fp, "%d\n", 1);

  for (j=NSIZE-1; j >=0; j--)
  {
    for (i=0; i < MSIZE; i++)
    {
      if (array[i][j] == 1)
      {
        fprintf(fp," 1");
      }
      else
      {
        fprintf(fp," 0");
      }
    }
    fprintf(fp,"\n");
  }
*/

  sprintf(filename, "cinput%07dx%07d.dat", MSIZE, NSIZE);

  printf("C filename is <%s>\n", filename);

  fp = fopen(filename, "w");
  
  /* Loop used as for very large arrays single write can fail */
  for (i=0; i < MSIZE; i++){
    fwrite(&array[i][0], sizeof(float), NSIZE, fp); /* C Style */
  }
  fclose(fp);

  sprintf(filename, "finput%07dx%07d.dat", MSIZE, NSIZE);

  printf("F filename is <%s>\n", filename);

  fp = fopen(filename, "w");
  
  for (j=0; j < NSIZE; j++)   /* F Style */
    {
      for (i=0; i < MSIZE; i++)
	{
	  fwrite(&array[i][j], sizeof(float), 1, fp);
	}
    }

  fclose(fp);
}
