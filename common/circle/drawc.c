/*=========================================================================
 * Daniel Greenhoe
 * Make LaTeX circle drawing commands
 * http://www.geocities.com/dgreenhoe/
 * National Chiao-Tung University
 *=========================================================================*/

/*-------------------------------------
 * Includes
 *-------------------------------------*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<time.h>
#include<sys\timeb.h>


/*-------------------------------------
 * Defines
 *-------------------------------------*/
#define PI 3.14159265358979323846

/*-------------------------------------
 * Prototypes
 *-------------------------------------*/
int      main(int argc, char *argv[] );
char *gen_title(char *title, int npts, char *prog);


/*-------------------------------------------------------------------------
 * Main
 *-------------------------------------------------------------------------*/
int main(int argc, char *argv[] )
{
   double const r = 100;

   double theta;
   int i;
   double x,y;
   int    npts;
   double delta;
   char   title[1024];


   if(argc>1) npts=atoi(argv[1]);
   else       npts=64; 
   delta = 2.*PI/(double)npts;
   gen_title(title, npts, argv[0]);
   printf("%s",title);
   printf("\n");

   for(i=0,theta=0; i<npts; i++, theta+=delta){
      x = r*cos(theta);
      y = r*sin(theta);
      printf("\\put(%+8.3lf, %+8.3lf)\{\\makebox(0,0)\{\$\\cdot\$\}\}  %% theta=PI*%6.4lf\n",x,y,theta/PI);
      }
   return 0;
}


/*-------------------------------------------------------------------------
 * Generate title
 *-------------------------------------------------------------------------*/
char *gen_title(char *title, int npts, char *prog)
{
   time_t  current_time;
   
   time(&current_time);
   
   sprintf(title,
          "%% -----------------------------------------\n"
          "%%| Daniel Greenhoe                         \n"
          "%%| http://www.geocities.com/dgreenhoe/     \n"
          "%%| LaTeX Drawing commands for circle       \n"
          "%%| number of points = %d                   \n"
          "%%| program: %s                             \n"
          "%%| %s"
          "%% -----------------------------------------\n"
          ,npts,prog,ctime(&current_time));

   return title;
}

