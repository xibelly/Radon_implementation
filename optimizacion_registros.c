/* Xibelly Eliseth Mosquera Escobar
 11) Se menciono en clase que hacer un buen uso de los registros puede ser conveniente para efectos de 
optimizar el usos del tiempo de computo.
Escriba un programa que a través de una secuencia de iteraciones muestre que en efecto hacer la optimización 

t1 = sin(alpha) * x + cos(alpha) * y;
t2 = ­cos(alsph) * x + sin(alpha) * y;
s = sin(alpha); c = cos(alpha);
t1 = s * x + c * y;
t2 = ­c * x + s * y

conduce a una ganancia de tiempo de cpu.
Haga una gráfica de tiempo de ejecución como función del numero de iteraciones y úsela para concluir 
para cuantas iteraciones se vuelve importante considerar la optimización mencionada
 */

/*
Analisis y diseño
Para resolver nuestro problema necesitamos:
-recivir el # en grados y pasarlo a radianes, para ello lo resiviremos por linea de comandos
-pasar el numero de iteraciones (N), de nuevo se hara por linea de comandos
-hacer uso de la funcion clock() para calcular el tiempo de ejecucion
-imprimir en disco el tiempo de ejecucion y el # de iteraciones
-realizar grafica de tiempo vs iteraciones para concluir cuando se vuelve importante la optimizacion
*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

///////////////////////Variables globales/////////////////////
#define pi 3.141592

///////////////////////Funciones/////////////////////////////
funcion1(double alpha, int Ntot, double x, double y)
{
  int i,j;
  double t1, t2;
  clock_t tini, tend, tacum;
  double cpu_time_used;
  
  FILE *pf = NULL;

  pf = fopen("iteraciones_modo1.dat","a");

  for(j=10; j<=Ntot; j*=10)
    {
      tini = clock();
      
      for(i=0; i<j; i++)
	{
	  t1 = sin(alpha)*x + cos(alpha)*y;
	  t2 = -cos(alpha)*x + sin(alpha)*y;  
	  
	}
      tend = clock();
      
      cpu_time_used = ((double) (tend - tini)) / CLOCKS_PER_SEC;
      
      fprintf(pf,"%d %16.8lf\n", j, cpu_time_used);
      
      
    }

    
  fclose(pf);
  
}  

funcion2(double alpha, int Ntot, double x, double y)
{
  double t1, t2, s, c;
  int i,j;
  clock_t tini, tend, tacum;
  double cpu_time_used;

  FILE *pf = NULL;

  pf = fopen("iteraciones_modo2.dat","a");

  s  = sin(alpha);
  c  = cos(alpha);

  for(j=10; j<=Ntot; j*=10)
    {
      
      tini = clock();
      
      for(i=0; i<j; i++)
	{
	  t1 = s*x + c*y;
	  t2 = -c*x + s*y;  
      
	}
      tend = clock();
      
      cpu_time_used = ((double) (tend - tini)) / CLOCKS_PER_SEC;
      fprintf(pf,"%d %16.8lf\n", j, cpu_time_used);
      
      
    }
  fclose(pf);
}  

//////////////////////////////////Programa principal/////////////////////
int main(int argc, char **argv){

  int N,i;
  double angulo, radianes, x, y;

  printf("%d\n",argc);

  if(argc != 5)
    {
      printf("ERROR--> use as:\n");
      printf("%s angulo #iteraciones x y\n",argv[0]);
      exit(0);  
    }
  
  angulo   = atof(argv[1]);
  radianes = angulo*pi/180.0;
  N        = atoi(argv[2]);
  x        = atoi(argv[3]);
  y        = atoi(argv[4]);
  
  printf("%lf %d %lf %lf\n", angulo, N, x, y);
  
  funcion1(radianes, N, x, y);

  funcion2(radianes, N, x, y);
  
  return 0;
  
}
