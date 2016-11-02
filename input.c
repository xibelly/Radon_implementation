/*Xibelly Eliseth Mosquera Escobar
 * 
 * programa: input.c
 * 
 * 
 * Se encarga de leer los datos de
 * un archivo dado le longitud arbitraria
 * usando malloc.
 * 
 * 
 * La funcion "read_file" recibe como entradas
 * la longitud del archivo y el nombre del archivo
 * a leer. Ambas entradas se hacen usando 
 * la variable ARGS en el makefile incluido con
 * el programa.
 *
 * 
 
 */



int read_file1(char *filename, int N)
{
  int i, nread;
  double X;

  FILE *pf = NULL;
  
  if(pf==NULL)
    printf("THE FILE CAN NOT BE OPENED\n");

  pf = fopen(filename,"r"); 
  
  data.ds = (double *) malloc(N *sizeof(double));
  
  printf(" ->In read_file:\n");
    
  for(i=0; i<N; i++)
    {
      nread = fscanf(pf,"%lf",&X);

      data.ds[i] = X ; //ray-path
      
    }

  printf("READ FILE STATE %s: SUCESSFUL\n", filename);
  
  return 0;
}

int read_file2(char *filename, int N)
{
  int i, nread;
  double Y;

  FILE *pf = NULL;
  
  if(pf==NULL)
    printf("THE FILE CAN NOT BE OPENED\n");

  pf = fopen(filename,"r"); 
  
  data.mm = (double *) malloc(N *sizeof(double));

  printf(" ->In read_file:\n");
    
  for(i=0; i<N; i++)
    {
      nread = fscanf(pf,"%lf",&Y);

      data.mm[i] = Y ; //travel time
      
    }

  printf("READ FILE STATE %s: SUCESSFUL\n", filename);
  
  return 0;
}

int read_file3(char *filename, int N)
{
  int i, nread;
  double Z;

  FILE *pf = NULL;
  
  if(pf==NULL)
    printf("THE FILE CAN NOT BE OPENED\n");

  pf = fopen(filename,"r"); 
  
  data.slowness = (double *) malloc(N *sizeof(double));

  printf(" ->In read_file:\n");
    
  for(i=0; i<N; i++)
    {
      nread = fscanf(pf,"%lf",&Z);

      data.slowness[i] = Z ; // 1/c -initial model-
      
    }

  printf("READ FILE STATE %s: SUCESSFUL\n", filename);
  
  return 0;
}
