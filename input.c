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

int read_file1(char *filename, int num_lineas)
{
  int i, k,nread ;
  double X, Z; //distance and depth               

  FILE *pf=NULL;
  
  data.ds_x = (double *) malloc(num_lineas*sizeof(double));

  data.ds_z = (double *) malloc(num_lineas*sizeof(double));
  
  data.ds = (double *) malloc(num_lineas*sizeof(double));

  printf("READING FILE: %s\n", filename);
  
  pf = fopen(filename,"r"); 

  
  if(pf==NULL)
    printf("THE FILE CAN NOT BE OPENED\n");
  
  printf(" ->In read_file:\n");
    
  for(i=0; i<num_lineas; i++)
    {
      nread = fscanf(pf,"%lf %lf",&X, &Z);

      data.ds_x[i] = X ; //ray coordinates

      data.ds_z[i] = Z ; 
                 
    }
   
  for(k=0; k<num_lineas; k++)
    {
      data.ds[k] = sqrt( ( (data.ds_x[k+1] - data.ds_x[k])*(data.ds_x[k+1] - data.ds_x[k]) ) + ( (data.ds_z[k+1] - data.ds_z[k])*(data.ds_z[k+1] - data.ds_z[k]) ) );//ray-path -> length of ray in a cell
      
          
      }

  printf("READ FILE STATE %s: SUCESSFUL\n", filename);

  fclose(pf);
  
  return num_lineas;
}

int read_file2(char *filename, int N)
{
  int i, nread;
  double t;

  FILE *pf2=NULL;
  
  printf("READING FILE: %s\n", filename);

  pf2 = fopen(filename,"r"); 
  
  if(pf2==NULL)
    printf("THE FILE CAN NOT BE OPENED\n");
   

  data.mm = (double *) malloc(N *sizeof(double));

  printf(" ->In read_file:\n");
    
  for(i=0; i<N; i++)
    {
      nread = fscanf(pf2,"%lf",&t);

      data.mm[i] = t ; //travel time
                
    }
  printf("READ FILE STATE %s: SUCESSFUL\n", filename);
  
  fclose(pf2);
  
  return 0;
}

int read_file3(char *filename, int N)
{
  int i, nread;
  double v;

  FILE *pf3=NULL;

  printf("READING FILE: %s\n", filename);
  
  pf3 = fopen(filename,"r"); 
  
  if(pf3==NULL)
    printf("THE FILE CAN NOT BE OPENED\n");
  

  data.slowness = (double *) malloc(N *sizeof(double));

  printf(" ->In read_file:\n");
    
  for(i=0; i<N; i++)
    {
      nread = fscanf(pf3,"%lf",&v);

      data.slowness[i] = v ; // 1/c -initial model-
      
    }

  printf("READ FILE STATE %s: SUCESSFUL\n", filename);

  fclose(pf3);
  
  return 0;
}
