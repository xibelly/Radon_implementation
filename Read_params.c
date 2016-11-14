int nread;
#define EXIT_ERROR printf("Error in parameter %s in parameter file\n",buf1); exit(0);


int Param_SPM(char *param_file)
{
  
  int aux_int;
  char buf[200],buf1[200],buf2[200];
  FILE *par_pf=NULL;
  
  printf("reading parameter file in %s\n",param_file);
  
  if(NULL==(par_pf=fopen(param_file,"r")))
    {
      printf("Parameterfile %s not found...\n",param_file);
      exit(0);
    }
  
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    gridx=atoi(buf2);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    gridz=atoi(buf2);
  }
  
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    NRAYS=atoi(buf2);
  }
  //SKIP;
    
  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    GRADIENTE=atof(buf2);
  }
  // SKIP;
      

 fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    MAX_ITERATIONS=atoi(buf2);
  }
  //SKIP;

 fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    Xini=atof(buf2);
  }

 fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    Zini=atof(buf2);
  }

 fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
   Xfin=atof(buf2);
  }

 fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    Zfin=atof(buf2);
  }


  fgets(buf,200,par_pf);
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    EXIT_ERROR;
  }
  else{
    r=atoi(buf2);
  }

    
  fclose(par_pf);

  return 0;
}
