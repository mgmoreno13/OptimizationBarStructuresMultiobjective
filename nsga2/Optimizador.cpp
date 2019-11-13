//*****************************
// Rutina Para Optimizar
//*****************************
#pragma hdrstop
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <limits>
//#include <vector>
#include <algorithm> 
#include "datos.h"
#include "estatico.h"
#include "ValoresPropios.h"
#include "fuerzas.h"
#include "main.h"
#include "posprogid.h"
#include "principal.h"
#include "raros.h"
#include "solver.h"
#include "memoria.h"
#include "Optimizador.h"

std::vector<Individuo> poblacion;
std::vector<Individuo> selec;
std::vector<Individuo> hijos;
std::vector<Individuo> newpob;
std::vector<Individuo> R_pob;
//std::vector<Individuo> NP;
//std::vector<Individuo> Elegibles;
//std::vector<Individuo> Penalizados;
// vecino;
int dimension;
int tam_pob=100;
std::vector<int> max_secciones;
int evaluaciones;

int secciones_max;
std::vector<std::vector<int> > Grafo;
std::vector<std::vector<int> > Grafobarras;
std::vector<std::vector<int> > MatrixAdy;
std::vector<bool> visitados;
std::vector<int> rama;
std::vector<int> bars;
std::vector<int> barras_tent;
Individuo vecino;
int cont_bfs=0;
int num1;
int objtosort;

//************************************************

int Optimizador( int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                  int nwktl,int prueba )
{
/**********************
Optimizador
**********************/
srand(prueba); //semilla 

Grafo.resize(npnod);
Grafobarras.resize(nelem);


MatrixAdy.assign(npnod,std::vector<int> (npnod,-1));

for(int ielem=1; ielem<=nelem; ielem++) {
    
    MatrixAdy[lnods[ielem][1]-1][lnods[ielem][2]-1]=ielem-1;//Llena matriz de adyacencias
    MatrixAdy[lnods[ielem][2]-1][lnods[ielem][1]-1]=ielem-1;
    Grafo[lnods[ielem][1]-1].push_back(lnods[ielem][2]-1);
    Grafo[lnods[ielem][2]-1].push_back(lnods[ielem][1]-1);

}

for(int ielem=1; ielem<=nelem; ielem++) {//elementos=cantidad de barras
	
	//printf("%d\n", (int) MatrixAdy[lnods[ielem][1]-1].size());

	for (int i = 0; i < MatrixAdy[lnods[ielem][1]-1].size(); ++i)
	{
		//printf(" m %d\n", MatrixAdy[lnods[ielem][1]-1][i]);
		if (MatrixAdy[lnods[ielem][1]-1][i]!=-1 && MatrixAdy[lnods[ielem][1]-1][i]!=ielem-1)
		{
			Grafobarras[ielem-1].push_back( MatrixAdy[lnods[ielem][1]-1][i]);
			//Grafobarras[MatrixAdy[lnods[ielem][1]-1][i]].push_back(ielem-1);
		}
	}

	for (int i = 0; i < MatrixAdy[lnods[ielem][2]-1].size(); ++i)
	{
		
		if (MatrixAdy[lnods[ielem][2]-1][i]!=-1 && MatrixAdy[lnods[ielem][2]-1][i]!=ielem-1)
		{
			Grafobarras[ielem-1].push_back( MatrixAdy[lnods[ielem][2]-1][i]);
			//Grafobarras[MatrixAdy[lnods[ielem][2]-1][i]].push_back(ielem-1);
		}

	}

 }


/////////////////////////////////////ARCHIVOS DE ESCRITURA/////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
    char buf1[150];
    char buf2[150];

    sprintf(buf1, "FRENTES_%d_%d_nsga_5000",nelem, prueba);
    sprintf(buf2, "INDIVIDUO_%d_%d_nsga_5000",nelem, prueba);

    FILE * p1File;
    FILE * p2File;

    p1File = fopen (buf1,"w");
    p2File = fopen (buf2,"w");
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
   

evaluaciones=0;
bool Optimiza = true;
bool generaMatRigidez;
int ival, casoCarga, jj, Num_iter, imats;
double beta, alfa, penaliza, dd, de, deltaf, fact;
//srand(time(NULL)); //semilla 
// Vector de fuerza de cada caso de carga

//cada caso de carga?
    for( casoCarga = 1; casoCarga <= ncarg; casoCarga++ ){
        npesp[casoCarga]=0;
    }

    
    //ncarg=1;
    
for( casoCarga = 1; casoCarga <= ncarg; casoCarga++ ){

// Lee caso de carga y al mismo tiempo ensambla el vector de fuerzas
    
    Fuerzas(nelem,npnod,ndime,ntipo,nevab,nnode,ngdln,iwrit,casoCarga,lnods,matnu,indfu,ntips,
            coord,xlong,props,fuerb,wcarg,carpp,fuepp,vectr,carga,fuemp,Optimiza,
            npesp[casoCarga]);

// Copia el vector de fuerzas del caso de carga, en caso de que se vaya a optimizar
    copia_vector_fuerza( casoCarga, nelem, nevab, carga, vectores_fuerzas );
    //mg vectores fuerzas es una matriz

}
    
    if( iwrit ==1 ) fprintf(fpLog,"Calculo de Fuerzas = %d \n\n",Prob());

// Inicializa valores de vector de materiales

    
Inicializa_Materiales( nmats, ndime, props, matva, cArmaduras.tamCat, cAceroRF.tamCat,// ----------> mg funcion que esta en principal.cpp
                       cAceroRC.tamCat, cConcr.tamCat, ubiCat );

	secciones_max=matva[1];//maxima cantidad de secciones

// Valor inicial de la funcion de optimizacion
//valor[4] = 1.e10;

// Primer iteracion
Actualiza(nmats,matva,matvm,valor,valom);

CargaPropCatOpt(matva,ndime,nmats,props,cArmaduras.arear,cAceroRF.arerf,cAceroRC.arerc,
                cConcr.arecr,ubiCat);
// Analiza los casos de carga
valor[2] = 0.0; // Para guardar el desplazamiento maximo
valor[3] = 0.0; // Para guardar la eficiencia maxima
generaMatRigidez = true;
    
for( casoCarga = 1; casoCarga <= ncarg; casoCarga++ ){
// Escoge el vector de fuerzas del caso de carga
    escoge_vector_fuerza(casoCarga,nelem,nevab,carga,vectores_fuerzas);
// Incluye las fuerzas de peso propio al vector de fuerzas
    incluye_peso_propio(npesp[casoCarga],nelem,ndime,ntipo,lnods,matnu,ntips,coord,
                        xlong,props,vectr,carga);
// Determina el vector de fuerzas de empotramiento perfecto
    vector_fuer_emp_perf(nelem,nevab,fuemp,carga);
    // Incluye las fuerzas de peso propio al vector de fuerzas
    incluye_peso_propio(npesp[casoCarga],nelem,ndime,ntipo,lnods,matnu,ntips,coord,
                        xlong,props,vectr,carga);
   // printf("indso = %d \n",indso);//mg indice del solucionador que se utilizara
    ival = Estatico( nelem,npnod,nevab,ndime,ntipo,ncaso,npres,false,indso,isale,ngdln,nnode,
                     ntotv,nincl,nreso,neqns,nwktl,ishot,lnods,matnu,inpre,iffix,ntips,maxad,
                     nodea,nodpr,leqns,lincl,lreso,react,coord,presc,fixed,despl,aslod,stiff,
                     srmat,props,xlong,vectr,rigid,carpp,fuepp,carga,fuemp,girom,girtm,tempr,
                     fuerc,valor,cArmaduras.arear,cAceroRF.arerf,cAceroRC.arerc,cConcr.arecr,
                     cConcr.dtcon,indfu,fuerb,fuefl,wcarg,xincl,resor,astif,vectp,vectu,vect1,
                     deslo,ubiCat,efi_max,generaMatRigidez,Optimiza);
// Ya se genero una vez la matriz de rigidez
generaMatRigidez = false;
}

    dimension=nmats;


    //@mg 
    max_secciones.resize(dimension);

    Genera_Poblacion_Inicial(ishot, npnod,nelem, npres,  ncarg, ntipo, ngdln, nmats,
                 ndime, iwrit, isale, nreso, nincl, nfami, ncaso,
                  nnode, mnode, nprop, ngaus, ntens, nevab, indso, ntotv,
                 lnods, inpre, nareas, narerf, narerc, narecon,
                 lreso, lincl, ntips, nodpr, matnu, matva, matvm,
                 iffix, indfu, nodea, maxad,mhigh, leqns,
                 coord, carga,carpp, presc, props,
                 fuemp, fuepp, rigid, xmasa, vmatr,
                 girom, girtm, tempr, fuefl, fuerb,
                 srmat, smmat, cArmaduras, cAceroRF,
                 cAceroRC,cConcr, resor, xincl,
                 fuerc, aslod, despl, fixed, react,
                 vecta, vect1, deslo, xlong, valor,
                 valom, angulo, angulg, wcarg, vectr,
                 astif, vectp, vectu, listaCatalogos, ubiCat,
                 vectores_fuerzas, npesp, efi_max, stiff, neqns,
                 nwktl);
    assig_rank_crow(poblacion);//asigna frente y distancia de crowding 


				for (int i = 0; i < poblacion.size(); ++i)
                {
                    fprintf(p1File,"%lf ",poblacion[i].objs[0]);
                    fflush(p1File);
                    fprintf(p1File,"%lf ",poblacion[i].objs[1]);
                    fflush(p1File);
                    fprintf(p1File,"%lf ",poblacion[i].peso);
                    fflush(p1File);
                    fprintf(p1File,"%lf ",poblacion[i].desplazamiento);
                    fflush(p1File);
                    fprintf(p1File,"%lf ",poblacion[i].efi);
                    fflush(p1File);
                    fprintf(p1File,"%lf ",poblacion[i].sum_pen);
                    fflush(p1File);
                    fprintf(p1File,"%d ",poblacion[i].rank);
                    fflush(p1File);



                    fprintf(p1File,"\n");
                    fflush(p1File);

                    for (int s = 0; s< poblacion[i].secciones.size(); ++s)
                    {
                           fprintf(p2File,"%d ",poblacion[i].secciones[s]);
                           fflush(p2File);
                    }  
                    fprintf(p2File,"\n");
                    fflush(p2File);
                }

                    double D=distancia_var(poblacion);
		            fprintf(p1File,"#Generacion: 1 \n");
		            fflush(p1File);
		            fprintf(p1File,"#Distancia: %lf \n",D);
		            fflush(p1File);


    for (int gen = 2; gen<=5000; ++gen)
    {
        
        Seleccion(poblacion);//Obtengo seleccionados

        Crossover_BCR(selec);//Obtengo hijos

        mutacion(ishot, npnod,nelem, npres,  ncarg, ntipo, ngdln, nmats,
                 ndime, iwrit, isale, nreso, nincl, nfami, ncaso,
                  nnode, mnode, nprop, ngaus, ntens, nevab, indso, ntotv,
                 lnods, inpre, nareas, narerf, narerc, narecon,
                 lreso, lincl, ntips, nodpr, matnu, matva, matvm,
                 iffix, indfu, nodea, maxad,mhigh, leqns,
                 coord, carga,carpp, presc, props,
                 fuemp, fuepp, rigid, xmasa, vmatr,
                 girom, girtm, tempr, fuefl, fuerb,
                 srmat, smmat, cArmaduras, cAceroRF,
                 cAceroRC,cConcr, resor, xincl,
                 fuerc, aslod, despl, fixed, react,
                 vecta, vect1, deslo, xlong, valor,
                 valom, angulo, angulg, wcarg, vectr,
                 astif, vectp, vectu, listaCatalogos, ubiCat,
                 vectores_fuerzas, npesp, efi_max, stiff, neqns,
                 nwktl);//modifica hijos y los evalua
                

        merge(poblacion,hijos,R_pob);
                
        poblacion.clear();
                
        hijos.clear();
               
        selec.clear();
               
        poblacion=fill_new_pob(R_pob);

        R_pob.clear();
       	Local_Search_fact(ishot, npnod,nelem, npres,  ncarg, ntipo, ngdln, nmats,
                 ndime, iwrit, isale, nreso, nincl, nfami, ncaso,
                  nnode, mnode, nprop, ngaus, ntens, nevab, indso, ntotv,
                 lnods, inpre, nareas, narerf, narerc, narecon,
                 lreso, lincl, ntips, nodpr, matnu, matva, matvm,
                 iffix, indfu, nodea, maxad,mhigh, leqns,
                 coord, carga,carpp, presc, props,
                 fuemp, fuepp, rigid, xmasa, vmatr,
                 girom, girtm, tempr, fuefl, fuerb,
                 srmat, smmat, cArmaduras, cAceroRF,
                 cAceroRC,cConcr, resor, xincl,
                 fuerc, aslod, despl, fixed, react,
                 vecta, vect1, deslo, xlong, valor,
                 valom, angulo, angulg, wcarg, vectr,
                 astif, vectp, vectu, listaCatalogos, ubiCat,
                 vectores_fuerzas, npesp, efi_max, stiff, neqns,
                 nwktl);
       	D=distancia_var(poblacion);


		if (gen%500==0)
       	{
       		    for (int i = 0; i < poblacion.size(); ++i)
                {
                    fprintf(p1File,"%lf ",poblacion[i].objs[0]);
                    fflush(p1File);
                    fprintf(p1File,"%lf ",poblacion[i].objs[1]);
                    fflush(p1File);
                    fprintf(p1File,"%lf ",poblacion[i].peso);
                    fflush(p1File);
                    fprintf(p1File,"%lf ",poblacion[i].desplazamiento);
                    fflush(p1File);
                    fprintf(p1File,"%lf ",poblacion[i].efi);
                    fflush(p1File);
                    fprintf(p1File,"%lf ",poblacion[i].sum_pen);
                    fflush(p1File);
                    fprintf(p1File,"%d ",poblacion[i].rank);
                    fflush(p1File);



                    fprintf(p1File,"\n");
                    fflush(p1File);

                    for (int s = 0; s< poblacion[i].secciones.size(); ++s)
                    {
                           fprintf(p2File,"%d ",poblacion[i].secciones[s]);
                           fflush(p2File);
                    }  
                    fprintf(p2File,"\n");
                    fflush(p2File);
                
                
                 }
                         
            D=distancia_var(poblacion);
            fprintf(p1File,"#Generacion: %d \n",gen);
            fflush(p1File);
            fprintf(p1File,"#Distancia: %lf \n",D);
            fflush(p1File);
       	}
       
              
    }

    
   
    fclose(p1File);
    fclose(p2File);



return ival;

}
//***********************************************************************************************************
//***********************************************************************************************************

double distancia_var(std::vector<Individuo> Pob){

	int cont=0,acum=0;

	for(int i=0; i<Pob.size(); i++){//ind
		for(int j=0; j<Pob.size(); j++){//itera con los demas ind
			if (i==j)continue;

			acum+=distancia_individuos(Pob[i],Pob[j]); 
			cont++;
			 
		}

	}
	return ((double)acum/(double)cont)*0.5;

}

double distancia_individuos(Individuo a, Individuo b){


	int sum=0;
	for (int i = 0; i < a.secciones.size(); ++i)
	{
		sum+=abs(a.secciones[i]-b.secciones[i]);
	}

	return sum;

}


void Local_Search_fact(int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                  int nwktl){


    double num,n;
    int secc,cont,j,jj,b=0,modif,mejor=0;
    //std::set<int> myset1;
    std::vector<double> v;

    for (int i = 0; i < tam_pob; ++i)//recorre  poblacion
    { 
    		if (poblacion[i].sum_pen == 0)//si el individuo no esta penalizado, no hacer la busqueda de factibilidad
    		{
    			continue;
    		}


    		b=0;
    		//CREA AL VECINO-ANTES DE MODIFICAR
    		for (int k= 0; k < poblacion[i].secciones.size(); ++k)//recorre secciones 
	        {

	            vecino.secciones.push_back(poblacion[i].secciones[k]);
	        }


	       for (int k= 0; k < poblacion[i].valores.size(); ++k)//recorre secciones 
	        {

	            vecino.valores.push_back(poblacion[i].valores[k]);
	        }


    	do{  

								//std::cout<<"vecino size "<<vecino.secciones.size()<<std::endl;
					    	modif=rand()%(poblacion[i].secciones.size());//cual posicion voy a modificar
					    	//printf("%d\n",poblacion[i].secciones.size() );


					    		//std::cout<<"modif "<<modif<<" secciones size "<<poblacion[i].secciones.size()<<" barras size "<<poblacion[i].barras.size()<<std::endl;
						        

						     //if (modif <= vecino.seccio

							for (int k= 0; k < poblacion[i].secciones.size(); ++k)//regreso al individuo original
					        {


					           matva[k+1]=poblacion[i].secciones[k];
                               vecino.secciones[k]=poblacion[i].secciones[k];
					                  // 	std::cout<<" "<<matva[k+1]<<" ";
					        }

					 //std::cout<<" matnu "<<std::endl;
					       /* for (int k= 0; k < poblacion[i].barras.size(); ++k)//recorre secciones 
					        {
					        	matnu[k+1]=vecino.barras[k];
					        	        //	std::cout<<" "<<matnu[k+1]<<" ";
					        }*/
                            do{
                                secc=1+rand()%(secciones_max+1-1);

                            }while(secc==vecino.secciones[modif]);


					        vecino.secciones[modif]=secc;

					        //if (modif+1 < poblacion[i].secciones.size())
					        //{
					        	matva[modif+1]=secc;
					        //}
					        
					          	
					        
					            
					                
					       
						//printf("secciones %d  %d \n",poblacion[i].secciones[0],poblacion[i].secciones[1]);


							v=Funcion_fitness(ishot, npnod,nelem,npres,ncarg,ntipo,ngdln,nmats,ndime,iwrit,isale,nreso,nincl,
					                                                  nfami,ncaso, nnode,mnode,nprop,ngaus,ntens,nevab,indso,ntotv,lnods,inpre,
					                                                  nareas,narerf,narerc,narecon,lreso,lincl,ntips,nodpr,matnu,matva,matvm,
					                                                  iffix,indfu,nodea,maxad,mhigh,leqns,coord,carga,carpp,presc,props,fuemp,
					                                                  fuepp,rigid,xmasa,vmatr, girom,girtm,tempr,fuefl,fuerb,srmat,smmat,
					                                                  cArmaduras,cAceroRF,cAceroRC,cConcr,resor,xincl,fuerc,aslod,despl,fixed,
					                                                  react,vecta,vect1,deslo,xlong,valor,valom,angulo,angulg,wcarg,vectr,astif,
					                                                  vectp,vectu,listaCatalogos,ubiCat,vectores_fuerzas,npesp,efi_max,stiff,
					                                                  neqns, nwktl);//@mg

							vecino.objs.push_back(v[0]);
							vecino.objs.push_back(v[1]);
							vecino.pena.push_back(v[2]);
							vecino.pena.push_back(v[3]);
							vecino.sum_pen=v[2]+v[3];



                            vecino.desplazamiento=valor[2];
                            vecino.peso=valor[4];
                            vecino.efi=valor[3];

					        vecino.valores.clear();
					        v.clear();

					       for (int j = 1; j < 5; ++j)//4 porque son 4 valores, peso propio,desplazamiento maximo,eficiencia maxima
					        {
					           vecino.valores.push_back(valor[j]);
					           // printf("valores %lf ",valor[j] );
					        }

					        //printf("\n");
					     


					        if (vecino.sum_pen < poblacion[i].sum_pen)
					        {
					        	//printf("vecino.sum_pen %lf, poblacion[i].sum_pen %lf\n",vecino.sum_pen,poblacion[i].sum_pen );
					        	actualizar_vecino_poblacion(vecino,i);
					        	//exit(0);
					        }


					        b++;
					       // printf("busqueda local\n");
					        //printf("%d\n",b);


					        vecino.valores.clear();
       						vecino.pena.clear();
       						vecino.objs.clear();
    	}while(b < 75);

        vecino.secciones.clear();
        

    
        
        
    }



}


void actualizar_vecino_poblacion(Individuo vecino,int i){

	poblacion.erase(poblacion.begin()+i);

	poblacion.insert(poblacion.begin()+i,vecino);



}


std::vector<Individuo> fill_new_pob(std::vector<Individuo> R_pob){

    std::vector< std::vector<int> > Frentes;
    //printf("R_pob size %d\n",R_pob.size() );
    Frentes=FNDS(R_pob);
    std::vector<Individuo> Pob_new;
    std::vector<Individuo> Last_front;
    int i=0;
    //printf("1.1\n");
    //printf("Pob_new %d  Frentes %d\n",Pob_new.size(),Frentes.size() );


    while(Pob_new.size()+Frentes[i].size() <= tam_pob){

        crow_distance(Frentes[i],R_pob);
        //printf("1.1.1\n");
        for (int j = 0; j < Frentes[i].size(); ++j)
        {
            //printf("Indice Frentes[i][j] %d y en i %d y j %d \n", Frentes[i][j],i,j);
            Pob_new.push_back(R_pob[Frentes[i][j]]);
        }
        //printf("1.1.2\n");
        i++;
    }
	if(i==0) crow_distance(Frentes[i],R_pob);
            //printf("2.1\n");

    for (int j = 0; j < Frentes[i].size(); ++j)
    {
            Last_front.push_back(R_pob[Frentes[i][j]]);
    }
            //printf("3.1\n");

    std::sort(Last_front.begin(),Last_front.end(),sortbycrow);
    int total=tam_pob-(int)Pob_new.size();
    for (int i = 0; i < total; ++i)
    {
        Pob_new.push_back(Last_front[i]);
    }
    return Pob_new;




}



void merge(std::vector<Individuo> poblacion,std::vector<Individuo> hijos,std::vector<Individuo>& R_pob){

    for (int i = 0; i < tam_pob; ++i)
    {   

        R_pob.push_back(poblacion[i]);
    }
    for (int i = 0; i < tam_pob; ++i)
    {   

        R_pob.push_back(hijos[i]);
    }

    //printf("R_pob size %d\n",R_pob.size() );

}


void mutacion(int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                  int nwktl){

    double num,n,secc;
    std::vector<double> v;

    for (int i = 0; i < tam_pob; ++i)//recorre  
    {   
        
        
        for (int j = 0; j < dimension; ++j)//recorre secciones segun sea viga o columna
        {
             num=rand()/(double)RAND_MAX;
            //cout<<num<<endl;
            if (num<0.055)//antes en 0.055-(1/(double)nelem)
            {
                secc=1+rand()%(max_secciones[j]+1-1);//necesito saber  cual es el numero maximo de secciones en el catalogo 1-59
                hijos[i].secciones[j]=(int)secc;
                    
            }


            matva[j+1]=hijos[i].secciones[j];
                
        }
    //printf("secciones %d  %d \n",hijos[i].secciones[0],hijos[i].secciones[1]);


//////////////////////////////////////////////////////////////////////////////////////
        v=Funcion_fitness(ishot, npnod,nelem,npres,ncarg,ntipo,ngdln,nmats,ndime,iwrit,isale,nreso,nincl,//QUEDE AQUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                                                  nfami,ncaso, nnode,mnode,nprop,ngaus,ntens,nevab,indso,ntotv,lnods,inpre,
                                                  nareas,narerf,narerc,narecon,lreso,lincl,ntips,nodpr,matnu,matva,matvm,
                                                  iffix,indfu,nodea,maxad,mhigh,leqns,coord,carga,carpp,presc,props,fuemp,
                                                  fuepp,rigid,xmasa,vmatr, girom,girtm,tempr,fuefl,fuerb,srmat,smmat,
                                                  cArmaduras,cAceroRF,cAceroRC,cConcr,resor,xincl,fuerc,aslod,despl,fixed,
                                                  react,vecta,vect1,deslo,xlong,valor,valom,angulo,angulg,wcarg,vectr,astif,
                                                  vectp,vectu,listaCatalogos,ubiCat,vectores_fuerzas,npesp,efi_max,stiff,
                                                  neqns, nwktl);//@mg'

          					 hijos[i].objs.push_back(v[0]);
							 hijos[i].objs.push_back(v[1]);
							 hijos[i].pena.push_back(v[2]);
							 hijos[i].pena.push_back(v[3]);
							 hijos[i].sum_pen=v[2]+v[3];

       /* for (int j = 0; j < nmats; ++j)
        {
             hijos[i].secciones[j]=matva[j+1];//haciendo individuo

            //poblacion[i].secciones.push_back(valor[j+3]);
        }*/
       // printf("pob secciones %d  %d\n", poblacion[i].secciones[0],poblacion[i].secciones[1]);
        for (int j = 1; j < 5; ++j)//4 porque son 4 valores, peso propio,desplazamiento maximo,eficiencia maxima
        {
             hijos[i].valores.push_back(valor[j]);//
            //poblacion[i].secciones.push_back(valor[j+3]);
        }
        hijos[i].desplazamiento=valor[2];
        hijos[i].peso=valor[4];
      	 hijos[i].efi=valor[3];

        
        v.clear();
    }




}



void DFS_DT(int v, int padre){//para recorrer grafo GRAFO
    
   // rama.push_back(v);
    if (padre != -1)
    {
        visitados[padre]=true;
        bars.push_back(MatrixAdy[v][padre]);
        cont_bfs++;
        
    }
    
    std::vector<int> aux;
    for (int n= 0; n< Grafo[v].size(); ++n)
    {
    aux.push_back(n);
    }
    
    random_shuffle(aux.begin(),aux.end());
    int j;

    for (int n= 0; n< Grafo[v].size(); ++n)
    {
    j=aux[n];
    
        if (cont_bfs>num1)
        {
            break;
        }
         if (visitados[ Grafo[v][j]] != true)
        {
            DFS_DT( Grafo[v][j],v);
            //break;
            
        }

                
            

    }

}

void BC_barras(int v){//para recorrer grafo GRAFO
	visitados.assign(Grafobarras.size(),false);
	//printf("%d\n",v);
   
	bars.push_back(v);
	cont_bfs++;
	for (int i = 0; i < Grafobarras[v].size(); ++i)
	{
		//printf("tentativos\n");
		barras_tent.push_back(Grafobarras[v][i]);
	}
	
	visitados[v]=true;
	double numrand;

	while(cont_bfs<num1){

		numrand=rand()%(barras_tent.size());
		if (!visitados[barras_tent[numrand]]){

			bars.push_back(barras_tent[numrand]);
			visitados[barras_tent[numrand]]=true;
			//printf("%d\n",barras_tent[numrand]);

			for (int i = 0; i < Grafobarras[barras_tent[numrand]].size(); ++i)
			{
				barras_tent.push_back(Grafobarras[barras_tent[numrand]][i]);
			}
			cont_bfs++;

		}




	}

}

void Crossover_BCR(std::vector<Individuo> selec){
    double point,num;
    

        int aux,unif,n;
        std::vector<int> hijo1;
        std::vector<int> hijo2;


        for (int k = 0; k < tam_pob;k=k+2)//OJO
        {

            num=rand()/(double)RAND_MAX;

            if(num<0.8){//PROBABILIDAD DE CROSSOVER

                hijo1=selec[k].secciones;
                hijo2=selec[k+1].secciones;

                n=rand()%(Grafo.size());//en que nodo empieza
                num1=rand()%(Grafo.size());//maximo numero de nodos  a recorrer
                bars.clear();//barras conexas
                cont_bfs=0;//cuenta los nodos que se van abriendo
                visitados.clear();//reinicia a los visitados.
                //visitados.assign(Grafo.size(),false);
                //DFS_DT(n,-1);
                barras_tent.clear();
				BC_barras(n);
                

                for (int i = 0; i < bars.size(); ++i)
                {
                    hijo1[bars[i]]=selec[k+1].secciones[bars[i]];
                    hijo2[bars[i]]=selec[k].secciones[bars[i]];
                }

                Individuo Indi1;
                Individuo Indi2;

                hijos.push_back(Indi1);
                hijos.push_back(Indi2);
                for (int l = 0; l < dimension; ++l)
                {
                    hijos[k].secciones.push_back(hijo1[l]);
                    hijos[k+1].secciones.push_back(hijo2[l]);
                }
                //hijos.push_back(hijo1);
                //hijos.push_back(hijo2);
                hijo1.clear();
                hijo2.clear();
            


            }else{
                Individuo Indi1;
                Individuo Indi2;


                hijos.push_back(Indi1);
                hijos.push_back(Indi2);
                for (int l = 0; l < dimension; ++l)
                {
                    hijos[k].secciones.push_back(selec[k].secciones[l]);
                    hijos[k+1].secciones.push_back(selec[k+1].secciones[l]);
                }

            }

        }   


}



void Seleccion(std::vector<Individuo> pob){
	double num;
    std::vector<int> pseudopob;
    for (int i = 0; i < tam_pob; ++i)
    {
        pseudopob.push_back(i);
    }

    std::random_shuffle(pseudopob.begin(),pseudopob.end() );//posiciciones aleatorias en la poblacion

    for (int i = 0; i < tam_pob; ++i)
    {
	num=rand()/(double)RAND_MAX;
        if (pob[i].rank < pob[pseudopob[i]].rank)
        {
            selec.push_back(pob[i]);

        }else if(pob[i].rank > pob[pseudopob[i]].rank){

            selec.push_back(pob[pseudopob[i]]);

        }else //if (pob[i].rank == pob[pseudopob[i]].rank)//si pertenecen al mismo frente comparar por distancia de crowding
        {
            if (pob[i].dist_crow<pob[pseudopob[i]].dist_crow)
            {
                selec.push_back(pob[pseudopob[i]]);
            }
            else if(pob[i].dist_crow>pob[pseudopob[i]].dist_crow){
                selec.push_back(pob[i]);
            }else if(num <= 0.5){
		selec.push_back(pob[pseudopob[i]]);
	    }else{
		selec.push_back(pob[i]);
		
                  }
        }
    }




}

void assig_rank_crow(std::vector<Individuo>& poblacionactual){


    std::vector< std::vector<int> > Frentes;
    Frentes=FNDS(poblacionactual);

    for (int i = 0; i < Frentes.size(); ++i)
    {
        crow_distance(Frentes[i],poblacionactual);
    }

    
}

bool sortbyobjs(Individuo a, Individuo b){


    return b.objs[objtosort] > a.objs[objtosort];

}


bool sortbycrow(Individuo a, Individuo b){


    return b.dist_crow < a.dist_crow;

}

//Fast Non-Dominated-Sort
std::vector< std::vector<int> > FNDS(std::vector<Individuo>& v){



			int bb=0;
			int cc=0;
			//cout<<Puntos2.size()<<endl;
			std::vector<int> aux;
			std::vector< std::vector<int> > S;
			std::vector< std::vector<int> > Frentes;
  			std::vector<int> n;
			std::vector<int> Q;
			S.resize(v.size());//puntos que el individuo domina
			n.resize(v.size());//cantidad de puntos que dominan al individuo
			Frentes.push_back(aux);

//			cout<<" puntos "<<v.size()<<endl;

	for (int i = 0; i < v.size(); ++i)
	{
		n[i]=0;
		S[i].clear();
		//cout<<v.size()<<endl;
		for (int j = 0; j< v.size(); ++j)
		{
			//cont++;
			if (i==j) continue;
			bb=0;
			cc=0;
			//cout<<v[i].coor.size()<<endl;
			for (int k = 0; k < v[i].objs.size(); ++k)
			{	
				
				if (v[i].objs[k] > v[j].objs[k])
				{
					bb=1;
				}else if (v[j].objs[k] > v[i].objs[k])
				{
					cc=1;
				}
			}
			if (bb == 0)
			{
				//cout<<"S i "<<i<<" j "<<j<<endl;
				S[i].push_back(j);


			}else if (cc == 0)
			{

				n[i]=n[i]+1;//numero de puntos q lo dominan
				//cout<<" n"<<n[i]<<" i "<<i<<endl;
			}



		}
		//cout<<v[0].S.size()<<" v"<<endl;
		//Frentes.push_back(aux);
		if (n[i] == 0)
		{
			v[i].rank=0;
			Frentes[0].push_back(i);
		}
	//cout<<" fre "<<Frentes[0][0].coor[0]<<" "<<Frentes[0][0].coor[1]<<endl;
	
	}
	//cout<<" fre "<<Frentes[0][0].S[0].n<<endl;
	//cout<<Frentes.size()<<endl;


	int k=0;

	//cout<<"Numero de frentes  VARIABLEEE "<<num_frentes<<endl;
	do{
	
		Q.clear();

//		cout<<" Frentes[k].size() "<<Frentes[k].size()<<endl;
		for (int p = 0; p < Frentes[k].size(); ++p)
		{
			//cont1++;
			for (int q = 0; q < S[Frentes[k][p]].size(); ++q)
			{
				//cont2++;
				n[S[Frentes[k][p]][q]]=n[S[Frentes[k][p]][q]]-1;
				if (n[S[Frentes[k][p]][q]] == 0)
				{
					v[S[Frentes[k][p]][q]].rank=k+1;
					
					Q.push_back(S[Frentes[k][p]][q]);
					//cout<<" Q size "<<Q.size()<<" para el frente "<<k<<endl;
				}
			}
		}

		if(Q.size()==0) break;
		k++;
		Frentes.push_back(aux);
		Frentes[k]=Q;

	
	}while(Frentes[k].size() != 0);
	//printf("Numero de frentes  %d \n",(int)Frentes.size());
	//cout<<"Numero de frentes "<<Frentes.size()<<endl;



        //printf("%d\n",(int)v[0].objs.size() );

	
		/*for (int i = 0; i <Frentes.size(); ++i)
		{
            for (int h = 0; h < Frentes[i].size(); ++h)
            {

			printf("Frente num %d \n",i);
			for (int j = 0; j < v[i].objs.size(); ++j)
			{
				//cout<<" Frente num "<<i<<" - ";
				printf("%lf ",v[Frentes[i][h]].objs[j]);
			}
			printf("\n");
            }
		}*/


return Frentes;//en el frente lo que se escribe es el numero donde se ubica el punto en la lista original V de puntos
//gguarda los indices de los puntos que se encuentran en cada frente

}

void crow_distance(std::vector<int> frente,std::vector<Individuo>& pob){
	//std::vector<int> frente- se refiere a que guarda los indices de los puntos que se encuentran en cada frente de la poblacion pob
    //std::vector<Individuo> pob -el frente pertenece a la poblacion pob
    float inf;
    double max,min;
    inf= std::numeric_limits<float>::infinity();
	int l=frente.size();
	std::vector<Individuo> frente_aux;
	for (int i = 0; i < frente.size(); ++i)
	{
		pob[frente[i]].dist_crow=0;//inicializo distancias de todos los puntos de mi frente actual
		frente_aux.push_back(pob[frente[i]]);
        	frente_aux[i].dist_crow=0;
        	frente_aux[i].indice=frente[i];
	}

   
	for (int i = 0; i < frente_aux[0].objs.size(); ++i)//recorre cada objetivo
	{

        objtosort=i;//objetivo por el que se ordenara
        sort(frente_aux.begin(),frente_aux.end(),sortbyobjs);
        frente_aux[0].dist_crow=inf;
        frente_aux[frente_aux.size()-1].dist_crow=inf;

        max= frente_aux[frente_aux.size()-1].objs[i];
        min=frente_aux[0].objs[i];

        for (int j = 1; j < frente_aux.size()-1; ++j)//recorre frente excepto los extremos
        {
            frente_aux[j].dist_crow=frente_aux[j].dist_crow+((frente_aux[j+1].objs[i]-frente_aux[j-1].objs[i])/(max-min));
        }

	}

    for (int i = 0; i < frente_aux.size(); ++i)
    {
        pob[frente_aux[i].indice].dist_crow=frente_aux[i].dist_crow;
    }


 

}






void Genera_Poblacion_Inicial(int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixed, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                 int nwktl){//@mg 

    //int tam_pob=10;
    poblacion.resize(tam_pob); 
    std::vector<double> v; 

    
    //poblacion.resize(1);
//
    for (int i = 0; i < tam_pob; ++i)
    //for (int i = 0; i < 1; ++i)
    {
       // printf("\nINDIVIDUO %d\n", i);
        Genera_Nuevo_Individuo(nmats, ndime, matva, props, cArmaduras.tamCat,
                           cAceroRF.tamCat, cAceroRC.tamCat, cConcr.tamCat, ubiCat);

        
        //matva[1]=8;
        //matva[2]=2;


       v=Funcion_fitness(ishot, npnod,nelem,npres,ncarg,ntipo,ngdln,nmats,ndime,iwrit,isale,nreso,nincl,//QUEDE AQUIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
                                                  nfami,ncaso, nnode,mnode,nprop,ngaus,ntens,nevab,indso,ntotv,lnods,inpre,
                                                  nareas,narerf,narerc,narecon,lreso,lincl,ntips,nodpr,matnu,matva,matvm,
                                                  iffix,indfu,nodea,maxad,mhigh,leqns,coord,carga,carpp,presc,props,fuemp,
                                                  fuepp,rigid,xmasa,vmatr, girom,girtm,tempr,fuefl,fuerb,srmat,smmat,
                                                  cArmaduras,cAceroRF,cAceroRC,cConcr,resor,xincl,fuerc,aslod,despl,fixed,
                                                  react,vecta,vect1,deslo,xlong,valor,valom,angulo,angulg,wcarg,vectr,astif,
                                                  vectp,vectu,listaCatalogos,ubiCat,vectores_fuerzas,npesp,efi_max,stiff,
                                                  neqns, nwktl);//@mg'

        poblacion[i].objs.push_back(v[0]);
		poblacion[i].objs.push_back(v[1]);
		poblacion[i].pena.push_back(v[2]);
		poblacion[i].pena.push_back(v[3]);
		poblacion[i].sum_pen=v[2]+v[3];


        for (int j = 0; j < nmats; ++j)
        {
            poblacion[i].secciones.push_back(matva[j+1]);//haciendo individuo

            //poblacion[i].secciones.push_back(valor[j+3]);
        }
       // printf("pob secciones %d  %d\n", poblacion[i].secciones[0],poblacion[i].secciones[1]);
        for (int j = 1; j < 5; ++j)//4 porque son 4 valores, peso propio,desplazamiento maximo,eficiencia maxima
        {
            poblacion[i].valores.push_back(valor[j]);//
            //poblacion[i].secciones.push_back(valor[j+3]);
        }
        poblacion[i].desplazamiento=valor[2];
        poblacion[i].peso=valor[4];
        poblacion[i].efi=valor[3];

        v.clear();
    }

}


std::vector<double> Funcion_fitness( int ishot, int npnod, int nelem, int npres, int ncarg, int ntipo, int ngdln, int nmats,
                 int ndime, int iwrit, int isale, int nreso, int nincl, int nfami, int ncaso,
                 int nnode, int mnode, int nprop, int ngaus, int ntens, int nevab, int indso, int ntotv,
                 int** lnods, int** inpre, int* nareas, int* narerf, int* narerc, int* narecon,
                 int* lreso, int* lincl,int* ntips, int* nodpr, int* matnu, int* matva, int* matvm,
                 int* iffix, int** indfu, int* nodea, int* maxad, int* mhigh, int** leqns,
                 double** coord, double** carga, double** carpp, double** presc, double** props,
                 double** fuemp, double** fuepp, double** rigid, double** xmasa, double** vmatr,
                 double** girom, double** girtm, double** tempr, double** fuefl,double** fuerb,
                 double*** srmat,double*** smmat, cat_armadura &cArmaduras, cat_acero_RF &cAceroRF,
                 cat_acero_RC &cAceroRC, cat_concreto &cConcr, double** resor, double* xincl,
                 double* fuerc, double* aslod, double* despl, double* fixd, double* react,
                 double* vecta, double* vect1, double* deslo, double* xlong, double* valor,
                 double* valom, double* angulo, double** angulg, double* wcarg, double** vectr,
                 double** astif, double* vectp, double* vectu, char **listaCatalogos, int* ubiCat,
                 double ***vectores_fuerzas, int *npesp, double *efi_max, double *stiff, int neqns,
                  int nwktl){//@mg 






    double dd,de,penalizacion1=0.0,penalizacion2=0.0;
    std::vector<double> fitness;
    double penaliza = 1.e7;
   // bool Optimiza = true;
    bool generaMatRigidez,Optimiza=true;
    int ival, casoCarga;

    
   
    CargaPropCatOpt(matva,ndime,nmats,props,cArmaduras.arear,cAceroRF.arerf,cAceroRC.arerc,cConcr.arecr,ubiCat);

        valor[2] = 0.0; // Para guardar el desplazamiento maximo
        valor[3] = 0.0; // Para guardar la eficiencia maxima
        generaMatRigidez = true;


        for( casoCarga = 1; casoCarga <= ncarg; casoCarga++ ){
            // Escoge el vector de fuerzas del caso de carga
          escoge_vector_fuerza(casoCarga,nelem,nevab,carga,vectores_fuerzas);
            // Incluye las fuerzas de peso propio al vector de fuerzas

          incluye_peso_propio(npesp[casoCarga],nelem,ndime,ntipo,lnods,matnu,ntips,coord,
                             xlong,props,vectr,carga);


            // Determina el vector de fueras de empotramiento perfecto
            //      vector_fuer_emp_perf(nelem,nevab,fuemp,carga);
          ival = Estatico( nelem,npnod,nevab,ndime,ntipo,ncaso,npres,false,indso,isale,ngdln,nnode,
                     ntotv,nincl,nreso,neqns,nwktl,ishot,lnods,matnu,inpre,iffix,ntips,maxad,
                     nodea,nodpr,leqns,lincl,lreso,react,coord,presc,fixd,despl,aslod,stiff,
                     srmat,props,xlong,vectr,rigid,carpp,fuepp,carga,fuemp,girom,girtm,tempr,
                     fuerc,valor,cArmaduras.arear,cAceroRF.arerf,cAceroRC.arerc,cConcr.arecr,
                     cConcr.dtcon,indfu,fuerb,fuefl,wcarg,xincl,resor,astif,vectp,vectu,vect1,
                     deslo,ubiCat,efi_max,generaMatRigidez,Optimiza);
            // Ya se genero una vez la matriz de rigidez
          generaMatRigidez = false;
        }

    dd = valor[2]-10.0;  // Castiga desplazamientos mayores a 10cm. Aparentemente aqui podemos poner claro/norma.
    de = valor[3]-1.0;   // Castiga eficiencia maxima. Ahorita la trata de llevar al 100%.
    valor[4]=valor[1];   // Peso de la estructura total.
            
            if (de >0) penalizacion1=de*penaliza;
            if (dd >0) penalizacion2=dd*penaliza;
   // fitness= valor[4];
    fitness.push_back(valor[4]+penalizacion1+penalizacion2);// mg se agrega al valor de la funcion //PESO 
    fitness.push_back(valor[2]+penalizacion1+penalizacion2);// mg se agrega al valor de la funcion  //DESPLAZAMIENTO
    fitness.push_back(penalizacion1);
    fitness.push_back(penalizacion2);
    fitness.push_back(valor[2]);
    fitness.push_back(valor[3]);
    fitness.push_back(valor[4]);
    evaluaciones++;
    
   /* printf(" En Funcion_fitness \n ",fitness);

    for(int i=1;i<=nelem;i++){
        printf(" %d ",matva[i]);
    }

     printf(" %lf \n ",fitness);*/
  

    return fitness;//DEVUELVE: obj1 obj2 penalizacion1 penalizacion 2


}



void Genera_Nuevo_Individuo( int nmats, int ndime, int* matva, double** props,//@mg 
                             int *Arm_tamCat, int *RF_tamCat, int *RC_tamCat,
                             int *Con_tamCat, int *ubiCat){


    int rolado, iCat;
    // Escoge un material al azar de la lista de materiales. Entre 1 y nmats.
    //mg otra opcion seria variable = limite_inferior + rand() % (limite_superior +1 - limite_inferior) ;
    //imats=(int)(((float) rand() / (float) RAND_MAX)*(float)nmats-.000000000001)+1;//lo que elige aqui es la columna o viga 
    /*printf("nMATS: %d\n",nmats ); //mg lo comento porque lo que se estaba generando era un vecino, ahora quiero es una poblacion
    
    printf("tam max_secciones %d\n", max_secciones.size() );*/
    //printf("ndime: %d\n",ndime );
    for (int imats = 1; imats <= nmats; ++imats)
    {
       // printf("iMATS: %d\n",imats );
        //printf("1 %lf\n",props[imats][5] );
        //printf("2 %lf\n",props[imats][8] );
        //if (ndime ==2 ) rolado= (int)(props[imats][5]+0.5);
    rolado= (int)(props[imats][8]+0.5);
        // printf("rolado: %d\n",rolado );

    // Del material antes escogido, se le cambia al azar la seccion de las
    // disponibles en el catalogo.
        iCat = ubiCat[imats];

        switch( rolado ){
        case 0://tamaño del catalogo para ese material "imats" en especifico.
        matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)Arm_tamCat[iCat]-.000000000001)+1;
            max_secciones[imats-1]=Arm_tamCat[iCat];
        //printf("------------- matva %d \n",matva[imats]);
            //printf(" Arm %d\n",Arm_tamCat[iCat]);
        break;
        case 1:
        matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)RF_tamCat[iCat]-.000000000001)+1;
             max_secciones[imats-1]=RF_tamCat[iCat];
        //printf("------------- matva %d \n",matva[imats]);
            //printf("RF %d\n",RF_tamCat[iCat]);
        break;
        case 2:
        matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)RC_tamCat[iCat]-.000000000001)+1;
            max_secciones[imats-1]=RC_tamCat[iCat];
            //printf(" RC %d\n",RC_tamCat[iCat]);
        //printf("------------- matva %d \n",matva[imats]);
        break;
        case 3:
        matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)Con_tamCat[iCat]-.000000000001)+1;
        //printf("------------- matva %d \n",matva[imats]);
            max_secciones[imats-1]=Con_tamCat[iCat];
            //printf("cat_concreto %d\n",Con_tamCat[iCat]);
        break;
        default:
        break;
        }

    }

    
    /*printf("matva1: %d matva2: %d \n",matva[1],matva[2]);
    printf("max_secciones: %d max_secciones: %d \n",max_secciones[0],max_secciones[1]);*/


}



/*void Genera_Nueva_Poblacion( int nmats, int ndime, int* matva, double** props,//mg LO QUE EN REALIDAD GENERA ES UN VECINO, escoje un numero al azar de cual tipo de estructura va a modificar
                             int *Arm_tamCat, int *RF_tamCat, int *RC_tamCat,
                             int *Con_tamCat, int *ubiCat )
{
int imats, rolado, iCat;
// Escoge un material al azar de la lista de materiales. Entre 1 y nmats.

//mg otra opcion seria variable = limite_inferior + rand() % (limite_superior +1 - limite_inferior) ;
imats=(int)(((float) rand() / (float) RAND_MAX)*(float)nmats-.000000000001)+1;//lo que elige aqui es cambiar la seccion de la columna o de la viga , imats es viga o columna
printf("IMATS: %d\n",imats );

if (ndime ==2 ) rolado= (int)props[imats][5]+0.5;
if (ndime ==3 ) rolado= (int)props[imats][8]+0.5;

// Del material antes escogido, se le cambia al azar la seccion de las
// disponibles en el catalogo.
iCat = ubiCat[imats];
switch( rolado ){
case 0://tamaño del catalogo para ese material "imats" en especifico.
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)Arm_tamCat[iCat]-.000000000001)+1;
break;
case 1:
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)RF_tamCat[iCat]-.000000000001)+1;
break;
case 2:
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)RC_tamCat[iCat]-.000000000001)+1;
break;
case 3:
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)Con_tamCat[iCat]-.000000000001)+1;
break;
default:
break;
}

printf("matva1: %d matva2: %d \n",matva[1],matva[2]);

}*/



/*void Genera_Nuevo_Individuo( int nmats, int ndime, int* matva, double** props,
                             int *Arm_tamCat, int *RF_tamCat, int *RC_tamCat,
                             int *Con_tamCat, int *ubiCat )
{
int imats, rolado, iCat;
// Escoge un material al azar de la lista de materiales. Entre 1 y nmats.

//mg otra opcion seria variable = limite_inferior + rand() % (limite_superior +1 - limite_inferior) ;
imats=(int)(((float) rand() / (float) RAND_MAX)*(float)nmats-.000000000001)+1;//lo que elige aqui 
printf("IMATS: %d\n",imats );

if (ndime ==2 ) rolado= (int)props[imats][5]+0.5;
if (ndime ==3 ) rolado= (int)props[imats][8]+0.5;

// Del material antes escogido, se le cambia al azar la seccion de las
// disponibles en el catalogo.
iCat = ubiCat[imats];//mg tama;o del catalogo?
switch( rolado ){
case 0://tamaño del catalogo para ese material "imats" en especifico.
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)Arm_tamCat[iCat]-.000000000001)+1;
break;
case 1:
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)RF_tamCat[iCat]-.000000000001)+1;
break;
case 2:
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)RC_tamCat[iCat]-.000000000001)+1;
break;
case 3:
matva[imats]= (int)(((float) rand() / (float) RAND_MAX)*(float)Con_tamCat[iCat]-.000000000001)+1;
break;
default:
break;
}

printf("matva1: %d matva2: %d \n",matva[1],matva[2]);

}*/


//******************************************************************************
void Recosido_Simulado(int Num_iter)
{  /*
int jj,imats, ival;
double dd,de,beta,deltaf,fact,penaliza;

penaliza=1.e07;
for(jj=1; jj<=Num_iter; jj++){
for(imats=1; imats<=nmats;imats++)
printf("matva[%d]= %d \t",imats,matva[imats]);
printf("\n");
Genera_Nueva_Poblacion(nmats,ndime,matva,valmax,props);
CargaPropCatOpt(matva);
ival =Estatico();
dd = valor[2]-10.0;
de = valor[3]-1.0;
valor[4]=valor[1];
if (dd >0)valor[4]+=dd*penaliza;
if (de >0)valor[4]+=de*penaliza;
if (beta < 1.e50) beta *= alfa ;
//	printf("caso de estudio = %d peso =%lf  desp. Max=%lf  eficiencia=%lf   funcion=%lf \n",jj,valor[1],valor[2],valor[3],valor[4]);
//	printf("Materiales 1=%d    2=%d    3=%d    4=%d \n",matva[1],matva[2],matva[3],matva[4]);
deltaf = valor[4] - valom[4];
if (deltaf < 0) Actualiza() ;
else{
//	printf(" %d deltaf =%lf valor=%lf\n",jj,deltaf,valor[4]);
deltaf *= beta ;
deltaf = exp(-deltaf) ;
fact=(float)rand()/(float)RAND_MAX;
if (fact < deltaf) {
Actualiza() ;
//printf("entro \n");
//printf("deltaf =%lf \n",deltaf);
//printf("fact=%lf deltaf=%lf beta=%lf  itera =%d\n",fact,deltaf,beta,jj);
}
else Recupera();
}
}    */
}


