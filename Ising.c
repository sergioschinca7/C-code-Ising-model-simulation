#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Declaro PARAMETROS CONSTANTES */
#define L 20
#define SEED 378520
#define J 1.0

// rango de temperaturas
#define T_max 5.0
#define T_min 0.01 
#define T_step 0.01

// sampleo
#define pasos_term (4*L*L) // pasos de metropolis para termalizar
#define n_sample 2000 // cantidad de samples que tomo para cada T 
#define dist_sample (2*L*L)// pasos de metropolis entre dos samples para descorrelacionar

#define SIMU 10 // cantidad de corridas simuladas desde configuraciones iniciales aleatorias distintas 

/* Declaro VARIABLES y ARRAYS */
static int red[L][L]; // red de spines
static int M; // magnetizacion
static int E; // energia

static double acept; // proporcion de aceptacion de pasos de MM

/* Declaro SUBRUTINAS */
void llenar(double p,int semilla);
void imprimir();
void metropolis(double T);
void flip(int fila, int col,double T);

int calcula_M_inicial();
int calcula_E_inicial();

// Output
FILE* file;

int main(){

	//Abro achivo de salida y les escribo el header
	char nom[100];
  	sprintf(nom,"Ising_statics_L%i.txt",L);   
	file=fopen(nom,"w");
	fprintf(file, "# ISING : SEED %i, J %f, pasos_term %i, n_sample %i, dist_sample %i \n", SEED, J, pasos_term, n_sample,dist_sample);
	fprintf(file,"L,T,simu,mean_M,mean_E,mean_M2, mean_E2,acept\n"); 
	
	int simu;
	for(simu=0; simu<SIMU; simu++){

		// Defino p = prob de spin up (+1) y pueblo la red
		float p = 0.5; 
		llenar(p, SEED*simu); 
		
		E = calcula_E_inicial();
		M =	calcula_M_inicial();

		// Considero distintas temperaturas T
		int cant_T = (T_max - T_min) / T_step + 1;
		int i_T;
		for (i_T=0; i_T<cant_T; i_T++){

			// Defino T
			double T = T_max - i_T * T_step; 

			// Termalizo el sistema (¡Aquí la "chanchada"!)
			int paso=0;
			for (paso=0; paso<pasos_term; paso++){
				metropolis(T);
			}

			//Sampleo
			double mean_M=0; // magnetizacion media
			double mean_E=0; // energia media
			double mean_M2=0; // magnetizacion cuadrada media
			double mean_E2=0; // energia cuadrada media
			acept=0; // prop de pasos de MM aceptados

			int sample;
			for(sample=0; sample < n_sample; sample++){
				
				//descorrelaciono las muestras
				for (paso=0; paso < dist_sample; paso++){
				metropolis(T);
			    }

			    mean_M = mean_M + M;
			    mean_E = mean_E + E;
			    mean_M2 = mean_M2 + (M*M);
			    mean_E2 = mean_E2 + (E*E);

			    printf("simu: %i, T: %f, sample:%i\n",simu, T, sample); 
			}

			mean_M = mean_M / (double) n_sample;
			mean_E = mean_E / (double) n_sample;
			mean_M2 = mean_M2 / (double) n_sample;
			mean_E2 = mean_E2 / (double) n_sample;
			acept = acept / (double)(n_sample*dist_sample);

			fprintf(file,"%i,%f,%i,%f,%f,%f,%f,%f\n",L,T,simu,mean_M,mean_E,mean_M2, mean_E2,acept); 
		}
	}

	fclose(file);
}

int calcula_E_inicial(){
	// E =  -J sum (si*sj) 
	float E =0;

	// Recorro sumando interaccion con vecino derecha y vecino abajo, para no contar doble
	int i,j, columna_derecha , fila_abajo;
	for(i=0; i<L; i++){
		for (j=0; j<L; j++){
			fila_abajo = (i+1)%L;
			columna_derecha = (j+1)%L;

			E  = E + red[i][j] * ( red[i][columna_derecha] + red[fila_abajo][j] );	
		}
	}
	E = -J * E;
return E;
}

int calcula_M_inicial(){
	// M = sum(s_i) 
	float res =0;
	int i,j;
	for(i=0; i<L; i++) {
		for (j=0; j<L; j++){
			res = res + red[i][j];
		}
	}
return res;
}

void metropolis(double T){
	// elijo un sitio al azar (i,j)=(fila,col)
	int fila =rand()%L;
	int col =rand()%L;
	flip(fila,col,T);
}


void flip(int fila, int col, double T){
	
	// condiciones periodicas de borde
	int fila_abajo = (fila+1)%L;
	int fila_arriba = (fila-1+L)%L;
	int columna_derecha = (col+1)%L;
	int columna_izquierda = (col-1+L)%L;
	

	int suma_vecinos = red[fila_abajo][col] + red[fila_arriba][col] + red[fila][columna_izquierda] + red[fila][columna_derecha];
	int dE= 2 * J * red [fila][col] * suma_vecinos; // dE= 2*J*s_i*(s1+s2+s3+s4) 
	
	//si dE < 0 flipeo, sino chequeo con la proba de aceptacion: exp(-dE/kT) con k=1 constante de Boltzmann
	if (dE<0){ 
		red[fila][col]= - red[fila][col];

		//actualizo valores
		M = M + 2 * red[fila][col];
		E = E + dE ;
		acept ++;

	}else{
		double exponente = (double)dE / T;
		double prob = exp(-exponente);
		double ratio = (double)rand()/(double)RAND_MAX;
		if (ratio <= prob){
			red[fila][col]= - red[fila][col];

			//actualizo valores
			M = M + 2 * red[fila][col];
			E = E + dE ;
			acept ++;
		} 
	}
} 

void llenar(double p,int semilla){
	int i,j;
	double ratio;
	srand(semilla);

	for(i=0;i<L;i++){
		for(j=0;j<L;j++){
			ratio = (double)rand()/(double)RAND_MAX;
			if (ratio<p){ 
				red[i][j] = 1;
			}else{
				red[i][j] = -1;
			}
		}
	}
}

void imprimir(){
	int i,j;
	printf("\n\n");

	for (i=0; i<L; i++){
		for(j=0; j<L; j++){
			printf("%i ", red[i][j]);
		}
		printf("\n");
	}
	printf("\n\n");
}
