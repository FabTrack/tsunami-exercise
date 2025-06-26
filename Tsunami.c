#include <stdio.h>
#include <math.h>

/*			 Tsunami 							
	 Fabiano Tracanna 10/09/20	*/

// Scala spaziale: sigma_x,	L = 3D  
// Scala temporale: 1 sec
int main(){
	double pi = 3.14159, g = 9.81;
	int Ttot = 3300;				// Durata evoluzione
	int N = 270, T = 330;
	int cut = 150;
	int i, j, k, l;
	double dfx, dfy;
	FILE *f, *ff;
	f = fopen("tsunami_data.txt","w+");
	ff = fopen("fondale.txt", "w+");
	
	//Altezza tsunami, epicentro, larghezza, distanza dalla costa a cui si innalza il fondo
	double r = 80000;
	double sigma_x = 80000/r, sigma_y = 200000/r;		// Larghezza onda
	double D = 286450/r, theta = pi/180.0;			// Rialzamento fondale
	double H_0 = 5000/r, L = 3*D, d = L-D;		// Fondale e griglia			
	double A = 1/r, x_c = d, y_c = L/2; 				// Altezza e coordinate onda
	g = g/r;
	printf("Paramentri adimensionali\n");
	printf("sigma_x=%f\tsigma_y=%f\nd=%f\tA=%f\nH_0=%f\tL=%f\nx_c=%f\ty_c=%f\n\n",sigma_x,sigma_y,d,A,H_0,L,x_c,y_c);
		
	double x[N+1], x2[2*N+2], h = (double)L/N, dt = (double)Ttot/T;
	double beta = dt*dt/(h*h);
	double u[N+1][N+1], v[N+1][N+1], w[N+1][N+1];		// w è allo step k-1, v è allo step k, u è allo step k+1
	
	double H[2*N+3], lambda[2*N+3];
	
	// Grid
	for(j = 0; j<=N; j++){
		x[j] = j*h;
		printf("%f\t%d\n", x[j], j);
	}
	for(j = 0; j<=2*N+1; j++){
		x2[j] = (j-1)*h/2;					// va da -h/2 ad L
		//printf("x=%f\t%d\n", x2[j], j);
	}	
	//printf("L = %f\n", L);
		
	// H e lambda
	// z = 0 sul fondo più profondo
	for(i = 0; i <= 2*N+1; i++){
		if(x2[i] <= d){
			H[i] = H_0;
		//	printf("H: %f\t%f\t%f\t%d\n",H[i],x[i],x[i],i);
		}else{
			H[i] = H_0 - (x2[i] - d)*tan(theta); 		
		//	printf("H: %f\t%f\t%f\t%d\n",H[i],x[i]-h/2,x[i],i);
		}
		lambda[i] = g*H[i];
		fprintf(ff,"%f\t%f\n", x2[i], H[i]);
		//printf("lambda: %f\t%d\n",lambda[i],i);
	}
	H[2*N+2] = 0; lambda[2*N+2] = 0;
	fprintf(ff,"%f\t%f", (x[N]+h/2), H[2*N+2]);

	//Control
	double dt_max = 1/sqrt(lambda[0])*1/sqrt(2/(h*h));
	printf("h=%f\tdt=%f\t dt_max=%f\tbeta=%f\n", h, dt, dt_max, beta);
	
	
	// Condizioni Iniziali
 	for(i = 0; i<=N; i++){
		for(j = 0; j<=N; j++){
			v[i][j] = A*exp( -pow( (x[i] - x_c )/sigma_x, 2.0) - pow( (x[j] - y_c )/sigma_y, 2.0) );
			//printf("%.20f\n", v[i][j]*r);
			//fprintf(f, "%f\t", v[i][j]*r);
		}
		//fprintf(f, "\n");
	}	
	
		
	// Evoluzione temporale		
	for(k = 1; k<=T; k++){
		// Primo step
		if(k == 1){
			for(i = 0; i<=N; i++){
				l = 2*i+1;
				//printf("lambda(i+1/2): %f\tlambda(i): %f  lambda(i-1/2): %f\tl: %d  i:%d\n", lambda[l+1], lambda[l], lambda[l-1], l, i);
				for(j = 0; j<=N; j++){
					//dfx			
					if(i == 0){
						dfx = lambda[l+1]*(v[i+1][j] - v[i][j]) - lambda[l-1]*(v[i][j]-v[i+1][j]);
					}else if(i == N){
						dfx = lambda[l+1]*(v[i-1][j] - v[i][j]) - lambda[l-1]*(v[i][j]-v[i-1][j]);
					}else{
						dfx = lambda[l+1]*(v[i+1][j] - v[i][j]) - lambda[l-1]*(v[i][j]-v[i-1][j]);
					}
					//dfy	
					if(j == 0){							
						dfy = lambda[l]*(v[i][j+1] - v[i][j]) - lambda[l]*(v[i][j]-v[i][j+1]);					
					} else if(j == N){
						dfy = lambda[l]*(v[i][j-1] - v[i][j]) - lambda[l]*(v[i][j]-v[i][j-1]);	
					}else{
						dfy = lambda[l]*(v[i][j+1] - v[i][j]) - lambda[l]*(v[i][j]-v[i][j-1]);
					}
					
					//step
					u[i][j] = v[i][j] + 0.5*beta*(dfx + dfy);	
				}
				
				//printf("%f\t%f\n",u[i][10],w[i][10]);
			}
			
		}else{
			for(i = 0; i<=N; i++){
				l = 2*i+1;
				for(j = 0; j<=N; j++){
					//dfx
					if(i == 0){
						dfx = lambda[l+1]*(v[i+1][j] - v[i][j]) - lambda[l-1]*(v[i][j]-v[i+1][j]);
					}else if(i == N){
						dfx = lambda[l+1]*(v[i-1][j] - v[i][j]) - lambda[l-1]*(v[i][j]-v[i-1][j]);
					}else{
						dfx = lambda[l+1]*(v[i+1][j] - v[i][j]) - lambda[l-1]*(v[i][j]-v[i-1][j]);
					}
					//dfy
					if(j == 0){							
						dfy = lambda[l]*(v[i][j+1] - v[i][j]) - lambda[l]*(v[i][j]-v[i][j+1]);					
					} else if(j == N){
						dfy = lambda[l]*(v[i][j-1] - v[i][j]) - lambda[l]*(v[i][j]-v[i][j-1]);	
					}else{
						dfy = lambda[l]*(v[i][j+1] - v[i][j]) - lambda[l]*(v[i][j]-v[i][j-1]);
					}

					//step
					u[i][j] = -w[i][j] + 2*v[i][j] + beta*(dfx + dfy);	
				}
				
				//printf("%f\t%f\n",u[i][10],w[i][10]);
			}

		}
		
			// Stampo
			if(k == 1){
				for(i = cut; i<=N; i++){
					//printf("i: %d\n",i);
					for(j = cut/2; j <= ( N - cut/2 + 1); j++){
						//printf("\tj: %d\n", j);
						fprintf(f, "%f\t", v[i][j]*r);
					}
					fprintf(f,"\n");
				}
					//printf("%d\n", (N/2 -cut/2));	
			}
			for(i = cut; i<=N; i++){
				for(j = cut/2; j <= ( N - cut/2 + 1); j++){
					fprintf(f, "%f\t", u[i][j]*r);
				}
				fprintf(f,"\n");
			}

		for(i = 0; i<=N; i++){
			for(j = 0; j<=N; j++){
				w[i][j] = v[i][j];			// k diventa k-1
				v[i][j] = u[i][j];			// k+1 diventa k
			}
		}
		
	}

	
	return 0;
}

