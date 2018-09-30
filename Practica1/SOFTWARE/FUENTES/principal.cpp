#include "arff_parser.h"
#include <vector>
#include <iostream>
#include <string>
#include <numeric>
#include <functional>
#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* Definicion del generador de numeros aleatorios */
#define MASK 2147483647
#define PRIME 65539
#define SCALE 0.4656612875e-9
#define Rand()  (( Seed = ( (Seed * PRIME) & MASK) ) * SCALE )
unsigned long Seed=2; /* semilla del generador de numeros aleatorios */

/* Parametro para el enfriamiento de la temperatura por
   descenso proporcional*/
#define TFACTR 0.9

using namespace std;



/**********************************************************/
/*            CALCULOS EUCLIDEOS                          */
/**********************************************************/
template<class Iter_T, class Iter2_T>
double vectorDistance(Iter_T first, Iter_T last, Iter2_T first2) {
  double ret = 0.0;
  while (first != last) {
    double dist = (*first++) - (*first2++);
    ret += dist * dist;
  }
  return ret > 0.0 ? sqrt(ret) : 0.0;
}

double calculaEuclidea(vector<double> &aux,vector<double> &ejemploElegido){
	return vectorDistance(aux.begin(), aux.end(), ejemploElegido.begin());
}

void distanciaEuclidea(vector<double> &ejemploElegido,vector<vector <double> > &vConMascara,vector<double> &vectorDistancias,int numeroFilas,int numeroAtributos){
	vector<double> aux;
		for(int i=0;i<numeroFilas;i++){
			aux=vConMascara[i];
			vectorDistancias.push_back(calculaEuclidea(aux,ejemploElegido));
		}
}

/**********************************************************/
/*            LECTOR ARFF                                 */
/**********************************************************/
void lectorARFF(int idArchivo,	std::vector<vector <double> > &v,int &numeroFilas,int &numeroAtributos){
	string nombreArchivo;
	switch(idArchivo){
		    case 1  :			nombreArchivo="./datos/libras00.arff";		    break;
		    case 2  :			nombreArchivo="./datos/libras01.arff";			break;
		    case 3  :			nombreArchivo="./datos/libras02.arff";			break;
		    case 4  :			nombreArchivo="./datos/libras03.arff";			break;
		    case 5  :			nombreArchivo="./datos/libras04.arff";			break;
		    case 6  :			nombreArchivo="./datos/libras05.arff";			break;
		    case 7  :			nombreArchivo="./datos/libras06.arff";			break;
		    case 8  :			nombreArchivo="./datos/libras07.arff";			break;
		    case 9  :			nombreArchivo="./datos/libras08.arff";			break;
		    case 10  :			nombreArchivo="./datos/libras09.arff";			break;

		    case 11  :			nombreArchivo="./datos/arrhythmia00.arff";		    break;
		    case 12  :			nombreArchivo="./datos/arrhythmia01.arff";			break;
		    case 13  :			nombreArchivo="./datos/arrhythmia02.arff";			break;
		    case 14  :			nombreArchivo="./datos/arrhythmia03.arff";			break;
		    case 15  :			nombreArchivo="./datos/arrhythmia04.arff";			break;
		    case 16  :			nombreArchivo="./datos/arrhythmia05.arff";			break;
		    case 17  :			nombreArchivo="./datos/arrhythmia06.arff";			break;
		    case 18  :			nombreArchivo="./datos/arrhythmia07.arff";			break;
		    case 19  :			nombreArchivo="./datos/arrhythmia08.arff";			break;
		    case 20  :			nombreArchivo="./datos/arrhythmia09.arff";			break;

		    case 21  :			nombreArchivo="./datos/wdbc00.arff";		    break;
		    case 22  :			nombreArchivo="./datos/wdbc01.arff";			break;
		    case 23  :			nombreArchivo="./datos/wdbc02.arff";			break;
		    case 24  :			nombreArchivo="./datos/wdbc03.arff";			break;
		    case 25  :			nombreArchivo="./datos/wdbc04.arff";			break;
		    case 26  :			nombreArchivo="./datos/wdbc05.arff";			break;
		    case 27  :			nombreArchivo="./datos/wdbc06.arff";			break;
		    case 28  :			nombreArchivo="./datos/wdbc07.arff";			break;
		    case 29  :			nombreArchivo="./datos/wdbc08.arff";			break;
		    case 30  :			nombreArchivo="./datos/wdbc09.arff";			break;
	}
	
	ArffParser parser(nombreArchivo);
	ArffData *data = parser.parse();
	
	numeroFilas=data->num_instances();
	numeroAtributos=data->num_attributes();

	//cout << data->num_attributes() << endl;
	//cout << data->num_instances() << endl;
	
	ArffValue *valor;
	ArffInstance *instancia;

	for(int j=0;j<numeroFilas;j++){
		instancia = data->get_instance(j);
		std::vector<double> aux;

		for (int i = 0; i < data->num_attributes(); i++) {
			aux.push_back(stod(instancia->get(i)->operator string()) );
		}
		v.push_back(aux);
		aux.clear();
	}

}

/**********************************************************/
/*            MAXIMOS Y MINIMOS POR COLUMNAS              */
/**********************************************************/
void calculaMinMax(vector<vector <double> > &v,vector<double> &minimoColumnas,vector<double>  &maximoColumnas,int numeroFilas,int numeroAtributos){
	maximoColumnas=minimoColumnas=v[0];
	for(int i=0;i<numeroFilas;i++){
		for(int j=0;j<numeroAtributos;j++){
			if(v[i][j]>maximoColumnas[j])
				maximoColumnas[j]=v[i][j];
			if(v[i][j]<minimoColumnas[j])
				minimoColumnas[j]=v[i][j];
		}
	}

}

/**********************************************************/
/*            NORMALIZANDO UNA MATRIZ                     */
/**********************************************************/
void normalizador(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado,int numeroFilas,int numeroAtributos,vector<double> &minimoColumnas,vector<double>  &maximoColumnas){
	std::vector<double> aux;
	for(int i=0;i<numeroFilas;i++){
		for(int j=0;j<numeroAtributos-1;j++){
			if(isnan((v[i][j]-minimoColumnas[j])/(maximoColumnas[j]-minimoColumnas[j])))
				aux.push_back(0);
			else
			aux.push_back((v[i][j]-minimoColumnas[j])/(maximoColumnas[j]-minimoColumnas[j]));
		}
		vNormalizado.push_back(aux);
		aux.clear();
	}
}

/**********************************************************/
/*            APLICAR MASCARA A UNA MATRIZ                */
/**********************************************************/
void aplicaMascara(vector<vector <double> > &vNormalizado,vector<vector <double> > &vConMascara,vector<int> &mask,int numeroFilas,int numeroAtributos){
	std::vector<double> aux;

	for(int i=0;i<numeroFilas;i++){
		for(int j=0;j<numeroAtributos-1;j++){
			if(mask[j]!=0)
				aux.push_back(vNormalizado[i][j]*mask[j]);
			else
				aux.push_back(0);
		}
		vConMascara.push_back(aux);
		aux.clear();
	}
	/*
		for(int i=0;i<numeroFilas;i++){
			for(int j=0;j<numeroAtributos-1;j++){
				cout << vConMascara[i][j] << " ";
			}
			cout << endl;
		}
	*/
}

/**********************************************************/
/*       DEVUELVE LOS TRES MENORES VALORES DE DISTANCIA  */
/**********************************************************/
void calculaMenores(vector<double> &vectorDistancias,int &valor1, int &valor2, int &valor3){
	valor1 = 0;
    for(int i=0; i < vectorDistancias.size(); ++i){
        if(vectorDistancias[i] < vectorDistancias[valor1])
            valor1 = i;        
    }

    if(valor1 > 0)   	valor2 = 0;
    else				valor2 = valor1+1;


    for(int i=0; i < vectorDistancias.size(); ++i){
        if(vectorDistancias[i] < vectorDistancias[valor2] && i != valor1)
            valor2 = i;        
    }

    if(valor2 <= vectorDistancias.size())	valor3 = vectorDistancias.size()-1;
    else									valor3=0;

    
    for(int i=0; i < vectorDistancias.size(); ++i){
        if(vectorDistancias[i] <= vectorDistancias[valor3] && i != valor1 && i != valor2)
            valor3 = i;        
    }
}

/**********************************************************/
/*      DEVUELVE EL INDICE CON MAYOR NUMERO DE ACIERTOS   */
/**********************************************************/
int indiceMaximo(vector<int> &resultado){
	int indice=0;
	for(int i=0;i<resultado.size();i++)
			if(resultado[i]>resultado[indice])
				indice=i;
	return indice;
}

/**********************************************************/
/*       CLASIFICADOR KNN-3				                  */
/**********************************************************/
int knnClas(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado,int ejemplo, vector<int> &mask,int numeroFilas,int numeroAtributos){
	vector<vector <double> > vConMascara;
	aplicaMascara(vNormalizado,vConMascara,mask,numeroFilas,numeroAtributos);

	vector<double> ejemploElegido=vConMascara[ejemplo];

	vector<double> vectorDistancias;
	distanciaEuclidea(ejemploElegido,vConMascara,vectorDistancias,numeroFilas,numeroAtributos-1);

	int valor1, valor2,  valor3;
	calculaMenores(vectorDistancias,valor1, valor2, valor3);
	//cout << valor1 << " "<< valor2 << " "<< valor3 << vectorDistancias[valor1] << " "<< vectorDistancias[valor2] << " "<< vectorDistancias[valor3] << endl;

	if(v[valor1][numeroAtributos-1]==v[valor2][numeroAtributos-1])
		return v[valor1][numeroAtributos-1];
	if(v[valor2][numeroAtributos-1]==v[valor3][numeroAtributos-1])
		return v[valor2][numeroAtributos-1];
	if(v[valor3][numeroAtributos-1]==v[valor1][numeroAtributos-1])
		return v[valor3][numeroAtributos-1];
	return v[valor1][numeroAtributos-1];

}

/**********************************************************/
/*            IMPLEMENTACION SFS                          */
/**********************************************************/
void sfs(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, vector<int> &mask,int numeroFilas,int numeroAtributos,vector<int> &resultado){
	bool seguir=true;
	vector<int> mascaraAuxiliar=resultado;
	int aciertosAnt=0;
	while(seguir){
			resultado.assign(resultado.size(),0);
			for(int i=0;i<numeroAtributos-1;i++){
				mask[i]=1;
				for(int j=0;j<numeroFilas;j++){
					if(knnClas(v,vNormalizado,j,mask,numeroFilas,numeroAtributos)==v[j][numeroAtributos-1])
					resultado[i]=resultado[i]+1;
				}
				if(mascaraAuxiliar[i]!=1)
				mask[i]=0;
			}
			
			if(resultado[indiceMaximo(resultado)]>aciertosAnt){
				mascaraAuxiliar[indiceMaximo(resultado)]=1;
				mask[indiceMaximo(resultado)]=1;
				aciertosAnt=resultado[indiceMaximo(resultado)];	
				//cout <<indiceMaximo(resultado) <<" " << aciertosAnt<<endl;

			}
			else
				seguir=false;

	}
	resultado=mascaraAuxiliar;
}

/**********************************************************/
/*            GENERADOR DE LLAMADAS PARA SFS              */
/**********************************************************/
void secuencial(int numSimulacion){
	cout << "LANZANDO SFS, Numero de simulacion:  "<< numSimulacion <<endl;
	std::vector<vector <double> > vectorEntero;
	std::vector<vector <double> > mitad1;
	std::vector<vector <double> > mitad2;

	int numeroFilasI;
	int numeroAtributos;
	lectorARFF(numSimulacion,vectorEntero,numeroFilasI,numeroAtributos);
	//Separando los datos en dos vectores, mitad en cada uno
	for(int i=0;i<numeroFilasI;i++){
		if(i<numeroFilasI/2)
			mitad1.push_back(vectorEntero[i]);
		else
			mitad2.push_back(vectorEntero[i]);
	}

	//Guardamos el numero de filas de cada una de las matrices
	int filasMitad1= numeroFilasI-numeroFilasI/2-1;
	int filasMitad2= numeroFilasI-filasMitad1;

	vector<double> minimoColumnas;
	vector<double> maximoColumnas;


	clock_t startTime = clock();
	std::vector<vector <double> > vNormalizado;
	std::vector<vector <double> > vNormalizado2;
	calculaMinMax(mitad1,minimoColumnas,maximoColumnas,filasMitad1,numeroAtributos);
	normalizador(mitad1,vNormalizado,filasMitad1,numeroAtributos,minimoColumnas,maximoColumnas);


	//	for(int i=0;i<numeroAtributos-1;i++)
		//cout << vNormalizado[0][i]<<" ";
		
	vector<int> mask;
	vector<int> resultado;
	for(int i=0;i<numeroAtributos-1;i++){
		mask.push_back(0);
		resultado.push_back(0);
	}

	 sfs(mitad1,vNormalizado,mask,filasMitad1,numeroAtributos,resultado);
	for(int i=0;i<numeroAtributos-1;i++)
		cout << resultado[i] << " ";
	cout << endl;
	double contadorCaracteristicas=0;
		for(int i=0;i<numeroAtributos-1;i++)
			if(resultado[i]==1)contadorCaracteristicas++;
	cout << "T_R:"<< 100*(((numeroAtributos-1)-(contadorCaracteristicas))/(numeroAtributos-1)) <<endl;


	double contadoraciertos=0;
	double contadorfallos=0;
	calculaMinMax(mitad2,minimoColumnas,maximoColumnas,mitad2.size(),numeroAtributos);
	normalizador(mitad2,vNormalizado2,mitad2.size(),numeroAtributos,minimoColumnas,maximoColumnas);
	for(int i=0;i<mitad2.size();i++)
		if(knnClas(mitad2,vNormalizado2,i,resultado,mitad2.size(),numeroAtributos)==mitad2[i][numeroAtributos-1])
			contadoraciertos++;
		else
			contadorfallos++;

	float porcentaje=	contadoraciertos*100/(contadoraciertos+contadorfallos);

	cout << "ACIERTOS: "<<contadoraciertos << "   FALLOS: "<<contadorfallos<< "   Porcentaje Acierto: "<< porcentaje <<endl;

	cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << endl;
	
}

/**********************************************************/
/*GENERADOR DE TODOS LOD VECINOS A PARTIR DE UNA MASCARA  */
/**********************************************************/
void generaVecinos(vector<int> mask,vector< vector<int> >&resultado){
	for(int i=0;i<mask.size();i++){
		if(mask[i]==0)
			mask[i]=1;
		else
			mask[i]=0;
		resultado.push_back(mask);
		
		if(mask[i]==0)
			mask[i]=1;
		else
			mask[i]=0;
	}
}

/**********************************************************/
/*   IMPLEMENTACION DE BUSQUEDA LOCAL LS                  */
/**********************************************************/
void ls(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, int numeroFilas,int numeroAtributos,vector<int> &salida){
	vector<int> mask;
	int resultado;
		for(int i=0;i<numeroAtributos-1;i++){
			if(Rand()>0.5)
				mask.push_back(1);
			else
				mask.push_back(0);
		}

	bool continua=true;
	int aciertosAnt=0;

		while(continua){
		vector<vector<int> > vecindario;
			generaVecinos(mask,vecindario);
			bool seguirfor=true;
			for(int i=0;i<vecindario.size()&&seguirfor;i++){
				resultado=0;
				for(int j=0;j<numeroFilas;j++){
					if(knnClas(v,vNormalizado,j,vecindario[i],numeroFilas,numeroAtributos)==v[j][numeroAtributos-1])
						resultado++;
				}
				if(aciertosAnt<resultado){
					aciertosAnt=resultado;
					seguirfor=false;
					mask=vecindario[i];
				}				
			}
			if(seguirfor)
				continua=false;
		}

		salida=mask;
}

/**********************************************************/
/*            GENERADOR DE LLAMADAS PARA LS               */
/**********************************************************/
void primeroMejor(int numSimulacion){
	std::vector<vector <double> > vectorEntero;
	std::vector<vector <double> > mitad1;
	std::vector<vector <double> > mitad2;

	int numeroFilasI;
	int numeroAtributos;
	lectorARFF(numSimulacion,vectorEntero,numeroFilasI,numeroAtributos);

	//Separando los datos en dos vectores, mitad en cada uno
	for(int i=0;i<numeroFilasI;i++){
		if(i<numeroFilasI/2)
			mitad1.push_back(vectorEntero[i]);
		else
			mitad2.push_back(vectorEntero[i]);
	}
	//Guardamos el numero de filas de cada una de las matrices
	int filasMitad1= numeroFilasI-numeroFilasI/2-1;
	int filasMitad2= numeroFilasI-filasMitad1;
	//Vectores que almacenan los maximos y minimos de las columnas para normalizar
	vector<double> minimoColumnas;
	vector<double> maximoColumnas;

	//Reloj de sistema
	clock_t startTime = clock();
	std::vector<vector <double> > vNormalizado;
	std::vector<vector <double> > vNormalizado2;
	vector<int> salida;

	//Normalizando la matriz 1
	calculaMinMax(mitad1,minimoColumnas,maximoColumnas,filasMitad1,numeroAtributos);
	normalizador(mitad1,vNormalizado,filasMitad1,numeroAtributos,minimoColumnas,maximoColumnas);

	//Llamando localSeach
	ls(mitad1,vNormalizado,filasMitad1,numeroAtributos,salida);
	for(int i=0;i<numeroAtributos-1;i++)
		cout << salida[i] << " ";
	cout << endl;
	double contadorCaracteristicas=0;
		for(int i=0;i<numeroAtributos-1;i++)
			if(salida[i]==1)contadorCaracteristicas++;
	cout << "T_R:"<< 100*(((numeroAtributos-1)-(contadorCaracteristicas))/(numeroAtributos-1)) <<endl;


	double contadoraciertos=0;
	double contadorfallos=0;
	//Normalizando la matriz 2
	calculaMinMax(mitad2,minimoColumnas,maximoColumnas,mitad2.size(),numeroAtributos);
	normalizador(mitad2,vNormalizado2,mitad2.size(),numeroAtributos,minimoColumnas,maximoColumnas);
	
	for(int i=0;i<mitad2.size();i++)
		if(knnClas(mitad2,vNormalizado2,i,salida,mitad2.size(),numeroAtributos)==mitad2[i][numeroAtributos-1])
			contadoraciertos++;
		else
			contadorfallos++;

	float porcentaje=	contadoraciertos*100/(contadoraciertos+contadorfallos);


	cout << "ACIERTOS: "<<contadoraciertos << "   FALLOS: "<<contadorfallos<< "   Porcentaje Acierto: "<< porcentaje <<endl;
	cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << endl;

}

/**********************************************************/
/*            CRITERIO DE METROPOLIS                      */
/**********************************************************/
int metrop (float de, float t){
 return de < 0.0 || Rand() < exp(-de/t);
 }

/**********************************************************/
/*            DADA UNA MASCARA LE HACE UN SWAP ALEATORIO  */
/**********************************************************/
void generaVecinoAleatorio(vector<int> &mask){
	int posicion= ((int)Rand())%mask.size();
	mask[posicion]=!mask[posicion];
}

/**********************************************************/
/*            IMPLEMENTACION ANNEAL (VERSION VIEJA)       */
/**********************************************************/
void annealV(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, int numeroFilas,int numeroAtributos,vector<int> &salida){
	//Calculo de una primera mascara aleatoria
	vector<int> mask;
	int resultado;
	int continua=0;
	int aciertosAnt=0;
		for(int i=0;i<numeroAtributos-1;i++){

			if(Rand()>0.5)
				mask.push_back(1);
			else
				mask.push_back(0);
		}
		
		//Calculo del coste usando KNN y lo almacenamos
		resultado=0;
		for(int j=0;j<numeroFilas;j++){
			if(knnClas(v,vNormalizado,j,mask,numeroFilas,numeroAtributos)==v[j][numeroAtributos-1])
				resultado++;
		}
		aciertosAnt=resultado;


	float tInicial=(0.3*aciertosAnt)/-log(0.3);
	float tActual=tInicial;
	float tFinal=0.001;
	cout << "temperatura inicial: "<<tInicial;

	int maxVecinos=numeroFilas;
	int maxExitos=0.1*maxVecinos;
	int evaluaciones=500;
	int M=evaluaciones/(maxVecinos*maxVecinos);
	float beta=(tInicial-tFinal)/M*tInicial*tFinal;

	for(int i=0;i<2;i++){
		int numeroVecinos=0;
		int exitos=0;
		bool enfriar=false;
		while(!enfriar){
			//Calculamos vecino aleatorio
			auto aux=mask;
			generaVecinoAleatorio(aux);
			numeroVecinos++;

			//Calculamos costo de la nueva solucion
			resultado=0;
			for(int j=0;j<numeroFilas;j++){
				if(knnClas(v,vNormalizado,j,aux,numeroFilas,numeroAtributos)==v[j][numeroAtributos-1])
					resultado++;
			}
			if(aciertosAnt<resultado){
				aciertosAnt=resultado;
				mask=aux;
				exitos++;
			}
			else{
				if(metrop(resultado-aciertosAnt,tActual)==1){
					aciertosAnt=resultado;
					mask=aux;
					exitos++;
				}				
			}

			printf("\nIteracion: %d.\n",i);
	     	
	     	//Si se han generado los vecinos o los exitos necesarios, se enfria
	     	if(numeroVecinos==maxVecinos||exitos==maxExitos)
	     		enfriar=true;
	    }
	    //Enfriamos temperatura
	    tActual =tActual/(1+beta*tActual);	

	}

	salida=mask;
}


/**********************************************************/
/*           IMPLEMENTACION anneal 		                  */
/**********************************************************/
void anneal(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, int numeroFilas,int numeroAtributos,vector<int> &salida){
	//Calculo de una primera mascara aleatoria
	vector<int> mask;
	int resultado;
	int continua=0;
	int aciertosAnt=0;
		for(int i=0;i<numeroAtributos-1;i++){

			if(Rand()>0.5)
				mask.push_back(1);
			else
				mask.push_back(0);

		}
		
		//Calculo del coste usando KNN y lo almacenamos
		resultado=0;
		for(int j=0;j<numeroFilas;j++){
			if(knnClas(v,vNormalizado,j,mask,numeroFilas,numeroAtributos)==v[j][numeroAtributos-1])
				resultado++;
		}
		aciertosAnt=resultado;


	float tInicial=(0.3*aciertosAnt)/-log(0.3);
	float tActual=tInicial;
	float tFinal=0.001;
	//cout << "temperatura inicial: "<<tInicial;

	int maxVecinos=10*numeroFilas;
	int maxExitos=0.1*maxVecinos;
	int evaluaciones=15000;
	int M=evaluaciones/(maxVecinos*maxVecinos);
	float beta=(tInicial-tFinal)/M*tInicial*tFinal;
	int exitos;
	do{
		int numeroVecinos=0;
		exitos=0;
		for(int i=0;i<maxVecinos;i++){
			//Generamos un vecino y lo almacenamos
			auto aux=mask;
			generaVecinoAleatorio(aux);
			numeroVecinos++;

			//Calculamos costo de la nueva solucion
				resultado=0;
				for(int j=0;j<numeroFilas;j++){
					if(knnClas(v,vNormalizado,j,aux,numeroFilas,numeroAtributos)==v[j][numeroAtributos-1])
						resultado++;
				}
				if(aciertosAnt<resultado){
					aciertosAnt=resultado;
					mask=aux;
					exitos++;
				}
				else{
					if(metrop(resultado-aciertosAnt,tActual)==1){
						aciertosAnt=resultado;
						mask=aux;
						exitos++;
					}				
				}
			//Si se han generado los vecinos o los exitos necesarios, se sale del bucle
		     if(exitos==maxExitos){
		     	i=evaluaciones;
		     }
		}

	    //Enfriamos temperatura
	    tActual =tActual/(1+beta*tActual);	
	    //cout <<"Temperadura: final:"<< tActual<<endl;

	}while(tActual>tFinal&&exitos>0);
	
	salida=mask;
}

/**********************************************************/
/*            GENERADOR LLAMADAS ANNEAL                   */
/**********************************************************/

void sa(int numSimulacion){
	std::vector<vector <double> > vectorEntero;
	std::vector<vector <double> > mitad1;
	std::vector<vector <double> > mitad2;

	int numeroFilasI;
	int numeroAtributos;
	lectorARFF(numSimulacion,vectorEntero,numeroFilasI,numeroAtributos);

//Separando los datos en dos vectores, mitad en cada uno
	for(int i=0;i<numeroFilasI;i++){
		if(i<numeroFilasI/2)
			mitad1.push_back(vectorEntero[i]);
		else
			mitad2.push_back(vectorEntero[i]);
	}
//Guardamos el numero de filas de cada una de las matrices
	int filasMitad1= numeroFilasI-numeroFilasI/2-1;
	int filasMitad2= numeroFilasI-filasMitad1;
//Vectores que almacenan los maximos y minimos de las columnas para normalizar
	vector<double> minimoColumnas;
	vector<double> maximoColumnas;

//Reloj de sistema
	clock_t startTime = clock();
	std::vector<vector <double> > vNormalizado;
	std::vector<vector <double> > vNormalizado2;
	vector<int> salida;

	//Normalizando la matriz 1
	calculaMinMax(mitad1,minimoColumnas,maximoColumnas,filasMitad1,numeroAtributos);
	normalizador(mitad1,vNormalizado,filasMitad1,numeroAtributos,minimoColumnas,maximoColumnas);

	//Llamando localSeach
	anneal(mitad1,vNormalizado,filasMitad1,numeroAtributos,salida);
	for(int i=0;i<numeroAtributos-1;i++)
		cout << salida[i] << " ";
	cout << endl;
	double contadorCaracteristicas=0;
		for(int i=0;i<numeroAtributos-1;i++)
			if(salida[i]==1)contadorCaracteristicas++;
	cout << "T_R:"<< 100*(((numeroAtributos-1)-(contadorCaracteristicas))/(numeroAtributos-1)) <<endl;


	double contadoraciertos=0;
	double contadorfallos=0;
	//Normalizando la matriz 2
	calculaMinMax(mitad2,minimoColumnas,maximoColumnas,mitad2.size(),numeroAtributos);
	normalizador(mitad2,vNormalizado2,mitad2.size(),numeroAtributos,minimoColumnas,maximoColumnas);
	
	for(int i=0;i<mitad2.size();i++)
		if(knnClas(mitad2,vNormalizado2,i,salida,mitad2.size(),numeroAtributos)==mitad2[i][numeroAtributos-1])
			contadoraciertos++;
		else
			contadorfallos++;

	float porcentaje=	contadoraciertos*100/(contadoraciertos+contadorfallos);

	cout << "ACIERTOS: "<<contadoraciertos << "   FALLOS: "<<contadorfallos<< "   Porcentaje Acierto: "<< porcentaje <<endl;

	cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << endl;
}



/**********************************************************/
/*           GENERADOR DE ENTORNOS PARA TABU              */
/**********************************************************/
void generaEntorno(vector<int> &mask,vector< vector<int> > &resultado,int numero){
	resultado.clear();
	for(int i=0;i<numero;i++){
		auto aux=mask;
		generaVecinoAleatorio(aux);
		resultado.push_back(aux);
	}
}

int indiceDistinto(vector<int> &mask,vector<int> &mask2){
	for(int i=0;i<mask.size();i++)
		if(mask[i]!=mask2[i])
			return i;
}


/**********************************************************/
/*           IMPLEMENTACION tabu 		                  */
/**********************************************************/
void tabu(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, int numeroFilas,int numeroAtributos,vector<int> &salida){
	//Calculo de una primera mascara aleatoria
	vector<int> mask;
	int resultado;
	int continua=0;
	int aciertosAnt=0;
	int mejorAcierto=0;
	int tamaTabu=numeroAtributos/30;
	vector<int> mascaraTabu(numeroAtributos-1,0);

	for(int i=0;i<numeroAtributos-1;i++){
			if(Rand()>0.5)
				mask.push_back(1);
			else
				mask.push_back(0);
	}
		
	//Calculo del coste usando KNN y lo almacenamos
	resultado=0;
	for(int j=0;j<numeroFilas;j++){
		if(knnClas(v,vNormalizado,j,mask,numeroFilas,numeroAtributos)==v[j][numeroAtributos-1])
			resultado++;
	}
	aciertosAnt=mejorAcierto=resultado;

	vector< vector<int> > entorno;

	for(int l=0;l<50;l++){
		//Generamosel entorno con la mascara actual
		generaEntorno(mask,entorno,30);
		int indiceMejorVecino;
		int resultadoMejorVecino=0;

		//Calculamos el mejor vecino posible
		for(int i=0;i<entorno.size();i++){
			int resultado=0;
			for(int j=0;j<numeroFilas;j++){
				if(knnClas(v,vNormalizado,j,entorno[i],numeroFilas,numeroAtributos)==v[j][numeroAtributos-1])
					resultado++;
				}
			if(resultado>resultadoMejorVecino){
				resultadoMejorVecino=resultado;
				indiceMejorVecino=i;
			}
		}

		//Miramos el indice que cambia con respecto a la mascara que tenemos actualmente
		int indice= indiceDistinto(mask,entorno[indiceMejorVecino]);
		//Si el indice no está marcadocomo Tabu, se marca como tabú, se actualiza la mascara y se disminuye el vlaor tabú de todos los elementos >0
		//Actualizamos la mascara actual ademas de actualizar los mejores aciertos que tenemos hasta ahora y anteriormente
		if(mascaraTabu[indice]==0){
			for(int k=0;k<mascaraTabu.size();k++){
				if(mascaraTabu[k]>0) mascaraTabu[k]--;
			}
			mascaraTabu[indice]=tamaTabu;
			mask=entorno[indiceMejorVecino];
			aciertosAnt=resultado;
			if(resultado>mejorAcierto)
				mejorAcierto=resultado;
		}
		else{ //Si existia el indice como tabú, miramos el criterio de aspiracion, si es mejor que el global, lo guardamos como mejor
			//Actualizamos lo necesario, la mascara y cambiamos los valores tabú de los elementos.
			if(resultado>mejorAcierto){
				mejorAcierto=resultado;
				for(int k=0;k<mascaraTabu.size();k++){
					if(mascaraTabu[k]>0) mascaraTabu[k]--;
				}
				mascaraTabu[indice]+=tamaTabu;
				mask=entorno[indiceMejorVecino];
				aciertosAnt=resultado;
			}
		}

	}

	salida=mask;
}



/**********************************************************/
/*            GENERADOR LLAMADAS TABU                   */
/**********************************************************/

void ts(int numSimulacion){
	std::vector<vector <double> > vectorEntero;
	std::vector<vector <double> > mitad1;
	std::vector<vector <double> > mitad2;

	int numeroFilasI;
	int numeroAtributos;
	lectorARFF(numSimulacion,vectorEntero,numeroFilasI,numeroAtributos);

//Separando los datos en dos vectores, mitad en cada uno
	for(int i=0;i<numeroFilasI;i++){
		if(i<numeroFilasI/2)
			mitad1.push_back(vectorEntero[i]);
		else
			mitad2.push_back(vectorEntero[i]);
	}
//Guardamos el numero de filas de cada una de las matrices
	int filasMitad1= numeroFilasI-numeroFilasI/2-1;
	int filasMitad2= numeroFilasI-filasMitad1;
//Vectores que almacenan los maximos y minimos de las columnas para normalizar
	vector<double> minimoColumnas;
	vector<double> maximoColumnas;

//Reloj de sistema
	std::vector<vector <double> > vNormalizado;
	std::vector<vector <double> > vNormalizado2;
	vector<int> salida;

	//Normalizando la matriz 1
	calculaMinMax(mitad1,minimoColumnas,maximoColumnas,filasMitad1,numeroAtributos);
	normalizador(mitad1,vNormalizado,filasMitad1,numeroAtributos,minimoColumnas,maximoColumnas);

	clock_t startTime = clock();
	//Llamando localSeach
	tabu(mitad1,vNormalizado,filasMitad1,numeroAtributos,salida);
	for(int i=0;i<numeroAtributos-1;i++)
		cout << salida[i] << " ";
	cout << endl;
	double contadorCaracteristicas=0;
		for(int i=0;i<numeroAtributos-1;i++)
			if(salida[i]==1)contadorCaracteristicas++;
	cout << "T_R:"<< 100*(((numeroAtributos-1)-(contadorCaracteristicas))/(numeroAtributos-1)) <<endl;


	double contadoraciertos=0;
	double contadorfallos=0;
	//Normalizando la matriz 2
	calculaMinMax(mitad2,minimoColumnas,maximoColumnas,mitad2.size(),numeroAtributos);
	normalizador(mitad2,vNormalizado2,mitad2.size(),numeroAtributos,minimoColumnas,maximoColumnas);
	
	for(int i=0;i<mitad2.size();i++)
		if(knnClas(mitad2,vNormalizado2,i,salida,mitad2.size(),numeroAtributos)==mitad2[i][numeroAtributos-1])
			contadoraciertos++;
		else
			contadorfallos++;

	float porcentaje=	contadoraciertos*100/(contadoraciertos+contadorfallos);

	cout << "ACIERTOS: "<<contadoraciertos << "   FALLOS: "<<contadorfallos<< "   Porcentaje Acierto: "<< porcentaje <<endl;

	cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds." << endl;
}


/*************************************************************************************/
/*           MAIN, PARA MENU Y LLAMADAS A LAS DIFERENTES FUNCIONES Y SIMULACIONES    */
/*************************************************************************************/

int main(){

	int opcionMenu;
	int numSimulacion;	
	do{
cout << "Menu de acciones:"<<endl;
cout << "1.- SFS:"<<endl;
cout << "2.- LS:"<<endl;
cout << "3.- SA:"<<endl;
cout << "4.- BT:"<<endl;
cout << "5.- SFS Completa:"<<endl;
cout << "6.- LS Completa:"<<endl;
cout << "7.- SA Completa:"<<endl;
cout << "8.- BT Completa:"<<endl;
cout << "9.- AYUDA"<<endl;
cout << "0.- SALIR:"<<endl;
cout << "---------:";
	cin >> opcionMenu;
	switch(opcionMenu){
	    case 1  :
			cout << "Numero de simulacion: 1-30 ";
			cin >> numSimulacion;
			secuencial(numSimulacion);
	    	break;
	    case 2  :
			cout << "Numero de simulacion: 1-30 ";
			cin >> numSimulacion; 
			primeroMejor(numSimulacion);
	    	break; 
	    case 3  :
			cout << "Numero de simulacion: 1-30 ";
			cin >> numSimulacion; 
			sa(numSimulacion);
	    	break; 	
	    case 4  :
			cout << "Numero de simulacion: 1-30 ";
			cin >> numSimulacion; 
			ts(numSimulacion);
	    	break; 	
	    case 5  :
			cout << "Simulando";
			for(int i=1;i<31;i++)
			secuencial(i);
	    	break;
	    case 6  :
			cout << "Simulando";
			for(int i=1;i<31;i++)
			primeroMejor(i);
	    	break; 
	    case 7  :
			cout << "Simulando";
			for(int i=1;i<31;i++)
			sa(i);
	    	break; 	
	    case 8  :
			cout << "Simulando";
			for(int i=1;i<31;i++)
			ts(i);
	    	break;
	    case 9  :
			cout << "Las simulaciones van de la 1 a la 30, divididas en 3 partes:"<<endl;
			cout << "1-10- Libras:"<<endl;
			cout << "11-20- Arrhythmia:"<<endl;
			cout << "21-30- Wdbc0:"<<endl;
			cout << "Las opciones 5-8 hacen todas las simulaciones en el orden anteriormente especificado"<<endl;
	    	break; 	
	    case 0  :
			cout << "SALIENDO..."<<endl;
	    	break; 

	    default : //Optional
	       cout << "OPCION INCORRECTA"<<endl;
	}
	}while(opcionMenu!=0);

}