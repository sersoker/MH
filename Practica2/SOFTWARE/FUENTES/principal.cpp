#include "arff_parser.h"
#include <vector>
#include <iostream>
#include <string>
#include <numeric>
#include <functional>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include "ordenaVector.h"

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
/*            DADA UNA MASCARA LE HACE UN SWAP ALEATORIO  */
/**********************************************************/
void generaVecinoAleatorio(vector<int> &mask){
	int posicion= (int)(Rand()*MASK)%mask.size();
	mask[posicion]=!mask[posicion];
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
/*   IMPLEMENTACION DE BUSQUEDA LOCAL multiarranque                  */
/**********************************************************/
void ma(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, int numeroFilas,int numeroAtributos,vector<int> &salida){
	vector<int> mask;
	int resultado;
	vector<int> mejorMascara;
	int mejorCosteActual=0;

	//Iniciamos la mascara aleatoria
		for(int i=0;i<numeroAtributos-1;i++){
			if(Rand()>0.5)
				mask.push_back(1);
			else
				mask.push_back(0);
			mejorMascara.push_back(0);
		}


	for(int k=0;k<25;k++){
		bool continua=true;
		int aciertosAnt=0;

		while(continua){
			vector<vector<int> > vecindario;
			generaVecinos(mask,vecindario);
			bool seguirfor=true;
			int contadorEvaluaciones=0;
			
			for(int i=0;i<vecindario.size()&&seguirfor;i++){
				resultado=0;
				for(int j=0;j<numeroFilas;j++){
					if(knnClas(v,vNormalizado,j,vecindario[i],numeroFilas,numeroAtributos)==v[j][numeroAtributos-1]){
						resultado++;
						contadorEvaluaciones++;
					}
				}
				if(aciertosAnt<resultado){
					aciertosAnt=resultado;
					seguirfor=false;
					mask=vecindario[i];
				}
				else if(contadorEvaluaciones>600)	//Para que no tarde tanto, cambiamos el maximo de evaluaciones de knn
					i=vecindario.size();			
			}
			if(seguirfor)
				continua=false;
		}
		if(resultado>mejorCosteActual){
			mejorCosteActual=resultado;
			mejorMascara=mask;
		}
	}

	salida=mejorMascara;
}

/**********************************************************/
/*            GENERADOR DE LLAMADAS PARA MA               */
/**********************************************************/
void multiarranque(int numSimulacion){
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
	ma(mitad1,vNormalizado,filasMitad1,numeroAtributos,salida);
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
/*   Operador de mutacion con t caracteristicas           */
/**********************************************************/
void mutador(int numeroMutaciones,vector<int> &mask){
	for(int i=0;i<numeroMutaciones;i++){
		int posicion= (int)(Rand()*MASK)%mask.size();
		mask[posicion]=!mask[posicion];
	}

}

/**********************************************************/
/*   IMPLEMENTACION DE BUSQUEDA LOCAL ILS                  */
/**********************************************************/
void iteradals(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, int numeroFilas,int numeroAtributos,vector<int> &salida){
	vector<int> mascaraInicial;
	int resultadoInicial=0;
	int resultado;
	
	//Generamos una mascara inicial
	for(int i=0;i<numeroAtributos-1;i++){
		if(Rand()>0.5)			mascaraInicial.push_back(1);
		else			mascaraInicial.push_back(0);
	}

	//Calculamos su valor
	for(int j=0;j<numeroFilas;j++){
		if(knnClas(v,vNormalizado,j,mascaraInicial,numeroFilas,numeroAtributos)==v[j][numeroAtributos-1])
			resultadoInicial++;
	}

	//Generamos una mascara como copia de la primera y la pasamos por LS
	vector<int>mascaraSecundaria(mascaraInicial);
	for(int k=0;k<25;k++){
		bool continua=true;
		int aciertosAnt=0;

		while(continua){
			vector<vector<int> > vecindario;
			generaVecinos(mascaraSecundaria,vecindario);
			bool seguirfor=true;
			int contadorEvaluaciones=0;
			
			for(int i=0;i<vecindario.size()&&seguirfor;i++){

				resultado=0;
				for(int j=0;j<numeroFilas;j++){
					if(knnClas(v,vNormalizado,j,vecindario[i],numeroFilas,numeroAtributos)==v[j][numeroAtributos-1]){
						resultado++;
						contadorEvaluaciones++;
					}
				}
				if(aciertosAnt<resultado){
					aciertosAnt=resultado;
					seguirfor=false;
					mascaraSecundaria=vecindario[i];
				}
				else if(contadorEvaluaciones>600)	//Para que no tarde tanto, cambiamos el maximo de evaluaciones de knn
					i=vecindario.size();			
			}
			if(seguirfor)
			continua=false;
		}

		//Si la salida de LS, no es mejor que la que teniamos anteriormente, mutamos la anterior y la guardamos para seguir con ella
		if(resultadoInicial>=resultado){
			mascaraSecundaria=mascaraInicial;
			mutador(mascaraInicial.size()*0.1,mascaraInicial);
		}
		else{ //Si mejora, mutamos esta y seguimos haciendo LS para esta guardandola como mejor
			mascaraInicial=mascaraSecundaria;
			resultadoInicial=resultado;
			mutador(mascaraSecundaria.size()*0.1,mascaraSecundaria);
		}
	}

	salida=mascaraSecundaria;
}

/**********************************************************/
/*            GENERADOR DE LLAMADAS PARA ILS               */
/**********************************************************/
void ils(int numSimulacion){
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
	iteradals(mitad1,vNormalizado,filasMitad1,numeroAtributos,salida);
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


int indiceMax(vector<int> &resultado){
	int mayor=0;
	int mayorvalor=resultado[0];

	for(int i=1;i<resultado.size();i++){
		if(resultado[i]>mayorvalor){
			mayor=i;
			mayorvalor=resultado[i];
		}
	}
	return mayor;
}
int indiceMin(vector<int> &resultado){
	int min=0;
	int minvalor=resultado[0];

	for(int i=1;i<resultado.size();i++){
		if(resultado[i]>minvalor){
			min=i;
			minvalor=resultado[i];
		}
	}
	return min;
}

int nMejorAleatorio(vector<int> resultado){

	int cMejor=resultado[indiceMax(resultado)];
	int cPeor=resultado[indiceMin(resultado)];
	vector<int> indicesMayores;

	for (int i=0;i<resultado.size();i++){
		int indice=indiceMax(resultado);
		if (resultado[indice]>=(cMejor-0.3*(cMejor-cPeor)))	{
			indicesMayores.push_back(indice);
		}
		resultado[indice]=-1;
	}

	return indicesMayores[(int)(Rand()*MASK)%indicesMayores.size()];

}

void generaGreedy(const vector<vector <double> > &v, vector<vector <double> > &vNormalizado, vector<int> &mascara, int numeroFilas, int numeroAtributos, vector<int> &s){
	vector<int> resultado;
	resultado.assign(numeroAtributos-1,0);

	//Inicializamos los valores del bucle a 0 y hacemos KNN para cada uno de los posibles
	resultado.assign(numeroAtributos,0);
	for(int i=0; i < numeroAtributos-1; ++i){
		if(mascara[i]==0){
			mascara[i] = !mascara[i];
			for(int j=0; j < numeroFilas; ++j){
				if(knnClas(v,vNormalizado,j,mascara,numeroFilas,numeroAtributos)==v[j][numeroAtributos-1])
					resultado[i]++;
			}

				mascara[i] = !mascara[i];
		}		
	}
    
    //Ya tenemos todos los aciertos de los generados por greedy, ahora llamamos a nMejorAleatorio que de los n mejores(0.3 de losposibles)devuelve un indice que invertimos el valor
    int numero=nMejorAleatorio(resultado);
	s[numero]=!s[numero];

}


/**********************************************************/
/*   IMPLEMENTACION DE BUSQUEDA GRASP                  */
/**********************************************************/
void grp(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, int numeroFilas,int numeroAtributos,vector<int> &salida){
	salida.assign(numeroAtributos-1,0);
	vector<int> mascaraGreedy;
	mascaraGreedy.assign(numeroAtributos-1,0);
	int resultadoInicial=0;
	int resultado=0;
	for(int i=0; i < 25; ++i){
		generaGreedy(v,vNormalizado,salida,numeroFilas,numeroAtributos,mascaraGreedy);

		bool continua=true;
		int aciertosAnt=0;

		while(continua){
			vector<vector<int> > vecindario;
			generaVecinos(mascaraGreedy,vecindario);
			bool seguirfor=true;
			int contadorEvaluaciones=0;
			
			for(int i=0;i<vecindario.size()&&seguirfor;i++){

				resultado=0;
				for(int j=0;j<numeroFilas;j++){
					if(knnClas(v,vNormalizado,j,vecindario[i],numeroFilas,numeroAtributos)==v[j][numeroAtributos-1]){
						resultado++;
						contadorEvaluaciones++;
					}
				}
				if(aciertosAnt<resultado){
					aciertosAnt=resultado;
					seguirfor=false;
					mascaraGreedy=vecindario[i];
				}
				else if(contadorEvaluaciones>600)	//Para que no tarde tanto, cambiamos el maximo de evaluaciones de knn
					i=vecindario.size();			
			}
			if(seguirfor)
			continua=false;
		}

		//Si la salida de LS, no es mejor que la que teniamos anteriormente, mutamos la anterior y la guardamos para seguir con ella
		if(resultado>resultadoInicial){
			salida=mascaraGreedy;
		}
	}

}

/**********************************************************/
/*            GENERADOR DE LLAMADAS PARA GRASP               */
/**********************************************************/
void grasp(int numSimulacion){
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
	grp(mitad1,vNormalizado,filasMitad1,numeroAtributos,salida);
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

/*	
	int opcionMenu;
	int numSimulacion;	
do{
cout << "Menu de acciones:"<<endl;
cout << "1.- SFS:"<<endl;
cout << "2.- MA:"<<endl;
cout << "3.- GRASP:"<<endl;
cout << "4.- ILS:"<<endl;
cout << "5.- SFS Completa:"<<endl;
cout << "6.- MA Completa:"<<endl;
cout << "7.- GRASP Completa:"<<endl;
cout << "8.- ILS:"<<endl;
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
			multiarranque(numSimulacion);
	    	break; 
	    case 3  :
			cout << "Numero de simulacion: 1-30 ";
			cin >> numSimulacion; 
			grasp(numSimulacion);
	    	break; 	
	    case 4  :
			cout << "Numero de simulacion: 1-30 ";
			cin >> numSimulacion; 
			ils(numSimulacion);
	    	break; 	
	    case 5  :
			cout << "Simulando";
			for(int i=1;i<31;i++)
			secuencial(i);
	    	break;
	    case 6  :
			cout << "Simulando";
			for(int i=1;i<31;i++)
			multiarranque(i);
	    	break; 
	    case 7  :
			cout << "Simulando";
			for(int i=1;i<31;i++)
			grasp(i);
	    	break; 	
	    case 8  :
			cout << "Simulando";
			for(int i=1;i<31;i++)
			ils(i);
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
*/

			for(int i=1;i<11;i++)
				grasp(i);
			for(int i=21;i<31;i++)
				grasp(i);
			grasp(11);
}