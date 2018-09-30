#include "arff_parser.h"
#include "random_ppio.h"
#include <vector>
#include <iostream>
#include <string>
#include <numeric>
#include <functional>
#include <math.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////
// Función para calcular la distancia euclídea entre dos vectores
/////////////////////////////////////////////////////////////////////////////////
template<class Iter_T, class Iter2_T>
double distanciaVectores(Iter_T first, Iter_T last, Iter2_T first2) {
  double ret = 0.0;
  while (first != last) {
    double dist = (*first++) - (*first2++);
    ret += dist * dist;
  }
  return ret > 0.0 ? sqrt(ret) : 0.0;
}


/////////////////////////////////////////////////////////////////////////////////
// Lector ARFF
/////////////////////////////////////////////////////////////////////////////////
void lectorARFF(int idArchivo,	std::vector<vector <double> > &v,int &numFilas,int &numAtributos){
	string nombreArchivo;
	switch(idArchivo){
	    case 1:			
	    	nombreArchivo = "./datos/libras00.arff";		    
	    	break;

	    case 2:			
	    	nombreArchivo = "./datos/libras01.arff";			
	    	break;

	    case 3:			
	    	nombreArchivo = "./datos/libras02.arff";			
	    	break;

	    case 4:			
	    	nombreArchivo = "./datos/libras03.arff";			
	    	break;
	    
	    case 5:			
	    	nombreArchivo = "./datos/libras04.arff";			
	    	break;
	    
	    case 6:			
	    	nombreArchivo = "./datos/libras05.arff";			
	    	break;
	    
	    case 7:			
	    	nombreArchivo = "./datos/libras06.arff";			
	    	break;
	    
	    case 8:			
	    	nombreArchivo = "./datos/libras07.arff";			
	    	break;
	    
	    case 9:			
	    	nombreArchivo = "./datos/libras08.arff";			
	    	break;
	    
	    case 10:			
	    	nombreArchivo = "./datos/libras09.arff";			
	    	break;

	    case 11:			
	    	nombreArchivo = "./datos/arrhythmia00.arff";		    
	    	break;
	    
	    case 12:			
	    	nombreArchivo = "./datos/arrhythmia01.arff";			
	    	break;
	    
	    case 13:			
	    	nombreArchivo = "./datos/arrhythmia02.arff";			
	    	break;
	    
	    case 14:			
	    	nombreArchivo = "./datos/arrhythmia03.arff";			
	    	break;
	    
	    case 15:			
	   		nombreArchivo = "./datos/arrhythmia04.arff";			
	    	break;
	    
	    case 16:			
	    	nombreArchivo = "./datos/arrhythmia05.arff";			
	    	break;
	    
	    case 17:			
	    	nombreArchivo = "./datos/arrhythmia06.arff";			
	    	break;
	    
	    case 18:			
	    	nombreArchivo = "./datos/arrhythmia07.arff";			
	    	break;
	    
	    case 19:			
	    	nombreArchivo = "./datos/arrhythmia08.arff";			
	    	break;
	    
	    case 20:			
	    	nombreArchivo = "./datos/arrhythmia09.arff";			
	    	break;

	    case 21:			
	    	nombreArchivo = "./datos/wdbc00.arff";		    
	    	break;
	    
	    case 22:			
	    	nombreArchivo = "./datos/wdbc01.arff";			
	    	break;
	    
	    case 23:			
	    	nombreArchivo = "./datos/wdbc02.arff";			
	    	break;
	    
	    case 24:			
	    	nombreArchivo = "./datos/wdbc03.arff";			
	    	break;
	    
	    case 25:			
	    	nombreArchivo = "./datos/wdbc04.arff";			
	    	break;
	    
	    case 26:			
	    	nombreArchivo = "./datos/wdbc05.arff";			
	    	break;
	    
	    case 27:			
	    	nombreArchivo = "./datos/wdbc06.arff";			
	    	break;
	    
	    case 28:			
	    	nombreArchivo = "./datos/wdbc07.arff";			
	    	break;
	    
	    case 29:			
	    	nombreArchivo = "./datos/wdbc08.arff";			
	    	break;
	    
	    case 30:			
	    	nombreArchivo = "./datos/wdbc09.arff";			
	    	break;
	}
	
	ArffParser parser(nombreArchivo);
	ArffData *data = parser.parse();
	
	numFilas = data->num_instances();
	numAtributos = data->num_attributes();

	ArffValue *valor;
	ArffInstance *instancia;

	for(int j=0;j<numFilas;++j){
		instancia = data->get_instance(j);
		std::vector<double> aux;

		for (int i = 0; i < data->num_attributes(); ++i) {
			aux.push_back(stod(instancia->get(i)->operator string()) );
		}
		v.push_back(aux);
		aux.clear();
	}

}

/////////////////////////////////////////////////////////////////////////////////
// Función para calcular el mínimo y el máximo de cada columna 
/////////////////////////////////////////////////////////////////////////////////
void calcularMinMax(vector<vector <double> > &v,vector<double> &minimoColumnas,vector<double>  &maximoColumnas,int numFilas,int numAtributos){
	maximoColumnas = minimoColumnas = v[0];
	for(int i=0; i < numFilas; ++i){
		for(int j=0; j < numAtributos; ++j){
			if(v[i][j] > maximoColumnas[j])
				maximoColumnas[j] = v[i][j];
			if(v[i][j] < minimoColumnas[j])
				minimoColumnas[j] = v[i][j];
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////
// Función para normalizar el 
/////////////////////////////////////////////////////////////////////////////////
void normalizar(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado,int numFilas,int numAtributos,vector<double> &minimoColumnas,vector<double>  &maximoColumnas){
	
	std::vector<double> aux;
	for(int i=0; i < numFilas; ++i){
		for(int j=0; j < numAtributos-1; ++j){
			if(isnan((v[i][j]-minimoColumnas[j])/(maximoColumnas[j]-minimoColumnas[j])))
				aux.push_back(0);
			else
			aux.push_back((v[i][j]-minimoColumnas[j])/(maximoColumnas[j]-minimoColumnas[j]));
		}

		vNormalizado.push_back(aux);
		aux.clear();
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Función para aplicar la máscara a los datos 
/////////////////////////////////////////////////////////////////////////////////
void aplicarMascara(vector<vector <double> > &vNormalizado,vector<vector <double> > &vConMascara,vector<int> &mascara,int numFilas,int numAtributos){
	
	std::vector<double> aux;

	for(int i=0; i < numFilas; ++i){
		for(int j=0; j < numAtributos; ++j){
			aux.push_back(vNormalizado[i][j]*mascara[j]);
		}

		vConMascara.push_back(aux);
		aux.clear();
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Función para calcular la distancia euclidea entre todas las filas
/////////////////////////////////////////////////////////////////////////////////
void distanciaEuclidea(vector<double> &ejemploElegido,vector<vector <double> > &vConMascara,vector<double> &vectorDistancias,int numFilas,int numAtributos){
	
	vector<double> aux;
	for(int i=0; i < numFilas; ++i){
		aux = vConMascara[i];
		vectorDistancias.push_back(distanciaVectores(aux.begin(), aux.end(), ejemploElegido.begin()));
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Función para calcular los 3 valores más cercanos (de menor distancia)
/////////////////////////////////////////////////////////////////////////////////
void calcularMenores(vector<double> &vectorDistancias,int &indice1, int &indice2, int &indice3){
	
	indice1 = 0;

	for(int i=0; i < vectorDistancias.size(); ++i){

		if(vectorDistancias[i] <= vectorDistancias[indice1])
			indice1 = i;		
	}

	if(indice1 != 0)
		indice2 = 0;

	else
		indice2 = 1;


	for(int i=0; i < vectorDistancias.size(); ++i){

		if(vectorDistancias[i] <= vectorDistancias[indice2] && i != indice1)
			indice2 = i;		
	}

	if(indice1 != 0 && indice2 != 0)
		indice3 = 0;
	
	else{

		if(indice1 != vectorDistancias.size()-1 && indice2 != vectorDistancias.size()-1)
			indice3 = vectorDistancias.size()-1;
		else if(indice1 == vectorDistancias.size()-1)
			indice3 = indice2+1;
		else
			indice3 = indice1+1;

	}

	
	for(int i=0; i < vectorDistancias.size(); ++i){

		if(vectorDistancias[i] <= vectorDistancias[indice3] && i != indice1 && i != indice2)
			indice3 = i;		
	}
	
}

/////////////////////////////////////////////////////////////////////////////////
// Función para calcular el índice de máximo valor
/////////////////////////////////////////////////////////////////////////////////
int indiceMax(vector<int> &resultado){

	int indice = 0;
	for(int i=0; i < resultado.size(); ++i)
		if(resultado[i] > resultado[indice])
			indice = i;
	return indice;
}

/////////////////////////////////////////////////////////////////////////////////
// Clasificador knn
/////////////////////////////////////////////////////////////////////////////////
int knn(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado,int ejemplo, vector<int> &mascara,int numFilas,int numAtributos){
	
	//Aplicamos la máscara
	vector<vector <double> > vConMascara;
	aplicarMascara(vNormalizado,vConMascara,mascara,numFilas,numAtributos);

	vector<double> ejemploElegido = vConMascara[ejemplo];

	//Calculamos las distancias 
	vector<double> vectorDistancias;
	distanciaEuclidea(ejemploElegido,vConMascara,vectorDistancias,numFilas,numAtributos);

	//Calculamos las 3 distancias menores
	int indice1, indice2,  indice3;
	calcularMenores(vectorDistancias,indice1, indice2, indice3);

	//Devolvemos la clase en la que más coincidan
	if(v[indice1][numAtributos-1] == v[indice2][numAtributos-1])
		return v[indice1][numAtributos-1];
	if(v[indice2][numAtributos-1] == v[indice3][numAtributos-1])
		return v[indice2][numAtributos-1];
	if(v[indice3][numAtributos-1] == v[indice1][numAtributos-1])
		return v[indice3][numAtributos-1];

	return v[indice1][numAtributos-1];
}

/////////////////////////////////////////////////////////////////////////////////
// Algoritmo SFS
/////////////////////////////////////////////////////////////////////////////////
void sfs(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, vector<int> &mascara,int numFilas,int numAtributos,vector<int> &resultado){

	bool fin = false;
	int aciertosAnt = 0;
	vector<int> mascaraAuxiliar = resultado;

	while(!fin){

		//Inicializamos los valores del bucle a 0
		resultado.assign(vNormalizado.size(),0);

		//Calculamos los aciertos de cada clase
		for(int i=0; i < numAtributos-1; ++i){
				mascara[i] = 1;
			for(int j=0; j < numFilas; ++j){
				if(knn(v,vNormalizado,j,mascara,numFilas,numAtributos) == v[j][numAtributos-1])
					resultado[i] += 1;
			}
			if(mascaraAuxiliar[i] != 1)
				mascara[i] = 0;
		}

		//Cogemos la característica más prometedora
		if(resultado[indiceMax(resultado)] > aciertosAnt){
	        mascaraAuxiliar[indiceMax(resultado)] = 1;
	        mascara[indiceMax(resultado)] = 1;
	        aciertosAnt = resultado[indiceMax(resultado)];  
      	}

      	//Si no mejoramos, terminamos
	    else
	    	fin = true;
		
	}

	resultado = mascaraAuxiliar;
}

/////////////////////////////////////////////////////////////////////////////////
// Función para generar vecinos
/////////////////////////////////////////////////////////////////////////////////
void generarVecinos(vector<int> mascara,vector< vector<int> >&resultado){

	for(int i=0; i < mascara.size(); ++i){
		if(mascara[i] == 0)
			mascara[i] = 1;
		else
			mascara[i] = 0;

		resultado.push_back(mascara);
		
		if(mascara[i] == 0)
			mascara[i] = 1;
		else
			mascara[i] = 0;
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Algoritmo LS
/////////////////////////////////////////////////////////////////////////////////
void ls(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, int numFilas,int numAtributos,vector<int> &resultado){
	
	vector<int> mascara;
	int coste;
	
	//Generamos una máscara inicial aleatoria
	for(int i=0;i < numAtributos-1;++i){
		mascara.push_back(Randint(0,1));
	}

	bool fin = false;
	int aciertosAnt = 0;

	while(!fin){

		//Generamos los vecinos
		vector<vector<int> > vecindario;
		generarVecinos(mascara,vecindario);
		bool continua = true;
		
		for(int i=0; i < vecindario.size() && continua; ++i){
			
			coste = 0;

			//Calculamos el coste de cada vecino usando KNN
			for(int j=0; j < numFilas; ++j){
				if(knn(v,vNormalizado,j,vecindario[i],numFilas,numAtributos) == v[j][numAtributos-1])
					coste++;
			}

			//Si encontramos una solución de mejor coste, dejamos de evaluar los vecinos y actualizamos la solución
			if(aciertosAnt < coste){
				aciertosAnt = coste;
				continua = false;
				mascara = vecindario[i];
			}				
		}

		//Si hemos evaluado todos los vecinos y no obtenemos ninguna mejora, terminamos
		if(continua)
			fin = true;
	}
}


/////////////////////////////////////////////////////////////////////////////////
// Función para comprobar la probabilidad de aceptación de una solución
/////////////////////////////////////////////////////////////////////////////////
int metrop (float coste, float t){
	
	return coste < 0.0 || Rand() < exp(-coste/t);
}

/////////////////////////////////////////////////////////////////////////////////
// Función para generar un vecino de forma aleatoria
/////////////////////////////////////////////////////////////////////////////////
void generarVecinoAleatorio(vector<int> &mascara){

	int posicion =  ((int)Rand())%mascara.size();
	mascara[posicion] = !mascara[posicion];
}

/////////////////////////////////////////////////////////////////////////////////
// Algoritmo SA
/////////////////////////////////////////////////////////////////////////////////
void sa(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, int numFilas,int numAtributos,vector<int> &resultado){
	
	vector<int> mascara;
	int coste;
	int aciertosAnt = 0;

	//Generamos una primera mascara aleatoria
	for(int i=0; i < numAtributos-1; ++i){
		cout << Randint(0,99);
		mascara.push_back(Randint(0,1));
	}

	//Calculamos el coste usando KNN
	coste = 0;
	
	for(int j=0; j < numFilas; ++j){
		if(knn(v,vNormalizado,j,mascara,numFilas,numAtributos) == v[j][numAtributos-1])
			coste++;
	}

	aciertosAnt = coste;

	float tInicial = (0.3*aciertosAnt)/-log(0.3);
	float tActual = tInicial;
	float tFinal = 0.001;

	int maxVecinos = numFilas;
	int maxExitos = 0.1*maxVecinos;
	int evaluaciones = 15000;
	int M = evaluaciones/(maxVecinos*maxVecinos);
	float beta = (tInicial-tFinal)/M*tInicial*tFinal;
	int exitos;

	do{
		int numVecinos = 0;
		exitos = 0;

		for(int i=0; i < evaluaciones; ++i){
			
			//Generamos un vecino
			auto aux = mascara;
			generarVecinoAleatorio(aux);
			numVecinos++;

			//Calculamos el coste de la nueva solución
			coste = 0;
			for(int j=0; j < numFilas; ++j){
				if(knn(v,vNormalizado,j,aux,numFilas,numAtributos) == v[j][numAtributos-1])
					coste++;
			}

			//Si la nueva solución es mejor, actualizamos la solución
			if(aciertosAnt < coste){
				aciertosAnt = coste;
				mascara = aux;
				exitos++;
			}
			
			//Si no es mejor, comprobamos la probabilidad de aceptación de la solución
			else if(metrop(coste-aciertosAnt,tActual) == 1){
				
				aciertosAnt = coste;
				mascara = aux;
				exitos++;		
			}

			//Si se han generado todos los vecinos o el número de éxitos necesarios, salimos del bucle
		    if(numVecinos == maxVecinos || exitos == maxExitos){
		    	i = evaluaciones;
		    }
		}

		cout << "Tactual: " << tActual << endl;

	    //Enfriamos la temperatura
	    tActual  = tActual/(1+beta*tActual);	

	}while(tActual > tFinal && exitos > 0);
		
	resultado = mascara;
}


/////////////////////////////////////////////////////////////////////////////////
// Función para generar entornos de tamaño tam para TB
/////////////////////////////////////////////////////////////////////////////////
void generarEntorno(vector<int> &mascara,vector< vector<int> > &resultado,int tam){
  
	resultado.clear();

	for(int i=0; i < tam; ++i){
		auto aux = mascara;
		generarVecinoAleatorio(aux);
		resultado.push_back(aux);
	}
}

/////////////////////////////////////////////////////////////////////////////////
// Función para encontrar el primer índice en el que los valores de dos máscaras dadas difieren
/////////////////////////////////////////////////////////////////////////////////
int indiceDistinto(vector<int> &mascara1,vector<int> &mascara2){

  for(int i=0; i < mascara1.size(); ++i)
    if(mascara1[i] != mascara2[i])
      return i;
}

/////////////////////////////////////////////////////////////////////////////////
// Algoritmo TB
/////////////////////////////////////////////////////////////////////////////////
void tb(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, int numFilas,int numAtributos,vector<int> &resultado){
	
	//Calculamos una primera mascara aleatoria
	vector<int> mascara;
	int coste=0;
	int aciertosAnt = 0;
	int mejorAcierto = 0;
	int tenencia = numAtributos/30;
	int iteraciones = 500;
	vector<int> listaTabu(numAtributos-1,0); 


	for(int i=0; i < numAtributos-1; ++i){
		mascara.push_back(Randint(0,1));
	}	

	//Calculamos el coste usando KNN 
	for(int j=0;j<numFilas;++j){
		if(knn(v,vNormalizado,j,mascara,numFilas,numAtributos) == v[j][numAtributos-1])
			coste++;
	}

	aciertosAnt = mejorAcierto = coste;
	vector< vector<int> > entorno;

	for(int l=0; l < iteraciones; ++l){

		//Generamos el entorno con la mascara actual
		generarEntorno(mascara,entorno,30);
		int indiceMejorVecino;
		int resultadoMejorVecino = 0;

		//Calculamos el mejor vecino del entorno generado
		for(int i=0; i < entorno.size(); ++i){

			//Calculamos el coste usando KNN
			coste = 0;
			for(int j=0; j < numFilas; ++j){
		   		if(knn(v,vNormalizado,j,entorno[i],numFilas,numAtributos) == v[j][numAtributos-1])
		      		coste++;
		    }
		  
		  	//Seleccionamos el mejor vecino
		  	if(coste > resultadoMejorVecino){
		    	resultadoMejorVecino = coste;
		    	indiceMejorVecino = i;
		  	}	
		}

		//Comprobamos el índice que cambia con respecto a la mascara que tenemos actualmente
		int indice =  indiceDistinto(mascara,entorno[indiceMejorVecino]);

		//Si el índice no está marcado como tabú
		if(listaTabu[indice] == 0){

			//Actualizamos los valores de tenencia de los elementos de la lista tabú
		  	for(int k=0; k < listaTabu.size(); ++k){
		    	if(listaTabu[k] > 0) 
		    		--listaTabu[k];
		  	}

		  	//Actualizamos la tenencia del índice y la máscara
		  	listaTabu[indice] = tenencia;
		  	mascara = entorno[indiceMejorVecino];
		  	aciertosAnt = coste;

		  	if(coste > mejorAcierto)
		   		mejorAcierto = coste;
		}

		//Si el índice esa marcado como tabú, comprobamos el criterio de aspiración
		else{
		  	
		  	//Si el coste es mejor que el del mejor acierto (global), lo actualizamos
		  	if(coste > mejorAcierto){
		    	mejorAcierto = coste;

		    	//Actualizamos los valores de tenencia de los elementos de la lista tabú
		    	for(int k=0; k < listaTabu.size(); ++k){
		      		if(listaTabu[k]>0) 
		      			listaTabu[k]--;
		    	}
		    	
		    	//Actualizamos la tenencia del índice y la máscara
		    	listaTabu[indice]+= tenencia;
		    	mascara = entorno[indiceMejorVecino];
		    	aciertosAnt = coste;
		  	}
		}
	}

	resultado = mascara;
}


int main(){

	int opcionMenu;
	int numSimulacion;

	cout << "-----------------" << endl;
	cout << "Menú:" << endl;
	cout << "Opción 1: SFS" << endl;
	cout << "Opción 2: LS:" << endl;
	cout << "Opción 3: SA:" << endl;
	cout << "Opción 4: TS:" << endl;
	cout << "Opción 5: salir" << endl;
	cout << "-----------------" << endl;
	cin >> opcionMenu;
	
	while(opcionMenu != 5){
	
		cout << "Introduzca el número de simulación(1-30): " << endl;
		cin >> numSimulacion;

		std::vector<vector <double> > vectorEntero;
		std::vector<vector <double> > mitad1;
		std::vector<vector <double> > mitad2;

		int numFilasI;
		int numAtributos;
		lectorARFF(numSimulacion,vectorEntero,numFilasI,numAtributos);
		//SepaRando los datos en dos vectores, mitad en cada uno
		for(int i=0;i<numFilasI;++i){
			if(i<numFilasI/2)
				mitad1.push_back(vectorEntero[i]);
			else
				mitad2.push_back(vectorEntero[i]);
		}

		//Guardamos el numero de filas de cada una de las matrices
		int filasMitad1 =  numFilasI-numFilasI/2-1;
		int filasMitad2 =  numFilasI-filasMitad1;

		vector<double> minimoColumnas;
		vector<double> maximoColumnas;

		clock_t startTime = clock();
		std::vector<vector <double> > vNormalizado;
		std::vector<vector <double> > vNormalizado2;
		calcularMinMax(mitad1,minimoColumnas,maximoColumnas,filasMitad1,numAtributos);
		normalizar(mitad1,vNormalizado,filasMitad1,numAtributos,minimoColumnas,maximoColumnas);
		vector<int> mascara;
		vector<int> resultado;
		mascara.assign(numAtributos,0);
		resultado.assign(numAtributos,0);

		switch(opcionMenu){
			case 1:
				cout << "Ejecutando SFS. Simulación: " << numSimulacion << endl;
				sfs(mitad1,vNormalizado,mascara,filasMitad1,numAtributos,resultado);
				break;

			case 2:
				cout << "Ejecutando LS. Simulación: " << numSimulacion << endl;
				ls(mitad1,vNormalizado,filasMitad1,numAtributos,resultado);
				break;

			case 3:
				cout << "Ejecutando SA. Simulación: " << numSimulacion << endl;
				sa(mitad1,vNormalizado,filasMitad1,numAtributos,resultado);
				break;

			case 4:
				cout << "Ejecutando TS. Simulación: " << numSimulacion << endl;
				tb(mitad1,vNormalizado,filasMitad1,numAtributos,resultado);
				break;
		}
 
		cout << "Resultado: " << endl;
		for(int i=0;i<numAtributos-1;++i)
			cout << resultado[i] << " ";
		cout << endl;

		int aciertos = 0;
		int fallos = 0;
		calcularMinMax(mitad2,minimoColumnas,maximoColumnas,mitad2.size(),numAtributos);
		normalizar(mitad2,vNormalizado2,mitad2.size(),numAtributos,minimoColumnas,maximoColumnas);
		for(int i=0;i<mitad2.size();++i)
			if(knn(mitad2,vNormalizado2,i,resultado,mitad2.size(),numAtributos) == mitad2[i][numAtributos-1])
				aciertos++;
			else
				fallos++;

		cout << "Número de aciertos: " << aciertos <<  "\tNúmero de fallos: " << fallos << endl << "Porcentaje de acierto: " << aciertos*100/(aciertos+fallos) << "%" << endl;
		cout << "Tiempo de ejecución: " << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " segundos." << endl;

		cout << "-----------------" << endl;
		cout << "Menú:" << endl;
		cout << "Opción 1: SFS" << endl;
		cout << "Opción 2: LS:" << endl;
		cout << "Opción 3: SA:" << endl;
		cout << "Opción 4: TS:" << endl;
		cout << "Opción 5: salir" << endl;
		cout << "-----------------" << endl;
		cin >> opcionMenu;

	}

	cout << "Saliendo..." << endl;
}
