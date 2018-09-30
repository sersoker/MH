#include "arff_parser.h"
#include <vector>
#include <iostream>
#include <string>
#include <numeric>
#include <functional>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Definicion del generador de numeros aleatorios */
#define MASK 2147483647
#define PRIME 65539
#define SCALE 0.4656612875e-9
#define Rand()  (( Seed = ( (Seed * PRIME) & MASK) ) * SCALE )

using namespace std;
unsigned long Seed=2; /* semilla del generador de numeros aleatorios */

int Randint(int low, int high){
	return (int)(low+(high-(low)+1)*Rand());
}


/**********************************************************/
/*            CALCULOS EUCLIDEOS                          */
/*********************************************************/
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
	vector<double> aux;

	for(int i=0;i<numeroFilas;i++){
		aux.clear();

		for(int j=0;j<numeroAtributos-1;j++){
			if(mask[j]!=0)
				aux.push_back(vNormalizado[i][j]*mask[j]);
			else
				aux.push_back(0);
		}
		vConMascara.push_back(aux);
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
int indiceMinimo(vector<int> &resultado){
	int indice=0;
	for(int i=0;i<resultado.size();i++)
			if(resultado[i]<resultado[indice])
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
	//vectorDistancias[ejemplo]=99999;
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



/********************************************************************************************************************************************************************************************/
/**********************************************************/
/* Calculo de un coste para una mascaara  */
/**********************************************************/
int calculoCoste(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, vector<int> &mask,int numeroFilas,int numeroAtributos){
	int resultado=0;
	for(int j=0;j<numeroFilas;j++){
		if(knnClas(v,vNormalizado,j,mask,numeroFilas,numeroAtributos)==v[j][numeroAtributos-1]){
			resultado++;
			}
	}	

	return resultado;
}
/**********************************************************/
/*            DADA UNA MASCARA LE HACE UN SWAP ALEATORIO  */
/**********************************************************/
void generaVecinoAleatorio(vector<int> &mask){

	int posicion= Randint(0,mask.size()-1);
	mask[posicion]=!mask[posicion];

}
/**********************************************************/
/*           GENERADOR DE POBLACION              */
/**********************************************************/
void generaEntorno(vector<int> &mask,vector< vector<int> > &resultado){
	resultado.clear();
	for(int i=0;i<30;i++){
		auto aux=mask;
		generaVecinoAleatorio(aux);
		resultado.push_back(aux);
		
	}
}

/**********************************************************/
/*           GENERADOR DE CRUcE                          */
/**********************************************************/
void cruzar(vector<int> &cruzado1,vector<int> &cruzado2){
	auto aux(cruzado1);
	int punto1;
	int punto2;

  /* Se generan los dos puntos de cruce*/
	punto1=Randint (0,cruzado1.size());
	if (punto1!=cruzado1.size())
		punto2=Randint (punto1+1,cruzado1.size());
	else
		punto2=cruzado1.size();

	for(int i=punto1;i<punto2;i++){
		cruzado1[i]=cruzado2[i];
		cruzado2[i]=aux[i];
	}


}


/**********************************************************/
/*           GENERADOR DE CRUcE                          */
/**********************************************************/
void mutador(vector<int> &mask){
	int punto1;
	punto1=Randint (0,mask.size());
	mask[punto1]=!mask[punto1];

}

/**********************************************************/
/*            IMPLEMENTACION AGG                          */
/**********************************************************/
void aggAlgotirmo(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, vector<int> &mask,int numeroFilas,int numeroAtributos,vector<int> &resultado){

	vector< vector<int> > poblacionInicial;
	vector< vector<int> > torneo;
	vector<int> costes;
	vector<int> costesTorneo;
	mask.assign(0,numeroAtributos-1);
	int num_eval=0;
	int numeroMutaciones=(30*numeroAtributos-1)*0.001;

	//Iniciamos la mascara aleatoria
		for(int i=0;i<numeroAtributos-1;i++){
			if(Rand()>0.5)		mask.push_back(0);
			else				mask.push_back(0);
		}

	generaEntorno(mask,poblacionInicial);

	for(int i=0;i<poblacionInicial.size();i++){
		costes.push_back(calculoCoste(v,vNormalizado,poblacionInicial[i],numeroFilas,numeroAtributos));
	}

	while(num_eval<15000){
		//REALIZAMOS UN TORNEO
		torneo.clear();
		costesTorneo.clear();
		
		int aleatorio1=0;
		int aleatorio2=0;
		for(int i=0;i<30;i++){
			do{
				aleatorio1=Randint(0,29);
				aleatorio2=Randint(0,29);
			}while(aleatorio2==aleatorio1);
			
			if(costes[aleatorio2]>costes[aleatorio1])				{torneo.push_back(poblacionInicial[aleatorio2]);costesTorneo.push_back(costes[aleatorio2]);}
			else													{torneo.push_back(poblacionInicial[aleatorio1]);costesTorneo.push_back(costes[aleatorio1]);}
		}

		//Hacemos el cruce del 0,7 de los ganadores del torneo
		for(int i=0;i<22;i+=2){
			cruzar(torneo[i],torneo[i+1]);
			costesTorneo[i]=calculoCoste(v,vNormalizado,torneo[i],numeroFilas,numeroAtributos);
			num_eval+=numeroFilas;
			costesTorneo[i+1]=calculoCoste(v,vNormalizado,torneo[i+1],numeroFilas,numeroAtributos);
			num_eval+=numeroFilas;
		}
		for(int i=0;i<numeroMutaciones;i++){

			aleatorio1=Randint(0,29);
			mutador(torneo[aleatorio1]);

			costesTorneo[aleatorio1]=calculoCoste(v,vNormalizado,torneo[aleatorio1],numeroFilas,numeroAtributos);
			num_eval+=numeroFilas;	
		}
		auto sustituir=poblacionInicial[indiceMaximo(costes)];
		auto costesustituir=costes[indiceMaximo(costes)];
		costes=costesTorneo;
		poblacionInicial=torneo;
		poblacionInicial[poblacionInicial.size()-1]=sustituir;
		costes[poblacionInicial.size()-1]=costesustituir;
	}

	resultado=poblacionInicial[indiceMaximo(costes)];

}



/**********************************************************/
/*            IMPLEMENTACION AGE                          */
/**********************************************************/
void ageAlgotirmo(const vector<vector <double> > &v,vector<vector <double> > &vNormalizado, vector<int> &mask,int numeroFilas,int numeroAtributos,vector<int> &resultado){

	vector< vector<int> > poblacionInicial;
	vector<int> costes;
	mask.assign(0,numeroAtributos-1);
	int num_eval=0;
	int numeroMutaciones=(30*numeroAtributos-1)*0.001;

	//Iniciamos la mascara aleatoria
		for(int i=0;i<numeroAtributos-1;i++){
			if(Rand()>0.5)		mask.push_back(0);
			else				mask.push_back(0);
		}

	generaEntorno(mask,poblacionInicial);

	for(int i=0;i<poblacionInicial.size();i++){
		costes.push_back(calculoCoste(v,vNormalizado,poblacionInicial[i],numeroFilas,numeroAtributos));
	}

	while(num_eval<15000){
		//REALIZAMOS UN TORNEO
		vector<int> ganador1,ganador2;
		int costeGanador1,costeGanador2;
		
		int aleatorio1=0;
		int aleatorio2=0;
			do{
				aleatorio1=Randint(0,29);
				aleatorio2=Randint(0,29);
			}while(aleatorio2==aleatorio1);
			
			if(costes[aleatorio2]>costes[aleatorio1])				{ganador1=poblacionInicial[aleatorio2];costeGanador1=costes[aleatorio2];}
			else													{ganador1=poblacionInicial[aleatorio1];costeGanador1=costes[aleatorio1];}

			do{
				aleatorio1=Randint(0,29);
				aleatorio2=Randint(0,29);
			}while(aleatorio2==aleatorio1);
			
			if(costes[aleatorio2]>costes[aleatorio1])				{ganador2=poblacionInicial[aleatorio2];costeGanador2=costes[aleatorio2];}
			else													{ganador2=poblacionInicial[aleatorio1];costeGanador2=costes[aleatorio1];}

		//Hacemos el cruce de los ganadores del torneo
			cruzar(ganador1,ganador2);

		if(Rand()<=0.001)
			mutador(ganador1);
		if(Rand()<=0.001)
			mutador(ganador2);

		costeGanador1=calculoCoste(v,vNormalizado,ganador1,numeroFilas,numeroAtributos);
		num_eval+=numeroFilas;	
		costeGanador2=calculoCoste(v,vNormalizado,ganador2,numeroFilas,numeroAtributos);
		num_eval+=numeroFilas;	

		if(costes[indiceMinimo(costes)]<costeGanador1){
			costes[indiceMinimo(costes)]=costeGanador1;
			poblacionInicial[indiceMinimo(costes)]=ganador1;
		}
		if(costes[indiceMinimo(costes)]<costeGanador2){
			costes[indiceMinimo(costes)]=costeGanador2;
			poblacionInicial[indiceMinimo(costes)]=ganador2;
		}
	}

	resultado=poblacionInicial[indiceMaximo(costes)];

}





/********************************************************************************************************************************************************************************************/
/**********************************************************/
/*            GENERADOR DE LLAMADAS PARA AGE              */
/**********************************************************/
void age(int numSimulacion){
	cout << "LANZANDO AGE, Numero de simulacion:  "<< numSimulacion <<endl;
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

	vector<int> mask;
	vector<int> resultado;

	ageAlgotirmo(mitad1,vNormalizado,mask,filasMitad1,numeroAtributos,resultado);
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
/*            GENERADOR DE LLAMADAS PARA AGG              */
/**********************************************************/
void agg(int numSimulacion){
	cout << "LANZANDO AGG, Numero de simulacion:  "<< numSimulacion <<endl;
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

	vector<int> mask;
	vector<int> resultado;

	aggAlgotirmo(mitad1,vNormalizado,mask,filasMitad1,numeroAtributos,resultado);
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








/*************************************************************************************/
/*           MAIN, PARA MENU Y LLAMADAS A LAS DIFERENTES FUNCIONES Y SIMULACIONES    */
/*************************************************************************************/

int main(){

	int opcionMenu;
	int numSimulacion;	

			for(int i=1;i<11;i++){
				agg(i);
				age(i);
			}


cout << "Menu de acciones:"<<endl;
cout << "1.- AGG:"<<endl;
cout << "2.- AGE:"<<endl;
cout << "3.- AYUDA"<<endl;
cout << "0.- SALIR:"<<endl;
cout << "---------:";

	cin >> opcionMenu;
	switch(opcionMenu){
	    case 1  :
			cout << "Numero de simulacion: 1-30 ";
			cin >> numSimulacion;
			agg(numSimulacion);
	    	break;
	    case 2  :
			cout << "Numero de simulacion: 1-30 ";
			cin >> numSimulacion; 
			age(numSimulacion);
	    	break; 
	    case 3  :
			cout << "Las simulaciones van de la 1 a la 30, divididas en 3 partes:"<<endl;
			cout << "1-10- Libras:"<<endl;
			cout << "11-20- Arrhythmia:"<<endl;
			cout << "21-30- Wdbc0:"<<endl;
	    	break; 	
	    case 0  :
			cout << "SALIENDO..."<<endl;
	    	break; 

	    default : //Optional
	       cout << "OPCION INCORRECTA"<<endl;
	}
}