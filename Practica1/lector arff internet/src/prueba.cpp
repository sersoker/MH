#include "arff_parser.h"
#include <iostream>
using namespace std;


int main(){
	ArffParser parser("ejemplo.arff");
	ArffData *data = parser.parse();
	
	int numeroFilas=data->num_instances();
	int numeroAtributos=data->num_attributes();

	cout << data->num_attributes() << endl;
	cout << data->num_instances() << endl;
	
	ArffValue *valor;
	ArffInstance *instancia;
	std::vector<vector <double> > v;

for(int j=0;j<numeroFilas;j++){
	instancia = data->get_instance(j);
	std::vector<double> aux;

	for (int i = 0; i < data->num_attributes(); i++) {
		aux.push_back(std::stod(instancia->get(i)->operator string()) );
	}
	v.push_back(aux);
	aux.clear();
}

	for(int i=0;i<numeroFilas;i++){
		for(int j=0;j<numeroAtributos;j++){
			cout << v[i][j] <<" ";
		}
		cout << endl;
	}

}

