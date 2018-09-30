#include "arff_parser.h"
#include <iostream>
using namespace std;


int main(){
	ArffParser parser("ejemplo.arff");
	ArffData *data = parser.parse();
	
	cout << data->num_attributes() << endl;
	cout << data->num_instances() << endl;
	
	ArffInstance *instancia = data->get_instance(0);
	ArffValue *valor;
	

	for (int i = 1; i < data->num_attributes(); i++) {
		valor = instancia->get(i);
		cout << valor->operator float() << " ";
	}

}