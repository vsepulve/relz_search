#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string.h>

#include <map>
#include <vector>

#include "RMQRMM64.h"

using namespace std;

int main(int argc, char* argv[]){
	
	cout << "Test - Start\n";
	
	
	
	
	
	
	RMQRMM64 *rmq = NULL;
	int *arr = NULL;
	unsigned int len = 10;
	unsigned int rmq_min;
	unsigned int beg, end;
	
	arr = new int[len];
	cout << "Test - Arreglo: ";
	for(unsigned int i = 0; i < len; ++i){
		arr[i] = rand() % 100;
		cout << arr[i] << " | ";
	}
	cout<<"\n";
	
	rmq = new RMQRMM64(arr, len);
	
	beg = 0;
	end = 10;
	rmq_min = rmq->queryRMQ(beg, end);
	cout << "Test - RMQ[" << beg << ", " << end << "]: " << rmq_min << " (" << arr[rmq_min] << ")\n";
	
	beg = 4;
	end = 8;
	rmq_min = rmq->queryRMQ(beg, end);
	cout << "Test - RMQ[" << beg << ", " << end << "]: " << rmq_min << " (" << arr[rmq_min] << ")\n";
	
	beg = 3;
	end = 6;
	rmq_min = rmq->queryRMQ(beg, end);
	cout << "Test - RMQ[" << beg << ", " << end << "]: " << rmq_min << " (" << arr[rmq_min] << ")\n";
	
	
	
	
	delete rmq;
	delete arr;
	
	cout << "Test - End\n";
	
}

















