#include<stdio.h>
#include<math.h>
#include<mpi.h>

double func(double);
double approximateArea(int, int);
double integrationPTP(int, int);
double integrationCC(int, int);


int main(void){
    const int REPEAT = 25;
    int myRank, commSz;

	MPI_Init(NULL, NULL);
	
    MPI_Comm_size(MPI_COMM_WORLD, &commSz);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	
    printf("Hello I am process %d of %d\n", myRank, commSz);

    double tStartPTP, tEndPTP, tStartCC, tEndCC;
	
    int count = 0;

	double resultPTP = -1, resultCC = -1; //if "-1" then the function has not been called 
    for(count; count < REPEAT; count++){
		
	
    	tStartPTP = MPI_Wtime();
    	resultPTP = integrationPTP(commSz, myRank);
		tEndPTP = MPI_Wtime();

		tStartCC = MPI_Wtime();
		//resultCC = integrationPTP(commSz, myRank);
		tEndCC = MPI_Wtime();

		if (myRank == 0){
			printf("====================\nCOUNT = %d\n\n", count);

			printf("PTP Result = %f\n", resultPTP);
			printf("Time taken PTP = %f milliseconds\n", 1000.0*(tEndPTP - tStartPTP));	//print wallTimeTaken

			printf("CC Result = %f\n", resultCC);
			printf("Time taken PTP = %f milliseconds\n", 1000.0*(tEndCC - tStartCC));
		}
	}
	MPI_Finalize();
}

double integrationPTP(int commSz, int myRank){
	
	double mySum = approximateArea(commSz, myRank);
	
	if(myRank != 0){
		MPI_Send(&mySum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	else{
		double recvSum;
		for(int i = 1; i < commSz; i++){
			MPI_Recv(&recvSum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			mySum += recvSum;
		}
	}

	return mySum;
}

/*TASK: Create a function here integrationCC where...
...each rank calls approximate area and finds their sum...
...using MPI_Reduce, find the final answer on Rank 0*/



double approximateArea(int commSz, int myRank){
	const double a=0.0, b=200.0;
	const int  quads = 50000000;
	double const width = (b-a) / (double) quads;
	
	//TASK: Include quadsPerRank, the number of quads to be for each rank
	int quadsPerRank = quads / commSz;

	//This is to be used for the final rank, if the number of quads is not divisible by the number of ranks
	if(myRank == commSz - 1){	
		quadsPerRank = quads - (quadsPerRank*myRank);
	}
	
	//TASK: Modify startIter and endIter to change based on the rank and quadsPerRank.
	int startIter = myRank * quadsPerRank;
	int endIter = startIter + quadsPerRank;
	
	double x, y;
    
	double sum;
	int i;
    for(i = startIter; i <= endIter; i++){
		x = a + i * width;
		y = x + width;
		sum += 0.5*width*(func(x) + func(y));
    }
	return sum;
}

double func(double x){
	return pow(x,1.5)/3.1 - x/log(3.0);
}

