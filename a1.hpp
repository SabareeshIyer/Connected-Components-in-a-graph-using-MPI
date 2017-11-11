/*  Sabareesh
 *  Meenakshisundaram Balasubramanian
 *  smeenaks
 */

#ifndef A1_HPP
#define A1_HPP

#include <vector>
#include <algorithm>
#include <fstream>
#include <mpi.h>

using namespace std;


int connected_components(std::vector<signed char>& A, int n, int q, const char* out, MPI_Comm comm) {
    
    int size = n/q;
    std::vector<int> P (size*size, 0);
    int rank;
    MPI_Comm_rank(comm, &rank);

    for(int i = 0; i < size; i++){
    	for(int j = 0; j < size; j++){
    		if(A[i*size + j] == 1)
    			P[i*size + j] = i+(rank/q)*size;
    	}
    }

    std::vector<int> max (size, -1);
    for(int i = 0; i < size; i++){
    	for(int j = 0; j < size; j++){
    		if(P[i*size + j] > max[j])
    			max[j] = P[i*size + j];
    	}
    }

    int col = rank % q;
    int row = rank / q;

    MPI_Comm col_comm;
    MPI_Comm row_comm;
    MPI_Comm_split(comm, col, rank, &col_comm);
    MPI_Comm_split(comm, row, rank, &row_comm);

    std::vector<int> recv(size, 0);
    MPI_Allreduce(max.data(), recv.data(), size, MPI_INT, MPI_MAX, col_comm);

    int new_rank;
    MPI_Comm_rank(col_comm, &new_rank);

    //int print_rank = 8;
    //print recv[] (after first Allreduce - col wise)
    /*
    if(rank == print_rank){
    	for(int i = 0; i< size; i++){
	    	cout<<"i is "<<i<<" and recv[i] is "<<recv[i]<<endl;
    	}
    }
	*/

	for(int i=0; i<size; i++){
		for(int j=0; j<size; j++){
			P[i*size + j] = recv[j];
		}
	}

	//print P
    /*
    if(rank == print_rank){
    	for(int i = 0; i< size; i++){
	    	for(int j = 0; j<size; j++){
	    		cout<<P[i*size +j]<<" ";
	    	}
	    	cout<<endl;
    	}
    }
    */
    

    std::vector<int> M (size*size, 0);
    for(int i = 0; i < size; i++){
    	for(int j = 0; j < size; j++){
    		if(A[i*size + j] == 1)
    			M[i*size + j] = P[i*size + j];
    	}
    }

    //print M
    /*
    if(rank == print_rank){
    	for(int i = 0; i< size; i++){
	    	for(int j = 0; j<size; j++){
	    		cout<<M[i*size +j]<<" ";
	    	}
	    	cout<<endl;
    	}
    }
    */
    
    //row-wise max, then all reduce
    std::vector<int> max_r(size, -1);
    for(int i = 0; i < size; i++){
    	for(int j = 0; j < size; j++){
    		if(M[i*size + j] > max_r[i])
    			max_r[i] = M[i*size + j];
    	}
    }

    MPI_Allreduce(max_r.data(), recv.data(), size, MPI_INT, MPI_MAX, row_comm);

    //print recv[] (after second Allreduce - row wise)
    /*
    if(rank == print_rank){
    	for(int i = 0; i< size; i++){
	    	cout<<"i is "<<i<<" and recv[i] is "<<recv[i]<<endl;
    	}
    }
    */

    std::vector<int> Q(size*size);
    for(int i=0; i<size; i++){
    	for(int j=0; j<size; j++){
    		Q[i*size + j] = recv[i];
    	}
    }

    //print Q
    /*
    if(rank == print_rank){
    	for(int i = 0; i< size; i++){
	    	for(int j = 0; j<size; j++){
	    		cout<<Q[i*size +j]<<" ";
	    	}
	    	cout<<endl;
    	}
    }
    */

	for(int i=0; i< size; i++){
		for(int j = 0; j< size; j++){
			if(Q[i*size + j] == j + (rank%q)*size)
				M[i*size + j] = P[i*size + j];
		}
	}

	//print M
	/*
    print_rank = 5;
    if(rank == print_rank){
    	for(int i = 0; i< size; i++){
	    	for(int j = 0; j<size; j++){
	    		cout<<M[i*size +j]<<" ";
	    	}
	    	cout<<endl;
    	}
    }
    */
    
	//row-wise max, then all reduce
	std::vector<int> P_Dash_Feeder(size, 0);
    std::fill(max_r.begin(), max_r.end(), -1);
    for(int i = 0; i < size; i++){
    	for(int j = 0; j < size; j++){
    		if(M[i*size + j] > max_r[i])
    			max_r[i] = M[i*size + j];
    	}
    }

    MPI_Allreduce(max_r.data(), P_Dash_Feeder.data(), size, MPI_INT, MPI_MAX, row_comm);

    std::vector<int> P_Dash(size*size, 0);
    for(int i = 0; i < size; i++){
    	for(int j = 0; j < size; j++){
    		P_Dash[i*size + j] = P_Dash_Feeder[i];
    	}
    }
    //row-wise ends

    /***************** Tree Hanging ***********************/

    for(int i = 0; i < size; i++){
    	for(int j = 0; j < size; j++){
    		if(P[i*size + j] == i+(rank/q)*size)
    			M[i*size + j] = P_Dash[i*size + j];
    	}
    }

    //row-wise max, then all reduce
    std::fill(max_r.begin(), max_r.end(), -1);
    for(int i = 0; i < size; i++){
    	for(int j = 0; j < size; j++){
    		if(M[i*size + j] > max_r[i])
    			max_r[i] = M[i*size + j];
    	}
    }

    MPI_Allreduce(max_r.data(), recv.data(), size, MPI_INT, MPI_MAX, row_comm);

     for(int i = 0; i < size; i++){
    	for(int j = 0; j < size; j++){
    		Q[i*size + j] = recv[i];
    	}
    }
    //row-wise ends
    
    for(int i = 0; i < size; i++){
    	for(int j = 0; j < size; j++){
    		int a = P_Dash[i*size + j];
    		int b = Q[i*size + j];
    		P[i*size + j] = a > b ? a : b;
    	}
    }

   	std::vector<int> retvec(size);
   	for(int i=0; i<size; i++)
   		retvec[i] = (P[i*size]);

   	std::vector<int> res(q*q*size);
   	MPI_Gather(retvec.data(), size, MPI_INT, res.data(), size, MPI_INT, 0, comm);
    MPI_Barrier(comm);

    if(rank == 0){
    	std::sort(res.begin(), res.end());
    	std::vector<int>::iterator it;
  		it = std::unique(res.begin(), res.end());
  		res.resize(std::distance(res.begin(), it));
	}	    

    return res.size();
} // connected_components

#endif // A1_HPP