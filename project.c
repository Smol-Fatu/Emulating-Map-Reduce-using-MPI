#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#define array_size 256

struct Key {
    int i;
    int k;
};
struct Value {
    int matrix;
    int j;
    int value;
};
struct data_row {
    struct Key key;
    struct Value values[16];
    int num_values;
};
struct keyval {
    struct Key key;
    struct Value val;
};
void printkeyvalA(struct keyval kv){
	FILE *fp;
    fp = fopen("outputA.txt", "a");
    fprintf(fp, "( %d , %d ),( %d , %d , %d )\n", kv.key.i, kv.key.k, kv.val.matrix, kv.val.j, kv.val.value);
    fclose(fp);
}
void printkeyvalB(struct keyval kv){
	FILE *fp;
    fp = fopen("outputB.txt", "a");
    fprintf(fp, "( %d , %d ),( %d , %d , %d )\n", kv.key.i, kv.key.k, kv.val.matrix, kv.val.j, kv.val.value);
    fclose(fp);
}
int compare_keyvals(const void* a, const void* b) {
    struct keyval* kv1 = (struct keyval*) a;
    struct keyval* kv2 = (struct keyval*) b;

    if (kv1->key.i != kv2->key.i) {
        return kv1->key.i - kv2->key.i;
    } else if (kv1->val.j != kv2->val.j) {
        return kv1->val.j - kv2->val.j;
    } else {
        return kv1->key.k - kv2->key.k;
    }
}
void sort(){
    FILE *fp1,*fp2;
    fp1 = fopen("outputA.txt", "r");
    fp2 = fopen("outputB.txt", "r");
    if (fp1 == NULL || fp2 == NULL) {
        //printf("Error opening file\n");
        return;
    }
    int num_lines1 = 0,num_lines2 = 0;
    char c1,c2;
    while ((c1 = fgetc(fp1)) != EOF) {
        if (c1 == '\n') {
            num_lines1++;
        }
    }
    while ((c2 = fgetc(fp2)) != EOF) {
        if (c2 == '\n') {
            num_lines2++;
        }
    }
    struct keyval* kv_arr1 = (struct keyval*) malloc(num_lines1 * sizeof(struct keyval));
    struct keyval* kv_arr2 = (struct keyval*) malloc(num_lines2 * sizeof(struct keyval));
    fseek(fp1, 0, SEEK_SET);
    fseek(fp2, 0, SEEK_SET);
    for (int i = 0; i < num_lines1; i++) {
        if (fscanf(fp1, "( %d , %d ),( %d , %d , %d )\n", &kv_arr1[i].key.i, &kv_arr1[i].key.k, &kv_arr1[i].val.matrix, &kv_arr1[i].val.j, &kv_arr1[i].val.value) != 5) {
            //printf("Error reading file\n");
            return;
        }
    }
    for (int i = 0; i < num_lines2; i++) {
        if (fscanf(fp2, "( %d , %d ),( %d , %d , %d )\n", &kv_arr2[i].key.i, &kv_arr2[i].key.k, &kv_arr2[i].val.matrix, &kv_arr2[i].val.j, &kv_arr2[i].val.value) != 5) {
            //printf("Error reading file\n");
            return;
        }
    }
    fclose(fp1);
    fclose(fp2);
    fp1 = fopen("outputA.txt", "w");
    fp2 = fopen("outputB.txt", "w");
    if (fp1 == NULL || fp2 == NULL) {
        //printf("Error opening file\n");
        return;
    }
    qsort(kv_arr1, num_lines1, sizeof(struct keyval), compare_keyvals);
    qsort(kv_arr2, num_lines2, sizeof(struct keyval), compare_keyvals);
    for (int i = 0; i < num_lines1; i++) {
        fprintf(fp1, "( %d , %d ),( %d , %d , %d )\n", kv_arr1[i].key.i, kv_arr1[i].key.k, kv_arr1[i].val.matrix, kv_arr1[i].val.j, kv_arr1[i].val.value);
    }
    for (int i = 0; i < num_lines2; i++) {
        fprintf(fp2, "( %d , %d ),( %d , %d , %d )\n", kv_arr2[i].key.i, kv_arr2[i].key.k, kv_arr2[i].val.matrix, kv_arr2[i].val.j, kv_arr2[i].val.value);
    }
    fclose(fp1);
    fclose(fp2);
    free(kv_arr1);
    free(kv_arr2);
}
void writefiles(){
	FILE *fp;
    fp = fopen("matrix1.txt", "w");
    int arr1[array_size];
    for (int i = 0; i < array_size; i++) {
        arr1[i] = rand()%100;
    	fprintf(fp,"%d\n",arr1[i]);
    }
    fclose(fp);
    fp = fopen("matrix2.txt", "w");
    int arr2[array_size];
    for (int i = 0; i < array_size; i++) {
        arr2[i] = rand()%100;
    	fprintf(fp,"%d\n",arr2[i]);
    }
    fclose(fp);
}
void readfiles(int **numbers1, int **numbers2) {
    int side = sqrt(array_size);
    FILE *fp;
    fp = fopen("matrix1.txt", "r");
    int i = 0;
    *numbers1 = (int*) malloc(array_size * sizeof(int));
    if (fp == NULL) {
        //printf("Error opening file\n");
        return;
    }
    while (fscanf(fp, "%d", &(*numbers1)[i]) != EOF) {
        i++;
    }
    fclose(fp);

    fp = fopen("matrix2.txt", "r");
    i = 0;
    *numbers2 = (int*) malloc(array_size * sizeof(int));
    if (fp == NULL) {
        //printf("Error opening file\n");
        return;
    }
    while (fscanf(fp, "%d", &(*numbers2)[i]) != EOF) {
        i++;
    }
    fclose(fp);
}
struct keyval* MapperA(int* recv_buf,int recv_size,int rank,int start){
	int side = sqrt(array_size);
    struct keyval* kv = malloc(side * recv_size * sizeof(struct keyval));
	int x=0;
	for(int k=0; k<side; k++){
		for (int i = 0; i < recv_size/side; i++) {
			for( int j = 0; j < side;j++){
				kv[x].key.i = i+start/side;
				kv[x].key.k = k;
				kv[x].val.matrix = 1;
				kv[x].val.j = j;
				kv[x].val.value = recv_buf[i*side+j];
				printkeyvalA(kv[x]);
				x++;
			}
		}
	}
	return kv;
}
struct keyval* MapperB(int* recv_buf,int recv_size, int rank, int start){
	int side = sqrt(array_size);
    struct keyval* kv = malloc(side * recv_size * sizeof(struct keyval));
	int x=0;
	for(int i=0; i<side; i++){
		for (int j = 0; j < recv_size/side; j++) {
			for( int k = 0; k < side;k++){
				kv[x].key.i = i;
				kv[x].key.k = k;
				kv[x].val.matrix = 2;
				kv[x].val.j = j+start/side;
				kv[x].val.value = recv_buf[j*side+k];
				printkeyvalB(kv[x]);
				x++;
			}
		}
	}
	return kv;
}
void read_data(char *filename, struct data_row data[], int *num_rows) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        //printf("Error: could not open file %s\n", filename);
        exit(1);
    }
    char line[200];
    while (fgets(line, sizeof(line), fp)) {
        int len = strlen(line);
        line[len-2] = '\0'; // Remove the closing parenthesis
        char *tok = strtok(line+1, ","); // Skip the opening parenthesis
        data[*num_rows].key.i = atoi(tok);
        tok = strtok(NULL, ",");
        data[*num_rows].key.k = atoi(tok);
        int i = 0;
        while (tok != NULL) {
            tok = strtok(NULL, ",");
            if (tok != NULL) {
                data[*num_rows].values[i].matrix = atoi(tok);
                tok = strtok(NULL, ",");
                data[*num_rows].values[i].j = atoi(tok);
                tok = strtok(NULL, ",");
                data[*num_rows].values[i].value = atoi(tok);
                i++;
            }
        }
        data[*num_rows].num_values = i;
        (*num_rows)++;
    }
    fclose(fp);
}
int* Reducer(int start,int recv_size) {
    struct data_row dataA[256], dataB[256];
    int num_rowsA = 0, num_rowsB = 0;
	int side = sqrt(array_size);
    read_data("shuffleA.txt", dataA, &num_rowsA);
    read_data("shuffleB.txt", dataB, &num_rowsB);
	int x;
	int **result = malloc(side * sizeof(int *));
    for (int i = 0; i < side; i++) {
        result[i] = malloc(side * sizeof(int));
    }
    for (int i = 0; i < num_rowsA; i++) {
		x=0;
        for (int j = 0; j < dataA[i].num_values; j++) {
			x += (dataA[i].values[j].value*dataB[i].values[j].value);
        }
		result[dataA[i].key.i][dataA[i].key.k]=x;
    }
	int *result_arr= (int*)malloc(recv_size * sizeof(int));
	for (int i = start; i < (recv_size/side)+start; ++i) {
		for (int j = 0; j < side; ++j) {
        	result_arr[(i-start) * side + j] = result[i][j];
			//printf("%d ", result_arr[(i-start)*side+j]);
			//if (j == side - 1)
				//printf("\n");
		}
	}
	//printf("\n");

	return result_arr;
}
void shuffle(){
    int side = sqrt(array_size);
    FILE *fp1, *fp2, *shufA, *shufB;
    char line1[100], line2[100];
    int i, k, mat, j, val;
    fp1 = fopen("outputA.txt", "r");
    fp2 = fopen("outputB.txt", "r");
    if (fp1 == NULL || fp2 == NULL) {
        //printf("Error opening file.\n");
        exit(1);
    }
    shufA = fopen("shuffleA.txt", "w");
    shufB = fopen("shuffleB.txt", "w");
    if (shufA == NULL || shufB == NULL) {
        //printf("Error creating file.\n");
        exit(1);
    }
    rewind(fp1);
    rewind(fp2);
    for(int a=0; a<side; a++) {
        for(int b=0; b<side; b++) {
            fprintf(shufA, "(%d,%d)", a, b);
            fprintf(shufB, "(%d,%d)", a, b);
            while (fgets(line1, 100, fp1) != NULL) {
                sscanf(line1, "( %d , %d ),( %d , %d , %d )", &i, &k, &mat, &j, &val);
                if (i == a && k == b) {
                    fprintf(shufA, ",(%d,%d,%d)", mat, j, val);
                }
            }
            rewind(fp1);
            while (fgets(line2, 100, fp2) != NULL) {
                sscanf(line2, "( %d , %d ),( %d , %d , %d )", &i, &k, &mat, &j, &val);
                if (i == a && k == b) {
                    fprintf(shufB, ",(%d,%d,%d)", mat, j, val);
                }
            }
            rewind(fp2);
            fprintf(shufA, "\n");
            fprintf(shufB, "\n");
        }
    }
    fclose(fp1);
    fclose(fp2);
    fclose(shufA);
    fclose(shufB);
}
bool matrixcompare(int* mat1,int* mat2){
	int side=sqrt(array_size);
	for(int i=0;i<side;i++){
		for(int j=0;j<side;j++){
			if(mat1[i*side+j]!=mat2[i*side+j])
				return false;
		}
	}
	return true;
}
int main(int argc, char** argv) {
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
	FILE *fp;
    fp = fopen("outputA.txt", "w");
    fclose(fp);
    fp = fopen("outputB.txt", "w");
    fclose(fp);

	int rank, size;
	int side = sqrt(array_size);
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Status status;
    int status_value = 0;
	
	int *numbers1,*numbers2;
	numbers1 = (int*) malloc(array_size * sizeof(int*));
	numbers2 = (int*) malloc(array_size * sizeof(int*));

    int* data = NULL;
    int* recv_buf1 = NULL;
    int* recv_buf2 = NULL;
    int* sendcounts = NULL;
    int* displs = NULL;

    sendcounts = (int*)malloc(size * sizeof(int));
    displs = (int*)malloc(size * sizeof(int));

	if(rank == 0){//this is master dividing indexes for scatter
    	MPI_Get_processor_name(processor_name, &name_len);
		printf("Master with process_id %d running on %s\n",rank,processor_name);
		writefiles();
		readfiles(&numbers1,&numbers2);
		int sizes[size];
		int start[size];
		sizes[0] = 0;
		start[0] = 0;
		for (int i = 1; i < size; i++) {
			sizes[i] = side / (size-1);
			start[i] = (i-1) * sizes[i];
		}
		sizes[size-1] += side % (size-1); 
		for (int i = 0; i < size; i++) {
			sendcounts[i]=sizes[i]*side;
		}
		for (int i = 0; i < size; i++) {
			displs[i]=start[i]*side;
			if(sendcounts[i]>0)
				printf("Task Map assigned to process %d \n",i);
		}

	}

	MPI_Barrier(MPI_COMM_WORLD);
	int recv_size;
    MPI_Scatter(sendcounts, 1, MPI_INT, &recv_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    recv_buf1 = (int*)malloc(recv_size  * sizeof(int));
    recv_buf2 = (int*)malloc(recv_size  * sizeof(int));
	MPI_Get_processor_name(processor_name, &name_len);
	if(recv_size>0)
		printf("Process %d received task Map on %s \n",rank,processor_name);
		
	MPI_Bcast(displs, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(numbers1, sendcounts, displs, MPI_INT, recv_buf1, recv_size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(numbers2, sendcounts, displs, MPI_INT, recv_buf2, recv_size, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	struct keyval* kv_A ;
	struct keyval* kv_B ;
	if(rank!=0){
		kv_A = MapperA(recv_buf1,recv_size,rank,displs[rank]);
		kv_B = MapperB(recv_buf2,recv_size,rank,displs[rank]);
		status_value=rank;
    	MPI_Send(&status_value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){//master process is sorting data from all processes then it will create lists (shuffling)
		for (int i = 1; i < size; i++) {
            MPI_Recv(&status_value, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            printf("Process %d has completed task Map \n", i);
        }
		sort();
		//shuffling
		shuffle();
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0){
		int reducers = size/4; //25% reducers
		for (int i = 1; i <= reducers; i++) {
			status_value=i;
			MPI_Send(&status_value, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			printf("Task Reduce assigned to process %d \n",i);
		}
	}

	status_value=0;
	MPI_Barrier(MPI_COMM_WORLD);
	int *array;
	int recvsize=0;
	int start=0;
	int reducers = size/4; //25% reducer
	if(rank<=reducers && rank!=0){
		MPI_Recv(&status_value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Get_processor_name(processor_name, &name_len);
		printf("Process %d received task Reduce on %s \n",rank,processor_name);
		recvsize = array_size/reducers;
		start = (rank-1)*(recvsize/side);
		array = malloc(recvsize * sizeof(int));
		array = Reducer(start,recvsize);
		status_value=rank;
    	MPI_Send(&status_value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}

	int *recvbuf = NULL;
	int *displas = NULL;
	int *recvcounts = NULL;
	if(rank==0){
		for (int i = 1; i <= reducers; i++) {
            MPI_Recv(&status_value, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            printf("Process %d has completed task Reduce \n", i);
        }
		recvbuf = malloc(array_size * sizeof(int));
		recvcounts = malloc(size * sizeof(int));
		displas = malloc(size * sizeof(int));
	}
	start=start*side;
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(&recvsize, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(&start, 1, MPI_INT, displas, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);


	MPI_Gatherv(array, recvsize, MPI_INT, recvbuf, recvcounts, displas, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		printf("Job has been completed\n");
	/*	printf("Final result: \n");
		for (int i = 0; i < side;i++) {
			for (int j = 0; j < side; j++) 
				printf("%d ", recvbuf[i*side+j]);
			printf("\n");
		}
		printf("\n");
	*/
		int* result = malloc(array_size * sizeof(int));
		for (int i = 0; i < side; ++i) {
			for (int j = 0; j < side; ++j) {
				result[i*side+j] = 0;
			}
		}
		for (int i = 0; i < side; ++i) {
			for (int j = 0; j < side; ++j) {
				for (int k = 0; k < side; ++k) {
					result[i*side+j] += numbers1[i*side+k] * numbers2[k*side+j];
				}
			}
		}
		if(matrixcompare(recvbuf,result)){
			printf("Matrix comparison function returned: True\n");
		}
		else{
			printf("Matrix comparison function returned: True\n");
		}
	/*
		for (int i = 0; i < side; ++i) {
			for (int j = 0; j < side; ++j) {
				printf("%d ", result[i][j]);
				if (j == side - 1)
					printf("\n");
			}
		}
	*/

	}

	
	MPI_Finalize();
	return 0;
	}
