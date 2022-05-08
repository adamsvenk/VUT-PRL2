/**
 * @file pro.cpp
 * @author Adam Švenk (xsvenk00@stud.fit.vutbr.cz)
 * @brief Binary-tree preorder traversal
 * @version 1.0
 * @date 2022-05-05
 */

#include <stdio.h>
#include <mpi.h>

#define MASTER_RANK 0

// Structure represeing node of the binary tree
typedef struct Node {
    int id = -1;
    int left = -1;
    int right = -1;
} node;

// Structure represeing edge of the binary tree
typedef struct Edge {
    int id = -1;
    int start = -1;
    int end = -1;
} edge;

// Structure represeing element of adjacency list
typedef struct Neighbour {
    int id = -1;
    int number = 0;
    int list[6] = {-1};
} neighbour;

// Assign edge to the neighbour structure list
void assignNeighbour(Neighbour *neighbour, Edge edge) {
    neighbour->list[neighbour->number++] = edge.id;
}

// Get largest integer from the array of integers
int largestIntFromArray(int array[], int length) {
    int largest = 0;
    
    for (int i = 0; i < length; i++) {
        if (largest < array[i]) {
                largest = array[i];
        }
    }

    return largest;
}

// Get reverse edge
Edge *getReverseEdge(Edge edges[], Edge edge, int edgeNumber) {
    
    for (int i = 0; i < edgeNumber; i++) {
        if (edges[i].start == edge.end && edges[i].end == edge.start) {
            return &edges[i];
        }
    }

    return NULL;
}

// Check if the edge is forwarding
int isForwardingEdge(Edge edge) {
    if (edge.start < edge.end) {
        return 1;
    }

    return 0;
}

// Main program body
int main(int argc, char *argv[]){
    MPI_Status status;
    int rank, size, data;

    // Initialize MPI values
    if (MPI_Init(&argc, &argv)) {
        fprintf(stderr, "Error during MPI initialization\n");
    }

    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank)) {
        fprintf(stderr, "Error during MPI initialization\n");
    }

    if (MPI_Comm_size(MPI_COMM_WORLD, &size)) {
        fprintf(stderr, "Error during MPI initialization\n");
    }

    // Read input string of nodes from program argument
    char *nodeString = argv[1];

    // Calculate number of nodes
    int nodeNumber = strlen(nodeString);
    // Calculate number of edges based od 2*n-2 evaluation
    int edgeNumber = 2 * nodeNumber - 2;
    
    // Array containing all nodes
    Node nodes[nodeNumber];
    // Array containing all edges
    Edge edges[edgeNumber];
    // Array containing neighbours of all nodes
    Neighbour neighbours[nodeNumber];

    // Init edges IDs
    for (int i = 0; i < edgeNumber; i++) {
        edges[i].id = i;
    }

    // Init nodes IDs (correct ID of character)
    // Init neighbours IDs
    for (int i = 0; i < nodeNumber; i++) {
        nodes[i].id = i + 1;
        neighbours[i].id = i + 1;
    }
    
    // Init nodes with their childrens
    for(int i = 0; i < nodeNumber / 2; i++) {
        nodes[i].left = 2 * (i + 1);
        nodes[i].right = 2 * (i + 1) + 1;
    }

    // Init all edges with their corresponding starting and ending nodes
    for (int i = 0; i < nodeNumber; i++) {
        if (nodes[i].left != - 1) {
            edges[(i + 1) * 4 - 4].start = nodes[i].id;
            edges[(i + 1) * 4 - 3].start = nodes[i].left;
            edges[(i + 1) * 4 - 4].end = nodes[i].left;
            edges[(i + 1) * 4 - 3].end = nodes[i].id;
        }

        if (nodes[i].right != -1) {
            edges[(i + 1) * 4 - 2].start = nodes[i].id;
            edges[(i + 1) * 4 - 1].start = nodes[i].right;
            edges[(i + 1) * 4 - 2].end = nodes[i].right;
            edges[(i + 1) * 4 - 1].end = nodes[i].id;
        }    
    }

    // Get neighbours of all nodes
    for (int i = 0; i < nodeNumber; i++) {
        for (int j = 0; j < edgeNumber; j++) {
            if (edges[j].start == i+1) {
                assignNeighbour(&neighbours[i], edges[j]);
                assignNeighbour(&neighbours[i], *getReverseEdge(edges, edges[j], edgeNumber));
            }
        }
    }

    Edge myEdge = edges[rank];
    Edge *myEdgeReverse = getReverseEdge(edges, myEdge, edgeNumber);
    int toPutInEtour = -1;

    // Create element for Eulerian Path
    // Parallel algorithm from the lecture
    for (int i = 0; i < nodeNumber; i++) {
        for (int j = 0; j < neighbours[i].number / 2; j++) {
            if (neighbours[i].list[j*2] == myEdgeReverse->id) {
                if ((j + 1) * 2 < neighbours[i].number) {
                    toPutInEtour = neighbours[i].list[(j + 1) * 2];
                } else {
                    toPutInEtour = neighbours[myEdge.end - 1].list[0];
                }
                break;
            }
        }
    }

    // Eulerian Path
    Edge eulerianPath[edgeNumber];

    // Receive Eulerian elements and build Eulerian Path
    if (rank == 0) {
        Edge eTour[edgeNumber];
        eTour[0] = edges[toPutInEtour];
        
        // Receive values from other ranks
        for (int i = 1; i < edgeNumber; i++) {
            MPI_Recv(&data, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
            
            if (status.MPI_ERROR != 0) {
                fprintf(stderr, "Error during MPI communication\n");
            }

            eTour[i] = edges[data];
        }

        int index = 0;
        eulerianPath[0] = edges[index];

        for (int i = 0; i < edgeNumber - 1; i++) {
            eulerianPath[i + 1] = eTour[index];
            index = eTour[index].id;
        }

    } else {
        // Send value to master rank
        MPI_Send(&toPutInEtour, 1, MPI_INT, MASTER_RANK, rank, MPI_COMM_WORLD);
    }

    // Sent Eulerian Path to all processes
    if (rank == 0) {
        for (int i = 1; i < edgeNumber; i++) {
            for (int j = 0; j < edgeNumber; j++) {
                MPI_Send(&eulerianPath[j].id, 1, MPI_INT, i, j, MPI_COMM_WORLD);
            }
        }
    } else {
        int edgeIndex = 0;
        for (int i = 0; i < edgeNumber; i++) {
            MPI_Recv(&edgeIndex, 1, MPI_INT, MASTER_RANK, i, MPI_COMM_WORLD, &status);

            if (status.MPI_ERROR != 0) {
                fprintf(stderr, "Error during MPI communication\n");
            }

            eulerianPath[i] = edges[edgeIndex];
        }
    }

    // Get weight of edge
    // Do parallel
    int myWeight = isForwardingEdge(myEdge);

    // Array of all weights
    int weights[edgeNumber];

    // Receive weights from all processes
    if (rank == 0) {
        weights[0] = myWeight;

        for (int i = 1; i < edgeNumber; i++) {
            MPI_Recv(&data, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);

            if (status.MPI_ERROR != 0) {
                fprintf(stderr, "Error during MPI communication\n");
            }

            weights[i] = data;
        }
    } else {
        MPI_Send(&myWeight, 1, MPI_INT, MASTER_RANK, rank, MPI_COMM_WORLD);
    }

    int suffixSums[edgeNumber] = {0};

    // Calculate Suffix Sums
    if (rank == 0) {
        int eulerianPathIndex[edgeNumber];

        for (int i = 0; i < edgeNumber; i++) {
            for (int j = 0; j < edgeNumber; j++) {
                if (i == eulerianPath[j].id) {
                    eulerianPathIndex[i] = j;
                    break;
                }
            }
        }

        int maxSum = 0;

        for (int i = edgeNumber - 1; i >= 0; i--) {
            for (int j = 0; j < edgeNumber; j++) {
                if (eulerianPathIndex[j] == i) {
                    maxSum += weights[j];
                    suffixSums[j] = maxSum;
                }
            }
        }
    }

    // Send Suffix Sums to all processes
    if (rank == 0) {
        for (int i = 1; i < edgeNumber; i++) {
            for (int j = 0; j < edgeNumber; j++) {
                MPI_Send(&suffixSums[j], 1, MPI_INT, i, j, MPI_COMM_WORLD);
            }
        }
    } else {
        int suffixReceived;
        
        for (int i = 0; i < edgeNumber; i++) {
            MPI_Recv(&suffixReceived, 1, MPI_INT, MASTER_RANK, i, MPI_COMM_WORLD, &status);

            if (status.MPI_ERROR != 0) {
                fprintf(stderr, "Error during MPI communication\n");
            }

            suffixSums[i] = suffixReceived;
        }
    }

    int preOrderValue = - 1;

    // If the edge is forwarding, do a correction
    // Parallel algorithm from the lecture
    if (rank != 0 && myWeight == 1) {
        preOrderValue = largestIntFromArray(suffixSums, edgeNumber) - suffixSums[myEdge.id] + 1;
    }

    int preOrder[nodeNumber];

    // Receive corrected data from processes
    if (rank == 0) {
        for (int i = 0; i < nodeNumber; i++) {
            preOrder[i] = i;
        }

        for (int i = 1; i < edgeNumber; i++) {
            data = - 1;
            
            MPI_Recv(&data, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);

            if (status.MPI_ERROR != 0) {
                fprintf(stderr, "Error during MPI communication\n");
            }

            if (data != - 1) {
                preOrder[data] = edges[i].end - 1;
            }
        }

        // Print final binary-tree preorder traversal
        for (int i = 0; i < nodeNumber; i++) {
            printf("%c", nodeString[preOrder[i]]);
        }
        
        printf("\n");
    } else {
        MPI_Send(&preOrderValue, 1, MPI_INT, MASTER_RANK, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}