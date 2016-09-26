#ifndef _PE_GRAPH_
#define _PE_GRAPH_

#define NODE_SIZE 4
#define ALPHA 0.015
#define SIGMA 0.0000016
#define EPSILON 5.5
#define T_ANG 0.966 //(cos(15))

#include <vector>
#include <set>
#include </Library/Eigen/include/eigen3/Eigen/Dense>

using namespace Eigen;

class param {
public:
	Vector3f normal;
	float d;
	param(Vector3f normal=Vector3f(0.0,0.0,0.0), float d=0.0) :normal(normal), d(d) {}
};



class Node {
public:
	Matrix<float, Dynamic, 3, RowMajor> pixels; //(row:(node*node) column:3)
	//Matrix<int, Dynamic, 2, RowMajor> pixelIndices;
	std::set<Node**> edges;
	param planeParam;

	struct {
		Matrix3f covarianceMatrix;
		Vector3f meanVector = Vector3f(0.0, 0.0, 0.0);
		int pixelNumber = 0;
	}extraData;


	float meanSquareError = 0;

	Node(float* arr, int length) {
		if (*arr == NULL || !(length % 3)) return;
		for (int i = 0; i < length; ++i) {
			pixels << arr[i];
			
		}
		initialization();
	}

	void initialization() {
		extraData.pixelNumber = pixels.rows();
		for (int i = 0; i < extraData.pixelNumber; i++) {
			extraData.meanVector += pixels.row(i);
		}
		extraData.meanVector /= pixels.rows();
		for (int i = 0; i <extraData.pixelNumber; i++) {
			pixels.row(i) -= extraData.meanVector;
		}
		extraData.covarianceMatrix = pixels.transpose()*pixels;
		for (int i = 0; i < extraData.pixelNumber; i++) {
			pixels.row(i) += extraData.meanVector;
		}
	}

	void merge(Node* node) {
		for (int i = 0; i < node->pixels.rows(); ++i)
			pixels << node->pixels.row(i);
		extraData.meanVector = extraData.meanVector*extraData.pixelNumber + node->extraData.meanVector*node->extraData.pixelNumber;
		extraData.pixelNumber += node->extraData.pixelNumber;
		extraData.meanVector /= extraData.pixelNumber;
		extraData.covarianceMatrix += node->extraData.covarianceMatrix;

	}

	void plane();

	bool rejectNode();

	bool operator>(const Node node) {
		return (meanSquareError < node.meanSquareError);
	}
};


//use Eigen for matrix operations ;OpenBLAS may be faster
void Node::plane() {
	JacobiSVD<MatrixXf> svd(extraData.covarianceMatrix, ComputeThinU | ComputeThinV);
	Vector3f eigenvalues = svd.singularValues();
	Matrix3f eigenvectors = svd.matrixU();

	int minIndex = 0;
	for (int i = 0; i < eigenvalues.rows(); ++i)
		if (eigenvalues(i) < eigenvalues(minIndex)) i = minIndex;
	switch (minIndex) {
	case 0:
		planeParam.normal = eigenvectors.col(1)*eigenvectors.col(2);
		break;
	case 1:
		planeParam.normal = eigenvectors.col(0)*eigenvectors.col(2);
		break;
	case 2:
		planeParam.normal = eigenvectors.col(0)*eigenvectors.col(1);
		break;
	}

	planeParam.d = planeParam.normal.dot(extraData.meanVector);
}

bool Node::rejectNode() {

	//depth discountinuity
	for (int i = 0; i < NODE_SIZE; ++i) 
		for (int j = 0; j < NODE_SIZE; ++j) {
			if (i != 0)
				if (abs(pixels(i*NODE_SIZE + j, 3) - pixels((i - 1)*NODE_SIZE + j, 3)) > 2 * ALPHA*(abs(pixels(i*NODE_SIZE + j, 3)) + 0.5)) return true;
			if (j != 0)
				if (abs(pixels(i*NODE_SIZE + j, 3) - pixels(i*NODE_SIZE + j - 1, 3)) > 2 * ALPHA*(abs(pixels(i*NODE_SIZE + j, 3)) + 0.5)) return true;
			if (i != NODE_SIZE - 1)
				if (abs(pixels(i*NODE_SIZE + j, 3) - pixels((i + 1)*NODE_SIZE + j, 3)) > 2 * ALPHA*(abs(pixels(i*NODE_SIZE + j, 3)) + 0.5)) return true;
			if (j != NODE_SIZE - 1)
				if (abs(pixels(i*NODE_SIZE + j, 3) - pixels(i *NODE_SIZE + j + 1, 3)) > 2 * ALPHA*(abs(pixels(i*NODE_SIZE + j, 3)) + 0.5)) return true;
		}

	//over-MSE
	plane();
	for (int i = 0; i < pixels.rows; i++) {
		meanSquareError += pow(planeParam.normal.dot(pixels.row(i)), 2);
	}
	if (meanSquareError > pow(((SIGMA*extraData.meanVector(2)*extraData.meanVector(2) + EPSILON)), 2))
		return true;

	return false;
}

class Graph {
public:
	int row, column;
	std::vector<Node*> graph;
	void connectEdge();
};

//¶þ¼¶Ö¸Õë

void Graph::connectEdge() {
	for (int i = 0; i < row; i++)
		for (int j = 0; j < column; j++) {
			if (graph[i*row + j] == NULL)
				continue;
			if (i != 0 && i != row-1) {
				if (graph[(i - 1)*row + j] != NULL&&graph[(i + 1)*row + j] != NULL)
					if (graph[(i - 1)*row + j]->planeParam.normal.dot(graph[(i + 1)*row + j]->planeParam.normal) > T_ANG) {
						graph[i*row + j]->edges.insert(&graph[(i - 1)*row + j]);
						graph[i*row + j]->edges.insert(&graph[(i + 1)*row + j]);
						graph[(i - 1)*row + j]->edges.insert(&graph[i*row + j]);
						graph[(i + 1)*row + j]->edges.insert(&graph[i*row + j]);
					}
			}
			if (j != 0 && j != column - 1) {
				if (graph[i*row + j - 1] != NULL&&graph[i*row + j + 1] != NULL)
					if (graph[i*row + j - 1]->planeParam.normal.dot(graph[i*row + (j + 1)]->planeParam.normal) > T_ANG)
					{
						graph[i*row + j]->edges.insert(&graph[i*row + j - 1]);
						graph[i*row + j]->edges.insert(&graph[i*row + j + 1]);
						graph[i*row + j - 1]->edges.insert(&graph[i*row + j]);
						graph[i*row + j + 1]->edges.insert(&graph[i*row + j]);
					}
			}
		}
}


#endif // !_PE_GRAPH_
