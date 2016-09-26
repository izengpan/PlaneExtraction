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



class PlaneParam {
public:
	Vector3f normal;
	float d;
	PlaneParam(Vector3f normal = Vector3f(0.0, 0.0, 0.0), float d = 0.0) :normal(normal), d(d) {}
};

class FittingParam {
public:
	Vector3f normal;
	float meanSquareError;
	FittingParam(Vector3f normal, float meanSqaureError) :normal(normal), meanSquareError(meanSqaureError) {}
};

FittingParam fitPlane(Matrix3f& const covarianceMatrix) {
	JacobiSVD<MatrixXf> svd(covarianceMatrix, ComputeThinU | ComputeThinV);
	Vector3f eigenvalues = svd.singularValues();
	Matrix3f eigenvectors = svd.matrixU();

	Vector3f normal;

	int minIndex = 0;
	for (int i = 0; i < eigenvalues.rows(); ++i)
		if (eigenvalues(i) < eigenvalues(minIndex)) i = minIndex;
	switch (minIndex) {
	case 0:
		normal = eigenvectors.col(1).cross(eigenvectors.col(2));
		break;
	case 1:
		normal = eigenvectors.col(0).cross(eigenvectors.col(2));
		break;
	case 2:
		normal = eigenvectors.col(0).cross(eigenvectors.col(1));
		break;
	}

	return FittingParam(normal, eigenvalues(minIndex));
}




class Node {
public:
	Matrix<float,Dynamic, 3, RowMajor> pixels; //(row:(node*node) column:3)
	
	std::set<Node**> edges;
	PlaneParam planeParam;
	float meanSquareError;

	struct {
		Matrix3f covarianceMatrix;
		Vector3f meanVector = Vector3f(0.0, 0.0, 0.0);
		int pixelNumber = 0;
	}extraData;


	float meanSquareError = 0;

	Node(int* arr=NULL, int x=0,int y=0,int nodesize=0,int row=0,int column=0) {
		pixels.resize(nodesize*nodesize, 3);
		if (arr == NULL) return;
		for (int i = 0; i < nodesize; ++i)
			for (int j = 0; j < nodesize; ++j) {
				pixels.row(i) << arr[x*nodesize*column + i*column + y*nodesize + j],arr[x*nodesize*column + i*column + y*nodesize + j + 1],	arr[x*nodesize*column + i*column + y*nodesize + j + 2];
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

	float merge(Node* node) {

		extraData.meanVector = extraData.meanVector*extraData.pixelNumber + node->extraData.meanVector*node->extraData.pixelNumber;
		extraData.pixelNumber += node->extraData.pixelNumber;
		extraData.meanVector /= extraData.pixelNumber;
		extraData.covarianceMatrix += node->extraData.covarianceMatrix;
	}

	void plane(Matrix3f& const covarianceMatrix) {
		FittingParam param = fitPlane(covarianceMatrix);
		meanSquareError = param.meanSquareError;
		planeParam = PlaneParam(param.normal,param.normal.dot(extraData.meanVector));
	}

	bool rejectNode();

	bool operator>(const Node node) {
		return (meanSquareError < node.meanSquareError);
	}
};

bool Node::rejectNode() {


	//missing data
	for (int i = 0; i < pixels.rows(); ++i) {
		if (pixels.row(i)(2) == 0)
			return true;
	}

	//depth discountinuity
	for (int i = 0; i < NODE_SIZE; ++i) 
		for (int j = 0; j < NODE_SIZE; ++j) {
			if (i != 0)
				if (abs(pixels(i*NODE_SIZE + j, 2) - pixels((i - 1)*NODE_SIZE + j, 2)) > 2 * ALPHA*(abs(pixels(i*NODE_SIZE + j, 2)) + 0.5)) return true;
			if (j != 0)
				if (abs(pixels(i*NODE_SIZE + j, 2) - pixels(i*NODE_SIZE + j - 1, 2)) > 2 * ALPHA*(abs(pixels(i*NODE_SIZE + j, 2)) + 0.5)) return true;
			if (i != NODE_SIZE - 1)
				if (abs(pixels(i*NODE_SIZE + j, 2) - pixels((i + 1)*NODE_SIZE + j, 2)) > 2 * ALPHA*(abs(pixels(i*NODE_SIZE + j, 2)) + 0.5)) return true;
			if (j != NODE_SIZE - 1)
				if (abs(pixels(i*NODE_SIZE + j, 2) - pixels(i *NODE_SIZE + j + 1, 2)) > 2 * ALPHA*(abs(pixels(i*NODE_SIZE + j, 2)) + 0.5)) return true;
		}

	//over-MSE
	plane(extraData.covarianceMatrix);
	
	if (meanSquareError > pow(((SIGMA*extraData.meanVector(2)*extraData.meanVector(2) + EPSILON)), 2))
		return true;
	
	pixels.resize(0, 0);

	return false;
}

class Graph {
public:
	int row, column;
	std::vector<Node*> graph;

	Graph(int row, int column) :row(row), column(column) {}

	void connectEdge();
	void addNode(Node* node) {
		graph.push_back(node);
	}
};

//¶þ¼¶Ö¸Õë

void Graph::connectEdge() {
	for (int i = 0; i < row; i++)
		for (int j = 0; j < column; j++) {
			if (graph[i*column + j] == NULL)
				continue;
			if (i != 0 && i != row-1) {
				if (graph[(i - 1)*column + j] != NULL&&graph[(i + 1)*column + j] != NULL)
					if (graph[(i - 1)*column + j]->planeParam.normal.dot(graph[(i + 1)*column + j]->planeParam.normal) > T_ANG) {
						graph[i*column + j]->edges.insert(&graph[(i - 1)*column + j]);
						graph[i*column + j]->edges.insert(&graph[(i + 1)*column + j]);
						graph[(i - 1)*column + j]->edges.insert(&graph[i*column + j]);
						graph[(i + 1)*column + j]->edges.insert(&graph[i*column + j]);
					}
			}
			if (j != 0 && j != column - 1) {
				if (graph[i*column + j - 1] != NULL&&graph[i*column + j + 1] != NULL)
					if (graph[i*column + j - 1]->planeParam.normal.dot(graph[i*column + (j + 1)]->planeParam.normal) > T_ANG)
					{
						graph[i*column + j]->edges.insert(&graph[i*column + j - 1]);
						graph[i*column + j]->edges.insert(&graph[i*column + j + 1]);
						graph[i*column + j - 1]->edges.insert(&graph[i*column + j]);
						graph[i*column + j + 1]->edges.insert(&graph[i*column + j]);
					}
			}
		}
}


#endif // !_PE_GRAPH_
