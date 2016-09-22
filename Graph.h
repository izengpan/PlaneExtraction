#ifndef _PE_GRAPH_
#define _PE_GRAPH_

#define NODE_SIZE 4
#define ALPHA 0.015
#define SIGMA 0.0000016
#define EPSILON 5.5

#include <vector>
#include <eigen3/Eigen/Dense>

using namespace Eigen;

class Node {
	Matrix<float,Dynamic,3> pixels; //(row:(node*node) column:3)
	std::vector<Node*> edges;
	std::vector<float> extraData;
	Vector4f param;

	Node(float* arr, int length) {
		if (!(length % 3)) return;
		for (int i = 0; i < length; ++i) {
			pixels << arr[i];
		}
	}

	bool rejectNode();

};

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
	//use Eigen for matrix operations ;OpenBLAS may be faster
	Vector3f mean(0.0, 0.0, 0.0);
	for (int i = 0; i < pixels.rows; i++) {
		mean += pixels.row(i);
	}
	for (int i = 0; i < pixels.rows; i++) {
		pixels.row(i) -= mean;
	} 
		//use PCA to fit the plane
	Matrix3f PCAMatrix = pixels.transpose()*pixels;
	JacobiSVD<MatrixXf> svd(PCAMatrix, ComputeThinU | ComputeThinV);
	Vector3f eigenvalues = svd.singularValues();
	Matrix3f eigenvectors = svd.matrixU();
	
	int minIndex=0;
	for (int i = 0; i < eigenvalues.rows(); ++i)
		if (eigenvalues(i) < eigenvalues(minIndex)) i = minIndex;
	Vector3f normalVector;
	switch (minIndex) {
		case 0:
			normalVector = eigenvectors.col(1)*eigenvectors.col(2);
			break;
		case 1:
			normalVector = eigenvectors.col(0)*eigenvectors.col(2);
			break;
		case 2:
			normalVector = eigenvectors.col(0)*eigenvectors.col(1);
			break;
	}

	float MSE;
	for (int i = 0; i < pixels.rows; i++) {
		MSE += pow(normalVector.dot(pixels.row(i)),2);
	}
	if (MSE > pow(((SIGMA*mean(2)*mean(2) + EPSILON)), 2))
		return true;
	
	for (int i = 0; i < pixels.rows; i++) {
		pixels.row(i) += mean;
	}


	return false;
}


class Graph {
	int row, column;
	std::vector<Node*> graph;
};



#endif // !_PE_GRAPH_
