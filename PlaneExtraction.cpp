#define T_NUM 100

#include <queue>
#include "Graph.h"
#include <fstream>
#include <iostream>
//#include <cv.hpp>
//#include <highgui.hpp>


//global variable
int row;
int column;
int* pixelNode;  //record which node does a pixel belong to
std::vector<Node*> coarseResult;

class cmp {
public:
	bool operator()(Node* n1, Node* n2) {
		return (*n1).meanSquareError > (*n2).meanSquareError;
	}
};


void AHCluster(std::vector<Node*>& graph) {
	std::priority_queue < Node*, std::vector<Node*>,cmp> queue;
	for (int i = 0; i < graph.size();++i) {
		if (graph[i] != NULL)
			queue.push(graph[i]);
	}
	while (!queue.empty())
	{
		Node* node = queue.top();
		queue.pop();
		if (node->ignore) continue;
		float minMSE=-1.0; 
		float z;
		Node* bestNode=NULL;
		Node* mergeNode=NULL;
		for (auto dest:node->edges) {
			float currentMSE;
			mergeNode = node;
			currentMSE = node->testMerge(dest);
			z = (node->extraData.meanVector[2] * node->extraData.pixelNumber + dest->extraData.meanVector[2] * dest->extraData.pixelNumber) / (node->extraData.pixelNumber + dest->extraData.pixelNumber);
			if ((minMSE < 0) || (currentMSE < minMSE)) {
				minMSE = currentMSE;
				bestNode = dest;
			}
		}
		if ((bestNode==NULL)||(minMSE >= pow(((SIGMA * z*z + EPSILON)), 2))) {
			if (node->extraData.pixelNumber > T_NUM) {
				coarseResult.push_back(node);
			}
			for (auto neighbouringNode : node->edges) {
				neighbouringNode->edges.erase(node);
			}
			node->ignore = true;
		}
		else {
			node->merge(bestNode);
			for (auto neighbouringNode : bestNode->edges) {
				neighbouringNode->edges.erase(bestNode);
				if (neighbouringNode != node) {
					neighbouringNode->edges.insert(node);
					node->edges.insert(neighbouringNode);
				}
			}
			bestNode->ignore = true;
			queue.push(node);
		}
	}

}

void refine(std::vector<Node*>& graph) {
	for (int i = 0; i < graph.size(); ++i) {
		for (auto index : graph[i]->pixelIndex) {
			pixelNode[index[0] * column + index[1]] = i;
		}
	}
}


int main() {
	std::ifstream os;
	os.open("depth.txt");
	float* data;
	float* testdata;
	os >> row>>column;
	data = new float[row*column*3];
	pixelNode = new int[row*column];
	int k = 0;
	for (int i = 0; i < row; ++i) 
		for (int j = 0; j < column; ++j) {
			int d;
			os >> d;
			data[k++] = (float)(i - row / 2) /528 * d;
			data[k++] = (float)(j - column / 2) / 528 * d;
			data[k++] = d;
		}

	std::vector<Node*> g;
	Graph rawdata(g);

	for (int i = 0; i < row / NODE_SIZE; ++i) 
		for (int j = 0; j < column / NODE_SIZE; ++j) {
			Node* currentNode = new Node(data, i, j, NODE_SIZE, row, column);
			if (currentNode->rejectNode()) {
				rawdata.addNode(NULL);
				continue;
			}
			else rawdata.addNode(currentNode);
		
		}
	
	rawdata.connectEdge(row/NODE_SIZE,column/NODE_SIZE);
	AHCluster(rawdata.graph);
	k = 0;
	/*cv::Mat segment(480,640,CV_32FC1);
	for (auto node : coarseResult) {
		k++;
		if (node->pixelIndex.size() > 800) {
			for (auto index : node->pixelIndex) {
				segment.at<float>(index[0], index[1]) = (float)k / 150;
			}
		}
	}
	cv::imshow("...", segment);
	cvWaitKey();*/
	return 0;
}