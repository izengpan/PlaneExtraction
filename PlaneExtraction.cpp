#define T_NUM 800

#include <queue>
#include "Graph.h"
#include <fstream>

std::vector<Node*> coarseResult;



void AHCluster(Graph* graph) {
	std::priority_queue<Node**> queue;
	for (auto node : graph->graph)
		if (node != NULL) queue.push(&node);
	while (!queue.empty())
	{
		Node** node = queue.top();
		queue.pop();
		if ((*node) == NULL) continue;
		float minMSE=-1.0; 
		Node** bestNode=NULL;
		Node* mergeNode=NULL;
		for (auto edge : (*node)->edges) {
			Node* currentNode = *edge;
			float currentMSE;
			mergeNode = *node;
			currentMSE = (*node)->testMerge(currentNode);
			if ((minMSE < 0) || (currentMSE < minMSE)) {
				minMSE = currentMSE;
				bestNode = edge;
			}
		}
		if (minMSE >= pow(((SIGMA* mergeNode->extraData.meanVector(2)* mergeNode->extraData.meanVector(2) + EPSILON)), 2)) {
			if ((*node)->extraData.pixelNumber > T_NUM) {
				coarseResult.push_back(*node);
			}
			for (auto neighbouringNode : (*node)->edges) {
				(*neighbouringNode)->edges.erase(node);
			}
			*node = NULL;
		}
		else {
			(*node)->merge(*bestNode);
			for (auto neighbouringNode : (*bestNode)->edges) {
				(*neighbouringNode)->edges.erase(bestNode);
				if (neighbouringNode!=node)
					(*neighbouringNode)->edges.insert(node);
			}
			delete *bestNode;
			*bestNode = NULL;
			queue.push(node);
		}
	}

}

void refine(Graph* graph) {


}


int main() {
	std::ifstream os;
	os.open("depth.txt");
	int row;
	int column;
	int* data;
	os >> row>>column;
	data = new int[row*column*3];
	
	int k = 0;
	for (int i = 0; i < row; ++i) 
		for (int j = 0; j < column; ++j) {
			int d;
			os >> d;
			data[k++] = (float)(i - row / 2) /528 * d;
			data[k++] = (float)(j - column / 2) / 528 * d;
			data[k++] = d;
		}
	Graph rawdata(row / NODE_SIZE, column / NODE_SIZE);

	for (int i = 0; i < row / NODE_SIZE; ++i) 
		for (int j = 0; j < column / NODE_SIZE; ++j) {
			Node* currentNode = new Node(data, i, j, NODE_SIZE, row, column);
			if (currentNode->rejectNode()) rawdata.addNode(NULL);
			else rawdata.addNode(currentNode);
		
		}
	rawdata.connectEdge();

	return 0;
}