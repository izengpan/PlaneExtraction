#define T_NUM 800
#include <queue>
#include "Graph.h"


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
		Node** bestNode;
		Node* mergeNode;
		for (auto edge : (*node)->edges) {
			Node* currentNode = *edge;
			float currentMSE;
			mergeNode = *node;
			mergeNode->merge(currentNode);
			mergeNode->plane();
			if ((minMSE < 0) || (mergeNode->meanSquareError < minMSE)) {
				minMSE = mergeNode->meanSquareError;
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

int main() {

	return 0;
}