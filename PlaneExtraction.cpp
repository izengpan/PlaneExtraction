#include <queue>
#include "Graph.h"


void AHCluster(Graph* graph) {
	std::priority_queue<Node*> queue;
	for (auto node : graph->graph)
		if (node != NULL) queue.push(node);

}

int main() {

	return 0;
}