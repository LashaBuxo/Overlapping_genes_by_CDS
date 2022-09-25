class AnalyzerGraph:
    class GraphEdge:
        def __init__(self, gene1_id, gene2_id):
            self.node1 = gene1_id
            self.node2 = gene2_id

    def __init__(self):
        self.neighbors = {}

    def add_edge(self, edge: GraphEdge):
        self.__add_neighbor(edge.node1, edge.node2)
        self.__add_neighbor(edge.node2, edge.node1)

    def __add_neighbor(self, node, neighbor_data):
        if not self.neighbors.__contains__(node):
            self.neighbors[node] = []
        self.neighbors[node].append(neighbor_data)

    _visited = {}

    def dfs(self, node):
        self._visited[node] = True
        if not self.neighbors.__contains__(node): return []
        cluster_members = [node]
        for neighbor in self.neighbors[node]:
            if not self._visited[neighbor]:
                cluster_members += self.dfs(neighbor)
        return cluster_members

    def get_connected_clusters(self):
        nodes = self.neighbors.keys()

        for node in nodes:
            self._visited[node] = False

        overlapped_gene_clusters = []
        for node in nodes:
            if not self._visited[node]:
                cluster_members = self.dfs(node)
                # we only needs clusters with >1 members, coz we are finding overlapped genes
                if len(cluster_members) > 1:
                    overlapped_gene_clusters.append(cluster_members)

        return overlapped_gene_clusters