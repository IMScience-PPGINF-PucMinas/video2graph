from grakel.graph import Graph
from sklearn.utils import Bunch
import numpy as np
import sys
import networkx as nx
import time
import os

def read_data(
        name,
        with_classes=True,
        as_graphs=False,
        is_symmetric=False):
    indicator_path = f"./datasets/{name}/{name}/raw/{name}_graph_indicator.txt"
    edges_path = f"./datasets/{name}/{name}/raw/{name}_A.txt"
    node_labels_path = f"./datasets/{name}/{name}/raw/{name}_node_labels.txt"
    edge_attributes_path = f"./datasets/{name}/{name}/raw/{name}_edge_attributes.txt"
    graph_classes_path = f"./datasets/{name}/{name}/raw/{name}_graph_labels.txt"

    # node graph correspondence
    ngc = dict()
    # edge line correspondence
    elc = dict()
    # dictionary that keeps sets of edges
    Graphs = dict()
    # dictionary of labels for nodes
    node_labels = dict()
    # dictionary of labels for edges
    edge_labels = dict()

    # Associate graphs nodes with indexes
    with open(indicator_path, "r") as f:
        for (i, line) in enumerate(f, 1):
            ngc[i] = int(line.strip())
            if int(line.strip()) not in Graphs:
                Graphs[int(line.strip())] = set()
            if int(line.strip()) not in node_labels:
                node_labels[int(line.strip())] = dict()
            if int(line.strip()) not in edge_labels:
                edge_labels[int(line.strip())] = dict()

    # Extract graph edges
    with open(edges_path, "r") as f:
        for (i, line) in enumerate(f, 1):
            edge = list(map(int, line.strip().split(",")))
            elc[i] = (edge[0], edge[1])
            Graphs[ngc[edge[0]]].add((edge[0], edge[1]))
            if is_symmetric:
                Graphs[ngc[edge[1]]].add((edge[1], edge[0]))

    # Extract node labels
    with open(node_labels_path, "r") as f:
        for (i, line) in enumerate(f, 1):
            labels = list(map(int, line.strip().split(",")))
            if len(labels) > 3:
                node_labels[ngc[i]][i] = (labels[2], labels[3])
            else:
                node_labels[ngc[i]][i] = (0, 0)

    # Extract edge attributes
    with open(edge_attributes_path, "r") as f:
        for (i, line) in enumerate(f, 1):
            label = float(line.strip())
            edge_labels[ngc[elc[i][0]]][elc[i]] = label
            if is_symmetric:
                edge_labels[ngc[elc[i][1]]][(elc[i][1], elc[i][0])] = label

    Gs = []
    if as_graphs:
        for i in range(1, len(Graphs)+1):
            Gs.append(Graph(Graphs[i], node_labels[i], edge_labels[i]))
    else:
        for i in range(1, len(Graphs)+1):
            Gs.append([Graphs[i], node_labels[i], edge_labels[i]])

    if with_classes:
        classes = []
        with open(graph_classes_path, "r") as f:
            for line in f:
                classes.append(int(line.strip()))

        classes = np.array(classes, dtype=np.int32)
        return Bunch(data=Gs, target=classes)
    else:
        return Bunch(data=Gs)

# A função de redução transitiva permanece a mesma, já que os grafos são agora DAGs
for dataset in os.listdir("./datasets"):
    edges_path = f"./output/{dataset}/{dataset}_A.txt"
    edge_attributes_path = f"./output/{dataset}/{dataset}_edge_attributes.txt"

    # Leitura dos dados
    data_att = read_data(dataset)
    data = data_att.data
    target = data_att.target

    # Garantir que os diretórios de saída existam
    os.makedirs(os.path.dirname(edges_path), exist_ok=True)

    with open(edges_path, "w") as edges_out, open(edge_attributes_path, "w") as edge_attrs_out:
        execution_time = 0
        for graph in data:
            # graph[0]: set of edges (tuplas)
            # graph[1]: node_labels (não utilizado aqui)
            # graph[2]: edge_labels (atributos das arestas)

            G = nx.DiGraph()
            G.add_edges_from(graph[0])

            # Como as arestas já são escritas apenas de menor para maior pixel_count, o grafo deve ser um DAG
            # Aplicar redução transitiva diretamente
            start_time = time.time()
            TR = nx.transitive_reduction(G)
            end_time = time.time()
            execution_time += (end_time - start_time) * 1000  # Em milissegundos

            # Ordenar as arestas de forma determinística (por nó de origem, depois nó de destino)
            edges_sorted = sorted(TR.edges(), key=lambda edge: (edge[0], edge[1]))
            for edge in edges_sorted:
                edges_out.write(f"{edge[0]}, {edge[1]}\n")
                edge_attrs_out.write(f"{graph[2].get(edge, 0)}\n")  # Usa 0 se o atributo não existir

        print(f"Execution time of {dataset}: {execution_time:.2f} ms")
