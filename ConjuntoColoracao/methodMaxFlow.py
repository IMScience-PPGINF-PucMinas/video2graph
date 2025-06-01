import numpy as np
import cv2
import os
import timeit
from adjImg import get_adj
from flowImg import parse_flow_dict, max_flow
import glob

edge_n = 0

# Escreve as arestas e pesos do frame
def writeGraph(adjMat, vals, count):
    global edge_n
    indicator_path = "./teste_graph_indicator.txt"
    edges_path = "./teste_A.txt"
    edge_attributes_path = "./teste_edge_attributes.txt"
    indicator = open(indicator_path, "a")
    edges = open(edges_path, "a")
    edge_attributes = open(edge_attributes_path, "a")
    for i, x in enumerate(adjMat):
        indicator.write(str(count)+ "\n")
        for j, y in enumerate(x):
            if(y > 0):
                edges.write(f'{i + 1 + edge_n}, {j + 1 + edge_n}\n')
                edge_attributes.write(str(0)+ "\n")
    edge_n+= len(vals)
    indicator.close()
    edges.close()
    edge_attributes.close()

# Escreve a classe
def writeClass(label:str):
    graph_labels_path = "./teste_graph_labels.txt"
    graph_labels = open(graph_labels_path, "a")
    graph_labels.write(label + "\n")
    graph_labels.close()

# Escreve o rótulo -> saída do fluxo
def getNodeLabels(string):
    node_labels_path = "./teste_node_labels.txt"
    node_labels = open(node_labels_path, "a")
    node_labels.write(string)
    node_labels.write("\n")
    node_labels.close()

# Escreve as arestas e pesos do par de frames -> relação de regiões
def relation(last_node, current_node, time):
    edges_path = "./teste_A.txt"
    edge_attributes_path = "./teste_edge_attributes.txt"
    edges = open(edges_path, "a")
    edge_attributes = open(edge_attributes_path, "a")
    edges.write(f'{last_node}, {current_node}\n')
    edge_attributes.write(str(time)+ "\n")
    edges.close()
    edge_attributes.close()

edge = c = 0

# Main
def writeForOneFile(fileDir: str, count: int, name, label):
    global edge, c

    videos = []

    for fileName in os.listdir("./videoClassification/dataset/" + name + "/" + fileDir + "/frame"):
        if "png" in fileName:
            img = cv2.imread('./videoClassification/dataset/' + name + '/' + fileDir + '/frame/' + fileName, 0)
            videos.append(img)
            size = len(np.unique(img))
            vals1 = {}
            adjMat1 = np.zeros((size, size))
            adjMat1, vals1 = get_adj(img)

            writeGraph(adjMat1, vals1, count)
            edge += len(vals1)

    writeClass(label)

    for i in range(len(videos)):
        if (i + 1) == len(videos):
            ultimo = len(vals1)
            for j in range(ultimo):
                #s = f"{i}, 0"
                s = f"0"
                getNodeLabels(s)
        else:
            img = videos[i]
            next_img = videos[i + 1]
            flow_value, flow_dict = max_flow(img, next_img)
            correspondences = parse_flow_dict(flow_dict)
            t = len(correspondences)

            for (source, target, flow) in correspondences:
                if target:
                    #s = f"{i}, 0, {i+1}, {int(flow * 100)}"
                    s = f"{int(flow * 100)}"
                    getNodeLabels(s)

                    src_region = source[1]
                    tgt_region = target[1]
                    t = len(correspondences)
                    relation((src_region + c + 1), (tgt_region + c + t + 1), (i + 1))
                else:
                    #s = f"{i}, 0"
                    s = f"0"
                    getNodeLabels(s)
            c += t

    c = edge

count = 1
for name in os.listdir("./videoClassification/dataset"):
    for folderName in os.listdir("./videoClassification/dataset/%s" %name):
        if os.path.isdir("./videoClassification/dataset/"+ name + "/" + folderName):
            if name == "walking":
                writeForOneFile(folderName, count, name, "1")
                count+=1
            if name == "handwaving":
                writeForOneFile(folderName, count, name, "5")
                count+=1

