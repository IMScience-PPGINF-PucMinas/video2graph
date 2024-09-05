import numpy as np
import cv2
import os
import timeit
from adjImg import get_adj
from flowImg import parse_flow_dict, max_flow
import glob

edge_n = 0

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

def writeClass(label:str):
    graph_labels_path = "./teste_graph_labels.txt"
    graph_labels = open(graph_labels_path, "a")
    graph_labels.write(label + "\n")
    graph_labels.close()

def addAdj1Method(adjMat1, vals1, adjMat2, vals2):
    vals1 = [*vals1]
    vals2 = [*vals2]
    newVals = np.unique(vals1 + vals2).tolist()
    size = len(newVals)
    adjMat=np.zeros((size, size))
    for i, x in enumerate(adjMat):
        for j, y in enumerate(x):
            v1 = 0
            v2 = 0
            label_i = newVals[i]
            label_j = newVals[j]
            # pegar o numero na matriz1
            if label_i in vals1 and label_j in vals1:
                index_i = vals1.index(label_i)
                index_j = vals1.index(label_j)
                v1=adjMat1[index_i][index_j]
            if label_i in vals2 and label_j in vals2:
                index_i = vals2.index(label_i)
                index_j = vals2.index(label_j)
                v2=adjMat2[index_i][index_j]
            adjMat[i][j] = v1 + v2
    return adjMat, newVals

def addAdj2Method(adjMat1, vals1, adjMat2, vals2, countMat):
    vals1 = [*vals1]
    vals2 = [*vals2]
    newVals = np.unique(vals1 + vals2).tolist()
    size = len(newVals)
    adjMat=np.zeros((size, size))
    for i in range (0, size):
        for j in range(i + 1, size):
            v1 = 0
            v2 = 0
            label_i = newVals[i]
            label_j = newVals[j]
            # pegar o numero na matriz1
            if label_i in vals1 and label_j in vals1:
                index_i = vals1.index(label_i)
                index_j = vals1.index(label_j)
                v1=adjMat1[index_i][index_j]
            if label_i in vals2 and label_j in vals2:
                index_i = vals2.index(label_i)
                index_j = vals2.index(label_j)
                v2=adjMat2[index_i][index_j]
            if v2 > 0 and v2 > 0:
                adjMat[i][j] = countMat
                adjMat[j][i] = countMat
            else:
                adjMat[i][j] = v1 + v2
                adjMat[j][i] = v1 + v2
    return adjMat, newVals

def getNodeLabels(string):
    node_labels_path = "./teste_node_labels.txt"
    node_labels = open(node_labels_path, "a")
    node_labels.write(string)
    node_labels.write("\n")
    node_labels.close()

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

def writeForOneFile(fileDir: str, count: int, name, label):
    global edge, c

    videos = []

    for fileName in os.listdir("./videoClassification/datasets-avi/" + name + "/" + fileDir + "/frameResult/00050"):
        if "ppm" in fileName:
            img = cv2.imread('./videoClassification/datasets-avi/' + name + '/' + fileDir + '/frameResult/00050/' + fileName, 0)
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
                #s = f"{i}, 100"
                s = f"0"
                getNodeLabels(s)
        else:
            img = videos[i]
            next_img = videos[i + 1]
            flow_value, flow_dict = max_flow(img, next_img)
            correspondences = parse_flow_dict(flow_dict)

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
                    #s = f"{i}, 100"
                    s = f"0"
                    getNodeLabels(s)
            c += t

    c = edge

count = 1
for name in os.listdir("./videoClassification/datasets-avi"):
    for folderName in os.listdir("./videoClassification/datasets-avi/%s" %name):
        if os.path.isdir("./videoClassification/datasets-avi/"+ name + "/" + folderName):
            # if name == "boxing":
            #     writeForOneFile(folderName, count, name, "0")
            #     count+=1
            if name == "walking":
                writeForOneFile(folderName, count, name, "1")
                count+=1
            # if name == "jogging":
            #     writeForOneFile(folderName, count, name, "2")
            #     count+=1
            # if name == "handclapping":
            #     writeForOneFile(folderName, count, name, "3")
            #     count+=1
            # if name == "running":
            #     writeForOneFile(folderName, count, name, "4")
            #     count+=1
            if name == "handwaving":
                writeForOneFile(folderName, count, name, "5")
                count+=1
