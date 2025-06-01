import numpy as np
import cv2
import os
from adjImg import get_adj
from flowImg import parse_flow_dict, max_flow

edge_n = 0
edge = c = 0

# Escreve as arestas e pesos do frame
def writeGraph(adjMat, vals, count, pixel_counts):
    global edge_n
    indicator_path = "./teste_graph_indicator.txt"
    edges_path = "./teste_A.txt"
    edge_attributes_path = "./teste_edge_attributes.txt"

    with open(indicator_path, "a") as indicator, \
         open(edges_path, "a") as edges, \
         open(edge_attributes_path, "a") as edge_attributes:

        for i, x in enumerate(adjMat):
            indicator.write(str(count) + "\n")
            for j, y in enumerate(x):
                if y > 0:
                    if pixel_counts[i] < pixel_counts[j]:
                        edges.write(f'{i + 1 + edge_n}, {j + 1 + edge_n}\n')
                        edge_attributes.write(str(0) + "\n")
    edge_n += len(vals)

# Escreve a classe
def writeClass(label: str):
    graph_labels_path = "./teste_graph_labels.txt"
    with open(graph_labels_path, "a") as graph_labels:
        graph_labels.write(label + "\n")

# Escreve o rótulo -> saída do fluxo
def getNodeLabels(string):
    node_labels_path = "./teste_node_labels.txt"
    with open(node_labels_path, "a") as node_labels:
        node_labels.write(string)
        node_labels.write("\n")

# Escreve as arestas e pesos do par de frames -> relação de regiões
def relation(last_node, current_node, time):
    edges_path = "./teste_A.txt"
    edge_attributes_path = "./teste_edge_attributes.txt"
    with open(edges_path, "a") as edges, \
         open(edge_attributes_path, "a") as edge_attributes:
        edges.write(f'{last_node}, {current_node}\n')
        edge_attributes.write(str(time) + "\n")

# Calcula a quantidade de pixels da região
def writeNodeInfo(img, node_info_path):
    height = len(img)
    width = len(img[0]) if height > 0 else 0
    labels = [[0 for _ in range(width)] for _ in range(height)]
    current_label = 1
    region_counts = []

    directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    #directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1), (-2, 0), (2, 0), (0, -2), (0, 2), (-2, -1), (-2, 1), (2, -1), (2, 1), (-1, -2), (-1, 2), (1, -2), (1, 2), (-2, -2), (-2, 2), (2, -2), (2, 2)]
        
    for y in range(height):
        for x in range(width):
            if labels[y][x] == 0:
                stack = [(y, x)]
                labels[y][x] = current_label
                color = img[y][x]
                count = 1

                while stack:
                    cy, cx = stack.pop()
                    for dy, dx in directions:
                        ny, nx = cy + dy, cx + dx
                        if 0 <= ny < height and 0 <= nx < width:
                            if labels[ny][nx] == 0 and img[ny][nx] == color:
                                labels[ny][nx] = current_label
                                stack.append((ny, nx))
                                count += 1
                region_counts.append(count)
                current_label += 1

    with open(node_info_path, "a") as f:
        for count in region_counts:
            f.write(f"{count}\n")

    return region_counts

# Main
def writeForOneFile(fileDir: str, count: int, name, label):
    global edge, c
    node_info_path = "./teste_node_info.txt"

    if count == 1:
        open(node_info_path, "w").close()

    videos = []

    for fileName in os.listdir(f"./videoClassification/dataset/{name}/{fileDir}/frameResult/00020"):
        if "ppm" in fileName:
            img = cv2.imread(f'./videoClassification/dataset/{name}/{fileDir}/frameResult/00020/{fileName}', 0)
            videos.append(img)
            size = len(np.unique(img))
            vals1 = {}
            adjMat1 = np.zeros((size, size))
            adjMat1, vals1 = get_adj(img)
            pixel_counts = writeNodeInfo(img, node_info_path)

            writeGraph(adjMat1, vals1, count, pixel_counts)
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
