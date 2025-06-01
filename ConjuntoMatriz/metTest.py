import numpy as np
from adjImg import get_adj
from flowFixo import parse_flow_dict, max_flow

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
def writeForOneFile(videos, count, label):
    global edge, c

    for img in videos:
        size = len(np.unique(img))
        vals1 = {}
        adjMat1=np.zeros((size, size))
        adjMat1, vals1 = get_adj(img)

        writeGraph(adjMat1, vals1, count)
        edge+= len(vals1)

    writeClass(label)

    for i in range(len(videos)):
        if((i+1) == len(videos)):
            ultimo = len(vals1)
            for j in range(ultimo):
                #s = f"{i}, 0"
                s = f"0"
                getNodeLabels(s)
        else :
            img = videos[i]
            next_img = videos[i + 1]
            flow_value, flow_dict = max_flow(img, next_img)
            correspondences = parse_flow_dict(flow_dict)

            for (source, target, flow) in correspondences:
                if(target) :
                    #s = f"{i}, 0, {i+1}, {int(flow * 100)}"
                    s = f"{int(flow * 100)}"
                    getNodeLabels(s)

                    src_region = source[1]
                    tgt_region = target[1]
                    t = len(correspondences)
                    relation((src_region+c+1), (tgt_region+c+t+1), (i+1))
                    #print((src_region+c+1), (tgt_region+c+t+1), (i+1))
                else : 
                    #s = f"{i}, 0"
                    s = f"0"
                    getNodeLabels(s)
            c += t
        
    c = edge

count = 1

# Objetos passando pela cena - 1

img0 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 4],
    [1, 1, 1, 1, 4],
    [2, 2, 2, 2, 4],
    [2, 2, 2, 2, 2],
]

img1 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 4, 1],
    [1, 1, 1, 4, 1],
    [2, 2, 2, 4, 2],
    [2, 2, 2, 2, 2],
]

img2 = [
    [1, 1, 1, 1, 1],
    [1, 3, 4, 1, 1],
    [1, 1, 4, 1, 1],
    [2, 2, 4, 2, 2],
    [2, 2, 2, 2, 2],
]

img3 = [
    [1, 1, 1, 1, 1],
    [1, 4, 1, 1, 1],
    [1, 4, 1, 1, 1],
    [2, 4, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img4 = [
    [1, 1, 1, 1, 1],
    [4, 3, 1, 1, 1],
    [4, 1, 1, 1, 1],
    [4, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

# Objetos passando pela cena - 2

img02 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 4],
    [1, 1, 1, 1, 4],
    [2, 2, 2, 2, 4],
    [2, 2, 2, 2, 2],
]

img12 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 4, 4],
    [1, 1, 1, 4, 4],
    [2, 2, 2, 4, 4],
    [2, 2, 2, 2, 2],
]

img22 = [
    [1, 1, 1, 1, 1],
    [1, 3, 4, 4, 1],
    [1, 1, 4, 4, 1],
    [2, 2, 4, 4, 2],
    [2, 2, 2, 2, 2],
]

img32 = [
    [1, 1, 1, 1, 1],
    [1, 4, 4, 1, 1],
    [1, 4, 4, 1, 1],
    [2, 4, 4, 2, 2],
    [2, 2, 2, 2, 2],
]

img42 = [
    [1, 1, 1, 1, 1],
    [4, 4, 1, 1, 1],
    [4, 4, 1, 1, 1],
    [4, 4, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

# Objetos passando pela cena - 3

img03 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 1, 1, 4],
    [2, 2, 2, 2, 4],
    [2, 2, 2, 2, 2],
]

img13 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 1, 4, 1],
    [2, 2, 2, 4, 2],
    [2, 2, 2, 2, 2],
]

img23 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 4, 1, 1],
    [2, 2, 4, 2, 2],
    [2, 2, 2, 2, 2],
]

img33 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 4, 1, 1, 1],
    [2, 4, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img43 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [4, 1, 1, 1, 1],
    [4, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

# Objetos passando pela cena - 4

img04 = [
    [1, 1, 1, 1, 1],
    [4, 3, 1, 1, 1],
    [4, 1, 1, 1, 1],
    [4, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img14 = [
    [1, 1, 1, 1, 1],
    [1, 4, 1, 1, 1],
    [1, 4, 1, 1, 1],
    [2, 4, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img24 = [
    [1, 1, 1, 1, 1],
    [1, 3, 4, 1, 1],
    [1, 1, 4, 1, 1],
    [2, 2, 4, 2, 2],
    [2, 2, 2, 2, 2],
]

img34 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 4, 1],
    [1, 1, 1, 4, 1],
    [2, 2, 2, 4, 2],
    [2, 2, 2, 2, 2],
]

img44 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 4],
    [1, 1, 1, 1, 4],
    [2, 2, 2, 2, 4],
    [2, 2, 2, 2, 2],
]

# Objetos passando pela cena - 5

img05 = [
    [1, 1, 1, 1, 1],
    [4, 1, 1, 1, 3],
    [4, 1, 1, 1, 1],
    [4, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img15 = [
    [1, 1, 1, 1, 1],
    [4, 1, 1, 3, 1],
    [4, 1, 1, 1, 1],
    [4, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img25 = [
    [1, 1, 1, 1, 1],
    [4, 1, 3, 1, 1],
    [4, 1, 1, 1, 1],
    [4, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img35 = [
    [1, 1, 1, 1, 1],
    [4, 3, 1, 1, 1],
    [4, 1, 1, 1, 1],
    [4, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img45 = [
    [1, 1, 1, 1, 1],
    [3, 1, 1, 1, 1],
    [4, 1, 1, 1, 1],
    [4, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

# Objeto pulando na cena

img5 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [2, 2, 5, 2, 2],
    [2, 2, 2, 2, 2],
]

img6 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 5, 1, 1],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img7 = [
    [1, 1, 1, 1, 1],
    [1, 3, 5, 1, 1],
    [1, 1, 1, 1, 1],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img8 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 5, 1, 1],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img9 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [2, 2, 5, 2, 2],
    [2, 2, 2, 2, 2],
]

# Objeto pulando na cena 2

img52 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [2, 5, 5, 2, 2],
    [2, 2, 2, 2, 2],
]

img62 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 5, 5, 1, 1],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img72 = [
    [1, 1, 1, 1, 1],
    [1, 5, 5, 1, 1],
    [1, 1, 1, 1, 1],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img82 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 5, 5, 1, 1],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img92 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [2, 5, 5, 2, 2],
    [2, 2, 2, 2, 2],
]


# Objeto pulando na cena 3

img53 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [2, 2, 5, 2, 2],
    [2, 2, 2, 2, 2],
]

img63 = [
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [2, 2, 5, 2, 2],
    [2, 2, 2, 2, 2],
]

img73 = [
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [2, 3, 5, 2, 2],
    [2, 2, 2, 2, 2],
]

img83 = [
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [2, 2, 5, 2, 2],
    [2, 2, 2, 2, 2],
]

img93 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [2, 2, 5, 2, 2],
    [2, 2, 2, 2, 2],
]

# Objeto pulando na cena 4

img54 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [2, 2, 2, 2, 5],
    [2, 2, 2, 2, 2],
]

img64 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 1, 1, 5],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img74 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 5],
    [1, 1, 1, 1, 1],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img84 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 1, 1, 5],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img94 = [
    [1, 1, 1, 1, 1],
    [1, 3, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [2, 2, 2, 2, 5],
    [2, 2, 2, 2, 2],
]


# Objeto pulando na cena 5

img55 = [
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 5, 5, 1],
    [2, 2, 5, 5, 2],
    [2, 2, 2, 2, 2],
]

img65 = [
    [1, 1, 1, 1, 1],
    [1, 1, 5, 5, 1],
    [1, 1, 5, 5, 1],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img75 = [
    [1, 1, 5, 5, 1],
    [1, 1, 5, 5, 1],
    [1, 1, 1, 1, 1],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img85 = [
    [1, 1, 1, 1, 1],
    [1, 1, 5, 5, 1],
    [1, 1, 5, 5, 1],
    [2, 2, 2, 2, 2],
    [2, 2, 2, 2, 2],
]

img95 = [
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 5, 5, 1],
    [2, 2, 5, 5, 2],
    [2, 2, 2, 2, 2],
]

video11 = [img0, img1, img2, img3, img4]
video12 = [img02, img12, img22, img32, img42]
video13 = [img03, img13, img23, img33, img43]
video14 = [img04, img14, img24, img34, img44]
video15 = [img05, img15, img25, img35, img45]

video21 = [img5, img6, img7, img8, img9]
video22 = [img52, img62, img72, img82, img92]
video23 = [img53, img63, img73, img83, img93]
video24 = [img54, img64, img74, img84, img94]
video25 = [img55, img65, img75, img85, img95]

allVideos = {
    'boxing': [video11, video12, video13, video14, video15],
    'walking': [video21, video22, video23, video24, video25],
}

for category, videos in allVideos.items():
    if category == "boxing":
        for video in videos:
            writeForOneFile(video, count, "0")
            count += 1
    if category == "walking":
        for video in videos:
            writeForOneFile(video, count, "1")
            count += 1