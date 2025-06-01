import networkx as nx
import cv2
import numpy as np

def flood_fill(matrix, x, y, visited):
    rows, cols = len(matrix), len(matrix[0])
    stack = [(x, y)]
    original_value = matrix[x][y]
    region = []

    while stack:
        cx, cy = stack.pop()
        if 0 <= cx < rows and 0 <= cy < cols and not visited[cx][cy] and matrix[cx][cy] == original_value:
            visited[cx][cy] = True
            region.append((cx, cy))
            # Utiliza adjacências de primeiro e segundo nível para segmentação
            stack.extend([
                (cx + dx, cy + dy) 
                for dx, dy in [
                    # Primeiro nível de adjacência
                    (-1, 0), (1, 0), (0, -1), (0, 1),
                    (-1, -1), (-1, 1), (1, -1), (1, 1)
                ]
            ])
    return region, original_value

'''
# Segundo nível de adjacência
(-2, 0), (2, 0), (0, -2), (0, 2),
(-2, -1), (-2, 1), (2, -1), (2, 1),
(-1, -2), (-1, 2), (1, -2), (1, 2),
(-2, -2), (-2, 2), (2, -2), (2, 2)
'''

def identify_regions(matrix):
    rows, cols = len(matrix), len(matrix[0])
    visited = [[False] * cols for _ in range(rows)]
    regions = []

    for i in range(rows):
        for j in range(cols):
            if not visited[i][j]:
                region, value = flood_fill(matrix, i, j, visited)
                if region:
                    regions.append((region, value))
    return regions

def compute_centroid(region):
    #Calcula o centroide (média das coordenadas) de uma região.
    x_sum = sum(x for (x, y) in region)
    y_sum = sum(y for (x, y) in region)
    n = len(region)
    return (x_sum / n, y_sum / n)

def get_weight_from_displacement(dx, dy):
    if dx == 0 and dy == 0:
        return 1.0
    # Se o movimento vertical é dominante:
    if abs(dx) >= abs(dy):
        return 0.8 if dx < 0 else 0.4
    else:  # movimento horizontal dominante
        return 0.6 if dy > 0 else 0.2

def build_flow_graph(img1, img2):
    regions_img1 = identify_regions(img1)
    regions_img2 = identify_regions(img2)
    G = nx.DiGraph()

    source = 'source'
    sink = 'sink'
    G.add_node(source)
    G.add_node(sink)

    # Pré-cálculo dos centroides das regiões
    centroids_img1 = [compute_centroid(region) for region, value in regions_img1]
    centroids_img2 = [compute_centroid(region) for region, value in regions_img2]

    used_img2_regions = set()

    for i, (region, value) in enumerate(regions_img1):
        node_img1 = ('img1', i, value)
        G.add_node(node_img1)
        best_match = None
        best_weight = -1.0
        centroid1 = centroids_img1[i]

        for j, (region2, value2) in enumerate(regions_img2):
            # Verifica se os valores (rótulos) são iguais e se a região candidata ainda não foi associada
            if value == value2 and j not in used_img2_regions:
                centroid2 = centroids_img2[j]
                dx = centroid2[0] - centroid1[0]
                dy = centroid2[1] - centroid1[1]
                weight = get_weight_from_displacement(dx, dy)
                if weight > best_weight:
                    best_weight = weight
                    best_match = j

        if best_match is not None:
            node_img2 = ('img2', best_match, value)
            # A capacidade da aresta será o peso calculado (entre 0.2 e 1.0)
            G.add_edge(node_img1, node_img2, capacity=best_weight)
            used_img2_regions.add(best_match)
        else:
            # Se não houver correspondência, conecta a região diretamente ao sink
            G.add_edge(node_img1, sink, capacity=len(region))

    # Conecta o source às regiões de img1 com capacidade igual ao tamanho da região
    for i, (region, value) in enumerate(regions_img1):
        node_img1 = ('img1', i, value)
        G.add_edge(source, node_img1, capacity=len(region))

    # Conecta as regiões de img2 ao sink com capacidade igual ao tamanho da região
    for j, (region, value) in enumerate(regions_img2):
        node_img2 = ('img2', j, value)
        G.add_edge(node_img2, sink, capacity=len(region))

    return G, source, sink

def max_flow(img1, img2):
    G, source, sink = build_flow_graph(img1, img2)
    flow_value, flow_dict = nx.maximum_flow(G, source, sink)
    return flow_value, flow_dict

def parse_flow_dict(flow_dict):
    correspondences = []
    for source_node, targets in flow_dict.items():
        if isinstance(source_node, tuple) and source_node[0] == 'img1':
            for target_node, flow_value in targets.items():
                if isinstance(target_node, tuple) and target_node[0] == 'img2' and flow_value > 0:
                    correspondences.append((source_node, target_node, flow_value))
                elif target_node == 'sink' and flow_value > 0:
                    correspondences.append((source_node, None, flow_value))
    return correspondences

'''
# Exemplo de uso com arquivos .ppm:

img1 = cv2.imread('./caminho/para/imagem1.ppm', 0)
img2 = cv2.imread('./caminho/para/imagem2.ppm', 0)

flow_value, flow_dict = max_flow(img1, img2)
correspondences = parse_flow_dict(flow_dict)

print("Flow value:", flow_value)
print("Correspondences:")
for (src, tgt, flow) in correspondences:
    src_region = src[1]
    src_value = src[2]
    if tgt:
        tgt_region = tgt[1]
        tgt_value = tgt[2]
        print(f"Região na img1: {src_region}, Valor: {src_value} -> Região na img2: {tgt_region}, Valor: {tgt_value}, Peso: {flow}")
    else:
        print(f"Região na img1: {src_region}, Valor: {src_value} -> Sem correspondência, Peso: {flow}")
'''
