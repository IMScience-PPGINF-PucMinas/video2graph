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
            stack.extend([(cx + dx, cy + dy) for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]])

    return region, original_value

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

def calculate_region_SAD(region, img1, img2, dx, dy):
    sad = 0
    for (x, y) in region:
        nx, ny = x + dx, y + dy
        if 0 <= nx < len(img2) and 0 <= ny < len(img2[0]):
            sad += abs(int(img1[x][y]) - int(img2[nx][ny]))
        else:
            sad += abs(int(img1[x][y]))
    return sad

def build_flow_graph(img1, img2, search_range=1):
    regions_img1 = identify_regions(img1)
    regions_img2 = identify_regions(img2)
    G = nx.DiGraph()

    source = 'source'
    sink = 'sink'
    G.add_node(source)
    G.add_node(sink)

    used_img2_regions = set()

    for i, (region, value) in enumerate(regions_img1):
        G.add_node(('img1', i, value))
        best_match = None
        best_sad = float('inf')
        for j, (region2, value2) in enumerate(regions_img2):
            if value == value2 and j not in used_img2_regions:
                for dx in range(-search_range, search_range + 1):
                    for dy in range(-search_range, search_range + 1):
                        SAD = calculate_region_SAD(region, img1, img2, dx, dy)
                        if SAD < best_sad:
                            best_sad = SAD
                            best_match = j
        if best_match is not None:
            capacity = 1 / (1 + best_sad)
            G.add_edge(('img1', i, value), ('img2', best_match, value), capacity=capacity)
            used_img2_regions.add(best_match)
        else:
            G.add_edge(('img1', i, value), sink, capacity=len(region))

    for i, (region, value) in enumerate(regions_img1):
        G.add_edge(source, ('img1', i, value), capacity=len(region))

    for j, (region, value) in enumerate(regions_img2):
        G.add_edge(('img2', j, value), sink, capacity=len(region))

    return G, source, sink

def max_flow(img1, img2, search_range=1):
    G, source, sink = build_flow_graph(img1, img2, search_range)
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
                    correspondences.append((source_node, None, 0))
    return correspondences

'''
img1 = cv2.imread('./caminho/para/imagem1.ppm', 0)
img2 = cv2.imread('./caminho/para/imagem2.ppm', 0)

flow_value, flow_dict = max_flow(img1, img2)
correspondences = parse_flow_dict(flow_dict)
print("Flow value:", flow_value)
print("Correspondences:", correspondences)
'''