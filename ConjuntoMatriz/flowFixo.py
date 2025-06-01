import networkx as nx

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
            stack.extend([(cx + dx, cy + dy) for dx, dy in 
                          [(-1, 0), (1, 0), (0, -1), (0, 1),
                           (-1, -1), (-1, 1), (1, -1), (1, 1)]])
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

def compute_centroid(region):
    x_sum = sum(x for (x, y) in region)
    y_sum = sum(y for (x, y) in region)
    n = len(region)
    return (x_sum / n, y_sum / n)

def get_weight_from_displacement(dx, dy):
    # Sem deslocamento:
    if dx == 0 and dy == 0:
        return 1

    # Decide qual variação é predominante
    if abs(dx) >= abs(dy):
        if dx < 0:
            return 0.8  # movimento para cima
        else:
            return 0.4  # movimento para baixo
    else:
        if dy > 0:
            return 0.6  # movimento para a direita
        else:
            return 0.2  # movimento para a esquerda

def build_flow_graph(img1, img2):
    regions_img1 = identify_regions(img1)
    regions_img2 = identify_regions(img2)
    G = nx.DiGraph()
    source = 'source'
    sink = 'sink'
    G.add_node(source)
    G.add_node(sink)
    
    # Pré-calcular os centroides das regiões de ambas as imagens
    centroids_img1 = [compute_centroid(region) for region, value in regions_img1]
    centroids_img2 = [compute_centroid(region) for region, value in regions_img2]
    
    used_img2_regions = set()
    
    # Para cada região em img1, procura a melhor correspondência em img2
    for i, (region, value) in enumerate(regions_img1):
        node_img1 = ('img1', i, value)
        G.add_node(node_img1)
        best_match = None
        best_weight = -1  # peso máximo encontrado até o momento
        centroid1 = centroids_img1[i]
        
        for j, (region2, value2) in enumerate(regions_img2):
            # Verifica se as regiões têm o mesmo valor e se a região candidata ainda não foi usada
            if value == value2 and j not in used_img2_regions:
                centroid2 = centroids_img2[j]
                dx = centroid2[0] - centroid1[0]
                dy = centroid2[1] - centroid1[1]
                weight = get_weight_from_displacement(dx, dy)
                # Seleciona a correspondência com maior peso (ou seja, a que apresenta menor deslocamento indesejado)
                if weight > best_weight:
                    best_weight = weight
                    best_match = j
                    
        if best_match is not None:
            node_img2 = ('img2', best_match, value)
            G.add_edge(node_img1, node_img2, capacity=best_weight)
            used_img2_regions.add(best_match)
        else:
            # Se não encontrar correspondência, conecta a região diretamente ao sink.
            G.add_edge(node_img1, sink, capacity=len(region))
    
    # Conecta o source às regiões de img1 e as regiões de img2 ao sink
    for i, (region, value) in enumerate(regions_img1):
        node_img1 = ('img1', i, value)
        G.add_edge(source, node_img1, capacity=len(region))
        
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
# Exemplo de imagens
img1 = [
    [1, 1, 1, 1, 1, 1],
    [1, 4, 1, 1, 1, 1],
    [1, 1, 2, 1, 1, 1],
    [1, 1, 2, 3, 1, 1],
    [1, 1, 2, 1, 3, 1],
    [1, 1, 1, 1, 1, 1],
]

img2 = [
    [1, 1, 1, 1, 1, 1],
    [1, 4, 1, 1, 1, 1],
    [1, 1, 2, 1, 1, 1],
    [1, 1, 2, 1, 1, 1],
    [1, 1, 2, 3, 3, 1],
    [1, 1, 1, 1, 1, 1],
]

flow_value, flow_dict = max_flow(img1, img2)
correspondences = parse_flow_dict(flow_dict)

for (source_node, target_node, flow) in correspondences:
    src_region = source_node[1]
    src_value = source_node[2]
    if target_node:
        tgt_region = target_node[1]
        tgt_value = target_node[2]
        print(f"Região na img1: {src_region}, Valor: {src_value} -> Região na img2: {tgt_region}, Valor: {tgt_value}, Fluxo (Peso): {flow}")
    else:
        print(f"Região na img1: {src_region}, Valor: {src_value} -> Sem correspondência, Fluxo (Peso): {flow}")
'''