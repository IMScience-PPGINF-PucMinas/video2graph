import numpy as np

def get_adj(img):
    def dfs(img, visited, i, j, current_label):
        stack = [(i, j)]
        directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]
        #directions = [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1), (-2, 0), (2, 0), (0, -2), (0, 2), (-2, -1), (-2, 1), (2, -1), (2, 1), (-1, -2), (-1, 2), (1, -2), (1, 2), (-2, -2), (-2, 2), (2, -2), (2, 2)]
        while stack:
            x, y = stack.pop()
            for direction in directions:
                nx, ny = x + direction[0], y + direction[1]
                if 0 <= nx < len(img) and 0 <= ny < len(img[0]) and not visited[nx][ny] and img[nx][ny] == img[i][j]:
                    visited[nx][ny] = current_label
                    stack.append((nx, ny))

    label_count = 0
    visited = [[0] * len(img[0]) for _ in range(len(img))]
    vals = {}

    # Label each connected region with a unique label
    for i in range(len(img)):
        for j in range(len(img[0])):
            if not visited[i][j]:
                label_count += 1
                visited[i][j] = label_count
                dfs(img, visited, i, j, label_count)
                if label_count not in vals:
                    vals[label_count] = img[i][j]

    # The size of the adjacency matrix should be based on the number of unique labels
    size = label_count
    adjMat = np.zeros((size, size))

    # Construct the adjacency matrix based on the new labels
    for i in range(len(img)):
        for j in range(len(img[0])):
            current_label = visited[i][j]
            if (i-1 >= 0):  # cima
                neighbor_label = visited[i-1][j]
                if current_label != neighbor_label:
                    adjMat[current_label-1][neighbor_label-1] = 1
                    adjMat[neighbor_label-1][current_label-1] = 1
            if (j-1 >= 0):  # esquerda
                neighbor_label = visited[i][j-1]
                if current_label != neighbor_label:
                    adjMat[current_label-1][neighbor_label-1] = 1
                    adjMat[neighbor_label-1][current_label-1] = 1
            if (j+1 < len(img[0])):  # direita
                neighbor_label = visited[i][j+1]
                if current_label != neighbor_label:
                    adjMat[current_label-1][neighbor_label-1] = 1
                    adjMat[neighbor_label-1][current_label-1] = 1
            if (i+1 < len(img)):  # abaixo
                neighbor_label = visited[i+1][j]
                if current_label != neighbor_label:
                    adjMat[current_label-1][neighbor_label-1] = 1
                    adjMat[neighbor_label-1][current_label-1] = 1

    return adjMat, vals
