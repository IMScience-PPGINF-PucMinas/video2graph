def ler_arquivo_grafo(nome_arquivo):
    grafo = {}
    with open(nome_arquivo, 'r') as arquivo:
        for linha in arquivo:
            no_origem, no_destino = map(int, linha.strip().split(','))
            if no_origem not in grafo:
                grafo[no_origem] = []
            grafo[no_origem].append(no_destino)
    return grafo

def ler_arquivo_pertencimento(nome_arquivo):
    pertencimento = {}
    with open(nome_arquivo, 'r') as arquivo:
        for indice, linha in enumerate(arquivo, start=1):
            grafo = int(linha.strip())
            pertencimento[indice] = grafo
    return pertencimento

def ler_arquivo_indicador(nome_arquivo):
    indicadores = {}
    with open(nome_arquivo, 'r') as arquivo:
        for indice, linha in enumerate(arquivo, start=1):
            id_video = int(linha.strip())
            indicadores[indice] = id_video
    return indicadores

def escrever_resultados_em_arquivo(resultados):
    with open("./graph/node_labels.txt", "a") as arquivo:
        if resultados:
            for par in resultados:
                for elemento in par:
                    arquivo.write(str(elemento) + ", ")
                arquivo.write("\n")

def escrever_pares_em_arquivo(no, outro_no):
    with open("./graph/A.txt", "a") as arquivo:
        arquivo.write(f"{no},{outro_no}\n")
        arquivo.write(f"{outro_no},{no}\n")

def escrever_tempo_em_arquivo(tempo):
    with open("./graph/edge_attributes.txt", "a") as arquivo:
        arquivo.write(f"{tempo}\n{tempo}\n")

def comparar_arestas_proximas(arquivo_grafo, arquivo_pertencimento, arquivo_indicador):
    grafo = ler_arquivo_grafo(arquivo_grafo)
    pertencimento = ler_arquivo_pertencimento(arquivo_pertencimento)
    indicadores = ler_arquivo_indicador(arquivo_indicador)

    resultados = []

    outro_nos_escolhidos = set()

    m = 3

    for no, vizinhos in grafo.items():
        nos_escolhidos = set()
        grafo_atual = pertencimento.get(no, 0)
        id_video_atual = indicadores.get(no, 0)
        for outro_no, outro_vizinhos in grafo.items():
            outro_grafo = pertencimento.get(outro_no, 0)
            id_video_outro = indicadores.get(outro_no, 0)
            if outro_grafo == grafo_atual + 1 and id_video_atual == id_video_outro:
                arestas = len(vizinhos)
                outro_arestas = len(outro_vizinhos)
                a = arestas - outro_arestas
                b = outro_arestas - arestas
                c = a if (a >= b) else b
                if (c <= m):
                    if (no in nos_escolhidos) or (outro_no in outro_nos_escolhidos):
                        continue
                    else:
                        nos_escolhidos.add(no)
                        outro_nos_escolhidos.add(outro_no)
                        resultados.append((grafo_atual, no, outro_grafo, outro_no))
                        escrever_pares_em_arquivo(no, outro_no)
                        escrever_tempo_em_arquivo(grafo_atual)
        if (no not in nos_escolhidos):
            resultados.append((grafo_atual, no))

    return resultados

resultados = comparar_arestas_proximas("./graph/A.txt", "./graph/frame_indicator.txt", "./graph/graph_indicator.txt")
escrever_resultados_em_arquivo(resultados)
