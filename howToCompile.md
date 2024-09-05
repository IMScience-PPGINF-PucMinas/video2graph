<!-- COMPILAÇÃO -->

Para compilar GFT:

1. Para executar o código:

   python3 methodMaxFlow.py

2. Serão gerados os seguintes arquivos:

   A.txt -> Aresta direcionada no estilo "Nó_origem, Nó_destino"
   edge_attributes.txt -> Peso da aresta. Indica se a aresta pertence a um frame ou se representa uma relação temporal
   node_labels.txt -> Apresenta o valor do fluxo máximo calculado entre um par de nós de dois subgrafos (frames) subsequentes
   graph_indicator.txt -> Atribui um id para cada grafo
   graph_labels.txt -> Apresenta a classe de cada grafo
