import os
from PIL import Image

def escala_de_cinza(imagem):
    # Convertendo para escala de cinza
    return imagem.convert('L')

def processar_imagem(nome_arquivo, pasta_origem, pasta_destino):
    # Abrindo a imagem
    imagem = Image.open(os.path.join(pasta_origem, nome_arquivo))
    
    # Aplicando escala de cinza
    imagem_cinza = escala_de_cinza(imagem)
    
    # Salvando a nova imagem na pasta de destino
    imagem_cinza.save(os.path.join(pasta_destino, nome_arquivo))

def main(pasta_origem, pasta_destino):
    # Verificando se a pasta de destino existe, se não, cria ela
    if not os.path.exists(pasta_destino):
        os.makedirs(pasta_destino)
    
    # Percorrendo os arquivos na pasta de origem
    for nome_arquivo in os.listdir(pasta_origem):
        if nome_arquivo.endswith('.ppm'):
            # Processando apenas arquivos do tipo .ppm
            processar_imagem(nome_arquivo, pasta_origem, pasta_destino)

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Uso: python script.py <pasta_origem> <pasta_destino>")
        sys.exit(1)
    
    pasta_origem = sys.argv[1]
    pasta_destino = sys.argv[2]
    main(pasta_origem, pasta_destino)
