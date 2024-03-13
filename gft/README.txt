Instruços de uso:

1) Copiar os arquivos da "gft" em alguma pasta do seu computador.
No meu caso eu coloquei em: /home/pmiranda/prog/lib/gft

2) Criar uma variável de ambiente GFT_DIR e colocar no .bashrc:
export GFT_DIR=/home/pmiranda/prog/lib/gft

3) Entrar na pasta da "gft" e fazer:
make clean
make

4) Se ocorrer o erro "fatal error: zlib.h: No such file or directory", então você tem que instalar o pacote do zlib: zlib1g-dev
Repita novamente o passo 3.

5) Entre na pasta do programa que usa a gft e faça:
make clean
make


export GFT_DIR=/home/carolina/Documentos/gft