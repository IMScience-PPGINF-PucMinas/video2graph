<!-- COMPILAÇÃO -->

Para compilar GFT:

1. Copiar os arquivos da "gft" em alguma pasta do seu computador.

2. Criar uma variável de ambiente GFT_DIR e colocar no .bashrc:

   export export GFT_DIR=/home/danielle/gft

3. Entrar na pasta da "gft" e fazer:
   
   make clean
   
   make

4. Se ocorrer o erro "fatal error: zlib.h: No such file or directory", então você tem que instalar o pacote do zlib: zlib1g-dev
   Repita novamente o passo 3.

5. Para compilar o Makefile:

   make -f MakefileDisf.make clean
   
   make -f MakefileDisf.make

6. Converter imagens para escala de cinza (opcional):
    
   python3 toGray.py _pastaOrigem_ _pastaDestino_
   
   (ex: python3 toGray.py dataset/girl dataset/grayGirl)

7. Para compilar e excutar o código:

   make clean
   
   make

8. Para criar as relações temporais dos videos:
   
   python3 node.py
