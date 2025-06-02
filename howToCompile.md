# 🧠 Instruções de Execução

Este repositório permite gerar e classificar grafos temporais a partir de vídeos.

---

## 📦 Conjuntos de Dados Disponibilizados

Os dados estão organizados nas seguintes pastas:

* **`ConjuntoMatriz/`**
  Vídeos sintéticos de baixa complexidade, construídos a partir da aplicação de valores de 1 a 4 (simbolizando cores) em uma matriz 5x5.

* **`ConjuntoKTH/`**
  Vídeos reais do dataset KTH. Foram selecionados os cinco primeiros vídeos das classes *handwaving* e *walking*.

* **`ConjuntoColoracao/`**
  Vídeos sintéticos a partir da coloração manual dos 14 primeiros frames do primeiro vídeo de cada classe do recorte do KTH. Cada frame foi replicado 5 vezes por classe.

---

## ⚙️ Gerar Grafos Temporais

1. Acesse a pasta desejada via terminal:

   ```bash
   cd video2graph/ConjuntoMatriz  # ou ConjuntoKTH, ou ConjuntoColoracao
   ```

2. Para gerar os grafos temporais completos, execute:

   ```bash
   python3 methodMaxFlow.py
   ```

   Ou, para gerar grafos temporais acíclicos:

   ```bash
   python3 methodPixel.py
   ```

3. Os seguintes arquivos serão gerados:

   * `teste_A.txt` – Arestas direcionadas no formato `nó_origem, nó_destino`.
   * `teste_edge_attributes.txt` – Pesos das arestas (identificam relações intra-frame ou temporais).
   * `teste_node_labels.txt` – Valor do fluxo máximo entre pares de nós de subgrafos consecutivos.
   * `teste_graph_indicator.txt` – Atribui um ID a cada grafo.
   * `teste_graph_labels.txt` – Classe de cada grafo.
   * `teste_node_info.txt` – Quantidade de pixels por região (apenas para grafos acíclicos).

---

## 🧪 Classificar Grafos Temporais

4. Clone o [TUDataset Benchmark](https://github.com/chrsmrrs/tudataset):

   ```bash
   git clone https://github.com/chrsmrrs/tudataset
   ```

5. Acesse o diretório do benchmark:

   ```bash
   cd tudataset/tud_benchmark
   ```

6. Copie os arquivos `applyTR.py` e `test.py` disponíveis em `video2graph/tudataset/tud_benchmark/` para esse diretório manualmente ou executando:

   ```bash
   cp [caminho_desde_a_raiz]/video2graph/tudataset/tud_benchmark/applyTR.py .
   cp [caminho_desde_a_raiz]/video2graph/tudataset/tud_benchmark/test.py .
   ```

7. Mova os arquivos dos grafos temporais (passo 3) para `tudataset/tud_benchmark/datasets/teste/teste/raw`.

   Obs: Pode ser necessário criar as pastas `teste/teste/raw`.

8. Execute o script de classificação:

   ```bash
   python3 test.py
   ```

---

## 🔁 Aplicar Redução Transitiva

9. Gere os grafos acíclicos (passo 2) e realize o passo 7.

10. Pelo terminal, acesse a pasta do benchmark:

```bash
cd tudataset/tud_benchmark
```

11. Execute o script de redução transitiva:

```bash
python3 applyTR.py
```

12. Serão gerados dois novos arquivos:

* `teste_A.txt`
* `teste_edge_attributes.txt`

Local: `tudataset/tud_benchmark/output/teste`

13. Substitua os arquivos de mesmo nome da pasta `tudataset/tud_benchmark/datasets/teste/teste/raw` por estes novos:

```bash
mv output/teste/A.txt datasets/teste/teste/raw/teste_A.txt
mv output/teste/edge_attributes.txt datasets/teste/teste/raw/teste_edge_attributes.txt
```

14. Execute o script de classificação:

```bash
python3 test.py
```
