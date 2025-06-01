import auxiliarymethods.auxiliary_methods as aux
import auxiliarymethods.datasets as dp
import kernel_baselines as kb
from auxiliarymethods.kernel_evaluation import kernel_svm_evaluation
import timeit
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# Download dataset.


def calculateGraphLet(name):
    all_matrices = []
    gm = kb.compute_graphlet_dense(name, True, False)
    gm = aux.normalize_gram_matrix(gm)
    all_matrices.append(gm)
    classes = dp.get_dataset(name)
    
    np.savetxt("GRAM", gm)
    
    print(gm.shape, gm.size)
    null_data = 100 * np.count_nonzero(gm) / gm.size
    print(f"Percentual de dados não nulos: {null_data:.2f}%")

    # Salvar o gráfico como JPEG
    plt.matshow(gm, cmap=matplotlib.cm.Greys)
    plt.savefig(f'{name}_graphlet.jpg')
    plt.close()

    # kernel_svm_evaluation(all_matrices, classes,
    #                         num_repetitions=10, all_std=False)
    print("GL: ")
    print(kernel_svm_evaluation(all_matrices, classes,
                            num_repetitions=10, all_std=False))
    
    print()


def calculateW1(name: str):
    all_matrices = []
    for i in range(1, 6):
        gm = kb.compute_wl_1_dense(name, i, True, False)
        gm = aux.normalize_gram_matrix(gm)
        all_matrices.append(gm)
        # print(gm.shape)
    classes = dp.get_dataset(name)
    # kernel_svm_evaluation(all_matrices, classes,
    #                             num_repetitions=10, all_std=False)

    print(gm.shape, gm.size)
    null_data = 100 * np.count_nonzero(gm) / gm.size
    print(f"Percentual de dados não nulos: {null_data:.2f}%")

    # Salvar o gráfico como JPEG
    plt.matshow(gm, cmap=matplotlib.cm.Greys)
    plt.savefig(f'{name}_wl.jpg')
    plt.close()

    print("WL: ")
    print(kernel_svm_evaluation(all_matrices, classes,
                                num_repetitions=10, all_std=False))
    
    print()


for name in os.listdir("./datasets"):
    print(name)
    timeGL = timeit.timeit('calculateGraphLet(name)', 'from __main__ import calculateGraphLet, name', number= 1)
    timeW1 = timeit.timeit('calculateW1(name)', 'from __main__ import calculateW1, name', number= 1)
    print(name + "\nTime GL: " + str(timeGL) + "    Time W1: " + str(timeW1) + "\n")
