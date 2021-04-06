import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn import preprocessing
from sklearn import decomposition

def main():
    dataset = datasets.load_wine()

    pca = decomposition.PCA(n_components=2)

    proj = pca.fit_transform(preprocessing.scale(dataset['data']))

    plt.figure(figsize=(4,4))
    plt.scatter(proj[:,0], proj[:,1], c=dataset['target'], marker='+')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
