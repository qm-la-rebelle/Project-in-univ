# SUN Qimeng, Mécanique quantique et topologie : modèle SSH, calculer le nombre d'enroulement
import numpy as np
import matplotlib.pyplot as plt


def up(k, fra):  # définir l'état propre + du modèle SSH
    v = fra + np.exp(1j * k)
    return 1 / np.sqrt(2) * np.array([1, v / np.abs(v)])


def um(k, fra):  # définir l'état propre - du modèle SSH
    v = fra + np.exp(1j * k)
    return 1 / np.sqrt(2) * np.array([1, -v / np.abs(v)])


def derup(k, fra, eps):  # définir la dérivée de l'état propre + du modèle SSH
    return (up(k + eps, fra) - up(k - eps, fra)) / (2 * eps)


def derum(k, fra, eps):  # définir la dérivée de l'état propre - du modèle SSH
    return (um(k + eps, fra) - um(k - eps, fra)) / (2 * eps)


Np = 10000
eps = 1e-15
k = np.linspace(-np.pi, np.pi, Np)
d = np.arange(0., 5., 0.1)
nm = np.zeros(2)
nnp = np.zeros(2)
for fra in d:  # calculer le nombre d'enroulement
    # calculer le nombre d'enroulement pour l'état propre + du modèle SSH
    N = 0
    for i in range(Np):
        v1 = up(k[i], fra)
        dv1 = derup(k[i], fra, eps)
        dk = 2 * np.pi / Np
        res = - (np.conj(v1) @ dv1) * 1j * dk / np.pi  # calculer l'intégral
        N = N + res
    # calculer le nombre d'enroulement pour l'état propre - du modèle SSH
    n = 0
    for j in range(Np):
        v1 = um(k[j], fra)
        dv1 = derum(k[j], fra, eps)
        dk = 2 * np.pi / Np
        res = - (np.conj(v1) @ dv1) * 1j * dk / np.pi  # calculer l'intégral
        n = n + res
    nmp = [np.real(n), np.imag(n)]
    npp = [np.real(N), np.imag(N)]
    nm = np.vstack((nm, nmp))
    nnp = np.vstack((nnp, npp))
print(nnp,nm)
nnp = np.delete(nnp, [0,11], axis=0)# supprimer la ligne non sens
nm = np.delete(nm, [0,11], axis=0)# supprimer la ligne non sens
print(nnp,nm)
# tracer le nombre d'enroulement
plt.plot(nnp[:, 0], nnp[:, 1], "ro", label="l'état propre + du modèle SSH")
plt.plot(nm[:, 0], nm[:, 1], "k", label="l'état propre - du modèle SSH")
plt.title('le nombre d enroulement avec ')
plt.xlabel('la partie réelle du nombre d enroulement ')
plt.ylabel('la partie imaginaire du nombre d enroulement ')
plt.legend()
plt.tight_layout()
plt.savefig("nbr enrou.png")
plt.show()
