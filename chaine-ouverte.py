#SUN Qimeng, Mécanique quantique et topologie : modèle SSH, calculer le spectre de la chaîne ouverte
import numpy as np
import pylab as pl

L = 50
hamil=np.zeros((2*L,2*L))
b=range(L)
c=range(L-1)
d=np.arange(0.,5.,0.1)
E=np.zeros(2*L)

for fra in d:
  # initializer la matrice de Hamiltonien
  for l in b:
    hamil[2*l,2*l+1]=fra
    hamil[2*l+1,2*l]=fra
  for l in c:
    hamil[2*l+1,2*l+2]=1
    hamil[2*l+2,2*l+1]=1
  # calculer les valeurs propres de H
  E1=np.linalg.eigvals(hamil)
  ord=np.argsort(E1)
  E1=E1[ord]
  # sauvegarder les valeurs propres calculées correspondant à la fraction de t1/t2
  E=np.vstack((E,E1))
E = np.delete(E,0,axis = 0)# supprimer la ligne non sens
# tracer le spectre
print(np.size(d))
print(np.size(E))
print(E)
e=range(np.size(d))
for n in e:
    xn=np.ones(2*L)*d[n]
    yn=E[n]
    pl.scatter(xn,yn)
pl.title('Spectre de chaîne ouverte avec {} cellules élémentaires '.format(L))
pl.xlabel('fraction de t1/t2')
pl.ylabel('énergie(valeur propre)/t2')
pl.grid()
pl.savefig("spe-ouv.png")
pl.show()
