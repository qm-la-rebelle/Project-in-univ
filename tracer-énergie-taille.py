#SUN Qimeng, Mécanique quantique et topologie : modèle SSH, tracer l'énergie dans le gap en fonction de la taille du système
import numpy as np
import pylab as pl


a=range(1,200,10)
d=np.arange(0.,1.,0.05)
Egapi=0
iEgap=0

fra=0.1
for L in a:
    # initializer la matrice de Hamiltonien
    hamil=np.zeros((2*L,2*L))
    b=range(L)
    c=range(L-1)
    for l in b:
        hamil[2*l,2*l+1]=fra
        hamil[2*l+1,2*l]=fra
    for l in c:
        hamil[2*l+1,2*l+2]=1
        hamil[2*l+2,2*l+1]=1
    # calculer les valeurs propres de H
    vap=np.linalg.eigvals(hamil)
    #La valeur d'index de la valeur du tableau de petite à grande
    ord=np.argsort(vap)
    #Utilisez la valeur de l'indice pour reorganiser toutes les énergies propres dans l'ordre croissant
    vapo=vap[ord]
    #sélectioner les deux énergies dans le gap
    Egap1=vapo[L]
    Egap2=vapo[L-1]
    Egapi=np.hstack((Egapi,Egap1))
    iEgap=np.hstack((iEgap,Egap2))

Egapi = np.delete(Egapi,0,axis = 0)# supprimer la ligne non sens
iEgap = np.delete(iEgap,0,axis = 0)# supprimer la ligne non sens
# convertir les énergies dans le gap en logarithme
Egapi=np.log(np.abs(Egapi[0:]))
iEgap=np.log(np.abs(iEgap[0:]))
# tracer l'énergie dans le gap en fonction de la taille du système
pl.plot(a,Egapi,"b", label="l'énergie propre du l'état de bord -")
pl.plot(a,iEgap,"ro", label="l'énergie propre du l'état de bord +")
pl.title('énergi des états dans le gap en fonction de\n la taille du système avec la fraction de t1/t2 = {}'.format(fra))
pl.xlabel('la taille du système')  
pl.ylabel('énergi des états dans le gap')
pl.legend()
pl.tight_layout()
pl.savefig("exp 0.1.png")
pl.show()

