Bonjour,

Lors de ce projet on a implémenter une classe de matrice nommée FredholmMatrix, puis on a essayé d'approcher des DenseMatrix par ces FredholmMatrix. On a ensuite regarder l'erreur en norme Frobenius de la difference entre la DenseMatrix et son approximation.

Après avoir ouvert un terminal et s'être placer dans le bon répertoire, il faut taper :
g++  projet.cpp -o test
./test

Comme ceci vous pourrez visualiser les tests effectues ainsi que les graphess demandés dans l'énoncé.

On remarquera que l'erreur stagne, ou varie doucement à partir d'un certains points, les courbes d'erreurs ressemble à des plateaux. Ainsi dès lors que l'erreur commence à stagner on peut se dire qu'il ne sert à rien de rajouter des vecteurs dans lru et lrv puisque la différence ne sera pas significative. Ainsi dans le graphe de la question 6, on peut se dire que s'arrêter à 10 vecteurs dans lru et lrv est une bonne idée.

Merci d'avance de m'avoir lu.
