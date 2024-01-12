from flask import Flask, request, jsonify
import sys 
import numpy as np
from RnaParser import *

 

class Predict_structure:
	
	# Liste de toutes les structures possibles
	all_structures = []
	
	# Dictionnaire avec les scores de liaison entre chaque paire de bases
	bases_scores = {'AU': 2, 'UA': 2, 'GC': 3, 'CG': 3, 'GU': 1, 'UG': 1}
	
	# Matrice pour stocker les scores de toutes les sous-structures de la séquence
	matrix = None
	
	# Longueur minimale des boucles non-appariées
	minimal_loop_length = 3
	
	# Séquence d'ARN
	rna = None
	
	# Structure prédite
	structure = None
	
	def __init__(self, rnaSeq, minloop):

		"""Initialisation de la classe avec la séquence d'ARN et la longueur minimale des boucles"""

		self.rna = rnaSeq
		self.minimal_loop_length = minloop
		self.matrix = self.fill_mat(self.rna, self.minimal_loop_length)
		
		##=> Programmation dynamique
	
	def fill_mat(self, rna, minimal_loop_length):
		"""Initialisation de la matrice avec les scores de toutes les sous-structures de la séquence"""
		
		n = len(rna)
		## on iniatilise la matrice M avec tout les scores de tout les structures de la sequence
		
		M = [[0] * n for i in range(n)]
	   
		## je parcours  les l valeurs  allant de minimal_loop_length à n 
		for l in range(minimal_loop_length,n):
			#pour chaque L
			#je parcours les indices i et j de  0 à n-l pour remplir les éléments de la matrice M.
			for i in range(n-l):
				j = i + l
			##La valeur de chaque élément M[i][j] est calculée en utilisant la récurrence de l'algorithme de programmation dynamique.
			#M[i][j] est égale à la valeur maximale parmi la somme des scores de deux sous-structures disjointes (M[i][k] + M[k+1][j]) pour chaque k allant de i à j
				M[i][j] = max(M[i][k] + M[k+1][j] for k in range(i, j))
			   
				if (j-i) > minimal_loop_length:
			#et la somme du score de l'appariement entre les nucléotides i et j et 
			#le score de la sous-structure interne (M[i+1][j-1]) pourvu que la taille de la boucle interne 
			#soit plus grande que minimal_loop_length
					M[i][j] = max(M[i][j], M[i+1][j-1] + self.bases_scores.get(rna[i]+rna[j], 0))
				self.all_structures.append((i, j, M[i][j]))

		return M
		
	
	def traceback_bis(self, M, rna, minimal_loop_length, fold=None, i=0, j=None):
		
		 
		"""
		
		Args:
			M (list of lists): Matrix of partition function values.
				rna (str): RNA sequence.
				minimal_loop_length (int): Minimal length of a loop.
				fold (set of tuples, optional): Set of base pairs forming the optimal structure.
				i (int, optional): Starting position of the subsequence.
				j (int, optional): Ending position of the subsequence.

		Returns: the optimal secondary structure in dot-bracket notation 

		"""
	   
		if j is None:
			j = len(rna)-1
		if i >= j:
			 
			return '.'
		##  si la valeur de la matrice de scores en position (i,j)est égale à la valeur en position (i,j-1) ==> CAS A
		if M[i][j] == M[i][j-1]:
			## appelle récursivement traceback_bis en conservant i et diminuant j  de -1 et on ajoute '.' au debut
			return self.traceback_bis(M, rna, minimal_loop_length, fold, i, j-1) + '.'
		#si la valeur de la matrice de scores en position (i,j)est égale à la valeur en position (i+1,j-1)
		##et  que le nombre d'hydrogene de rna[i] + rna[j] est different de none ===> CAS B
		elif self.bases_scores.get(rna[i]+rna[j]) is not None and M[i][j] == M[i+1][j-1]+self.bases_scores.get(rna[i]+rna[j]):
			#alors on appelle récursivement traceback en diminuant j de 1 et en augmentant  i de 1  de plus on ajoute de par 
			##et d'autre respectivement une parenthese ouvrante et fermante pour reconstruire notre structure  
	
			return '(' + self.traceback_bis(M, rna, minimal_loop_length, fold, i+1, j-1) +')'
		else: #==> CAS C 
			### pour chaque k allant de i+1 à j-minimal_loop_length
			for k in range(i+1, j-minimal_loop_length):
				##Soit j est apparié à une position k avec i < k < j, j − k > θ, et des structures se
				##forment alors dans les régions [i , k − 1] et [k + 1, j − 1].
				if self.bases_scores.get(rna[k]+rna[j]) is not None and M[i][j] == self.bases_scores.get(rna[k]+rna[j]) + M[i][k-1] + M[k+1][j-1]:
					#return self.traceback_bis(M, rna, minimal_loop_length, fold, i, k-1)+'('+self.traceback_bis(M, rna, minimal_loop_length, fold, k+1, j-1)+')'
					res=self.traceback_bis(M, rna, minimal_loop_length, fold, i, k-1)+'('+self.traceback_bis(M, rna, minimal_loop_length, fold, k+1, j-1)+')'
					
					
	   
					print("score" ,"=",max(max(M)))
					
					return res
	
	
	def traceback_all(self, M, rna, minimal_loop_length, fold=None, i=0, j=None):
		
		"""
		
		Args:
			M (list of lists): Matrix of partition function values.
				rna (str): RNA sequence.
				minimal_loop_length (int): Minimal length of a loop.
				fold (set of tuples, optional): Set of base pairs forming the optimal structure.
				i (int, optional): Starting position of the subsequence.
				j (int, optional): Ending position of the subsequence.

		Returns:
			list of str: List of optimal secondary structures in dot-bracket notation.
		"""
		if j is None:
			j = len(rna)-1
		## liste qui va prendre toutes les structures optimales 
		res = []
		## i - j est inferieur ou egal minimal loop length  donc pas d'appariement donc on met autant de point que (j-i+1)
		if j-i<=minimal_loop_length:
			return(['.'*(j-i+1)])
		
		##  si la valeur de la matrice de scores en position (i,j)est égale à la valeur en position (i,j-1) cas A 
		## donc il ya pas de liaison entre i et j==> libre 
		if M[i][j] == M[i][j-1]:
			## appelle récursivement traceback en augmentant i de 1 et en conservant j 
			## pour tous s dans notre fonction recursive  on ajoute dans notre liste avec un point car il ya pas d'appariement 
			for s in self.traceback_all(M, rna, minimal_loop_length, fold, i, j-1):
				res.append(s+'.')
				
			
		##CAS  B si la valeur de la matrice de scores en position (i,j)est égale à la valeur en position (i+1,j-1) plus nombre de liaison hydrogene de  rna[i]+rna[j]
		## L(Si , j)=L(Si+1 , j−1)+α(ri,rj)
		if self.bases_scores.get(rna[i]+rna[j]) is not None and M[i][j] == M[i+1][j-1]+self.bases_scores.get(rna[i]+rna[j]):
			#pour tous s dans notre fonction recursive  on ajoute dans notre liste une parenthese ouvrante ,  s et une  parenthese fermante  car il ya  appariement entre i et j  
			for s in self.traceback_all(M, rna, minimal_loop_length, fold, i+1, j-1):
				res.append('(' + s +')')
		##### pour chaque k allant de i+1 à j-minimal_loop_length	CAS C
		for k in range(i+1, j-minimal_loop_length):
			# si j est  apparié à une position k avec i < k < j et , j − k > θ  et si la valeur de la matrice de scores en position (i,j) est egal à la valeur en position (k+1,j-1) + la nb liaison d'hydro de rna[k]+rna[j]
			if self.bases_scores.get(rna[k]+rna[j]) is not None and M[i][j] == self.bases_scores.get(rna[k]+rna[j]) + M[i][k-1] + M[k+1][j-1]:
				# on stocke dans subres1 les structures se forment  alors dans les régions [i , k − 1]
				subres1 = self.traceback_all(M, rna, minimal_loop_length, fold, i, k-1)
				# on stocke dans subres2 les structures qui se forment  alors dans les régions  [k + 1, j − 1]
				subres2 = self.traceback_all(M, rna, minimal_loop_length, fold, k+1, j-1)
				## pour chaque S1 dans subres1 
				for s1 in subres1:
					## pour chaque S2 dans subres2 
					for s2 in subres2:
					#### si fold n'est pas  et que (i j) est pas  dans fold et   et que tous les u et v ne sont pas fold aussi  
						if not fold or (i, j) in fold or all((u, v) not in fold 
															for u in range(i, k+1) for v in range(k+1, j+1)):
							# ajoute  dans notre liste encore des parenthese 
		

							res.append(s1 + '(' + s2 + ')')
		return res
app = Flask(__name__)

@app.route('/traitement', methods=['POST'])


def traitement():
	if 'fichier' not in request.files:
		return jsonify({'message': 'Aucun fichier trouvé dans la requête'}), 400
	
	fichier = request.files['fichier']
	if fichier.filename == '':
		return jsonify({'message': 'Aucun fichier sélectionné'}), 400
	
	rnaSeq = fichier.read().decode('utf-8')  # Lecture du fichier
	
	minimal_loop_length = 3
	a = Predict_structure(rnaSeq, minimal_loop_length)
	b = a.fill_mat(rnaSeq, minimal_loop_length)
	test = a.traceback_all(b, rnaSeq, minimal_loop_length)
	
	# Écriture des résultats dans un fichier de sortie
	output = 'output.txt'
	with open(output, 'w') as file:
		out = 'Le nombre de structures optimales possibles est de : {} \n'.format(len(test))
		for i in range(len(test)):
			out += test[i] + "	  score = " + str(max(max(b))) + " liaisons" + '\n'
			file.write(out)
	contenu = ''
		
	with open(output, 'r') as file:
		contenu = file.readlines()
		
	return jsonify({'output': contenu})


@app.route('/validation', methods=['POST'])

def validation():
    if 'fichier' not in request.files:
        return jsonify({'message': 'Aucun fichier trouvé dans la requête'}), 400
    
    fichier = request.files['fichier']
    if fichier.filename == '':
        return jsonify({'message': 'Aucun fichier sélectionné'}), 400
    
    minimal_loop_length = 3
    uploaded_filepath = 'input.txt'
    fichier.save(uploaded_filepath)

    a = Rnaparser.parse_connect_file(uploaded_filepath, minimal_loop_length)
    
    output_filepath = 'output.txt'
    with open(output_filepath, 'w') as file:
        file.write('Number of hydrogen bonds : {} \n'.format(a.Hydrogen_bonds()))
        file.write(str(a))
    
    contenu = ''
    with open(output_filepath, 'r') as file:
        contenu = file.readlines()
    
    return jsonify({'output': contenu})



if __name__ == "__main__":
	app.run(debug=True)
