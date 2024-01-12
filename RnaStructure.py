
class Rnastructure:
	"""
	A class to represent an RNA, its identifier and its parenthetic structure

	Attributes
	----------
	
	seq_ID : str
	inrna: str
	structure : str 

	Methods
	-------
	get_structure:
		return the parenthesis structure of the RNA 
	check_rna:
		Raise Error if the sequence is not a RNA sequence
	check_structure :
		Raise Error if the structure is not acceptable
	Hydrogen bonds :
		Compute the number of hydrogen bonds in the RNA
	


	"""
	
	def __init__(self,seq_id,inrna,structure):
	   
		self.__seq_id=seq_id
		self.__inrna=inrna.upper()
		self.__structure=structure
		self.listop=[]
		## dictonnaire qui va prendre nos piles 
		self.pairs_dic={}
		
		self.check_structure()
		
		## on intialise notre pile 
		stack=[]
		## on intilaise la liste des paires 
		pairs=[]

		## on enumere notre structure de sorte qu'on peut obtenir l'indice ete l'element  
		for i,c in enumerate(self.__structure):
			##  si c est une parenthese ouvrante donc 
			if c == '(':
				#on ajoute l'indice de  c dans la pile 
				stack.append(i)
			## si c est une parenthese fermante  donc 
			elif c == ')':
				## on verifie si notre liste est n'est pas  si c'est le cas oon retourne zero  
				if len(stack)==0:
					return False
				## sinon   
				else :
					## on pop notre pile pour avoir l'indice des parenthese fermante 
					stock=stack.pop()
					## on stock nos positions appariés  dans la dictionnaire 
					self.pairs_dic[i]=stock
					### et on ajoute ici nos tuples de pairs appariés 
					pairs.append((i,stock))
		### pour chaque i dans pairs on garde les minimum et les maximum pour eviter les doublons		   
		for i in pairs :
			self.listop.append((min(i),max(i)))

		## on trie notre liste de tuples 
		self.listop=sorted(set(self.listop))

	def get_structure(self):
		
		return self.__structure


	def check_rna(self):
		
		rna_alphabet='AUCGaucg'

		for i in self.__inrna:
			if i not in self.__rna_alphabet:
				return('Invalid RNA character in sequence')
			
	def check_structure(self):

		#Number of parentesis

		## intialise  les  nombre parenthese ouvrant et fermante  à zero 
		po=0
		pf=0
		# pour tout i dans la la structure 
		for i in range(len(self.__structure)):
			## si i est une parnethese ouvrante 

			if self.__structure[i] == '(':
				## on incremente  de 1  po 
				po+=1
			## si i est une parnethese est fermante 
			elif self.__structure[i]==')':
				## on incremente  de 1  po 
				pf+=1
		## si  le nombre po est different de  pf 
		if pf != po :
			# donc on retourne une erreur 
			raise ValueError('The number of opening parenthesis and closing parenthesis are not equal')
		
		#checking for the minimal loop length  

		###  pour chaque postion (pos1 et pos2) dans notre dictionnaire appariement 
		
		for pos1,pos2 in self.pairs_dic.items():
			## on verifie qui il y a au moins trois base libre dans boucle  si ce n'est pas le cas 
			if abs(pos1-pos2)-1 < 3:
				### on retourne une erreur 
				raise ValueError("La distance entre les positions doit être supérieure à 3.")
			## on verifie aussi si nos deux positions ne sont egale 
			if pos1 == pos2:
				 raise ValueError("Les positions ne peuvent pas être les mêmes.")
	
		#Checking the bases pairing
		## on s'assure qu'il ya pas d'appariement entre deux base simillaire  
		for i in self.listop:

			if (self.__inrna[i[0]] == 'G' and self.__inrna[i[1]] == 'C') or (self.__inrna[i[0]] == 'C' and self.__inrna[i[1]] == 'G'):
				pass
			elif (self.__inrna[i[0]] == 'A' and self.__inrna[i[1]] == 'U') or (self.__inrna[i[0]] == 'U' and self.__inrna[i[1]] == 'A'):
				pass
			elif (self.__inrna[i[0]] == 'G' and self.__inrna[i[1]] == 'U') or (self.__inrna[i[0]] == 'U' and self.__inrna[i[1]] == 'G'):
				pass
			else:
				raise ValueError('The bases pairing is not respected')

		
	### calcule des nombres de liaisons d'hydrogenes	
  
	def Hydrogen_bonds(self):

		h_bond=0
  
		for i in range(len(self.listop)): 
			if (self.__inrna[self.listop[i][0]] == 'C' and self.__inrna[self.listop[i][1]] == 'G') or (self.__inrna[self.listop[i][0]] == 'G' and self.__inrna[self.listop[i][1]]=='C'):
				h_bond+=3

			elif (self.__inrna[self.listop[i][0]] == 'A' and self.__inrna[self.listop[i][1] ]== 'U') or (self.__inrna[self.listop[i][0]] == 'U' and self.__inrna[self.listop[i][1] ]=='A'):
				h_bond+=2

			elif (self.__inrna[self.listop[i][0]] == 'G' and self.__inrna[self.listop[i][1]]== 'U') or (self.__inrna[self.listop[i][0]] == 'U' and self.__inrna[self.listop[i][1] ]=='A'):
				h_bond+=1

		return h_bond

	def __str__(self):
		
		return "{}\n{}\n{}".format(self.__seq_id , self.__inrna ,self.__structure)
	
