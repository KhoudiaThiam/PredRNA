from RnaStructure import *
class Rnaparser:
	
	'''
	A class for parsing a file in connect format or parenthestic format

	Methods
	-------
	parse_connect_file:
		A method for parsing connect files
	parse_dotbrackets_file:
		A method for parsing parenthesis files 



	'''
				  
	def parse_connect_file(filename,minimal_loop_length):

		with open (filename,'r') as fileOut:
			lines=fileOut.readlines()
			# initialisation de la sequence 
			arnseq=''
			## liste des pairs 
			List_pair=[]
			list2=[]
			List3=[]
			structure=[]
			# Parcourir les lignes
			for line in lines:
				# Séparer les valeurs
				line=line.split()
				# Récupérer les positions
				pos1=int(line[0])
				pos2=int(line[2])
				#Récuperer la sequence
				arnseq+=line[1]
				
				arnseq=arnseq.upper()
				### eliminer les positions 2 avec 0  
				if pos2 != 0:
					### on ajoute les positions correctes
					List_pair.append((pos1,pos2))
					
			 #### enlever les doublons	  
			for element in List_pair :
				list2.append((min(element),max(element)))
			tuples=sorted(set(list2))
			### on transforme notre sequence obtenue grace au pairs appariement en structure en format parenthese
			# on parcours la  sequence 
			for i in range(len(arnseq)):
				## on ajoute des points de la taille de la sequence  dans notre nouvelle liste 
				List3.append('.')
			## pour chaque couple dans nos tuples 
			for a in tuples:
				# la premier est le debut 
				start=a[0]
				## la deuxieme position est la fin 
				end=a[1]
				## pour tous les "start" on mets une parenthese ouvrante 
				List3[start-1]='('
				### pour tous les "end" on met une parenthese fermante 
				List3[end-1]=')'
			## on convertit  la liste en chaine caractere ==> on obtient notre structure 
			ok=''.join(List3)
			
			## on parcours notre liste de tuples 
			for i in range(len(tuples)):
				# Vérifier que la distance entre les positions est supérieure à 3
				if abs(tuples[i][0]-tuples[i][1])-1 < minimal_loop_length:
					raise ValueError('La distance entre les 2 positions doit être supérieur à la distance minimale')
				if abs(tuples[i][0]-tuples[i][1]) == 0:
					raise ValueError('La distance ne doit pas être la même')
				
				
		### pour eliminer les Croisement
				# avec deux boucles imbrique  on parcours nos tuples de telle sorte
				## qu'on obtient quatre positions ab cd 
			for i,(a,b) in enumerate(tuples):
				for j,(c,d) in enumerate(tuples):
					if i == j :
						continue
					#Croisements interdits : Si (i , j),(k ,l) telles que i<k
					#i<k<l< j ou i< j<k<l
					if (a<c<b<d) or (c<a<d<b):
						raise ValueError('Les appariements se croisent')
			## lien entre Rnastructure et Arnparser		 
			rna=Rnastructure(filename,arnseq,ok)
			#structure.append(rna)   
		return rna
				   
	def parse_dotbrackets_file(filename):
		## ouvrir le fichier 
		with open (filename,'r') as fileOut:
			## on lit tous les lignes 
			lines=fileOut.readlines()
			### on intialise la liste
			List_arn=[]
			## on lit les lignes par 3 
			for i in range(0,len(lines),3):
				## ligne 1 identifiant 
				ID=lines[i].strip()
				## ligne 2 sequence 
				Seq=lines[i+1].strip()
				## ligne 3 la structure en formant parenthese 
				Dot=lines[i+2].strip()
				### lien entre RNastructure et Rnaparser 
				Rna=Rnastructure(ID,Seq,Dot)
				List_arn.append(Rna)
		# pour lire ligne par ligne		 
		for i in List_arn:
			i.check_structure()
			
		return List_arn


		