"""Generateur de séquences ayant un score ForSA ameliore a partir de sequences de Blocs Proteiques

Utilise la matrice d'occurence de FoRSA pour determiner la probabilite qu un acide amine soit positionne
a un une position donne.

===============================================================================

Usage :
	python3 Generateur.py matrice_score.format maseq.pbseq maseq.aaseq

"""
import sys, os, copy, random
import subprocess
import numpy as np
import forsaGlobal as fGpy
import multiprocessing
import time

#==============================================================================================================================#
# Utilitaires 

def reformate(file):
	"""Modifie la matrice d occurence du code forsa_global.c en un fichier csv classique
	avec l extension .format
	"""
	with open(file, 'r') as fil, open('matrice_score.format', 'w') as out:
		for line in fil:
			line = line.strip()
			out.write("{}\n".format(line[1:-2])) #retire la premiere et derniere accolade ainsi que la virgule de fin


def lit_Matrice(file):
	"""Lit la matrice d'occurence .format et renvoie un dicctionnaire qui la contient
		-Args:
			_file .format file qui contient la matrice d'occurence
	Le dictionnaire est de la forme mat[PB] = {}
	clés = Blocs proteiques (lettres de a->p + x et z)
	valeurs = dictionaire dont les clés sont A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
	et les valeurs sont des listes, chaque liste contient 15 valeurs correspondants aux probabilité de l'aa pour chaque 
	"""
	if file.split('.')[-1] != 'format':
		sys.exit("Erreur de format de fichier\nIndiquer le chemin vers un fichier .format qui contient la matrice d'occurences")
	PB = 'ABCDEFGHIJKLMNOPXZ' #blocs proteiques de a -> p
	AA = 'ACDEFGHIKLMNPQRSTVWY' #acides amines code une lettre
	indice = 0 # variable qui garde en memoire la progression dans le fichier
	pointePB = 0 # compteur qui donne la lettre PB correspondant à ce qui est lu dans le fichier
	pointeAA = 0 # compteur acide amine
	matPB = {} # dictionnaire qui contient la matrice
	dAA = {} # dictionnaire qui contient les cles AA pour les listes de 15 valeurs 
	with open(file,'r') as fmt:
		for i,line in enumerate(fmt):
			if i >= indice + 20: # toutes les 20 lignes 
				matPB[PB[pointePB]] = copy.deepcopy(dAA) # enregistre le dictionnaire precedent
				indice += 20 
				pointePB +=1 # pointe vers le PB suivant
				pointeAA = 0 # remet le compteur acide amine a 0
				dAA = {} # efface le contenu (pas obligatoire)
				matPB[PB[pointePB]] = dAA # assigne au nouveau PB un dictionnaire d'acide amine vide
			line = line.strip().split(',')
			for i,elem in enumerate(line):
				line[i] = float(elem)
			dAA[AA[pointeAA]] = line
			pointeAA +=1
	return matPB

def lit_seq_AA(file):
	"""Lit une ou plusieurs sequences dans un fichier .aaseq et les renvoient sous forme d'un dictionnaire 
		-Args:
			_file : .aaseq file qui contient les sequences en format aaseq
		le dictionnaire sera de la forme aaSeq[identifiant] = 'ACDEF..'
	"""
	if file.split('.')[-1] != 'aaseq':
		sys.exit("Erreur de format de fichier\nIndiquer le chemin vers un fichier .aaseq qui contient une ou plusieurs sequences d'acide aminés")
	aaSeq = {}
	with open(file, 'r') as seq:
		for line in seq:
			if line.startswith('>'):
				identifiant = line.strip()[1:].split('_')[0]
				continue
			line = line.strip()
			aaSeq[identifiant] = line
	return aaSeq

def lit_seq_PB(file):
	"""Lit une ou plusieurs sequences dans un fichier .pbseq et les renvoient sous forme d'un dictionnaire
		-Args:
			_file : .pbseq file qui contient les sequences en format pbseq
		le dictionnaire sera de la forme pbSeq[identifiant] = 'PBSEQ..'
	"""
	if file.split('.')[-1] != 'pbseq':
		sys.exit("Erreur de format de fichier\nIndiquer le chemin vers un fichier .pbseq qui contient une ou plusieurs sequences de PB")
	pbSeq = {}
	with open(file, 'r') as seq:
		for line in seq:
			if line.startswith('>'):
				identifiant = line.strip()[1:].split('_')[0]
				continue
			line = line.strip()
			pbs = ''
			for i,elem in enumerate(line): # modifie les Z en X puisqu'ils aprtagent les mêms valeurs
				if elem == 'Z':
					pbs += 'X'
				else:
					pbs += elem
			pbSeq[identifiant] = pbs
	return pbSeq

def writeSeq(seqAA):
	"""Ecrit dans un fichier aaseq la nouvelle sequence
		-Args:
			_seqAA : chaine de caractères contenant la nouvelle sequence
	"""
	with open('newseq.aaseq' , 'w') as out:
		out.write('>newseq_by_Generateur.py\n')
		out.write(seqAA)
		out.write('\n')

def ReadFoRSA(file):
	"""Lit un fichier generé par FoRSA et renvoie le Zscore 
		-Args :
			_file: nom du fichier qui contient les informations
	"""
	with open(file,'r') as src:
		for line in src:
			if line.startswith('raw'):
				Zscore = float(line.strip().split(':')[-1])
				return Zscore

def writeFileOfNewSeq(memory,nSeqPop):
	"""Ecrit dans un fichier .txt les sequences améliorées et leurs Zscores associés
		-Args :
			_memory: dictionnaire contenant la sequence et le Zscore associé 
	"""
	with open('./out/MesSequencesAmeliorees.txt','w') as out:
		ecart = len(list( memory.values() )[0][0]) -8
		out.write('ID   Sequence-I{}Sequence-F{}Z-score\n'.format(' '*ecart,' '*ecart))
		for i,key in enumerate(memory.keys()):
			out.write('{}  {}  {}  {:.3f}\n'.format(memory[key][-1],key,memory[key][0],memory[key][1]))
			# if i == nSeqPop:
			# 	out.write('Population ci dessus : {} sequences\n'.format(nSeqPop))

def writeDataMutation(memory):
	"""Ecrit dans un fichier les differences entre le Zscore initial d une sequence et les Z scores obtenu en mutant chaque position par chaque acide aminé
		-Args :
			_memory: liste de listes contenant 20 DeltaZscores pour chaque position
	"""
	fileName = 'DataMutation.txt'
	if not os.path.isfile('DataMutation.txt'):
		mode = 'w'
	else : 
		mode = 'a'
	with open(fileName, mode) as out:
		# out.write('itération : {}\n'.format(iteration))
		if mode == 'w':
			out.write('A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y\n')
		for i in range(len(memory)):
			trns = [str(elem) for elem in memory[i]]
			line = ",".join(trns)
			out.write(line + '\n')

def vAbs(a,b):
	"""Calcule la valeur absolue de la difference entre a et b
	"""
	if a > b:
		return a - b
	return b - a

def is_equal(a,b):
	"""Teste l'égalité entre deux floats
	"""
	d = vAbs(a,b)
	if d > 10**(-6):
		return False
	return True

def writeEvolScores(keepscore):
	with open('./out/Evolution_Scores.csv', 'w') as out:
		for key in keepscore.keys():
			strlist = [str(elem) for elem in keepscore[key]]
			out.write("{},{}\n".format(key,",".join(strlist)))

def writeDiffScoreMute(memory,idSeq):
	"""Ecrit dans un fichier correspondant a la bonne séquence les differences entre le Zscore initial d une sequence et les Z scores obtenu en mutant chaque position par chaque acide aminé
		-Args :
			_memory: liste de listes contenant 20 DeltaZscores pour chaque position
			_idseq: l'indice de la sequence
	"""
	fileName = './out/{}.csv'.format(idSeq)
	if not os.path.isfile(fileName):
		mode = 'w'
	else : 
		mode = 'a'
	with open(fileName, mode) as out:
		if mode == 'w':
			out.write('mutation,score\n')
		else:
			out.write('//,//\n')
		for pos in memory:
			for mut in pos.keys():
				out.write('{},{}\n'.format(mut,pos[mut]))

#==============================================================================================================================#
# Approche Graine

def fixeAA(PB,occ,index,Rand):
	"""Donne pour un PB la sequence optimale d'acide aminé d'après la matrice d'occurence
		-Args:
			_PB : un Bloc Proteique sous forme de char
			_occ : la matrice d'occurence transformée
			_index : permet de savoir si l'on s'occupe d'une (0,14) ou quinze positions (15) 
	La 'graine' donnera une suite de 15 aa sous forme de char, ensuite une seule position a la fois
	"""
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	if index != 15:
		# maxi = -100
		lmax = []
		for i in range(20):
			aa = AA[i]
			lmax.append(occ[PB][aa][index])
			# if occ[PB][aa][index] > maxi:
			# 	maxi = occ[PB][aa][index]
			# 	keep = aa
		k = random.choice(np.argsort(lmax)[19-Rand:])
		return AA[k]
	seqOpt = '' # sequence dont les acides amines sont optimaux a tt les positions pour PB
	for i in range(15):
		# maxi = -100
		lmax = []
		for j in range(20):
			aa = AA[j]
			lmax.append(occ[PB][aa][i])
			# if occ[PB][aa][i] > maxi:
			# 	maxi = occ[PB][aa][i]
			# 	keep = aa
		k = random.choice(np.argsort(lmax)[19-Rand:])
		seqOpt += AA[k]
	return seqOpt

def graine(pbseq,occ,Rand):
	"""Determine la position de la graine de depart sur la sequence de PB
	"""
	if len(pbseq) < 15:
		sys.exit('Impossible pour une sequence de moins de 15 acides aminés...')
	pGraines = [] # comporte toutes les positions possibles pour la graine
	ori = (len(pbseq) / 2)
	npositions = len(pbseq) - 14
	newseqs = []
	for g in range(npositions):
		nseq = fixeAA(pbseq[g+7],occ,15,Rand)
		pref,suff = '',''
		for p in range(npositions):
			if p == g:
				continue
			if p < g:
				pref += fixeAA(pbseq[p+7],occ,0,Rand)
			if p > g:
				suff += fixeAA(pbseq[p+7],occ,14,Rand)
		newseqs.append(pref + nseq + suff)
	return newseqs




#==============================================================================================================================#
# Somme

def matrice_sequence(occ,pbs):
	"""Construit et retourne une matrice qui contient les informations de la sequence en PB et de la matrice d'occurence
		-Args:
			_occ : anciennement matPB, la matrice d'occurence transformée
			_pbs : une sequence de PB sous forme de chaine de caractères
	Cette matrice aura une taille n*20. n pour le nombre de positions, et 20 valeurs associées aux 20 acides aminés
	Ces valeurs sont obtenues en moyennant sur les 15 positions de la fenetre.
	Exemple : une alanine au centre d'une sequence de 15 PB de long aura une valeur associée issue de la somme de 
	15 valeurs (de la matrice d'occurence) correspondant a l'occurence pour chaque PB de l'alanine à une position donnée
	"""
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	m_seq = []
	for i in range(len(pbs)):
		m_seq.append([])
		for j in range(20):
			m_seq[i].append([0,0])

	for i in range(len(pbs)):
		pbCurr = pbs[i]
		for j in range(20):
			aa = AA[j]
			m_seq[i][j][0] += occ[pbCurr][aa][7]
			m_seq[i][j][1] += 1
			iterateur = 0 # permet d'incrementer sur la fenetre (de 0 à 15 ) sauf si on est sur les bords
			while iterateur < 7 and (i+iterateur) < len(pbs):
				if iterateur == 0:
					iterateur += 1
					continue
				m_seq[i+iterateur][j][0] += occ[pbCurr][aa][7+iterateur]
				m_seq[i][j][1] += 1
				iterateur += 1
			iterateur = 0
			while iterateur < 7 and (i-iterateur) > 0:
				m_seq[i-iterateur-1][j][0] += occ[pbCurr][aa][7-iterateur-1]
				m_seq[i][j][1] += 1
				iterateur += 1
	return m_seq

def moyenne(m_seq):
	"""Moyenne la somme faite sur la matrice d'occurence et renvoie une matrice 
		-Args:
			_m_seq : matrice derivee de la matrice d'occurence construite dans matrice_sequence
	"""
	for i in range (len(m_seq)):
		for j in range(20):
			m_seq[i][j] = m_seq[i][j][0] / m_seq[i][j][1]
	return m_seq

def retourneSequence(mseq,naa,Rand):
	"""Retourne une sequence qui correspond aux valeurs les plus elevees dans mseq
		-Args:
			_mseq : matrice derivee de la matrice d'occurence construite dans matrice_sequence et modifiee dans moyenne
			_naa : taille de la sequence d'origine
	"""
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	newseq = ''
	npb = len(mseq)
	for i in range (npb-naa,npb): # doit faire une taille de naa. On saute les premières positions au cas ou
		#k = mseq[i].index(max(mseq[i]))# l'index des acides aminés de mseq (mseq[x]) correspond a l'index de AA
		k = random.choice(np.argsort(mseq[i])[19-Rand:])#[:0+Rand+1])#[19-Rand:]) 
		newseq += AA[k]
	return newseq

#==============================================================================================================================#

def mutagenese(seq,n):
	"""Prend une sequence en acide aminé et renvoie une sequence mutée sur n positions aléatoires
		-Args:
			_seq : la sequence en acide aminé sous forme de chaine de char
			_n : le nombre de mutations désirée 
	Pas de mutations identiques, et les mutations se font sur n positions distinctes
	"""
	intervalle = list(range(len(seq)))
	positions = []
	for i in range(n):
		rd = random.randint(0,len(intervalle)-1)
		positions.append(intervalle[rd])
		del intervalle[rd]
	
	seqMut = seq
	for pos in positions:
		AA = 'ACDEFGHIKLMNPQRSTVWY'
		aa = seq[pos]
		if AA.index(aa) == 19:	
			exclu = AA[:AA.index(aa)]
		else:
			exclu = AA[:AA.index(aa)] + AA[AA.index(aa)+1:]
		rd = random.randint(0,len(exclu)-1)
		mut = exclu[rd]
		if pos != len(seq)-1:
			seqMut = seqMut[:pos] + mut + seqMut[pos+1:]
		else:
			seqMut = seqMut[:pos] + mut
	return seqMut

def computeFoRSA(seqaa,seqpb,gapP):
	"""Compute FoRSA global on two sequences (PB and aa)
		-Args:
			_seqaa : Fichier de sequence d'acide aminé
			_seqpb : Fichier de sequence de Protein Blocs
			_gapP : gap penalty (int)
	Sauve l'info dans un fichier A MODIFIER POUR NE PAS ECRIRE PUIS LIRE LE FICHIER
	"""
	cmd = 'forsa_global '+seqaa+ ' '+ seqpb+' '+gapP+' >'+' '+ 'ephScore.txt'
	subprocess.run( ['/bin/bash', '-i', '-c', cmd])#.decode('utf-8')
	Zscore = ReadFoRSA('ephScore.txt')
	return Zscore
	# bashOutput = bashOutput.decode('utf-8')
	# memory = [] # contient en 0 la sequence aa, en 1 la sequence PB et en 2 une liste [rawscore,score normal, longueur d alignement, Z-score]
	# memory.append(bashOutput.split('\n')[1])
	# memory.append(bashOutput.split('\n')[2])
	# memory.append(bashOutput.split('\n')[3].split(':')[1:])
	# return memory

#==============================================================================================================================#
#  MONTE CARLO

def transforme(liste):
	"""Transforme une liste de chiffre en liste de probabilité correspondante
		-Args:
			_liste : une liste de double 
	"""
	m = min(liste)
	if m < 0:
		m = -m
		for i in range(len(liste)):
			liste[i] += m
	somme = 0
	for i in range(len(liste)):
		somme += liste[i]
	for i in range(len(liste)):
		liste[i] = liste[i] / somme
	return liste

def generateurStochastique(m_seq,naa):
	"""Genère une sequence aléatoire sur la base de la matrice de somme d'occurence spécifique de sequence
		-Args:
			_mseq : matrice derivee de la matrice d'occurence construite dans matrice_sequence et modifiee dans moyenne
			_naa : taille de la sequence d'origine
	Les sommes d'occurences sont transformées en poids/probabilité puis on determine un aa par position selon cette distribution
	"""
	repartition = [0.4,0.3,0.15,0.12,0.001875,0.001875,0.001875,0.001875,0.001875,0.001875,0.001875,0.001875,0.001875,0.001875,0.001875,0.001875,0.001875,0.001875,0.001875,0.0018750000000006262]
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	
	seq = ''
	for i in range(naa):
		Arg = np.argsort(m_seq[i]) # trier les indices de m_seq par ordre croissant de valeur
		AaTrie = [] # va contenir les acides aminés par ordre décroissant de valeur
		for i in range(len(Arg)-1,-1,-1):
			AaTrie.append(AA[Arg[i]])
		choix = random.choices(AaTrie,repartition)
		seq += choix[0]
	return seq

def probabilite_mseq(colonne):
	"""Prend en arg une colonne de la mseq et renvoie une liste de probabilité associée à chaque valeur
	"""
	# ecartM = vAbs(max(colonne), min(colonne)) # ecart maximum entre deux valeurs de la colonne

	Somme = sum(colonne)
	Probas = []
	for j in range(len(colonne)):
		Probas.append(colonne[j] / Somme)
	return Probas

def generateurStochastique2(m_seq,naa):
	"""Genère une sequence aléatoire sur la base de la matrice de somme d'occurence spécifique de sequence
		-Args:
			_mseq : matrice derivee de la matrice d'occurence construite dans matrice_sequence et modifiee dans moyenne
			_naa : taille de la sequence d'origine
	Les sommes d'occurences sont transformées en poids/probabilité puis on determine un aa par position selon cette distribution
	"""
	repartition = []
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	
	seq = ''
	for i in range(naa):
		probas = probabilite_mseq(m_seq[i])
		# print(probas)
		choix = random.choices(AA,probas)
		seq += choix[0]
	return seq

def muteSelonMatrice(mseq,posi,aapred):
	"""Fonctionne comme generateurStochastique mais pour une position
		-Args:
			_mseq : matrice derivee de la matrice d'occurence construite dans matrice_sequence et modifiee dans moyenne
			_posi : position sur la sequence de la mutation a generer
			_aapred : acide amine a la posi avant la mutation
	"""
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	nliste = [elem for i,elem in enumerate(mseq[posi]) if i != AA.index(aapred)]
	AA = [elem for elem in AA if elem != aapred]
	lproba = transforme(nliste)
	choix = random.choices(AA,lproba)
	return choix[0]

def mutOpti(seq,sPred,pbs):
	"""Mute la sequence de maniere a ce que la mutation augmente le scoe au maximum
		-Args :
			_seq: sequence en aa (chaine de caractères) 
			_sPred: Zscore de la sequence
			_pbs: sequence de protein blocks en chaine de char
	"""
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	maximum = -100
	memory = []
	g = 0
	for pos in range(len(seq)):
		memory.append([])
		for aa in AA:
			cseq = list(seq)
			if aa == cseq[pos]: # on ne mute pas par le meme aa
				memory[pos].append(0)
				continue
			cseq[pos] = aa
			cseq = "".join(cseq)
			score = fGpy.launchFoRSAfromPy(cseq,pbs,-5)
			diff = score - sPred
			memory[pos].append(diff)
			if diff > maximum:
				maximum = diff
				keep = [aa,pos,score]
	writeDataMutation(memory)
	newseq = list(seq)
	newseq[keep[1]] = keep[0]
	newseq = "".join(newseq)
	print("Sequences pred&mut :\n{} Score: {:.3f}\n{} Score: {:.3f}".format(seq,sPred,newseq,keep[2]))
	return [newseq,keep[2]]

def mutOpti2(seq,sPred,pbs):
	"""Mute la sequence de maniere a ce que la mutation augmente le scoe au maximum
		-Args :
			_seq: sequence en aa (chaine de caractères) 
			_sPred: Zscore de la sequence
			_pbs: sequence de protein blocks en chaine de char
	"""
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	maximum = -100
	memory = []
	order = [] # permet de conserver les meilleures mutation a chaque position
	g = 0
	for pos in range(len(seq)):
		memory.append([])
		maxLocal = -100
		for aa in AA:
			cseq = list(seq)
			if aa == cseq[pos]: # on ne mute pas par le meme aa
				memory[pos].append(0)
				continue
			cseq[pos] = aa
			cseq = "".join(cseq)
			score = fGpy.launchFoRSAfromPy(cseq,pbs,-5)
			diff = score - sPred
			memory[pos].append(diff)
			if diff > maxLocal:
				maxLocal = diff
				keepLocal = [aa,pos,score]
			if diff > maximum:
				maximum = diff
				keep = [aa,pos,score]
		order.append(keepLocal)
	ordArgS = np.argsort(order,axis=0)



	# print('Selection globale : {}\nSelection locale : {}\nSelection suivante : {}'.format(keep,order[ordArgS[-1][2]],order[ordArgS[-2][2]]))
	lastTokeep = [keep] # permet de muter sur les positions les plus favorables a conditions qu'elles soient à + de 15 positions d'ecarts
	for i in range(len(ordArgS)-1,-1,-1):
		flag = 1
		for j in range(len(lastTokeep)):
			if vAbs(order[ordArgS[i][2]][1], lastTokeep[j][1]) < 14:
				flag = 0
				break
		if flag:
			# print(i,order[ordArgS[i][2]][1],lastTokeep[j][1])
			lastTokeep.append(order[ordArgS[i][2]])
	# sys.exit()
	# print(lastTokeep)

	writeDataMutation(memory)
	newseq = list(seq)
	for mutation in lastTokeep:
		newseq[mutation[1]] = mutation[0]
	newseq = "".join(newseq)
	print("Sequences pred&mut :\n{} Score: {}\n{} Score: {}".format(seq,sPred,newseq,keep[2]))
	return [newseq,keep[2]]

def mutScc(pos,seq,sPred,pbs,return_mutIn):
	"""Mute une position par chaque autre acide amine possible et calcule un FoRSA score 
	Garde dans la memoire partagee return_mutIn le meilleur score obtenu
		-Args :
			_pos: position dans la sequence
			_seq: sequence en aa (chaine de caractères) 
			_sPred: Zscore de la sequence
			_pbs: sequence de protein blocks en chaine de char
			_return_mutIn: memoire partagée (permet de recuperer les donnee issues de pls processus)
	"""
	AA = 'ACDEFGHIKLMNPQRSTVWY'
	memory = {}
	maxLocal = -100
	for aa in AA:
		cseq = list(seq)
		if aa == cseq[pos]: # on ne mute pas par le meme aa
			continue
		cseq[pos] = aa
		cseq = "".join(cseq)
		score = fGpy.launchFoRSAfromPy(cseq,pbs,-5)
		diff = score - sPred
		mut = '{}{}{}'.format(seq[pos],pos,aa)
		memory[mut] = diff
		if diff > maxLocal:
			maxLocal = diff
			keepLocal = [aa,pos,diff,score]
	return_mutIn[pos] = (keepLocal,memory)


def mutOptiParallel(seq,sPred,pbs,idSeq):
	"""Mute la sequence de maniere a ce que la mutation augmente le scoe au maximum
		-Args :
			_seq: sequence en aa (chaine de caractères) 
			_sPred: Zscore de la sequence
			_pbs: sequence de protein blocks en chaine de char
	Le processus est parallélisé sur le nombre de positions cad pour chaque position sur la sequence on lance 
	mutScc qui va muter par chaque autre acide amine et calculer un score FoRSA.
	"""
	maximum = -100
	memory = []
	g = 0
	manager2 = multiprocessing.Manager()
	return_mutIn = manager2.dict()
	Sprocesses = []
	for pos in range(len(seq)): # un processus par position
		Sprocessus = multiprocessing.Process(target=mutScc,args=[pos,seq,sPred,pbs,return_mutIn])
		Sprocessus.start()
		Sprocesses.append(Sprocessus)

	for Sprocess in Sprocesses:
		Sprocess.join()

	order = [] # permet de conserver les meilleures mutation a chaque position
	for pos in return_mutIn.keys():
		order.append(return_mutIn[pos][0]) # ajouter à order [aa,pos,score] 
		memory.append(return_mutIn[pos][1])
	ordArgS = np.argsort(order,axis=0) # ordArgs : liste des indicies qui trieraient order par valeur de Zscore (et trie les autres aussi)
	lastTokeep = []
	if order[ordArgS[-1][-1]][-2] > 0:
		lastTokeep = [order[ordArgS[-1][-1]]] # permet de muter sur les positions les plus favorables a conditions qu'elles soient à + de 15 positions d'ecarts
		# la première position de lastTokeep est la liste [aa,pos,score] qui comporte la valeur max de score dans order soit order[dernier indice de ordArgs]

	for i in range(len(ordArgS)-2,-1,-1): # en partant de la fin, i parcours ordArgS (du maximum au minimum) pour prioriser les scores max en premiers
		flag = 1
		if order[ordArgS[i][-1]][-2] < 0: # si la valeur de Z score est en dessous de 0 alors on s'arrête, ces mut ne seront pas séléctionnées
			break
		# continue
		for j in range(len(lastTokeep)): # j parcours l'ensemble des [aa,pos,score] conservées pour la mutation
		# noter qu'au départ lastToKeep ne contient que la valeur maximum
		# ainsi on teste pour chaque mutation dans l'ordre décroissant de Z score si sa position est à une distance inférieure à 7,
		# si oui alors on n'inclus pas  la mutation dans la séquence
			if vAbs(order[ordArgS[i][2]][1], lastTokeep[j][1]) < 7: 
				flag = 0
				break
		if flag: # uniquement si la mutation est à distance de toutes les autres (précédentes)
			lastTokeep.append(order[ordArgS[i][-1]]) 

	writeDiffScoreMute(memory,idSeq)
	# writeDataMutation(memory)
	newseq = list(seq)
	for mutation in lastTokeep:
		newseq[mutation[1]] = mutation[0]
	newseq = "".join(newseq)
	if lastTokeep:
		print("Sequences {} pred&mut :\n{} Score: {:.4f}\n{} Score: {:.4f}".format(idSeq,seq,sPred,newseq,lastTokeep[0][-1]))
		return [newseq,lastTokeep[0][-1]]
	else: 
		print("Sequences {} reached Ostate :\n{} Score: {:.4f}\n{} Score: {:.4f}".format(idSeq,seq,sPred,newseq,sPred))
		return [newseq,sPred]

def SSM2(seq,sPred,pbs,N,return_dict,keep_Scores):
	"""Site saturation mutagenesis done multiple time on one sequence
		-Args :
			_seq: sequence en aa (chaine de caractères) 
			_sPred: Zscore de la sequence
			_pbs: sequence de protein blocks en chaine de char
			_N: number of iteration on sequence
			_return_dict: permet d'acceder a la valeur retournée par le processus
	"""
	for _ in range(N):
		seqM = mutOptiParallel(seq,sPred,pbs)
		seq = seqM[0]
		sPred = seqM[1]
		keep_Scores[seq].append(seqM[1])
	return_dict[seq] = seqM

def SSM(seq,sPred,pbs,N,return_dict,keep_Scores):
	"""Site saturation mutagenesis done multiple time on one sequence
		-Args :
			_seq: sequence en aa (chaine de caractères) 
			_sPred: Zscore de la sequence
			_pbs: sequence de protein blocks en chaine de char
			_N: number of iteration on sequence
			_return_dict: permet d'acceder a la valeur retournée par le processus
	"""
	L = keep_Scores[seq]
	seqPred = seq
	finish = False
	for _ in range(N):
		seqM = mutOptiParallel(seq,sPred,pbs,return_dict[seqPred][-1])
		if is_equal(seqM[1],sPred):
			finish = True
		seq = seqM[0]
		sPred = seqM[1]
		L.append(seqM[1])
		if finish:
			break

	keep_Scores[seqPred] = copy.deepcopy(L)
	seqM.append(return_dict[seqPred][-1])
	return_dict[seqPred] = seqM

#==============================================================================================================================#

if __name__ == '__main__':
	Randomiseur = 0

	if len(sys.argv) == 2 and sys.argv[1] == 'help':
		print(__doc__)
		sys.exit()
	occ = sys.argv[1]
	pbs = sys.argv[2]
	aa = sys.argv[3]

	P = lit_seq_PB(pbs)
	M = lit_Matrice(occ)
	A = lit_seq_AA(aa)

	if not os.path.isdir('./out'):
		subprocess.run(['mkdir','out']) # crée un répertoire out pour stocker les fichiers de sortie

########################################################################## SOMME
	# for key in P.keys():
	# 	a = A[key]
	# 	p = P[key]
	# 	m = moyenne(matrice_sequence(M,p))
	# 	nseq = retourneSequence(m,len(a),Randomiseur)
	# 	print(nseq)
	# 	mut = muteSelonMatrice(m,0,nseq[0])
	# 	nseq = mut + nseq[1:]
	# 	print(nseq)
	# 	writeSeq(nseq)
	# 	spred = computeFoRSA('newseq.aaseq',pbs,'-5')
	# 	 
	# 	print("le score de FoRSA initial est {}".format(spred))
##########################################################################

########################################################################## MONTE CARLO
	for key in P.keys(): # Ne tourne que pour une sequence 
		a = A[key]
		p = P[key]
		m = moyenne(matrice_sequence(M,p))

		# Score de reference
		nseq = retourneSequence(m,len(a),Randomiseur)
		# writeSeq(nseq)
		# ZscoreI = computeFoRSA('newseq.aaseq',pbs,'-5')
		ZscoreI = fGpy.launchFoRSAfromPy(nseq,p,-5)
		print("Le score initial de la sequence {} est : {:.3f}".format(nseq,ZscoreI))

		# Generer la population de sequences
		nSeqPop = 0 
		manager = multiprocessing.Manager()
		return_dict = manager.dict()
		return_dict[nseq] = [nseq,ZscoreI,'ID0']
		keep_Scores = manager.dict()
		keep_Scores[nseq] = [ZscoreI]
		while nSeqPop < 4:
			nseq = generateurStochastique(m,len(a))
			Zscore = fGpy.launchFoRSAfromPy(nseq,p,-5)
			if Zscore > ZscoreI*0.65:
				if nseq not in return_dict.keys():
					return_dict[nseq] = [nseq,Zscore,'ID{}'.format(nSeqPop+1)]
					keep_Scores[nseq] = [Zscore]
					nSeqPop += 1
					print("On ajoute une sequence {} de score {:.3f} à la population...".format(nseq,Zscore))
		print('Nombre de séquence dans la population : {}. Début de l\'optimisation...'.format(nSeqPop))
		# Mutagenese complete SSM

		start = time.perf_counter()
		processes = []
		N = 20
		for elem in return_dict.keys():
			processus = multiprocessing.Process(target=SSM,args=[return_dict[elem][0],return_dict[elem][1],p,N,return_dict,keep_Scores])
			processus.start()
			processes.append(processus)

		for process in processes:
			process.join()

		finish = time.perf_counter()
		print("TIME start to finish = {:.2f}".format(finish-start))

		# Sortie

		writeEvolScores(keep_Scores) # ecrit l'evolution des scores pour chaque séquence

		if return_dict:
			writeFileOfNewSeq(return_dict,nSeqPop) # nouvelles sequences obtenues dans un fichier
		else:
			print(return_dict)
##########################################################################

########################################################################## GRAINE
	# for key in P.keys():
	# 	a = A[key]
	# 	p = P[key]
	# 	seqs = graine(p,M,3)
	# for seq in seqs:
	# 	writeSeq(seq)
	# 	memory = computeFoRSA('newseq.aaseq',pbs,'-5')
	# 	spred = memory[2][3]
	# 	print("Le score de FoRSA pour cette sequence est {}".format(spred))
########################################################################## 

########################################################################## MUTAGENESE ET AMELIORATION DU SCORE
	# for i in range(1000):
	# 	seqm = mutagenese(nseq,4)
	# 	writeSeq(seqm)
	# 	memory = computeFoRSA('newseq.aaseq',pbs,'-5')
	# 	if memory[2][3] > spred:
	# 		spred  = memory[2][3]
	# 		print("le score de FoRSA a été amélioré, il vaut : {}".format(spred))
	# 		print('La sequence correspondante est : {}'.format(seqm))
	# 		nseq = seqm
##########################################################################
