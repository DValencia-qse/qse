


import numpy as np
import wolframalpha as wa
import sympy as sp
from dwave.system import DWaveSampler,EmbeddingComposite
import dwave.inspector
import sys
numero=15
if len(sys.argv)==2:
	numero=int(sys.argv[1])

tam_auxiliares=0
binario=bin(numero)
tam=len(bin(numero))-2
qbits=np.floor(np.log2(numero)**2/4).astype(dtype='int')
tam_p=np.ceil(qbits/2).astype(dtype='int')
tam_q=qbits-tam_p

MM=np.zeros(shape=(tam+2,tam),dtype='U16')
MM[:,:]="0"
MM[0,0]="PP"
MM[0,tam-1]="1"
MM[1,0]="QQ"
MM[1,tam-1]="1"
for i in range(0,tam):
	MM[tam+1,i]=binario[i+2]
cadSimbolos=''
for i in range(1,tam_p+1):
	MM[0,tam-i-1]="p_"+str(i)
	cadSimbolos=cadSimbolos+' '+"p_"+str(i)
for i in range(1,tam_q+1):
	MM[1,tam-i-1]="q_"+str(i)
	cadSimbolos=cadSimbolos+' '+"q_"+str(i)

for i in range(0,tam_q+1):
	for j in range(0,tam_p+1):
		MM[2+i,tam-j-i-1]=MM[1,tam-i-1]+"*"+MM[0,tam-j-1]

print(MM)

ecuaciones=np.zeros(shape=(tam),dtype='U64')
for i in range(0,tam):
	for j in range(0,tam):
		if (i==tam-1):
			ecuaciones[j]=ecuaciones[j]+"-"+MM[2+i,j]
		else:
			ecuaciones[j]=ecuaciones[j]+MM[2+i,j]+"+"

print(ecuaciones)
print(cadSimbolos)
listaSimb=sp.symbols(cadSimbolos)
exp=""

consulta="expand: "
for i in  range(len(ecuaciones)-1):
	consulta=consulta+"("+ecuaciones[i]+")^2+"
	exp=exp+"("+ecuaciones[i]+")**2+"

consulta=consulta+"("+ecuaciones[len(ecuaciones)-1]+")^2"
exp=exp+"("+ecuaciones[len(ecuaciones)-1]+")**2"
#print(consulta)
print(exp)
resultado=sp.expand(exp)
print(resultado)
print(listaSimb)
print(len(listaSimb))
for i in range(0,len(listaSimb)):
	resultado=resultado.subs(listaSimb[i]**2,listaSimb[i])

#print(resultado)
auxiliares=[]
contAux=0
resultadoFinal=resultado
for t in resultado.args:
#	print(t)
	suma=0
	for l in range(0,len(listaSimb)):
		suma=suma+t.count(listaSimb[l])	
	if (suma==3):
		contAux=contAux+1 
		auxi=sp.Symbol('a_'+str(contAux))
		auxiliares.append(auxi)
#	Esto habra que cambiarrlo
		mm4=sp.Mul(3,auxi)
		mm3=sp.Mul(-2,sp.Mul(t.args[len(t.args)-2],auxi))
		mm2=sp.Mul(-2,sp.Mul(t.args[len(t.args)-3],auxi))
		mm1=sp.Mul(t.args[len(t.args)-3],t.args[len(t.args)-2])
		mm=sp.Mul(2,sp.Add(sp.Add(mm1,mm2),sp.Add(mm3,mm4)))
		#Aqui hay un error
		mm=sp.Add(sp.Mul(auxi,t.args[len(t.args)-1]))
		if (len(t.args)==4):
			mm=sp.Mul(t.args[0],mm)

		#print(mm)
		#input()
			
		resultadoFinal=resultadoFinal.subs(t,mm)
#print(resultado.type())
print("Sin mul*2")
print(resultadoFinal)
#Ahora sobre el resultado final hago una sustitucion
resultadoFinal=sp.Mul(2,resultadoFinal)
listaSimbolosFinal=list(listaSimb)
for a in auxiliares:
	listaSimbolosFinal.append(a)
listaSimb=tuple(listaSimbolosFinal)
#print(listaSimb)
sustitutos=[]
contS=0
for t in listaSimb:
	contS=contS+1
	s=sp.Symbol('s_'+str(contS))
	sustitutos.append(s)
#	resultadoFinal.subs(t,(
#print(sustitutos)
contS=0
print("Antes de sustitucion")
print(resultadoFinal)
for t in listaSimb:
	resultadoFinal=resultadoFinal.subs(t,sp.Add(1,-sustitutos[contS]))
	contS=contS+1

print("Antes de expansion")
print(resultadoFinal)
resultadoFinal=sp.expand(resultadoFinal)
print("Tras expansion")
print(resultadoFinal)

#Ahoraa se crea el diccionario
h={}
J={}
for t in resultadoFinal.args:
	suma=0
	for s in sustitutos:
		suma=suma+t.count(s)
	if (suma>0):
		print(t)
		if (suma==1):
			h[str(t.args[1])]=t.args[0]
		else:
			J[(str(t.args[1]),str(t.args[2]))]=t.args[0]

for s in sustitutos:
	if not (str(s) in h):
		h[str(s)]=0


print(h)
print(J)

sampler=EmbeddingComposite(DWaveSampler())

sampleset=sampler.sample_ising(h,J,num_reads=100)
for sample in sampleset.samples(sorted_by='energy'):
	print(sample)
#dwave.inspector.show(sampleset)
result=sampleset.samples(sorted_by='energy')[0]
print(result)
#print(sustitutos)
cadenas=list(map(str,sustitutos))
res_p='0b'+'0'*tam_p
res_q='0b'+'0'*tam_q
res_p=list(res_p)
res_p[len(res_p)-1]='1'
res_q=list(res_q)
res_q[len(res_q)-1]='1'

for k in result:
	pos=cadenas.index(k)
#	print(pos)
#	print(str(listaSimb[pos])+"==>"+str(int((1-result[k])/2)))
	l=str(listaSimb[pos]).split('_')
	if l.count('q')>0:
		res_q[len(res_q)-1-int(l[1])]=str(int((1-result[k])/2))
	if l.count('p')>0:
		res_p[len(res_p)-1-int(l[1])]=str(int((1-result[k])/2))
print(listaSimb)
print(sustitutos)
#print(res_q)
#print(res_p)
print(int(''.join(res_q),2))
print(int(''.join(res_p),2))

"""
cliente=wa.Client('UX5V4Y-7X5VXU8RQK')
result=cliente.query(consulta)
for pod in result.results:
	print(pod.text)
	cadenaResultado=pod.text

#	for sub in pod.subpods:
#		print(sub.plaintext)
simplificados=cadenaResultado.replace('^2','')
print(simplificados)
res2=cliente.query("expand: "+simplificados)
for pod in res2.results:
	print(pod.text)
	reducidos=pod.text
recortados=reducidos.split('+')
final=[]
for t in recortados:
	if (t.count("-")==0):
		final.append(t)
	else:
		menos=t.split("-")
		primero=0
		for m in menos:
			if (primero==0):
				final.append(m)
				primero=1
			else:
				final.append("-"+m)


print(recortados)
print(final)
#Comprobar si hay que introducir variables auxiliares
print("Comprobar si hay que introuir variables auxiliares")
definitiva=[]
for t in final:
	if ((t.count('q')+t.count('p'))==3):
		print(t.split())
		valores=t.split()
		nuevo=[]
		tam_auxiliares=tam_auxiliares+1
		if (len(valores)>3):
			nuevo.append(valores[0])
		nuevo.append("(")
#		Hay que ver si es positiva o negativa
		nuevo.append("a_"+str(tam_auxiliares)+valores[len(valores)-1])
		nuevo.append("+2(3a_"+str(tam_auxiliares)+"+")
		nuevo.append("-2a_"+str(tam_auxiliares)+valores[len(valores)-2])
		nuevo.append("-2a_"+str(tam_auxiliares)+valores[len(valores)-3])
		nuevo.append("+"+valores[len(valores)-3]+valores[len(valores)-2]+")")
		nuevo.append(")")
		definitiva.append("".join(nuevo))
		
	else:
		definitiva.append(t)


print(definitiva)
consultaDef=definitiva[0]
for s in definitiva[1:]:
	consultaDef=consultaDef+" + "+s
print(consultaDef)
res3=cliente.query("expand: "+consultaDef)
for pod in res3.results:
	#print(pod.text)
	hyJ=pod.text

#positivos=cadenaResultado.split('+')
#print(positivos)
h={}
J={}
finalisimo=[]
rec=hyJ.split('+')
for r in rec:
	if (r.count('-')==0):
		finalisimo.append(r)
	else:
		menos=r.split("-")
		primero=0
		for m in menos:
			if (primero==0):
				finalisimo.append(m)
				primero=1
			else:
				finalisimo.append("-"+m)

print(finalisimo)
#Aqui ha que realizar la conversion de x1 a (1-s1)
dictConversion={}
contadorS=0
for d in finalisimo:
	if ((d.count('q')+d.count('p')+d.count('a'))==1):
		contadorS=contadorS+1
		datos=d.split()
		dictConversion[datos[1]]='s_'+str(contadorS)

print("=====Cambio de variables=====")
print(dictConversion)
print("=====Fin cambio variables=====")
#Aactualizacion de variables
consultaConversion=[]
for s in finalisimo:
	n=s
	for k in dictConversion:
#		print(k+"->"+dictConversion[k])
		n=n.replace(k,"(1 - "+dictConversion[k]+")")
	consultaConversion.append(n)
inicio=1
if (consultaConversion[0]==''):
	cc="2( "+consultaConversion[1]
	inicio=2
else:
	cc="2( "+consultaConversion[0]
for s in consultaConversion[inicio:]:
	cc=cc+" + "+s
cc=cc+" )"
print(cc)
res4=cliente.query("expand: "+cc)
hyJ2=''
for pod in res4.results:
	print(pod.text)
	hyJ2=pod.text
print(hyJ2)

for d in finalisimo:
#	print(d)
	if ((d.count('q')+d.count('p')+d.count('a'))==1):
		print(d)
		datos=d.split()
		print(datos)
		if len(datos)==2:
			if (len(datos[0])==1) and (datos[0].count('-')==1):
				h[datos[1]]=-1
			else:
				h[datos[1]]=int(datos[0])
		else:
			h[datos[0]]=1

print(h)
for d in finalisimo:
#	print(d)
	if ((d.count('q')+d.count('p')+d.count('a'))==2):
		print(d)
		datos=d.split()
		print(datos)
		if (len(datos)==3):
			if (len(datos[0])==1) and (datos[0].count('-')==1):
				J[(datos[1],datos[2])]=-1
			else:
				J[(datos[1],datos[2])]=int(datos[0])
		else:
			if len(datos)==2:
				J[(datos[0],datos[1])]=1
			else:
				if (len(datos[0])==1) and (datos[0].count('-')==1):
					J[(datos[2],datos[3])]=int(datos[0]+datos[1])
				
				

print(J)
"""
#		if len(datos)==2:
#				h[datos[1]]=-1
#			else:
#				h[datos[1]]=int(datos[0])
#		else:
#			h[datos[0]]=1


#h={'s1': 580, 's2': 420, 's3': 144, 's4': 128}
#J={('s1','s2'): 152, ('s1','s3'): -144, ('s1','s4'): -512, ('s2','s3'): 16,('s2','s4'):  -512, ('s3','s4'): 128}
