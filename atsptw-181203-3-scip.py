#!/usr/bin/env python
from __future__ import print_function
#from google.apputils import app
from ortools.linear_solver import pywraplp
import json, sys
#from pulp import LpProblem, LpMinimize, LpVariable, LpStatus, lpSum
from subprocess import call
import os
import csv
import numpy
import random
import os.path
import goto
from goto import goto, comefrom, label
from math import exp
import math

# mememp (0=Def/Normal) (1): MemoryEmphasis conserve memory where possible
# NodeFile: (1) Default on Memory (2) Node file on disk
# NodeSel: (1) Best-bound search default, (2) Best-estimate search  
# MIPEmphasis:(1) Emphasize feasibility over optimality
# ParallelMode: (0) Automatic/default: CPLEX decide deterministic or opportunistic, (-1) Enable opportunistic parallel search mode
# VarSel:(0) Automatic/default CPLEX choose variable to branch, (3) Strong branching

################################################################################################################################################################
#### PARAMETERS  (Hernan = 18 Rutas)
################################################################################################################################################################
ksize    = 87 + 1   # Nodo DEPOT
xtol     = 20  		# 15
ktour    = 13	   	# 29,12		# First Route lenght attempt 
################################################################################################################################################################
klastrun  = -6					# Use klastrun  = 0 to start from scratch. Use Negative to start from any given run number.
xstartlen = 0					# Use xstartlen = 0 to start from scratch. 
xtypesol = ""					# Use "" for regular / "opt" (pre-solved) 
wantopt  = True
xquit    = True				
################################################################################################################################################################
runopt   = False
xendopt  = False
xhardend = False
xthold = 60
bigm = 1000
kmaxlate = 100
worsttw = 0
worstrun = -99
kruns = 20
xrouteinc = 0
xcom = ""
################################################################################################################################################################
#### STORES DEFINITION
################################################################################################################################################################
kstore0 = 0
kstore1 = 85
kstore2 = 86
kstore3 = 87

xstores = [False for i in range(ksize)]
xstores[0]  = True
xstores[85] = True
xstores[86] = True
xstores[87] = True

################################################################################################################################################################
#### LOAD DISTANCE/COST TABLE
################################################################################################################################################################
xgeo   = numpy.zeros(ksize)
ygeo   = numpy.zeros(ksize)
xdep   = numpy.zeros(ksize, dtype = int)
twa   = [0 for i in range(ksize)] 
twb   = [0 for i in range(ksize)] 
tserv = [0 for i in range(ksize)] 

infile = open("heb-instance.csv", "r")
mybuf = csv.reader(infile, delimiter=',')
next(mybuf)
kp = 0

for myrow in mybuf:
	kp		  = int(myrow[0])
	xdep[kp]  = int(myrow[1])
	xgeo[kp]  = float(myrow[2])
	ygeo[kp]  = float(myrow[3])
	twa[kp]   = int(myrow[4])
	twb[kp]   = int(myrow[5])
	tserv[kp] = int(myrow[6])
	if xstores[kp] :	
		twa[kp]   = 0
		twb[kp]   = bigm
	kp = kp + 1

tserv[0] = 0
infile.close()

################################################################################################################################################################
#### GENERATE TIMEs / COST TABLE
################################################################################################################################################################
cost = [[0 for j in range(ksize)] for i in range(ksize)]
for i in range(ksize):
	for j in range(ksize):
		if i <> j :
			lon_1, lat_1, lon_2, lat_2 = map(math.radians, [xgeo[i], ygeo[i], xgeo[j], ygeo[j]])
			Delta_lat = lat_2 - lat_1
			Delta_lon = lon_2 - lon_1
			a = (math.sin(Delta_lat / 2))**2
			c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
			lat_distance = 3959 * c
			a = (math.sin(Delta_lon / 2))**2
			c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
			lon_distance = 3959 * c
			#xdist = round(150*((xgeo[i] - xgeo[j])**2 + (ygeo[i] - ygeo[j])**2)**0.5 ,0)
			xdist = abs(lat_distance) + abs(lon_distance)		# Distance miles		
			cost[i][j] = xdist  * 1.81 							# Travel_Time

################################################################################################################################################################
#### INICIALIZAR ESTADO DE LOS NODOS Y VECTORES DE CAMBIOS DE ORDENES (ADD / DELETES)
################################################################################################################################################################
xsol  = [[[False for j in range(ksize)] for i in range(ksize)] for k in range(kruns)]	
xlate = [[0 for i in range(ksize)] for k in range(kruns)]	
xchaddord = [999 for i in range(ksize)] 
xchdelord = [999 for i in range(ksize)] 
xcntord   = [0 for k in range(kruns)]  

################################################################################################################################################################
#### INICIA BUCLE DE OPTIMIZACION
################################################################################################################################################################
label .mystart

################################################################################################################################################################
#### INICIALIZAR ESTADO DE LOS NODOS
################################################################################################################################################################
xactive = [True for i in range(ksize)]
xdif = 0

################################################################################################################################################################
#### PRE-CARGA DE LAS CORRIDAS
################################################################################################################################################################
### Use klastrun = 0 to start from scratch
################################################################################################################################################################
if klastrun <> 0 :
	tstamp = [[0 for i in range(ksize)] for k in range(kruns)]
	xsol  = [[[False for j in range(ksize)] for i in range(ksize)] for k in range(kruns)]	
	xlate = [[0 for i in range(ksize)] for k in range(kruns)]	
	maxorders = 0
	worsttw = 0
	worstord = 0

	for krun in range(abs(klastrun)+1):
		kfilesol = "sol" + xtypesol + "-" + str(krun) + ".sol"	
		if os.path.isfile(kfilesol) : 
			if krun == 0 :
				initime = os.path.getctime(kfilesol) 
			infile = open(kfilesol, "r")
			mybuf = csv.reader(infile, delimiter='\t')    
			klastour = 0
			qtyst = 0

			if len(list(mybuf)) <= 2 :
				print("General Error. No file found...")
				quit()
			infile.seek(0)		# Incia en la 2da linea 
			next(mybuf)
			next(mybuf)			
			
			for myrow in mybuf:
				wvariable = myrow[0][0:25].strip()
				wvalue	  = float(myrow[0][26:])	
				if wvariable[0:4] == "pend" and wvalue >= 0.99 :			
					kind = int(wvariable[5:wvariable.find(']')])
					klastour += 1
					ktour = klastour
					if not xstores[kind] :  
						xactive[kind]= True
						xdif += 1
					else :
						qtyst += 1
				if wvariable[0:4] == "arcs" and wvalue >= 0.99 :
					kind1 = int(wvariable[5:wvariable.find(']')])
					kcoma = wvariable.find(',')
					kind2 = int(wvariable[kcoma+2:wvariable.find(']',kcoma)])
					xsol[krun][kind1][kind2] = True
					if not xstores[kind1] :
						xactive[kind1] = False	
					if not xstores[kind2] :
						xactive[kind2] = False	
				if wvariable[0:5] == "arrhr" and wvalue > 0 :
					kind = int(wvariable[6:wvariable.find(']')])
					tstamp[krun][kind] = wvalue
				if wvariable[0:4] == "late" and wvalue >= 1 :
					kind = int(wvariable[5:wvariable.find(']')])
					xlate[krun][kind] = wvalue
					print("Late: ", kind," = ",round(wvalue,0))
					if xlate[krun][kind] > worsttw :
						worsttw  = xlate[krun][kind]
						worstord = kind
						worstrun = krun
						
		infile.close()
		print("**************** Route Built #: ", krun, " = ", klastour - qtyst)
		if klastour - qtyst > maxorders :
			maxorders = klastour - qtyst 

	endtime = os.path.getctime(kfilesol) 
	
	#### OPTIMIZACION PARA INCLUIR CLIENTES ESPECIFICOS EN LA SIGUIENTE RUTA A CONSTRUIR
	if runopt and worstrun <> -99 and not xhardend :
		xtypesol = "opt"
		print("*** Looking (1) to optimize order: ",worstord, " from route # ", worstrun, " Late = ", worsttw)
	elif klastrun >= 1 :
		##########################################################################################################################
		#### MATRICS PART
		##########################################################################################################################
		totdur = 0
		totdist = 0
		maxdur = 0
		maxdist = 0
		for k in range(klastrun+1) :
			routedist = 0
			xcom   = "Route [" + str(k) + "] = "
			totcom = "            "
			twcom  = "            "
			i = 0
			j = 0
			xcom   = xcom   + str(i) + ", "
			totcom = totcom + str(i) + "[" + str(0) + "], "
			twcom = twcom + str(0) + ", "
			while (j<>999) :
				for j in range(ksize) :
					if xsol[k][i][j] :
						xcom    = xcom   + str(j) + ", "
						totcom  = totcom + str(j) + "[" + str(int(tstamp[k][j])) + "], "
						twcom = twcom + str(int(tstamp[k][j])) + ", "
						
						if j <> 0 : 
							last = int(tstamp[k][j] + tserv[j])
							if last > maxdur :
								maxdur = last
							totdist = totdist + cost[i][j]
							routedist = routedist + cost[i][j]
						if j == 0 :
							print(xcom)
							print(twcom)
							print(totcom)
							totdur = totdur + last
							if routedist > maxdist :
								maxdist = routedist
							j = 999
							break 
						i = j
						break
		totords = ksize-4.0
		totroutes = k + 1
		totdist = totdist 
		maxdist = maxdist
		print("")
		print("")
		print("**** Results & Metrics **** : ")
		print("       Total duration: ", totdur )
		print("       Total distance: ", round(totdist,2))
		print("     Number of routes: ", totroutes)
		print("     Number of orders: ", int(totords))
		print("Avg. distance / order: ", round(totdist/totords,2))
		print("  Avg. orders / route: ", round(totords/totroutes,2))
		print("Avg. distance / route: ", round(totdist/totroutes,2))
		print("  Max. route duration: ", maxdur)
		print("  Max. route distance: ", round(maxdist,2))
		print("     Max. route oders: ", maxorders)
		print("  Elapsed time (mins): ", round((endtime - initime)/60,2))
		print("============ OPTIMIZATION ENDS =================")
		print("")
		print("")
		if xquit :
			quit()
	else :
		klastrun = abs(klastrun) + 1
		print("Restart Routing....", klastrun)
		
nrun = klastrun 

################################################################################################################################################################
#### MAIN LOOP 
################################################################################################################################################################
while True:
	###############################################################################################################################
	# *** DEFINE INPUT VARIABLES
	###############################################################################################################################
	arcs = [[0 for j in range(ksize)] for i in range(ksize)]
	xbolarc = [[False for j in range(ksize)] for i in range(ksize)]	
	arrhr = [0 for i in range(ksize)] 
	pend  = [0 for i in range(ksize)] 
	late  = [0 for i in range(ksize)] 
	early = [0 for i in range(ksize)] 
	###############################################################################################################################	
			
	###############################################################################################################################
	# Instantiate a mixed-integer solver, naming it MYPROB.
	# **** Se Inicializan  las Variables Binarias ****
	###############################################################################################################################
	if nrun <> 0 and xcom <> "" :  
		del solver
		solver = 0
	
	solver = pywraplp.Solver('MYPROB', pywraplp.Solver.CBC_MIXED_INTEGER_PROGRAMMING)
	kvars = 0

	if xendopt or runopt :  		
		xactive = [False for i in range(ksize)]
		if runopt and not xendopt:
			xactive[worstord] = True  
		for i in range(ksize):
			for j in range(ksize):
				if xsol[nrun][i][j] :
					if xchdelord[i] <> nrun :
						xactive[i] = True
					if xchdelord[j] <> nrun :
						xactive[j] = True
					j = ksize
			if xchaddord[i] == nrun : 
				xactive[i] = True
				
	#### ***** REVISION DE LAS TIENDAS ACTIVAS ******
	kcnt1 = 0
	kcnt2 = 0
	kcnt3 = 0
	for i in range(ksize):
		if xactive[i] and not xstores[i] :
			if xdep[i] == kstore1 :
				kcnt1 += 1
			if xdep[i] == kstore2 :
				kcnt2 += 1
			if xdep[i] == kstore3 :
				kcnt3 += 1
	if kcnt1 == 0 :
		xactive[kstore1] = False
	else :
		xactive[kstore1] = True			
	if kcnt2 == 0 :
		xactive[kstore2] = False	
	else :
		xactive[kstore2] = True			
	if kcnt3 == 0 :
		xactive[kstore3] = False	
	else :
		xactive[kstore3] = True	

	if xendopt :  
		ktour = 0
		for i in range(ksize):
			if xactive[i] :
				ktour += 1
		xdif = ksize - ktour    
	else :	
		if runopt :
			xcount = 0
			for i in range(ksize) :
				if not xactive[i] :
					xcount += 1	
			xdif = xcount
			ktour = ksize - xdif  		
		else :
			xcount = 0
			for i in range(ksize) :
				if xactive[i] and not xstores[i] :
					xcount += 1	
			xdif = ksize - xcount
			ktour = min(ktour, ksize - xdif)   

	if xstartlen <> 0 :
		ktour = xstartlen
		xstartlen = 0
	print("***********************************************************  Route: ", nrun, " size: ", ksize - xdif, " Tour Length attempt = ", ktour) 

	# *** ACTUALIZAR LOS TIEMPOS DE PROCESO DE LAS TIENDAS
	if xactive[kstore1] :
		tserv[kstore1] = 5 * min(kcnt1,6) 
	if xactive[kstore2] :
		tserv[kstore2] = 5 * min(kcnt2,6) 
	if xactive[kstore3] :
		tserv[kstore3] = 5 * min(kcnt3,6) 
	
	for i in range(ksize):
		if xactive[i] :
			arrhr[i]   = solver.NumVar(0, bigm, 'arrhr[%i]' % i)
			pend[i]    = solver.BoolVar('pend[%i]' % i)	
			if xendopt : 
				solver.Add(pend[i] == 0)
			if not xstores[i] :
				late[i]  = solver.NumVar(0, xtol, 'late[%i]' % i)
				early[i] = solver.NumVar(0, xtol, 'early[%i]' % i)
			
			for j in range(ksize):
				if xactive[j] and cost[i][j] <= xthold and i<>j :
					if twa[i] + tserv[i] + cost[i][j] < twb[j] :
						kvars += 1
						xbolarc[i][j] = True
						arcs[i][j] = solver.BoolVar('arcs[%i],[%i]' % (i, j))
						if xendopt : 
							solver.Add(pend[j] == 0)

	print("Density = ", kvars, ": %", round(100*kvars/((ksize-xdif)**2),2)) 

	xprocst1 = solver.NumVar(0, bigm, 'xprocst1')	
	xprocst2 = solver.NumVar(0, bigm, 'xprocst2')	
	xprocst3 = solver.NumVar(0, bigm, 'xprocst3')	

	if runopt :
		xmaxlate  = solver.NumVar(0, xtol, 'xmaxlate')	
	
	###############################################################################################################################
	##### PROCESS TIME CONSTRAINTS AT STORES
	###############################################################################################################################
	#solver.Add(solver.Sum([ (1 - pend[i]) for i in (ii for ii in range(ksize) if xdep[ii] == kstore1 and ii<>kstore1 and xactive[ii]) ]) == xprocst1)	
	#solver.Add(solver.Sum([ (1 - pend[i]) for i in (ii for ii in range(ksize) if xdep[ii] == kstore2 and ii<>kstore2 and xactive[ii]) ]) == xprocst2)	
	#solver.Add(solver.Sum([ (1 - pend[i]) for i in (ii for ii in range(ksize) if xdep[ii] == kstore3 and ii<>kstore3 and xactive[ii]) ]) == xprocst3)	
	
	solver.Add( xprocst1 == tserv[kstore1] * (1 - pend[kstore1]) )
	solver.Add( xprocst2 == tserv[kstore2] * (1 - pend[kstore2]) )
	solver.Add( xprocst3 == tserv[kstore3] * (1 - pend[kstore3]) )

	if (ksize - xdif) >= 2*ktour :
		xheu = True
	else :
		xheu = False

	###############################################################################################################################
	##### HEURISTICS 
	###############################################################################################################################
	### HEU-1: Del DEPOT a las Tiendas
	if xheu or True :
		xcount = 0
		for i in range(1,ksize) :
			if xstores[i] and xactive[i] :
				xcount += 1	
		solver.Add(solver.Sum([ arcs[0][j] for j in (jj for jj in range(1,ksize) if xstores[jj] and xactive[jj]) ]) == 1)
		if ksize - xdif > ktour :     
			solver.Add(solver.Sum([ (1-pend[i]) for i in (ii for ii in range(1,ksize) if xstores[ii] and xactive[ii]) ]) <= min(2,xcount))
			solver.Add(solver.Sum([ arcs[i][j] for i in (ii for ii in range(1,ksize) if xstores[ii] and xactive[ii]) for j in (jj for jj in range(1,ksize) if xstores[jj] and xactive[jj]) ])  <= min(2,xcount))
		else :
			solver.Add(solver.Sum([ (1-pend[i]) for i in (ii for ii in range(1,ksize) if xstores[ii] and xactive[ii]) ]) <= xcount)
			solver.Add(solver.Sum([ arcs[i][j] for i in (ii for ii in range(1,ksize) if xstores[ii] and xactive[ii]) for j in (jj for jj in range(1,ksize) if xstores[jj] and xactive[jj]) ])  <= xcount)
			
	if xheu or True :
		### HEU-2: De los Clientes al DEPOT 
		xcount = 0
		for i in range(1,ksize) :
			if not xstores[i] and xactive[i] and twa[i]>=180 :
				xcount += 1	
		if xcount >= 1 :
			solver.Add(solver.Sum([ arcs[i][0] for i in (ii for ii in range(1,ksize) if not xstores[ii] and xactive[ii] and twa[ii]>=180) ]) == 1)
			solver.Add(solver.Sum([ arcs[i][0] for i in (ii for ii in range(1,ksize) if xactive[ii] and twa[ii]<180) ]) == 0)
		else :
			xcount = 0
			for i in range(1,ksize) :
				if not xstores[i] and xactive[i] and twa[i]>=120 :
					xcount += 1	
			if xcount >= 1 :
				solver.Add(solver.Sum([ arcs[i][0] for i in (ii for ii in range(1,ksize) if not xstores[ii] and xactive[ii] and twa[ii]>=120) ]) == 1)
				solver.Add(solver.Sum([ arcs[i][0] for i in (ii for ii in range(1,ksize) if xactive[ii] and twa[ii]<120) ]) == 0)

	xcount = 0
	for i in range(1,ksize) :
		if not xstores[i] and xactive[i] and twb[i]<=120 :
			xcount += 1	
	if xcount >= 1 and xheu :
		### HEU-3: De las Tiendas a los 1eros Clientes (Des-activado: Muchos Problemas)
		if True :   # Cuidado este es el nuevo 
			solver.Add( solver.Sum([ arcs[i][j] for i in (ii for ii in range(1,ksize) if xstores[ii]) for j in (jj for jj in range(1,ksize) if not xstores[jj] and twb[jj]<=120) ]) >= 1)
			solver.Add( solver.Sum([ arcs[i][j] for i in (ii for ii in range(1,ksize) if xstores[ii]) for j in (jj for jj in range(1,ksize) if not xstores[jj] and twb[jj]<=120) ]) <= 3 - pend[kstore1] - pend[kstore2] - pend[kstore3] )
			solver.Add( solver.Sum([ arcs[i][j] for i in (ii for ii in range(1,ksize) if xstores[ii]) for j in (jj for jj in range(1,ksize) if not xstores[jj] and twb[jj]>=180) ]) == 0)	
		for i in range(1,ksize) :
			if xstores[i] and xactive[i] :
				xcount = 0
				for j in range(1,ksize) :
					if not xstores[j] and xactive[j] and twb[j]<=120 :
						xcount += 1
				if xcount >= 1 :
					### Aqui abajo des-activado NO funciono
					#solver.Add( solver.Sum([ arcs[i][j] for j in (jj for jj in range(1,ksize) if not xstores[jj] and xactive[jj] and twb[jj]<=120) ]) == 1 - pend[i])
					### Cuidado aqui abajo es el nuevo 
					solver.Add( solver.Sum([ arcs[i][j] for j in (jj for jj in range(1,ksize) if not xstores[jj] and xactive[jj] and twb[jj]<=120 or xstores[jj] and xactive[jj]) ]) <= 1 - pend[i] )
					solver.Add(solver.Sum([ arcs[i][j] for j in (jj for jj in range(1,ksize) if not xstores[jj] and xactive[jj] and twb[jj]>120) ]) == 0)

	###############################################################################################################################
	##### RELAXED CONSTRAINTS
	##### OBJECTIVE: total cost, to be minimized
	###############################################################################################################################
	for i in range(ksize):
		if xactive[i] :
			#print(i,",",xdep[i],",",twa[i],",",twb[i])   

			# Time Windows Contraints
			solver.Add(arrhr[i] + early[i] >= twa[i] * (1-pend[i]))
			solver.Add(arrhr[i] - late[i]  <= twb[i] * (1-pend[i]))
			
			#  *** PRECAUCION 
			if runopt :
				solver.Add(late[i]  <= xmaxlate)
			
			# Time Windows Precedence Constraints  (Over dependent Store)
			if i <> xdep[i] :
				solver.Add(arrhr[i] + bigm*pend[i] >= arrhr[xdep[i]])
							
			solver.Add(solver.Sum([ arcs[i][j] for j in (jj for jj in range(ksize) if xbolarc[i][jj]) ]) + pend[i] == 1)
			solver.Add(solver.Sum([ arcs[j][i] for j in (jj for jj in range(ksize) if xbolarc[jj][i]) ]) + pend[i] == 1)
			
			for j in range(ksize):
				if xbolarc[i][j] :	
					solver.Add(arcs[i][j] + arcs[j][i] <= 1)
					solver.Add(arrhr[i] >= cost[0][i] * arcs[i][j])
					solver.Add(arrhr[0] >= cost[i][0] * arcs[i][0])
					
					if i <> 0 : 
						if not xstores[i] : 
							solver.Add(arrhr[i] + tserv[i] + cost[i][j] <= arrhr[j] + bigm*(1 - arcs[i][j]))    
						elif i == kstore1 : 
							solver.Add(arrhr[i] + xprocst1 + cost[i][j] <= arrhr[j] + bigm*(1 - arcs[i][j]))		
						elif i == kstore2 : 
							solver.Add(arrhr[i] + xprocst2 + cost[i][j] <= arrhr[j] + bigm*(1 - arcs[i][j]))		
						elif i == kstore3 :
							solver.Add(arrhr[i] + xprocst3 + cost[i][j] <= arrhr[j] + bigm*(1 - arcs[i][j]))
					else:
						solver.Add(cost[i][j] <= arrhr[j] + bigm*(1 - arcs[i][j]))

	xtour = min(ksize - xdif, ktour)  
	solver.Add( solver.Sum([ pend[i] for i in (ii for ii in range(1,ksize) if not xstores[ii] and xactive[ii]) ]) == ksize - xdif - xtour)   
	solver.Add(pend[kstore0] == 0)
	
	if xactive[kstore1] :
		solver.Add(solver.Sum([ (1 - pend[i]) for i in (ii for ii in range(ksize) if xdep[ii] == kstore1 and ii<>kstore1 and xactive[ii]) ]) >= 1 - pend[kstore1])	
		solver.Add(solver.Sum([ (1 - pend[i]) for i in (ii for ii in range(ksize) if xdep[ii] == kstore1 and ii<>kstore1 and xactive[ii]) ]) <= bigm*(1 - pend[kstore1]))

	if xactive[kstore2] :
		solver.Add(solver.Sum([ (1 - pend[i]) for i in (ii for ii in range(ksize) if xdep[ii] == kstore2 and ii<>kstore2 and xactive[ii]) ]) >= 1 - pend[kstore2])	
		solver.Add(solver.Sum([ (1 - pend[i]) for i in (ii for ii in range(ksize) if xdep[ii] == kstore2 and ii<>kstore2 and xactive[ii]) ]) <= bigm*(1 - pend[kstore2]))

	if xactive[kstore3] :
		solver.Add(solver.Sum([ (1 - pend[i]) for i in (ii for ii in range(ksize) if xdep[ii] == kstore3 and ii<>kstore3 and xactive[ii]) ]) >= 1 - pend[kstore3])	
		solver.Add(solver.Sum([ (1 - pend[i]) for i in (ii for ii in range(ksize) if xdep[ii] == kstore3 and ii<>kstore3 and xactive[ii]) ]) <= bigm*(1 - pend[kstore3]))

	#solver.Add(arrhr[???] == 0)    (1er cliente a visitar)								

	if runopt :
		z = solver.Sum([ cost[i][j] * arcs[i][j]  for i in range(ksize)  for j in (jj for jj in range(ksize) if xbolarc[i][jj]) ]) + \
				kmaxlate * solver.Sum([ late[i] + early[i] for i in range(ksize) ]) + kmaxlate * xmaxlate
	else :
		z = solver.Sum([ cost[i][j] * arcs[i][j]  for i in range(ksize)  for j in (jj for jj in range(ksize) if xbolarc[i][jj]) ]) + \
				kmaxlate * solver.Sum([ late[i] + early[i] for i in range(ksize) ])
		
	objective = solver.Minimize(z)
	
	###############################################################################################################################
	##### SOLVER PARAMETERS
	###############################################################################################################################
	
	mymodel = solver.ExportModelAsMpsFormat(False,False)	
	wfile    = "myprob"  + xtypesol + "-" + str(nrun) + ".mps"
	kfilesol = "sol"     + xtypesol + "-" + str(nrun) + ".sol"
	
	#scip -c "set presolving emphasis aggressive set emphasis feasibility set heuristics emphasis aggressive set limits time 200 read myprob-0.mps optimize write solution sol-0.sol quit"
	#xcom    = 'scip -c "set presolving emphasis aggressive set node selection estimate set emphasis feasibility set limits time 120 read ' + wfile + ' optimize write solution ' + kfilesol + ' quit"'
	xcom     = 'scip -c "set presolving emphasis aggressive set emphasis feasibility set heuristics emphasis aggressive set limits time 60 read ' + wfile + ' optimize write solution ' + kfilesol + ' quit"'
	
	if runopt :
		if os.path.isfile(kfilesol) :
			os.remove(kfilesol)

	fileout = open(wfile, 'w')
	print(mymodel, file=fileout)
	fileout.close()
		
	os.system(xcom)
	
	###############################################################################################################################
	##### OPEN SOLUTION FILE  
	###############################################################################################################################
	infile = open(kfilesol, "r")
	mybuf = csv.reader(infile, delimiter='\t')    
	if len(list(mybuf)) > 2 :
		infile.seek(0)		# Incia en la 2da linea 
		next(mybuf)
		next(mybuf)
		xerr = False
		for myrow in mybuf:			
			wvariable = myrow[0][0:25].strip()
			wvalue	  = float(myrow[0][26:])
			if wvariable[0:4] == "pend" and wvalue >= 0.99 :		
				kind = int(wvariable[5:wvariable.find(']')])
				if not xstores[kind] and not runopt:
					xactive[kind] = True	
			if wvariable[0:4] == "arcs" :
				kind1 = int(wvariable[5:wvariable.find(']')])
				kcoma = wvariable.find(',')
				kind2 = int(wvariable[kcoma+2:wvariable.find(']',kcoma)])
				xsol[nrun][kind1][kind2] = False
				if wvalue >= 0.99 : 
					xsol[nrun][kind1][kind2] = True
					print(kind1," : ",kind2)
					if not xstores[kind1] and not runopt:
						xactive[kind1] = False	
					if not xstores[kind2] and not runopt:
						xactive[kind2] = False	
			if wvariable[0:4] == "late" :
				kind1 = int(wvariable[5:wvariable.find(']')])
				xlate[nrun][kind1] = wvalue
				if wvalue >= 1 :
					print("Late: ", kind1," = ", wvalue)  
					if wvalue > worsttw :
						if runopt and not xendopt :   
							xerr = True	
							xchdelord[worstord] = nrun   
							xcntord[nrun] -= 1 
						else :
							worstrun = nrun			# ???? Se puede Borrar?
			if wvariable[0:5] == "Early" and wvalue >= 1 :
				kind1 = int(wvariable[6:wvariable.find(']')])
				print("Early: ", kind1," = ",wvalue)
				#raw_input("Press Enter to continue...")
				
		infile.close()
		
		if xerr :  
			### EXEPCION ANTICIPADA: TERMINAR DE ITERAR EN OPTIMIZACION. APLICAR TODOS LOS CAMBIOS Y EJECUTAR LAS RUTAS FINALES. GO TO LOOP  
			xcount = 20
			xrouteinc = nrun
			for krun in range(klastrun+1):  
				xsum = numpy.sum(xsol[krun]) 
				print(krun,",",numpy.sum(xsol[krun]),",",xcntord[krun],"... Not store dependent")
				if xsum < xcount : 
					xcount = xsum
					xrouteinc = krun
			print("*** Looking (4) to optimize order: ",worstord, " from route # ", worstrun, " Late = ", worsttw, ". Incumbent Route: ", xrouteinc) 
			nrun = xrouteinc
		else :	
			print()
			print("***********************************************************  Route #: ", nrun, " = ", xtour)   
			if not runopt or xendopt :
				if xendopt :  
					if nrun == klastrun :   
						print("***************************************** Total Termination *****************************************")
						###  *** AQUI TERMINA TODO (VAMOS PARA ABAJO / REPORT METRICS ONLY)
						xtypesol = "opt" 
						xhardend = True
						break
					else :
						### IR A LA SIGUIENTE ITERACION DE RUTAS YA OPTIMIZADAS. GO TO LOOP 
						nrun += 1
				else :
					xcount = 0
					for i in range(ksize) :
						if xactive[i] and not xstores[i] :
							xcount += 1	
					if xcount == 0 :
						# ****** AQUI EMPIEZA LA FASE DE OPTIMIZACION (VAMOS PARA ABAJO / GO TO START-OVER)  *****
						break
					else :
						# *** TODO OK: PROCESAR SIGUIENTE RUTA USANDO MISMA LONGITUD. GO TO LOOP 
						nrun += 1
			else :
				# *** ITERACION: PROCESO DE OPTIMIZACION. ACTUALIZAR VECTORES DE ADDS / DELETES
				xchaddord[worstord] = nrun
				xcntord[nrun] += 1 
				xchdelord[worstord] = worstrun
				xcntord[worstrun] -= 1 
				print("Route Updates: ADDORD(",worstord,":",nrun,"). DELORD(",worstord,":",worstrun,")")
				# *** IMPORTANTE: ACTUALIZAR RUTAS AFECTADAS. ELIMINAR EL TIEMPO DE RETRASO DE LA ORDEN INCUMBENTE
				xlate[worstrun][worstord] = 0  
								
				# *** ITERACION: PROCESO DE OPTIMIZACION. IDENTIFICAR LA ORDEN CON MAS RETRASO  
				worsttw = 0
				for krun in range(klastrun+1):  
					for i in range(ksize):
						if xlate[krun][i] > worsttw :  
							worsttw  = xlate[krun][i]
							worstord = i
							worstrun = krun
							
				# *** ITERACION: PROCESO DE OPTIMIZACION. IDENTIFICAR LA RUTA CON MENOS ORDENES QUE INCLUYA LA TIENDA DEPENDIENTE DE LA ORDEN CON MAS RETRASO  
				xcount = 20
				xrouteinc = nrun
				for krun in range(klastrun+1):  
					if numpy.sum(xsol[krun][xdep[worstord]]) >= 1 :  
						xsum = numpy.sum(xsol[krun]) 
						print(krun,",",numpy.sum(xsol[krun]),",",xcntord[krun],"... store dependent")
						if xsum < xcount : 
							xcount = xsum
							xrouteinc = krun
						
				if worsttw >= 0.5 :
					print("*** Looking (2) to optimize order: ",worstord, " from route # ", worstrun, " Late = ", worsttw, ". Incumbent Route: ", xrouteinc) 
					nrun = xrouteinc
				else :
					### OPTIMO GENERAL: TERMINAR DE ITERAR EN OPTIMIZACION. APLICAR TODOS LOS CAMBIOS Y EJECUTAR LAS RUTAS FINALES. GO TO LOOP
					xendopt = True   
					xtypesol = "opt"
					nrun = 0
					print("***********************************************************  Optimization ENDING ... ")				
	else:
		# *** PROCESAR MISMA RUTA PERO USANDO LONGITUD = (N-1). GO TO LOOP
		# raw_input("No solution found. Press Enter to continue...")
		if not runopt :
			ktour = ktour - 1
		else :
			xchdelord[worstord] = nrun   
			xcntord[nrun] -= 1 
			
			# *** ITERACION: PROCESO DE OPTIMIZACION. IDENTIFICAR LA RUTA CON MENOS ORDENES DE CUALQUIER TIPO
			xcount = 20
			xrouteinc = nrun
			for krun in range(xrouteinc+1,klastrun+1):  
				if numpy.sum(xsol[krun][xdep[worstord]]) >= 1 :  
					xsum = numpy.sum(xsol[krun]) 
					print(krun,",",numpy.sum(xsol[krun]),",",xcntord[krun],"...store dependent (2)")
					if xsum < xcount : 
						xcount = xsum
						xrouteinc = krun
			if xcount == 20 :
				for krun in range(klastrun+1):  
					xsum = numpy.sum(xsol[krun]) 
					print(krun,",",numpy.sum(xsol[krun]),",",xcntord[krun],"... Not store dependent")
					if xsum < xcount : 
						xcount = xsum
						xrouteinc = krun
			if xrouteinc <> nrun : 
				print("*** Looking (3) to optimize order: ",worstord, " from route # ", worstrun, " Late = ", worsttw, ". Incumbent Route: ", xrouteinc) 
				nrun = xrouteinc
			else :
				### OPTIMO GENERAL: TERMINAR DE ITERAR EN OPTIMIZACION. APLICAR TODOS LOS CAMBIOS Y EJECUTAR LAS RUTAS FINALES. GO TO LOOP
				xendopt = True   
				xtypesol = "opt"
				nrun = 0
				print("***********************************************************  Optimization ENDING ... ")				
			
			
# ****** AQUI EMPIEZA LA FASE DE OPTIMIZACION  (VAMOS DE PARA ARRIBA)  *******
klastrun = nrun
if wantopt :
	if not runopt :  
		print("*************************** ROUTE RE-BUILDING STARTING ******************************************")
		runopt  = True   
		
goto .mystart
