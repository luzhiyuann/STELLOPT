##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This library provides a python class for calculating
quantities related to collisions
"""

# Libraries
#from libstell import libstell

# Constants
EC = 1.602176634E-19 # Electron charge [C]
HBAR = 1.054571817E-34 # Planks reduced constant [J.s]
DA = 1.66053906660E-27 # Dalton
ME = 9.1093837E-31 # Electron mass [kg]
C  = 299792458 # Speed of light [m/s]
EPS0 = 8.8541878188E-12 # Vacuum permittivity [F/m]

# VMEC Class
class COLLISIONS():
	"""Class for working with FUSION data

	"""
	def __init__(self):
		# Electron Mass
		self.ME = ME
		# Proton Mass
		self.MP = 1.007276466621*DA 
		# Deuturium Mass
		self.MD = 2.01410177811*DA
		# Tritium Mass
		self.MT = 3.01604928*DA
		# Helium Mass
		self.MHe3 = 3.0160293*DA
		self.MHe4 = 4.002603254*DA

	def coullog_ee(self,ne,te):
		"""Computes the Coulomb e-e logarithm

		This routine calculates the Coulomb Logarithm acording to the
		NRL plasma formulary general deffinition.

		Parameters
		----------
		ne : real
			Electron density [m^-3]
		te : real
			Electron Temperature [eV]
		Returns
		----------
		clog : real
			Coulomb logarithm
		"""
		import numpy as np
		ne_cm = ne*1E-6
		a1 = -np.log(np.sqrt(ne_cm)*te**(-1.25))
		a2 = -np.sqrt(1E-5 + ((np.log(te)-2)**2)/16.0)
		return 23.5 + a1 + a2

	def coullog_ei(self,ne,te,mi,Z,ni,ti):
		"""Computes the Coulomb e-i logarithm

		This routine calculates the Coulomb Logarithm acording to the
		NRL plasma formulary general deffinition.

		Parameters
		----------
		ne : real
			Electron density [m^-3]
		te : real
			Electron Temperature [eV]
		mi : real
			Ion mass [kg]
		Z : real
			Ion Charge number
		ni : real
			Ion density [m^-3]
		ti : real
			Ion Temperature [eV]
		Returns
		----------
		clog : real
			Coulomb logarithm
		"""
		import numpy as np
		ne_cm = ne*1E-6
		mu = mi/self.MP
		# No need to convert masses from kg to g
		if ti/mi < (te/self.ME):
			if (te < (10 * Z * Z)):
				clog = 23 - np.log(np.sqrt(ne_cm)*Z*te**(-1.5))
			else:
				clog = 24 - np.log(np.sqrt(ne_cm)/te)
		else:
			mu = (self.ME*mi)/(self.ME+mi)
			clog  = 16 - np.log(mu*Z*Z*np.sqrt(ni)*ti**-1.5)
		return clog

	def coullog_ii(self,mi1,Z1,ni1,ti1,mi2,Z2,ni2,ti2):
		"""Computes the Coulomb i-i logarithm

		This routine calculates the Coulomb Logarithm acording to the
		NRL plasma formulary general deffinition.

		Parameters
		----------
		mi1 : real
			Ion #1 mass [kg]
		Z1  : real
			Ion #1 Charge number
		ni1 : real
			Ion #1 density [m^-3]
		ti1 : real
			Ion #1 Temperature [eV]
		mi2 : real
			Ion #2 mass [kg]
		Z2  : real
			Ion #2 Charge number
		ni2 : real
			Ion #2 density [m^-3]
		ti2 : real
			Ion #2 Temperature [eV]
		Returns
		----------
		clog : real
			Coulomb logarithm
		"""
		import numpy as np
		ni1_cm = ni1*1E-6
		ni2_cm = ni2*1E-6
		mu1 = mi1/self.MP
		mu2 = mi2/self.MP
		# No need to convert masses from kg to g
		clog = (ni1_cm*Z1*Z1/ti1) + (ni2_cm*Z2*Z2/ti2)
		clog = Z1*Z2*(mu1+mu2)*np.sqrt(clog)/(mu1*ti2+mu2*ti1)
		clog = 23 - np.log(clog)
		return clog

	def coullog_iifast(self,ne,te,mi1,Z1,mi2,Z2,vd):
		"""Computes the Coulomb counterstreaming i-i logarithm

		This routine calculates the Coulomb Logarithm acording to the
		NRL plasma formulary general deffinition.

		Parameters
		----------
		ne : real
			Electron density [m^-3]
		te : real
			Electron Temperature [eV]
		mi1 : real
			Ion #1 mass [kg]
		Z1  : real
			Ion #1 Charge number
		mi2 : real
			Ion #2 mass [kg]
		Z2  : real
			Ion #2 Charge number
		vd  : real
			Ion relative velocity [m/s]
		Returns
		----------
		clog : real
			Coulomb logarithm
		"""
		import numpy as np
		ne_cm = ne*1E-6
		mu1 = mi1/self.MP
		mu2 = mi2/self.MP
		# No need to convert masses from kg to g
		beta = vd/C
		clog = ne_cm/te
		clog = Z1*Z2*(mu1+mu2) * np.sqrt(clog) / (mu1*mu2*beta*beta)
		clog = 43-np.log(clog)
		return clog

	def collisionfreq_thermal_equilibration(self,m1,Z1,T1,m2,Z2,n2,T2,clog):
		"""Computes the collision frequency for two species

		This routine calculates the collision frequency acording to the
		NRL plasma formulary definition of Thermal Equilibration.

		Parameters
		----------
		m1 : real
			Particle #1 mass [kg]
		Z1 : real
			Particle #1 Charge number
		T1 : real
			Particle #1 Temperature [eV]
		m2 : real
			Particle #2 mass [kg]
		Z2 : real
			Particle #2 Charge number
		n2 : real
			Particle #2 Density [m^-3]
		T2 : real
			Particle #2 Temperature [eV]
		clog : real
			Coulomb logarithm
		Returns
		----------
		freq : real
			Collision frequency
		"""
		import numpy as np
		n2_cm = n2*1E-6
		# factor of 10^(-3/2) comes from mass in kg to g
		freq = 1.8E-19*np.sqrt(m1*m2)*Z1*Z1*Z2*Z2*n2_cm*clog*10**(-1.5)
		freq = freq / (m1*T2+m2*T1)**1.5
		return freq

	def psi(self,x):
		# psi function as defined in the NRL plasma formulary
		from scipy.special import gammainc, gamma
		from numpy import sqrt, pi
  
		return (2.0/sqrt(pi))*gammainc(3./2,x)*gamma(3./2)

	def psip(self,x):
		# derivative of psi function
		from numpy import sqrt, pi, exp
  
		return (2.0/sqrt(pi)) * sqrt(x) * exp(-x)
		
	def collisionfreq_slowingdown(self,m1,Z1,v1,m2,Z2,T2,n2,clog):
		"""Computes the slowing-down collision frequency of particle 1
  		inside thermal bath of species 2
		according to NRL plasma formulary

		Parameters
		----------
		m1 : real
			Particle #1 mass [kg]
		Z1 : real
			Particle #1 Charge number
		v1 : real
			Particle #1 Velocity [m/s]
		m2 : real
			Particle #2 mass [kg]
		Z2 : real
			Particle #2 Charge number
		T2 : real
			Particle #2 Temperature [eV]
		n2 : real
			Particle #2 Density [m^-3]
		clog : real
			Coulomb logarithm
		Returns
		----------
		freq : real
			Collision frequency
		"""
  
		# Note that quantities are given in SI units. In order to recover the formula
		# of the NRL formulary, which is in gaussian units, one only needs to
		# make the substitution: eps0 --> 1/4pi
     
		import numpy as np
		
		x = m2*v1*v1 / (2*EC*T2)

		freq = Z1**2 * Z2**2 * EC**4 * clog * n2 / (m1**2 * v1**3 * 4*np.pi * EPS0**2)
		freq = freq * (1+m1/m2) * self.psi(x)

		return freq

	def collisionfreq_transverse(self,m1,Z1,v1,m2,Z2,T2,n2,clog):
		"""Computes the transverse diffusion collision frequency 
  		of particle 1 inside thermal bath of species 2
		according to NRL plasma formulary

		Parameters
		----------
		m1 : real
			Particle #1 mass [kg]
		Z1 : real
			Particle #1 Charge number
		v1 : real
			Particle #1 Velocity [m/s]
		m2 : real
			Particle #2 mass [kg]
		Z2 : real
			Particle #2 Charge number
		T2 : real
			Particle #2 Temperature [eV]
		n2 : real
			Particle #2 Density [m^-3]
		clog : real
			Coulomb logarithm
		Returns
		----------
		freq : real
			Collision frequency
		"""
  
		# Note that quantities are given in SI units. In order to recover the formula
		# of the NRL formulary, which is in gaussian units, one only needs to
		# make the substitution: eps0 --> 1/4pi
     
		import numpy as np
  
		x = m2*v1*v1 / (2*EC*T2)

		freq = Z1**2 * Z2**2 * EC**4 * clog * n2 / (m1**2 * v1**3 * 4*np.pi * EPS0**2)
		freq = freq * 2 * ( (1-0.5/x)*self.psi(x) + self.psip(x) )
  
		return freq

	def collisionfreq_PENTA(self,v1,m,Z,T,n,clog):
		"""Computes the perpendicular collision frequency 
  		of particle 1 inside thermal bath of species [1,2,...,n]
		as in the PENTA code

		Parameters
		----------
		v1 : real
			Particle #1 Velocity [m/s]
		m : list
			Particles [#1,#2,...,#n] mass [kg]
		Z : list
			Particles [#1,#2,...,#n] Charge number
		T : list
			Particles [#1,#2,...,#n] Temperature [eV]
		n : list
			Particles [#1,#2,...,#n] Density [m^-3]
		clog : list
			Coulomb logarithm [lambda_11,lambda_12,...lambda_1n]
		Returns
		----------
		freq : list
			Collision frequency [nu_11,nu_12,...,nu_1n]
		"""
     
		import numpy as np
		from scipy.special import erf
  
		m1 = m[0]
		Z1 = Z[0]

		freq = []
		for m2,Z2,T2,n2,clog2 in zip(m,Z,T,n,clog):
  
			x = m2*v1*v1 / (2*EC*T2)

			nu = Z1**2 * Z2**2 * EC**4 * clog2 * n2 / (m1**2 * v1**3 * 4*np.pi * EPS0**2)
			nu = nu * ( (1-0.5/x)*erf(np.sqrt(x)) + np.exp(-x)/np.sqrt(x*np.pi) )
   
			freq.append( nu )
  
		return freq

	def collisionfreq_energy(self,m1,Z1,v1,m2,Z2,T2,n2,clog):
		"""Computes the energy loss collision frequency 
  		of particle 1 inside thermal bath of species 2
		according to NRL plasma formulary

		Parameters
		----------
		m1 : real
			Particle #1 mass [kg]
		Z1 : real
			Particle #1 Charge number
		v1 : real
			Particle #1 Velocity [m/s]
		m2 : real
			Particle #2 mass [kg]
		Z2 : real
			Particle #2 Charge number
		T2 : real
			Particle #2 Temperature [eV]
		n2 : real
			Particle #2 Density [m^-3]
		clog : real
			Coulomb logarithm
		Returns
		----------
		freq : real
			Collision frequency
		"""
  
		# Note that quantities are given in SI units. In order to recover the formula
		# of the NRL formulary, which is in gaussian units, one only needs to
		# make the substitution: eps0 --> 1/4pi
     
		import numpy as np

		x = m2*v1*v1 / (2*EC*T2)

		freq = Z1**2 * Z2**2 * EC**4 * clog * n2 / (m1**2 * v1**3 * 4*np.pi * EPS0**2)
		freq = freq * 2 * ( (m1/m2)*self.psi(x) - self.psip(x) )
  
		return freq

	def tauspitzer(self,mi,Zi,ne,Te,clog):
		"""Computes the Spitzer ion-electron exchange time

		This routine calculates the Spitzer ion-electron momentum
		exchange time.

		Parameters
		----------
		mi : real
			Ion Mass [kg]
		Zi : real
			Ion Charge number
		ne : real
			Electron Density [m^-3]
		Te : real
			Electron Temperature [eV]
		clog : real
			Coulomb logarithm
		Returns
		----------
		tau : real
			Spitzer ion-electron momentum exchange time
		"""
		import numpy as np
		return 3.777183E41*mi*np.sqrt(Te*Te*Te)/(ne*Zi*Zi*clog)

	def criticalvelocity(self,mp,Te):
		"""Computes the critical velocity

		This routine calculates the critical velocity for an energetic
		particle.

		Parameters
		----------
		mi : real
			Plasma Ion Mass [kg]
		Te : real
			Electron Temperature [eV]
		Returns
		----------
		vcrit : real
			Critical velocity [m/s]
		"""
		import numpy as np
		return (0.75*np.sqrt(np.pi*mp/self.ME))**(1.0/3.0) * np.sqrt(2*EC*Te/mp)

if __name__=="__main__":
	import sys
	sys.exit(0)
