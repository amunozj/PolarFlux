__authors__ = ["Zach Werginz", "Andres Munoz-Jaramillo"]
__email__ = ["zachary.werginz@snc.edu", "amunozj@gsu.edu"]

import numpy as np
import sunpy.map
from sunpy.sun import constants
from sunpy.sun import sun
import astropy.units as u

#np.set_printoptions(threshold=np.nan)

RSUN_METERS = sun.constants.radius.si.value
DSUN_METERS = sun.constants.au.si.value


class CRD:
	"""Calculates various magnetogram coordinate information.
	Can calculate heliographic coordinate information,
	line of sight (LOS) corrections for the magnetic field,
	area elements for each pixel, and magnetic flux. This can
	be done for one pixel, or the whole data map. If the whole
	data map is given as a parameter, it will save the information
	as an instance attribute for the object.

	"""
	
	def __init__(self, filename):
		"""Reads magnetogram as a sunpy.map object."""
		self.im_raw = sunpy.map.Map(filename)

		if self.im_raw.detector == '512':

			# Define center of sun and location of detector.
			self.X0 = self.im_raw.meta['CRPIX1A']
			self.Y0 = self.im_raw.meta['CRPIX2A']
			self.B0 = self.im_raw.meta['B0']
			self.L0 = self.im_raw.meta['L0']
			self.rsun = self.im_raw.rsun_obs.value
			self.dsun = DSUN_METERS

		elif self.im_raw.detector == 'SPMG':

			# Define center of sun and location of detector.
			self.X0 = self.im_raw.meta['CRPIX1A']
			self.Y0 = self.im_raw.meta['CRPIX2A']
			self.B0 = self.im_raw.meta['B0']
			self.L0 = self.im_raw.meta['L0']
			self.rsun = self.im_raw.rsun_obs.value / self.im_raw.meta['SCALE']
			self.dsun = DSUN_METERS

		elif self.im_raw.detector == 'MDI':

			# Define center of sun and location of detector.
			self.X0 = self.im_raw.meta['X0']
			self.Y0 = self.im_raw.meta['Y0']
			self.B0 = self.im_raw.meta['B0']
			self.L0 = self.im_raw.meta['L0']
			self.rsun = self.im_raw.meta['R_SUN']
			self.dsun = self.im_raw.meta['OBS_DIST']/.0046491

		elif self.im_raw.detector == 'HMI':
			self.X0 = self.im_raw.meta['CRPIX1']
			self.Y0 = self.im_raw.meta['CRPIX2']
			self.B0 = self.im_raw.meta['CRLT_OBS']
			self.L0 = self.im_raw.meta['CRLN_OBS']
			self.rsun = self.im_raw.rsun_obs.value
			self.dsun = self.im_raw.dsun.value

		# try:
		# 	self.B0 = self.im_raw.meta['B0']
		# except KeyError:
		# 	self.B0 = self.im_raw.meta['OBS_B0']
		# try:
		# 	self.L0 = self.im_raw.meta['L0']
		# except KeyError:
		# 	self.L0 = self.im_raw.meta['OBS_L0']
		# try:
		# 	self.X0 = self.im_raw.meta['CRPIX1']
		# except KeyError:
		# 	self.X0 = self.im_raw.meta['CRPIX1A']
		# except KeyError:
		# 	self.X0 = self.im_raw.meta['X0']
		# except KeyError:
		# 	self.X0 = self.im_raw.meta['IMG_X0']
		# try:
		# 	self.Y0 = self.im_raw.meta['Y0']
		# except KeyError:
		# 	self.Y0 = self.im_raw.meta['IMG_Y0']
		# if self.im_raw.detector == 'SPMG':
		# 	self.rsun = self.im_raw.rsun_obs.value / self.im_raw.meta['SCALE']	
		
	def __repr__(self):
		#TODO
		return None

	def heliographic(self, *args, corners=False):
		"""Calculate hg coordinates from hpc and returns it.

		Can accept either a coordinate pair (x, y) or an entire 2D array.
		This pair corresponds the the pixel you want information on.

		The function will return a coordinate pair if given a coordinate
		pair or two arrays (latitude and longitude) if given a data array.
		Those arrays are indexed such that [0,0] is the top left pixel.

		Use standard python indexing conventions for both the single
		coordinate and array calculations [row, column].

		Examples: 
		lath, lonh = kpvt.heliographic(kpvt.im_raw.data)
		aia.heliographic(320, 288)
		"""

		xScl = self.im_raw.scale[0].value
		yScl = self.im_raw.scale[1].value
		
		# Check for single coordinate or ndarray object.
		if isinstance(args[0], np.ndarray):
			# Retrieve integer dimensions and create arrays holding
			# x and y coordinates of each pixel
			xdim = np.int(np.floor(self.im_raw.dimensions[0].value))
			ydim = np.int(np.floor(self.im_raw.dimensions[1].value))
			try:
				xrow = (np.arange(0, xdim) - self.X0 + args[1])*xScl
				yrow = (np.arange(0, ydim) - self.Y0 + args[2])*yScl
				self.xg, self.yg = np.meshgrid(xrow, yrow, indexing='xy')
				self.rg = np.sqrt(self.xg**2 + self.yg**2)
				x = self.xg
				y = -self.yg
			except IndexError:
				if corners:
					xrow = (np.arange(0, xdim + 1) - self.X0 - 0.5)*xScl
					yrow = (np.arange(0, ydim + 1) - self.Y0 - 0.5)*yScl
				else:
					xrow = (np.arange(0, xdim) - self.X0)*xScl
					yrow = (np.arange(0, ydim) - self.Y0)*yScl
				self.xg, self.yg = np.meshgrid(xrow, yrow, indexing='xy')
				self.rg = np.sqrt(self.xg**2 + self.yg**2)
				x = self.xg
				y = -self.yg
		else:
			# Have to switch coordinate conventions because calculations
			# assume standard cartesian whereas python indexing is 
			# [row, column]
			x = (args[1] - self.X0)*xScl
			y = (self.Y0 - args[0])*yScl

		# First convert to heliocentric cartesian coordinates.
		# Calculations taken from sunpy.wcs.
		x *= np.deg2rad(1)/3600.0
		y *= np.deg2rad(1)/3600.0

		# cosx = np.cos(x*np.deg2rad(1)/3600.0)
		# sinx = np.sin(x*np.deg2rad(1)/3600.0)
		# cosy = np.cos(y*np.deg2rad(1)/3600.0)
		# siny = np.sin(y*np.deg2rad(1)/3600.0)

		q = self.dsun * np.cos(y) * np.cos(x)
		distance = q**2 - self.dsun**2 + RSUN_METERS**2
		distance = q - np.sqrt(distance)

		rx = distance * np.cos(y) * np.sin(x)
		ry = distance * np.sin(y)
		rz = np.sqrt(RSUN_METERS**2 - rx**2 - ry**2)

		# Now convert to heliographic coordinates.
		cosb = np.cos(np.deg2rad(self.B0))
		sinb = np.sin(np.deg2rad(self.B0))

		hecr = np.sqrt(rx**2 + ry**2 + rz**2)
		hgln = np.arctan2(rx, rz*cosb - ry*sinb) \
		       + np.deg2rad(self.L0)
		hglt = np.arcsin((ry * cosb + rz * sinb)/hecr)
		# Only add the instance attribute if it doesn't exist.
		if isinstance(args[0], np.ndarray) and not hasattr(self, 'lonh'):
			self.lonh = np.rad2deg(hgln)
			self.lath = np.rad2deg(hglt)

		return np.rad2deg(hgln), np.rad2deg(hglt)

	def los_corr(self, *args):
		"""Takes in coordinates and returns corrected magnetic field.

		Applies the dot product between the observers unit vector and
		the heliographic radial vector to get the true magnitude of 
		the magnetic field vector. See geometric projection for
		calulations.
		"""

		# TODO optimize with saved lath, lonh data
		if isinstance(args[0], np.ndarray):
			lonh, lath = np.deg2rad(self.heliographic(args[0]))
		else:
			lonh, lath = np.deg2rad(self.heliographic(args[0], args[1]))

		Xobs = np.cos(np.deg2rad(self.B0))*np.cos(np.deg2rad(self.L0))
		Yobs = np.cos(np.deg2rad(self.B0))*np.sin(np.deg2rad(self.L0))
		Zobs = np.sin(np.deg2rad(self.B0))

		corr_factor = (np.cos(lath)*np.cos(lonh)*Xobs
					   + np.cos(lath)*np.sin(lonh)*Yobs
					   + np.sin(lath)*Zobs)
		if isinstance(args[0], np.ndarray):
			self.im_corr = self.im_raw.data/corr_factor
			return self.im_corr
		else:
			return self.im_raw.data[args[0], args[1]]/corr_factor

	def eoa(self, *args):
		""" Takes in coordinates and returns the area of pixels on sun.

		Each pixel is projected onto the sun, and therefore pixels close to
		the limbs have vastly greater areas. This function uses a closed form
		solution to a spherical area integral to calulate the area based on
		the heliographic coordinate unit vectors of each corner of the pixel.
		We use these to calculate a solid angle of a pyramid with its apex
		at the center of the sun.
		"""

		#Assume coordinate is in center of pixel.
		#Information on pixel standard is in this article.
		#http://www.aanda.org/component/article?access=bibcode&bibcode=&bibcode=2002A%2526A...395.1061GFUL
		if isinstance(args[0], np.ndarray):
			lon, lat = self.heliographic(args[0], corners=True)
			# Calculating unit vectors of pixel corners for solid angle.
			r1 = np.array([np.cos(np.deg2rad(lat[0:len(lat) - 2, 0: len(lat) - 2]))*np.cos(np.deg2rad(lon[0:len(lat) - 2, 0: len(lat) - 2])),
							np.cos(np.deg2rad(lat[0:len(lat) - 2, 0: len(lat) - 2]))*np.sin(np.deg2rad(lon[0:len(lat) - 2, 0: len(lat) - 2])),
							np.sin(np.deg2rad(lat[0:len(lat) - 2, 0: len(lat) - 2]))])

			r2 = np.array([np.cos(np.deg2rad(lat[1:len(lat) - 1, 0: len(lat) - 2]))*np.cos(np.deg2rad(lon[1:len(lat) - 1, 0: len(lat) - 2])),
							np.cos(np.deg2rad(lat[1:len(lat) - 1, 0: len(lat) - 2]))*np.sin(np.deg2rad(lon[1:len(lat) - 1, 0: len(lat) - 2])),
							np.sin(np.deg2rad(lat[1:len(lat) - 1, 0: len(lat) - 2]))])

			r3 = np.array([np.cos(np.deg2rad(lat[1:len(lat) - 1, 1: len(lat) - 1]))*np.cos(np.deg2rad(lon[1:len(lat) - 1, 1: len(lat) - 1])),
							np.cos(np.deg2rad(lat[1:len(lat) - 1, 1: len(lat) - 1]))*np.sin(np.deg2rad(lon[1:len(lat) - 1, 1: len(lat) - 1])),
							np.sin(np.deg2rad(lat[1:len(lat) - 1, 1: len(lat) - 1]))])

			r4 = np.array([np.cos(np.deg2rad(lat[0:len(lat) - 2, 1: len(lat) - 1]))*np.cos(np.deg2rad(lon[0:len(lat) - 2, 1: len(lat) - 1])),
							np.cos(np.deg2rad(lat[0:len(lat) - 2, 1: len(lat) - 1]))*np.sin(np.deg2rad(lon[0:len(lat) - 2, 1: len(lat) - 1])),
							np.sin(np.deg2rad(lat[0:len(lat) - 2, 1: len(lat) - 1]))])

		else:
			x = args[0]
			y = args[1]
			lonUL, latUL = self.heliographic(x - .5, y - .5)
			lonLL, latLL = self.heliographic(x + .5, y - .5)
			lonLR, latLR = self.heliographic(x + .5, y + .5)
			lonUR, latUR = self.heliographic(x - .5, y + .5)

			# Calculating unit vectors of pixel corners for solid angle.
			r1 = np.array([np.cos(np.deg2rad(latUL))*np.cos(np.deg2rad(lonUL)),
							np.cos(np.deg2rad(latUL))*np.sin(np.deg2rad(lonUL)),
							np.sin(np.deg2rad(latUL))])

			r2 = np.array([np.cos(np.deg2rad(latLL))*np.cos(np.deg2rad(lonLL)),
							np.cos(np.deg2rad(latLL))*np.sin(np.deg2rad(lonLL)),
							np.sin(np.deg2rad(latLL))])

			r3 = np.array([np.cos(np.deg2rad(latLR))*np.cos(np.deg2rad(lonLR)),
							np.cos(np.deg2rad(latLR))*np.sin(np.deg2rad(lonLR)),
							np.sin(np.deg2rad(latLR))])

			r4 = np.array([np.cos(np.deg2rad(latUR))*np.cos(np.deg2rad(lonUR)),
							np.cos(np.deg2rad(latUR))*np.sin(np.deg2rad(lonUR)),
							np.sin(np.deg2rad(latUR))])

		# Calculate solid angle of pixel based on a pyrimid shaped polygon.
		# See 
		cross1 = np.cross(r1, r2, axis=0)
		cross2 = np.cross(r3, r4, axis=0)
		numerator1 = dot(cross1, r3)
		numerator2 = dot(cross2, r1)
		solid_angle1 = 2*np.arctan2(numerator1,
						(dot(r1, r2) + dot(r2, r3) + dot(r3, r1) + 1))
		solid_angle2 = 2*np.arctan2(numerator2, 
						(dot(r3, r4) + dot(r4, r1) + dot(r3, r1) + 1))
		solid_angle = solid_angle1 + solid_angle2

		r = RSUN_METERS*100 # Convert to centimeters
		if isinstance(args[0], np.ndarray):
			self.area = np.abs((r**2)*solid_angle)
			ind = np.where(self.rg[1:len(self.rg)-1, 1:len(self.rg)-1] > self.rsun)
			self.area[ind] = np.nan
			return self.area
		else:
			return np.abs((r**2)*solid_angle)

	def magnetic_flux(self, *args, raw_field=False):
		""" Takes in coordinates and returns magnetic flux of pixel."""
		# Use existing attributes if saved.
		if hasattr(self, 'area'):
			area = self.area
		else:	
			area = self.eoa(*args)

		if raw_field:
			field = self.im_raw.data
			if isinstance(args[0], np.ndarray) and not hasattr(self, 'mgnt_flux_raw'):
				self.mgnt_flux_raw = area*field
			return self.mgnt_flux_raw
		else:
			if not hasattr(self, 'im_corr'):
				field = self.los_corr
			else:
				field = self.im_corr
			if isinstance(args[0], np.ndarray) and not hasattr(self, 'mgnt_flux'):
				self.mgnt_flux_corr = area*field
			return self.mgnt_flux_corr

def dot(a, b):
	return  np.array(a[0]*b[0] + a[1]*b[1] + a[2]*b[2])