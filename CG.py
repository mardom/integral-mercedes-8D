import numpy as  np
import astropy.units as u
import astropy.coordinates as coorde
from astropy.coordinates import SkyCoord
from astropy.coordinates import cartesian_to_spherical
import gala.coordinates as gc
import pandas as pd

#datos de la LMC Sacados de van der Marel, ApJL 832, L23 (2016)
#Tabla 2 columna (2) PMSTGAS
def CG_datos(cte):
    #DATOS
    ra_LMC=78.76*np.pi/180
    dec_LMC=-69.19*np.pi/180
    D_LMC=50.1 #en kpc
    vl_LMC=262.2 #en km/s
    mual_LMC=1.872 #en mas/yr
    mude_LMC=0.224#en mas/yr
    #Datosde SMC
    ra_SMC=12.80*np.pi/180
    dec_SMC=-73.15*np.pi/180
    D_SMC=62.80 #en kpc
    vl_SMC=145.6   #en km/s
    mual_SMC=0.79 #en mas/yr
    mude_SMC=-1.256 #en mas/yr
    
    #Seteamos los valores para el sistema galactico ya que los que vienen por default estan desactualizados
    n_frame = coorde.ICRS()
    v_sun = [-12.9, 245.6, 7.78] * (u.km / u.s)  # [vx, vy, vz]
    gc_frame = coorde.Galactocentric(galcen_distance=8.122*u.kpc,galcen_v_sun=v_sun,z_sun=20.8*u.pc)
    
    #inicializamos el objeto skycoord y hacemos la corrección por el mov solar de la LMC
    LMC_car=SkyCoord(ra=ra_LMC*u.rad,dec=dec_LMC*u.rad,distance=D_LMC*u.kpc,frame='icrs')
    c_LMC=coorde.SkyCoord(ra=ra_LMC*u.rad,dec=dec_LMC*u.rad,distance=D_LMC*u.kpc,pm_ra_cosdec=mual_LMC*u.mas/u.yr,pm_dec=mude_LMC*u.mas/u.yr,radial_velocity=vl_LMC*u.km/u.s,frame='icrs')
    corr_LMC=gc.reflex_correct(c_LMC,galactocentric_frame = gc_frame)
    #inicializamos el objeto skycoord y hacemos la corrección por el mov solar de la SMC
    SMC_car=SkyCoord(ra=ra_SMC*u.rad,dec=dec_SMC*u.rad,distance=D_SMC*u.kpc,frame='icrs')
    c_SMC=coorde.SkyCoord(ra=ra_SMC*u.rad,dec=dec_SMC*u.rad,distance=D_SMC*u.kpc,pm_ra_cosdec=mual_SMC*u.mas/u.yr,pm_dec=mude_SMC*u.mas/u.yr,radial_velocity=vl_SMC*u.km/u.s,frame='icrs')
    corr_SMC=gc.reflex_correct(c_SMC,galactocentric_frame = gc_frame)
    
    vcmx=(cte*corr_LMC.velocity.d_x.value+corr_SMC.velocity.d_x.value)/(cte+1)
    vcmy=(cte*corr_LMC.velocity.d_y.value+corr_SMC.velocity.d_y.value)/(cte+1)
    vcmz=(cte*corr_LMC.velocity.d_z.value+corr_SMC.velocity.d_z.value)/(cte+1)
    vsh=np.sqrt(vcmx**2+vcmy**2+vcmz**2)
    print('salida de CG')
    print('velocidad del CM',vcmx*u.km/u.s,vcmy*u.km/u.s,vcmz*u.km/u.s)
    print('módulo de velocidad',vsh*u.km/u.s)
    
    #obtengo el primer ángulo para rotar
    ang=cartesian_to_spherical(vcmx,vcmy,vcmz)
    t1=np.arctan2(vcmy,vcmx)
    t2=np.arctan2(vcmz,np.sqrt(vcmx**2+vcmy**2))
    vcm=np.sqrt(vcmx**2+vcmy**2+vcmz**2)
    print('chequeo de velocidad con ángulos',vcmx,vcmy,vcmz)
    print(vcm*np.cos(t2)*np.cos(t1),vcm*np.cos(t2)*np.sin(t1),vcm*np.sin(t2))
    print('chequeo primer ángulos',ang)# chequeo de rotaciones
    print(vcm,t2,t1)
    print('chequeo primer rotación',np.cos(t1)*np.cos(t2)*vcmx+np.sin(t1)*np.cos(t2)*vcmy+np.sin(t2)*vcmz,-np.sin(t1)*vcmx+np.cos(t1)*vcmy,-np.cos(t1)*np.sin(t2)*vcmx-np.sin(t1)*np.sin(t2)*vcmy+np.cos(t2)*vcmz)
    rcm=(cte*LMC_car.cartesian+SMC_car.cartesian)/(cte+1)
    #hago la primer roto-translacion
    CG=SkyCoord(l=0*u.deg,b=0*u.deg,frame='galactic')
    CG1=SkyCoord(ra=CG.icrs.ra.deg*u.deg,dec=CG.icrs.dec.deg*u.deg,distance=8.122*u.kpc,frame='icrs')
    LMCboost_cg=CG1.cartesian-rcm
    LMCcg_x=np.cos(t1)*np.cos(t2)*LMCboost_cg.x.value+np.sin(t1)*np.cos(t2)*LMCboost_cg.y.value+np.sin(t2)*LMCboost_cg.z.value
    LMCcg_y=-np.sin(t1)*LMCboost_cg.x.value+np.cos(t1)*LMCboost_cg.y.value
    LMCcg_z=-np.cos(t1)*np.sin(t2)*LMCboost_cg.x.value-np.sin(t1)*np.sin(t2)*LMCboost_cg.y.value+np.cos(t2)*LMCboost_cg.z.value

    return vsh,[LMCcg_x,LMCcg_y,LMCcg_z]
