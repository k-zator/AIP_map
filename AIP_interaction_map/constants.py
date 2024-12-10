import numpy as np
import pandas as pd
import pathlib

PATH = pathlib.Path(__file__).parents[0]
FF_VS_PATH = ((PATH / "data" / "protein_aip_ff.xml").absolute().as_posix())
FF_VS_PATH_DUAL = ((PATH / "data" / "protein_aip_ff_dual.xml").absolute().as_posix())
CML_NS = "http://www.xml-cml.org/schema"
SSIP_NS = "http://www-hunter.ch.cam.ac.uk/SSIP"
LIG = "MOL"
AA_LIST = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS',
           'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ASX', 'GLH', 'CSS']

LIG = "MOL"
TEMPERATURE = 298.0
CONVERSION_FROM_J_TO_KJ = 0.001
GAS_CONSTANT = 8.3144598 #kJK^-1mol^-1,
#value of the van der Waals interaction of 2 SSIPs interacting in a
#pairwise manner. As detailed in DOI: http://dx.doi.org/10.1039/c3sc22124e.
EVDW = -5.6 #kJmol^-1.
#value of maximum concentration of SSIPs in a phase.
CMAX = 300.0

POLY_PATH = ((PATH / "data" / "SSIMPLEpolynomials.csv").absolute().as_posix())
POLY_COEFFS = {}
df = pd.read_csv(POLY_PATH, header=0, sep="\t")
for i, row in df.iterrows():
    POLY_COEFFS[row[0]] = [row[2], {'positive': np.array(row[12:21]), 'negative': np.array(row[3:12])}] 


poly_coeffs_wat = {'positive': np.array([0.838237, -0.51907000, -7.714317E-02, -8.1380530E-01, 3.524312E-01, -7.0570040E-02, 7.6254820E-03, -4.2989670E-04, 9.9333680E-06]),
                   'negative': np.array([0.7751790, 1.25897301E-02, -3.7236050E-01, 4.1098470E-02, 1.7334470E-02, 2.097013000E-03, 1.2569490E-04, 3.79972E-06, 4.62349E-08])}
poly_coeffs_hex = {'positive': np.array([-1.257667, -1.343249E-02, -1.0086660E-03, 2.992371E-06, 9.5324600E-05, -2.9965860E-05, 4.0427420E-06, -2.6054980E-07, -6.573138E-09]),
                   'negative': np.array([-1.256501, 1.30147200E-01, -3.56959E-03, 2.6307740E-04, 3.59959700E-05, 2.771214000E-06, 1.3473240E-07, 3.622555E-09, 4.0709470E-11])}
poly_coeffs_chl = {'positive': np.array([-1.220126, -1.662607E-01, -9.94286E-03, -4.672309E-04, -2.678255E-06, 1.475957E-05, -2.411077E-06, 1.74012E-07, -4.822481E-09]),
                   'negative': np.array([-1.187073, 4.507678E-01, 1.42246E-01, 9.4427E-02, 1.417198E-02, 1.072633E-03, 4.489395E-05, 9.877984E-07, 8.891901E-09])}
poly_coeffs_cs2 = {'positive': np.array([-1.21343295, -2.49821200E-14, 2.6091539E-14, -1.3373634E-14, 3.8865428E-15, -6.7489199E-16, 6.9181780E-17, -3.8478277E-18, 8.9347593E-20]),
                   'negative': np.array([-1.21293227, 1.9037009E-01, -1.6858256E-02, 7.7350536E-04, 2.6699505E-05, -6.8706906E-06, -5.9836524E-07, -1.9282643E-08, -2.2981511E-10])}
poly_coeffs_dcm = {'positive': np.array([-1.10099626, -0.29869698, -4.6780011E-02, -3.1416267E-03, -2.39800E-04, 1.544126E-04, -1.6726002E-05, 7.59042177E-07, -1.2701268E-08]),
                   'negative': np.array([-1.08795265, 0.30898577, 5.2884521E-02, 7.8507850E-02, 1.4260985E-02, 1.2737424E-03, 6.2943009E-05, 1.6485356E-06, 1.7901279E-08])}
poly_coeffs_toluene = {'positive': np.array([-1.063877347, -0.48268113, -0.10558987, -4.7627433E-02, 1.1666038E-02, -8.867099E-04, -9.2032543E-06, 4.5104983E-06, -1.6817837E-07]),
                       'negative': np.array([-1.064333460, 0.29451487, -2.8469254E-02, 2.8269265E-03, 4.799617E-04, 2.8344249E-05, 7.84524129E-07, 7.9749664E-09, -1.9539591E-11])}
poly_coeffs_mesitylene = {'positive': np.array([-1.103849090, -0.37903855, -8.751164E-02, -4.6345369E-02, 5.2364119E-03, 9.862940E-04, -2.466683E-04, 1.90397017E-05, -5.1946617E-07]),
                          'negative': np.array([-1.104818021, 0.30446835, -2.4126399E-02, 1.2238449E-03, 2.552200E-04, 1.5827083E-05, 4.8391690E-07, 6.8743554E-09, 2.8392763E-11])}
poly_coeffs_cyclohexane = {'positive': np.array([-1.20321812, 1.1656654E-14, -1.6844975E-14, 1.032598E-14, -3.3840555E-15, 6.3696204E-16, -6.9015864E-17, 3.9989851E-18, -9.5811505E-20]),
                           'negative': np.array([-1.20204813, 2.8445054E-01, -1.0259071E-02, 2.6637015E-04, 8.0451629E-05, 6.5457825E-06, 3.0242337E-07, 7.7914070E-09, 8.6030397E-11])}
poly_coeffs_acetonitrile = {'positive': np.array([-0.6769533416, -1.1734686, 1.0260458, -0.96749574, 0.24525482, -2.9545913E-02, 1.6865190E-03, -2.9952386E-05, -5.0541168E-07]),
                            'negative': np.array([-0.7280925810, 0.34542032, -3.9179104E-02, 4.141877E-02, 8.5852890E-03, 8.078711E-04, 4.1276042E-05, 1.109405E-06, 1.2313114E-08])}
poly_coeffs_noble = {'positive': np.array([-1.2566373, -6.87104537E-15, 2.582209E-15, 1.7683018E-16, -2.796134E-16, 5.9802896E-17, -5.453352E-18, 2.1687955E-19, -2.6261837E-21]),
                     'negative': np.array([-1.2566373, 1.30168203E-14, 1.101443E-14, 3.9705275E-15, 7.39536634E-16, 7.710874E-17, 4.5474337E-18, 1.41760874E-19, 1.8151619E-21])}

theta_wat = 0.73848929
theta_hex = 0.41933333
theta_chl = 0.45466667
theta_cs2 = 0.35943116
theta_dcm = 0.369778336
theta_toluene = 0.45665283
theta_mesitylene = 0.4483533
theta_cyclohexane = 0.34772212
theta_acetonitrile = 0.4768397
theta_noble = 0.41491364 #i.e. perfectly neutral chloroform footprint
