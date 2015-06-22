
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import re
from operator import itemgetter
import sys

import iotbx.cif
import iotbx.cif.builders
from iotbx.cif.builders import CifBuilderError
import cctbx.eltbx.xray_scattering
import cctbx.eltbx.sasaki
import cctbx.sgtbx
from cctbx.sgtbx import space_group, space_group_symbols
from cctbx.array_family import flex

import json

class General:
    def __init__(self,parent,reference=None,temperature_range=None,experiment_temperature=None,cifstring=None):
        self.reference = reference
        self.temperature_range = temperature_range
        self.experiment_temperature = experiment_temperature
        self.cifstring = cifstring
        self.parent = parent

class Lattice:
    def __init__(self,parent,parameters=None):
        self._parameters = parameters      # intended for access from outside the class with the setter/getter method self.parameters
        self.parent = parent
    def parameters(self,parameters=None):
        # As a setter method: e.g. mycrys.lattice.parameters([5,5,5,90,90,90])
        # As a getter method: e.g. params = mycrys.lattice.parameters()
        if parameters is not None:
            parameters = list(parameters)
            parameters = map(float,parameters)
            self._parameters = tuple(parameters)
        else:
            return self._parameters
    def reciprocal_parameters(self):
        cell_volume = mt.sqrt(np.linalg.det(self.g_matrix()))
        lp = list(self.parameters())    # casting as list creates a deepcopy
        lp[3:] = map(mt.radians,lp[3:])
        rp = []
        rp.append(lp[1]*lp[2]*np.sin(lp[3])/cell_volume)
        rp.append(lp[2]*lp[0]*np.sin(lp[4])/cell_volume)
        rp.append(lp[0]*lp[1]*np.sin(lp[5])/cell_volume)
        rp.append(np.arccos( (np.cos(lp[4])*np.cos(lp[5])-np.cos(lp[3])) / (np.sin(lp[4])*np.sin(lp[5])) ))
        rp.append(np.arccos( (np.cos(lp[3])*np.cos(lp[5])-np.cos(lp[4])) / (np.sin(lp[3])*np.sin(lp[5])) ))
        rp.append(np.arccos( (np.cos(lp[3])*np.cos(lp[4])-np.cos(lp[5])) / (np.sin(lp[3])*np.sin(lp[4])) ))
        rp[3:] = map(mt.degrees,rp[3:])
        rp[3:] = map(round,rp[3:])
        return list(rp)
    def b_matrix(self):
        # RecipCrys to CartRecipCrys
        lp = list(self.parameters())    # casting as list creates a deepcopy
        lp[3:] = map(mt.radians,lp[3:])
        rp = list(self.reciprocal_parameters())
        rp[3:] = map(mt.radians,rp[3:])
        b_matrix = np.array([[ rp[0], rp[1]*np.cos(rp[5]), rp[2]*np.cos(rp[4]) ], 
                             [ 0, rp[1]*np.sin(rp[5]), -rp[2]*np.sin(rp[4])*np.cos(lp[3]) ], 
                             [ 0, 0, 1/lp[2] ]])
        return b_matrix
    def g_matrix(self):
        lp = list(self.parameters())    # casting as list creates a deepcopy
        lp[3:] = map(mt.radians,lp[3:])
        g_matrix = np.array([[ lp[0]**2, lp[0]*lp[1]*np.cos(lp[5]), lp[0]*lp[2]*np.cos(lp[4]) ],
                             [ lp[1]*lp[0]*np.cos(lp[5]), lp[1]**2, lp[1]*lp[2]*np.cos(lp[3]) ],
                             [ lp[2]*lp[0]*np.cos(lp[4]), lp[2]*lp[1]*np.cos(lp[3]), lp[2]**2 ]])
        return g_matrix
    def g_star_matrix(self):
        rp = self.reciprocal_parameters()
        rp[3:] = map(mt.radians,rp[3:])
        g_star_matrix = np.array([[ rp[0]**2, rp[0]*rp[1]*np.cos(rp[5]), rp[0]*rp[2]*np.cos(rp[4]) ],
                                  [ rp[1]*rp[0]*np.cos(rp[5]), rp[1]**2, rp[1]*rp[2]*np.cos(rp[3]) ],
                                  [ rp[2]*rp[0]*np.cos(rp[4]), rp[2]*rp[1]*np.cos(rp[3]), rp[2]**2 ]])
        return g_star_matrix
    def d_spacing(self,refl):
        refl_vec = np.array(refl)
        g_star_matrix = self.g_star_matrix()
        refl_norm = np.sqrt(np.dot(refl_vec.transpose(),np.dot(g_star_matrix,refl_vec)))
        d_spacing = 1/refl_norm
        return d_spacing
    def theta(self,refl,energy):
        d_spacing = self.d_spacing(refl)
        sth = (12.3984/energy)/(2*d_spacing)
        if 0 <= sth <= 1:
            theta = np.degrees(np.arcsin(sth))
        else:
            theta = None
            print 'reflection '+str(refl)+' not reachable at Energy = '+str(energy)+' keV'
        return theta
    def blume_matrix(self,refl,azimuth_ref):
        # UAziRef to CartRecipCrys
        # Crystallographic convention: reflection // -u3
        # Vector No 1
        # RecipCrys
        refl_vec = np.array(refl)
        g_star_matrix = self.g_star_matrix()
        refl_norm = np.sqrt(np.dot(refl_vec.transpose(),np.dot(g_star_matrix,refl_vec)))
        refl_normalised = refl_vec/refl_norm
        # RecipCrys to CartRecipCrys
        b_matrix = self.b_matrix()
        refl_normalised_transformed = np.dot(b_matrix,refl_normalised)
        # Vector No 2, not perpendicular
        # RecipCrys
        azi_vec = np.array(azimuth_ref)
        # RecipCrys to CartRecipCrys
        azi_transformed = np.dot(b_matrix,azi_vec)
        # Vector No 2, perpendicular
        azi_perp_transformed = azi_transformed - np.dot(refl_normalised_transformed.transpose(),azi_transformed)*refl_normalised_transformed
        azi_perp__normalised_transformed = azi_perp_transformed/np.linalg.norm(azi_perp_transformed)
        # Vector No 3
        v1 = refl_normalised_transformed
        v2 = azi_perp__normalised_transformed
        v1.shape = (3,)
        v2.shape = (3,)
        v3 = np.cross(v1,v2)
        blume_matrix = np.transpose(np.vstack((v2,
                                               -v3,
                                               -v1)))
        return blume_matrix
    def azimuth_matrix(self,azimuth_angle):
        # UAziCurrent_to_UAziRef
        azimuth_matrix = np.array([[np.cos(azimuth_angle),-np.sin(azimuth_angle),0],
                                   [np.sin(azimuth_angle),np.cos(azimuth_angle),0],
                                   [0,0,1]])
        return azimuth_matrix
    def blume_polin_matrix(self,refl,energy):
        theta = np.radians(self.theta(refl,energy))
        blume_polin_matrix = np.array([[0,np.sin(theta),np.cos(theta)],[-1,0,0],[0,-np.cos(theta),np.sin(theta)]])
        return blume_polin_matrix
    def blume_polout_matrix(self,refl,energy):
        theta = np.radians(self.theta(refl,energy))
        blume_polout_matrix = np.array([[0,-np.sin(theta),np.cos(theta)],[-1,0,0],[0,-np.cos(theta),-np.sin(theta)]])
        return blume_polout_matrix
    def direct_matrix(self):
        # NormDirectCrys to CartDirectCrys1
        lp_norm = [1,1,1] + list(self.parameters()[3:6])
        lattice_norm = Lattice(None,parameters=lp_norm)
        rp_norm = list(lattice_norm.reciprocal_parameters())
        direct_matrix = np.array([[ lp_norm[0], lp_norm[1]*np.cos(np.radians(lp_norm[5])), lp_norm[2]*np.cos(np.radians(lp_norm[4])) ],
                                    [ 0,        lp_norm[1]*np.sin(np.radians(lp_norm[5])), -lp_norm[2]*np.sin(np.radians(lp_norm[4]))*np.cos(np.radians(rp_norm[3])) ],
                                    [ 0,        0,                                         1/rp_norm[2] ]])
        return direct_matrix
    def cartesian_matrix(self):
        direct_matrix = self.direct_matrix()
        # CartDirectCrys1
        b = np.dot(direct_matrix,np.array([0,1,0]))
        c = np.dot(direct_matrix,np.array([0,0,1]))
        # CartDirectCrys1 to CartDirectCrys2 <-> CartRecCrys
        a_star_cart = np.cross(b,c)/np.linalg.norm(np.cross(b,c))
        b_star_cart = np.cross(c,np.cross(b,c))/np.linalg.norm(np.cross(c,np.cross(b,c)))
        c_star_cart = np.cross(a_star_cart,b_star_cart)
        cartesian_matrix_inv = np.transpose(np.vstack((a_star_cart,
                                                       b_star_cart,
                                                       c_star_cart)))
        cartesian_matrix = np.linalg.inv(cartesian_matrix_inv)
        return cartesian_matrix


class Structure:
    def __init__(self,parent,site_labels=None,space_group=None):
        self._site = {}    # intended for access from outside the class with the setter/getter method self.site
        self._space_group = None    # intended for access from outside the class with the setter/getter method self.space_group
        if site_labels is not None:
            self.site(label=site_labels)
        if space_group is not None:
            self.space_group(space_group)
        self.parent = parent
    def site(self,label=None,coordinates=None,atom_type=None,occupancy=None,u_isotropic=None,u_anisotropic=None):
        if label is not None:
            if not hasattr(label,'__iter__'):
                label = [label]
                if coordinates is not None:
                    coordinates = [coordinates]
                if atom_type is not None:
                    atom_type = [atom_type]
                if occupancy is not None:
                    occupancy = [occupancy]
                if u_isotropic is not None:
                    u_isotropic = [u_isotropic]
                if u_anisotropic is not None:
                    u_anisotropic = [u_anisotropic]
            for jj in range(0,len(label)):
                if label[jj] not in self._site:
                    self._site[label[jj]] = {'coordinates':None,'occupancy':1,'atom_type':None,
                                              'u_isotropic':0,'u_anisotropic':0}
                if coordinates is not None:
                    self._site[label[jj]]['coordinates'] = coordinates[jj]
                if atom_type is not None:
                    self._site[label[jj]]['atom_type'] = atom_type[jj]
                if occupancy is not None:
                    self._site[label[jj]]['occupancy'] = occupancy[jj]
                if u_isotropic is not None:
                    self._site[label[jj]]['u_isotropic'] = u_isotropic[jj]
                if u_anisotropic is not None:
                    self._site[label[jj]]['u_anisotropic'] = u_anisotropic[jj]
        else:
            return self._site
    def space_group(self,space_group=None,format='symbol'):
        # format can be 'number', 'symbol' (Hermann-Mauguin)
        # As a setter method, no need to specify the format; always converts to symbol to store in self._space_group
        # and presumes standard setting if number is given
        # As a getter method, 'symbol' is the default
        if space_group is not None:
                sgs = space_group_symbols(space_group)
                self._space_group = sgs.universal_hermann_mauguin()
                # self._space_group = sgs.universal_hall() needs reading with sgs = space_group_symbols(self._space_group,table_id = "Hall")
        else:
            sgs = space_group_symbols(self._space_group)
            if format=='number':
                return sgs.number()
            elif format=='symbol':
                return sgs.universal_hermann_mauguin()
    def point_group(self):
        sg = cctbx.sgtbx.space_group_info(self.space_group()).group()      # reconstruct object space_group from number
        return sg.point_group_type()
    def equivalent_points(self):
        sg = cctbx.sgtbx.space_group_info(self.space_group()).group()      # reconstruct object space_group from number
        equivalent_points = []
        for ep in sg.info().group():
            equivalent_points.append(str(ep).split(','))
        return equivalent_points
    def space_group_operators(self):
        sg = cctbx.sgtbx.space_group_info(self.space_group()).group()      # reconstruct object space_group from number
        ops=[]
        for op in sg.all_ops():
            rot = np.array(op.r().as_double())
            rot.shape = (3,3)
            trasl = np.array(op.t().as_double())
            trasl.shape = (3,)
            ops.append([rot,trasl])
        return ops
    def point_group_operators(self):
        sg = cctbx.sgtbx.space_group_info(self.space_group()).group()      # reconstruct object space_group from number
        pg=sg.build_derived_point_group()
        ops=[]
        for op in pg.all_ops():
            rot = np.array(op.r().as_double())
            rot.shape = (3,3)
            ops.append(rot)
        return ops
    def site_coordinates_exact(self,site_label):
        locsym = self.parent.get_locsym(site_label)
        coord_exact = locsym.exact_site()
        return coord_exact
    def site_equivalent_coordinates(self,site_label):
        cctbx_crystal = self.parent.get_cctbx_structure()
        cctbx_crystal = cctbx_crystal.expand_to_p1(append_number_to_labels=False,sites_mod_positive=True)
        cctbx_atoms = cctbx_crystal.scatterers()
        eq_coord = []
        for jj in cctbx_atoms:
            if jj.label==site_label:
               eq_coord.append(jj.site)
        return eq_coord
    def site_symmetry_points(self,site_label):
        locsym = self.parent.get_locsym(site_label)
        symm_pos = []
        for ep in locsym.matrices():
            symm_pos.append(str(ep).split(','))
        return symm_pos
    def site_multiplicity(self,site_label):
        locsym = self.parent.get_locsym(site_label)
        return locsym.multiplicity()
    def site_wyckoff(self,site_label):
        locsym = self.parent.get_locsym(site_label)
        sg = cctbx.sgtbx.space_group_info(self.space_group()).group()      # reconstruct object space_group from number
        wyckoff_position = sg.info().wyckoff_table().mapping(locsym)
        wyckoff_letter = wyckoff_position.position().letter()
        return wyckoff_letter
    def site_point_group(self,site_label):
        locsym = self.parent.get_locsym(site_label)
        return locsym.point_group_type()
    def site_point_group_operators(self,site_label):
        locsym = self.parent.get_locsym(site_label)
        ops=[]
        for op in locsym.matrices():
            rot = np.array(op.r().as_double())
            rot.shape = (3,3)
            ops.append(rot)
        return ops


class Crystal:
    def __init__(self,reference=None,temperature_range=None,experiment_temperature=None,parameters=None,site_labels=None,cifstring=None):
        self.general = General(self,reference,temperature_range,experiment_temperature,cifstring)
        self.lattice = Lattice(self,parameters)
        self.structure = Structure(self,site_labels)
    def load_cif(self,filename,data_code=None):
        # Sets: 
        # cifstring;
        # lattice_parameters; 
        # space_group (= only the space group number - cctbx object to be reconstructed from that);
        # equivalent_positions;
        # for every site: label; coordinates; atom_type; occupancy; u_iso; equivalent_positions
        # 1) read file as string
        with open(filename, "r") as ciffile:
            #cifstring = ciffile.readlines()
            #cifstring = ciffile.read().replace('\n', '')
            cifstring = ciffile.read()
        self.general.cifstring = cifstring
        # 2) read file with cctbx
        # e.g. filename = 'D:\test.cif'; icsd_code = '151772-ICSD' from 1st uncommented line in cif file = data_151772-ICSD
        cif_crystal = iotbx.cif.reader(file_path = filename)
        try:
            cctbx_crystal = cif_crystal.build_crystal_structures()
            if data_code is not None:
                cctbx_crystal = cctbx_crystal[data_code]
            else:
                cctbx_crystal = cctbx_crystal.values()[0]
            lp = cctbx_crystal.crystal_symmetry().unit_cell().parameters()     # type = cctbc.uctbx.ext.unit_cell or tuple(lattice parameters)
            self.lattice.parameters(lp)
            sg = cctbx_crystal.crystal_symmetry().space_group()    # type = cctbx.sgtbx.space_group (defined in cctbx_build\lib\cctbx_sgtbx_ext.pyd - not visible)
            sgsymbol = sg.info().symbol_and_number().rpartition(' (No.')[0]
            self.structure.space_group(sgsymbol)
            cctbx_atoms = cctbx_crystal.scatterers()
            for jj in cctbx_atoms:
                self.structure.site(label=jj.label)
                self.structure.site()[jj.label]['coordinates'] = jj.site
                self.structure.site()[jj.label]['atom_type'] = jj.scattering_type
                self.structure.site()[jj.label]['occupancy'] = jj.occupancy
                self.structure.site()[jj.label]['u_isotropic'] = jj.u_iso     # cctbx automatically converts from B to U = B/(8*np.pi**2) - checked;
                                                                              # automatically converts 0 to -1 if anisotropic parameters are present
        except CifBuilderError:
            if data_code is not None:
                cif_dict = cif_crystal.model()[data_code]
            else:
                cif_dict = cif_crystal.model().values()[0]
            lp = [iotbx.cif.builders.float_from_string(cif_dict['_cell_length_a']),
                  iotbx.cif.builders.float_from_string(cif_dict['_cell_length_b']),
                  iotbx.cif.builders.float_from_string(cif_dict['_cell_length_c']),
                  iotbx.cif.builders.float_from_string(cif_dict['_cell_angle_alpha']),
                  iotbx.cif.builders.float_from_string(cif_dict['_cell_angle_beta']),
                  iotbx.cif.builders.float_from_string(cif_dict['_cell_angle_gamma'])]
            self.lattice.parameters(lp)
            try:
                msgsymbol = cif_dict['_space_group.magn_name_BNS']
                sgsymbol = msgsymbol.replace("'","")
                self.structure.space_group(sgsymbol)
            except KeyError:
                sgsymbol = cif_dict['_symmetry_space_group_name_H-M']
                self.structure.space_group(sgsymbol)
            for jj in range(0,len(cif_dict['_atom_site_label'])):
                site_label = cif_dict['_atom_site_label'][jj]
                self.structure.site(label=site_label)
                self.structure.site()[site_label]['coordinates'] = [iotbx.cif.builders.float_from_string(cif_dict['_atom_site_fract_x'][jj]),
                                                                    iotbx.cif.builders.float_from_string(cif_dict['_atom_site_fract_y'][jj]),
                                                                    iotbx.cif.builders.float_from_string(cif_dict['_atom_site_fract_z'][jj])]
                self.structure.site()[site_label]['atom_type'] = cif_dict['_atom_site_type_symbol'][jj]
                try:
                    self.structure.site()[site_label]['occupancy'] = iotbx.cif.builders.float_from_string(cif_dict['_atom_site_occupancy'][jj])
                except KeyError:
                    pass
                try:
                    self.structure.site()[site_label]['u_isotropic'] = iotbx.cif.builders.float_from_string(cif_dict['_atom_site_U_iso_or_equiv'][jj])     # ??? to check
                except KeyError:
                    pass
    def save(self,filename):
        jdata = {'reference':self.general.reference,'temperature_range':self.general.temperature_range,'experiment_temperature':self.general.experiment_temperature,
                 'parameters':self.lattice.parameters(),'site':self.structure.site(),'space_group':self.structure.space_group(),
                 'cifstring':self.general.cifstring}
        with open(filename,'w') as json_outfile:
            json.dump(jdata,json_outfile,sort_keys=True,indent=4,ensure_ascii=False)
    def load(self,filename):
        with open(filename) as json_inputfile:
            jdata = json.load(json_inputfile)
        self.general.reference = str(jdata['reference'])
        self.general.temperature_range = jdata['temperature_range']
        self.general.experiment_temperature = jdata['experiment_temperature']
        self.general.cifstring = str(jdata['cifstring'])
        self.lattice.parameters(jdata['parameters'])
        self.structure.space_group(str(jdata['space_group']))
        sites = jdata['site']
        for lab in sites.keys():
            self.structure.site(label=str(lab),atom_type=str(sites[lab]['atom_type']),occupancy=sites[lab]['occupancy'],coordinates=sites[lab]['coordinates'],
                                u_isotropic=sites[lab]['u_isotropic'],u_anisotropic=sites[lab]['u_anisotropic'])
    def save_cif(self,filename):
        with open(filename, "w") as ciffile:
            ciffile.write(self.general.cifstring)
    def reflection_list(self,energy,refl='sym',anomalous_flag=False,sort='two_theta',print_list=True):
        # refl = 'all'; 'allowed'; 'sym
        # sort = 'int'; 'two_theta'; 'd_spacing'; 'index'
        
        wavelength = 12.39842/energy
        d_min = wavelength/2
        cctbx_atoms = self.validate_cctbx_atoms()
        cctbx_xray = self.get_cctbx_structure(cctbx_atoms)
        cctbx_xray.scattering_type_registry(table="wk1995")
        if anomalous_flag is True:
            cctbx_xray.set_inelastic_form_factors(wavelength,'sasaki')
        if refl in ['all','allowed']:
            cctbx_xray = cctbx_xray.expand_to_p1(sites_mod_positive=True)
        scattering = cctbx_xray.structure_factors(anomalous_flag,algorithm='direct',d_min=d_min).f_calc()
        ind = scattering.indices()
        inten = (scattering.amplitudes().data())**2 #|F(hkl)|**2
        amp = list(scattering.data()) # F(hkl) complex
        two_theta = scattering.two_theta(wavelength,deg=True).data()
        d_spacing = scattering.d_spacings().data()
        reflist = []
        if refl=='all':
            sg = cctbx.sgtbx.space_group_info(self.structure.space_group()).group()
            for jj in range(0,len(ind)):
                if not sg.is_sys_absent(ind[jj]):
                    reflist.append([ind[jj], inten[jj], amp[jj], two_theta[jj], d_spacing[jj]])
                else:
                    reflist.append([ind[jj], 0, 0, two_theta[jj], d_spacing[jj]])
        elif refl=='allowed':
            sg = cctbx.sgtbx.space_group_info(self.structure.space_group()).group()
            for jj in range(0,len(ind)):
                if not sg.is_sys_absent(ind[jj]):
                    reflist.append([ind[jj], inten[jj], amp[jj], two_theta[jj], d_spacing[jj]])
        elif refl=='sym':
            for jj in range(0,len(ind)):
                reflist.append([ind[jj], inten[jj], amp[jj], two_theta[jj], d_spacing[jj]])
        maxint = max([row[1] for row in reflist])
        for row in reflist:
            row.insert(2,row[1]/maxint*100)
        if sort=='two_theta':
            item = 4
        elif sort=='d_spacing':
            item = 5
        elif sort=='int':
            item = 1
        elif sort=='index':
            item = 0
        reflist = sorted(reflist, key=itemgetter(item))
        if print_list is True:
            self.print_list(reflist)
        return reflist
    def print_list(self,reflist):
        sys.stdout.write("{:<5}{:<18}{:<20}{:<18}{:<40}{:<30}{:<20}\n".format('#','(h,k,l)','Intensity','Norm. Int.','Amplitude','TwoTheta (deg)','d_spacing (A)'))
        for jj in range(0,len(reflist)):
            #print jj, reflist[jj][0], '  ', reflist[jj][1], '  ', reflist[jj][2], '  ', reflist[jj][3], '  ', reflist[jj][4]
            sys.stdout.write("{:<5}{:<18}{:<20}{:<18}{:<40}{:<30}{:<20}\n".format(jj+1,reflist[jj][0],reflist[jj][1],reflist[jj][2],reflist[jj][3],reflist[jj][4],reflist[jj][5]))
    def get_cctbx_structure(self,cctbx_atoms=None,cctbx_symm=None):
        if cctbx_atoms is None:
            cctbx_atoms = []
            for site_label in self.structure.site().keys():
                cctbx_atoms.append(cctbx.xray.scatterer(label=str(site_label),
                                                        site=(self.structure.site()[site_label]['coordinates']),
                                                        u=self.structure.site()[site_label]['u_isotropic'],
                                                        occupancy=self.structure.site()[site_label]['occupancy'],
                                                        scattering_type=str(self.structure.site()[site_label]['atom_type'])))
        if cctbx_symm is None:
            cctbx_symm = cctbx.crystal.symmetry(unit_cell=self.lattice.parameters(),
                                                space_group_symbol=self.structure.space_group())
        settings = cctbx.crystal.special_position_settings(crystal_symmetry=cctbx_symm)
        cctbx_struct = cctbx.xray.structure(special_position_settings=settings,
                                            scatterers = flex.xray_scatterer(cctbx_atoms))
        return cctbx_struct
    def validate_cctbx_atoms(self):
        validated_atoms = {}
        for site_label in self.structure.site().keys():
            atom_type = self.structure.site()[site_label]['atom_type']
            try:
                atom_type_validated = cctbx.eltbx.xray_scattering.get_standard_label(atom_type,exact=True)
            except RuntimeError:
                atom_type_validated = cctbx.eltbx.xray_scattering.get_standard_label(atom_type,exact=False)
                print 'Warning: scattering atom type %s replaced with %s' %(atom_type, atom_type_validated)
            validated_atoms[site_label] = atom_type_validated
        cctbx_atoms = []
        for site_label in self.structure.site().keys():
            cctbx_atoms.append(cctbx.xray.scatterer(label=str(site_label),
                                                    site=(self.structure.site()[site_label]['coordinates']),
                                                    u=self.structure.site()[site_label]['u_isotropic'],
                                                    occupancy=self.structure.site()[site_label]['occupancy'],
                                                    scattering_type=validated_atoms[site_label]))
        return cctbx_atoms
    def get_locsym(self,site_label):
        cctbx_struct = self.get_cctbx_structure()
        locsym = cctbx_struct.site_symmetry(self.structure.site()[site_label]['coordinates'])
        return locsym

# class ReflectionList:
#     def __init__:
#         self.reflections = None
#     def printlist(self):
#         pass




