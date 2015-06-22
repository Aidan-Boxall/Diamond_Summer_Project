

from Crystal import *

import copy

class MagneticStructure(Structure):
    def __init__(self,parent,site_labels=None,site_structure=None,space_group=None,magnetic_space_group=None):
        self._site = {}    # intended for access from outside the class with the setter/getter method self.site
        self._general_space_group = None    # Can be magnetic or crystallographic space group symbol
        self._magnetic_space_group_primings = None
        if site_labels is not None or site_structure is not None:
            self.site(site_labels=site_labels,site_structure=site_structure)
        self.parent = parent
    def site(self,label=None,coordinates=None,atom_type=None,occupancy=None,u_isotropic=None,u_anisotropic=None,magnetic_moment=None):
        if label is not None:
            Structure.site(self,label,coordinates,atom_type,occupancy,u_isotropic,u_anisotropic)
            if not hasattr(label,'__iter__'):
                label = [label]
                if magnetic_moment is not None:
                    magnetic_moment = [magnetic_moment]
            for jj in range(0,len(label)):
                if magnetic_moment is not None:
                    self._site[label[jj]]['magnetic_moment'] = magnetic_moment[jj]
                else:
                    self._site[label[jj]]['magnetic_moment'] = None
        else:
            return self._site
    def magnetic_space_group(self,magnetic_space_group=None):
        if magnetic_space_group is not None:
            self._general_space_group = magnetic_space_group
        else:
            if "'" in self._general_space_group:
                return self._general_space_group
            else:
                return None
    def space_group(self,space_group=None,format='symbol'):
        if space_group is not None:
            sgs = space_group_symbols(space_group)
            self._general_space_group = sgs.hermann_mauguin()
        else:
            msg = self._general_space_group
            sgsymbol = msg.replace("'","")
            sgs = space_group_symbols(sgsymbol)
            if format=='number':
                return sgs.number()
            elif format=='symbol':
                return sgs.hermann_mauguin()
    def magnetic_space_group_primings(self,magnetic_space_group_primings=None):
        # Must be ordered like the space_group_operators
        if magnetic_space_group_primings is not None:
            self._magnetic_space_group_primings = magnetic_space_group_primings
        else:
            return self._magnetic_space_group_primings
    def magnetic_space_group_operators(self):
        mops = []
        for jj in range(0,len(self.space_group_operators())):
            mops.append([self.space_group_operators()[jj][0],self.space_group_operators()[jj][1],self.magnetic_space_group_primings()[jj]])
        return mops
    def site_magnetic_point_group_operators(self,site_label):
        cpg_ops = self.site_point_group_operators(site_label)
        msg_ops = self.magnetic_space_group_operators()
        mpg_ops = []
        for jj in range(0,len(cpg_ops)):
            kk_skip = False
            for kk in range(0,len(msg_ops)):
                if (cpg_ops[jj]==msg_ops[kk][0]).all() and kk_skip is False:
                    mpg_ops.append([cpg_ops[jj],self.magnetic_space_group_primings()[kk]])
                    kk_skip = True
        return mpg_ops
    def site_equivalent_magnetic_operators(self,site_label):
        # Not univocal choice
        equiv_coord_oper = []
        coordinates = self.site_equivalent_coordinates(site_label)
        in_coord = self.site_coordinates_exact(site_label)
        magnetic_operators = self.magnetic_space_group_operators()
        transf_coordinates = []
        for mop in magnetic_operators:
            transf_coordinates.append(self.apply_to_coordinates(mop,in_coord))
        for coord in coordinates:
            distances = []
            for transf_coord in transf_coordinates:
                dist = np.linalg.norm([np.abs(coord[0]-transf_coord[0])%1,
                                       np.abs(coord[1]-transf_coord[1])%1,
                                       np.abs(coord[2]-transf_coord[2])%1])
                distances.append(dist)
            jj_min = np.argmin(distances)
            equiv_coord_oper.append([coord,magnetic_operators[jj_min]])
        return equiv_coord_oper
    def apply_to_coordinates(self,operator,coords):
        if len(operator)>1 and len(operator[1])==3:
            # operator = space group operator or magnetic space group operator
            
            new_coords = np.dot(operator[0],np.array(coords)) + np.array(operator[1])
        else:
            # operator = point group operator or magnetic point group operator
            new_coords = coords 
        return new_coords
    def apply_to_magnetic_moment(self,magnetic_operator,magnetic_moment):
        if len(magnetic_operator)==3:
            # operator = magnetic space group operator
            priming = magnetic_operator[2]
        elif len(magnetic_operator)==2:
            # operator = magnetic point group operator
            priming = magnetic_operator[1]
        rot = magnetic_operator[0]
        new_magnetic_moment = priming*np.linalg.det(rot)*np.dot(rot,np.array(magnetic_moment))
        return list(new_magnetic_moment)


class ParentPhase:
    def __init__(self,parent,coord_to_pp_operator=None,space_group=None):
        self._coord_to_pp_operator = coord_to_pp_operator
        self._space_group = space_group
        self.parent = parent
    def coord_to_pp_operator(self,coord_to_pp_operator=None):
        if coord_to_pp_operator is not None:
            self._coord_to_pp_operator = coord_to_pp_operator
        else:
            return self._coord_to_pp_operator
    def space_group(self,space_group=None):
        if space_group is not None:
            sgs = space_group_symbols(space_group)
            self._space_group = sgs.universal_hermann_mauguin()
        else:
            return self._space_group
    def find_cell(self):
        transf_matrix = np.linalg.inv(self.coord_to_pp_operator()[0])
        mag_g_matrix = self.parent.lattice.g_matrix()
        pp_g_matrix = np.dot(np.dot(np.transpose(transf_matrix),mag_g_matrix),transf_matrix)
        a_cell = np.sqrt(pp_g_matrix[0,0])
        b_cell = np.sqrt(pp_g_matrix[1,1])
        c_cell = np.sqrt(pp_g_matrix[2,2])
        alpha_cell = np.degrees(np.arccos(np.average([pp_g_matrix[1,2],pp_g_matrix[2,1]])/(b_cell*c_cell)))
        beta_cell = np.degrees(np.arccos(np.average([pp_g_matrix[0,2],pp_g_matrix[2,0]])/(a_cell*c_cell)))
        gamma_cell = np.degrees(np.arccos(np.average([pp_g_matrix[0,1],pp_g_matrix[1,0]])/(a_cell*b_cell)))
        pp_cell = self.validate_cell([a_cell,b_cell,c_cell,alpha_cell,beta_cell,gamma_cell])
        return pp_cell
    def validate_cell(self,pp_cell):
        sg = cctbx.sgtbx.space_group_info(self.space_group()).group()
        cctbx_cell = cctbx.uctbx.unit_cell(pp_cell)
        validated_cctbx_cell = sg.average_unit_cell(cctbx_cell)
        return validated_cctbx_cell.parameters()
    def validate_coord_to_pp_operator(self):
        op = self.coord_to_pp_operator()[0]
        pp_lattice = Lattice(None,parameters=self.find_cell())
        pp_g_matrix = pp_lattice.g_matrix()
        mag_g_matrix = self.parent.lattice.g_matrix()
        for jj in range(0,100):
            new_op = np.dot(np.dot(np.linalg.inv(mag_g_matrix),np.linalg.inv(np.transpose(op))),pp_g_matrix)
            op = new_op
        return [new_op, self.coord_to_pp_operator()[1]]
    def find_atoms(self):
        validated = {}
        oper = self.validate_coord_to_pp_operator()
        pp_cell = self.find_cell()
        pp_space_group =self.space_group()
        for site_label in self.parent.structure.site().keys():
            mag_coords = self.parent.structure.site()[site_label]['coordinates']
            pp_coords = self.parent.structure.apply_to_coordinates(oper,mag_coords)
            pp_coords_rounded = [round(x,10) for x in pp_coords]
            frac_pp_coords = list(flex.fmod_positive(flex.double(pp_coords_rounded),1))
            partial_pp = Crystal(parameters=pp_cell,site_labels=site_label)
            partial_pp.structure.space_group(pp_space_group)
            partial_pp.structure.site()[site_label]['coordinates'] = frac_pp_coords
            distances = []
            for val_label in validated.keys():
                tmp_pp = Crystal(parameters=pp_cell,site_labels=val_label)
                tmp_pp.structure.space_group(pp_space_group)
                tmp_pp.structure.site()[val_label]['coordinates'] = validated[val_label]
                validated_coordinates = tmp_pp.structure.site_equivalent_coordinates(val_label)
                for vc in validated_coordinates:
                    distances.append(np.linalg.norm(np.array(partial_pp.structure.site_coordinates_exact(site_label))-np.array(vc)))
            if all(distances)>0.001:    # True for empty list
                validated[site_label] = frac_pp_coords
        return validated
    def get_parentphase_crystal(self):
        pp_cell = self.find_cell()
        pp_space_group = self.space_group()
        pp_atoms = self.find_atoms()
        pphase_crystal = Crystal(reference=self.parent.general.reference,
                                 temperature_range=self.parent.general.temperature_range,
                                 experiment_temperature=self.parent.general.experiment_temperature,
                                 cifstring=self.parent.general.cifstring,
                                 parameters=pp_cell)
        for site_label in pp_atoms.keys():
            pphase_crystal.structure.site(label=site_label,coordinates=pp_atoms[site_label],
                                          atom_type=self.parent.structure.site()[site_label]['atom_type'],
                                          occupancy=self.parent.structure.site()[site_label]['occupancy'],
                                          u_isotropic=self.parent.structure.site()[site_label]['u_isotropic'])
        pphase_crystal.structure.space_group(pp_space_group)
        return pphase_crystal
    def hkl_to_pp(self,hkl):
        if not hasattr(hkl[0],'__iter__'):
            hkl_list = [hkl]
        else:
            hkl_list = hkl
        pp_hkl_list = []
        transf_matrix = np.transpose(np.linalg.inv(self.coord_to_pp_operator()[0]))
        for refl in hkl_list:
            pp_hkl_list.append(np.dot(transf_matrix,np.array(refl)))
        if not hasattr(hkl[0],'__iter__'):
            pp_hkl = pp_hkl_list[0]
        else:
            pp_hkl = pp_hkl_list
        return pp_hkl
    def hkl_from_pp(self,hkl):
        if not hasattr(hkl[0],'__iter__'):
            hkl_list = [hkl]
        else:
            hkl_list = hkl
        pp_hkl_list = []
        transf_matrix = np.transpose(self.coord_to_pp_operator()[0])
        for refl in hkl_list:
            pp_hkl_list.append(np.dot(transf_matrix,np.array(refl)))
        if not hasattr(hkl[0],'__iter__'):
            pp_hkl = pp_hkl_list[0]
        else:
            pp_hkl = pp_hkl_list
        return pp_hkl


class MagneticSpacegroup(Crystal):
    def __init__(self,reference=None,temperature_range=None,experiment_temperature=None,parameters=None,site_labels=None,cifstring=None,
                 coord_to_pp_operator=None,pp_space_group=None):
        self.general = General(self,reference,temperature_range,experiment_temperature,cifstring)
        self.lattice = Lattice(self,parameters)
        self.structure = MagneticStructure(self,site_labels)
        self.parentphase = ParentPhase(self,coord_to_pp_operator,pp_space_group)
    def load_cif(self,filename,data_code=None):
        Crystal.load_cif(self,filename,data_code)
        cif_crystal = iotbx.cif.reader(file_path = filename)
        if data_code is not None:
            cif_dict = cif_crystal.model()[data_code]
        else:
            cif_dict = cif_crystal.model().values()[0]
        # magnetic moments info
        for jj in range(0,len(cif_dict['_atom_site_moment_label'])):
            self.structure.site()[cif_dict['_atom_site_moment_label'][jj]]['magnetic_moment'] = [iotbx.cif.builders.float_from_string(cif_dict['_atom_site_moment_crystalaxis_x'][jj]),
                                                                                                 iotbx.cif.builders.float_from_string(cif_dict['_atom_site_moment_crystalaxis_y'][jj]),
                                                                                                 iotbx.cif.builders.float_from_string(cif_dict['_atom_site_moment_crystalaxis_z'][jj])]
        # magnetic space group info
        self.structure.magnetic_space_group(cif_dict['_space_group.magn_name_BNS'])
        position_p = []
        priming_p = []
        position_t = []
        priming_t = []
        for jj in range(0,len(cif_dict['_space_group_symop.magn_id'])):
            trasf = cif_dict['_space_group_symop.magn_operation_xyz'][jj].split(',')
            position_p.append(trasf[0:3])
            priming_p.append(float(trasf[3]))
        operators_p = self.find_operators(position_p)
        for jj in range(0,len(cif_dict['_space_group_symop.magn_centering_id'])):
            trasf = cif_dict['_space_group_symop.magn_centering_xyz'][jj].split(',')
            position_t.append(trasf[0:3])
            priming_t.append(float(trasf[3]))
        operators_t = self.find_operators(position_t)
        operators = []
        primings = []
        for jj in range(0,len(operators_t)):
            for kk in range(0,len(operators_p)):
                operators.append(self.multiply_operators(operators_t[jj],operators_p[kk]))
                primings.append(priming_t[jj]*priming_p[kk])
        sorted_primings = []
        for jj in range(0,len(self.structure.space_group_operators())):
            sg_op = self.structure.space_group_operators()[jj]
            for kk in range(0,len(operators)):
                mg_op = operators[kk]
                if (sg_op[0]==mg_op[0]).all() and (sg_op[1]==mg_op[1]).all():
                    sorted_primings.append(primings[kk])
        self.structure.magnetic_space_group_primings(sorted_primings)
        # parent phase info
        # add handling of nonstandard settings!!
        pp_space_group = cif_dict['_parent_space_group.name_H-M']
        self.parentphase.space_group(pp_space_group)
        coord_to_pp_operator = self.find_coord_to_pp_operator(cif_dict['_magnetic_space_group.transform_from_parent_Pp_abc'])
        self.parentphase.coord_to_pp_operator(coord_to_pp_operator)
    def save(self,filename):
        jdata = {'reference':self.general.reference,'temperature_range':self.general.temperature_range,'experiment_temperature':self.general.experiment_temperature,
                 'parameters':self.lattice.parameters(),'site':self.structure.site(),
                 'cifstring':self.general.cifstring,'magnetic_space_group_primings':self.structure.magnetic_space_group_primings(),'general_space_group':self.structure._general_space_group}
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
        self.structure.magnetic_space_group(str(jdata['general_space_group']))
        self.structure.magnetic_space_group_primings(jdata['magnetic_space_group_primings'])
        sites = jdata['site']
        for lab in sites.keys():
            self.structure.site(label=str(lab),atom_type=str(sites[lab]['atom_type']),occupancy=sites[lab]['occupancy'],coordinates=sites[lab]['coordinates'],
                                u_isotropic=sites[lab]['u_isotropic'],u_anisotropic=sites[lab]['u_anisotropic'],magnetic_moment=sites[lab]['magnetic_moment'])
    def find_operators(self,positions):
        operators = []
        for ep in positions:
            rot = np.zeros((3,3))
            transl = np.zeros(3)
            for jj in range(0,3):
                pos = filter(None, ep[jj].replace('-','+-').split('+'))
                for kk in range(0,len(pos)):
                    if 'x' in pos[kk]:
                        rot[jj,0] = eval(pos[kk].replace('x','1'))
                    elif 'y' in pos[kk]:
                        rot[jj,1] = eval(pos[kk].replace('y','1'))
                    elif 'z' in pos[kk]:
                        rot[jj,2] = eval(pos[kk].replace('z','1'))
                    else:
                        transl[jj] = eval(pos[kk])
            oper = [rot,transl]
            operators.append(oper)
        return operators
    def multiply_operators(self,op_first,op_second):
        rot = np.dot(op_second[0],op_first[0])
        transl = np.dot(op_first[0],op_second[1]) + op_first[1]
        oper = [rot,transl]
        return oper
    def find_coord_to_pp_operator(self,coord_to_pp_string):
        rot_string, transl_string = coord_to_pp_string.split(';')
        rot_vectors_string = rot_string.split(',')
        rot = np.zeros((3,3))
        for jj in range(0,3):
            column_vector = filter(None, rot_vectors_string[jj].replace('-','+-').split('+'))
            for kk in range(0,len(column_vector)):
                    if 'a' in column_vector[kk]:
                        try:
                            rot[0,jj] = eval(column_vector[kk].replace('a',''))
                        except SyntaxError:
                            rot[0,jj] = eval(column_vector[kk].replace('a','1'))
                    elif 'b' in column_vector[kk]:
                        try:
                            rot[1,jj] = eval(column_vector[kk].replace('b',''))
                        except SyntaxError:
                            rot[1,jj] = eval(column_vector[kk].replace('b','1'))
                    elif 'c' in column_vector[kk]:
                        try:
                            rot[2,jj] = eval(column_vector[kk].replace('c',''))
                        except SyntaxError:
                            rot[2,jj] = eval(column_vector[kk].replace('c','1'))
        transl = [eval(transl_string.split(',')[0]),eval(transl_string.split(',')[1]),eval(transl_string.split(',')[2])]
        operator = [rot,transl]
        return operator
    def reflection_list(self,energy,anomalous_flag=False,sort='two_theta',print_list=True,azimuth_ref=None,azimuth_angle=0,polin=[1,0],polout=[0,1]):
        pp_crystal = self.parentphase.get_parentphase_crystal()
        pp_bragg_list = pp_crystal.reflection_list(energy,refl='all',anomalous_flag=anomalous_flag,sort=sort,print_list=False)
        bragg_list = Crystal.reflection_list(self,energy,refl='all',anomalous_flag=False,sort='two_theta',print_list=False)
        # Everything is consistent except for:
        # 1) bragg_list has stronger intensity than pp_bragg_list because more atoms
        #     e.g. 10 atoms rather than 1 -> 10**2 times the structure factor
        # 2) bragg_list has reflections (not in pp_bragg_list) that are not exactly 0:
        #     because the real structure is slightly distorted OR because the space group unprimed is not the real space group??
#         Crystal.print_list(pp_crystal,pp_bragg_list)
#         Crystal.print_list(self,bragg_list)
        for row in bragg_list:
            row.insert(1,self.parentphase.hkl_to_pp(row[0]))
        mag_list = copy.deepcopy(bragg_list)
        for row in mag_list:
            refl = row[0]
            if azimuth_ref is None:
                mag_ampl = np.linalg.norm(self.magnetic_structure_factor(refl))
            else:
                mag_ampl = self.magnetic_amplitude(refl,energy,azimuth_ref=azimuth_ref,azimuth_angle=azimuth_angle,
                                                   polin=polin,polout=polout)
            mag_inten = np.abs(mag_ampl)**2
            row.append(mag_ampl)
            row.append(mag_inten)
        if print_list is True:
            self.print_list(mag_list)
        return mag_list
    def print_list(self,reflist):
        sys.stdout.write("{:<5}{:<18}{:<22}{:<22}{:<22}{:<22}{:<30}{:<20}\n".format('#','(h,k,l)','(h,k,l) parent','Intensity','Norm. Int.','Mag. Int.','TwoTheta (deg)','d_spacing (A)'))
        for jj in range(0,len(reflist)):
            #print jj, reflist[jj][0], '  ', reflist[jj][1], '  ', reflist[jj][2], '  ', reflist[jj][3], '  ', reflist[jj][4]
            sys.stdout.write("{:<5}{:<18}{:<22}{:<22}{:<22}{:<22}{:<30}{:<20}\n".format(jj+1,reflist[jj][0],reflist[jj][1],reflist[jj][2],reflist[jj][3],reflist[jj][8],reflist[jj][5],reflist[jj][6]))
    def magnetic_amplitude(self,refl,energy,azimuth_ref,azimuth_angle=0,polin=[1,0],polout=[0,1],process='e1e1',lsratio=1):
        polin = np.append(np.array(polin),0)
        polout = np.append(np.array(polout),0)
        # UAziCurrent
        e1 = np.dot(self.lattice.blume_polin_matrix(refl,energy),polin)
        e2 = np.dot(self.lattice.blume_polout_matrix(refl,energy),polout)
        e1.shape = (3,)
        e1.shape = (3,)
        if process=='e1e1':
            polfactor = -1j*np.cross(e2,e1)
            # UAziCurrent to UAziRef
            azimuth_matrix = self.lattice.azimuth_matrix(azimuth_angle)
            polfactor = np.dot(azimuth_matrix,polfactor)
            # UAziRef to CartRecipCrys
            blume_matrix = self.lattice.blume_matrix(refl,azimuth_ref)
            polfactor = np.dot(blume_matrix,polfactor)
        elif process=='nonresonant':
            theta = np.radians(self.lattice.theta(refl,energy))
            # UAziCurrent
            k1 = np.dot(self.lattice.blume_polin_matrix(refl,energy),np.array([0,0,1]))
            k2 = np.dot(self.lattice.blume_polout_matrix(refl,energy),np.array([0,0,1]))
            k1.shape = (3,)
            k2.shape = (3,)
            polfactor_s = np.cross(e2,e1) + \
                          np.dot(k2,e1)*np.cross(k2,e2) - \
                          np.dot(k1,e2)*np.cross(k1,e1) - \
                          np.cross(np.cross(k2,e2),np.cross(k1,e1))
            K = np.array([0,0,1])
            K.shape = (3,)
            polfactor_l = 1/2*(-4)*(np.sin(theta))**2*\
                          (np.cross(e2,e1) - \
                          np.dot(K,np.cross(e2,e1))*K)
            # UAziCurrent to UAziRef
            azimuth_matrix = self.lattice.azimuth_matrix(azimuth_angle)
            polfactor_s = np.dot(azimuth_matrix,polfactor_s)
            polfactor_l = np.dot(azimuth_matrix,polfactor_l)
            # UAziRef to CartRecipCrys
            blume_matrix = self.lattice.blume_matrix(refl,azimuth_ref)
            polfactor_s = np.dot(blume_matrix,polfactor_s)
            polfactor_l = np.dot(blume_matrix,polfactor_l)
        # NormDirectCrys
        magfactor = self.magnetic_structure_factor(refl)
        # NormDirectCrys to CartDirectCrys1
        direct_matrix = self.lattice.direct_matrix()
        magfactor = np.dot(direct_matrix,magfactor)
        # CartDirectCrys1 to CartDirectCrys2 <=> CartRecCrys
        cartesian_matrix = self.lattice.cartesian_matrix()
        magfactor = np.dot(cartesian_matrix,magfactor)
        # CartRecCrys <=> CartDirectCrys2
        if process=='e1e1':
            ampl = np.dot(magfactor,polfactor)
        elif process=='nonresonant':
            if lsratio == 0:
                l = [0,0,0]
            else:
                l = magfactor/(1+2/lsratio)
            s = magfactor/(2+lsratio)
            ampl = np.dot(l,polfactor_l) + np.dot(s,polfactor_s)
        return ampl
    def magnetic_structure_factor(self,refl):
        # NormDirectCrys
        refl = np.array(refl)
        magfactor = 0
        for site_label in self.structure.site().keys():
            if self.structure.site()[site_label]['magnetic_moment'] is not None:
                equiv_coord_oper = self.structure.site_equivalent_magnetic_operators(site_label)
                in_magnetic_moment = self.structure.site()[site_label]['magnetic_moment']
                for coord_oper in equiv_coord_oper:
                    magnetic_operator = coord_oper[1]
                    magnetic_moment = np.array(self.structure.apply_to_magnetic_moment(magnetic_operator,in_magnetic_moment))
                    coordinates = coord_oper[0]
                    magfactor = magfactor + magnetic_moment*np.exp(1j*2*np.pi*np.dot(coordinates,refl))
        return magfactor


