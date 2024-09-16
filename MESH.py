import numpy as np
from scipy.spatial import Delaunay
from scipy.sparse.linalg import spsolve
from scipy.linalg import solve

class MKE_solver:

    def compute(self, file, Delaunay_mesh=False, sparse=True, assume_a='pos',scale=1,circle='half'):
        K,b = self.meshing(file, Delaunay_mesh,scale)
        if sparse:
            u = spsolve(K,b)
        else:
            u = solve(K,b, assume_a=assume_a)
        self.u = u
        self.C = self.C_MKE(u, self.elements[:,1:], self.nodes[:,1:], zamik=True,circle=circle)
        return u
    
    def C_MKE(self, u, elements, nodes, circle='half',zamik=False):
        '''integracija po vsakem KE'''
        phi=0 #pretok
        if zamik: zamik=-1 #zakaj!?
        else: zamik=0
        for i in range(self.N_el):
            voz_index = elements[i]+zamik
            (x1,y1),(x2,y2),(x3,y3) = nodes[voz_index]
            A = 0.5 * np.abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
            U = u[voz_index]
            phi+= np.sum(U*A/3)
        R=eta=dp=1
        if circle=='full': A=np.pi*R**2 #krog
        elif circle=='half': A=0.5*np.pi*R**2 #polkrog
        elif circle=='hiska': A=5.
        else: print(' -> wrong area type!')
        S=A**2
        C = phi * 8*np.pi*eta / (S*dp)
        return C

    def meshing(self, file, Delaunay_mesh=False,scale=1):
        nodes, elements, boundary, mesher = self.read(file)
        n_node,x,y = np.transpose(nodes)
        boundary_nodes = boundary['boundary_nodes']
        self.boundary_nodes = boundary_nodes

        #scale
        x=x*scale
        y=y*scale
        nodes = np.transpose([n_node,x,y])
        self.nodes=nodes

        print('št. vozlišč:', self.N_voz)
        #print('št. elementov:',self.N_el)

        if Delaunay_mesh: #Delaunay meshing on same nodes as in program
            triang = Delaunay(np.transpose([x,y]))
            elementi = triang.simplices #indeksi vozlišč po KE, shape(št.el, št.vozlišč na KE=3)
            self.elementi = elementi
            if elementi.shape[0] == elements.shape[0] and (elementi==elements[:,-3:]-1).any():
                print('enaka mreža')
            #K,b = build_delaunay(elementi, nodes, boundary_nodes) #TODO
        else: 
            K,b = self.build(elements, nodes, boundary_nodes)
        return K,b

    def K_e3(self, vozlisca,A):
        '''trivozliščni KE
        
        (x1,y1),(x2,y2),(x3,y3) = vozlisca'''
        (x1,y1),(x2,y2),(x3,y3) = vozlisca
        #A = 0.5 * np.abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
        K = np.array([[(x2-x3)**2+(y2-y3)**2, -((x1-x3)*(x2-x3)+(y1-y3)*(y2-y3)), (x1-x2)*(x2-x3)+(y1-y2)*(y2-y3)],
                    [-((x1-x3)*(x2-x3)+(y1-y3)*(y2-y3)), (x1-x3)**2+(y1-y3)**2, -((x1-x2)*(x1-x3)+(y1-y2)*(y1-y3))],
                    [(x1-x2)*(x2-x3)+(y1-y2)*(y2-y3), -((x1-x2)*(x1-x3)+(y1-y2)*(y1-y3)), (x1-x2)**2+(y1-y2)**2]
                    ])
        return K * 1/(4*A)

    def build(self, elements, nodes, boundary_nodes, dtype=np.float64):
    
        #zlaganje
        K=np.zeros((self.N_voz,self.N_voz),dtype=dtype)
        b=np.zeros(self.N_voz)
        f = 1 #podatki - desna stran DE

        for i in range(self.N_el):
            voz_index = np.array(elements[i,1:]-1)
            voz = np.array(nodes[voz_index,1:])
            (x1,y1),(x2,y2),(x3,y3) = voz
            A = 0.5 * np.abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
            K[np.transpose([voz_index]), voz_index] += self.K_e3(voz,A)
            b[voz_index] += f*A/3

        ########### RP ###########
        #Dirichlet
        u_rp=0 #non_slip
        for j in boundary_nodes:
            i=j-1 #ansys -> boundary nodes
            K[i,:] = 0
            K[i,i] = 1
            b[i] = u_rp
            #plt.plot(nodes[i,1],nodes[i,2],'k.') #boundary nodes?
        return K,b

    def read(self, file):
        filetype = file.split('.')[-1]
        if filetype == 'dat':
            mesher = 'ANSYS'
            data = self.read_ansys(file)
        elif filetype == 'inp':
            mesher = 'ABAQUS'
            data = self.read_abaqus(file)
        else: 
            raise Exception('ERROR: Wrong file type')
        nodes, elements, boundary = data
        self.N_el  = elements.shape[0] #št. elementov
        self.N_voz = nodes.shape[0] #št. vozlišč
        self.nodes = nodes
        self.elements = elements
        self.boundary = boundary
        self.mesher = mesher
        return nodes, elements, boundary, mesher

    def read_abaqus(self, file):
        start_node = start_element = start_nset = start_elset = False
        nodes = []
        elements=[]
        nset = {}
        elset = {}
        with open(file) as dat:
            for line in dat:
                line = line.replace('\n', '').split(',')
                if '*NODE' in line[0].upper():
                    start_node = True
                elif '*ELEMENT' in line[0].upper():
                    start_node = False
                    start_element = True
                elif '*ELEMENT' in line[0].upper():
                    start_node = False
                    start_element = True
                elif '*END PART' in line[0].upper():
                    start_element = False
                elif '*NSET' in line[0].upper():
                    start_element = False
                    start_nset = True
                    for el in line:
                        el = el.split('=')
                        if 'NSET' in el[0].upper() and len(el)>1:
                            name = str(el[1])
                            nset[name] = np.array([])
                elif '*ELSET' in line[0].upper():
                    start_element = False
                    start_elset = True
                    for el in line:
                        el = el.split('=')
                        if 'ELSET' in el[0].upper() and len(el)>1:
                            name = str(el[1])
                            elset[name] = np.array([])
                elif '** SECTION' in line[0].upper():
                    start_nset = False
                    start_elset = False
                    #print(line)
                    #ne zanima
                elif '*END ASSEMBLY' in line[0].upper():
                    start_nset = False
                    start_elset = False
                elif start_node:
                    nodes.append(np.double(line)) #double precision
                elif start_element:
                    elements.append(np.int_(line))
                elif start_nset:
                    nset[name] = np.concatenate((nset[name].astype(int), np.int_(line).astype(int)))
                elif start_elset:
                    elset[name] = np.concatenate((nset[name], np.int_(line)), dtype=int)
        nodes = np.array(nodes)
        elements = np.array(elements, dtype=int)
        return nodes, elements, nset

    def read_ansys(self, file):
        start_node = start_element = start_boundary = np.inf
        nodes = []
        elements=[]
        boundary = {}
        boundary_nodes = []
        i=0
        with open(file) as dat:
            for line in dat:
                i+=1
                line = line.replace('\n', '').split(',')
                if '/com' in line[0] and len(line)>1:
                    if '*********** NODES FOR THE WHOLE ASSEMBLY ***********' in line[1].upper():
                        start_node = i
                    elif '*********** ELEMENTS FOR' in line[1].upper():
                        start_element = i
                    elif '*********** DEFINE TEMPERATURE CONSTRAINT ***********' in line[1].upper():
                        start_boundary = i
                elif '-1' == line[0]:
                    start_node = np.inf
                    start_element = np.inf
                    start_boundary = np.inf
                elif '' == line[0]:
                    start_boundary = np.inf
                elif i>start_node+2:
                    nodes.append(np.double(line[0].split('    ')[-4:-1])) #double precision, without z-axis
                elif i>start_element+4:
                    el = line[0].split(' ') #5 presledkov
                    el = [el_i for el_i in el if el_i != ''][-5:]
                    if el[-1] == el[-2]:
                        elements.append(np.int_(el[:-1])) #without fourth node
                    else:
                        elements.append(np.int_(el)) #all nodes
                elif i>start_boundary+2:
                    el = line[0].split('     ')[1:]
                    boundary_nodes.extend(np.int_(el))
                    if len(el) != 8:
                        start_boundary = np.inf
        nodes = np.array(nodes)
        elements = np.array(elements,dtype=int)
        boundary_nodes = np.array(boundary_nodes,dtype=int)
        boundary['boundary_nodes'] = boundary_nodes
        return nodes, elements, boundary
    
    def manual_meshing(self, Nr, N_phi, bias=1, dr_min=None, sparse=True, circle='half'):
        '''BIAS NE DELA!'''
        R=1#polmer
        dr=R/(Nr-1) #ni biasa
        Nr=Nr-1#sredina ni vsakič zraven
        # if dr_min==None:
        #     dr_min=R/(Nr-1)
        # if bias == 1:
        #     dr_max=R/(Nr-1)
        #     dr=np.ones(Nr)*dr_max
        #     if dr_max != dr_min:
        #         print('Warning: min korak in dimenzija nista usklajena!')
        # else: #bias!=1
        #     #dr_min=R/(Nr-1)
        #     dr_max = bias * dr_min
        #     dr = np.linspace(dr_max, dr_min, Nr)
        if circle == 'half':
            d_phi=np.pi/(N_phi-1)
        elif circle == 'full':
            d_phi=2*np.pi/(N_phi-1)
            N_phi=N_phi-1
        f = 1 #desna stran DE 'volumski izvor/ponor'
        #print('bias = ', bias, ', dr_max = ', round(dr_max,3), ', dr_min =', round(dr_min,3))

        #diskretizacija
        pointsall=np.zeros((Nr*N_phi+1,2))
        for j in range(N_phi):
            for i in range(1,Nr+1):
                index = (j)*(Nr)+i
                #edge_factor = 1 - abs(i - (Nr - 1) / 2) / ((Nr - 1) / 2)  # This will be 1 at the edge and decrease towards 0 towards the center
                #dr[i] = dr_min + (dr_max - dr_min) * edge_factor
                # xi = np.sum(dr[:i]) * np.cos(j*d_phi)
                # yi = -np.sum(dr[:i]) * np.sin(j*d_phi)
                xi = (dr*i) * np.cos(j*d_phi)
                yi = -(dr*i) * np.sin(j*d_phi)
                pointsall[index] = [xi,yi]
        pointsall[0] = [0,0]
        #triangulacija
        triang = Delaunay(pointsall)
        elementi = triang.simplices #indeksi vozlišč po KE, shape(št.el, št.vozlišč na KE=3)
        vozlisca = pointsall[triang.simplices] #koordinate vozlišč po KE, shape(št.el, št.vozlišč na KE=3, dimenzije=2)

        N_el = elementi.shape[0] #št. elementov
        N_voz = pointsall.shape[0] #št. vozlišč
        print('št. vozlišč:',N_voz)
        #print('št. elementov:',N_el)
        self.N_el = N_el
        self.N_voz = N_voz

        #zlaganje
        K=np.zeros((N_voz,N_voz))
        b=np.zeros(N_voz)

        for i in range(N_el):
            voz = vozlisca[i]
            voz_index = elementi[i]
            (x1,y1),(x2,y2),(x3,y3) = voz
            A = 0.5 * np.abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))
            K[np.transpose([voz_index]), voz_index] += self.K_e3(voz,A)
            b[voz_index] += f*A/3

        ########### RP ###########
        #Dirichlet
        u_rp=0 #non_slip
        if circle=='half':
            for i in range(1,Nr+1):
                #spodaj
                K[i,:] = 0
                K[i,i] = 1
                b[i] = u_rp
                #zgoraj
                K[-i,:] = 0
                K[-i,-i] = 1
                b[-i] = u_rp
            #središče
            K[0,:] = 0
            K[0,0] = 1
            b[0] = u_rp
        for j in range(0,N_phi):
            # if circle=='half':
            #     #znotraj
            #     K[j*Nr,:] = 0
            #     K[j*Nr,j*Nr] = 1
            #     b[j*Nr] = u_rp
            #zunaj
            K[j*Nr+Nr,:] = 0
            K[j*Nr+Nr,j*Nr+Nr] = 1
            b[j*Nr+Nr] = u_rp
        if sparse:
            u = spsolve(K,b)
        else:
            u = solve(K,b, assume_a='pos')
        self.u = u
        #self.u_2d = u.reshape((N_phi,Nr))
        self.pointsall = pointsall
        #self.nodes = vozlisca #/=
        self.elements = elementi
        if circle == 'full':
            self.C = self.C_MKE(u, elementi, pointsall,circle=circle)
        else:
            self.C = self.C_MKE(u, elementi, pointsall)
        return u