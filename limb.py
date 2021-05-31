import numpy as np
from scipy import interpolate
import os

#################################################################################
class Limb:

    def __init__(self, source, author="unknown"):

        self.Line = None
        self.LineReg = None
        self.author = author.lower()
        self.datapoints = []
        self.side = ''
        self.name = ''
        self.filename = source
        self.scale = 1.
        self.icp_score = 0
        self.day=0
        self.hour = 0
        self.ageAsString=''
        self.litterID=''
        self.embryoID = ''
        self.age=0
        self.ageIndex = None

        if isinstance(source, dict):
            self.author = source["author"]
            self.name = source["name"]
            self.filename = source["filename"]
            self.side = source["side"]
            self.datapoints = source["datapoints"]
            self.day = source["day"]
            self.hour = source["hour"]
            self.age = source["age"]
            self.ageAsString = source["ageAsString"]
            self.icp_score = source["icp_score"]
            self.litterID = source["litterID"]
            self.embryoID = source["embryoID"]
            return ################################### to be filled from npy file

        if 'RH' in source:
            self.side = 'R'
        elif 'LH' in source:
            self.side = 'L'
        elif 'right' in source:
            self.side = 'R'
        elif 'left' in source:
            self.side = 'L'
        else:
            raise RuntimeError()

        f = open(source, "r")
        lines = f.readlines()
        f.close()

        self.name = source.split('/')[-1]

        for i in range(len(lines)):
            if i<1: continue
            if author == "heura":
                line = lines[i].split()
            else:
                line = lines[i].split(',')
            self.datapoints.append([float(line[1]), float(line[2]), 0.0])

        self.datapoints = np.asarray(self.datapoints)

        if self.side == 'L':
            self.datapoints[:,0] *= -1

        if self.author == "heura":
            sfn = self.name.split('.')
            self.day = int(sfn[0].replace("E",""))
            self.hour = int(sfn[1].split('_')[0])
            self.ageAsString = str(self.day)+"."+sfn[1].split('_')[0]
            self.litterID = sfn[1].split('_')[1]
            self.embryoID = ""
        else:
            sfn = self.name.split('_')
            self.age = 0
            self.day = int(sfn[0].split(';')[0].replace("E",""))
            self.hour = int(sfn[0].split(';')[1])
            self.ageAsString = "E"+str(self.day)+"."+sfn[0].split(';')[1]
            self.litterID = sfn[1]
            self.embryoID = sfn[2]

        self.age = 24*self.day + self.hour


#################################################################################
def load_heura_limbs(stages=[14960, 15660, 16380, 17100, 17820],
                      mydir='/home/musy/Projects/limb_remeshing/coarse/',
                      c='k'):
    from vedo import load
    tclimbs = []
    for i,s in enumerate(stages):
        limb = load(mydir+'/*_'+str(s)+'min_coarse.vtk')
        limb.cutWithPlane(origin=(5,0,0))
        limb.scale(0.22) #############
        limb.lw(0.1).c(i).lighting('off')
        tclimbs.append(limb)
    return tclimbs

def load_heura_limbs2(source='/home/musy/Projects/stagingsystem/data/heura/'):
    limbs = []
    if source.endswith(".npy"):
        llist = np.load(source, allow_pickle=True)
        for dlimb in llist:
            if "heura" not in dlimb["author"]: continue
            lm = Limb(dlimb)
            limbs.append(lm)
    else:
        for fn in sorted(os.listdir(source)):
            if 'txt' not in fn: continue
            lm = Limb(source+fn, author="heura")
            limbs.append(lm)
    return limbs, ['10.00', '10.03', '10.06', '10.11', '10.14', '10.19']

#################################################################################
def load_james_limbs(source='/home/musy/Projects/stagingsystem/data/james/'):
    if source.endswith(".npy"):
        llist = np.load(source, allow_pickle=True)
        limbs = []
        for dlimb in llist:
            if "james" not in dlimb["author"]: continue
            lm = Limb(dlimb)
            limbs.append(lm)
    else:
        limbs = []
        for fn in sorted(os.listdir(source)):
            if '.csv' not in fn: continue
            lm = Limb(source+fn, author="james")
            limbs.append(lm)
    return limbs, ['10.09','10.21','11.09','11.21','12.09','12.21']

#################################################################################
def load_kevin_limbs(source='/home/musy/Projects/stagingsystem/data/kevin/'):
    limbs = []
    if source.endswith(".npy"):
        llist = np.load(source, allow_pickle=True)
        for dlimb in llist:
            if "kevin" not in dlimb["author"]: continue
            lm = Limb(dlimb)
            limbs.append(lm)
    else:
        for fn in sorted(os.listdir(source)):
            if '.csv' not in fn: continue
            lm = Limb(source+fn, author="kevin")
            limbs.append(lm)
    return limbs, ['12.09','13.09','14.09','14.21','15.09','15.21']

#################################################################################
def load_james_kevin_limbs(source='/home/musy/Projects/stagingsystem/data/james_kevin/'):
    limbs = []
    for fn in sorted(os.listdir(source)):
        if '.csv' not in fn: continue
        lm = Limb(source+fn, author="james_kevin")
        limbs.append(lm)
    return limbs, ['10.09','10.21','11.09','11.21','12.09','12.21',
                   '13.09','13.21','14.09','14.21','15.09','15.21']


#################################################################################
def load_marco_refs(stages=[249, 261, 273, 285, 297],
                    # mydir='/home/musy/Projects/tissuesim/data/timecourse1d/',
                    mydir='/home/musy/Projects/limb_remeshing/timecourse1d/',
                    c='k',
                    triangulate=True,
    ):
    from vedo import load
    tclimbs = []
    for s in stages:
        limb = load(mydir+'/*_'+str(s)+'.vtk').cutWithPlane()
        limb.scale(88) #############
        limb.c(c).alpha(0.5)
        pts = limb.points()
        limb.addPos(-(pts[0]+pts[-1])/2)
        if triangulate:
            ll = limb.triangulate()
        else:
            ll = limb
        tclimbs.append(ll.lighting('off'))
    return tclimbs

#################################################################################
def load_luciano_refs(stages=[14960, 15660, 16380, 17100, 17820],
                      mydir='/home/musy/Projects/limb_remeshing/coarse/',
                      c='k',
    ):
    from vedo import load
    tclimbs = []
    for i,s in enumerate(stages):
        limb = load(mydir+'/*_'+str(s)+'min_coarse.vtk')
        limb.cutWithPlane(origin=(5,0,0))
        limb.scale(0.22) #############
        limb.lw(0.1).c(i).lighting('off')
        tclimbs.append(limb)
    return tclimbs


####################################################
def getDist(source, target):
    a,b = source, target
    apts = a.points()
    bpts =[]
    for apt in apts:
        bpt = b.closestPoint(apt)
        bpts.append(bpt)
    bpts = np.array(bpts)
    dists = np.square(apts-bpts).sum(axis=1)
    return np.mean(dists)

def handedness(line):
    line = line.points()
    p0 = line[0]
    p1 = line[-1]
    pm = line[int(len(line)/2)]
    handedness = np.sign(np.cross(pm-p0, p1-p0)[2])
    return handedness

def ageAsString(agehour):
    day = int(agehour/24.0)
    h   = int(agehour) - 24*day
    if h<10: s = "0"
    else: s = ""
    return "E"+str(day)+"."+s+str(h)

def ageInHours(age_string, separator='.'):
    d, h = age_string.replace("E","").split(separator)
    d = int(d)
    h = int(h)
    return d*24+h

def spline_xy(lines):
    allx, ally, t = [], [], []
    NPT = lines[0].N()

    for i in range(NPT):
        for l in lines:
            l.c('k', 0.5).lw(2)
            p = l.points(i)
            allx.append(p[0])
            ally.append(p[1])
            t.append(float(i))

    ts = np.linspace(0, NPT, NPT, endpoint=False)
    knotstep=2
    knots = np.array(range(knotstep*2, NPT-knotstep*2, knotstep))
    #interpolate X
    tck   = interpolate.splrep(t, allx, s=0, k=3, t=knots, task=-1)
    allxs = interpolate.splev(ts, tck, der=0)
    #interpolate Y
    tck   = interpolate.splrep(t, ally, s=0, k=3, t=knots, task=-1)
    allys = interpolate.splev(ts, tck, der=0)
    return t, ts, allx, allxs, ally, allys

def reposition(lines, angle=0):
    # rigidly put lines contour about horizontal and their base at (0,0,0) (inplace)
    # without modifying their relative positioning
    p0s = []
    p1s = []
    pvs = []
    for l in lines:
        pts = l.points()
        p0s.append(pts[0])
        p1s.append(pts[-1])
        pvs.append(pts[int(l.N()/2)])

    p0 = np.mean(np.array(p0s), axis=0)
    p1 = np.mean(np.array(p1s), axis=0)
    pm = (p0+p1)/2
    pv = np.mean(np.array(pvs), axis=0)
    v = (pv-pm)/ np.linalg.norm(pv-pm)

    cr = np.cross(v, [0,1,0])[2]
    a = np.arccos(np.dot([0,1,0], v)) * 57.3 * cr/abs(cr)
    for l in lines:
        l.shift(-pm).rotateZ(-90 + a + angle)
    return


def distPD(line, isHeura=False):
    # proximo-distal distance of a limb shape
    if isinstance(line, np.ndarray):
        pts = line
    else:
        pts = line.points()
    p0 = pts[0]
    p1 = pts[-1]
    pm = (p0+p1)/2
    n = len(pts)
    if isHeura:
        return np.linalg.norm(pts[int(n/2)]-pm)
    maxd = 0
    for i in range(n):
        v = pts[i]-pm
        d = np.linalg.norm(v)
        if d>maxd:
            maxd=d
    return maxd
















