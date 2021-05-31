from limb import load_james_limbs, ageInHours
from vedo import Spline, Line, interactive, write, load, mag2
from vedo import Ribbon, merge, Axes, printc, precision, Video
from vedo import Text2D, DashedLine, Plotter, getColor, linInterpolate
from vedo.pyplot import plot
import numpy as np
import scipy.optimize as opt
import time


########################################################### load data
zscale = 10
limbs, stages = load_james_limbs("data/limbs_db.npy")
averages      = load("data/E*_james.vtk")
assert len(stages) == len(averages)

########################################################### limb tower
ave_sizes = []
for a in averages:
    a.lw(4)
    ave_sizes.append(a.averageSize())

averages_times = np.array([ageInHours(s)   for s in stages])
averages_sizes = np.array([a.averageSize() for a in averages])

growth = Line(averages_times, averages_sizes)

zlabs = []
for i,l in enumerate(averages):
    l.color(i).rotateZ(0)
    age = ageInHours(stages[i])
    l.z(age*zscale)
    zlabs.append([age*zscale, "E"+stages[i]+"\n("+str(age)+"h)"])

rbs = [] #### create the mesh by merging the ribbon strips
for i in range(len(averages) - 1):
    rb = Ribbon(averages[i], averages[i+1], closed=False, res=(200,12))
    rbs.append(rb)
mesh = merge(rbs).clean().color('k').alpha(0.2).pickable(True)

axs = Axes(mesh, zValuesAndLabels=zlabs, zLabelRotation=1, zTitleOffset=-0.03,
           xtitle='x \mum', ytitle='y \mum', ztitle='age', yzGridTransparent=True)

plt = Plotter(size=(1100,1350), interactive=0, title="Mouse Limb Fitter")
plt.camera.SetPosition( [534.088, -1335.679, 3762.633] )
plt.camera.SetFocalPoint( [153.867, -30.44, 2780.391] )
plt.camera.SetViewUp( [-0.108, 0.577, 0.809] )
plt.camera.SetDistance( 1677.205 )

# instructions message
def keycallb(keypressed):
    global skip
    if keypressed == 'i':
        if plt.interactive:
            plt.interactor.ExitCallback()
        plt.interactive = not(plt.interactive)
        printc('\t*** Interaction mode set to:', plt.interactive, '***', italic=1)
    elif keypressed == 's':
        skip = current_limb_age

# def fffunc(evt): # bug
#     print(evt)
# plt.addCallback("LeftButtonPress", fffunc)

plt.keyPressFunction = keycallb
plt.show(averages, mesh, axs, resetcam=False)

skip = ''
tinf = Text2D("Press:\ni to toggle interaction mode"+
              "\ns to jump to next stage"+
              "\nF1 to abort", font="Quikhand", pos='bottom-right', s=0.8)


########################################################### fitting
def fit_limb(limb):

    global ncalls, prev_score, pmin

    ###############
    def func(pars):

        global ncalls, prev_score, pmin, ln_output

        ncalls += 1

        ln = limb.Line.clone(transformed=True).c('k3').lw(5)
        ln.rotateZ(pars[3], locally=True).shift(pars[0], pars[1], pars[2])
        cpts = np.array([mesh.closestPoint(p) for p in ln.points()])
        lpts = np.array([ln.closestPoint(p)   for p in cpts]) # go back and forth
        d = cpts - lpts

        age_z = limb.age + pars[2]/zscale
        x = linInterpolate(age_z, [averages_times[0], averages_times[-1]], [0,1])
        gr_factor = growth.eval(x)[1]

        score = np.sum(mag2(d)) / gr_factor / ln.N() #*1000/ gr_factor
        ln_output = ln
        # print(d)
        # print(d*d)
        # exit()

        if score < prev_score or viz:
            prev_score = score
            minplot_data.append([ncalls, score])

            # --------- visualization:
            if (len(minplot_data) > 1 and score < minplot_data[-2][1]*0.95)  or viz:
                pmin = None
                if ncalls>9:
                    pmin = plot(minplot_data,
                                xtitle='#calls', ytitle='\epsilon',
                                lc='lb',
                                ylim=(0, minplot_data[0][1]*1.1),
                                aspect=16/9,
                                axes=dict(xyGridColor='k4',
                                          xyPlaneColor='k3',
                                          gridLineWidth=0.2,
                                          xyAlpha=1,
                                          textScale=2,
                                          htitle='SLSQP minimization',
                                          hTitleSize=0.04,
                                ),
                    )
                    pmin += DashedLine([0,15/gr_factor], [ncalls*1.1,15/gr_factor],
                                        c='r', lw=1, spacing=0.3)
                    pmin.rotateX(60, locally=True).scale(180/ncalls).pos(250,200,249*zscale)

                dage = precision(ln.z()/zscale, 2, vrange=20)+" h"
                pscale = "1.0 (fixed)"
                if len(pars)>4: pscale = precision(pars[4], 4)
                msg = ( "----  Fit parameters: ----"+
                        "\nCall #"+str(ncalls)+
                        "\nage = "+precision(age_z, 4)+" h, \Delta= "+dage+
                        "\nshift: "+precision(pars[:2], 2, 50)+
                        "\nangle: "+precision(pars[3], 2)+"\circ"+
                        "\nscale: "+pscale+
                        "\n\epsilon = "+precision(score, 3)+
                        "\n-------------------------"
                )
                plt.clear().show(averages, mesh, ln, *outlines, t2dl.text(msg), t2dr, tinf, pmin,
                                 axs, resetcam=False, interactive=False)
                # if plt.escaped: # ESC was pressed
                #     plt.close()
                #     exit(0)
                if vd: vd.addFrame().pause(0.1)
                if viz: # make the line flash in white
                    ln.c('w').lw(7); plt.render()
                    if vd: vd.addFrame().pause(0.1)
                    time.sleep(0.1)
                    ln.c('k3').lw(5); plt.render()
                    if vd: vd.addFrame().pause(0.1)

        return score

    # ----------------
    limb.Line.transformWithLandmarks(limb.Line,
                                     averages[limb.ageIndex].clone().z(limb.age*zscale),
                                     rigid=True)

    t2dr.text(limb.name.replace('.csv', '')+"\nAge group: "+limb.ageAsString)
    t2dr.background(limb.ageIndex, 0.9)

    # param are: x,y,z shifts, 1 angle, 1 scale
    x0    = [0,0,0, 0]                # initial guess, central values
    bnds  = [(-50, 50), (-50, 50)]    # x,y ranges
    bnds += [(-20*zscale, 20*zscale)] # z range
    bnds += [(-15, 15)]               # angle range

    ncalls = 0
    prev_score = 1e10
    viz, pmin= False, None
    minplot_data = []

    res = opt.minimize(func, x0, bounds=bnds, method="SLSQP", tol=0.0001)
    viz = True # rerun func for all pts:
    fscore = func(res["x"])
    limb.Line = ln_output
    printc(" Final fit score", precision(fscore,3), 'ncalls=', ncalls, c='k8')
    return fscore

# ##################################################################### FIT LOOP
vd = None # Video('slidingline.mp4')
pcloud, pcloud_cols, pcloud_ages, pcloud_scores = [],[],[],[]
t2dl = Text2D(font="VictorMono", c='k1', pos='top-left', s=0.8)
t2dr = Text2D(font="Calco", pos='top-right', c='w', s=1.05)

outlines = []
for i in range(110,len(limbs)):
    limb = limbs[i]
    limb.Line = Spline(limb.datapoints, res=200)
    limb.ageIndex = stages.index(limb.ageAsString.replace("E",""))
    current_limb_age = limb.ageAsString
    if skip == current_limb_age:
        continue

    col = getColor(limb.ageIndex)
    printc(i, "LIMB:", limb.ageAsString, str(limb.age)+'h', c=col, invert=1, end=' ')

    score = fit_limb(limb)

    limb.Line.c(col)
    limb.Line.info['score'] = score
    limb.Line.info['age'] = limb.age
    limb.Line.info['ageAsString'] = limb.ageAsString
    limb.Line.info['ageIndex'] = limb.ageIndex
    limb.Line.info['side'] = limb.side
    limb.Line.info['litterID'] = limb.litterID
    limb.Line.info['embryoID'] = limb.embryoID
    outlines.append(limb.Line.lw(1).alpha(0.5))

write(outlines, "outlines.npy")
if vd: vd.close()
interactive()







